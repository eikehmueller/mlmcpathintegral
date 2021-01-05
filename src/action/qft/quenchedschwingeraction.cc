#include "quenchedschwingeraction.hh"
/** @file quenchedschwingeraction.cc
 * @brief Implementation of quenchedschwingeraction.hh
 */

/* Value of action for a given configuration */
const double QuenchedSchwingerAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta =  phi_state->data[lattice->link_cart2lin(i  ,j  ,0)]
                          + phi_state->data[lattice->link_cart2lin(i+1,j  ,1)]
                          - phi_state->data[lattice->link_cart2lin(i  ,j+1,0)]
                          - phi_state->data[lattice->link_cart2lin(i  ,j  ,1)];
            S += (1.-cos(theta));
        }
    }
    return beta*S;
}

/* Compute staple angles */
void QuenchedSchwingerAction::compute_staple_angles(std::shared_ptr<SampleState> phi_state,
                                                    const int i,
                                                    const int j,
                                                    const int mu,
                                                    double& theta_p,
                                                    double& theta_m) {
    if (mu==0) {
        theta_p = phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                + phi_state->data[lattice->link_cart2lin(i,  j  ,1)]
                - phi_state->data[lattice->link_cart2lin(i+1,j  ,1)];
        theta_m = phi_state->data[lattice->link_cart2lin(i,  j-1,0)]
                + phi_state->data[lattice->link_cart2lin(i+1,j-1,1)]
                - phi_state->data[lattice->link_cart2lin(i,  j-1,1)];
    } else {
        theta_p = phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                - phi_state->data[lattice->link_cart2lin(i,  j+1,0)];
        theta_m = phi_state->data[lattice->link_cart2lin(i-1,j+1,0)]
                + phi_state->data[lattice->link_cart2lin(i-1,j  ,1)]
                - phi_state->data[lattice->link_cart2lin(i-1,j,  0)];
    }
}

/* local heat bath update */
void QuenchedSchwingerAction::heatbath_update(std::shared_ptr<SampleState> phi_state,
                                              const unsigned int ell) {
    int i, j, mu;
    double theta_p;
    double theta_m;
    lattice->link_lin2cart(ell,i,j,mu);
    compute_staple_angles(phi_state,i,j,mu,theta_p,theta_m);
    double theta_0 = 0.5*(theta_p+theta_m)+(fabs((theta_p-theta_m))>M_PI)*M_PI;
    double sigma = 2.*beta*fabs(cos(0.5*(theta_p-theta_m)));
    phi_state->data[ell] = mod_2pi(theta_0 + exp_sin2_dist.draw(engine,sigma));
}

/* local overrelaxation update */
void QuenchedSchwingerAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                                    const unsigned int ell) {
    int i, j, mu;
    double theta_p;
    double theta_m;
    lattice->link_lin2cart(ell,i,j,mu);
    compute_staple_angles(phi_state,i,j,mu,theta_p,theta_m);
    double theta_0 = 0.5*(theta_p+theta_m)+(fabs((theta_p-theta_m))>M_PI)*M_PI;
    phi_state->data[ell] = mod_2pi(2.0*theta_0 - phi_state->data[ell]);
}

/* Force for HMC integrator */
void QuenchedSchwingerAction::force(const std::shared_ptr<SampleState> phi_state,
                                    std::shared_ptr<SampleState> p_state) const {
    for (unsigned int ell=0;ell<p_state->data.size();++ell) {
        p_state->data[ell] = 0.0;
    }
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta = phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                         + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                         - phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                         - phi_state->data[lattice->link_cart2lin(i,  j,  1)];
            double F = beta*sin(theta);
            p_state->data[lattice->link_cart2lin(i,  j,  0)] += F;
            p_state->data[lattice->link_cart2lin(i+1,j,  1)] += F;
            p_state->data[lattice->link_cart2lin(i,  j+1,0)] -= F;
            p_state->data[lattice->link_cart2lin(i,  j  ,1)] -= F;
        }
    }
}

/* Copy coarse data from path on coarser level */
void QuenchedSchwingerAction::copy_from_coarse(const std::shared_ptr<SampleState> x_coarse,
                                               std::shared_ptr<SampleState> x_path) {
    mpi_parallel::cerr << "ERROR: QuenchedSchwingerAction::copy_from_coarse not implemented yet" << std::endl;
    mpi_exit(EXIT_FAILURE);
}

/* Copy coarse data from path on finer level */
void QuenchedSchwingerAction::copy_from_fine(const std::shared_ptr<SampleState> x_fine,
                                             std::shared_ptr<SampleState> x_path) {
    mpi_parallel::cerr << "ERROR: QuenchedSchwingerAction::copy_from_fine not implemented yet" << std::endl;
    mpi_exit(EXIT_FAILURE);
}

/* Initialise state with random entries */
void QuenchedSchwingerAction::initialise_state(std::shared_ptr<SampleState> phi_state) const {
    std::uniform_real_distribution<double> uniform(-M_PI,M_PI);
    std::generate(phi_state->data.data(),
                  phi_state->data.data()+phi_state->data.size(),
    [this,&uniform]() {
        return uniform(engine);
    });
}

/* Check whether coarsening is permitted */
void QuenchedSchwingerAction::check_coarsening_is_permitted(const unsigned int n_level) {
    if ( (Mt_lat>>n_level)<<n_level == Mt_lat) {
        mpi_parallel::cout << "M_{t,lat} = " << Mt_lat << " = 2^{" << n_level << "-1} * " << (Mt_lat>>(n_level-1)) << std::endl;
    } else {
        mpi_parallel::cout << "ERROR: M_{t,lat} = " << Mt_lat << " is not a multiple of 2^{n_level} = 2^{"<<n_level << "}" << std::endl;
        mpi_exit(-1);
    }
    if ( (Mx_lat>>n_level)<<n_level == Mx_lat) {
        mpi_parallel::cout << "M_{x,lat} = " << Mx_lat << " = 2^{" << n_level << "-1} * " << (Mx_lat>>(n_level-1)) << std::endl;
    } else {
        mpi_parallel::cout << "ERROR: M_{x,lat} = " << Mx_lat << " is not a multiple of 2^{n_level} = 2^{"<<n_level << "}" << std::endl;
        mpi_exit(-1);
    }
}
