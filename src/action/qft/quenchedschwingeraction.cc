#include "quenchedschwingeraction.hh"
/** @file quenchedschwingeraction.cc
 * @brief Implementation of quenchedschwingeraction.hh
 */

/* Value of action for a given configuration */
const double QuenchedSchwingerAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0;
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
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
        theta_p = mod_2pi(phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                        + phi_state->data[lattice->link_cart2lin(i,  j  ,1)]
                        - phi_state->data[lattice->link_cart2lin(i+1,j  ,1)]);
        theta_m = mod_2pi(phi_state->data[lattice->link_cart2lin(i,  j-1,0)]
                        + phi_state->data[lattice->link_cart2lin(i+1,j-1,1)]
                        - phi_state->data[lattice->link_cart2lin(i,  j-1,1)]);
    } else {
        theta_p = mod_2pi(phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                        + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                        - phi_state->data[lattice->link_cart2lin(i,  j+1,0)]);
        theta_m = mod_2pi(phi_state->data[lattice->link_cart2lin(i-1,j+1,0)]
                        + phi_state->data[lattice->link_cart2lin(i-1,j  ,1)]
                        - phi_state->data[lattice->link_cart2lin(i-1,j,  0)]);
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
    phi_state->data[ell] = exp_cos_dist.draw(engine,theta_p,theta_m);
}

/* local overrelaxation update */
void QuenchedSchwingerAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                                    const unsigned int ell) {
    int i, j, mu;
    double theta_p;
    double theta_m;
    lattice->link_lin2cart(ell,i,j,mu);
    compute_staple_angles(phi_state,i,j,mu,theta_p,theta_m);
    phi_state->data[ell] = mod_2pi((theta_p+theta_m) - phi_state->data[ell]);
}

/* Force for HMC integrator */
void QuenchedSchwingerAction::force(const std::shared_ptr<SampleState> phi_state,
                                    std::shared_ptr<SampleState> p_state) const {
    for (unsigned int ell=0;ell<p_state->data.size();++ell) {
        p_state->data[ell] = 0.0;
    }
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
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

/* Copy coarse links from state on coarser level */
void QuenchedSchwingerAction::copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                               std::shared_ptr<SampleState> phi_state) {
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    const unsigned int Mt_c_lat = coarse_lattice->getMt_lat();
    const unsigned int Mx_c_lat = coarse_lattice->getMx_lat();
    double theta_c;
    if ( ( Mt_c_lat == Mt_lat/2 ) and ( Mx_c_lat == Mx_lat/2 ) ) {
        /* case 1: coarsened in both directions */
        for (unsigned int i=0;i<Mt_c_lat;++i) {
            for (unsigned int j=0;j<Mx_c_lat;++j) {
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,0)];
                phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,0)] = 0.5*theta_c;
                phi_state->data[lattice->link_cart2lin(2*i+1,2*j  ,0)] = 0.5*theta_c;
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,1)];
                phi_state->data[lattice->link_cart2lin(2*i  ,2*j  ,1)] = 0.5*theta_c;
                phi_state->data[lattice->link_cart2lin(2*i  ,2*j+1,1)] = 0.5*theta_c;
            }
        }
    } else if ( ( Mt_c_lat == Mt_lat/2 ) and ( Mx_c_lat == Mx_lat ) ) {
        /* case 2: coarsened in temporal direction only */
        for (unsigned int i=0;i<Mt_c_lat;++i) {
            for (unsigned int j=0;j<Mx_lat;++j) {
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,0)];
                phi_state->data[lattice->link_cart2lin(2*i  ,j  ,0)] = 0.5*theta_c;
                phi_state->data[lattice->link_cart2lin(2*i+1,j  ,0)] = 0.5*theta_c;
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,1)];
                phi_state->data[lattice->link_cart2lin(2*i  ,j  ,1)] = theta_c;
            }
        }
    } else if ( ( Mt_c_lat == Mt_lat ) and ( Mx_c_lat == Mx_lat/2 ) ) {
        /* case 3: coarsened in spatial direction only */
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_c_lat;++j) {
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,0)];
                phi_state->data[lattice->link_cart2lin(i  ,2*j  ,0)] = theta_c;
                theta_c = phi_coarse->data[coarse_lattice->link_cart2lin(i  ,j  ,1)];
                phi_state->data[lattice->link_cart2lin(i  ,2*j  ,1)] = 0.5*theta_c;
                phi_state->data[lattice->link_cart2lin(i  ,2*j+1,1)] = 0.5*theta_c;
            }
        }
    } else {
        mpi_parallel::cerr << "ERROR: cannot copy from coarse lattice." << std::endl;
        mpi_exit(EXIT_FAILURE);
        throw std::runtime_error("...");
    }
}

/* Copy coarse links from state on finer level */
void QuenchedSchwingerAction::copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                                             std::shared_ptr<SampleState> phi_state) {
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    const unsigned int Mt_f_lat = fine_lattice->getMt_lat();
    const unsigned int Mx_f_lat = fine_lattice->getMx_lat();
    if ( ( Mt_f_lat == 2*Mt_lat ) and ( Mx_f_lat == 2*Mx_lat ) ) {
        /* case 1: refined in both directions */
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_lat;++j) {
                phi_state->data[lattice->link_cart2lin(i  ,j  ,0)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(2*i  ,2*j  ,0)]
                            + phi_fine->data[fine_lattice->link_cart2lin(2*i+1,2*j  ,0)]);
                phi_state->data[lattice->link_cart2lin(i  ,j  ,1)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(2*i  ,2*j  ,1)]
                            + phi_fine->data[fine_lattice->link_cart2lin(2*i  ,2*j+1,1)]);
            }
        }
    } else if ( ( Mt_f_lat == 2*Mt_lat ) and ( Mx_f_lat == Mx_lat ) ) {
        /* case 2: refined in temporal direction only */
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_lat;++j) {
                phi_state->data[lattice->link_cart2lin(i  ,j  ,0)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(2*i  ,j  ,0)]
                            + phi_fine->data[fine_lattice->link_cart2lin(2*i+1,j  ,0)]);
                phi_state->data[lattice->link_cart2lin(i  ,j  ,1)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(2*i  ,j  ,1)]);
            }
        }
    } else if ( ( Mt_f_lat == Mt_lat ) and ( Mx_f_lat == 2*Mx_lat ) ) {
        /* case 3: refined in spatial direction only */
        for (unsigned int i=0;i<Mt_lat;++i) {
            for (unsigned int j=0;j<Mx_lat;++j) {
                phi_state->data[lattice->link_cart2lin(i  ,j  ,0)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(i  ,2*j  ,0)]);
                phi_state->data[lattice->link_cart2lin(i  ,j  ,1)]
                    = mod_2pi(phi_fine->data[fine_lattice->link_cart2lin(i  ,2*j  ,1)]
                            + phi_fine->data[fine_lattice->link_cart2lin(i  ,2*j+1,1)]);
            }
        }
    } else {
        mpi_parallel::cerr << "ERROR: cannot copy from fine lattice." << std::endl;
        mpi_exit(EXIT_FAILURE);
        throw std::runtime_error("...");
    }
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

/* Return lattice information */
std::string QuenchedSchwingerAction::info_string() const {
    std::stringstream sstr;
    sstr << QFTAction::info_string() << ", beta = " << beta;
    return sstr.str();
}
