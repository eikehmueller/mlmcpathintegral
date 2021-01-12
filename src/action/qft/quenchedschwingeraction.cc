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
    double theta0 = theta0_expsin2(theta_p,theta_m);
    double sigma = sigma_expsin2(theta_p,theta_m);
    double dtheta = exp_sin2_dist.draw(engine,beta*sigma);
    phi_state->data[ell] = mod_2pi(theta0 + dtheta);
}

/* local overrelaxation update */
void QuenchedSchwingerAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                                    const unsigned int ell) {
    int i, j, mu;
    double theta_p;
    double theta_m;
    lattice->link_lin2cart(ell,i,j,mu);
    compute_staple_angles(phi_state,i,j,mu,theta_p,theta_m);
    double theta0 = theta0_expsin2(theta_p,theta_m);
    phi_state->data[ell] = mod_2pi(2.0*theta0 - phi_state->data[ell]);
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

/* Exact result for topological susceptibility */
double QuenchedSchwingerAction::chit_exact() const {
    /* Number of plaquettes */
    int n_plaq = Mt_lat*Mx_lat;
    /* Truncation index of sums */
    const int nmax = 10;
    /* Compute functions I_n(\beta), I'_n(\beta), I''_n(\beta) */
    std::vector<double> In(nmax);
    std::vector<double> dIn(nmax);
    std::vector<double> ddIn(nmax);
    compute_In(beta,In,dIn,ddIn);
    /* Un-normalised weights */
    std::vector<double> weight(nmax);
    /* Sum of weights needed for normalisation */
    double weight_sum = 0.0;
    /* Compute weights w_n(\beta,P) */
    for (int n=0;n<nmax;++n) {
        /* Terms with n>0 are counted twice to account for
           corresponding negative entries in original sum */
        double duplicity = (1 + (n > 0));
        weight[n] = duplicity*pow(In[n],n_plaq);
        weight_sum += weight[n];
    }
    double chit=0.0;
    for (int n=0;n<nmax;++n) {
        chit += weight[n]/weight_sum * (ddIn[n]/In[n] -(n_plaq-1)*(dIn[n]*dIn[n])/(In[n]*In[n]));
    }
    return chit;
}

/* Compute functions I'_n(x) and I''_n(x) required in calculation of topological susceptibility */
void QuenchedSchwingerAction::compute_In(const double x,
                                         std::vector<double>& In,
                                         std::vector<double>& dIn,
                                         std::vector<double>& ddIn) const {
    
    /* Integrand -1/(4*pi^2)*phi*exp(-x*cos(phi)) */
    auto integrand_phi1 = [](double phi, void * p) -> double {
        double x = *((double*) p);
        return  -1./(4.*M_PI*M_PI)*phi*exp(x*cos(phi));
    };

    /* Integrand 1/(8*pi^3)*phi^2*exp(-x*cos(phi)) */
    auto integrand_phi2 = [](double phi, void * p) -> double {
        double x = *((double*) p);
        return  1./(8.*M_PI*M_PI*M_PI)*phi*phi*exp(x*cos(phi));
    };

    /* GSL function required for integration */
    gsl_function integrand;
    integrand.params = (void*) &x;
    
    /* GSL workspace required for numerical integration with QAWO routine */
    const size_t n_workspace = 1024;
    const size_t n_level = 10; // 2^{n_level} must not exceed n_workspace
    gsl_integration_workspace* workspace;
    workspace = gsl_integration_workspace_alloc(n_workspace);
    gsl_integration_qawo_table* qawo_table;
    qawo_table = gsl_integration_qawo_table_alloc(1.0,2.0*M_PI,GSL_INTEG_COSINE,n_level);
    
    /* Tolerances for numerical integration */
    const double epsabs = 1.E-10;
    const double epsrel = 1.E-12;
    
    /* Evaluate all functions*/
    for (int n=0;n<In.size();++n) {
        /* --- I_n(x) --- */
        In[n] = gsl_sf_bessel_In(n,x);
        /* Evaluate the integrals */
        double abserr; // Absolute error of numerical integration
        /* --- I'_n(x), weight function is sin(n*phi) --- */
        gsl_integration_qawo_table_set(qawo_table,n,2.0*M_PI,GSL_INTEG_SINE);
        integrand.function = integrand_phi1;
        gsl_integration_qawo(&integrand, -M_PI, epsabs, epsrel, n_workspace, workspace, qawo_table, &dIn[n], &abserr);
        /* --- I''_n(x), weight function is cos(n*phi) --- */
        gsl_integration_qawo_table_set(qawo_table,n,2.0*M_PI,GSL_INTEG_COSINE);
        integrand.function = integrand_phi2;
        gsl_integration_qawo(&integrand, -M_PI, epsabs, epsrel, n_workspace, workspace, qawo_table, &ddIn[n], &abserr);
    }
    /* Free woskspace memory */
    gsl_integration_qawo_table_free(qawo_table);
    gsl_integration_workspace_free(workspace);
}

