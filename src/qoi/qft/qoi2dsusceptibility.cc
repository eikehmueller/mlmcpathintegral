#include "qoi2dsusceptibility.hh"
/** @file qoi2dsusceptibility.cc
 * @brief Implementation of qoi2dsusceptibility.hh
 */

/* Evaluate QoI */
const double QoI2DSusceptibility::evaluate(const std::shared_ptr<SampleState> phi_state) {
    if (phi_state->data.size() != 2*Mt_lat*Mx_lat) {
        mpi_parallel::cout << "ERROR: Evaluating QoI2DSusceptibility on state of wrong size." << std::endl;
        mpi_exit(EXIT_FAILURE);
    }
    // lambda function for working out linear index of link
    double Q = 0.0;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            double theta = phi_state->data[lattice->link_cart2lin(i,  j,  0)]
                         + phi_state->data[lattice->link_cart2lin(i+1,j,  1)]
                         - phi_state->data[lattice->link_cart2lin(i,  j+1,0)]
                         - phi_state->data[lattice->link_cart2lin(i,  j  ,1)];
            Q += mod_2pi(theta);
        }
    }
    return four_pi2_inv*Q*Q;
}

/* Analytical expression for topological susceptibility */
double quenchedschwinger_chit_analytical(const double beta,
                                         const unsigned int n_plaq) {
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
        double rho = In[n]/In[0];
        weight[n] = duplicity*pow(rho,n_plaq);
        weight_sum += weight[n];
    }
    double chit=0.0;
    for (int n=0;n<nmax;++n) {
        chit += n_plaq*weight[n]/weight_sum * (ddIn[n]/In[n] -(n_plaq-1)*(dIn[n]*dIn[n])/(In[n]*In[n]));
    }
    return chit;
}

/* Compute functions I'_n(x) and I''_n(x) required in calculation of topological susceptibility */
void compute_In(const double x,
                std::vector<double>& In,
                std::vector<double>& dIn,
                std::vector<double>& ddIn) {
    
    /* Integrand -1/(4*pi^2)*phi*exp(-x*cos(phi)) */
    auto integrand_phi1 = [](double phi, void * p) -> double {
        double x = *((double*) p);
        return  -1./(4.*M_PI*M_PI)*phi*exp(x*(cos(phi)-1.0));
    };

    /* Integrand 1/(8*pi^3)*phi^2*exp(-x*cos(phi)) */
    auto integrand_phi2 = [](double phi, void * p) -> double {
        double x = *((double*) p);
        return  1./(8.*M_PI*M_PI*M_PI)*phi*phi*exp(x*(cos(phi)-1.0));
    };

    /* GSL function required for integration */
    gsl_function integrand;
    integrand.params = (void*) &x;
    
    /* GSL workspace required for numerical integration with QAWO routine */
    const size_t n_level = 10; // 2^{n_level} must not exceed n_workspace
    const size_t n_workspace = 1<<n_level;
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
        In[n] = gsl_sf_bessel_In_scaled(n,x);
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
