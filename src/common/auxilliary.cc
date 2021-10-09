#include "auxilliary.hh"
/** @file auxilliary.cc
 * @brief Implementation of auxilliary.hh
 */

/* Implementation of Sigma_hat */
double Sigma_hat(const double xi,
                 const unsigned int p) {
    if (p % 2 == 0) {
        if (p == 0) {
            // Ratio is one if p==0
            return 1.0;
        } else {
            unsigned int mmax = 100; // Number of terms in sum
            double num = 0.0;
            double denom = 1.0;
            for (unsigned int m=1; m<mmax; ++m) {
                double exp_factor = exp(-0.5*xi*m*m);
                num += 2.*pow(m,p)*exp_factor;
                denom += 2.*exp_factor;
            }
            return num/denom;
        }
    } else {
        // Sum is zero if p is odd (symmetry)
        return 0.0;
    }
}

/* Implementation of log_factorial */
double log_factorial(unsigned int n) {
  double s=0.0;
  for (unsigned int k=2;k<=n;++k) {
    s += log(k);
  }
  return s;
}

/* Implementation of log_nCk */
double log_nCk(unsigned int n, unsigned int k) {
 return log_factorial(n)-log_factorial(k)-log_factorial(n-k);
}

/* Analytical expression for susceptibility function Phi_chit(beta,P) */
double Phi_chit(const double beta,
                const unsigned int n_plaq) {
    if (beta > 2000.0) {
        // The computation becomes unstable for large values of beta
        mpi_parallel::cerr << "ERROR: Phi_chit(beta,P) unstable for beta>2000. Use Phi_chit_perturbative instead." << std::endl;
        mpi_exit(EXIT_FAILURE);
    }
    /* Truncation index of sums */
    const int nmax = 20;
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
    double phi_chit=0.0;
    for (int n=0;n<nmax;++n) {
        phi_chit += beta*weight[n]/weight_sum * (ddIn[n]/In[n] -(n_plaq-1)*(dIn[n]*dIn[n])/(In[n]*In[n]));
    }
    return phi_chit;
}

/* Perturbative approximation for susceptibility function Phi_chit(beta,P) */
double Phi_chit_perturbative(const double beta,
                             const unsigned int n_plaq) {
    double xi = n_plaq/beta;
    double z = 1./beta;
    double Sigma_hat2 = Sigma_hat(xi,2);
    double Sigma_hat4 = Sigma_hat(xi,4);
    // Leading order
    double Phi_LO = 1.0-xi*Sigma_hat2;
    // Next-to-leading order coefficient
    double Phi_NLO = 0.5-xi*Sigma_hat2+0.25*xi*xi*(Sigma_hat4-Sigma_hat2*Sigma_hat2);
    // Total result
    return (Phi_LO + z*Phi_NLO)/(4.*M_PI*M_PI);
}


/* Compute functions I'_n(x) and I''_n(x) required in calculation of topological 
 * susceptibility function \f$\Phi_{\chi_t(\beta,P)}\f$ */
void compute_In(const double x,
                std::vector<double>& In,
                std::vector<double>& dIn,
                std::vector<double>& ddIn) {
        
    /* Structures for parameters to be passed to the integrand */
    struct ParamType {
        double x;      // argument x
        int n;         // order n
        bool use_qawo; // Use the QAWO integration method for trig-functions?
    };
    
    /* Integrand -1/(4*pi^2)*phi*exp(-x*cos(phi)) */
    auto integrand_phi1 = [](double phi, void * p) -> double {
        struct ParamType *params = (struct ParamType *) p;
        double x = params->x;
        int n = params->n;
        bool include_trig = not params->use_qawo;
        double f = -1./(4.*M_PI*M_PI)*phi*exp(x*(cos(phi)-1.0));
        if (include_trig) f *= sin(n*phi);
        return f;
    };

    /* Integrand 1/(8*pi^3)*phi^2*exp(-x*cos(phi)) */
    auto integrand_phi2 = [](double phi, void * p) -> double {
        struct ParamType *params = (struct ParamType *) p;
        double x = params->x;
        int n = params->n;
        bool include_trig = not params->use_qawo;
        double f = 1./(8.*M_PI*M_PI*M_PI)*phi*phi*exp(x*(cos(phi)-1.0));
        if (include_trig) f *= cos(n*phi);
        return f;
    };

    /* GSL function required for integration */
    gsl_function integrand;
    ParamType params;
    params.x = x;
    /*
     * Only use QAWO integration for not too large values of x.
     * For large values of x the function is concentrated near zero and Gaussian
     * quadrature is more robust.
     */
    params.use_qawo = (x < 32.0);
    integrand.params = &params;
    
    /* GSL workspace required for numerical integration with QAWO routine */
    const size_t n_level = 20; // 2^{n_level} must not exceed n_workspace
    const size_t n_workspace = 1<<n_level;
    gsl_integration_workspace* workspace;
    workspace = gsl_integration_workspace_alloc(n_workspace);
    gsl_integration_qawo_table* qawo_table;
    if (params.use_qawo) {
        qawo_table = gsl_integration_qawo_table_alloc(1.0,2.0*M_PI,GSL_INTEG_COSINE,n_level);
    }
    
    /* Tolerances for numerical integration */
    const double epsabs = 1.E-15;
    const double epsrel = 1.E-12;
    
    /* Evaluate all functions*/
    for (int n=0;n<In.size();++n) {
        /* --- I_n(x) --- */
        In[n] = gsl_sf_bessel_In_scaled(n,x);
        /* Evaluate the integrals */
        double abserr; // Absolute error of numerical integration
        params.n = n;
        /* --- I'_n(x), weight function is sin(n*phi) --- */
        integrand.function = integrand_phi1;
        if (params.use_qawo) {
            gsl_integration_qawo_table_set(qawo_table,n,2.0*M_PI,GSL_INTEG_SINE);
            gsl_integration_qawo(&integrand, -M_PI, epsabs, epsrel, n_workspace, workspace, qawo_table, &dIn[n], &abserr);
        } else {
            gsl_integration_qag(&integrand, -M_PI, +M_PI, epsabs, epsrel, n_workspace, 1, workspace,&dIn[n], &abserr);
        }
        /* --- I''_n(x), weight function is cos(n*phi) --- */
        integrand.function = integrand_phi2;
        if (params.use_qawo) {
            gsl_integration_qawo_table_set(qawo_table,n,2.0*M_PI,GSL_INTEG_COSINE);
            gsl_integration_qawo(&integrand, -M_PI, epsabs, epsrel, n_workspace, workspace, qawo_table, &ddIn[n], &abserr);
        } else {
            gsl_integration_qag(&integrand, -M_PI, +M_PI, epsabs, epsrel, n_workspace, 1, workspace,&ddIn[n], &abserr);
        }
    }
    /* Free woskspace memory */
    if (params.use_qawo) {
        gsl_integration_qawo_table_free(qawo_table);
    }
    gsl_integration_workspace_free(workspace);
}

/* exact expectation value of phi^2 for GFF */
double gff_phi_squared_analytical(const double mu2,
                                  const double Mt_lat,
                                  const double Mx_lat) {
    double sum_spectral = 0.0;
    for (unsigned int k1=0;k1<Mt_lat;++k1) {
        for (unsigned int k2=0;k2<Mx_lat;++k2) {
            double sin_k1 = sin(0.5*M_PI*k1/Mt_lat);
            double sin_k2 = sin(0.5*M_PI*k2/Mx_lat);
            sum_spectral += 1./(4.*(sin_k1*sin_k1+sin_k2*sin_k2)+mu2);
        }
    }
    return sum_spectral/(Mt_lat*Mx_lat);
}
