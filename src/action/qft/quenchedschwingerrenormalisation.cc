#include "quenchedschwingerrenormalisation.hh"

/** @file quenchedschwingerrenormalisation.cc
 * @brief Implementation of quenchedschwingerrenormalisation.hh
 */

double RenormalisedQuenchedSchwingerParameters::betacoarse_nonperturbative() {
    int rho_refine = (coarsening_type==CoarsenBoth)?4:2;
    struct ParamType params = {beta,lattice->getNcells(),rho_refine};
    
    gsl_function f_gsl;
    f_gsl.function = &f_root;
    f_gsl.params = &params;
    
    // Status of root finder
    int status;
    // Maximal number of iterations
    unsigned int max_iter = 100;
    // Select a simple bisection solver
    const gsl_root_fsolver_type *SolverType;
    gsl_root_fsolver *solver;
    SolverType = gsl_root_fsolver_bisection;
    // Allocate memory
    solver = gsl_root_fsolver_alloc (SolverType);
    // Set lower and upper bound of search interval
    double x_lo = 0.01;
    double x_hi = 2.0;
    // relative tolerance
    double rel_tol = 1.E-12;
    
    // Root
    double x;
    
    // Initialise root finder
    double f_lo = f_root (x_lo, &params);
    double f_hi = f_root (x_hi, &params);
    if ( ( (f_lo>0) and (f_hi>0) ) or ( (f_lo<0) and (f_hi<0) )) {
        // Use fallback value if there is no root in the interval [x_lo,x_hi]
        x = coarsening_type==CoarsenBoth?0.25:0.5;
    } else {
        gsl_root_fsolver_set (solver, &f_gsl, x_lo, x_hi);
    
        // Set this to true to print output for debugging
        bool verbose = false;
        for (unsigned int k=0;k<max_iter;++k) {
            status = gsl_root_fsolver_iterate (solver);
            x = gsl_root_fsolver_root (solver);
            x_lo = gsl_root_fsolver_x_lower (solver);
            x_hi = gsl_root_fsolver_x_upper (solver);
            // Check for convergence
            status = gsl_root_test_interval (x_lo, x_hi, 0, rel_tol);
            if (verbose) {
                if (status == GSL_SUCCESS)
                    printf ("Converged:\n");

                printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                        k, x_lo, x_hi,x, x_hi - x_lo);
            }
            if (status != GSL_CONTINUE) break;
        }
    }
    gsl_root_fsolver_free(solver);
    return x*beta;
}
