#ifndef AUXILLIARY_HH
#define AUXILLIARY_HH AUXILLIARY_HH
/** @file auxilliary.hh
 *
 * @brief Several auxilliary functions
 */
#include "config.h"
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include "mpi/mpi_wrapper.hh"


/* Generate coloured output? This is pretty, but will add escape
 * sequences if the output is precessed with a program which does not
 * support this
 */
#ifdef USECOLOR
#define RST "\x1B[0m"
#define CBOLD "\x1B[1m"
#define CRED "\x1B[31m"
#define CBLUE "\x1B[34m"
#define CGREEN "\x1B[32m"
#define CMAGENTA "\x1B[35m"
#define FBOLD(X) CBOLD X RST
#define FRED(X) CRED X RST
#define FBLUE(X) CBLUE X RST
#define FGREEN(X) CGREEN X RST
#define FMAGENTA(X) CMAGENTA X RST
#else
#define FBOLD(X) X
#define FRED(X) X
#define FBLUE(X) X
#define FGREEN(X) X
#define FMAGENTA(X) X
#endif // USECOLOR

/** @brief Calculate \f$ x mod [-pi,pi) \f$
 *
 * @param[in] x Value of \f$x\f$
 */
double inline mod_2pi(const double x) {
    return x - 2.*M_PI*floor(0.5*(x+M_PI)/M_PI);
}

/** @brief Calculate \f$ x mod [-pi/2,pi/2) \f$
 *
 * @param[in] x Value of \f$x\f$
 */
double inline mod_pi(const double x) {
    return x - M_PI*floor((x+0.5*M_PI)/M_PI);
}


/** @brief Calculate Sum \f$\hat{\Sigma_p(\xi)}\f$
 *
 * Return value of sum
 *
 \f[
 \hat{\Sigma}_p(\xi) = \frac{\sum_{m\in\mathbb{Z}}m^p \exp\left[-\frac{1}{2}\xi m^2\right]}{\sum_{m\in\mathbb{Z}}\exp\left[-\frac{1}{2}\xi m^2\right]}
 \f]
 *
 * @param[in] xi value of \f$\xi\f$
 * @param[in] p power \f$p\f$
 */
double Sigma_hat(const double xi,
                 const unsigned int p);

/** @brief Compute logarithm of factorial
 *
 * Returns \f$\log(n!)\f$
 *
 * @param[in] n Parameter \f$n\f$
 */
double log_factorial(unsigned int n);

/** @brief Compute logarithm of n choose k
 *
 * Returns \f$\log(n!)-\log(k!)-\log((n-k)!)\f$
 *
 * @param[in] n Parameter \f$n\f$
 * @param[in] k Parameter \f$k\f$
 */
double log_nCk(unsigned int n, unsigned int k);

/** @brief Analytical result for function \f$\Phi_{\chi_t}\f$ required in topological
 * susceptibility calculation.
 *
 * Computes the analytical value of the function \f$\Phi_{\chi_t}\f$ which is required in the
 * calculation of the analytical expression for the topological susceptibility.
 * If the number of plaquettes is \f$P\f$, then this is given by
 *
 * \f[
 *   \Phi_{\chi_t}(\beta,P)
  *    = \beta \sum_{n=-\infty}^{\infty} w_n(\beta,P)\left[
 *      (P-1) \left(\frac{I'_n(\beta)}{I_n(\beta)}\right)^2 - \frac{I''_n(\beta)}{I_n(\beta)}
 *   \right]
 * \f]
 *  with the weights
 * \f[
 *  w_n(\beta,P) = \frac{(I_n(\beta)^P)}{\sum_{n=-\infty}^{\infty} (I_n(\beta))^P}
 * \f]
 *  and the functions
 *  \f[
 *  \begin{aligned}
 *    I_n(x) &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= \frac{1}{2\pi} \int_{-\pi}^{+\pi} e^{x(\cos(\phi)-1)} \cos(n\phi) \; d\phi \\
 *    I'_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \frac{\phi}{2\pi }e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= -\frac{1}{4\pi^2} \int_{-\pi}^{+\pi} \phi e^{x(\cos(\phi)-1)} \sin(n\phi) \; d\phi \\
 *    I''_n(x) &= \frac{i}{2\pi} \int_{-\pi}^{+\pi} \left(\frac{\phi}{2\pi }\right)^2 e^{i n\phi + x(\cos(\phi)-1)} \; d\phi \\
 *        &= -\frac{1}{8\pi^3} \int_{-\pi}^{+\pi} \phi^2 e^{x(\cos(\phi)-1)} \cos(n\phi) \; d\phi
 *  \end{aligned}
 *  \f]
 *  
 *  The topological susceptibility (scaled by the lattice volume)of the quenched Schwinger
 *  model is then given by
 * 
 * \f[
 *    V \chi_t(\beta,P) = P/\beta \Phi_{\chi_t}(\beta,P)
 * \f]
 * 
 *  All functions are scaled by a factor \f$e^{-x}\f$ for numerical stability (this factor will cancel out, since we
 *  only every compute ratios of the above functions). Note that \f$I_n(x)\f$ is the (rescaled) modified Bessel
 *  function of the first kind. For more details see the following two papers:
 *
 *  - Bonati, C. and Rossi, P., 2019. "Topological susceptibility of two-dimensional U (N) gauge theories."
 *   Physical Review D, 99(5), p.054503 https://doi.org/10.1103/PhysRevD.99.054503
 *
 *  - Kiskis, J., Narayanan, R. and Sigdel, D., 2014. "Correlation between Polyakov loops oriented in two different
 *   directions in SU(N) gauge theory on a two-dimensional torus. Physical Review D, 89(8), p.085031.
 *   https://doi.org/10.1103/PhysRevD.89.085031
 *
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double Phi_chit(const double beta,
                const unsigned int n_plaq);

/** @brief Perturbative approximation to \f$\Phi_{\chi_t}\f$ up to (and including) terms of 
 * \f$O(1/\beta)\f$
 * 
 *   @param[in] beta Coupling constant \f$beta\f$
 *   @param[in] n_plaq Number of plaquettes \f$P\f$
 */
double Phi_chit_perturbative(const double beta,
                             const unsigned int n_plaq);
                             
                             
/** @brief Compute functions required in calculation of topological susceptibility
 * function \f$\Phi_{\chi_t}\f$
 *
 * Evaluate functions \f$I_n(x)\f$, \f$I'_n(x)\f$ and \f$I''_n(x)\f$
 * required in calculation of \f$\chi_t\f$.
 *
 * @param[in] x Value at which to evaluate functions
 * @param[inout] In Vector which will contain \f$I_n(\beta)\f$
 * @param[inout] dIn Vector which will contain \f$I'_n(\beta)\f$
 * @param[inout] ddIn Vector which will contain \f$I''_n(\beta)\f$
 */
void compute_In(const double x,
                std::vector<double>& In,
                std::vector<double>& dIn,
                std::vector<double>& ddIn);

#endif // AUXILLIARY_HH
