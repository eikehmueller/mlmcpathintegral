#ifndef AUXILLIARY_HH
#define AUXILLIARY_HH AUXILLIARY_HH
/** @file auxilliary.hh
 *
 * @brief Several auxilliary functions
 */
#include "config.h"
#include <cmath>

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

/** @brief Calculate Sum \f$\hat{\Sigma_p(\xi)}\f$
 *
 * Return value of sum
 *
 \f[
 \hat{\Sigma}_p(\xi) = \frac{\sum_{m\in\mathbb{Z}}m^p \exp\left[-\frac{1}{2}\xi m^2}\right]}{\sum_{m\in\mathbb{Z}}\exp\left[-\frac{1}{2}\xi m^2}\right]}
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

#endif // AUXILLIARY_HH
