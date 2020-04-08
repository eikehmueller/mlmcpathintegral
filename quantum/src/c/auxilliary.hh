#ifndef AUXILLIARY_HH
#define AUXILLIARY_HH AUXILLIARY_HH
/** @file auxilliary.hh
 *
 * @brief Several auxilliary functions
 */
#include <cmath>

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

#endif // AUXILLIARY_HH
