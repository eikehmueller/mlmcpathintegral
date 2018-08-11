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

#endif // AUXILLIARY_HH
