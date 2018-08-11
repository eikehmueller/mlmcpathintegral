#include "expsin2distribution.hh"
/** @file expsin2distribution.cc
 * @brief Implementation of expsin2distribution.hh
 */

/** Evaluate distribution */
double ExpSin2Distribution::evaluate(const double x) const {
  double sin_x_half = sin(0.5*x);
  return Znorm_inv*exp(-sigma*sin_x_half*sin_x_half);
}
