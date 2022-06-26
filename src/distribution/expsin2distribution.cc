#include "expsin2distribution.hh"
/** @file expsin2distribution.cc
 * @brief Implementation of expsin2distribution.hh
 */

/* Fast bessel function */
double ExpSin2Distribution::fast_2pi_I0_scaled(const double z) const {
  if (z > 100.) {
    /* Use Taylor expansion for large arguments, see https://dlmf.nist.gov/10.40
     */
    double z_inv = 1. / z;
    return sqrt(2. * M_PI * z_inv) *
           (1. + 0.125 * z_inv + 0.0703125 * z_inv * z_inv);
  } else {
    return 2. * M_PI * gsl_sf_bessel_I0_scaled(z);
  }
}

/** Evaluate distribution */
double ExpSin2Distribution::evaluate(const double x, const double sigma) const {
  double sin_x_half = sin(0.5 * x);
  double besselI0 = fast_2pi_I0_scaled(0.5 * sigma);
  return exp(-sigma * sin_x_half * sin_x_half) / besselI0;
}
