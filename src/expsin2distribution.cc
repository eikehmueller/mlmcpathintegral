#include "expsin2distribution.hh"
/** @file expsin2distribution.cc
 * @brief Implementation of expsin2distribution.hh
 */

/** Evaluate distribution */
double ExpSin2Distribution::evaluate(const double x, const double sigma) const {
  double sin_x_half = sin(0.5*x);
  double besselI0;
  if (sigma>100.) {
    double sigma_inv = 1./sigma;
    besselI0 = 2.*sqrt(M_PI*sigma_inv)*(1.+0.25*sigma_inv+1.125*sigma_inv*sigma_inv);
  } else {
    besselI0 = 2*M_PI*gsl_sf_bessel_I0_scaled(0.5*sigma);
  }
  return exp(-sigma*sin_x_half*sin_x_half)/besselI0;
}
