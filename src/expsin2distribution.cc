#include "expsin2distribution.hh"
/** @file expsin2distribution.cc
 * @brief Implementation of expsin2distribution.hh
 */

/** Evaluate distribution */
double ExpSin2Distribution::evaluate(const double x) const {
  double sin_x_half = sin(0.5*x);
  return Znorm_inv*exp(-sigma*sin_x_half*sin_x_half);
}

/** Evaluate distribution */
double ExpSin2Distribution::evaluate(const double x, const double sigma_) const {
  double sin_x_half = sin(0.5*x);
  double Znorm_inv_ = 1./(2.0*M_PI*exp(-0.5*sigma_)*gsl_sf_bessel_I0(0.5*sigma_));
  return Znorm_inv_*exp(-sigma_*sin_x_half*sin_x_half);
}
