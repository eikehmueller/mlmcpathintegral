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
  double besselI0;
  if (sigma_>100.) {
    besselI0 = 1./sqrt(M_PI*sigma_)*(1.+0.25/sigma_+9.0/(32.*sigma_*sigma_));
  } else {
    besselI0 = exp(-0.5*sigma_)*gsl_sf_bessel_I0(0.5*sigma_);
  }
  double Znorm_inv_ = 1./(2.0*M_PI*besselI0);
  return Znorm_inv_*exp(-sigma_*sin_x_half*sin_x_half);
}
