#include "expcosdistribution.hh"
/** @file expcosdistribution.cc
 * @brief Implementation of expcosdistribution.hh
 */

/* Evaluate at a given point */
double ExpCosDistribution::evaluate(const double x, const double x_p,
                                    const double x_m) const {
  double dx = x_p - x_m;
  double z = x - x_m;
  int sign_flip = (dx < 0.0) ? -1 : +1;
  dx *= sign_flip;
  if (dx > M_PI) {
    sign_flip *= -1;
    dx = 2. * M_PI - dx;
  }
  z *= sign_flip;
  double sigma = 2. * beta * fabs(cos(0.5 * dx));
  double Z_norm = 2. * M_PI * fast_bessel_I0_scaled(sigma);
  return 1. / Z_norm * exp(sigma * (cos(z - 0.5 * dx) - 1.0));
}
