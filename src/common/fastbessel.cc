#include "fastbessel.hh"
/** @file fastbessel.cc
 * @brief Implementation of fastbessel.hh
 */

/* Evaluate scaled modified bessel function I_0 */
const double fast_bessel_I0_scaled(const double z) {
  if (z > 1100.) {
    double z_inv = 1. / z;
    double p = ModifiedBesselCoefficient<4>::value;
    p = z_inv * p + ModifiedBesselCoefficient<3>::value;
    p = z_inv * p + ModifiedBesselCoefficient<2>::value;
    p = z_inv * p + ModifiedBesselCoefficient<1>::value;
    p = z_inv * p + ModifiedBesselCoefficient<0>::value;
    return p / sqrt(2. * M_PI * z);
  } else if (z > 400.) {
    double z_inv = 1. / z;
    double p = ModifiedBesselCoefficient<5>::value;
    p = z_inv * p + ModifiedBesselCoefficient<4>::value;
    p = z_inv * p + ModifiedBesselCoefficient<3>::value;
    p = z_inv * p + ModifiedBesselCoefficient<2>::value;
    p = z_inv * p + ModifiedBesselCoefficient<1>::value;
    p = z_inv * p + ModifiedBesselCoefficient<0>::value;
    return p / sqrt(2. * M_PI * z);
  } else if (z > 200.) {
    double z_inv = 1. / z;
    double p = ModifiedBesselCoefficient<6>::value;
    p = z_inv * p + ModifiedBesselCoefficient<5>::value;
    p = z_inv * p + ModifiedBesselCoefficient<4>::value;
    p = z_inv * p + ModifiedBesselCoefficient<3>::value;
    p = z_inv * p + ModifiedBesselCoefficient<2>::value;
    p = z_inv * p + ModifiedBesselCoefficient<1>::value;
    p = z_inv * p + ModifiedBesselCoefficient<0>::value;
    return p / sqrt(2. * M_PI * z);
  } else if (z > 100.) {
    double z_inv = 1. / z;
    double p = ModifiedBesselCoefficient<7>::value;
    p = z_inv * p + ModifiedBesselCoefficient<6>::value;
    p = z_inv * p + ModifiedBesselCoefficient<5>::value;
    p = z_inv * p + ModifiedBesselCoefficient<4>::value;
    p = z_inv * p + ModifiedBesselCoefficient<3>::value;
    p = z_inv * p + ModifiedBesselCoefficient<2>::value;
    p = z_inv * p + ModifiedBesselCoefficient<1>::value;
    p = z_inv * p + ModifiedBesselCoefficient<0>::value;
    return p / sqrt(2. * M_PI * z);
  } else {
    return gsl_sf_bessel_I0_scaled(z);
  }
}
