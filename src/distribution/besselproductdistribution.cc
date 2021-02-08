#include "besselproductdistribution.hh"
/** @file besselproductdistribution.cc
 * @brief Implementation of besselproductdistribution.hh
 */

/* Evaluate at a given point */
double BesselProductDistribution::evaluate(const double x,
                                           const double x_p,
                                           const double x_m) const {
    double I0_p = gsl_sf_bessel_I0(2*beta*cos(0.5*(x-x_p)));
    double I0_m = gsl_sf_bessel_I0(2*beta*cos(0.5*(x-x_m)));
    return Znorm_inv(x_p-x_m)*I0_p*I0_m;
}

/* Compute inverse normalisation constant */
const double BesselProductDistribution::Znorm_inv(double phi,
                                                  const bool rescaled) const {
    double s=1.0;
    for (unsigned int k=1;k<=kmax;++k) {
        s += alphaZ[k]*cos(k*phi);
    }
    if (not rescaled) {
        s*=alphaZ[0];
    }
    return 1.0/s;
}
