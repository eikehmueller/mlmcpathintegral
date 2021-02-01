#include "expcosdistribution.hh"
/** @file expcosdistribution.cc
 * @brief Implementation of expcosdistribution.hh
 */

/* Evaluate at a given point */
double ExpCosDistribution::evaluate(const double x,
                                    const double x_p,
                                    const double x_m) const {
    double Z_norm = 2.*M_PI*gsl_sf_bessel_I0(2.*beta*cos(0.5*(x_p-x_m)));
    return 1./Z_norm*exp(beta*(cos(x-x_p)+cos(x-x_m)));
}
