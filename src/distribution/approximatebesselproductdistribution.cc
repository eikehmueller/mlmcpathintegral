#include "approximatebesselproductdistribution.hh"
/** @file approximatebesselproductdistribution.cc
 * @brief Implementation of approximatebesselproductdistribution.hh
 */

/* Evaluate at a given point */
double ApproximateBesselProductDistribution::evaluate(const double x,
                                                      const double x_p,
                                                      const double x_m) const {
    double x0 = x_p-x_m; // x_+ - x_- needs to be in [-2*pi,+2*pi]
    double z = x-x_m;
    // Map x_0 to range [0,pi]
    double sign_flip = (x0<0)?-1:+1; // Flip sign if difference is negative
    x0 *= sign_flip;
    if (x0 > M_PI) {
        x0 = 2.*M_PI - x0;
        sign_flip *= -1;
    }
    z *= sign_flip;
    double sigma2_inv = beta*cos(0.25*x0);
    double sigma2_inv_tilde = beta*sin(0.25*x0);
    double rho_r = gsl_sf_bessel_I0(2.*sigma2_inv_tilde)/gsl_sf_bessel_I0(2.*sigma2_inv);
    double N_p = 1./(1.+rho_r*rho_r);
    double N_m = 1.-N_p;
    double s_p=0.0;
    double s_m=0.0;
    for (int k=-kmax;k<=kmax;++k) {
        double z_shifted = z-0.5*x0+2*k*M_PI;
        s_p += exp(-0.5*sigma2_inv*z_shifted*z_shifted);
        z_shifted += M_PI;
        s_m += exp(-0.5*sigma2_inv*z_shifted*z_shifted);
    }
    return sqrt(0.5*sigma2_inv/M_PI)*(N_p*s_p + N_m*s_m);
}
