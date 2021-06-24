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
    double N_p;
    double sigma2_p_inv,sigma2_m_inv;
    compute_N_p_sigma2inv(beta,x0,N_p,sigma2_p_inv,sigma2_m_inv);
    double N_m = 1.-N_p;
    double s_p=0.0;
    double s_m=0.0;
    for (int k=-kmax;k<=kmax;++k) {
        double z_shifted = z-0.5*x0+2*k*M_PI;
        s_p += sqrt(sigma2_p_inv)*exp(-0.5*sigma2_p_inv*z_shifted*z_shifted);
        z_shifted += M_PI;
        s_m += sqrt(sigma2_m_inv)*exp(-0.5*sigma2_m_inv*z_shifted*z_shifted);
    }
    return sqrt(0.5/M_PI)*(N_p*s_p + N_m*s_m);
}

/* Compute probability N_p of drawing from the main mode of the distribution and
 * widths of Gaussians */
void ApproximateBesselProductDistribution::compute_N_p_sigma2inv(const double beta,
                                                                 const double x0,
                                                                 double& N_p,
                                                                 double& sigma2_p_inv,
                                                                 double& sigma2_m_inv) const {
    double epsilon = 0.125*M_PI;
    if (x0 < epsilon) {
        sigma2_p_inv = beta;
        sigma2_m_inv = 0.0;
        N_p = 1.0;
    } else {
        sigma2_p_inv = beta*cos(0.25*x0);
        sigma2_m_inv = beta*sin(0.25*x0);
        double rho = pow(sigma2_p_inv/sigma2_m_inv,1.5)*exp(-4.0*(sigma2_p_inv-sigma2_m_inv));
        N_p = 1.0/(1.0+rho);
    }
}

