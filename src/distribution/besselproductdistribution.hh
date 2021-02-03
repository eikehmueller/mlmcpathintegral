#ifndef BESSELPRODUCTDISTRIBUTION_HH
#define BESSELPRODUCTDISTRIBUTION_HH BESSELPRODUCTDISTRIBUTIONHH
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include "common/auxilliary.hh"
#include "mpi/mpi_wrapper.hh"

/** @file besselproductdistribution.hh
 * @brief Header file for distribution given by the product of two modified Bessel functions
 */


/** @class BesselProductDistribution
 *
 * @brief Class for sampling from the distribution \f$p(x|x_+,x_-) = Z_^{-1}I_0(2\beta \cos((x-x_+)/2))I_0(2\beta \cos((x-x_-)/2))\f$
 *
 * where \f$I_0\f$ is the modified Bessel function of the first kind and \f$Z\f$
 * is a normalisation constant.
 *
 * To sample from this distribution we use rejection sampling with a Gaussian envelope.
 * Due to the symmetries of the distribution is is sufficient to be able to sample from \f$p(x|y,0)\f$ with
 * \f$y\ge 0\f$. In this case it is possible to construct two separate Gaussian envelopes for the two
 * intervals \f$[-\pi,-\pi+y]\f$ and \f$[-\pi+y,+\pi]\f$. To generate a sample, first work out which of
 * those two intervals the envelope sample is draw from (using the appropriate integrals over the envelope
 * distribution), then draw a sample from the envelope distribution and finally accept/reject with the
 * appropriate ratio.
 *
 * Note that the normalisation constant of the distribution \f$p(x|x_+,x_-)\f$ can be written as a Fourier-
 * cosine expansion of the angle \f$x_+-x_-\f$, where the expansion coefficients depend on the value of the
 * parameter \f$\beta\f$. Some care needs to be taken when truncating this sum, the code currently only reliably
 * supports value of \f$\beta\f$ not exceeding 8. For larger values, kmax and nmax would need to be adjusted.
 */

class BesselProductDistribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance for a given value of \f$\beta\f$
   *
   * @param[in] beta_ Parameter \f$\beta\f$
   */
  BesselProductDistribution(double beta_) :
    beta(beta_),
    distribution(0.0,1.0),
    normal_distribution(0.0,1.0),
    kmax(20), nmax(40),
    I0_twobeta(gsl_sf_bessel_I0(2*beta)),
    sigma_beta(M_PI/sqrt(2*log(I0_twobeta))) {
        if (beta > 16.0) {
            mpi_parallel::cerr << "ERROR beta must not exceed 16.0" << std::endl;
            exit(-1);
        }
        /* Compute expansion coefficients for normalisation constant */
        for (unsigned int k=0;k<=kmax;++k) {
            double s = 0.0;
            for (unsigned int n=k;n<=nmax;++n) {
                for (unsigned int m=k;m<=nmax;++m) {
                    double log_comb = log_nCk(2*n,n-k)+log_nCk(2*m,m-k)-2*(log_factorial(n)+log_factorial(m));
                    s += pow(0.5*beta,2*(n+m))*exp(log_comb);
                }
            }
            alphaZ.push_back(((k==0)?2:4)*M_PI*s);
        }
    }

  /** @brief Draw number from distribution for different \f$x_+,x_-\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] x_p Value of parameter \f$x_+\f$
   * @param[in] x_m Value of parameter \f$x_-\f$
   */
  template <class URNG>
    const double draw(URNG& engine, const double x_p, const double x_m) const {
        double a_min, a_max, mu;
        double dx = x_m-x_p;
        // Flip sign if necessary to ensure dx is always non-negative
        // (this will be compensated for later)
        double sign_flip = (dx<0)?-1:+1; // Flip sign if difference is negative
        dx *= sign_flip;
        // Work out (scaled) normalisation constants of both sides of the
        // envelope density
        double N_p = gsl_sf_erf((M_PI-0.5*dx)/sigma_beta);
        double N_m = gsl_sf_erf(0.5*dx/sigma_beta)*pow(I0_twobeta,2.*(dx/M_PI-1.));
        double C_gauss_p = pow(I0_twobeta,2.*(1.-dx*dx/(4.*M_PI*M_PI)));
        double C_gauss_m = pow(I0_twobeta,2.*(1.-(dx-2.*M_PI)*(dx-2.*M_PI)/(4.*M_PI*M_PI)));
        double sigma = sigma_beta/sqrt(2.);
        bool accepted = false;
        double xi, x, C_gauss;
        while (not accepted) {
            xi = distribution(engine);
            // Work out which part of the envelope distribution will be used
            if (xi >= N_m/(N_p+N_m)) {
                // right side
                a_min = -M_PI+dx;
                a_max = +M_PI;
                mu = 0.5*dx;
                C_gauss = C_gauss_p;
            } else {
                // left side
                a_min = -M_PI;
                a_max = -M_PI+dx;
                mu = 0.5*(dx-2.*M_PI);
                C_gauss = C_gauss_m;
            }
            // Draw from envelope distribution
            // (keep trying until we find a sample in the
            // desired interval)
            while (not accepted) {
                x = sigma*normal_distribution(engine)+mu;
                accepted = ( (x >= a_min) and (x < a_max) );
            }
            // Check by comparing envelope to actual distribution
            double I0 = gsl_sf_bessel_I0(2.*beta*cos(0.5*x));
            double I0_dx = gsl_sf_bessel_I0(2.*beta*cos(0.5*(x-dx)));
            double x_shifted = (x-mu)/sigma_beta;
            double rho_accept = I0*I0_dx/C_gauss*exp(x_shifted*x_shifted);
            xi = distribution(engine);
            accepted = (xi<=rho_accept);
        }
        // flip sign if necessary and add x_p back on
        return mod_2pi(sign_flip*x+x_p);
    }
  
  /** @brief Evaluate distribution for a given value of \f$x\f$ and parameters \f$x_+,x_-\f$
   *
   * Calculate the value of the distribution \f$p(x|x_+,x_-)\f$ at a point
   * \f$x\in[-\pi,\pi]\f$.
   *
   * @param[in] x Point \f$x\f$ at which to evaluate the distribution
   * @param[in] x_p Value of parameter \f$x_+\f$
   * @param[in] x_m Value of parameter \f$x_-\f$
   */
  double evaluate(const double x, const double x_p, const double x_m) const;
  
    /** @brief Return parameter \f$\beta\f$ */
  double get_beta() const { return beta; }

  /** @brief Inverse of normalisation constant
   *
   * Returns \f$Z^{-1}\f$
   *
   * @param[in] phi angle \f$\phi\f$
   */
   const double Znorm_inv(double phi) const;

private:
    
  /** @brief parameter \f$\beta\f$*/
  const double beta;
  /** @brief Uniform distribution for sampling */
  mutable std::uniform_real_distribution<double> distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** @brief Number  of expansion terms for normalisation constant */
  const unsigned int kmax;
  /** @brief Number  of expansion terms in computation of normalisation constant */
  const unsigned int nmax;
  /** @brief Expansion coefficients for normalisation constant */
  std::vector<double> alphaZ;
  /** Precomputed value of \f$I_0(2\beta)\f$ */
  const double I0_twobeta;
  /** Width of Gaussian envelope \f$\sigma_\beta=\frac{\pi}{\sqrt{2\log(I_0(2\beta))}}\f$ */
  const double sigma_beta;
};

#endif // BESSELPRODUCTDISTRIBUTION
