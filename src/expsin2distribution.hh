#ifndef EXPSIN2DISTRIBUTION_HH
#define EXPSIN2DISTRIBUTION_HH EXPSIN2DISTRIBUTIONHH
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "auxilliary.hh"

/** @file expsin2distribution.hh
 * @brief Header file for exponential sin-squared distribution 
 */


/** @class ExpSin2Distribution 
 *
 * @brief Class for sampling from the distribution \f$p(x) = Z_{\sigma}^{-1}\exp\left[-\sigma \sin^2\left(\frac{x}{2}\right)\right]\f$ with \f$x\in[-\pi,\pi]\f$
 *
 * The normalisation constant is
 * \f$Z_{\sigma}=2\pi e^{-\frac{\sigma}{2}}I_0\left(\frac{\sigma}{2}\right)\f$
 * where \f$I_0\f$ is the modified Bessel function of the first kind.
 * 
 * To sample from this distribution we use rejection sampling.
 * For this, consider the distribution on \f$[0,\pi]\f$ defined by
 * \f$\hat{p}(x) = 2p(x)\f$. Construct a piecewise constant distribution
 * \f$q(x)=\left(Z'_{\sigma}\right)^{-1} \exp\left[-\sigma f(x)\right]\f$ where the function
 * \f$f(x)\le \sin^2\left(\frac{x}{2}\right)\f$ is piecewise constant on
 * the intervals \f$[x_i,x_{i+1}]\f$ with \f$x_0=0\f$, \f$x_n=\pi\f$. Draw
 * a point \f$X\f$ from \f$q(x)\f$ using the cumulative distribution function,
 * and then accept with probability \f$r(x) = C\hat{p}(X)/p(X)=\exp\left[-\sigma\left(\sin^2\left(\frac{X}{2}\right)-f(x)\right)\right]\f$. Note that the 
 * constant \f$C\f$ is chosen such that the ratio \f$r\f$ does not exceed
 * \f$1\f$.
 */

class ExpSin2Distribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance
   *
   * @param[in] sigma_ Parameter \f$\sigma\f$ of distribution
   * @param[in] nint_ Number of intervals \f$n\f$ used for piecewise
   *                  constant distribution
   */
  ExpSin2Distribution(const double sigma_,
                      const unsigned int nint_=0) : sigma(sigma_),
                                                    nint(nint_),
                                                    Znorm_inv(1./(2.0*M_PI*exp(-0.5*sigma)*gsl_sf_bessel_I0(0.5*sigma))),
                                                    distribution(0.0,1.0),
                                                    normal_distribution(0.0,1.0),
                                                    sigma_threshold(16.) {
    /* Calculate cumulative probability density for piecewise constant
     * distribution */
    // Width of intervales
   if (nint>0) {
      double h = M_PI/(1.0*nint);
      p_cdf.push_back(0.0);
      for (unsigned int i=0;i<=nint;++i) {
        double x_tmp = i*h;
        x.push_back(x_tmp);
        if (i<nint) {
          double sin_psi_half = sin(0.5*x_tmp);
          double p = exp(-sigma*sin_psi_half*sin_psi_half);
          p_pdf.push_back(p);
          p_cdf.push_back(p_cdf.back()+p);
        }
      }
      double p_cdf_total_inv = 1.0/p_cdf[nint];
      for (unsigned int i=0;i<=nint;++i) {
        p_cdf[i] *= p_cdf_total_inv;
      }
    }
  }

  /** @brief Draw number from distribution
   *
   * @param[in] engine Random number generator engine
   */
  template <class URNG>
  const double draw(URNG& engine) const {
    // Repeat until a number is accepted (and return in this case)
    if (nint==0) {
      return draw(engine,sigma);
    } else {
      while (true) {
        double r_interval = distribution(engine);
        // Find interval
        unsigned int i=0;
        while (r_interval > p_cdf[i+1]) i++;
        // Draw uniform number from interval
        double r_x = x[i] + distribution(engine)*(x[i+1]-x[i]);
        // accept or reject with given probability
        double r_accept = distribution(engine);
        double sin_psi_half = sin(0.5*r_x);
        if ( r_accept < exp(-sigma*sin_psi_half*sin_psi_half)/p_pdf[i] ) {
          // Find sign
          double r_sign = 2.0*(distribution(engine)>0.5)-1.0;
          return r_sign*r_x;
        }
      }
    }
  }

  /** @brief Draw number from distribution for different \f$sigma\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] sigma_ Value of \f$\sigma\f$
   */
  template <class URNG>
  const double draw(URNG& engine, const double sigma_) const {
    if (sigma_>sigma_threshold) {
      // Sample from approximate distribution
      return sqrt(2./sigma_)*normal_distribution(engine);
    } else {
      // Repeat until a number is accepted (and return in this case)
      while (true) {
        double r_x = M_PI*(2.0*distribution(engine)-1.0);
        double sin_psi_half = sin(0.5*r_x);
        double r_accept = distribution(engine);
        if ( r_accept < exp(-sigma_*sin_psi_half*sin_psi_half)) {
          return r_x;
        }
      }
    }
  }

  
  /** @brief Evaluate distribution 
  *
   * Calculate the value of the distribution \f$p(x)\f$ at a point
   * \f$x\in[-\pi,\pi]\f$.
   *
   * @param[in] x Point \f$x\f$ at which to evaluate the distribution
   */
  double evaluate(const double x) const;

  /** @brief Evaluate distribution for a different value of \f$\sigma\f$
  *
   * Calculate the value of the distribution \f$p(x)\f$ at a point
   * \f$x\in[-\pi,\pi]\f$.
   *
   * @param[in] x Point \f$x\f$ at which to evaluate the distribution
   * @param[in] sigma_ Value of sigma
   */
  double evaluate(const double x, const double sigma_) const;
  
private:
  /** @brief Parameter \f$\sigma\f$ of distribution */
  const double sigma;
  /** @brief Number of intervals of piecewise constant distribution */
  const unsigned int nint;
  /** @brief Vector with points \f$x_i\f$ defininit intervals for piecewise constant distribution */
  std::vector<double> x;
  /** @brief Probability in \f$i\f$-th interval of piecewise constant distribution */
  std::vector<double> p_pdf;
  /** @brief Cumulative probability in \f$i\f$-th interval of piecewise constant distribution */
  std::vector<double> p_cdf;
  /** @brief Inverse normalisation constant \f$Z_{\sigma}^{-1}\f$ */
  const double Znorm_inv;
  /** @brief Uniform distribution for sampling */
  mutable std::uniform_real_distribution<double> distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** @brief Threshold for using approximate distribution */
  const double sigma_threshold;
};

#endif // EXPSIN2DISTRIBUTIONHH
