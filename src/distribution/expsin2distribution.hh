#ifndef EXPSIN2DISTRIBUTION_HH
#define EXPSIN2DISTRIBUTION_HH EXPSIN2DISTRIBUTIONHH
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "common/auxilliary.hh"

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
 * To sample from this distribution we use rejection sampling with a Gaussian envelope.
 */

class ExpSin2Distribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance
   *
   */
  ExpSin2Distribution() : distribution(0.0,1.0),
                          normal_distribution(0.0,1.0) {}

  /** @brief Draw number from distribution for different \f$sigma\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] sigma Value of \f$\sigma\f$
   */
  template <class URNG>
  const double draw(URNG& engine, const double sigma) const {
    // Repeat until a number is accepted (and return in this case)
    while (true) {
      double r_x = M_PI/sqrt(2.*sigma)*normal_distribution(engine);
      if (fabs(r_x)<M_PI) {
        double sin_psi_half = sin(0.5*r_x);
        double r_accept = distribution(engine);
        if ( r_accept < exp(-sigma*(sin_psi_half*sin_psi_half-r_x*r_x/(M_PI*M_PI)))) {
          return r_x;
        }
      }
    }
  }
  
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
  /** @brief Fast bessel function
   * Compute \f$2\pi e^{-z}I_0(z)\f$, using a Taylor approximation for large values of z
   *
   *  @param[in] z Value at which the function is evaluated
   */
  double fast_2pi_I0_scaled(const double z) const;
  
  /** @brief Uniform distribution for sampling */
  mutable std::uniform_real_distribution<double> distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
};

#endif // EXPSIN2DISTRIBUTIONHH
