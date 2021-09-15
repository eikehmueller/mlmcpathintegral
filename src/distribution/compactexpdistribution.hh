#ifndef COMPACTEXPDISTRIBUTION_HH
#define COMPACTEXPDISTRIBUTION_HH COMPACTEXPDISTRIBUTION_HH
#include <random>
#include <cmath>

/** @file compactexpdistribution.hh
 * @brief Header file for exponential distribution on finite interval
 */


/** @class CompactExpDistribution 
 *
 * @brief Class for sampling from the distribution \f$p_\sigma(x) = Z_{\sigma}^{-1}\exp\left[\sigma x\right]\f$ with \f$x\in[-1,+1]\f$
 *
 * The normalisation constant is
 * \f$Z_{\sigma}=\frac{e^\sigma-e^{-\sigma}}{\sigma}\f$ and the cumulative distribution 
 * is \f$F_{\sigma}(x) = \frac{e^{\sigma x}-e^{-\sigma}}{e^{\sigma}-e^{-\sigma}}\f$.
 * 
 *
 * To sample from this distribution we use the inverse transform method
 * (see e.g. https://en.wikipedia.org/wiki/Inverse_transform_sampling). For this, we
 * draw a number \f$u\f$ uniformly from the interval \f$[0,1)\f$ and compute
 * 
 * \f[
 *    x = \frac{1}{\sigma}\log\left[u e^{\sigma}+(1-u)e^{-\sigma}\right]
 * \f]
 *
 * which is distributed according to \f$p_\sigma\f$.
 */

class CompactExpDistribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance
   *
   */
  CompactExpDistribution() : uniform_distribution(0.0,1.0) {}

  /** @brief Draw number from distribution for different \f$sigma\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] sigma Value of \f$\sigma\f$
   */
  template <class URNG>
  const double draw(URNG& engine, const double sigma) const {
    // Repeat until a number is accepted (and return in this case)
    double u = uniform_distribution(engine);
    double exp_sigma = exp(sigma);
    return 1./sigma*log(u*exp_sigma + (1.-u)/exp_sigma);
  }
  
  /** @brief Evaluate distribution for a different value of \f$\sigma\f$
  *
   * Calculate the value of the distribution \f$p_\sigma(x)\f$ at a point
   * \f$x\in[-1,+1]\f$.
   *
   * @param[in] x Point \f$x\f$ at which to evaluate the distribution
   * @param[in] sigma_ Value of sigma
   */
  double evaluate(const double x, const double sigma_) const;
  
private:
  
  /** @brief Uniform distribution for sampling */
  mutable std::uniform_real_distribution<double> uniform_distribution;
};

#endif // COMPACTEXPDISTRIBUTION_HH
