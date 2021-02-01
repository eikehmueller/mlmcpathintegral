#ifndef EXPCOSDISTRIBUTION_HH
#define EXPCOSDISTRIBUTION_HH EXPCOSDISTRIBUTIONHH
#include <algorithm>
#include <random>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "common/auxilliary.hh"

/** @file expcosdistribution.hh
 *
 * @brief Header file for distribution given by the product of the exponentials of two cosines
 */

/** @class ExpCosDistribution
 *
 * @brief Class for sampling from the distribution \f$p(x|x_+,x_-) = Z_^{-1}\exp\left[\beta(\cos(x-x_+)+\cos(x-x_-))\right]\f$
 *
 * where \f$Z=2\pi I_0(2\beta\cos((x_+-x_-)/2))\f$ is the normalisation constant and \f$I_0\f$ the
 * modified Bessel function of first kind.
 *
 * To sample from this distribution we use rejection sampling with a Gaussian envelope.
 * Due to the symmetries of the distribution is is sufficient to be able to sample from \f$p(x|y,0)\f$.
 *
 * Note that this distribution can be reduced to a ExpSin2Distribution.
 */

class ExpCosDistribution {
public:
  /** @brief Constructor
   * 
   * Create a new instance for a given value of \f$\beta\f$
   *
   * @param[in] beta_ Parameter \f$\beta\f$
   */
  ExpCosDistribution(double beta_) :
    beta(beta_),
    distribution(0.0,1.0),
    normal_distribution(0.0,1.0),
    fourpi2_inv(1./(4.*M_PI*M_PI)) {}

  /** @brief Draw number from distribution for different \f$x_+,x_-\f$
   *
   * @param[in] engine Random number generator engine
   * @param[in] x_p Value of parameter \f$x_+\f$
   * @param[in] x_m Value of parameter \f$x_-\f$
   */
  template <class URNG>
    const double draw(URNG& engine, const double x_p, const double x_m) const {
        double dx = x_m-x_p;
        double tau = 2.*beta*fabs(cos(0.5*dx));
        double sigma = M_PI*sqrt(2./tau);
        bool accepted = false;
        double x;
        while (not accepted) {
            x = sigma*normal_distribution(engine);
            if ( (-M_PI <= x) and (x < M_PI) ) {
                double xi = distribution(engine);
                accepted = (xi <= exp(tau*(cos(x)-1.+fourpi2_inv*x*x)));
            }
        }
        return mod_2pi(x+0.5*(x_p+x_m)+(fabs(dx)>M_PI)*M_PI);
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
private:
    
  /** @brief parameter \f$\beta\f$*/
  const double beta;
  /** @brief Uniform distribution for rejection sampling */
  mutable std::uniform_real_distribution<double> distribution;
  /** @brief Normal distribution for approximate sampling */
  mutable std::normal_distribution<double> normal_distribution;
  /** Constant \f$1/(4\pi^2)\f$ used in rejection sampling*/
  const double fourpi2_inv;
};

#endif // EXPCOSDISTRIBUTION
