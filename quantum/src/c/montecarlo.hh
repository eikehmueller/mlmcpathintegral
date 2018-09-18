#ifndef MONTECARLO_HH
#define MONTECARLO_HH MONTECARLO_HH
#include <utility>
#include <cmath>
#include <iostream>
#include "sampler.hh"
#include "action.hh"
#include "quantityofinterest.hh"
#include "statistics.hh"

/** @file montecarlo.hh
 * @brief Header file for Monte Carlo base class
 */

/** @class MonteCarlo
 * 
 * @brief Monte Carlo base class
 */

class MonteCarlo {
public:
  /** @brief Create new instance
   *
   * @param[in] n_samples_ Number of samples
   * @param[in] n_burnin_ Number of burn-in steps
   */
  MonteCarlo(const unsigned int n_samples_,
             const unsigned int n_burnin_) :
    n_samples(n_samples_), 
    n_burnin(n_burnin_) {}
protected:
  /** @brief Number of samples */
  const unsigned int n_samples;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
};

#endif // MONTECARLO_HH
