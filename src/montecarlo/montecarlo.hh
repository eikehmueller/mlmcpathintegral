#ifndef MONTECARLO_HH
#define MONTECARLO_HH MONTECARLO_HH
#include "action/action.hh"
#include "common/statistics.hh"
#include "qoi/quantityofinterest.hh"
#include "sampler/sampler.hh"
#include <cmath>
#include <iostream>
#include <utility>

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
   * @param[in] n_burnin_ Number of burn-in steps
   */
  MonteCarlo(const unsigned int n_burnin_) : n_burnin(n_burnin_) {}

protected:
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
};

#endif // MONTECARLO_HH
