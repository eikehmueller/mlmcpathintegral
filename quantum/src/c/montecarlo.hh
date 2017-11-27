#ifndef MONTECARLO_HH
#define MONTECARLO_HH MONTECARLO_HH
#include <utility>
#include <cmath>
#include "sampler.hh"
#include "action.hh"
#include "quantityofinterest.hh"
#include "twolevelmetropolissampler.hh"

/** @class MonteCarlo
 * 
 * @brief Monte Carlo base class
 */

class MonteCarlo {
public:
  MonteCarlo(unsigned int n_samples_,
             unsigned int n_burnin_) : 
    n_samples(n_samples_), 
    n_burnin(n_burnin_) {}
protected:
  const unsigned int n_samples;
  const unsigned int n_burnin;
};

/** @class MonteCarloSingleLevel
 *
 * @brief Single level Monte Carlo sampler
 *
 * Calculates MCMC estimator by drawing samples and evaluating the 
 * quantity of interest on those
 */
class MonteCarloSingleLevel : public MonteCarlo {
public:
  /** @brief Create new instance
   *
   * @param[in] sampler Sampler to draw from
   * @param[in] qoi Quantity of interest to evaluate on samples
   * @param[in] n_samples Number of samples to evaluate
   * @param[in] n_burnin Number of burnin samples (QoI not evaluated on those)
   */
  MonteCarloSingleLevel(Sampler& sampler_,
                        QoI& qoi_,
                        unsigned int n_samples_,
                        unsigned int n_burnin_) :
    MonteCarlo(n_samples_,n_burnin_), sampler(sampler_), qoi(qoi_)
  {}

  /** @brief Calculate QoI
   * 
   * Calculate the Quantity of interest by Monte Carlo sampling. Return
   * estimator for the mean and standard deviation.
   */
  std::pair<double,double> evaluate();
private:
  Sampler& sampler;
  QoI& qoi;
};

/** @class MonteCarloTwoLevel
 * 
 * @brief Two level Monte Carlo method
 * 
 * 
 */
class MonteCarloTwoLevel : public MonteCarlo {
public:
  /** \brief Create new instance 
   *
   * @param[in] coarse_sample_ Sampler on coarse level
   * @param[in] coarse_action_ Action on coarse level
   * @param[in] fine_action_ Action on fine level
   * @param[in] coarse_qoi_ Quantity of interest on coarse level
   * @param[in] fine_qoi_ Quantity of interest on fine level
   * @param[in] n_samples Number of samples to evaluate
   * @param[in] n_burnin Number of burnin samples (QoI not evaluated on those)
   */
  MonteCarloTwoLevel(Sampler& coarse_sampler_,
                     Action& coarse_action_,
                     Action& fine_action_,
                     QoI& qoi_,
                     unsigned int n_samples_,
                     unsigned int n_burnin_) : 
    MonteCarlo(n_samples_,n_burnin_),
    coarse_sampler(coarse_sampler_),
    coarse_action(coarse_action_),
    fine_action(fine_action_),
    qoi(qoi_),
    twolevel_sampler(coarse_sampler_,
                     coarse_action_,
                     fine_action_) {}

  /** @brief Calculate mean and variance of difference
   *
   * Calculate the mean and variance of the difference in the QoI 
   * evaluated at two subsequent levels.
   */
  std::pair<double,double> evaluate_difference();

private:
  Sampler& coarse_sampler;
  Action& coarse_action;
  Action& fine_action;
  QoI& qoi;
  TwoLevelMetropolisSampler twolevel_sampler;
};

#endif // MONTECARLO_HH
