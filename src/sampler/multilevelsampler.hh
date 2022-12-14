#ifndef MULTILEVELSAMPLER_HH
#define MULTILEVELSAMPLER_HH MULTILEVELSAMPLER_HH
#include "action/action.hh"
#include "action/conditionedfineaction.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "common/statistics.hh"
#include "common/timer.hh"
#include "montecarlo/twolevelmetropolisstep.hh"
#include "mpi/mpi_wrapper.hh"
#include "qoi/quantityofinterest.hh"
#include "sampler/hierarchicalsampler.hh"
#include "sampler/sampler.hh"
#include <memory>
#include <sstream>
#include <string>
#include <vector>

/** @file multilevelsampler.hh
 * @brief Header file for Hierarchical sampler class
 */

/** @class MultilevelSampler
 * @brief Multilevel sampler
 *
 * Implementation of multilevel sampling. Generate samples by successively
 * screening samplers from coarser levels
 */
class MultilevelSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] fine_action_ Fine level Action to sample from
   * @param[in] qoi_factory_ Quantiy of interest factory
   * @param[in] coarse_sampler_factory Factory for sampler on coarsest level
   * @param[in] conditioned_fine_action_factory Factory for conditioned fine
   * actions
   * @param[in] param_general General parameters
   * @param[in] param_statistics Statistics parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  MultilevelSampler(
      const std::shared_ptr<Action> fine_action,
      const std::shared_ptr<QoIFactory> qoi_factory_,
      const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
      const std::shared_ptr<ConditionedFineActionFactory>
          conditioned_fine_action_factory,
      const StatisticsParameters param_stats,
      const HierarchicalParameters param_hierarchical);
  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~MultilevelSampler() {}

  /** @brief Draw a sample
   *
   * returns a sample state \f$\phi\f$
   *
   * @param[out] phi_state State \f$\phi\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<SampleState> phi_state);

  /** @brief Show statistics
   *
   * Print out statistics
   */
  virtual void show_stats();

  /** @brief Set current state to particular value
   *
   * @param[in] phi_state
   */
  virtual void set_state(std::shared_ptr<SampleState> phi_state);

  /** Return cost per sample */
  virtual double cost_per_sample() { return cost_per_sample_; }

private:
  /** @brief Action to sample from */
  std::vector<std::shared_ptr<Action>> action;
  /** @brief Quantity of interest */
  std::vector<std::shared_ptr<QoI>> qoi;
  /** @brief Number of levels in hierarchy */
  const unsigned int n_level;
  /** @brief Two level step on all levels of the multigrid hierarchy */
  std::vector<std::shared_ptr<TwoLevelMetropolisStep>> twolevel_step;
  /** @brief Sampler on coarsest level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Sampler state on each level  */
  std::vector<std::shared_ptr<SampleState>> phi_sampler_state;
  /** @brief number of skipped samples between independent samples on all levels
   */
  std::vector<double> t_indep;
  /** @brief number of independent samples on all levels */
  std::vector<int> n_indep;
  /** @brief number of samples generated since last independent sample */
  std::vector<int> t_sampler;
  /** @brief vector with statistics of sampler states */
  std::vector<std::shared_ptr<Statistics>> stats_sampler;
  /** @brief size of window for computing autocorrleations */
  unsigned int n_autocorr_window;
  /** @brief cost per sample */
  mutable double cost_per_sample_;
};

class MultilevelSamplerFactory : public SamplerFactory {
public:
  /** @brief Create new instance
   *
   * @param[in] qoi_factory_ Factory for quantity of interest
   * @param[in] coarse_sampler_factory_ Factory for coarse level sampler
   * @param[in] conditioned_fine_action_factory_ Factory for conditioned fine
   * action
   * @param[in] param_stats Statistics parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  MultilevelSamplerFactory(
      const std::shared_ptr<QoIFactory> qoi_factory_,
      const std::shared_ptr<SamplerFactory> coarse_sampler_factory_,
      const std::shared_ptr<ConditionedFineActionFactory>
          conditioned_fine_action_factory_,
      const StatisticsParameters param_stats_,
      const HierarchicalParameters param_hierarchical_)
      : qoi_factory(qoi_factory_),
        coarse_sampler_factory(coarse_sampler_factory_),
        conditioned_fine_action_factory(conditioned_fine_action_factory_),
        param_stats(param_stats_), param_hierarchical(param_hierarchical_) {}

  /** @brief Destructor */
  virtual ~MultilevelSamplerFactory() {}

  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::make_shared<MultilevelSampler>(
        action, qoi_factory, coarse_sampler_factory,
        conditioned_fine_action_factory, param_stats, param_hierarchical);
  }

private:
  /** Quantity of interest */
  const std::shared_ptr<QoIFactory> qoi_factory;
  /** Factory for coarsest level sampler*/
  const std::shared_ptr<SamplerFactory> coarse_sampler_factory;
  /** Factory for conditioned fine action */
  const std::shared_ptr<ConditionedFineActionFactory>
      conditioned_fine_action_factory;
  /** Statistics parameters */
  const StatisticsParameters param_stats;
  /** Hierarchical sampler parameters */
  const HierarchicalParameters param_hierarchical;
};

#endif // MULTILEVELSAMPLER
