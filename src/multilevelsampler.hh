#ifndef MULTILEVELSAMPLER_HH
#define MULTILEVELSAMPLER_HH MULTILEVELSAMPLER_HH
#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include "mpi_wrapper.hh"
#include "path.hh"
#include "action.hh"
#include "parameters.hh"
#include "sampler.hh"
#include "conditionedfineaction.hh"
#include "twolevelmetropolisstep.hh"
#include "statistics.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"
#include "statistics.hh"
#include "quantityofinterest.hh"
#include "timer.hh"
#include "hierarchicalsampler.hh"

/** @file multilevelsampler.hh
* @brief Header file for Hierarchical sampler class
*/

/** @class MultilevelSampler
 * @brief Multilevel sampler
 *
 * Implementation of multilevel sampling. Generate samples by successively screening samplers
 * from coarser levels
 */
class MultilevelSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] action_ Fine level Action to sample from
   * @param[in] param_general General parameters
   * @param[in] param_hmc HMC sampler parameters
   * @param[in] param_cluster Cluster sampler parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  MultilevelSampler(const std::shared_ptr<Action> fine_action,
                    const std::shared_ptr<QoI> qoi_,
                    const GeneralParameters param_general,
                    const StatisticsParameters param_stats,
                    const HMCParameters param_hmc,
                    const ClusterParameters param_cluster,
                    const HierarchicalParameters param_hierarchical);
  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~MultilevelSampler() {}

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<Path> x_path);
  
/** @brief Show statistics
     *
     * Print out statistics
     */
    virtual void show_stats();
  
  /** @brief Return cost per sample (in microseconds) */
  double cost_per_sample() {
    return cost_per_sample_;
  }
  
  /** @brief Set current state to particular value
   *
   * @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<Path> x_path);
  
private:
  /** @brief Action to sample from */
  std::vector<std::shared_ptr<Action>> action;
  /** @brief Quantity of interest */
  const std::shared_ptr<QoI> qoi;
  /** @brief Number of levels in hierarchy */
  const unsigned int n_level;
  /** @brief Two level step on all levels of the multigrid hierarchy */
  std::vector<std::shared_ptr<TwoLevelMetropolisStep> > twolevel_step;
  /** @brief Sampler on coarsest level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Sampler path on each level  */
  std::vector<std::shared_ptr<Path> > x_sampler_path;
  /** @brief cost per sample */
  double cost_per_sample_;
  /** @brief number of skipped samples between independent samples on all levels */
  std::vector<double> t_indep;
  /** @brief number of independent samples on all levels */
  std::vector<int> n_indep;
  /** @brief number of samples generated since last independent sample */
  std::vector<int> t_sampler;
  /** @brief vector with statistics of sampler paths */
  std::vector<std::shared_ptr<Statistics> > stats_sampler;
  /** @brief size of window for computing autocorrleations */
  unsigned int n_autocorr_window;
};

#endif // MULTILEVELSAMPLER
