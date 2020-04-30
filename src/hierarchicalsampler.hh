#ifndef HIERARCHICALSAMPLER_HH
#define HIERARCHICALSAMPLER_HH HIERARCHICALSAMPLER_HH
#include <vector>
#include <memory>
#include "mpi_wrapper.hh"
#include "path.hh"
#include "action.hh"
#include "parameters.hh"
#include "sampler.hh"
#include "montecarlomultilevel.hh"
#include "conditionedfineaction.hh"
#include "twolevelmetropolisstep.hh"
#include "statistics.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"

/** @file hierarchicalsampler.hh
 * @brief Header file for Hierarchical sampler class
 */

/** @class HierarchicalSampler
 * @brief Hierarchical sampler
 *
 * Implementation of hierarchical sampling. Generate samples by successively screening samplers
 * from coarser levels
 */
class HierarchicalSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] action_ Fine level Action to sample from
   * @param[in] param_general General parameters
   * @param[in] param_hmc HMC sampler parameters
   * @param[in] param_cluster Cluster sampler parameters
   * @param[in] param_multilevelmc Multilevel parameters
   */
  HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                      const GeneralParameters param_general,
                      const HMCParameters param_hmc,
                      const ClusterParameters param_cluster,
                      const MultiLevelMCParameters param_multilevelmc);
  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~HierarchicalSampler() {}

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<Path> x_path);
  
private:
  /** @brief Action to sample from */
  std::vector<std::shared_ptr<Action>> action;
  /** @brief Number of levels in hierarchy */
  const unsigned int n_level;
  /** @brief Two level step on all levels of the multigrid hierarchy */
  std::vector<std::shared_ptr<TwoLevelMetropolisStep> > twolevel_step;
  /** @brief Sampler on coarsest level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Sampler path on each level  */
  std::vector<std::shared_ptr<Path> > x_sampler_path;

};

#endif // HIERARCHICALSAMPLER
