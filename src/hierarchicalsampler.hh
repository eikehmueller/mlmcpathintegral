#ifndef HIERARCHICALSAMPLER_HH
#define HIERARCHICALSAMPLER_HH HIERARCHICALSAMPLER_HH
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
#include "timer.hh"

/** @class HierarchicalParameters
 *
 * @brief Class for storing parameters of Hierarchical sampler.
 */
class HierarchicalParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  HierarchicalParameters() :
    Parameters("hierarchical"),
    n_level_(2),
  coarsesampler_(SamplerHMC) {
    addKey("n_level",Integer,Positive);
    addKey("coarsesampler",String);
  }
  
  /** @brief Constructor for setting specific values
   *  @param[in] n_lvl Number of levels
   *  @param[in] coarse_spl Coarse level sampler
   */
  HierarchicalParameters(const unsigned int n_lvl,
                         const SamplerType coarse_spl) :
  Parameters("hierarchical") {
    n_level_ = n_lvl;
    coarsesampler_ = coarse_spl;
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_level_ = getContents("n_level")->getInt();
      std::string sampler_str = getContents("coarsesampler")->getString();
      if (sampler_str == "HMC") {
        coarsesampler_ = SamplerHMC;
      } else if (sampler_str == "cluster") {
        coarsesampler_ = SamplerCluster;
      } else if (sampler_str == "exact") {
        coarsesampler_ = SamplerExact;
      } else  {
        mpi_parallel::cerr << " ERROR: Unknown coarse sampler: " << sampler_str;
        mpi_parallel::cerr << std::endl;
        mpi_parallel::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
        mpi_exit(EXIT_FAILURE);
      }
    }
    return readSuccess;
  }

  /** @brief Return number of levels */
  unsigned int n_level() const { return n_level_; }
  /** @brief Return sampler type */
  SamplerType coarsesampler() const { return coarsesampler_; }
private:
  /** @brief Number of levels */
  unsigned int n_level_;
  /** @brief tolerance epsilon */
  SamplerType coarsesampler_;
};

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
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                      const GeneralParameters param_general,
                      const StatisticsParameters param_stats,
                      const HMCParameters param_hmc,
                      const ClusterParameters param_cluster,
                      const HierarchicalParameters param_hierarchical);
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
  
/** @brief Show statistics
     *
     * Print out statistics
     */
    virtual void show_stats();
  
  /** @brief Set current state to particular value
   *
   * @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<Path> x_path);
  
  /** Return cost per sample */
  virtual double cost_per_sample() { return cost_per_sample_; }
  
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
  /** @brief cost per sample */
  double cost_per_sample_;
};

#endif // HIERARCHICALSAMPLER
