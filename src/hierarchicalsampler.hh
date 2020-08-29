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
    n_max_level_(2),
  coarsesampler_(SamplerHMC) {
    addKey("n_max_level",Integer,Positive);
    addKey("coarsesampler",String);
  }
  
  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_max_level_ = getContents("n_max_level")->getInt();
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

  /** @brief Return maximal number of levels */
  unsigned int n_max_level() const { return n_max_level_; }
  /** @brief Return sampler type */
  SamplerType coarsesampler() const { return coarsesampler_; }
private:
  /** @brief Number of levels */
  unsigned int n_max_level_;
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
   *
   * @param[in] fine_action_ Fine level Action to sample from
   * @param[in] coarse_sampler_factory Sampler factory for generating coarse level sampler
   * @param[in] param_general General parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                      const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
                      const GeneralParameters param_general,
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

class HierarchicalSamplerFactory : public SamplerFactory {
  public:
  /** @brief Create new instance
   *
   * @param[in] coarse_sampler_factory_ Factory for coarse level sampler
   * @param[in] param_hierarchical General parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
   */
  HierarchicalSamplerFactory(const std::shared_ptr<SamplerFactory> coarse_sampler_factory_,
                             const GeneralParameters param_general_,
                             const HierarchicalParameters param_hierarchical_) :
    coarse_sampler_factory(coarse_sampler_factory_),
    param_general(param_general_),
    param_hierarchical(param_hierarchical_) {}
  
  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::make_shared<HierarchicalSampler>(action,
                                                 coarse_sampler_factory,
                                                 param_general,
                                                 param_hierarchical);
  }
private:
  /** Factory for coarsest level sampler*/
  const std::shared_ptr<SamplerFactory> coarse_sampler_factory;
  /** General parameters */
  const GeneralParameters param_general;
  /** Hierarchical sampler parameters */
  const HierarchicalParameters param_hierarchical;
};

#endif // HIERARCHICALSAMPLER
