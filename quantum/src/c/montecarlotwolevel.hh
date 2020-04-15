#ifndef MONTECARLOTWOLEVEL_HH
#define MONTECARLOTWOLEVEL_HH MONTECARLOTWOLEVEL_HH
#include <utility>
#include <cmath>
#include <iostream>
#include "sampler.hh"
#include "action.hh"
#include "conditionedfineaction.hh"
#include "quantityofinterest.hh"
#include "twolevelmetropolisstep.hh"
#include "statistics.hh"
#include "parameters.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"
#include "montecarlo.hh"
#include "mpi_wrapper.hh"

/** @file montecarlotwolevel.hh
 * @brief Header file for two level Monte Carlo classes
 */

/** @class TwoLevelMCParameters
 *
 * @brief Class for storing parameters of two level Monte Carlo integrator.
 */
class TwoLevelMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  TwoLevelMCParameters() :
    Parameters("twolevelmc"),
    n_burnin_(100),
    n_samples_(100),
    coarsesampler_(SamplerHMC) {
    addKey("n_burnin",Integer,Positive);
    addKey("n_samples",Integer,Positive);
    addKey("coarsesampler",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_burnin_ = getContents("n_burnin")->getInt();
      n_samples_ = getContents("n_samples")->getInt();
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

  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return number of samples */
  unsigned int n_samples() const { return n_samples_; }
  /** @brief Return sampler type */
  SamplerType coarsesampler() const { return coarsesampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
  /** @brief Sampler type */
  SamplerType coarsesampler_;
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
   * @param[in] fine_action_ Action on fine level
   * @param[in] qoi_ Quantity of interest
   * @param[in] param_general General parameters
   * @param[in] param_hmc HMC sampler parameters
   * @param[in] param_cluster Cluster sampler parameters
   * @param[in] param_twolevelmc Two level sampler parameters
   */
  MonteCarloTwoLevel(std::shared_ptr<Action> fine_action_,
                     std::shared_ptr<QoI> qoi_,
                     const GeneralParameters param_general,
                     const HMCParameters param_hmc,
                     const ClusterParameters param_cluster,
                     const TwoLevelMCParameters param_twolevelmc);
  
  /** @brief Calculate mean and variance of difference
   *
   * Calculate the mean and variance of the difference in the QoI 
   * evaluated at two subsequent levels.
   *
   * @param[inout] stats_fine Statistics for fine level
   * @param[inout] stats_coarse Statistics for coarse level
   * @param[inout] stats_diff Statistics for difference
   */
  void evaluate_difference(Statistics& stats_fine,
                           Statistics& stats_coarse,
                           Statistics& stats_diff);

  /** @brief Return coarse sampler */
  std::shared_ptr<Sampler> get_coarsesampler() { return coarse_sampler; }
  
  /** @brief Return reference to two-level sampler 
   */
  std::shared_ptr<TwoLevelMetropolisStep> get_twolevelstep() {
    return twolevel_step;
  }
    
private:
  /** @brief Number of samples */
  const unsigned int n_samples;
  /** @brief Sampler on coarse level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Action on coarse level */
  std::shared_ptr<Action> coarse_action;
  /** @brief Action on fine level */
  std::shared_ptr<Action> fine_action;
  /** @brief Conditioned fine action */
  std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
  /** @brief Quantity of interest */
  std::shared_ptr<QoI> qoi;
  /** Two-level sampler */
  std::shared_ptr<TwoLevelMetropolisStep> twolevel_step;
};

#endif // MONTECARLOTWOLEVEL_HH
