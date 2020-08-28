#ifndef MONTECARLOSINGLELEVEL_HH
#define MONTECARLOSINGLELEVEL_HH MONTECARLOSINGLELEVEL_HH
#include <utility>
#include <cmath>
#include <iostream>
#include "config.h"
#include "timer.hh"
#include "sampler.hh"
#include "action.hh"
#include "quantityofinterest.hh"
#include "statistics.hh"
#include "parameters.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"
#include "montecarlo.hh"
#include "mpi_wrapper.hh"
#include "hierarchicalsampler.hh"
#include "multilevelsampler.hh"

/** @file montecarlosinglelevel.hh
 * @brief Header file for single level Monte Carlo classes
 */

/** @class SingleLevelMCParameters
 *
 * @brief Class for storing parameters of single level Monte Carlo integrator.
 */
class SingleLevelMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  SingleLevelMCParameters() :
    Parameters("singlelevelmc"),
    n_burnin_(100),
    epsilon_(1.0),
    sampler_(SamplerHMC) {
    addKey("n_burnin",Integer,Positive);
    addKey("n_samples",Integer,NonNegative);
    addKey("epsilon",Double,Positive);
    addKey("sampler",String);
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
      epsilon_ = getContents("epsilon")->getDouble();
      std::string sampler_str = getContents("sampler")->getString();
      if (sampler_str == "HMC") {
        sampler_ = SamplerHMC;
      } else if (sampler_str == "cluster") {
        sampler_ = SamplerCluster;
      } else if (sampler_str == "exact") {
        sampler_ = SamplerExact;
      } else if (sampler_str == "hierarchical") {
        sampler_ = SamplerHierarchical;
      } else if (sampler_str == "multilevel") {
          sampler_ = SamplerMultilevel;
      } else  {
        mpi_parallel::cerr << " ERROR: Unknown sampler: " << sampler_str << std::endl;
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
  /** @brief Return tolerance epsilon */
  double epsilon() const { return epsilon_; }
  /** @brief Return sampler type */
  SamplerType sampler() const { return sampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
  /** @brief tolerance epsilon */
  double epsilon_;
  /** @brief Sampler type */
  SamplerType sampler_;
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
   * @param[in] action_ Action to use
   * @param[in] qoi_ Quantity of interest to evaluate on samples
   * @param[in] param_general General parameters
   * @param[in] param_hmc HMC sampler parameters
   * @param[in] param_cluster Cluster sampler parameters
   * @param[in] param_singlelevelmc Single level sampler parameters
   * @param[in] param_hierarchical Hierarchical sampler parameters
  */ 
  MonteCarloSingleLevel(std::shared_ptr<Action> action_,
                        std::shared_ptr<QoI> qoi_,
                        std::shared_ptr<Sampler> sampler_,
                        const StatisticsParameters param_stats,
                        const SingleLevelMCParameters param_singlelevelmc);
  
  /** @brief Calculate QoI
   * 
   * Calculate the Quantity of interest by Monte Carlo sampling. Return
   * estimator for the mean and variance
   *
   * @param[inout] stats Object for recording statistics
   */
  void evaluate();

  /** @brief Show statistics */
  void show_statistics();

  /** @brief Return sampler */
  std::shared_ptr<Sampler> get_sampler() { return sampler; }
  
  /** @brief Return numerical result */
  double numerical_result() const {
    return stats_Q->average();
  }
  
  /** @brief Return statistical error */
  double statistical_error() const {
    return stats_Q->error();
  }

private:
  /** @brief Action action to use */
  std::shared_ptr<Action> action;
  /** @brief Sampler class for creating samples */
  std::shared_ptr<Sampler> sampler;
  /** @brief Quantity of interest */
  std::shared_ptr<QoI> qoi;
  /** @brief Statistics for measuring QoI */
  std::shared_ptr<Statistics> stats_Q;
  /** @brief Size of autocorrelation window */
  unsigned int n_autocorr_window;
  /** @brief Minimal number of samples for qoi */
  unsigned int n_min_samples_qoi;
  /** @brief Exact number of samples to use (ignored if zero) */
  unsigned int n_samples;
  /** @brief Tolerance epsilon */
  const double epsilon;
  /** @brief time */
  Timer timer;
};

#endif // MONTECARLOSINGLELEVEL_HH
