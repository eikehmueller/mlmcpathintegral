#ifndef MONTECARLOSINGLELEVEL_HH
#define MONTECARLOSINGLELEVEL_HH MONTECARLOSINGLELEVEL_HH
#include <utility>
#include <cmath>
#include <iostream>
#include "sampler.hh"
#include "action.hh"
#include "quantityofinterest.hh"
#include "statistics.hh"
#include "parameters.hh"
#include "parameters.hh"
#include "hmcsampler.hh"
#include "clustersampler.hh"
#include "montecarlo.hh"

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
    n_samples_(100),
    sampler_(SamplerHMC) {
    addKey("n_burnin",Integer,Positive);
    addKey("n_samples",Integer,Positive);
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
      std::string sampler_str = getContents("sampler")->getString();
      if (sampler_str == "HMC") {
        sampler_ = SamplerHMC;
      } else if (sampler_str == "cluster") {
        sampler_ = SamplerCluster;
      } else if (sampler_str == "exact") {
        sampler_ = SamplerExact;
      } else  {
        std::cerr << " ERROR: Unknown sampler: " << sampler_str << std::endl;
        std::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
        exit(-1);
      }
    }
    return readSuccess;
  }

  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return number of samples */
  unsigned int n_samples() const { return n_samples_; }
  /** @brief Return sampler type */
  SamplerType sampler() const { return sampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
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
   */
  MonteCarloSingleLevel(std::shared_ptr<Action> action_,
                        std::shared_ptr<QoI> qoi_,
                        const GeneralParameters param_general,
                        const HMCParameters param_hmc,
                        const ClusterParameters param_cluster,
                        const SingleLevelMCParameters param_singlelevelmc);
  
  /** @brief Calculate QoI
   * 
   * Calculate the Quantity of interest by Monte Carlo sampling. Return
   * estimator for the mean and variance
   *
   * @param[inout] stats Object for recording statistics
   */
  void evaluate(Statistics& stats);

  /** @brief Return sampler */
  std::shared_ptr<Sampler> get_sampler() { return sampler; }
  

private:
  /** @brief Action action to use */
  std::shared_ptr<Action> action;
  /** @brief Sampler class for creating samples */
  std::shared_ptr<Sampler> sampler;
  /** @brief Quantity of interest */
  std::shared_ptr<QoI> qoi;
};

#endif // MONTECARLOSINGLELEVEL_HH