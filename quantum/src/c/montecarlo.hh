#ifndef MONTECARLO_HH
#define MONTECARLO_HH MONTECARLO_HH
#include <utility>
#include <cmath>
#include <iostream>
#include "sampler.hh"
#include "action.hh"
#include "conditionedfineaction.hh"
#include "quantityofinterest.hh"
#include "twolevelmetropolissampler.hh"
#include "statistics.hh"
#include "parameters.hh"

/** @file montecarlo.hh
 * @brief Header file for Monte Carlo classes
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
        std::cerr << " ERROR: Unknown coarse sampler: " << sampler_str;
        std::cerr << std::endl;
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
  SamplerType coarsesampler() const { return coarsesampler_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of samples */
  unsigned int n_samples_;
  /** @brief Sampler type */
  SamplerType coarsesampler_;
};

/** @class MonteCarlo
 * 
 * @brief Monte Carlo base class
 */

class MonteCarlo {
public:
  /** @brief Create new instance
   *
   * @param[in] n_samples_ Number of samples
   * @param[in] n_burnin_ Number of burn-in steps
   * @param[in] record_stats_ Record statistics of sampler
   */
  MonteCarlo(const unsigned int n_samples_,
             const unsigned int n_burnin_,
             const bool record_stats_=false) :
    n_samples(n_samples_), 
    n_burnin(n_burnin_),
    record_stats(record_stats_) {}
protected:
  /** @brief Number of samples */
  const unsigned int n_samples;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
  /** @brief Record statistics */
  const bool record_stats;
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
   * @param[in] sampler_ Sampler to draw from
   * @param[in] qoi_ Quantity of interest to evaluate on samples
   * @param[in] n_samples_ Number of samples to evaluate
   * @param[in] n_burnin_ Number of burnin samples (QoI not evaluated on those)
   * @param[in] record_stats_ Record statistics of sampler
   */
  MonteCarloSingleLevel(Action& action_,
                        Sampler& sampler_,
                        QoI& qoi_,
                        const unsigned int n_samples_,
                        const unsigned int n_burnin_,
                        const bool record_stats_=false) :
    MonteCarlo(n_samples_,n_burnin_,record_stats_), 
    action(action_), 
    sampler(sampler_), 
    qoi(qoi_)
  {}

  /** @brief Calculate QoI
   * 
   * Calculate the Quantity of interest by Monte Carlo sampling. Return
   * estimator for the mean and variance
   *
   * @param[inout] stats Object for recording statistics
   */
  void evaluate(Statistics& stats);
private:
  /** @brief Action action to use */
  Action& action;
  /** @brief Sampler class for creating samples */
  Sampler& sampler;
  /** @brief Quantity of interest */
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
   * @param[in] coarse_action_ Action on coarse level
   * @param[in] coarse_sampler_ Sampler on coarse level
   * @param[in] fine_action_ Action on fine level
   * @param[in] conditioned_fine_action_ Conditioned fine action
   * @param[in] qoi_ Quantity of interest
   * @param[in] n_samples_ Number of samples to evaluate
   * @param[in] n_burnin_ Number of burnin samples (QoI not evaluated on those)
   * @param[in] record_stats_ Record statistics of sampler
   */
  MonteCarloTwoLevel(Action& coarse_action_,
                     Sampler& coarse_sampler_,
                     Action& fine_action_,
                     ConditionedFineAction& conditioned_fine_action_,
                     QoI& qoi_,
                     const unsigned int n_samples_,
                     const unsigned int n_burnin_,
                     const bool record_stats_=false) :
    MonteCarlo(n_samples_,n_burnin_,record_stats_),
    coarse_sampler(coarse_sampler_),
    coarse_action(coarse_action_),
    fine_action(fine_action_),
    conditioned_fine_action(conditioned_fine_action_),
    qoi(qoi_),
    twolevel_sampler(coarse_sampler_,
                     coarse_action_,
                     fine_action_,
                     conditioned_fine_action_,
                     record_stats_) {}

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

  /** @brief Return reference to two-level sampler 
   */
  TwoLevelMetropolisSampler& get_twolevelsampler() {
    return twolevel_sampler;
  }
    
private:
  /** @brief Sampler on coarse level */
  Sampler& coarse_sampler;
  /** @brief Action on coarse level */
  Action& coarse_action;
  /** @brief Action on fine level */
  Action& fine_action;
  /** @brief Conditioned fine action */
  ConditionedFineAction& conditioned_fine_action;
  /** @brief Quantity of interest */
  QoI& qoi;
  /** Two-level sampler */
  TwoLevelMetropolisSampler twolevel_sampler;
};

#endif // MONTECARLO_HH
