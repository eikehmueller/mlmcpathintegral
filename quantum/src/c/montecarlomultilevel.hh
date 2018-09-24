#ifndef MONTECARLOMULTILEVEL_HH
#define MONTECARLOMULTILEVEL_HH MONTECARLOMULTILEVEL_HH
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include "timer.hh"
#include "path.hh"
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

/** @file montecarlomultilevel.hh
 * @brief Header file for multilevel Monte Carlo classes
 */

/** @class MultiLevelMCParameters
 *
 * @brief Class for storing parameters of multilevel Monte Carlo integrator.
 */
class MultiLevelMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  MultiLevelMCParameters() :
    Parameters("multilevelmc"),
    n_level_(2),
    epsilon_(1.0),
    coarsesampler_(SamplerHMC),
    show_detailed_stats_(false) {
    addKey("n_level",Integer,Positive);
    addKey("epsilon",Double,Positive);
    addKey("coarsesampler",String);
    addKey("show_detailed_stats",Bool);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_level_ = getContents("n_level")->getInt();
      epsilon_ = getContents("epsilon")->getDouble();
      show_detailed_stats_ = getContents("show_detailed_stats")->getBool();
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

  /** @brief Return number of levels */
  unsigned int n_level() const { return n_level_; }
  /** @brief Return tolerance epsilon */
  double epsilon() const { return epsilon_; }
  /** @brief Return sampler type */
  SamplerType coarsesampler() const { return coarsesampler_; }
  /** @brief Show detailed statistics? */
  bool show_detailed_stats() const { return show_detailed_stats_; }
private:
  /** @brief Number of levels */
  unsigned int n_level_;
  /** @brief tolerance epsilon */
  double epsilon_;
  /** @brief Sampler type */
  SamplerType coarsesampler_;
  /** @brief Show detailed statistics? */
  bool show_detailed_stats_;
};

/** @class MonteCarloMultiLevel
 * 
 * @brief Multilevel Monte Carlo method
 * 
 * 
 */
class MonteCarloMultiLevel : public MonteCarlo {
public:
  /** \brief Create new instance 
   *
   * @param[in] fine_action_ Action on fine level
   * @param[in] qoi_ Quantity of interest
   * @param[in] param_general General parameters
   * @param[in] param_hmc HMC sampler parameters
   * @param[in] param_cluster Cluster sampler parameters
   * @param[in] param_multilevelmc Multilevel parameters
   */
  MonteCarloMultiLevel(std::shared_ptr<Action> fine_action_,
                       std::shared_ptr<QoI> qoi_,
                       const GeneralParameters param_general,
                       const StatisticsParameters param_stats,
                       const HMCParameters param_hmc,
                       const ClusterParameters param_cluster,
                       const MultiLevelMCParameters param_multilevelmc);

  /** @brief Run multilevel method */
  void evaluate();

  /** @brief Print out detailed statistics on all levels */
  void show_detailed_statistics();
  
  /** @brief Print out statistics */
  void show_statistics();
  
private:
  /** @brief Sampler on coarsest level */
  std::shared_ptr<Sampler> coarse_sampler;
  /** @brief Action on fine level */
  std::shared_ptr<Action> fine_action;
  /** @brief Action on all levels of the multigrid hierarchy */
  std::vector<std::shared_ptr<Action> > action;
  /** @brief Two level step on all levels of the multigrid hierarchy */
  std::vector<std::shared_ptr<TwoLevelMetropolisStep> > twolevel_step;
  /** @brief Quantity of interest */
  std::shared_ptr<QoI> qoi;
  /** @brief Number of levels */
  const unsigned int n_level;
  /** @brief Tolerance epsilon */
  const double epsilon;
  /** @brief Path on a particular level */
  std::vector<std::shared_ptr<Path> > x_path;
  /** @brief Coarse on a particular level */
  std::vector<std::shared_ptr<Path> > x_coarse_path;
  /** @brief vector with statistics of (correlated) Q_fine */
  std::vector<std::shared_ptr<Statistics> > stats_corr;
  /** @brief vector with statistics of uncorrelated Y's*/
  std::vector<std::shared_ptr<Statistics> > stats_qoi;
  /** @brief Size of autocorrelation window */
  unsigned int n_autocorr_window;
  /** @brief Minimal number of samples for correlated quantities */
  unsigned int n_min_samples_corr;
  /** @brief Minimal number of samples for qoi */
  unsigned int n_min_samples_qoi;
  /** @brief Timer class */
  Timer timer;
};

#endif // MONTECARLOMULTILEVEL_HH
