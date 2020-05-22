#ifndef MONTECARLOMULTILEVEL_HH
#define MONTECARLOMULTILEVEL_HH MONTECARLOMULTILEVEL_HH
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <memory>
#include <typeinfo>
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
#include "multilevelsampler.hh"
#include "montecarlo.hh"
#include "mpi_wrapper.hh"

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
    n_burnin_(100),
    epsilon_(1.0),
    show_detailed_stats_(false),
    sampler_(SamplerMultilevel) {
    addKey("n_level",Integer,Positive);
    addKey("n_burnin",Integer,Positive);
    addKey("epsilon",Double,Positive);
    addKey("sampler",String);
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
      n_burnin_ = getContents("n_burnin")->getInt();
      epsilon_ = getContents("epsilon")->getDouble();
      show_detailed_stats_ = getContents("show_detailed_stats")->getBool();
      std::string sampler_str = getContents("sampler")->getString();
      if (sampler_str == "hierarchical") {
        sampler_ = SamplerHierarchical;
      } else if (sampler_str == "multilevel") {
          sampler_ = SamplerMultilevel;
      } else {
        mpi_parallel::cerr << " ERROR: Unknown sampler: " << sampler_str << std::endl;
        mpi_parallel::cerr << "        allowed values are \'HMC\', \'cluster\', \'exact\'" << std::endl;
        mpi_exit(EXIT_FAILURE);
      }
    }
    return readSuccess;
  }

  /** @brief Return number of levels */
  unsigned int n_level() const { return n_level_; }
  /** @brief Return number burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return tolerance epsilon */
  double epsilon() const { return epsilon_; }
  /** @brief Show detailed statistics? */
  bool show_detailed_stats() const { return show_detailed_stats_; }
  /** @brief Return sampler type */
  SamplerType sampler() const { return sampler_; }

private:
  /** @brief Number of levels */
  unsigned int n_level_;
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief tolerance epsilon */
  double epsilon_;
  /** @brief Show detailed statistics? */
  bool show_detailed_stats_;
  /** @brief Sampler type */
  SamplerType sampler_;
};

/** @class MonteCarloMultiLevel
 * 
 * @brief Multilevel Monte Carlo method
 * 
 * 
 */
class MonteCarloMultiLevel : public MonteCarlo {
public:
  /** @brief Create new instance 
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
                       const MultiLevelMCParameters param_multilevelmc,
                       const HierarchicalParameters param_hierarchical);

  /** @brief Run multilevel method */
  void evaluate();

  /** @brief Print out detailed statistics on all levels */
  void show_detailed_statistics();
  
  /** @brief Print out statistics */
  void show_statistics();
  
  /** @brief Return numerical result */
  double numerical_result() const;
  
  /** @brief Return statistical error */
  double statistical_error() const;
  
private:

  /** Draw (independent) coarse level sample
   * @param[in] level Level on which to draw sample (i.e. the coarse level)
   * @param[out] x_path Path which will contain the sample
   */
  void draw_coarse_sample(const unsigned int level,
                          std::shared_ptr<Path> x_path);
  
  /** @brief Calculate effective cost \f$C_{\ell}^{eff}\f$ on a level
   * 
   * Calculate
   * \f[
   *   C^{eff}_{\ell} = \lceil \tau^{int}_{\ell} \rceil (C_\ell + C_{\ell+1}^{indep})
   * \f]
   *
   * @param[in] ell Level \f$\ell\f$
   */
  double cost_eff(const int ell) const;
  
private:
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
  /** @brief Coarse path on a particular level */
  std::vector<std::shared_ptr<Path> > x_coarse_path;
  /** @brief vector with statistics of uncorrelated Y's*/
  std::vector<std::shared_ptr<Statistics> > stats_qoi;
  /** @brief hierarchical sampler on each level */
  std::vector<std::shared_ptr<Sampler>> coarse_sampler;
  /** @brief target number of samples on each level */
  std::vector<int> n_target;
  /** @brief Size of autocorrelation window */
  unsigned int n_autocorr_window;
  /** @brief Minimal number of samples for qoi */
  unsigned int n_min_samples_qoi;
  /** @brief Timer class */
  Timer timer;
  /** @brief number of skipped samples between independent coarse samples on all levels */
  std::vector<double> t_indep;
  /** @brief number of independent coarse samples on all levels */
  std::vector<int> n_indep;
  /** @brief number of samples generated since last independent coarse sample */
  std::vector<int> t_sampler;
  /** @brief vector with statistics of coarse sampler paths */
  std::vector<std::shared_ptr<Statistics> > stats_coarse_sampler;
  /** @brief Sub-sample coarse level sampler? */
  bool sub_sample_coarse;
};

#endif // MONTECARLOMULTILEVEL_HH
