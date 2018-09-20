#include "montecarlosinglelevel.hh"

/** @file montecarlosinglelevel.cc
 * @brief Implementation of montecarlosinglelevel.hh 
 */

/* Constructor */
MonteCarloSingleLevel::MonteCarloSingleLevel(std::shared_ptr<Action> action_,
                                             std::shared_ptr<QoI> qoi_,
                                             const GeneralParameters param_general,
                                             const HMCParameters param_hmc,
                                             const ClusterParameters param_cluster,
                                             const SingleLevelMCParameters param_singlelevelmc) :
  MonteCarlo(param_singlelevelmc.n_burnin()),
  action(action_), 
  qoi(qoi_) {
  if (param_singlelevelmc.sampler() == SamplerHMC) {
    sampler = std::make_shared<HMCSampler>(action,
                                           param_hmc.T(),
                                           param_hmc.dt(),
                                           param_hmc.n_burnin());
  } else if (param_singlelevelmc.sampler() == SamplerCluster) {
    if (param_general.action() != ActionRotor) {
      std::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
      exit(-1);
    }
    sampler = std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(action),
                                               param_cluster.n_burnin(),
                                               param_cluster.n_updates());
  } else if (param_singlelevelmc.sampler() == SamplerExact) {
    if (param_general.action() != ActionHarmonicOscillator) {
      std::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
      exit(-1);
    }
    sampler = std::dynamic_pointer_cast<Sampler>(action);
  }
  int k_max = 20;
  stats_corr = std::make_shared<Statistics>("corr",k_max);
  stats_Q = std::make_shared<Statistics>("Q",k_max);
  }

/** Calculate Monte Carlo estimate with single level method */
void MonteCarloSingleLevel::evaluate() {
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(action->getM_lat(),
                           action->getT_final());
  // Window over which autocorrelation is measured
  double epsilon=1.E-2;
  unsigned int n_min_samples_autocorr = 10;
  unsigned int n_min_samples_qoi = 10;

  for (int i=i;i<n_burnin;++i)
    sampler->draw(x_path);
  
  double two_epsilon_inv2 = 2./(epsilon*epsilon);
  stats_corr->reset();
  stats_Q->reset();
  int t=0;
  bool sufficient_stats = false;
  while (not sufficient_stats) {
    sampler->draw(x_path);
    /* Save (some) paths to disk? Edit file config.h */
#ifdef SAVE_PATHS
    if ( (SAVE_FIRST_PATH<=k) and (k<=SAVE_LAST_PATH) ) {
      std::stringstream filename;      
      filename << "path_";
      filename << std::setw(8) << std::setfill('0');
      filename << k << ".dat";
      x_path->save_to_disk(filename.str());
    }
#endif
    // Quantity of interest
    double qoi_Q = qoi->evaluate(x_path);
    stats_corr->record_sample(qoi_Q);
    if ( (t > stats_corr->tau_int()) and
         (stats_corr->samples() > n_min_samples_autocorr) ) {
      t = 0;
      stats_Q->record_sample(qoi_Q);
      int n_samples = stats_Q->samples();
      if (n_samples > n_min_samples_qoi) {
        int n_target = ceil(two_epsilon_inv2*stats_Q->variance());
        sufficient_stats = (n_samples > n_target);
      }
    }
    t++;
  }
}

/* Print out statistics */
void MonteCarloSingleLevel::show_statistics() {
  std::cout << *stats_corr;
  std::cout << *stats_Q;
  std::cout << std::endl;
}
