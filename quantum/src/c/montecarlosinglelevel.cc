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
  MonteCarlo(param_singlelevelmc.n_samples(),
             param_singlelevelmc.n_burnin()),
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
  }

/** Calculate Monte Carlo estimate with single level method */
void MonteCarloSingleLevel::evaluate(Statistics& stats) {
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(action->getM_lat(),
                           action->getT_final());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    sampler->draw(x_path);
  }
  stats.reset();
  for (unsigned int k=0;k<n_samples;++k) {
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
    stats.record_sample(qoi->evaluate(x_path));
  }
}
