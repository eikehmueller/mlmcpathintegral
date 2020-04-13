#include "montecarlotwolevel.hh"
#include "config.h"
#include <sstream>
#include <iomanip>

/** @file montecarlotwolevel.cc
 * @brief Implementation of montecarlotwolevel.hh 
 */
/* Constructor */
MonteCarloTwoLevel::MonteCarloTwoLevel(std::shared_ptr<Action> fine_action_,
                                       std::shared_ptr<QoI> qoi_,
                                       const GeneralParameters param_general,
                                       const HMCParameters param_hmc,
                                       const ClusterParameters param_cluster,
                                       const TwoLevelMCParameters param_twolevelmc) :
  MonteCarlo(param_twolevelmc.n_burnin()),
  n_samples(param_twolevelmc.n_samples()),
  fine_action(fine_action_),
  qoi(qoi_) {
  coarse_action = fine_action->coarse_action();
  if (param_twolevelmc.coarsesampler() == SamplerHMC) {
    coarse_sampler = std::make_shared<HMCSampler>(coarse_action,
                                                  param_hmc.T(),
                                                  param_hmc.dt(),
                                                  param_hmc.n_burnin());
  } else if (param_twolevelmc.coarsesampler() == SamplerCluster) {
    if (param_general.action() != ActionRotor) {
      mpi_parallel::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
      exit(-1);
    }
    coarse_sampler =
      std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(coarse_action),
                                       param_cluster.n_burnin(),
                                       param_cluster.n_updates());
  } else if (param_twolevelmc.coarsesampler() == SamplerExact) {
    if (param_general.action() != ActionHarmonicOscillator) {
      mpi_parallel::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
      exit(-1);
    }
    coarse_sampler = std::dynamic_pointer_cast<Sampler>(coarse_action);
  }
  if (param_general.action() == ActionRotor) {
    conditioned_fine_action =
      std::make_shared<RotorConditionedFineAction>(std::dynamic_pointer_cast<RotorAction>(fine_action));
  } else {
    conditioned_fine_action =
      std::make_shared<GaussianConditionedFineAction>(fine_action);
  }
  twolevel_step = std::make_shared<TwoLevelMetropolisStep>(coarse_action,
                                                           fine_action,
                                                           conditioned_fine_action);
}

/** Calculate mean and variance of difference in QoI */
void MonteCarloTwoLevel::evaluate_difference(Statistics& stats_fine,
                                             Statistics& stats_coarse,
                                             Statistics& stats_diff) {
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(fine_action->getM_lat(),
                           fine_action->getT_final());
  std::shared_ptr<Path> x_coarse_path =
    std::make_shared<Path>(coarse_action->getM_lat(),
                           coarse_action->getT_final());

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step->draw(x_coarse_path,x_path);
  }

  // Sampling phase
  stats_diff.reset();
  for (unsigned int k=0;k<n_samples;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step->draw(x_coarse_path,x_path);
    double qoi_fine = qoi->evaluate(x_path);
    double qoi_coarse = qoi->evaluate(x_coarse_path);
    stats_fine.record_sample(qoi_fine);
    stats_coarse.record_sample(qoi_coarse);
    stats_diff.record_sample(qoi_fine-qoi_coarse);
  }
}
