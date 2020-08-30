#include "montecarlotwolevel.hh"

/** @file montecarlotwolevel.cc
 * @brief Implementation of montecarlotwolevel.hh 
 */
/* Constructor */
MonteCarloTwoLevel::MonteCarloTwoLevel(const std::shared_ptr<Action> fine_action_,
                                       const std::shared_ptr<QoI> qoi_,
                                       const std::shared_ptr<SamplerFactory> sampler_factory,
                                       const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                       const TwoLevelMCParameters param_twolevelmc) :
  MonteCarlo(param_twolevelmc.n_burnin()),
  n_samples(param_twolevelmc.n_samples()),
  fine_action(fine_action_),
  qoi(qoi_) {
  coarse_action = fine_action->coarse_action();
  coarse_sampler = sampler_factory->get(coarse_action);
  conditioned_fine_action = conditioned_fine_action_factory->get(fine_action);
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
  stats_coarse.hard_reset();
  stats_fine.hard_reset();
  stats_diff.hard_reset();

  // Burn-in phase
  for (unsigned int k=0;k<n_burnin;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step->draw(x_coarse_path,x_path);
    double qoi_fine = qoi->evaluate(x_path);
    double qoi_coarse = qoi->evaluate(x_coarse_path);
    stats_fine.record_sample(qoi_fine);
    stats_coarse.record_sample(qoi_coarse);
    stats_diff.record_sample(qoi_fine-qoi_coarse);
  }
  mpi_parallel::cout << "Burnin completed" << std::endl;

  unsigned int n_local_samples = distribute_n(n_samples);
  // Sampling phase
  // Do a hard reset since we are interested in the variance
  stats_coarse.hard_reset();
  stats_fine.hard_reset();
  stats_diff.hard_reset();
  for (unsigned int k=0;k<n_local_samples;++k) {
    coarse_sampler->draw(x_coarse_path);
    twolevel_step->draw(x_coarse_path,x_path);
    double qoi_fine = qoi->evaluate(x_path);
    double qoi_coarse = qoi->evaluate(x_coarse_path);
    stats_fine.record_sample(qoi_fine);
    stats_coarse.record_sample(qoi_coarse);
    stats_diff.record_sample(qoi_fine-qoi_coarse);
  }
}
