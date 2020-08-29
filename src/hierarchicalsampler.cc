#include "hierarchicalsampler.hh"

/** @file hierarchicalsampler.cc
 * @brief Implementation of hierarchicalsampler.hh
 */

/* Construct new instance */
HierarchicalSampler::HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                                         const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
                                         const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                         const HierarchicalParameters param_hierarchical) :
  Sampler(),
  n_level(param_hierarchical.n_max_level()-fine_action->get_coarsening_level()),
  cost_per_sample_(0.0) {
  // Check that Number of lattice points permits number of levels
  unsigned int M_lat = fine_action->getM_lat();
  if ( (M_lat>>n_level)<<n_level == M_lat) {
    mpi_parallel::cout << " Hierarchical sampler: M_lat = " << M_lat << " = 2^{" << n_level << "-1} * " << (M_lat>>(n_level-1)) << std::endl;
  } else {
    mpi_parallel::cout << "ERROR: M_lat = " << M_lat << " is not a multiple of 2^{n_level} = 2^{"<<n_level << "}" << std::endl;
  }
  action.push_back(fine_action);
  // Construct action and two-level MCMC step on all levels
  for (unsigned int ell=0;ell<n_level-1;++ell) {
    std::shared_ptr<Action> action_tmp = action[ell];
    std::shared_ptr<Action> coarse_action_tmp = action[ell]->coarse_action();
    action.push_back(coarse_action_tmp);
    std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
    conditioned_fine_action = conditioned_fine_action_factory->get(action_tmp);
    std::shared_ptr<TwoLevelMetropolisStep> twolevel_step_tmp =
      std::make_shared<TwoLevelMetropolisStep>(coarse_action_tmp,
                                               action_tmp,
                                               conditioned_fine_action);
    twolevel_step.push_back(twolevel_step_tmp);
  }
  // Action on coarsest level
  std::shared_ptr<Action> coarse_action = action[n_level-1];
  double T_final = fine_action->getT_final();
  for (unsigned int ell=0;ell<n_level;++ell) {
    unsigned int M_lat = action[ell]->getM_lat();    
    x_sampler_path.push_back(std::make_shared<Path>(M_lat,T_final));
  }
  // Construct sampler on coarsest level
  coarse_sampler = coarse_sampler_factory->get(coarse_action);
  std::shared_ptr<Path> meas_path=std::make_shared<Path>(fine_action->getM_lat(),fine_action->getT_final());
  Timer timer_meas;
  unsigned int n_meas = 10000;
    timer_meas.start();
  for (unsigned int k=0;k<n_meas;++k) {
    draw(meas_path);
  }
  timer_meas.stop();
  cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
}

void HierarchicalSampler::draw(std::shared_ptr<Path> x_path) {
  accept = true;
  for (int ell=1;ell<n_level;++ell) {
    x_sampler_path[ell]->copy_strided(x_sampler_path[ell-1],2);
  }
  for (int ell=n_level-1;ell>=0;--ell) {
    if (ell == (n_level-1)) {
      /* Sample directly on coarsest level */
      coarse_sampler->set_state(x_sampler_path[ell]);
      coarse_sampler->draw(x_sampler_path[ell]);
      accept = accept and coarse_sampler->accepted();
    } else {
      twolevel_step[ell]->set_state(x_sampler_path[ell]);
      twolevel_step[ell]->draw(x_sampler_path[ell+1],
                               x_sampler_path[ell]);
      accept = accept and twolevel_step[ell]->accepted();
    }
    if (not accept) break;
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  if (accept or copy_if_rejected) {
    x_path->copy(x_sampler_path[0]);
  }
}

/* Set current state */
void HierarchicalSampler::set_state(std::shared_ptr<Path> x_path) {
  const unsigned int M_lat = x_path->M_lat;
  x_sampler_path[0]->copy(x_path);
}

/* Show statistics on all levels */
void HierarchicalSampler::show_stats() {
  mpi_parallel::cout << std::setprecision(3) << std::fixed;
  mpi_parallel::cout << "  cost per sample = " << cost_per_sample() << " mu s" << std::endl << std::endl;
  mpi_parallel::cout << "  acceptance rate = " << p_accept() << std::endl;
  for (unsigned int ell=0;ell<n_level;++ell) {
    double p_acc;
    if (ell==n_level-1) {
      p_acc = coarse_sampler->p_accept();
    } else {
      p_acc = twolevel_step[ell]->p_accept();
    }
    std::string level_str;
    if (ell == 0) {
      level_str = "[finest]  ";
    } else if (ell == n_level-1) {
      level_str = "[coarsest]";
    } else {
      level_str = "          ";
    }
    std::stringstream sstream;
    sstream.setf(std::ios::fixed);
    sstream.width(2);
    sstream.precision(3);
    sstream << "  level " << ell << " " << level_str <<" : ";
    sstream << " p = " << p_acc << std::endl;
    mpi_parallel::cout << sstream.str();
  }
}
