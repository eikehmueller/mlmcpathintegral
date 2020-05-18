#include "multilevelsampler.hh"

/** @file multilevelsampler.cc
 * @brief Implementation of multilevelsampler.hh
 */

/* Construct new instance */
MultilevelSampler::MultilevelSampler(const std::shared_ptr<Action> fine_action,
                                     const std::shared_ptr<QoI> qoi_,
                                     const GeneralParameters param_general,
                                     const StatisticsParameters param_stats,
                                     const HMCParameters param_hmc,
                                     const ClusterParameters param_cluster,
                                     const HierarchicalParameters param_hierarchical) :
  Sampler(),
  qoi(qoi_),
  n_level(param_hierarchical.n_level()),
  t_indep(param_hierarchical.n_level(),0.0),
  n_indep(param_hierarchical.n_level(),0),
  t_sampler(param_hierarchical.n_level(),0),
  n_autocorr_window(param_stats.n_autocorr_window()),
  cost_per_sample_(0.0){
  
  // Check that Number of lattice points permits number of levels
  unsigned int M_lat = fine_action->getM_lat();
  if ( (M_lat>>n_level)<<n_level == M_lat) {
    mpi_parallel::cout << " Multilevel sampler: M_lat = " << M_lat << " = 2^{" << n_level << "-1} * " << (M_lat>>(n_level-1)) << std::endl;
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
    if (param_general.action() == ActionRotor) {
      conditioned_fine_action =
        std::make_shared<RotorConditionedFineAction>(std::dynamic_pointer_cast<RotorAction>(action_tmp));
    } else {
      conditioned_fine_action =
        std::make_shared<GaussianConditionedFineAction>(action_tmp);
    } 
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
  if (param_hierarchical.coarsesampler() == SamplerHMC) {
    coarse_sampler = std::make_shared<HMCSampler>(coarse_action,
                                                  param_hmc);
  } else if (param_hierarchical.coarsesampler() == SamplerCluster) {
    if (param_general.action() != ActionRotor) {
      mpi_parallel::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    coarse_sampler =
      std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(coarse_action),
                                       param_cluster.n_burnin(),
                                       param_cluster.n_updates());
  } else if (param_hierarchical.coarsesampler() == SamplerExact) {
    if (param_general.action() != ActionHarmonicOscillator) {
      mpi_parallel::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    coarse_sampler = std::dynamic_pointer_cast<Sampler>(coarse_action);
  }
  // Statistics on all levels
  for (unsigned int level=0;level<n_level;++level) {
    std::stringstream stats_sampler_label;
    stats_sampler_label << "   Q_{sampler}[" << level << "]";
    stats_sampler.push_back(std::make_shared<Statistics>(stats_sampler_label.str(),n_autocorr_window));
  }
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

/* Draw next sample */
void MultilevelSampler::draw(std::shared_ptr<Path> x_path) {
  accept = true;
  int level = n_level-1;
  do {
    if (level == (n_level-1)) {
      /* Sample directly on coarsest level */
      coarse_sampler->draw(x_sampler_path[level]);
    } else {
      /*
       * On all other levels, sample by using the two level MCMC process
       * We know that the coarse level samples are decorrelated, since the
       * algorithm only ever proceeds to the next finer level if this is the
       * case.
       */
      twolevel_step[level]->draw(x_sampler_path[level+1],
                                 x_sampler_path[level]);
    }
    // The QoI of the independent sampler, Q_{ell}
    double qoi_sampler = qoi->evaluate(x_sampler_path[level]);
    stats_sampler[level]->record_sample(qoi_sampler);
    t_sampler[level]++;
    if (t_sampler[level] >= ceil(stats_sampler[level]->tau_int())) {
      // t_sampler[level] is the number of x_sampler_path samples
      // generated on this level since the last independent sample was used
      t_indep[level] = (n_indep[level]*t_indep[level]+t_sampler[level])/(1.0+n_indep[level]);
      // t_indep[level] is the average number of samples between indpendent
      // samples
      n_indep[level]++;
      // n_indep is the number of independent samples of x_sampler_path on
      // this level
      t_sampler[level] = 0; // Reset number of independent samples
      // Move to next-finer level since we have obtained a new independent
      // sample
      level--;
    } else {
      // Return to coarsest level
      level = n_level-1;
    }
  } while (level>=0);
  // Copy path
  x_path->copy(x_sampler_path[0]);
}

/* Set current state */
void MultilevelSampler::set_state(std::shared_ptr<Path> x_path) {
  x_sampler_path[0]->copy(x_path);
}

/* Show statistics on all levels */
void MultilevelSampler::show_stats() {
  mpi_parallel::cout << std::setprecision(3) << std::fixed;
  mpi_parallel::cout << "  cost per sample = " << cost_per_sample() << " mu s" << std::endl << std::endl;
  for (unsigned int ell=0;ell<n_level;++ell) {
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
    sstream << " level " << ell << " " << level_str <<" :";
    mpi_parallel::cout << sstream .str() << std::endl;
    if (ell==n_level-1) {
      coarse_sampler->show_stats();
    } else {
      twolevel_step[ell]->show_stats();
    }
    mpi_parallel::cout << "    average spacing between sample " << t_indep[ell] << std::endl;
    mpi_parallel::cout << *stats_sampler[ell] << std::endl;
  }
}
