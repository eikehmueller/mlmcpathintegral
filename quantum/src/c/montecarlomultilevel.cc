#include "montecarlomultilevel.hh"

/** @file montecarlomultilevel.cc
 * @brief Implementation of montecarlomultilevel.hh 
 */

MonteCarloMultiLevel::MonteCarloMultiLevel(std::shared_ptr<Action> fine_action_,
                                           std::shared_ptr<QoI> qoi_,
                                           const GeneralParameters param_general,
                                           const StatisticsParameters param_stats,
                                           const HMCParameters param_hmc,
                                           const ClusterParameters param_cluster,
                                           const MultiLevelMCParameters param_multilevelmc) :
  MonteCarlo(param_multilevelmc.n_burnin()),
  fine_action(fine_action_),
  qoi(qoi_),
  n_level(param_multilevelmc.n_level()),
  epsilon(param_multilevelmc.epsilon()),
  n_target(param_multilevelmc.n_level(),0),
  t_indep(param_multilevelmc.n_level(),0.0),
  n_indep(param_multilevelmc.n_level(),0),
  t_sampler(param_multilevelmc.n_level(),0),
  n_autocorr_window(param_stats.n_autocorr_window()),
  n_min_samples_qoi(param_stats.n_min_samples_qoi()),
  n_min_samples_sampler(param_multilevelmc.n_min_samples_sampler()),
  timer("MultilevelMC") {
  // Check that Number of lattice points permits number of levels
  unsigned int M_lat = fine_action->getM_lat();
  if ( (M_lat>>n_level)<<n_level == M_lat) {
    std::cout << " M_lat = " << M_lat << " = 2^{" << n_level << "} * " << (M_lat>>n_level) << std::endl;
  } else {
    std::cout << "ERROR: M_lat = " << M_lat << " is not a multiple of 2^{n_level} = 2^{"<<n_level << "}" << std::endl;
  }
  action.push_back(fine_action);
  // Construct action and two-level MCMC step on all levels
  for (int ell=0;ell<n_level-1;++ell) {
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
  std::shared_ptr<Action> coarse_action = action[n_level-1];
  // Construct paths on all levels
  double T_final = fine_action->getT_final();
  for (int ell=0;ell<n_level;++ell) {
    unsigned int M_lat = action[ell]->getM_lat();    
    x_path.push_back(std::make_shared<Path>(M_lat,T_final));
    x_coarse_path.push_back(std::make_shared<Path>(M_lat,T_final));
    x_sampler_path.push_back(std::make_shared<Path>(M_lat,T_final));
  }
  // Construct sampler on coarsest level
  if (param_multilevelmc.coarsesampler() == SamplerHMC) {
    coarse_sampler = std::make_shared<HMCSampler>(coarse_action,
                                                  param_hmc.T(),
                                                  param_hmc.dt(),
                                                  param_hmc.n_burnin());
  } else if (param_multilevelmc.coarsesampler() == SamplerCluster) {
    if (param_general.action() != ActionRotor) {
      std::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
      exit(-1);
    }
    coarse_sampler =
      std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(coarse_action),
                                       param_cluster.n_burnin(),
                                       param_cluster.n_updates());
  } else if (param_multilevelmc.coarsesampler() == SamplerExact) {
    if (param_general.action() != ActionHarmonicOscillator) {
      std::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
      exit(-1);
    }
    coarse_sampler = std::dynamic_pointer_cast<Sampler>(coarse_action);
  }
  // Construct statistics on all levels
  for (int level=0;level<n_level;++level) {
    std::stringstream stats_sampler_label;
    stats_sampler_label << "Q_{sampler}[" << level << "]";
    stats_sampler.push_back(std::make_shared<Statistics>(stats_sampler_label.str(),n_autocorr_window));
    std::stringstream stats_qoi_label;
    stats_qoi_label << "Y[" << level << "]";
    stats_qoi.push_back(std::make_shared<Statistics>(stats_qoi_label.str(),n_autocorr_window));
  }
}    

void MonteCarloMultiLevel::evaluate() {
  /* Vector recording time since taking last sample from fine level chain
   * on a particular level. The samples are decorrelated of this time is
   * larger than the autocorrelation time of the fine level process.
   */
  for (int level=0;level<n_level;++level) {
    stats_sampler[level]->reset();
    t_sampler[level] = 0;
    n_indep[level] = 0;
    t_indep[level] = 0;
  }
  double two_epsilon_inv2 = 2./(epsilon*epsilon);

  // Burn in phase
  // Loop over all levels and create n_burnin samples
  for (int level=n_level-1;level>=0;level--) {
    for (int j=0;j<n_burnin;++j) {
      if (level==n_level-1) {
        draw_independent_sample(level,x_path[level]);
      } else {
        draw_independent_sample(level+1,x_coarse_path[level+1]);
        twolevel_step[level]->draw(x_coarse_path[level+1],x_path[level]);
      }
    }
  }
  std::cout << "Burnin completed" << std::endl;

  // Reset everything before actual run
  for (int level=0;level<n_level;++level) {
    stats_qoi[level]->reset();
    n_target[level] = n_min_samples_qoi;
  }

  timer.reset();
  timer.start();
  bool sufficient_stats=false;
  do {
    // Loop over all levels and measure QoI
    for (int level=n_level-1;level>=0;level--) {
      double qoi_Y;
      unsigned int j_start = stats_qoi[level]->samples();
      for (int j=j_start;j<n_target[level];++j) {
        if (level==n_level-1) {
          draw_independent_sample(level,x_path[level]);
          qoi_Y = qoi->evaluate(x_path[level]);
        } else {
          draw_independent_sample(level+1,x_coarse_path[level+1]);
          twolevel_step[level]->draw(x_coarse_path[level+1],x_path[level]);
          double qoi_fine = qoi->evaluate(x_path[level]);
          double qoi_coarse = qoi->evaluate(x_sampler_path[level+1]);
          qoi_Y = qoi_fine-qoi_coarse;
        }
        stats_qoi[level]->record_sample(qoi_Y);
      }
    }
    // Now recompute target numbers on all levels
    sufficient_stats = true;
    double sum_s_ell2_C_ell_eff = 0;
    for (int ell=0;ell<n_level;++ell) {
      double V_ell = stats_qoi[ell]->variance();
      double C_ell_eff = cost_eff(ell);
      sum_s_ell2_C_ell_eff += sqrt(V_ell*C_ell_eff);
    }
    for (int ell=0;ell<n_level;++ell) {
      int n_samples = stats_qoi[ell]->samples();
      // Calculate variance and predict new number of target samples
      double V_ell = stats_qoi[ell]->variance();
      double tau_int = stats_qoi[ell]->tau_int();
      double C_ell_eff = cost_eff(ell);
      n_target[ell] = ceil(two_epsilon_inv2*sum_s_ell2_C_ell_eff*sqrt(V_ell/C_ell_eff)*tau_int);
      sufficient_stats = sufficient_stats and (n_samples > n_target[ell]);
    }
  } while (not sufficient_stats); 
  timer.stop();
}

/* Draw an independent sample on a particular level */
void MonteCarloMultiLevel::draw_independent_sample(const int ell,
                                                   std::shared_ptr<Path> x_path) {
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
    if ( (stats_sampler[level]->samples() > n_min_samples_sampler) and
         (t_sampler[level] >= ceil(stats_sampler[level]->tau_int())) ) {
      // t_sampler[level] is the number of x_sampler_path samples
      // generated on this level since the last independent sample was used
      t_indep[level] = (n_indep[level]*t_indep[level]+t_sampler[level])/(1.0+n_indep[level]);
      // t_indep[level] is the average number of samples between indpendent
      // samples
      n_indep[level]++;
      // n_indep is the number of independent samples of x_sampler_path on
      // this level
      t_sampler[level] = 0; // Reset number of independent samples
      // Move to next-finer level if we have obtained a new independent sample
      level--;
    } else {
      // Return to coarsest level
      level = n_level-1;
    }
  } while (level>=ell);
  // Copy path
  std::copy(x_sampler_path[ell]->data,
            x_sampler_path[ell]->data+x_sampler_path[ell]->M_lat,
            x_path->data);
}

/* Calculate effective cost on a particular level ell */
unsigned int MonteCarloMultiLevel::cost_eff(const int ell) const {
  // Cost on current level
  int cost = action[ell]->evaluation_cost();
  // Add cost on coarser levels
  int T_k = 1;
  for (int k=ell+1;k<n_level;++k) {
    T_k *= ceil(t_indep[k]);
    cost += T_k*action[k]->evaluation_cost();
  }
  return cost*ceil(stats_qoi[ell]->tau_int());
}
  
/* Show detailed statistics on all levels */
void MonteCarloMultiLevel::show_detailed_statistics() {
  std::cout << "=== Statistics of coarse level samplers ===" << std::endl;
  for (int level=1;level<n_level;++level) {
    std::cout << "level = " << level << std::endl;
    std::cout << *stats_sampler[level];
    std::cout << " spacing between samples " << t_indep[level] << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "=== Statistics of QoI ===" << std::endl;
  for (int level=0;level<n_level;++level) {
    std::cout << "level = " << level << std::endl;
    std::cout << *stats_qoi[level];
    std::cout << " target number of samples = " << n_target[level] << std::endl;
    std::cout << " cost " << cost_eff(level) << std::endl;
    std::cout << "------------------------------------" << std::endl;
  }
  std::cout << std::endl;
}

/* Show summary statistics and timing */
void MonteCarloMultiLevel::show_statistics() {
  double qoi_value = 0.0;
  double qoi_error = 0.0;
  for (int level=0;level<n_level;++level) {
    qoi_value += stats_qoi[level]->average();
    double error = stats_qoi[level]->error();
    qoi_error += error*error;
  }
  qoi_error = sqrt(qoi_error);
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << std::endl;
  std::cout << " Q: Avg +/- Err = " << qoi_value << " +/- " << qoi_error << std::endl;
  std::cout << std::endl;
  std::cout << timer << std::endl;
  std::cout << std::endl;
}
