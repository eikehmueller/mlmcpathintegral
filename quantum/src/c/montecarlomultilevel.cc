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
  n_autocorr_window(param_stats.n_autocorr_window()),
  n_min_samples_corr(param_stats.n_min_samples_corr()),
  n_min_samples_qoi(param_stats.n_min_samples_qoi()),
  t_indep(param_multilevelmc.n_level(),0),
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
  std::vector<int> t_sampler(n_level,0);
  std::vector<int> n_indep(n_level,0);
  // Vector with target samples on each level. Record at least 100 samples.
  for (int level=0;level<n_level;++level) {
    stats_sampler[level]->reset();
    stats_qoi[level]->reset();
  }
  double two_epsilon_inv2 = 2./(epsilon*epsilon);
  int level = n_level-1; // Current level
  // Samples on a given level
  std::vector<int> n_samples_burnin(n_level,0);
  // Burnin phase
  do {
    // The QoI for the coarse path sampler on the current level
    double qoi_sampler;
    if (level == (n_level-1)) {
      coarse_sampler->draw(x_path[level]);
      coarse_sampler->draw(x_sampler_path[level]);
      qoi_sampler = qoi->evaluate(x_sampler_path[level]);
    } else {
      twolevel_step[level]->draw(x_sampler_path[level+1],x_path[level]);
      if (level > 0) {
        twolevel_step[level]->draw(x_sampler_path[level+1],
                                   x_sampler_path[level]);
      }
      qoi_sampler = qoi->evaluate(x_sampler_path[level]);
    }
    t_sampler[level]++;
    stats_sampler[level]->record_sample(qoi_sampler);
    n_samples_burnin[level]++;
    if (level > 0) {
      if ( (stats_sampler[level]->samples()>n_min_samples_corr) and
           (t_sampler[level]>=ceil(stats_sampler[level]->tau_int())) ) {
        t_sampler[level] = 0;
        // If we haven't reached the finest level, pass down to next level
        level--;
      } else {
        // Cycle back to coarsest level to pull down further samples
        level = n_level-1;
      } 
    } else {
      // Check that we've collected sufficient samples
      if (std::all_of(n_samples_burnin.begin(),
                      n_samples_burnin.end(),
                      [=](int n_samples) { return n_samples > n_burnin; })) {
        break;
      }
      level = n_level-1;
    }      
  } while(true);
  
  std::cout << "Burnin completed" << std::endl;

  double sum_s_ell2_C_ell_eff = 1E9;
  // Reset everything before actual run
  for (int level=0;level<n_level;++level) {
    stats_sampler[level]->reset();
    stats_qoi[level]->reset();
    t_sampler[level] = 0;
    t_indep[level] = 0;
    n_indep[level] = 0;
  }
  timer.reset();
  timer.start();
  level = n_level-1;
  do {
    // The QoI which is recorded on this level. This is Q_{L-1} on the
    // coarsest level and Y_{ell} = Q_{ell+1}-Q_{ell} on all othe levels
    double qoi_Y;
    // The QoI of the independent sampler, Q_{ell}
    double qoi_sampler;
    if (level == (n_level-1)) {
      /* Sample directly on coarsest level */
      coarse_sampler->draw(x_path[level]); // Path used for measuring QoI
      coarse_sampler->draw(x_sampler_path[level]); // Path used on next level
      qoi_sampler = qoi->evaluate(x_sampler_path[level]);
      qoi_Y = qoi->evaluate(x_path[level]);
    } else {
      /* 
       * On all other levels, sample by using the two level MCMC process
       * We know that the coarse level samples are decorrelated, since the
       * algorithm only ever proceeds to the next finer level if this is the
       * case.
       */
      twolevel_step[level]->draw(x_sampler_path[level+1],x_path[level]);
      if (level > 0) {
        // Generate a new sample which can be used on the next finer level
        twolevel_step[level]->draw(x_sampler_path[level+1],
                                   x_sampler_path[level]);
      }
      qoi_sampler = qoi->evaluate(x_sampler_path[level]);
      double qoi_fine = qoi->evaluate(x_path[level]);
      double qoi_coarse = qoi->evaluate(x_sampler_path[level+1]);
      qoi_Y = qoi_fine-qoi_coarse;
    }
    stats_sampler[level]->record_sample(qoi_sampler);
    t_sampler[level]++;
    stats_qoi[level]->record_sample(qoi_Y);
    if (level > 0) {
      if ( (stats_sampler[level]->samples() > n_min_samples_corr) and
           (t_sampler[level] >= ceil(stats_sampler[level]->tau_int())) ) {
        t_indep[level] = (n_indep[level]*t_indep[level]+t_sampler[level])/(1.0+n_indep[level]);
        n_indep[level]++;
        t_sampler[level] = 0;
        level--;
      } else {
        level = n_level-1;
      }
    } else {
      sum_s_ell2_C_ell_eff = 0;
      for (int ell=0;ell<n_level;++ell) {
        double V_ell = stats_qoi[ell]->variance();
        double C_ell_eff = cost_eff(ell);
        sum_s_ell2_C_ell_eff += sqrt(V_ell*C_ell_eff);
      }
      bool sufficient_stats = true;
      for (int ell=0;ell<n_level;++ell) {
        int n_samples = stats_qoi[ell]->samples();
        // Calculate variance and predict new number of target samples
        double V_ell = stats_qoi[ell]->variance();
        double tau_int = stats_qoi[ell]->tau_int();
        double C_ell_eff = cost_eff(ell);
        int n_target = ceil(two_epsilon_inv2*sum_s_ell2_C_ell_eff*sqrt(V_ell/C_ell_eff)*tau_int);
        sufficient_stats = sufficient_stats and (n_samples > std::max(int(n_min_samples_qoi),n_target));
      }
      if (sufficient_stats) break;
      level = n_level-1;
    }
    // Abort if we have collected sufficient statistics, this can
    // only be judged on the coarsest level
  } while (true);
  timer.stop();
}

/* Calculate effective cost on a particular level ell */
int MonteCarloMultiLevel::cost_eff(const int ell) const {
  int M_lat = fine_action->getM_lat();
  // Cost on current level
  int cost = (M_lat >> ell);
  // Add cost on coarser levels
  int T_k = ceil(t_indep[n_level-1]);
  int C_k = M_lat >> (n_level-1);
  for (int k=n_level-1;k>ell;--k) {
    cost += T_k*C_k;
    T_k *= ceil(t_indep[k]);
    C_k *= 2;
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
    qoi_error += stats_qoi[level]->error();
  }
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << std::endl;
  std::cout << " Q: Avg +/- Err = " << qoi_value << " +/- " << qoi_error << std::endl;
  std::cout << std::endl;
  std::cout << timer << std::endl;
  std::cout << std::endl;
}
