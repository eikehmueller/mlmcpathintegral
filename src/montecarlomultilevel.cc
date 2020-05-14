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
                                           const MultiLevelMCParameters param_multilevelmc,
                                           const HierarchicalParameters param_hierarchical) :
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
  timer("MultilevelMC") {
  // Check that Number of lattice points permits number of levels
  unsigned int M_lat = fine_action->getM_lat();
  if ( (M_lat>>n_level)<<n_level == M_lat) {
    mpi_parallel::cout << " MultilevelMC: M_lat = " << M_lat << " = 2^{" << n_level << "-1} * " << (M_lat>>(n_level-1)) << std::endl << std::endl;
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
    HierarchicalParameters param_hierarchical_tmp = HierarchicalParameters(param_hierarchical.n_level()-ell,
                               param_hierarchical.coarsesampler());
    std::shared_ptr<HierarchicalSampler> sampler_tmp = std::make_shared<HierarchicalSampler>(coarse_action_tmp,
                                              param_general,
                                              param_hmc,
                                              param_cluster,
                                              param_hierarchical_tmp);
    hierarchical_sampler.push_back(sampler_tmp);
  }
  std::shared_ptr<Action> coarse_action = action[n_level-1];
  // Construct paths on all levels
  double T_final = fine_action->getT_final();
  for (unsigned int ell=0;ell<n_level;++ell) {
    unsigned int M_lat = action[ell]->getM_lat();    
    x_path.push_back(std::make_shared<Path>(M_lat,T_final));
    x_coarse_path.push_back(std::make_shared<Path>(M_lat,T_final));
  }
  // Construct statistics on all levels
  for (unsigned int level=0;level<n_level;++level) {
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
  for (unsigned int level=0;level<n_level;++level) {
    stats_sampler[level]->hard_reset();
    stats_qoi[level]->hard_reset();
    t_sampler[level] = 0;
    n_indep[level] = 0;
    t_indep[level] = 0;
  }
  double two_epsilon_inv2 = 2./(epsilon*epsilon);

  // Burn in phase
  // Loop over all levels and create n_burnin samples
  for (int level=n_level-1;level>=0;level--) {
    double qoi_Y;
    for (unsigned int j=0;j<n_burnin;++j) {
      if (level==n_level-1) {
        draw_independent_sample(level,x_path[level]);
        qoi_Y = qoi->evaluate(x_path[level]);
      } else {
        draw_independent_sample(level+1,x_coarse_path[level+1]);
        twolevel_step[level]->draw(x_coarse_path[level+1],x_path[level]);
        double qoi_fine = qoi->evaluate(x_path[level]);
        double qoi_coarse = qoi->evaluate(x_coarse_path[level+1]);
        qoi_Y = qoi_fine-qoi_coarse;
      }
      stats_qoi[level]->record_sample(qoi_Y);
    }
  }
  mpi_parallel::cout << "Burnin completed" << std::endl;

  // Reset everything before actual run
  for (unsigned int level=0;level<n_level;++level) {
    stats_qoi[level]->reset();
    stats_sampler[level]->reset();
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
          double qoi_coarse = qoi->evaluate(x_coarse_path[level+1]);
          qoi_Y = qoi_fine-qoi_coarse;
        }
        stats_qoi[level]->record_sample(qoi_Y);
      }
    }
    // Now recompute target numbers on all levels \ell=0...L-1
    //
    // N_\ell^{eff} = 2/\epsilon^2*S*\sqrt{V_\ell/C_\ell^{eff}}
    //
    // with:
    //   * C_\ell^{eff} = cost for generating a sample on level \ell
    //   * V_\ell = variance on level \ell
    //   * S = \sum_{\ell=0}^{L-1} \sqrt{V_\ell*C_\ell^{eff}}
    //
    // Taking into account autocorrelations, the target number on level \ell
    // is then N_\ell = \ceil(\tau_\ell^{int}*N_\ell^{eff})
    sufficient_stats = true;
    // Compute sum S defined above
    double sum_s_ell2_C_ell_eff = 0;
    for (unsigned int ell=0;ell<n_level;++ell) {
      double V_ell = stats_qoi[ell]->variance();
      double C_ell_eff = cost_eff(ell);
      sum_s_ell2_C_ell_eff += sqrt(V_ell*C_ell_eff);
    }
    for (unsigned int ell=0;ell<n_level;++ell) {
      int n_samples = stats_qoi[ell]->samples();
      // Calculate variance and predict new number of target samples
      double V_ell = stats_qoi[ell]->variance();
      double tau_int = stats_qoi[ell]->tau_int();
      double C_ell_eff = cost_eff(ell);
      n_target[ell] = ceil(two_epsilon_inv2*sum_s_ell2_C_ell_eff*sqrt(V_ell/C_ell_eff)*tau_int);
      sufficient_stats = sufficient_stats and (n_samples >= n_target[ell]);
    }
  } while (not sufficient_stats); 
  timer.stop();
}

/* Draw an independent sample on a particular level */
void MonteCarloMultiLevel::draw_independent_sample(const unsigned int ell,
                                                   std::shared_ptr<Path> x_path) {
  while (true) {
    hierarchical_sampler[ell-1]->draw(x_path);
    double qoi_sampler = qoi->evaluate(x_path);
    stats_sampler[ell]->record_sample(qoi_sampler);
    t_sampler[ell]++;
    if (t_sampler[ell] >= ceil(stats_sampler[ell]->tau_int())) {
      t_indep[ell] = (n_indep[ell]*t_indep[ell]+t_sampler[ell])/(1.0+n_indep[ell]);
      n_indep[ell]++;
      t_sampler[ell] = 0;
      break;
    }
  }
}

/* Calculate effective cost on a particular level ell */
double MonteCarloMultiLevel::cost_eff(const int ell) const {
  // Cost on current level
  double cost;
  if (ell==n_level-1) {
    cost = hierarchical_sampler[ell-1]->cost_per_sample()*t_indep[ell];
  } else {
    double cost_twolevel = twolevel_step[ell]->cost_per_sample();
    double cost_coarse = hierarchical_sampler[ell]->cost_per_sample();
    cost = cost_twolevel + t_indep[ell+1]*cost_coarse;
  }
  return ceil(stats_qoi[ell]->tau_int())*cost;
}
  
/* Show detailed statistics on all levels */
void MonteCarloMultiLevel::show_detailed_statistics() {
  mpi_parallel::cout << "=== Statistics of coarse level samplers ===" << std::endl;
  for (unsigned int level=1;level<n_level;++level) {
    mpi_parallel::cout << "level = " << level << std::endl;
    mpi_parallel::cout << *stats_sampler[level];
    mpi_parallel::cout << " spacing between samples " << t_indep[level] << std::endl;
    mpi_parallel::cout << " hierarchical sampler stats " << std::endl;
    hierarchical_sampler[level-1]->show_stats();
    mpi_parallel::cout << "------------------------------------" << std::endl;    
  }
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << "=== Statistics of QoI ===" << std::endl;
  for (unsigned int level=0;level<n_level;++level) {
    mpi_parallel::cout << "level = " << level << std::endl;
    mpi_parallel::cout << *stats_qoi[level];
    mpi_parallel::cout << " target number of samples = " << n_target[level] << std::endl;
    mpi_parallel::cout << " cost [per indep. sample]              " << cost_eff(level) << " mu s" << std::endl;
    if (level < n_level-1) {
      mpi_parallel::cout << " cost [per sample in two-level step] " << twolevel_step[level]->cost_per_sample() << " mu s" << std::endl;
    }
    mpi_parallel::cout << "------------------------------------" << std::endl;
  }
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << "=== Breakdown of cost ===" << std::endl;
  std::vector<double> cost_level;
  double cost_total=0.0;
  for (unsigned int level=0;level<n_level;++level) {
    double cost_per_level = 1.E-6*cost_eff(level)/ceil(stats_qoi[level]->tau_int())*stats_qoi[level]->samples();
    cost_level.push_back(cost_per_level);
    cost_total+=cost_per_level;
  }
  for (unsigned int level=0;level<n_level;++level) {
    mpi_parallel::cout << " level " << level << " : " << cost_level[level] << " s [ " << 100*cost_level[level]/cost_total << " ] %" << std::endl;
  }
  mpi_parallel::cout << " total = " << cost_total << " s" << std::endl;
  mpi_parallel::cout << std::endl;
}

/* Calculate numerical value */
double MonteCarloMultiLevel::numerical_result() const {
  double qoi_value = 0.0;
  for (unsigned int level=0;level<n_level;++level) {
    qoi_value += stats_qoi[level]->average();
  }
  return qoi_value;
}

/* Calculate statistical error */
double MonteCarloMultiLevel::statistical_error() const {
  double qoi_error = 0.0;
  for (unsigned int level=0;level<n_level;++level) {
    double error = stats_qoi[level]->error();
    qoi_error += error*error;
  }
  return sqrt(qoi_error);
}

/* Show summary statistics and timing */
void MonteCarloMultiLevel::show_statistics() {
  mpi_parallel::cout << std::setprecision(6) << std::fixed;
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << " Q: Avg +/- Err = " << numerical_result() << " +/- " << statistical_error() << std::endl;
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << timer << std::endl;
  mpi_parallel::cout << std::endl;
}
