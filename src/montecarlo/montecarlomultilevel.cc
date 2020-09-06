#include "montecarlomultilevel.hh"

/** @file montecarlomultilevel.cc
 * @brief Implementation of montecarlomultilevel.hh 
 */

MonteCarloMultiLevel::MonteCarloMultiLevel(std::shared_ptr<Action> fine_action_,
                                           std::shared_ptr<QoI> qoi_,
                                           std::shared_ptr<SamplerFactory> sampler_factory,                                           
                                           std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                           const StatisticsParameters param_stats,
                                           const MultiLevelMCParameters param_multilevelmc) :
  MonteCarlo(param_multilevelmc.n_burnin()),
  fine_action(fine_action_),
  qoi(qoi_),
  n_level(param_multilevelmc.n_level()),
  epsilon(param_multilevelmc.epsilon()),
  n_target(param_multilevelmc.n_level(),0),
  n_autocorr_window(param_stats.n_autocorr_window()),
  n_min_samples_qoi(param_stats.n_min_samples_qoi()),
  timer("MultilevelMC"),
  t_indep(param_multilevelmc.n_level(),0.0),
  n_indep(param_multilevelmc.n_level(),0),
  t_sampler(param_multilevelmc.n_level(),0) {
  mpi_parallel::cout << " MultilevelMC: ";
  fine_action->check_coarsening_is_permitted(n_level);
  action.push_back(fine_action);
  // Construct action and two-level MCMC step on all levels
  for (unsigned int ell=0;ell<n_level-1;++ell) {
    std::shared_ptr<Action> action_tmp = action[ell];
    std::shared_ptr<Action> coarse_action_tmp = action[ell]->coarse_action();
    action.push_back(coarse_action_tmp);
    std::shared_ptr<ConditionedFineAction> conditioned_fine_action=conditioned_fine_action_factory->get(action_tmp);
    std::shared_ptr<TwoLevelMetropolisStep> twolevel_step_tmp =
      std::make_shared<TwoLevelMetropolisStep>(coarse_action_tmp,
                                               action_tmp,
                                               conditioned_fine_action);
    twolevel_step.push_back(twolevel_step_tmp);
    std::shared_ptr<Sampler> sampler_tmp=sampler_factory->get(coarse_action_tmp);
    sub_sample_coarse = true;
    coarse_sampler.push_back(sampler_tmp);
    std::stringstream stats_sampler_label;
    stats_sampler_label << "   Q_{sampler}[" << ell << "]";
    stats_coarse_sampler.push_back(std::make_shared<Statistics>(stats_sampler_label.str(),
                                                                n_autocorr_window));
  }
  std::shared_ptr<Action> coarse_action = action[n_level-1];
  // Construct samples on all levels
  for (unsigned int ell=0;ell<n_level;++ell) {
    x_path.push_back(std::make_shared<SampleState>(action[ell]->sample_size()));
    x_coarse_path.push_back(std::make_shared<SampleState>(action[ell]->sample_size()));
  }
  // Construct statistics on all levels
  for (unsigned int level=0;level<n_level;++level) {
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
    stats_qoi[level]->hard_reset();
  }
  double two_epsilon_inv2 = 2./(epsilon*epsilon);

  // Burn in phase
  // Loop over all levels and create n_burnin samples
  for (int level=n_level-1;level>=0;level--) {
    double qoi_Y;
    for (unsigned int j=0;j<n_burnin;++j) {
      if (level==n_level-1) {
        coarse_sampler[level-1]->draw(x_path[level]);
        qoi_Y = qoi->evaluate(x_path[level]);
      } else {
        coarse_sampler[level]->draw(x_coarse_path[level+1]);
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
    if (level < n_level-1) {
      stats_coarse_sampler[level]->reset();
    }
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
          draw_coarse_sample(level,x_path[level]);
          qoi_Y = qoi->evaluate(x_path[level]);
        } else {
          draw_coarse_sample(level+1,x_coarse_path[level+1]);
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

/* Draw (independent) coarse level sample */
void MonteCarloMultiLevel::draw_coarse_sample(const unsigned int level,
                                              std::shared_ptr<SampleState> x_path) {
  // sub-sample if we are using the hierarchical sampler
  if (sub_sample_coarse) {
    while (t_sampler[level-1] < ceil(2.*stats_coarse_sampler[level-1]->tau_int())) {
      coarse_sampler[level-1]->draw(x_path);
      double qoi_sampler = qoi->evaluate(x_path);
      stats_coarse_sampler[level-1]->record_sample(qoi_sampler);
      t_sampler[level-1]++;
    }
  } else {
    coarse_sampler[level-1]->draw(x_path);
    t_sampler[level-1] = 1;
  }
  t_indep[level-1] = (n_indep[level-1]*t_indep[level-1]+t_sampler[level-1])/(1.0+n_indep[level-1]);
  n_indep[level-1]++;
  t_sampler[level-1] = 0; // Reset number of independent samples
}

/* Calculate effective cost on a particular level ell */
double MonteCarloMultiLevel::cost_eff(const int ell) const {
  // Cost on current level
  double cost;
  if (ell==n_level-1) {
    cost = t_indep[ell-1]*coarse_sampler[ell-1]->cost_per_sample();
  } else {
    double cost_twolevel = twolevel_step[ell]->cost_per_sample();
    double cost_coarse = t_indep[ell]*coarse_sampler[ell]->cost_per_sample();
    cost = cost_twolevel + cost_coarse;
  }
  return ceil(stats_qoi[ell]->tau_int())*cost;
}
  
/* Show detailed statistics on all levels */
void MonteCarloMultiLevel::show_detailed_statistics() {
  mpi_parallel::cout << "=== Statistics of coarse level samplers ===" << std::endl;
  for (unsigned int level=1;level<n_level;++level) {
    mpi_parallel::cout << "level = " << level << std::endl;
    mpi_parallel::cout << " coarse level sampler stats " << std::endl;
    coarse_sampler[level-1]->show_stats();
    mpi_parallel::cout << " sub-sampling rate = " << t_indep[level-1]<< std::endl;
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
