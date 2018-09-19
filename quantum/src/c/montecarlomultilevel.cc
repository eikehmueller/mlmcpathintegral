#include "montecarlomultilevel.hh"

/** @file montecarlomultilevel.cc
 * @brief Implementation of montecarlomultilevel.hh 
 */

MonteCarloMultiLevel::MonteCarloMultiLevel(std::shared_ptr<Action> fine_action_,
                                           std::shared_ptr<QoI> qoi_,
                                           const GeneralParameters param_general,
                                           const HMCParameters param_hmc,
                                           const ClusterParameters param_cluster,
                                           const MultiLevelMCParameters param_multilevelmc) :
  MonteCarlo(1,1),
  fine_action(fine_action_),
  qoi(qoi_),
  n_level(param_multilevelmc.n_level()),
  epsilon(param_multilevelmc.epsilon()) {
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
    if (ell < n_level-1) {
      x_coarse_path.push_back(std::make_shared<Path>(M_lat/2,T_final));
    }
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
  int k_max = 20; // Record autocorrelations over a window of 20 steps
  for (int level=0;level<n_level;++level) {
    std::stringstream stats_label;
    stats_label << "Q[" << level << "]";
    stats.push_back(std::make_shared<Statistics>(stats_label.str(),k_max));
    std::stringstream stats_Y_label;
    stats_Y_label << "Y[" << level << "]";
    stats_Y.push_back(std::make_shared<Statistics>(stats_Y_label.str(),k_max));
  }
}    

void MonteCarloMultiLevel::evaluate() {
  /* Vector recording time since taking last sample from fine level chain
   * on a particular level. The samples are decorrelated of this time is
   * larger than the autocorrelation time of the fine level process.
   */
  std::vector<int> t(n_level,0);
  // Vector with target samples on each level. Record at least 100 samples.
  for (int level=0;level<n_level;++level) {
    stats[level]->reset();
    stats_Y[level]->reset();
  }
  double two_epsilon_inv2 = 2./(epsilon*epsilon);
  // Minimal number of samples to measure autocorrelation time
  int n_min_samples_autocorr = 100;
  // Minimal number of samples to measure QoI
  int n_min_samples_qoi = 100; 
  // Array which records whether we collected sufficient (uncorrelated)
  // samples on a particular level
  std::vector<bool> sufficient_stats(n_level,false);
  int level = n_level-1; // Current level
  double sum_V_ell_over_h_ell; // sum_{ell=0}^{L-1} V_{ell}/h_{ell}
  std::vector<double> V_ell(n_level,0.0); // Variance on a particular level
  do {
    double h_ell = action[level]->geta_lat(); // Lattice spacing on cur. level
    double qoi_fine; // The QoI for the fine level samples
    // The QoI which is recorded on this level. This is Q_{L-1} on the
    // coarsest level and Y_{ell} = Q_{ell+1}-Q_{ell} on all othe levels
    double qoi_Y; 
    if (level == (n_level-1)) {
      /* Sample directly on coarsest level */
      coarse_sampler->draw(x_path[level]);
      qoi_fine = qoi->evaluate(x_path[level]);
      qoi_Y = qoi_fine;
    } else {
      /* 
       * On all other levels, sample by using the two level MCMC process
       * We know that the coarse level samples are decorrelation, since the
       * algorithm only ever proceeds to the next finer level if this is the
       * case.
       */
      twolevel_step[level]->draw(x_path[level+1],x_path[level]);
      qoi_fine = qoi->evaluate(x_path[level]);
      double qoi_coarse = qoi->evaluate(x_path[level+1]);      
      qoi_Y = qoi_fine-qoi_coarse;
    }
    stats[level]->record_sample(qoi_fine);
    t[level] += 1;
    /* If the current path is sufficiently well decorrelated, use it to
     * calculate the estimator and down to the next-finer level.
     * Requiring the number of samples to be larger than n_min_samples ensures
     * that the measured autocorrelation time has been measured to 
     * sufficient accuracy.
     */
    if ( (stats[level]->samples() > n_min_samples_autocorr) and
         (t[level] > stats[level]->tau_int()) ) {
      stats_Y[level]->record_sample(qoi_Y);
      // Calculate variance and predict new number of target samples
      V_ell[level] = stats_Y[level]->variance();
      int n_samples = stats_Y[level]->samples();
      if (n_samples > n_min_samples_qoi) {
        int n_target = ceil(two_epsilon_inv2*sqrt(V_ell[level]*h_ell)*sum_V_ell_over_h_ell);
        sufficient_stats[level] = (n_samples > n_target);
      }
      if (level > 0) {
        // If we haven't reached the finest level, pass down to next level
        t[level] = 0;
        level--;
      } else {
        level = n_level-1;
        // Calculate the sum \sum_{ell=0}^{L-1} V_ell/h_ell 
        sum_V_ell_over_h_ell = 0.0;
        for (int i=0;i<n_level;++i) {
          double h_ell = action[i]->geta_lat();
          sum_V_ell_over_h_ell += V_ell[i]/h_ell;
        }
        if (std::all_of(sufficient_stats.begin(),
                        sufficient_stats.end(),
                        [](bool v) { return v; })) {
          break;
        }
      }      
    } else {
      // Cycle back to coarsest level
      level = n_level-1;
    }
    // Abort if we have collected sufficient statistics, this can
    // only be judged on the coarsest level
  } while (true);
  for (int level=0;level<n_level;++level) {
    std::cout << *stats[level];
    std::cout << *stats_Y[level];
    std::cout << "------------------------------------" << std::endl;
    std::cout << std::endl;
  }  
}