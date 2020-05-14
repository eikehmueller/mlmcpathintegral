#include "hierarchicalsampler.hh"

/** @file hierarchicalsampler.cc
 * @brief Implementation of hierarchicalsampler.hh
 */

/* Construct new instance */
HierarchicalSampler::HierarchicalSampler(const std::shared_ptr<Action> fine_action,
                                         const GeneralParameters param_general,
                                         const HMCParameters param_hmc,
                                         const ClusterParameters param_cluster,
                                         const HierarchicalParameters param_hierarchical) :
  Sampler(),
  n_level(param_hierarchical.n_level()),
  cost_per_sample_(0.0){
  
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
                                                  param_hmc.nt(),
                                                  param_hmc.dt(),
                                                  param_hmc.n_burnin());
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
void HierarchicalSampler::draw(std::shared_ptr<Path> x_path) {
  accept = true;
  for (int ell=n_level-1;ell>=0;--ell) {
    if (ell == (n_level-1)) {
      /* Sample directly on coarsest level */
      for (int i=0;i<1;++i) {
        coarse_sampler->draw(x_sampler_path[ell]);
      }
    } else {
      /* 
       * On all other levels, sample by using the two level MCMC process
       * We know that the coarse level samples are decorrelated, since the
       * algorithm only ever proceeds to the next finer level if this is the
       * case.
       */
      twolevel_step[ell]->draw(x_sampler_path[ell+1],
                               x_sampler_path[ell]);
    }
  }
  if (n_level > 1) {
    accept = twolevel_step[0]->accepted();
  } else {
    accept = coarse_sampler->accepted();
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  // Copy path
  if (copy_if_rejected or accept) {
    std::copy(x_sampler_path[0]->data,
              x_sampler_path[0]->data+x_sampler_path[0]->M_lat,
              x_path->data);
  }
}

/* Show statistics on all levels */
void HierarchicalSampler::show_stats() {
  mpi_parallel::cout << std::setprecision(3) << std::fixed;
  mpi_parallel::cout << "   cost per sample = " << cost_per_sample() << " mu s" << std::endl;
    mpi_parallel::cout << "   acceptance/rejection probability   p      1-p" << std::endl;
    for (unsigned int ell=0;ell<n_level;++ell) {
        double p_acc;
        if (ell == n_level-1) {
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
        sstream << "   level " << ell << " " << level_str <<" : ";
        sstream << "              " << p_acc << "  " << 1.-p_acc;
        mpi_parallel::cout << sstream.str() << std::endl;
    }
}
