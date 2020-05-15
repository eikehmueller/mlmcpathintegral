#include "montecarlosinglelevel.hh"

/** @file montecarlosinglelevel.cc
 * @brief Implementation of montecarlosinglelevel.hh 
 */

/* Constructor */
MonteCarloSingleLevel::MonteCarloSingleLevel(std::shared_ptr<Action> action_,
                                             std::shared_ptr<QoI> qoi_,
                                             const GeneralParameters param_general,
                                             const StatisticsParameters param_stats,
                                             const HMCParameters param_hmc,
                                             const ClusterParameters param_cluster,
                                             const SingleLevelMCParameters param_singlelevelmc,
                                             const HierarchicalParameters param_hierarchical) :
  MonteCarlo(param_singlelevelmc.n_burnin()),
  action(action_), 
  qoi(qoi_),
  n_autocorr_window(param_stats.n_autocorr_window()),
  n_min_samples_qoi(param_stats.n_min_samples_qoi()),
  n_samples(param_singlelevelmc.n_samples()),
  epsilon(param_singlelevelmc.epsilon()),
  timer("SinglevelMC") {
  if (param_singlelevelmc.sampler() == SamplerHMC) {
    sampler = std::make_shared<HMCSampler>(action,
                                           param_hmc);
  } else if (param_singlelevelmc.sampler() == SamplerCluster) {
    if (param_general.action() != ActionRotor) {
      mpi_parallel::cerr << " ERROR: can only use cluster sampler for QM rotor action." << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    sampler = std::make_shared<ClusterSampler>(std::dynamic_pointer_cast<ClusterAction>(action),
                                               param_cluster.n_burnin(),
                                               param_cluster.n_updates());
  } else if (param_singlelevelmc.sampler() == SamplerExact) {
    if (param_general.action() != ActionHarmonicOscillator) {
      mpi_parallel::cerr << " ERROR: can only sample exactly from harmonic oscillator action." << std::endl;
      mpi_exit(EXIT_FAILURE);
    }
    sampler = std::dynamic_pointer_cast<Sampler>(action);
  } else if (param_singlelevelmc.sampler() == SamplerHierarchical) {
    sampler = std::make_shared<HierarchicalSampler>(action,
                                                    qoi,
                                                    param_general,
                                                    param_stats,
                                                    param_hmc,
                                                    param_cluster,
                                                    param_hierarchical);
  } else {
    mpi_parallel::cerr << " ERROR: Unknown sampler." << std::endl;
    mpi_exit(EXIT_FAILURE);
  }
  stats_Q = std::make_shared<Statistics>("Q",n_autocorr_window);
  }

/** Calculate Monte Carlo estimate with single level method */
void MonteCarloSingleLevel::evaluate() {
  std::ofstream qoi_file; // File to write QoI to
  std::shared_ptr<Path> x_path =
    std::make_shared<Path>(action->getM_lat(),
                           action->getT_final());
  stats_Q->hard_reset();
  for (unsigned int i=0;i<n_burnin;++i) {
    sampler->draw(x_path);
    double qoi_Q = qoi->evaluate(x_path);
    stats_Q->record_sample(qoi_Q);
  }

  mpi_parallel::cout << "Burnin completed" << std::endl;
  
  double two_epsilon_inv2 = 2./(epsilon*epsilon);
  stats_Q->reset();
  bool sufficient_stats = false;
  unsigned int n_target;
  unsigned int n_local_target;
  if (n_samples > 0) {
    n_target = n_samples;
  } else {
    n_target = n_min_samples_qoi;
  }
#ifdef LOG_QOI
  if (mpi_comm_size() > 1) {
    mpi_parallel::cerr << "ERROR: can only log QoI in sequential runs." << std::endl;
    mpi_exit(EXIT_FAILURE);
  }
  qoi_file.open("qoi.dat");
#endif // LOG_QOI
  n_local_target = distribute_n(n_target);
  timer.reset();
  timer.start();
  do {
    unsigned int k_start = stats_Q->local_samples();
    for (unsigned int k=k_start;k<n_local_target;++k) {
      sampler->draw(x_path);
#ifdef SAVE_PATHS
      if ( (SAVE_FIRST_PATH<=k) and (k<=SAVE_LAST_PATH) and (mpi_master())) {
        std::stringstream filename;      
        filename << "path_";
        filename << std::setw(8) << std::setfill('0');
        filename << k << ".dat";
        x_path->save_to_disk(filename.str());
      }
#endif // SAVE_PATHS
      // Quantity of interest
      double qoi_Q = qoi->evaluate(x_path);
      stats_Q->record_sample(qoi_Q);
#ifdef LOG_QOI
      qoi_file.write( reinterpret_cast<char*>( &qoi_Q ), sizeof(double));
#endif // LOG_QOI
    }
    if (n_samples > 0) {
      n_target = n_samples;
    } else {
      n_target = ceil(stats_Q->tau_int()*two_epsilon_inv2*stats_Q->variance());
    }
    n_local_target = distribute_n(n_target);
    sufficient_stats = mpi_allreduce_and(stats_Q->local_samples() >= n_local_target);
    // If the target number of samples is given explicitly, generate exactly
    // the number of requested samples
  } while (not sufficient_stats);
  timer.stop();
#ifdef LOG_QOI
  qoi_file.close();
#endif // LOG_QOI
}

/* Print out statistics */
void MonteCarloSingleLevel::show_statistics() {
  mpi_parallel::cout << *stats_Q;
  mpi_parallel::cout << std::endl;
  mpi_parallel::cout << timer << std::endl;
  mpi_parallel::cout << std::endl;
}
