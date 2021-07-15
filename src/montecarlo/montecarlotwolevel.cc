#include "montecarlotwolevel.hh"

/** @file montecarlotwolevel.cc
 * @brief Implementation of montecarlotwolevel.hh
 */
/* Constructor */
MonteCarloTwoLevel::MonteCarloTwoLevel(const std::shared_ptr<Action> fine_action_,
                                       const std::shared_ptr<QoIFactory> qoi_factory_,
                                       const std::shared_ptr<SamplerFactory> sampler_factory,
                                       const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                       const StatisticsParameters param_stats,
                                       const TwoLevelMCParameters param_twolevelmc) :
    MonteCarlo(param_twolevelmc.n_burnin()),
    n_samples(param_twolevelmc.n_samples()),
    fine_action(fine_action_),
    qoi_fine(qoi_factory_->get(fine_action_)),
    qoi_coarse(qoi_factory_->get(fine_action_->coarse_action())),
    t_indep(0.0),
    n_indep(0),
    t_sampler(0),
    stats_fine("QoI[fine]",param_twolevelmc.n_fine_autocorr_window()),
    stats_coarse("QoI[coarse]",param_twolevelmc.n_coarse_autocorr_window()),
    stats_diff("delta QoI",param_twolevelmc.n_delta_autocorr_window()),
    stats_coarse_sampler("QoI[coarsesampler]",param_stats.n_autocorr_window()) {
        mpi_parallel::cout << "Twolevel Monte Carlo:" << std::endl;
        coarse_action = fine_action->coarse_action();
        mpi_parallel::cout << "  fine action   : " << fine_action->info_string() << std::endl;
        mpi_parallel::cout << "  coarse action : " << coarse_action->info_string() << std::endl;
        coarse_sampler = sampler_factory->get(coarse_action);
        conditioned_fine_action = conditioned_fine_action_factory->get(fine_action);
        twolevel_step = std::make_shared<TwoLevelMetropolisStep>(coarse_action,
                                                                 fine_action,
                                                                 conditioned_fine_action);
}

/** Calculate mean and variance of difference in QoI */
void MonteCarloTwoLevel::evaluate_difference() {
    std::shared_ptr<SampleState> phi_state =
        std::make_shared<SampleState>(fine_action->sample_size());
    std::shared_ptr<SampleState> phi_coarse_state =
        std::make_shared<SampleState>(coarse_action->sample_size());
    stats_coarse.hard_reset();
    stats_coarse_sampler.hard_reset();
    stats_fine.hard_reset();
    stats_diff.hard_reset();

    // Burn-in phase
    for (unsigned int k=0; k<n_burnin; ++k) {
        draw_coarse_sample(phi_coarse_state);
        twolevel_step->draw(phi_coarse_state,phi_state);
        double qoi_fine_value = qoi_fine->evaluate(phi_state);
        double qoi_coarse_value = qoi_coarse->evaluate(phi_coarse_state);
        stats_fine.record_sample(qoi_fine_value);
        stats_coarse.record_sample(qoi_coarse_value);
        stats_diff.record_sample(qoi_fine_value-qoi_coarse_value);
    }
    mpi_parallel::cout << "Burnin completed" << std::endl;
    stats_coarse_sampler.reset();

    // Work out number of samples
    int n_proc = mpi_comm_size();
    unsigned int n_local_samples = (unsigned int) ceil(n_samples/(1.0*n_proc));
    
    // Sampling phase
    // Do a hard reset since we are interested in the variance
    stats_coarse.hard_reset();
    stats_fine.hard_reset();
    stats_diff.hard_reset();
    for (unsigned int k=0; k<n_local_samples; ++k) {
        draw_coarse_sample(phi_coarse_state);
        twolevel_step->draw(phi_coarse_state,phi_state);
        double qoi_fine_value = qoi_fine->evaluate(phi_state);
        double qoi_coarse_value = qoi_coarse->evaluate(phi_coarse_state);
        stats_fine.record_sample(qoi_fine_value);
        stats_coarse.record_sample(qoi_coarse_value);
        stats_diff.record_sample(qoi_fine_value-qoi_coarse_value);
    }
}

/* Draw (independent) coarse level sample */
void MonteCarloTwoLevel::draw_coarse_sample(std::shared_ptr<SampleState> phi_state) {
    double two_tau_int = fmin(100,ceil(2.*stats_coarse_sampler.tau_int()));
    while (t_sampler < two_tau_int) {
        coarse_sampler->draw(phi_state);
        double qoi_sampler = qoi_coarse->evaluate(phi_state);
        stats_coarse_sampler.record_sample(qoi_sampler);
        t_sampler++;
    }
    t_indep = (n_indep*t_indep+t_sampler)/(1.0+n_indep);
    n_indep++;
    t_sampler = 0; // Reset number of independent samples
}

/* Print out statistics */
void MonteCarloTwoLevel::show_statistics() const {
    mpi_parallel::cout << stats_fine << std::endl;
    mpi_parallel::cout << stats_coarse << std::endl;
    mpi_parallel::cout << stats_diff << std::endl;
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << "=== Coarse level sampler statistics === " << std::endl;
    mpi_parallel::cout << stats_coarse_sampler << std::endl;
    coarse_sampler->show_stats();
    mpi_parallel::cout << std::endl;
    mpi_parallel::cout << "=== Two level sampler statistics === " << std::endl;
    twolevel_step->show_stats();
}
