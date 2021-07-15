#include "multilevelsampler.hh"

/** @file multilevelsampler.cc
 * @brief Implementation of multilevelsampler.hh
 */

/* Construct new instance */
MultilevelSampler::MultilevelSampler(const std::shared_ptr<Action> fine_action,
                                     const std::shared_ptr<QoIFactory> qoi_factory_,
                                     const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
                                     const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
                                     const StatisticsParameters param_stats,
                                     const HierarchicalParameters param_hierarchical) :
    Sampler(),
    n_level(param_hierarchical.n_max_level()-fine_action->get_coarsening_level()),
    t_indep(n_level,0.0),
    n_indep(n_level,0),
    t_sampler(n_level,0),
    n_autocorr_window(param_stats.n_autocorr_window()),
    cost_per_sample_(0.0) {
    mpi_parallel::cout << " Multilevel sampler: ";
    action.push_back(fine_action);
    // Construct action and two-level MCMC step on all levels
    for (unsigned int ell=0; ell<n_level-1; ++ell) {
        std::shared_ptr<Action> action_tmp = action[ell];
        std::shared_ptr<Action> coarse_action_tmp = action[ell]->coarse_action();
        action.push_back(coarse_action_tmp);
        std::shared_ptr<ConditionedFineAction> conditioned_fine_action=conditioned_fine_action_factory->get(action_tmp);
        std::shared_ptr<TwoLevelMetropolisStep> twolevel_step_tmp =
            std::make_shared<TwoLevelMetropolisStep>(coarse_action_tmp,
                    action_tmp,
                    conditioned_fine_action);
        twolevel_step.push_back(twolevel_step_tmp);
    }
    for (unsigned int ell=0; ell<n_level; ++ell) {
        mpi_parallel::cout << "  action on level " << ell << " : " << action[ell]->info_string() << std::endl;
    }
    // Construct QoI sampler factory
    for (unsigned int ell=0; ell<n_level; ++ell) {
        qoi.push_back(qoi_factory_->get(action[ell]));
    }
    // Action on coarsest level
    std::shared_ptr<Action> coarse_action = action[n_level-1];
    for (unsigned int ell=0; ell<n_level; ++ell) {
        phi_sampler_state.push_back(std::make_shared<SampleState>(action[ell]->sample_size()));
    }
    // Construct sampler on coarsest level
    coarse_sampler = coarse_sampler_factory->get(coarse_action);
    // Statistics on all levels
    for (unsigned int level=0; level<n_level; ++level) {
        std::stringstream stats_sampler_label;
        stats_sampler_label << "   Q_{sampler}[" << level << "]";
        stats_sampler.push_back(std::make_shared<Statistics>(stats_sampler_label.str(),n_autocorr_window));
    }
    std::shared_ptr<SampleState> meas_state=std::make_shared<SampleState>(fine_action->sample_size());
    Timer timer_meas;
    unsigned int n_meas = 10000;
    timer_meas.start();
    for (unsigned int k=0; k<n_meas; ++k) {
        draw(meas_state);
    }
    timer_meas.stop();
    cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
}

/* Draw next sample */
void MultilevelSampler::draw(std::shared_ptr<SampleState> phi_state) {
    accept = true;
    int level = n_level-1;
    do {
        if (level == (n_level-1)) {
            /* Sample directly on coarsest level */
            coarse_sampler->draw(phi_sampler_state[level]);
        } else {
            /*
             * On all other levels, sample by using the two level MCMC process
             * We know that the coarse level samples are decorrelated, since the
             * algorithm only ever proceeds to the next finer level if this is the
             * case.
             */
            twolevel_step[level]->draw(phi_sampler_state[level+1],
                                       phi_sampler_state[level]);
        }
        // The QoI of the independent sampler, Q_{ell}
        double qoi_sampler = qoi[level]->evaluate(phi_sampler_state[level]);
        stats_sampler[level]->record_sample(qoi_sampler);
        t_sampler[level]++;
        if (t_sampler[level] >= ceil(stats_sampler[level]->tau_int())) {
            // t_sampler[level] is the number of phi_sampler_state samples
            // generated on this level since the last independent sample was used
            t_indep[level] = (n_indep[level]*t_indep[level]+t_sampler[level])/(1.0+n_indep[level]);
            // t_indep[level] is the average number of samples between indpendent
            // samples
            n_indep[level]++;
            // n_indep is the number of independent samples of phi_sampler_state on
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
    // Copy state
    phi_state->data = phi_sampler_state[0]->data;
}

/* Set current state */
void MultilevelSampler::set_state(std::shared_ptr<SampleState> phi_state) {
    phi_sampler_state[0]->data = phi_state->data;
}

/* Show statistics on all levels */
void MultilevelSampler::show_stats() {
    mpi_parallel::cout << std::setprecision(3) << std::fixed;
    mpi_parallel::cout << "  cost per sample = " << cost_per_sample() << " mu s" << std::endl << std::endl;
    for (unsigned int ell=0; ell<n_level; ++ell) {
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
