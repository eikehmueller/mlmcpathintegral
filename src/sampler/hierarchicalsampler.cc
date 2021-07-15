#include "hierarchicalsampler.hh"

/** @file hierarchicalsampler.cc
 * @brief Implementation of hierarchicalsampler.hh
 */

/* Construct new instance */
HierarchicalSampler::HierarchicalSampler(const std::shared_ptr<Action> fine_action,
        const std::shared_ptr<SamplerFactory> coarse_sampler_factory,
        const std::shared_ptr<ConditionedFineActionFactory> conditioned_fine_action_factory,
        const HierarchicalParameters param_hierarchical) :
    Sampler(),
    n_level(param_hierarchical.n_max_level()-fine_action->get_coarsening_level()),
    cost_per_sample_(0.0) {
    mpi_parallel::cout << "Hierarchical sampler: " << std::endl;
    action.push_back(fine_action);
    // Construct action and two-level MCMC step on all levels
    for (unsigned int ell=0; ell<n_level-1; ++ell) {
        std::shared_ptr<Action> action_tmp = action[ell];
        std::shared_ptr<Action> coarse_action_tmp = action[ell]->coarse_action();
        action.push_back(coarse_action_tmp);
        std::shared_ptr<ConditionedFineAction> conditioned_fine_action;
        conditioned_fine_action = conditioned_fine_action_factory->get(action_tmp);
        std::shared_ptr<TwoLevelMetropolisStep> twolevel_step_tmp =
            std::make_shared<TwoLevelMetropolisStep>(coarse_action_tmp,
                    action_tmp,
                    conditioned_fine_action);
        twolevel_step.push_back(twolevel_step_tmp);
    }
    for (unsigned int ell=0; ell<n_level; ++ell) {
        mpi_parallel::cout << "  action on level " << ell << " : " << action[ell]->info_string() << std::endl;
    }
    // Action on coarsest level
    std::shared_ptr<Action> coarse_action = action[n_level-1];
    for (unsigned int ell=0; ell<n_level; ++ell) {
        phi_sampler_state.push_back(std::make_shared<SampleState>(action[ell]->sample_size()));
    }
    // Construct sampler on coarsest level
    coarse_sampler = coarse_sampler_factory->get(coarse_action);
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

void HierarchicalSampler::draw(std::shared_ptr<SampleState> phi_state) {
    accept = true;
    for (int ell=1; ell<n_level; ++ell) {
        action[ell]->copy_from_fine(phi_sampler_state[ell-1],phi_sampler_state[ell]);
    }
    for (int ell=n_level-1; ell>=0; --ell) {
        if (ell == (n_level-1)) {
            /* Sample directly on coarsest level */
            coarse_sampler->set_state(phi_sampler_state[ell]);
            coarse_sampler->draw(phi_sampler_state[ell]);
            accept = accept and coarse_sampler->accepted();
        } else {
            twolevel_step[ell]->set_state(phi_sampler_state[ell]);
            twolevel_step[ell]->draw(phi_sampler_state[ell+1],
                                     phi_sampler_state[ell]);
            accept = accept and twolevel_step[ell]->accepted();
        }
        if (not accept) break;
    }
    n_total_samples++;
    n_accepted_samples += (int) accept;
    if (accept or copy_if_rejected) {
        phi_state->data = phi_sampler_state[0]->data;
    }
}

/* Set current state */
void HierarchicalSampler::set_state(std::shared_ptr<SampleState> phi_state) {
    phi_sampler_state[0]->data = phi_state->data;
}

/* Show statistics on all levels */
void HierarchicalSampler::show_stats() {
    mpi_parallel::cout << std::setprecision(4) << std::fixed;
    mpi_parallel::cout << "  cost per sample = " << cost_per_sample() << " mu s" << std::endl << std::endl;
    mpi_parallel::cout << "  acceptance rate = " << p_accept() << std::endl;
    for (unsigned int ell=0; ell<n_level; ++ell) {
        double p_acc;
        if (ell==n_level-1) {
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
        sstream << "  level " << ell << " " << level_str <<" : ";
        sstream << " p = " << p_acc << std::endl;
        mpi_parallel::cout << sstream.str();
    }
}
