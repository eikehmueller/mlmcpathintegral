#include "overrelaxedheatbathsampler.hh"

/** @file overrelaxedheatbathsampler.cc
 * @brief Implementation of overrelaxedheatbathsampler.hh
 */

/** Draw next sample */
void OverrelaxedHeatBathSampler::draw(std::shared_ptr<SampleState> phi_state) {
    // Overrelaxation sweeps
    for (unsigned int sweep=0;sweep<n_sweep_overrelax;++sweep) {
        if (random_order) {
            shuffle (index_map.begin(), index_map.end(), engine);
        }
        for (unsigned int i_dof=0; i_dof<action->sample_size(); ++i_dof) {
            action->overrelaxation_update(phi_state_cur,index_map[i_dof]);
        }
    }
    // Heat bath sweeps
    for (unsigned int sweep=0;sweep<n_sweep_heatbath;++sweep) {
        if (random_order) {
            shuffle (index_map.begin(), index_map.end(), engine);
        }
        for (unsigned int i_dof=0; i_dof<action->sample_size(); ++i_dof) {
            action->heatbath_update(phi_state_cur,index_map[i_dof]);
        }
    }
    accept = true;
    n_total_samples++;
    n_accepted_samples++;
    phi_state->data = phi_state_cur->data;
}

/* Set current state */
void OverrelaxedHeatBathSampler::set_state(std::shared_ptr<SampleState> phi_state) {
    phi_state_cur->data = phi_state->data;
}

