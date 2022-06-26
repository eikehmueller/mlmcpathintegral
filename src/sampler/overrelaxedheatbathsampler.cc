#include "overrelaxedheatbathsampler.hh"

/** @file overrelaxedheatbathsampler.cc
 * @brief Implementation of overrelaxedheatbathsampler.hh
 */

/** Draw next sample */
void OverrelaxedHeatBathSampler::draw(std::shared_ptr<SampleState> phi_state) {
  // Overrelaxation sweeps
  for (unsigned int sweep = 0; sweep < n_sweep_overrelax; ++sweep) {
    if (random_order) {
      shuffle(index_map.begin(), index_map.end(), engine);
    }
    for (auto it = index_map.begin(); it != index_map.end(); ++it) {
      action->overrelaxation_update(phi_state_cur, *it);
    }
  }
  // Heat bath sweeps
  for (unsigned int sweep = 0; sweep < n_sweep_heatbath; ++sweep) {
    if (random_order) {
      shuffle(index_map.begin(), index_map.end(), engine);
    }
    for (auto it = index_map.begin(); it != index_map.end(); ++it) {
      action->heatbath_update(phi_state_cur, *it);
    }
  }
  accept = true;
  n_total_samples++;
  n_accepted_samples++;
  phi_state->data = phi_state_cur->data;
}

/* Set current state */
void OverrelaxedHeatBathSampler::set_state(
    std::shared_ptr<SampleState> phi_state) {
  phi_state_cur->data = phi_state->data;
}
