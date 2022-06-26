#include "twolevelmetropolisstep.hh"
/** @file twolevelmetropolissampler.cc
 * @brief Implementation of twolevelmetropolissampler.hh
 */
TwoLevelMetropolisStep::TwoLevelMetropolisStep(
    const std::shared_ptr<Action> coarse_action_,
    const std::shared_ptr<Action> fine_action_,
    const std::shared_ptr<ConditionedFineAction> conditioned_fine_action_)
    : MCMCStep(), coarse_action(coarse_action_), fine_action(fine_action_),
      conditioned_fine_action(conditioned_fine_action_), cost_per_sample_(0.0) {
  theta_fine = std::make_shared<SampleState>(fine_action->sample_size());
  theta_fine_C = std::make_shared<SampleState>(coarse_action->sample_size());
  theta_prime = std::make_shared<SampleState>(fine_action->sample_size());
  engine.seed(89216491);
  std::shared_ptr<SampleState> phi_fine =
      std::make_shared<SampleState>(fine_action->sample_size());
  std::shared_ptr<SampleState> phi_coarse =
      std::make_shared<SampleState>(coarse_action->sample_size());

  fine_action_theta_fine = fine_action->evaluate(theta_fine);
  conditioned_fine_action_theta_fine =
      conditioned_fine_action->evaluate(theta_fine);
  Timer timer_meas;
  unsigned int n_meas = 10000;
  timer_meas.start();
  for (unsigned int k = 0; k < n_meas; ++k) {
    draw(phi_coarse, phi_fine);
  }
  timer_meas.stop();
  cost_per_sample_ = 1.E6 * timer_meas.elapsed() / n_meas;
  reset_stats();
}

/** Draw new sample pair */
void TwoLevelMetropolisStep::draw(
    const std::shared_ptr<SampleState> phi_coarse_state,
    std::shared_ptr<SampleState> phi_state) {
  // Populate the fine level trial state
  // Step 1: coarse level points
  fine_action->copy_from_coarse(phi_coarse_state, theta_prime);
  // Step 2: fine level points from conditioned action
  conditioned_fine_action->fill_fine_points(theta_prime);
  /*
   * Calculate the difference in level-\ell actions,
   * required to calculate the ratio
   * \pi^{\ell}(\theta'_\ell)/\pi^{\ell}(\theta_\ell^{n})
   */
  double fine_action_theta_prime = fine_action->evaluate(theta_prime);
  double deltaS_fine = fine_action_theta_prime - fine_action_theta_fine;
  /*
   * Calculate the difference in level-(\ell-1) actions,
   * required to calculate the ratio
   * \pi^{\ell-1}(\theta_{\ell,C}^{n})/\pi^{\ell-1}(\theta'_{\ell,C})
   */
  coarse_action->copy_from_fine(theta_fine, theta_fine_C);

  double deltaS_coarse = coarse_action->evaluate(theta_fine_C) -
                         coarse_action->evaluate(phi_coarse_state);
  /*
   * Calculate the difference in free level-\ell actions,
   * required to calculate the ratio
   * q_{ML}^{\ell,F}(\theta_{\ell,F}^{n}|\theta'_{\ell,C}) /
   * q_{ML}^{\ell,F}(\theta'_{\ell,F}|\theta_{\ell,C}^{n})
   */
  double conditioned_fine_action_theta_prime =
      conditioned_fine_action->evaluate(theta_prime);
  double deltaS_trial =
      conditioned_fine_action_theta_fine - conditioned_fine_action_theta_prime;

  double deltaS = deltaS_fine + deltaS_coarse + deltaS_trial;
  accept = false;
  if (deltaS < 0.0) {
    accept = true;
  } else {
    double threshold = exp(-deltaS);
    accept = (uniform_dist(engine) < threshold);
  }
  if (accept) {
    theta_fine->data = theta_prime->data;
    fine_action_theta_fine = fine_action_theta_prime;
    conditioned_fine_action_theta_fine = conditioned_fine_action_theta_prime;
  }
  n_total_samples++;
  n_accepted_samples += (int)accept;
  // Copy back to state
  if (copy_if_rejected or accept) {
    phi_state->data = theta_fine->data;
  }
}

/* Set current state */
void TwoLevelMetropolisStep::set_state(std::shared_ptr<SampleState> phi_state) {
  theta_fine->data = phi_state->data;
  fine_action_theta_fine = fine_action->evaluate(theta_fine);
  conditioned_fine_action_theta_fine =
      conditioned_fine_action->evaluate(theta_fine);
}
