#include "twolevelmetropolisstep.hh"
/** @file twolevelmetropolissampler.cc
 * @brief Implementation of twolevelmetropolissampler.hh
 */
TwoLevelMetropolisStep::TwoLevelMetropolisStep(const std::shared_ptr<Action> coarse_action_, const std::shared_ptr<Action> fine_action_, const std::shared_ptr<ConditionedFineAction> conditioned_fine_action_) :
   MCMCStep(),
   coarse_action(coarse_action_),
   fine_action(fine_action_),
   conditioned_fine_action(conditioned_fine_action_),
   cost_per_sample_(0.0) {
   assert(2*coarse_action->getM_lat()==fine_action->getM_lat());
   theta_fine = std::make_shared<Path>(fine_action->getM_lat(),
                                       coarse_action->getT_final());
   theta_fine_C = std::make_shared<Path>(coarse_action->getM_lat(),
                                         coarse_action->getT_final());
   theta_prime = std::make_shared<Path>(fine_action->getM_lat(),
                                        coarse_action->getT_final());
   engine.seed(89216491);
   std::shared_ptr<Path> x_fine = std::make_shared<Path>(fine_action->getM_lat(),
                                   coarse_action->getT_final());
   std::shared_ptr<Path> x_coarse = std::make_shared<Path>(coarse_action->getM_lat(),
                                     coarse_action->getT_final());

   Timer timer_meas;
   unsigned int n_meas = 10000;
   timer_meas.start();
   for (unsigned int k=0;k<n_meas;++k) {
       draw(x_coarse,x_fine);
   }
   timer_meas.stop();
   cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
   reset_stats();
 }

/** Draw new sample pair */
void TwoLevelMetropolisStep::draw(const std::shared_ptr<Path> x_coarse_path,
                                  std::shared_ptr<Path> x_path) {
  unsigned int M_lat = fine_action->getM_lat();
  // Populate the fine level trial state
  // Step 1: coarse level points
  for (unsigned int j=0; j<M_lat/2; ++j) {
    theta_prime->data[2*j] = x_coarse_path->data[j];    
  }
  // Step 2: fine level points from conditioned action
  conditioned_fine_action->fill_fine_points(theta_prime);
  /*
   * Calculate the difference in level-\ell actions,
   * required to calculate the ratio
   * \pi^{\ell}(\theta'_\ell)/\pi^{\ell}(\theta_\ell^{n}) 
   */
  double deltaS_fine = fine_action->evaluate(theta_prime)
    - fine_action->evaluate(theta_fine);
  /*
   * Calculate the difference in level-(\ell-1) actions,
   * required to calculate the ratio   
   * \pi^{\ell-1}(\theta_{\ell,C}^{n})/\pi^{\ell-1}(\theta'_{\ell,C}) 
   */
  for (unsigned int j=0; j<M_lat/2; ++j) {
    theta_fine_C->data[j] = theta_fine->data[2*j];
  }
  double deltaS_coarse = coarse_action->evaluate(theta_fine_C)
    - coarse_action->evaluate(x_coarse_path);
  /*
   * Calculate the difference in free level-\ell actions,
   * required to calculate the ratio                                       
   * q_{ML}^{\ell,F}(\theta_{\ell,F}^{n}|\theta'_{\ell,C}) /
   * q_{ML}^{\ell,F}(\theta'_{\ell,F}|\theta_{\ell,C}^{n})
   */
  double deltaS_trial = conditioned_fine_action->evaluate(theta_fine) \
    - conditioned_fine_action->evaluate(theta_prime);

  double deltaS = deltaS_fine + deltaS_coarse + deltaS_trial;
  accept = false;
  if (deltaS < 0.0) {
    accept = true;
  } else {
    double threshold = exp(-deltaS);
    accept = (uniform_dist(engine) < threshold);
  }
  if (accept) {
    std::copy(theta_prime->data,
              theta_prime->data+M_lat,
              theta_fine->data);
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  // Copy back to path
  if (copy_if_rejected or accept) {
    std::copy(theta_fine->data,
              theta_fine->data+M_lat,
              x_path->data);
  }
}

/* Set current state */
void TwoLevelMetropolisStep::set_state(std::shared_ptr<Path> x_path) {
  const unsigned int M_lat = x_path->M_lat;
  std::copy(x_path->data,
            x_path->data+M_lat,
            theta_fine->data);
}
