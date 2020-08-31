#include "twolevelmetropolisstep.hh"
/** @file twolevelmetropolissampler.cc
 * @brief Implementation of twolevelmetropolissampler.hh
 */
TwoLevelMetropolisStep::TwoLevelMetropolisStep(const std::shared_ptr<Action> coarse_action_,
                                               const std::shared_ptr<Action> fine_action_,
                                               const std::shared_ptr<ConditionedFineAction> conditioned_fine_action_) :
   MCMCStep(),
   coarse_action(coarse_action_),
   fine_action(fine_action_),
   conditioned_fine_action(conditioned_fine_action_),
   cost_per_sample_(0.0) {
   std::shared_ptr<Lattice1D> fine_lattice = fine_action->get_lattice();
   std::shared_ptr<Lattice1D> coarse_lattice = coarse_action->get_lattice();
   theta_fine = std::make_shared<Path>(fine_lattice);
   theta_fine_C = std::make_shared<Path>(coarse_lattice);
   theta_prime = std::make_shared<Path>(fine_lattice);
   engine.seed(89216491);
   std::shared_ptr<Path> x_fine = std::make_shared<Path>(fine_lattice);
   std::shared_ptr<Path> x_coarse = std::make_shared<Path>(coarse_lattice);

   fine_action_theta_fine = fine_action->evaluate(theta_fine);
   conditioned_fine_action_theta_fine = conditioned_fine_action->evaluate(theta_fine);
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
  // Populate the fine level trial state
  // Step 1: coarse level points
  theta_prime->copy_from_coarse(x_coarse_path);
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
  theta_fine_C->copy_from_fine(theta_fine);
  
  double deltaS_coarse = coarse_action->evaluate(theta_fine_C)
    - coarse_action->evaluate(x_coarse_path);
  /*
   * Calculate the difference in free level-\ell actions,
   * required to calculate the ratio                                       
   * q_{ML}^{\ell,F}(\theta_{\ell,F}^{n}|\theta'_{\ell,C}) /
   * q_{ML}^{\ell,F}(\theta'_{\ell,F}|\theta_{\ell,C}^{n})
   */
  double conditioned_fine_action_theta_prime = conditioned_fine_action->evaluate(theta_prime);
  double deltaS_trial = conditioned_fine_action_theta_fine - conditioned_fine_action_theta_prime;

  double deltaS = deltaS_fine + deltaS_coarse + deltaS_trial;
  accept = false;
  if (deltaS < 0.0) {
    accept = true;
  } else {
    double threshold = exp(-deltaS);
    accept = (uniform_dist(engine) < threshold);
  }
  if (accept) {
    theta_fine->copy(theta_prime);
    fine_action_theta_fine = fine_action_theta_prime;
    conditioned_fine_action_theta_fine = conditioned_fine_action_theta_prime;
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  // Copy back to path
  if (copy_if_rejected or accept) {
    x_path->copy(theta_fine);
  }
}

/* Set current state */
void TwoLevelMetropolisStep::set_state(std::shared_ptr<Path> x_path) {
  theta_fine->copy(x_path);
  fine_action_theta_fine = fine_action->evaluate(theta_fine);
  conditioned_fine_action_theta_fine = conditioned_fine_action->evaluate(theta_fine);
}
