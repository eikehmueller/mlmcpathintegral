#include "twolevelmetropolissampler.hh"
/** @file twolevelmetropolissampler.cc
 * @brief Implementation of twolevelmetropolissampler.hh
 */

/** Draw new sample pair */
void TwoLevelMetropolisSampler::draw(std::vector<std::shared_ptr<Path>> x_path) {
  unsigned int M_lat = fine_action->getM_lat();
  // Draw new coarse level sample
  std::vector<std::shared_ptr<Path>> coarse_path;
  coarse_path.push_back(theta_coarse);
  coarse_sampler->draw(coarse_path);
  // Populate the fine level trial state
  // Step 1: coarse level points
  for (unsigned int j=0; j<M_lat/2; ++j) {
    theta_prime->data[2*j] = theta_coarse->data[j];    
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
    - coarse_action->evaluate(theta_coarse);
  /*
   * Calculate the difference in free level-\ell actions,
   * required to calculate the ratio                                       
   * q_{ML}^{\ell,F}(\theta_{\ell,F}^{n}|\theta_{\ell,C}^{n}) /
   * q_{ML}^{\ell,F}(\theta'_{\ell,F}|\theta'_{\ell,C})
   */
  double deltaS_trial = conditioned_fine_action->evaluate(theta_fine) \
    - conditioned_fine_action->evaluate(theta_prime);

  double deltaS = deltaS_fine + deltaS_coarse + deltaS_trial;
  bool accept = false;
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
  if (record_stats) {
    n_total_samples++;
    n_accepted_samples += (int) accept;
  }
  // Copy back to path
  std::copy(theta_fine->data,
            theta_fine->data+M_lat,
            x_path[0]->data);
  std::copy(theta_coarse->data,
            theta_coarse->data+M_lat/2,
            x_path[1]->data);
}

