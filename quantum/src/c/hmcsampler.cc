#include "hmcsampler.hh"
#include <iostream>

/** @file hmcsampler.cc
 * @brief Implementation of hmcsampler.hh
 */

/** Draw next sample */
void HMCSampler::draw(std::vector<Path*> x_path) {
  const unsigned int M_lat = action.getM_lat();
  const double T_final = action.getT_final();
  // Trial state
  Path* x_path_trial = new Path(M_lat,T_final);
  // Momentum change from force term
  Path* dp_path = new Path(M_lat,T_final);
  // Draw random momentum from normal distribution
  for (unsigned int j=0;j<M_lat;++j) {
    p_path_cur->data[j] = normal_dist(engine);
  }
  // STEP 1: Integrate deterministic trajectories with symplectic Euler method
  // Copy current state to trial state
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path_trial->data);
  const unsigned int hmc_steps = floor(T_hmc/dt_hmc);
  for (unsigned int k=0;k<hmc_steps;++k) {
    // Momentum update P -> P - dt_{hmc}*dS(X)/dX
    action.force(x_path_trial,dp_path);
    for(unsigned int j=0;j<M_lat;++j) {
      p_path_cur->data[j] -= dt_hmc*dp_path->data[j];
      // Position update X -> X + dt_{hmc}*P
      x_path_trial->data[j] += dt_hmc*p_path_cur->data[j];
    }
  }
  // STEP 2: Accept-reject step
  bool accept = false;
  double deltaS = action.evaluate(x_path_trial) - action.evaluate(x_path_cur);
  if (deltaS < 0.0) {
    accept = true;
  } else {
    double threshold = exp(-deltaS);
    accept = (uniform_dist(engine) < threshold);
  }
  // If accepted, copy state
  if (accept) {
    std::copy(x_path_trial->data,
              x_path_trial->data+M_lat,
              x_path_cur->data);
  }
  if (record_stats) {
    n_total_samples++;
    n_accepted_samples += (int) accept;
  }
  
  // Copy to output vector
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path[0]->data);

  // Deallocate memory
  delete x_path_trial;
  delete dp_path;
}
