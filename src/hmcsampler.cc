#include "hmcsampler.hh"
#include <iostream>

/** @file hmcsampler.cc
 * @brief Implementation of hmcsampler.hh
 */

/** Draw next sample */
void HMCSampler::draw(std::shared_ptr<Path> x_path) {
  const unsigned int M_lat = action->getM_lat();
  // Initial kinetic energy
  double T_kin_cur = 0.0;
  // Draw random momentum from normal distribution
  for (unsigned int j=0;j<M_lat;++j) {
    double tmp = normal_dist(engine);
    p_path_cur->data[j] = tmp;
    T_kin_cur += 0.5*tmp*tmp;
  }
  // STEP 1: Integrate deterministic trajectories with symplectic Euler method
  // Copy current state to trial state
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path_trial->data);
  for(unsigned int k =0;k<=nt_hmc;++k) {
    double dt_p = dt_hmc;
    double dt_x = dt_hmc;
    if (k==0) {
      dt_p = 0.5*dt_hmc;
    }
    if (k==nt_hmc) {
      dt_p = 0.5*dt_hmc;
      dt_x = 0.0;
    }
    // Calculate force
    action->force(x_path_trial,dp_path);
    for(unsigned int j=0;j<M_lat;++j) {
      // Momentum update P -> P - 0.5*dt_{hmc}*dS(X)/dX
      p_path_cur->data[j] -= dt_p*dp_path->data[j];
      // Position update X -> X + dt_{hmc}*P
      x_path_trial->data[j] += dt_x*p_path_cur->data[j];
    }
  }
  // Calculate kinetic energy from trial state at end of trajectory
  double T_kin_trial = 0.0;
  for (unsigned int j=0;j<M_lat;++j) {
    double tmp = p_path_cur->data[j];
    T_kin_trial += 0.5*tmp*tmp;
  }
  // STEP 2: Accept-reject step
  accept = false;
  // Change in action S(X)
  double deltaS = action->evaluate(x_path_trial) - action->evaluate(x_path_cur);
  // Change in kinetic energy T(P)
  double deltaT = T_kin_trial - T_kin_cur;
  // Change in H(X,P) = T(P) + S(X)
  double deltaH = deltaS + deltaT;
  if (deltaH < 0.0) {
    accept = true;
  } else {
    double threshold = exp(-deltaH);
    accept = (uniform_dist(engine) < threshold);
  }
  // If accepted, copy state
  if (accept) {
    std::copy(x_path_trial->data,
              x_path_trial->data+M_lat,
              x_path_cur->data);
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  
  // Copy to output vector
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path->data);

}
