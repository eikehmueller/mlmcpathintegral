#include "hmcsampler.hh"

/** @file hmcsampler.cc
 * @brief Implementation of hmcsampler.hh
 */

/** Draw next sample */
void HMCSampler::draw(std::shared_ptr<SampleState> x_path) {
  accept = false;
  for (unsigned int r=0;r<n_rep;++r) {
    accept = accept or single_step();
  }
  n_total_samples++;
  n_accepted_samples += (int) accept;
  // Copy to output vector
  if (copy_if_rejected or accept) {
    x_path->data = x_path_cur->data;
  }
}

/* Do a single HMC accept/reject step */
bool HMCSampler::single_step() {
  // Initial kinetic energy
  std::generate(p_path_cur->data.data(),
                p_path_cur->data.data()+p_path_cur->data.size(),
                [this]() {return normal_dist(engine);});
  
  double T_kin_cur = 0.5*p_path_cur->data.squaredNorm();
  // STEP 1: Integrate deterministic trajectories with symplectic Euler method
  // Copy current state to trial state
  x_path_trial->data = x_path_cur->data;
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
    p_path_cur->data -= dt_p*dp_path->data;
    x_path_trial->data += dt_x*p_path_cur->data;    
  }
  // Calculate kinetic energy from trial state at end of trajectory
  double T_kin_trial = 0.5*p_path_cur->data.squaredNorm();
  // STEP 2: Accept-reject step
  bool accept_step = false;
  // Change in action S(X)
  double deltaS = action->evaluate(x_path_trial) - action->evaluate(x_path_cur);
  // Change in kinetic energy T(P)
  double deltaT = T_kin_trial - T_kin_cur;
  // Change in H(X,P) = T(P) + S(X)
  double deltaH = deltaS + deltaT;
  if (deltaH < 0.0) {
    accept_step = true;
  } else {
    double threshold = exp(-deltaH);
    accept_step = (uniform_dist(engine) < threshold);
  }
  // If accepted, copy state
  if (accept_step) {
    x_path_cur->data = x_path_trial->data;
  }
  return accept_step;
}

/* Set current state */
void HMCSampler::set_state(std::shared_ptr<SampleState> x_path) {
  x_path_cur->data = x_path->data;
}

/* automatically tune stepsize */
void HMCSampler::autotune_stepsize(const double p_accept_target) {
  unsigned int n_autotune_samples=1000;
  double tolerance = 1.E-2; // Tolerance
  double dt_hmc_original = dt_hmc;
  double dt_hmc_min = 0.5*dt_hmc;
  double dt_hmc_max = 2.*dt_hmc;
  bool converged=false;
  mpi_parallel::cout << std::setprecision(4) << std::fixed;
  mpi_parallel::cout << " Auto-tuning HMC step size to achieve acceptance rate of " << p_accept_target << " ..." << std::endl;
  mpi_parallel::cout << "  Starting with dt_{HMC} = " << dt_hmc << std::endl;
  for (unsigned int k=0;k<100;++k) {
    reset_stats();
    dt_hmc = 0.5*(dt_hmc_min+dt_hmc_max);
    for (unsigned int j=0;j<n_autotune_samples;++j) {
      n_accepted_samples += (int) single_step();
      n_total_samples++;
    }
    if (p_accept() > p_accept_target) {
      dt_hmc_min = dt_hmc;
    } else {
      dt_hmc_max = dt_hmc;
    }
    if (fabs(p_accept()-p_accept_target) < tolerance) converged=true;
  }
  if (converged) {
    mpi_parallel::cout << "  Tuned         dt_{HMC} = " << dt_hmc << std::endl;
  } else {
    dt_hmc = dt_hmc_original;
    mpi_parallel::cout << "  FAILED to tune, reverting to " << dt_hmc << std::endl;
  }
  mpi_parallel::cout << std::endl;
  reset_stats();
}
