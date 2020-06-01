#include "clustersampler.hh"

/** @file clustersampler.cc
 * @brief Implementation of clustersampler.hh
 */

/* Constructor */
ClusterSampler::ClusterSampler(const std::shared_ptr<ClusterAction> action_,
                               const unsigned int n_burnin_,
                               const unsigned int n_updates_) :
  Sampler(),
  action(action_),
  n_burnin(n_burnin_),
  n_updates(n_updates_),
  uniform_dist(0.0,1.0),
  uniform_int_dist(0,action_->getM_lat()-1)
{
  engine.seed(2141517);
  const unsigned int M_lat = action->getM_lat();
  const double T_final = action->getT_final();
  // Create temporary workspace
  x_path_cur = std::make_shared<Path>(M_lat,T_final);
  action->initialise_path(x_path_cur);
  // Measure cost per sample
  Timer timer_meas;
  unsigned int n_meas = 10000;
    timer_meas.start();
  for (unsigned int k=0;k<n_meas;++k) {
    draw(x_path_cur);
  }
  timer_meas.stop();
  cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
  // Burn in
  std::shared_ptr<Path> x_path_tmp = std::make_shared<Path>(M_lat,T_final);
  for (unsigned int i=0;i<n_burnin;++i) {
    draw(x_path_tmp);
  }
}

/** Process next link */
std::pair<bool,int> ClusterSampler::process_link(const int i,
                                                 const int direction) {
  // Neighbouring site
  int i_neighbour = (i+direction)%x_path_cur->M_lat;
  // Check if neighbouring site is bonded
  double Sell = action->S_ell(x_path_cur->data[i],
                             x_path_cur->data[i_neighbour]);
  // Connection probability
  double p_connect = 1.-exp(fmin(0,-Sell));
  // Flip neighbouring site if it is bonded
  bool bonded = (uniform_dist(engine)<p_connect);
  if (bonded) {
    x_path_cur->data[i_neighbour] = action->flip(x_path_cur->data[i_neighbour]);
  }
  return std::make_pair(bonded,i_neighbour);
}

/** Draw next sample */
void ClusterSampler::draw(std::shared_ptr<Path> x_path) {
  for (unsigned int i=0;i<n_updates;i++) {
    // Pick new subgroup
    action->new_angle();
    // Pick a random site and flip it
    int i0 = uniform_int_dist(engine);
    x_path_cur->data[i0] = action->flip(x_path_cur->data[i0]);
    // Process cluster forward
    int i_p=i0;
    int i_last_p;
    bool bonded;
    do {
      i_last_p = i_p;
      std::pair<bool,int> processed_link = process_link(i_p,+1);
      bonded = processed_link.first;
      i_p = processed_link.second;
    } while ( (i_p != i0) and bonded );
    // Process cluster backwards
    int i_m = i0;
    do {
      std::pair<bool,int> processed_link = process_link(i_m,-1);
      bonded = processed_link.first;
      i_m = processed_link.second;
    } while ( (i_m != i_last_p) and bonded );
    n_total_samples++;
    n_accepted_samples++;
  }
  // Copy to output vector
  x_path->copy(x_path_cur);
  accept = true;
}
