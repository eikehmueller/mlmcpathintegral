#include "clustersampler.hh"

/** @file clustersampler.cc
 * @brief Implementation of clustersampler.hh
 */

/** Process next link */
std::pair<bool,int> ClusterSampler::process_link(const int i,
                                                 const int direction) {
  // Neighbouring site
  int i_neighbour = (i+direction)%x_path_cur->M_lat;
  // Check if neighbouring site is bonded
  double Sell = action.S_ell(x_path_cur->data[i],
                             x_path_cur->data[i_neighbour]);
  // Connection probability
  double p_connect = 1.-exp(fmin(0,-Sell));
  // Flip neighbouring site if it is bonded
  bool bonded = (uniform_dist(engine)<p_connect);
  if (bonded) {
    x_path_cur->data[i_neighbour] = action.flip(x_path_cur->data[i_neighbour]);
  }
  return std::make_pair(bonded,i_neighbour);
}

/** Draw next sample */
void ClusterSampler::draw(std::vector<std::shared_ptr<Path>> x_path) {
  const unsigned int M_lat = action.getM_lat();
  for (int i=0;i<n_updates;i++) {
    // Pick new subgroup
    action.new_angle();
    // Pick a random site and flip it
    int i0 = uniform_int_dist(engine);
    x_path_cur->data[i0] = action.flip(x_path_cur->data[i0]);
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
    if (record_stats) {
      n_total_samples++;
      n_accepted_samples++;
    }
  }
  // Copy to output vector
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path[0]->data);
}
