#include "clustersampler.hh"

/** @file clustersampler.cc
 * @brief Implementation of clustersampler.hh
 */

/** Draw next sample */
void ClusterSampler::draw(std::vector<Path*> x_path) {
  const unsigned int M_lat = action.getM_lat();
  // Pick new subgroup
  action.new_subgroup();
  for (int i=0;i<M_lat;++i) {
    double dE = action.link_E(x_path_cur->data[i],
                              x_path_cur->data[(i+1)%M_lat])
      - action.link_Q(x_path_cur->data[i],
                      x_path_cur->data[(i+1)%M_lat]);
    // Connection probability
    double p_connect = 1.-exp(dE);
    bonds[i] = (uniform_dist(engine)<p_connect);
  }
  // Find first cluster
  int i0=0;
  while (bonds[(i0-1)%M_lat]) {
    i0 = (i0-1)%M_lat;
    if (i0==0) break;
  }
  int i = i0;
  action.new_subgroup_element();
  while (not ( (i+1)%M_lat == i0 ) ) {
    x_path_cur->data[i] = action.apply_current_subgroup_element(x_path_cur->data[i]);
      if (not bonds[i]) {
        action.new_subgroup_element();
      }
    i = (i+1)%M_lat;
  }
  // Copy to output vector
  std::copy(x_path_cur->data,
            x_path_cur->data+M_lat,
            x_path[0]->data);
}
