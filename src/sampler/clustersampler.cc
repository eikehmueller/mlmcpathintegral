#include "clustersampler.hh"

/** @file clustersampler.cc
 * @brief Implementation of clustersampler.hh
 */

/* Constructor */
ClusterSampler::ClusterSampler(const std::shared_ptr<ClusterAction> action_,                                   
                               const ClusterParameters cluster_param) :
    Sampler(),
    action(action_),
    lattice(action->get_generic_lattice()),
    n_vertices(lattice->getNvertices()),
    n_burnin(cluster_param.n_burnin()),
    n_updates(cluster_param.n_updates()),
    uniform_dist(0.0,1.0),
    uniform_int_dist(0,n_vertices-1) {
    engine.seed(2141517);
    // Create temporary workspace
    x_path_cur = std::make_shared<SampleState>(action->sample_size());
    action->initialise_state(x_path_cur);
    // Measure cost per sample
    Timer timer_meas;
    unsigned int n_meas = 10000;
    timer_meas.start();
    for (unsigned int k=0; k<n_meas; ++k) {
        draw(x_path_cur);
    }
    timer_meas.stop();
    cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
    // Burn in
    std::shared_ptr<SampleState> x_path_tmp = std::make_shared<SampleState>(action->sample_size());
    for (unsigned int i=0; i<n_burnin; ++i) {
        draw(x_path_tmp);
    }
}

/** Draw next sample */
void ClusterSampler::draw(std::shared_ptr<SampleState> x_path) {
    for (unsigned int i=0; i<n_updates; i++) {
        if (lattice->get_dim()==1) {
            single_cluster_update1d();
        } else {
            single_cluster_update();
        }
    }
    n_total_samples+=n_updates;
    n_accepted_samples+=n_updates;
    // Copy to output vector
    x_path->data = x_path_cur->data;
    accept = true;
}

/* Carry out a single generic cluster update */
void ClusterSampler::single_cluster_update() {
    // Pick new subgroup
    action->new_reflection();
    // Get neighbour list
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    std::set<unsigned int> cluster;
    // Pick a random site and flip it
    unsigned int ell0 = uniform_int_dist(engine);
    action->flip(x_path_cur,ell0);
    cluster.insert(ell0);
    // Queue with candidate vertices
    std::queue<unsigned int> active_vertices;
    active_vertices.push(ell0);
    while (not active_vertices.empty()) {
        // Get next vertex
        unsigned int ell = active_vertices.front();
        active_vertices.pop();
        // Loop over all neighbours of currently considered vertex
        for (auto it=neighbour_vertices[ell].begin();it!=neighbour_vertices[ell].end();++it) {
            unsigned int ell_neighbour = *it;
            // If neighbour is not in cluster yet, process link
            if (cluster.count(ell_neighbour)==0) {
                double Sell = action->S_ell(x_path_cur,ell,ell_neighbour);
                // Connection probability
                double p_connect = 1.-exp(fmin(0,-Sell));
                // Flip neighbouring site if it is bonded and add to lists
                bool bonded = (uniform_dist(engine)<p_connect);
                if (bonded) {
                    action->flip(x_path_cur,ell_neighbour);
                    cluster.insert(ell_neighbour);
                    active_vertices.push(ell_neighbour);
                }
            }
        }
    }
}

/* Carry out a single 1d cluster update */
void ClusterSampler::single_cluster_update1d() {
    // Pick new subgroup
    action->new_reflection();
    // Pick a random site and flip it
    int i0 = uniform_int_dist(engine);
    action->flip(x_path_cur,i0);
    // Process cluster forward
    int i_p=i0;
    int i_last_p;
    bool bonded;
    do {
        i_last_p = i_p;
        std::pair<bool,int> processed_link = process_link1d(i_p,+1);
        bonded = processed_link.first;
        i_p = processed_link.second;
    } while ( (i_p != i0) and bonded );
    // Process cluster backwards
    int i_m = i0;
    do {
        std::pair<bool,int> processed_link = process_link1d(i_m,-1);
        bonded = processed_link.first;
        i_m = processed_link.second;
    } while ( (i_m != i_last_p) and bonded );
}

/** Process next link in 1d algorithm */
std::pair<bool,int> ClusterSampler::process_link1d(const int i,
                                                   const int direction) {
    // Neighbouring site
    int i_neighbour = (i+direction + n_vertices) % n_vertices;
    // Check if neighbouring site is bonded
    double Sell = action->S_ell(x_path_cur,i,i_neighbour);
    // Connection probability
    double p_connect = 1.-exp(fmin(0,-Sell));
    // Flip neighbouring site if it is bonded
    bool bonded = (uniform_dist(engine)<p_connect);
    if (bonded) {
        action->flip(x_path_cur,i_neighbour);
    }
    return std::make_pair(bonded,i_neighbour);
}
