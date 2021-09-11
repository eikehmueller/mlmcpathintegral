#include "clustersampler.hh"

/** @file clustersampler.cc
 * @brief Implementation of clustersampler.hh
 */

/* Constructor */
QMClusterSampler::QMClusterSampler(const std::shared_ptr<QMClusterAction> action_,                                   
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
    std::cout << "HERE before" << std::endl;
    action->initialise_state(x_path_cur);
    std::cout << "HERE after" << std::endl;
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

/** Process next link */
std::pair<bool,int> QMClusterSampler::process_link(const int i,
        const int direction) {
    // Neighbouring site
    int i_neighbour = (i+direction + n_vertices) % n_vertices;
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
void QMClusterSampler::draw(std::shared_ptr<SampleState> x_path) {
    for (unsigned int i=0; i<n_updates; i++) {
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
    x_path->data = x_path_cur->data;
    accept = true;
}
