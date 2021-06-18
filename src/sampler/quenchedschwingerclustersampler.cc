#include "quenchedschwingerclustersampler.hh"

/** @file quenchedschwingerclustersampler.cc
 * @brief Implementation of quenchedschwingerclustersampler.hh
 */

/* Constructor */
QuenchedSchwingerClusterSampler::QuenchedSchwingerClusterSampler(const std::shared_ptr<QuenchedSchwingerAction> action_,
                                                                 const ClusterParameters cluster_param) :
    Sampler(),
    action(action_),
    engine(723851),
    uniform_distribution(-M_PI,+M_PI) {
    // Create temporary workspace
    unsigned int ndof = action->sample_size();
    psi_cluster_state = std::make_shared<SampleState>(ndof);
    
    // Create rotor action and cluster sampler
    double beta = action->getbeta();    
    std::shared_ptr<Lattice1D> lattice1d = std::make_shared<Lattice1D>(ndof/2,1.0);
    rotor_action = std::make_shared<RotorAction>(lattice1d,
                                                 RenormalisationNone,
                                                 beta*lattice1d->geta_lat());
    cluster_sampler = std::make_shared<ClusterSampler>(rotor_action,cluster_param);
    rotor_action->initialise_state(psi_cluster_state);
    // Measure cost per sample
    std::shared_ptr<SampleState> phi_state_tmp = std::make_shared<SampleState>(ndof);
    Timer timer_meas;
    unsigned int n_meas = 10000;
    timer_meas.start();
    for (unsigned int k=0; k<n_meas; ++k) {
        draw(phi_state_tmp);
    }
    timer_meas.stop();
    cost_per_sample_ = 1.E6*timer_meas.elapsed()/n_meas;
}
/** Draw next sample */
void QuenchedSchwingerClusterSampler::draw(std::shared_ptr<SampleState> phi_state) {
    std::shared_ptr<Lattice2D> lattice2d = action->get_lattice();
    unsigned int Mt_lat = lattice2d->getMt_lat();
    unsigned int Mx_lat = lattice2d->getMx_lat();
    unsigned int ndof = action->sample_size();
    // Draw path with cluster sampler
    cluster_sampler->draw(psi_cluster_state); 
    for (unsigned int j=0; j<ndof; ++j) phi_state->data[j] = 0.0;
    // Set vertical links
    int i_lin = 0;
    for (int i=0;i<Mt_lat-1;++i) {
        for (int j=0;j<Mx_lat;++j) {            
            phi_state->data[lattice2d->link_cart2lin(i+1,j,1)] \
              = phi_state->data[lattice2d->link_cart2lin(i,j,1)] \
              + psi_cluster_state->data[i_lin+1] - psi_cluster_state->data[i_lin];
            i_lin++;
        }
    }
    // Set horizontal links
    for (int j=0;j<Mx_lat-1;++j) {
        phi_state->data[lattice2d->link_cart2lin(Mt_lat-1,j+1,0)] \
          = phi_state->data[lattice2d->link_cart2lin(Mt_lat-1,j,0)] \
          - phi_state->data[lattice2d->link_cart2lin(Mt_lat-1,j,1)] \
          - psi_cluster_state->data[i_lin+1] + psi_cluster_state->data[i_lin];
        i_lin++;
    }
    // Gauge transformation
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {            
            double theta = uniform_distribution(engine);
            phi_state->data[lattice2d->link_cart2lin(i,j,0)] \
              = mod_2pi(phi_state->data[lattice2d->link_cart2lin(i,j,0)]+theta);
            phi_state->data[lattice2d->link_cart2lin(i-1,j,0)] \
              = mod_2pi(phi_state->data[lattice2d->link_cart2lin(i-1,j,0)]-theta);
            phi_state->data[lattice2d->link_cart2lin(i,j,1)] \
              = mod_2pi(phi_state->data[lattice2d->link_cart2lin(i,j,1)]+theta);
            phi_state->data[lattice2d->link_cart2lin(i,j-1,1)] \
              = mod_2pi(phi_state->data[lattice2d->link_cart2lin(i,j-1,1)]-theta);
        }
    }
    accept = true;
    n_total_samples++;
    n_accepted_samples++;
}
