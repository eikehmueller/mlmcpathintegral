#include "nonlinearsigmaaction.hh"
/** @file nonlinearsigmaaction.cc
 * @brief Implementation of nonlinearsigmaaction.hh
 */

/* Value of action for a given configuration */
const double NonlinearSigmaAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0;
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n;
    for (unsigned int ell=0;ell<lattice->getNvertices();++ell) {
        sigma_n.setZero();
        add_sigma(phi_state,ell,sigma_n);
        // Sum of neighbouring fields at point n
        Delta_n = delta_neighbours(phi_state,ell);
        // Add dot-product of sigma_n and Delta_n to action
        S += sigma_n.dot(Delta_n);
    }
    return -0.5*beta*S;
}

/* local heat bath update */
void NonlinearSigmaAction::heatbath_update(std::shared_ptr<SampleState> phi_state,
                                           const unsigned int ell) {
    // temporary vectors
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n; // sum of nearest neighbour vectors
    Eigen::Vector3d Delta_n_hat; // unit vector pointing in direction Delta_n
    Eigen::Vector3d Delta_n_perp; // vector perpendicular to Delta_n
    double Delta_n_nrm;
    Delta_n = delta_neighbours(phi_state,ell);
    Delta_n_nrm = Delta_n.norm();
    // Unit vector pointing in same direction as Delta_n
    Delta_n_hat = Delta_n;
    Delta_n_hat.normalize();
    // Work out the 'best' vector which is perpendicular to Delta_n
    // ('best' = largest angle to Delta_n)
    double tmp[3];
    tmp[0] = abs(Delta_n_hat[0]);
    tmp[1] = abs(Delta_n_hat[1]);
    tmp[2] = abs(Delta_n_hat[2]);
    int index = std::distance(tmp,std::min_element(tmp,tmp+3));
    double rho_inv = 1./sqrt(1.-tmp[index]*tmp[index]);
    switch(index) {
        case(0):
            Delta_n_perp[0] = 0.0;
            Delta_n_perp[1] = -Delta_n_hat[2]*rho_inv;
            Delta_n_perp[2] = +Delta_n_hat[1]*rho_inv;
            break;
        case(1):
            Delta_n_perp[0] = -Delta_n_hat[2]*rho_inv;
            Delta_n_perp[1] = 0.0;
            Delta_n_perp[2] = +Delta_n_hat[0]*rho_inv;
            break;
        case(2):
            Delta_n_perp[0] = +Delta_n_hat[1]*rho_inv;
            Delta_n_perp[1] = -Delta_n_hat[0]*rho_inv;
            Delta_n_perp[2] = 0.0;
            break;
    }
    // Draw altitude and azimuth rotation angles
    double theta = exp_sin2_dist.draw(engine,2.*beta*Delta_n_nrm);
    double phi = uniform_dist(engine);
    sigma_n = Eigen::AngleAxisd(phi, Delta_n_hat)
            * Eigen::AngleAxisd(theta, Delta_n_perp)
            * Delta_n_hat;
    phi = atan2(sigma_n[1],sigma_n[0]);
    theta = atan2(sqrt(sigma_n[0]*sigma_n[0]+sigma_n[1]*sigma_n[1]),sigma_n[2]);
    phi_state->data[2*ell] = theta;
    phi_state->data[2*ell+1] = phi;
}

/* local overrelaxation update */
void NonlinearSigmaAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                                 const unsigned int ell) {
    // temporary vectors
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n; // sum of nearest neighbour vectors
    double Delta_n_nrm;
    // Field at point n
    sigma_n.setZero();
    add_sigma(phi_state,ell,sigma_n);
    // Sum of nearest neihbours
    Delta_n = delta_neighbours(phi_state,ell);
    Delta_n.normalize();
    // Rotate around vector Delta_n (this does not change the action)
    double phi = uniform_dist(engine);
    sigma_n = Eigen::AngleAxisd(phi, Delta_n) * sigma_n;
    phi = atan2(sigma_n[1],sigma_n[0]);
    double theta = atan2(sqrt(sigma_n[0]*sigma_n[0]+sigma_n[1]*sigma_n[1]),sigma_n[2]);
    phi_state->data[2*ell] = theta;
    phi_state->data[2*ell+1] = phi;
}

/* Force for HMC integrator */
void NonlinearSigmaAction::force(const std::shared_ptr<SampleState> phi_state,
                                 std::shared_ptr<SampleState> p_state) const {
    for (unsigned int ell=0;ell<p_state->data.size();++ell) {
        p_state->data[ell] = 0.0;
    }
    Eigen::Vector3d Delta_n;
    for (unsigned int ell=0;ell<lattice->getNvertices();++ell) {
        double theta = phi_state->data[2*ell];
        double phi = phi_state->data[2*ell+1];
        Delta_n = delta_neighbours(phi_state,ell);
        double dS_dtheta = -beta*((Delta_n[0]*cos(phi)+Delta_n[1]*sin(phi))*cos(theta)-Delta_n[2]*sin(theta));
        double dS_dphi = -beta*(-Delta_n[0]*sin(phi)+Delta_n[1]*cos(phi))*sin(theta);
        p_state->data[2*ell] = dS_dtheta;
        p_state->data[2*ell+1] = dS_dphi;
    }
}

/* Copy coarse links from state on coarser level */
void NonlinearSigmaAction::copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                            std::shared_ptr<SampleState> phi_state) {
    const std::map<unsigned int, unsigned int>& fine2coarse_map = lattice->get_fine2coarse_map();
    for (auto it=fine2coarse_map.begin();it!=fine2coarse_map.end();++it) {
        unsigned int ell = it->first;
        unsigned int ell_coarse = it->second;
        phi_state->data[2*ell] = phi_coarse->data[2*ell_coarse];
        phi_state->data[2*ell+1] = phi_coarse->data[2*ell_coarse+1];
    }
}

/* Copy coarse links from state on finer level */
void NonlinearSigmaAction::copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                                          std::shared_ptr<SampleState> phi_state) {
    const std::map<unsigned int, unsigned int>& fine2coarse_map = fine_lattice->get_fine2coarse_map();
    for (auto it=fine2coarse_map.begin();it!=fine2coarse_map.end();++it) {
        unsigned int ell_fine = it->first;
        unsigned int ell = it->second;
        phi_state->data[2*ell] = phi_fine->data[2*ell_fine];
        phi_state->data[2*ell+1] = phi_fine->data[2*ell_fine+1];
    }
}

/* Initialise state with random entries */
void NonlinearSigmaAction::initialise_state(std::shared_ptr<SampleState> phi_state) const {
    std::uniform_real_distribution<double> uniform(-1.0,1.0);
    for (unsigned int ell=0;ell<sample_size()/2;++ell) {
        // Draw a point inside shell with inner radius 0.1 and outer radius 1.0
        double sigma_n[3];
        double nrm_sq;
        do {
            sigma_n[0] = uniform(engine);
            sigma_n[1] = uniform(engine);
            sigma_n[2] = uniform(engine);
            nrm_sq = sigma_n[0]*sigma_n[0]+sigma_n[1]*sigma_n[1]+sigma_n[2]*sigma_n[2];
        } while ( (nrm_sq > 1.0) or (nrm_sq < 0.01) );
        // convert to angles
        double phi = atan2(sigma_n[1],sigma_n[0]);
        double theta = atan2(sqrt(sigma_n[0]*sigma_n[0]+sigma_n[1]*sigma_n[1]),sigma_n[2]);
        phi_state->data[2*ell] = theta;
        phi_state->data[2*ell+1] = phi;
    }
}

/* Return lattice information */
std::string NonlinearSigmaAction::info_string() const {
    std::stringstream sstr;
    sstr << QFTAction::info_string() << ", beta = " << beta;
    return sstr.str();
}
