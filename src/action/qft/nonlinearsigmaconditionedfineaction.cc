#include "nonlinearsigmaconditionedfineaction.hh"
/** @file nonlinearsigmaconditionedfineaction.cc
 * @brief Implementation of nonlinearsigmaconditionedfineaction.hh
 */

/* Fill in fine links */
void NonlinearSigmaConditionedFineAction::fill_fine_points(std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            // If action is formulated on rotated lattice, only consider points
            // at double-odd vertices. Otherwise, consider all points not on
            // the coarse lattice
            if (((action->is_rotated()) and (i&2) and (j&1)) or ((i+j)&1)) {
                action->heatbath_ij_update(phi_state,i,j);                
            }
        }
    }
}

/* Evaluate conditioned action at fine links */
double NonlinearSigmaConditionedFineAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const unsigned int Mt_lat = lattice->getMt_lat();
    const unsigned int Mx_lat = lattice->getMx_lat();
    // TO DO
    double S=0;
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            // If action is formulated on rotated lattice, only consider points
            // at double-odd vertices. Otherwise, consider all points not on
            // the coarse lattice
            if (((action->is_rotated()) and (i&2) and (j&1)) or ((i+j)&1)) {
                // Field at point n
                sigma_n.setZero();
                action->add_sigma(phi_state,i,j,sigma_n);
                // Sum of neighbouring fields at point n
                Delta_n = action->delta_neighbours(phi_state,i,j);
                // Add dot-product of sigma_n and Delta_n to action
                S += sigma_n.dot(Delta_n);
            }
        }
    }
    return -0.5*beta*S;
}
