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
            // the rotated coarse lattice
            if (is_fillin_point(i,j)) {
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
    double S=0;
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n;
    for (int i=0;i<Mt_lat;++i) {
        for (int j=0;j<Mx_lat;++j) {
            // If action is formulated on rotated lattice, only consider points
            // at double-odd vertices. Otherwise, consider all points not on
            // the rotated coarse lattice
            if (is_fillin_point(i,j)) {
                // Field at point n
                sigma_n.setZero();
                action->add_sigma(phi_state,i,j,sigma_n);
                // Sum of neighbouring fields at point n
                Delta_n = action->delta_neighbours(phi_state,i,j);
                // Work out length of Delta_n and angle between Delta_n and sigma_n
                double Delta_n_nrm = Delta_n.norm();
                double theta = acos(sigma_n.dot(Delta_n)/Delta_n_nrm);
                S -= log(exp_sin2_dist.evaluate(theta,2.*beta*Delta_n_nrm));
            }
        }
    }
    return S;
}