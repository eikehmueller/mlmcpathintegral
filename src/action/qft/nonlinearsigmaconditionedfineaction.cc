#include "nonlinearsigmaconditionedfineaction.hh"
/** @file nonlinearsigmaconditionedfineaction.cc
 * @brief Implementation of nonlinearsigmaconditionedfineaction.hh
 */

/* Fill in fine links */
void NonlinearSigmaConditionedFineAction::fill_fine_points(std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const std::vector<unsigned int>& fineonly_vertices = lattice->get_fineonly_vertices();
    // Iterate only over the fine-only vertices
    for (auto it=fineonly_vertices.begin();it!=fineonly_vertices.end();++it) {
        unsigned int ell = *it;
        // we need to use the _dof_-index in the call to the heatbath update. This
        // changes both dof 2*ell and 2*ell+1
        action->heatbath_update(phi_state,2*ell);
    }
}

/* Evaluate conditioned action at fine links */
double NonlinearSigmaConditionedFineAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    double S=0;
    Eigen::Vector3d sigma_n;
    Eigen::Vector3d Delta_n;
    const std::vector<unsigned int>& fineonly_vertices = lattice->get_fineonly_vertices();
    // Iterate only over the fine-only vertices
    for (auto it=fineonly_vertices.begin();it!=fineonly_vertices.end();++it) {
        unsigned int ell = *it;
        // Field at point n
        sigma_n.setZero();
        action->add_sigma(phi_state,ell,sigma_n);
        // Sum of neighbouring fields at point n
        Delta_n = action->delta_neighbours(phi_state,ell);
        // Work out length of Delta_n and angle between Delta_n and sigma_n
        double Delta_n_nrm = Delta_n.norm();
        double theta = acos(sigma_n.dot(Delta_n)/Delta_n_nrm);
        S -= log(exp_sin2_dist.evaluate(theta,2.*beta*Delta_n_nrm));
    }
    return S;
}
