#include "nonlinearsigmaconditionedfineaction.hh"
/** @file nonlinearsigmaconditionedfineaction.cc
 * @brief Implementation of nonlinearsigmaconditionedfineaction.hh
 */

/* Fill in unknowns at fine vertices */
void NonlinearSigmaConditionedFineAction::fill_fine_points(
    std::shared_ptr<SampleState> phi_state) const {
  std::shared_ptr<Lattice2D> lattice = action->get_lattice();
  const std::vector<unsigned int> &fineonly_vertices =
      lattice->get_fineonly_vertices();
  // Iterate only over the fine-only vertices
  for (auto it = fineonly_vertices.begin(); it != fineonly_vertices.end();
       ++it) {
    unsigned int ell = *it;
    action->heatbath_update(phi_state, ell);
  }
}

/* Evaluate conditioned action at fine vertices */
double NonlinearSigmaConditionedFineAction::evaluate(
    const std::shared_ptr<SampleState> phi_state) const {
  std::shared_ptr<Lattice2D> lattice = action->get_lattice();
  double S = 0;
  Eigen::Vector3d Delta_n;
  const std::vector<unsigned int> &fineonly_vertices =
      lattice->get_fineonly_vertices();
  // Iterate only over the fine-only vertices
  for (auto it = fineonly_vertices.begin(); it != fineonly_vertices.end();
       ++it) {
    unsigned int ell = *it;
    // Field at point n
    Eigen::Vector3d sigma_n = action->get_sigma(phi_state, ell);
    // Sum of neighbouring fields at point n
    Delta_n = action->delta_neighbours(phi_state, ell);
    // Work out length of Delta_n and angle between Delta_n and sigma_n
    double Delta_n_nrm = Delta_n.norm();
    Delta_n.normalize();
    double z = sigma_n.dot(
        Delta_n); // projection of spin onto direction of neighbour sum
    S -= log(compact_exp_dist.evaluate(z, beta * Delta_n_nrm));
  }
  return S;
}
