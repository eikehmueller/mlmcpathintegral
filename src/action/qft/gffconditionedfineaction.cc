#include "gffconditionedfineaction.hh"
/** @file gffconditionedfineaction.cc
 * @brief Implementation of gffconditionedfineaction.hh
 */

/* Fill in unknowns at fine vertices */
void GFFConditionedFineAction::fill_fine_points(std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    const std::vector<unsigned int>& fineonly_vertices = lattice->get_fineonly_vertices();
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    double sigma = 1./sqrt(4.+mu2);
    // Iterate only over the fine-only vertices
    for (auto it=fineonly_vertices.begin();it!=fineonly_vertices.end();++it) {
        unsigned int ell = *it;
        double Delta = 0.0;
        for (int k=0;k<4;++k) {
            Delta += phi_state->data[neighbour_vertices[ell][k]];
        }        
        phi_state->data[ell] = sigma*(normal_dist(engine) + sigma*Delta);
    }
}

/* Evaluate conditioned action at fine vertices */
double GFFConditionedFineAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    std::shared_ptr<Lattice2D> lattice = action->get_lattice();
    double S=0;
    const std::vector<unsigned int>& fineonly_vertices = lattice->get_fineonly_vertices();
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    double sigma2 = 1./(4.+mu2);
    double sigma2_inv = 1./sigma2;
    // Iterate only over the fine-only vertices
    for (auto it=fineonly_vertices.begin();it!=fineonly_vertices.end();++it) {
        unsigned int ell = *it;
        double Delta = 0.0;
        for (int k=0;k<4;++k) {
            Delta += phi_state->data[neighbour_vertices[ell][k]];
        }
        double dphi = phi_state->data[ell] - sigma2*Delta;
        S += 0.5*sigma2_inv*dphi*dphi;
    }
    return S;
}