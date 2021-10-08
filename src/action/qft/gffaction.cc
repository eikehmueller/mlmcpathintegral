#include "gffaction.hh"
/** @file gffaction.cc
 * @brief Implementation of gffaction.hh
 */

/* Value of action for a given configuration */
const double GFFAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0;
    double kappa = 2. + 0.5*mu2;
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    unsigned int Nvertices = lattice->getNvertices();
    for (unsigned int ell=0;ell<Nvertices;++ell) {
        double phi_n = phi_state->data[ell];
        // zero order term
        S += kappa*phi_n*phi_n;
        // Finite difference terms
        for (int k=0;k<4;++k) {
            S -= 0.5*phi_n * phi_state->data[neighbour_vertices[ell][k]];
        }
    }
    return S;
}

/* local heat bath update */
void GFFAction::heatbath_update(std::shared_ptr<SampleState> phi_state,
                                const unsigned int ell) {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();    
    double Delta =0.0;
    for (int k=0;k<4;++k) {
        Delta += 0.5*phi_state->data[neighbour_vertices[ell][k]];
    }
    phi_state->data[ell] = sigma*normal_dist(engine)+Delta/(4.+mu2);
}

/* local overrelaxation update */
void GFFAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                      const unsigned int ell) {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();    
    double Delta =0.0;
    for (int k=0;k<4;++k) {
        Delta += 0.5*phi_state->data[neighbour_vertices[ell][k]];
    }
    phi_state->data[ell] = 2.*Delta/(4.+mu2)-phi_state->data[ell];
}

/* Force for HMC integrator */
void GFFAction::force(const std::shared_ptr<SampleState> phi_state,
                      std::shared_ptr<SampleState> p_state) const {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    unsigned int Nvertices = lattice->getNvertices();
    double kappa = 4. + mu2;
    for (unsigned int ell=0;ell<Nvertices;++ell) {
        double phi_n = phi_state->data[ell];
        double momentum = kappa*phi_n;
        for (int k=0;k<4;++k) {
            momentum -= phi_state->data[neighbour_vertices[ell][k]];
        }
        p_state->data[ell] = momentum;
    }
}

/* Copy coarse unknowns from state on coarser level */
void GFFAction::copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                 std::shared_ptr<SampleState> phi_state) {
    const std::map<unsigned int, unsigned int>& fine2coarse_map = lattice->get_fine2coarse_map();
    for (auto it=fine2coarse_map.begin();it!=fine2coarse_map.end();++it) {
        unsigned int ell = it->first;
        unsigned int ell_coarse = it->second;
        phi_state->data[ell] = phi_coarse->data[ell_coarse];
    }
}

/* Copy unknowns from state on finer level */
void GFFAction::copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                               std::shared_ptr<SampleState> phi_state) {
    const std::map<unsigned int, unsigned int>& fine2coarse_map = fine_lattice->get_fine2coarse_map();
    for (auto it=fine2coarse_map.begin();it!=fine2coarse_map.end();++it) {
        unsigned int ell_fine = it->first;
        unsigned int ell = it->second;
        phi_state->data[ell] = phi_fine->data[ell_fine];
    }
}

/* Initialise state with random entries */
void GFFAction::initialise_state(std::shared_ptr<SampleState> phi_state) const {
    // IMPLEMENT THIS
}

/* Return lattice information */
std::string GFFAction::info_string() const {
    std::stringstream sstr;
    sstr << QFTAction::info_string() << ", mu2 = " << mu2;
    return sstr.str();
}
