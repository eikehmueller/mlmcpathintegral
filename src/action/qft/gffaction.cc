#include "gffaction.hh"
/** @file gffaction.cc
 * @brief Implementation of gffaction.hh
 */

/* Value of action for a given configuration */
const double GFFAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0.0;
    double kappa = 4. + mu2;
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    unsigned int Nvertices = lattice->getNvertices();
    for (unsigned int ell=0;ell<Nvertices;++ell) {        
        double phi_n = phi_state->data[ell];
        double S_local = kappa*phi_n;
        // nearest neighbour terms
        for (int k=0;k<4;++k) {
            S_local -= phi_state->data[neighbour_vertices[ell][k]];
        }
        S += phi_n*S_local;
    }
    return 0.5*S;
}

/* local heat bath update */
void GFFAction::heatbath_update(std::shared_ptr<SampleState> phi_state,
                                const unsigned int ell) {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();    
    double Delta =0.0;
    for (int k=0;k<4;++k) {
        Delta += phi_state->data[neighbour_vertices[ell][k]];
    }
    phi_state->data[ell] = sigma*normal_dist(engine)+Delta/(4.+mu2);
}

/* local overrelaxation update */
void GFFAction::overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                      const unsigned int ell) {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();    
    double Delta =0.0;
    for (int k=0;k<4;++k) {
        Delta += phi_state->data[neighbour_vertices[ell][k]];
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
    const_cast<GFFAction *>(this)->draw(phi_state);
}

/* Return lattice information */
std::string GFFAction::info_string() const {
    std::stringstream sstr;
    sstr << QFTAction::info_string() << ", mu2 = " << mu2;
    return sstr.str();
}

/* Build Cholesky factorisation for direct sampling */ 
void GFFAction::buildCholesky() {
    unsigned int Nvertices = lattice->getNvertices();
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    // Construct the precision matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletlist(5*Nvertices);
    double kappa = 4. + mu2;
    for (unsigned int ell=0;ell<Nvertices;++ell) {        
        tripletlist[5*ell]   = T(ell,ell,kappa);
        for (int k=0;k<4;++k) {
            tripletlist[5*ell+1+k] = T(ell,neighbour_vertices[ell][k],-1.0);
        }
    }    
    Eigen::SparseMatrix<double> Q_precision(Nvertices,Nvertices);
    Q_precision.setFromTriplets(tripletlist.begin(),tripletlist.end());
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > sparse_cholesky(Q_precision);
    choleskyLT = sparse_cholesky.matrixU();
}

/* Draw sample from true distribution */
void GFFAction::draw(std::shared_ptr<SampleState> phi_state) {
    unsigned int Nvertices = lattice->getNvertices();
    // Draw uncorrelated sample vector psi
    for (unsigned int ell=0;ell<Nvertices;++ell) {
        rhs_sample[ell] = normal_dist(engine);
    }
    // solve L^T.phi = psi to obtain correlated sample phi
    phi_state->data = choleskyLT.triangularView<Eigen::Upper>().solve(rhs_sample);
    n_total_samples++;
    n_accepted_samples++;
}