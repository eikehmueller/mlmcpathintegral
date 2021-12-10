#include "gffaction.hh"
/** @file gffaction.cc
 * @brief Implementation of gffaction.hh
 */

/* Value of action for a given configuration */
const double GFFAction::evaluate(const std::shared_ptr<SampleState> phi_state) const {
    double S=0.0;
    if (n_gibbs_smooth==0) {        
        // No Gibbs smoothing
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
    } else {
        // Compute \phi^T.\hat{Q}.\phi
        S = phi_state->data.dot(Q_precision_hat*phi_state->data);
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

/* global heat bath update with effective action */
void GFFAction::global_heatbath_update_eff(std::shared_ptr<SampleState> phi_state) const {
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    unsigned int Nvertices = lattice->getNvertices();
    double sigma_eff = 1./sqrt(4.+0.5*mu2 - 4./(4.+0.5*mu2));
    double kappa = 1./(4.+0.5*mu2);
    for (unsigned int ell=0;ell<Nvertices;++ell) {
        double Delta = 0.0;
        for (int k=0;k<4;++k) {
            Delta += 2.*kappa*phi_state->data[neighbour_vertices[ell][k]];
        }
        for (int k=4;k<8;++k) {
            Delta += kappa*phi_state->data[neighbour_vertices[ell][k]];
        }
        phi_state->data[ell] = sigma_eff*(normal_dist(engine)+Delta*sigma_eff);
    }
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

/* Build matrices required for direct sampling and Gibbs acceleration */ 
void GFFAction::buildMatrices() {
    unsigned int Nvertices = lattice->getNvertices();
    
    // Precision- and covariance- matrices
    std::vector<double> stencil {4.+mu2,-1.};
    Eigen::SparseMatrix<double> Q_precision;
    Q_precision = buildPrecisionMatrix(stencil);
    Eigen::MatrixXd Sigma = Eigen::MatrixXd(Q_precision).inverse();

    // Effective precision- and covariance- matrices
    std::vector<double> stencil_eff {4.+0.5*mu2 - 4./(4.+0.5*mu2),
                                     -2./(4.+0.5*mu2),
                                     -1./(4.+0.5*mu2)};
    Eigen::SparseMatrix<double> Q_precision_eff;
    Q_precision_eff = buildPrecisionMatrix(stencil_eff);
    Eigen::MatrixXd Sigma_eff = Eigen::MatrixXd(Q_precision_eff).inverse();
    
    Eigen::MatrixXd M_mat = Q_precision_eff.triangularView<Eigen::Lower>();
    Eigen::MatrixXd G_mat = Eigen::MatrixXd::Identity(Nvertices, Nvertices);
    if (n_gibbs_smooth>0) {
        Eigen::MatrixXd G_mat_tmp = (M_mat.inverse()*Q_precision_eff
                                    - Eigen::MatrixXd::Identity(Nvertices, Nvertices));
        for (int k=0;k<n_gibbs_smooth;++k) {
            G_mat = G_mat_tmp*G_mat;
        }
    }
    
    // Precision matrix after Gibbs smoothing
    Q_precision_hat = (Sigma_eff + G_mat*(Sigma-Sigma_eff)*G_mat.transpose()).inverse();
    
    // Cholesky factorisation for exact sampling
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Lower,Eigen::NaturalOrdering<int>> sparse_cholesky(Q_precision);
    choleskyLT = sparse_cholesky.matrixU();
    choleskyL = sparse_cholesky.matrixL();
}

/* Build precision matrix based on a given stencil */
Eigen::SparseMatrix<double> GFFAction::buildPrecisionMatrix(std::vector<double> stencil) {        
    unsigned int Nvertices = lattice->getNvertices();
    const std::vector<std::vector<unsigned int> >& neighbour_vertices = lattice->get_neighbour_vertices();
    size_t stencil_size = 1+4*(stencil.size()-1);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletlist(stencil_size*Nvertices);
    for (unsigned int ell=0;ell<Nvertices;++ell) {        
        tripletlist[stencil_size*ell] = T(ell,ell,stencil[0]);
        for (int j=0;j<stencil.size()-1;++j) {
            for (int k=0;k<4;++k) {
                tripletlist[stencil_size*ell+4*j+k+1] = T(ell,neighbour_vertices[ell][k],stencil[j+1]);
            }
        }
    }    
    Eigen::SparseMatrix<double> Q_prec(Nvertices,Nvertices);
    Q_prec.setFromTriplets(tripletlist.begin(),tripletlist.end());
    return Q_prec;
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
    for (int k=0;k<n_gibbs_smooth;++k) {
        global_heatbath_update_eff(phi_state);
    }
    n_total_samples++;
    n_accepted_samples++;
}

void GFFAction::print_correlation_function() const {
    unsigned int Nvertices = lattice->getNvertices();
    Eigen::VectorXd rhs;
    rhs.resize(Nvertices);
    // Draw uncorrelated sample vector psi
    for (unsigned int ell=0;ell<Nvertices;++ell) {
        rhs[ell] = 0.0;
    }    
    rhs[0] = 1.0;
    Eigen::VectorXd u = choleskyLT.triangularView<Eigen::Upper>().solve(choleskyL.triangularView<Eigen::Lower>().solve(rhs));
    std::vector<double> X;
    std::vector<double> Y;
    for (int i=0;i<lattice->getMt_lat();++i) {
        unsigned int ell = lattice->vertex_cart2lin(i,i);
        X.push_back(sqrt(2.)*i/(1.*lattice->getMt_lat()));
        Y.push_back(u[ell]);
    }
    printf("=== correlation function ===\n");
    for (int i=0;i<lattice->getMt_lat();++i) {
        printf("%8.4f %12.6e\n",X[i],Y[i]);
    }
}