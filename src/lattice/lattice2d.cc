#include "lattice/lattice2d.hh"

/** Implementation of lattice2d.hh */

/* constructor */
Lattice2D::Lattice2D(const unsigned int Mt_lat_,
                     const unsigned int Mx_lat_,
                     const CoarseningType coarsening_type_,
                     const int coarsening_level_) : 
        Lattice(coarsening_level_,2),
        Mt_lat(Mt_lat_),
        Mx_lat(Mx_lat_),
        coarsening_type(coarsening_type_),
        rotated( (coarsening_type_==CoarsenRotate) and (coarsening_level_%2) ) {            
    // Construct coarse lattice
    if ( (rotated) and ( ( Mx_lat%2) or (Mt_lat%2) ) ) {
        mpi_parallel::cerr << "ERROR: Both Mx_lat and Mt_lat have to be even for rotated lattices." << std::endl;
        mpi_exit(EXIT_FAILURE);
        throw std::runtime_error("...");
    }
    int rho_coarsen_t;   // temporal coarsening factor
    int rho_coarsen_x;   // spatial coarsening factor
    bool coarsening_allowed = true;    
    switch (coarsening_type) {
        case CoarsenBoth:
            // Coarsen in both directions
            rho_coarsen_t = 2;
            rho_coarsen_x = 2;
        break;
        case CoarsenTemporal:
            // Coarsen in temporal direction only
            rho_coarsen_t = 2;
            rho_coarsen_x = 1;
        break;
        case CoarsenSpatial:
            // Coarsen in spatial direction only
            rho_coarsen_t = 1;
            rho_coarsen_x = 2;
        break;
        case CoarsenAlternate:
            // Coarsen in alternate directions
            if (coarsening_level%2==0) {
                rho_coarsen_t = 2;
                rho_coarsen_x = 1;
            } else {
                rho_coarsen_t = 1;
                rho_coarsen_x = 2;
            }
        break;
        case CoarsenRotate:
            if (rotated) {
                if ( (Mt_lat%2) or (Mx_lat%2) ) {
                    coarsening_allowed = false;                    
                }                
                rho_coarsen_t = 2;
                rho_coarsen_x = 2;
            } else {
                rho_coarsen_t = 1;
                rho_coarsen_x = 1;
            }
        break;
    }
    // Check that lattice can be coarsened
    unsigned int Mt_lat_coarse = Mt_lat;
    unsigned int Mx_lat_coarse = Mx_lat;        
    if ( rho_coarsen_t > 1 ) {
        if (Mt_lat%rho_coarsen_t) {
            coarsening_allowed = false;
        }
        Mt_lat_coarse = Mt_lat/rho_coarsen_t;
    }
    if ( rho_coarsen_x > 1 ) {
        if (Mx_lat%rho_coarsen_x) {
            coarsening_allowed = false;
        }
        Mx_lat_coarse = Mx_lat/rho_coarsen_x;
    }
    coarsening_allowed = coarsening_allowed and (Mt_lat_coarse>1) and (Mx_lat_coarse>1);
    if (coarsening_allowed) {
        coarse_lattice = std::make_shared<Lattice2D>(Mt_lat_coarse,
                                                     Mx_lat_coarse,
                                                     coarsening_type,
                                                     coarsening_level+1);
        // list of coarse and fine-only vertices        
        if (coarsening_type == CoarsenRotate) {
            if (rotated) {
                for (int i=0;i<Mt_lat;++i) {
                    for (int j=0;j<Mx_lat;++j) {
                        if ( (i+j)%2 == 0 ) {
                            unsigned int ell = vertex_cart2lin(i,j);                       
                            if ( (i%2==0) and (j%2==0) ) {
                                coarse_vertices.push_back(ell);
                            } else  {
                                fineonly_vertices.push_back(ell);
                            }
                        }
                    }
                }
            } else {
                for (int i=0;i<Mt_lat;++i) {
                    for (int j=0;j<Mx_lat;++j) {
                        unsigned int ell = vertex_cart2lin(i,j);                       
                        if ( (i+j)%2 == 0 ) {
                            coarse_vertices.push_back(ell);
                        } else  {
                            fineonly_vertices.push_back(ell);
                        }
                    }
                }
            }
        } else {
            for (int i=0;i<Mt_lat;++i) {
                for (int j=0;j<Mx_lat;++j) {
                    unsigned int ell = vertex_cart2lin(i,j);
                    if ( (i%rho_coarsen_t==0) and (j%rho_coarsen_x==0) ) {
                        coarse_vertices.push_back(ell);
                    } else  {
                        fineonly_vertices.push_back(ell);
                    }
                }
            }
        }
        std::sort(coarse_vertices.begin(),coarse_vertices.end());
        std::sort(fineonly_vertices.begin(),fineonly_vertices.end());
        // list of coarse dofs corresponding to a particular fine vertex
        for (auto it=coarse_vertices.begin();it<coarse_vertices.end();++it) {
            unsigned int ell = *it;
            int i,j;
            vertex_lin2cart(ell,i,j);
            unsigned int ell_coarse = coarse_lattice->vertex_cart2lin(i/rho_coarsen_t,
                                                                      j/rho_coarsen_x);            
            fine2coarse_map[ell] = ell_coarse;
            
        }
    } else {
        coarse_lattice = nullptr;
    }
    // list of direct neighbour vertices
    std::array<int,4> offset_i;
    std::array<int,4> offset_j;
    if (rotated) {
        offset_i = {+1,+1,-1,-1};
        offset_j = {+1,-1,+1,-1};
    } else {
        offset_i = {+1,-1, 0, 0};
        offset_j = { 0, 0,+1,-1};
    }
    for (unsigned int ell=0;ell<getNvertices();++ell) {
        int i,j;
        vertex_lin2cart(ell,i,j);
        std::vector<unsigned int> local_nb;
        for (int k=0;k<4;++k) {
            local_nb.push_back(vertex_cart2lin(i+offset_i[k],j+offset_j[k]));
        }
        neighbour_vertices.push_back(local_nb);
    }
}
        