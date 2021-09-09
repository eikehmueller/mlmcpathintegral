#include "lattice/lattice2d.hh"

/** Implementation of lattice2d.hh */

Lattice2D::Lattice2D(const unsigned int Mt_lat_,
                     const unsigned int Mx_lat_,
                     const CoarseningType coarsening_type_,
                     const int coarsening_level_) : 
        Lattice(coarsening_level_),
        Mt_lat(Mt_lat_),
        Mx_lat(Mx_lat_),
        coarsening_type(coarsening_type_),
        rotated(false) {
    // Construct coarse lattice
    int rho_coarsen_t;   // temporal coarsening factor
    int rho_coarsen_x;   // spatial coarsening factor
    bool coarse_rotated; // Is the coarse lattice rotated?
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
                rho_coarsen_t = 2;
                rho_coarsen_x = 2;
                coarse_rotated = false;
            } else {
                rho_coarsen_t = 1;
                rho_coarsen_x = 1;
                coarse_rotated = true;
            }
        break;
    }
    // Check that lattice can be coarsened
    bool coarsening_allowed = true;
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
    if (coarsening_allowed) {
        coarse_lattice = std::make_shared<Lattice2D>(Mt_lat_coarse,
                                                     Mx_lat_coarse,
                                                     coarsening_type,
                                                     coarsening_level+1);
        // Rotate coarse lattice if necessary
        if (coarsening_type == CoarsenRotate) {
            coarse_lattice->rotated = coarse_rotated;
        }
    } else {
        coarse_lattice = nullptr;
    }
}
        