#include "lattice/lattice1d.hh"

/** Implementation of lattice1d.hh */

/* Constructor */
Lattice1D::Lattice1D(const unsigned int M_lat_, const double T_final_,
                     const int coarsening_level_)
    : Lattice(coarsening_level_, 1), M_lat(M_lat_), T_final(T_final_),
      a_lat(T_final_ / M_lat_) {
  assert(T_final > 0.0);
  // Construct neighbour list
  for (unsigned int ell = 0; ell < M_lat; ++ell) {
    std::vector<unsigned int> local_nb;
    for (int offset = -1; offset < 2; offset += 2) {
      local_nb.push_back((ell + offset + M_lat) % M_lat);
    }
    neighbour_vertices.push_back(local_nb);
  }
}
