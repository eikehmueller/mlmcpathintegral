#ifndef LATTICE1D_HH
#define LATTICE1D_HH LATTICE1D_HH
#include <memory>
#include <stdexcept>
#include <cassert>
#include "mpi/mpi_wrapper.hh"

/** @file lattice1d.hh
 * @brief Header file for one-dimensional lattice class
 */

/** @class Lattice1D
 *
 * @brief Class for 1d lattice
 *
 */
class Lattice1D {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   */
  Lattice1D(const unsigned int M_lat_,
            const double T_final_,
            const int coarsening_level_=0) :
  M_lat(M_lat_),
  T_final(T_final_),
  a_lat(T_final_/M_lat_),
  coarsening_level(coarsening_level_) {
    assert(T_final>0.0);
  }

  /** @brief Return number of timeslices \f$M\f$ */
  unsigned int getM_lat() const { return M_lat; }

  /** @return final time \f$T\f$ */
  double getT_final() const { return T_final; }

  /** @brief Return lattice spacing \f$a\f$ */
  double geta_lat() const { return a_lat;}

  std::shared_ptr<Lattice1D> coarse_lattice() {
    if (M_lat%2) {
      mpi_parallel::cerr << "ERROR: cannot coarsen 1d lattice with M = " << M_lat << " points." << std::endl;
      mpi_exit(EXIT_FAILURE);
      throw std::runtime_error("...");
    }
    return std::make_shared<Lattice1D>(M_lat/2,T_final,coarsening_level+1);
  };

  /** @brief Return coarsening level of action */
  const int get_coarsening_level() const { return coarsening_level; }
  
protected:
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief lattice spacing */
  const double a_lat;
  /** @brief coarsening level (0=finest level) */
  mutable int coarsening_level;
};

#endif // LATTICE1D_HH
