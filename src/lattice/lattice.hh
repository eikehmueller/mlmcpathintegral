#ifndef LATTICE_HH
#define LATTICE_HH LATTICE_HH
#include <memory>
#include <stdexcept>
#include <cassert>
#include "mpi/mpi_wrapper.hh"

/** @file lattice.hh
 * @brief Header file for abstract lattice base class
 */

/** @class Lattice
 *
 * @brief Abstract lattice base class
 *
 */
class Lattice {
public:
    /** @brief Initialise class
     *
     *
     * @param[in] coarsening_level_ Coarsening level (0=finest)
     * @param[in] T_final_ Final time \f$T\f$
     */
    Lattice(const int coarsening_level_=0) : coarsening_level(coarsening_level_) {}

    /** @brief Return coarsening level of action */
    const int get_coarsening_level() const {
        return coarsening_level;
    }

protected:
    /** @brief coarsening level (0=finest level) */
    mutable int coarsening_level;
};

#endif // LATTICE_HH
