#ifndef RENORMALISATION_HH
#define RENORMALISATION_HH RENORMALISATION_HH
#include <math.h>
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "lattice/lattice1d.hh"
#include "mpi/mpi_wrapper.hh"

/** @file renormalisation.hh
 * @brief Methods for calculation renormalised parameters
 */

/** @brief Enum for renormalisations:
 *  - 0: No renormalisation
 *  - 1: Perturbative renormalisation
 *  - 2: Exact renormalisation
*/
enum RenormalisationType {
    RenormalisationNone = 0,
    RenormalisationPerturbative = 1,
    RenormalisationExact = 2
};

/** @class RenormalisedParameters
 * @brief Base class for renormalised parameters
 *
 */

class RenormalisedParameters {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Underlying lattice
     * @param[in] renormalisation_ Type of renormalisation to use
     *              (0: none, 1: perturbative, 2: exact)
     */
    RenormalisedParameters(const std::shared_ptr<Lattice1D> lattice_,
                           const RenormalisationType renormalisation_) :
        lattice(lattice_), renormalisation(renormalisation_) {}

protected:
    /** @brief Number of time slices */
    const std::shared_ptr<Lattice1D> lattice;
    /** @brief Type of renormalisation */
    const RenormalisationType renormalisation;
};

#endif // RENORMALISATION_HH
