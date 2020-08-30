#ifndef RENORMALISATION_HH
#define RENORMALISATION_HH RENORMALISATION_HH
#include <math.h>
#include "auxilliary/auxilliary.hh"
#include "auxilliary/parameters.hh"
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
   * @param[in] M_lat_ Number of time slices
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] renormalisation_ Type of renormalisation to use 
   *              (0: none, 1: perturbative, 2: exact)
   */
  RenormalisedParameters(const unsigned int M_lat_,
                         const double T_final_,
                         const RenormalisationType renormalisation_) :
    M_lat(M_lat_), T_final(T_final_), 
    a_lat(T_final_/M_lat_), renormalisation(renormalisation_) {}
  
protected:
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief Lattice spacing */
  const double a_lat;
  /** @brief Type of renormalisation */
  const RenormalisationType renormalisation;
};

#endif // RENORMALISATION_HH
