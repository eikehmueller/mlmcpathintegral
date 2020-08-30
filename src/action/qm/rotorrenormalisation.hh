#ifndef ROTORRENORMALISATION_HH
#define ROTORRENORMALISATION_HH ROTORRENORMALISATION_HH
#include "config.h"
#include <memory>
#include "action/renormalisation.hh"

/** @file rotoraction.hh
 * @brief Header file for renormalisation of quantum mechanical rotor class
 */

/** @class RenormalisedRotorParameters 
 * @brief Renormalised coarse grid parameters for the topological rotor
 * 
 * Calculate the renormalised coarse grid mass (= moment of inertia) 
 *
 \f[
  m_0^{(c)} = \left(1+\delta_I(T/m_0)\frac{a}{m_0}\right)m_0
 \f]
 * \f$\delta_I(T/m_0)\f$ depends on the ratio \f$T/m_0\f$. 
 */
class RenormalisedRotorParameters : public RenormalisedParameters {
public:
  /** @brief Create new instance
   *
   * @param[in] M_lat_ Number of time slices
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass (i.e. moment of inertia) \f$m_0\f$
   * @param[in] renormalisation_ Type of renormalisation to use 
   *              (0: none, 1: perturbative, 2: exact [not implemented])
   */
  RenormalisedRotorParameters(const unsigned int M_lat_,
                              const double T_final_,
                              const double m0_,
                              const RenormalisationType renormalisation_) :
    RenormalisedParameters(M_lat_, T_final_, renormalisation_),
    m0(m0_) {}
  
  /** @brief Renormalised coarse level mass \f$m_0^{(c)}\f$*/
  double m0_coarse() {
    double m0coarse;
    switch (renormalisation) {
    case RenormalisationNone:
      m0coarse = m0;
      break;
    case RenormalisationPerturbative:
      m0coarse = (1.+deltaI(T_final/m0)*a_lat/m0)*m0;
      break;
    case RenormalisationExact:
      m0coarse = 1.0;
      mpi_parallel::cerr << "ERROR: exact renormalisation not implemented for rotor action " << std::endl;
      mpi_exit(EXIT_FAILURE);
      break;
    }
    return m0coarse;
  }

  /** @brief Function \f$\delta (\xi)\f$
   *
   * @param[in] xi Argument \f$\xi\f$of the function
   */
  double deltaI(const double xi);

private:
  /** @brief Mass (moment of inertia) \f$m_0\f$ */
  const double m0;
};

#endif // ROTORRENORMALISATION_HH
