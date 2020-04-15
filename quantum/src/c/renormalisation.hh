#ifndef RENORMALISATION_HH
#define RENORMALISATION_HH RENORMALISATION_HH
#include <math.h>
#include "auxilliary.hh"
#include "parameters.hh"
#include "mpi_wrapper.hh"

/** @file renormalisation.hh
 * @brief Methods for calculation renormalised parameters
 */

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

/** @class RenormalisedHOParameters 
 * @brief Renormalised coarse grid parameters for the harmonic oscillator
 * 
 * Calculate the renormalised coarse grid mass and harmonic oscillator
 * potential parameter \f$\mu^2\f$
 *
 \f[
  m_0^{(c)} = m_0\left(1+\frac{a^2\mu^2}{2}\right)^{-1} 
 \f]
 * 
 \f[
  \left(\mu^{(c)}\right)^2 = \mu^2\left(1+\frac{a^2\mu^2}{4}\right)^{-1} 
 \f]
 *
 * Instead of those exact formulae, the perturbative expansion in \f$\mu\f$
 * can also be used.
 */
class RenormalisedHOParameters : public RenormalisedParameters {
public:
  /** @brief Create new instance
   *
   * @param[in] M_lat_ Number of time slices
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass \f$m_0\f$
   * @param[in] mu2_ Harmonic oscillator potential parameter \f$\mu^2\f$
   * @param[in] renormalisation_ Type of renormalisation to use 
   *              (0: none, 1: perturbative, 2: exact)
   */
  RenormalisedHOParameters(const unsigned int M_lat_,
                           const double T_final_,
                           const double m0_,
                           const double mu2_,
                           const RenormalisationType renormalisation_) :
    RenormalisedParameters(M_lat_, T_final_, renormalisation_),
    m0(m0_), mu2(mu2_) {}

  /** @brief Renormalised coarse level mass \f$m_0^{(c)}\f$*/
  double m0_coarse() {
    double m0coarse;
    switch (renormalisation) {
    case RenormalisationNone:
      m0coarse = m0;
      break;
    case RenormalisationPerturbative:
      m0coarse = m0*(1.-0.5*a_lat*a_lat*mu2);
      break;
    case RenormalisationExact:
      m0coarse = m0/(1.+0.5*a_lat*a_lat*mu2);
      break;
    }
    return m0coarse;
  }
  /** @brief Renormalised coarse level oscillator parameter
   * \f$\left(\mu^{(c)}\right)^2\f$ */
  double mu2_coarse() {
    double mu2coarse;
    switch (renormalisation) {
    case RenormalisationNone:
      mu2coarse = mu2;
      break;
    case RenormalisationPerturbative:
      mu2coarse = mu2*(1.+0.25*a_lat*a_lat*mu2); 
      break;
    case RenormalisationExact:
      mu2coarse = mu2*(1.+0.25*a_lat*a_lat*mu2);
      break;
    }
    return mu2coarse;
  }
  
private:
  /** @brief Mass \f$m_0\f$ */
  const double m0;
  /** @brief Harmonic oscillator potential parameter \f$\mu^2\f$*/
  const double mu2;
};

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

#endif // RENORMALISATION_HH
