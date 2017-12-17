#ifndef RENORMALISATION_HH
#define RENORMALISATION_HH RENORMALISATION_HH
/** @file renormalisation.hh
 * @brief Methods for calculation renormalised parameters
 */

/** @class RenormalisedParameters 
 * @brief Renormalised coarse grid parameters
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
class RenormalisedHOParameters {
public:
  /** @brief Create new instance
   *
   * @param[in] M_lat_ Number of time slices
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass \f$m_0\f$
   * @param[in] mu2_ Harmonic oscillator potential parameter \f$\mu^2\f$
   * @param[in] perturbative_ Use perturbative expansion?
   */
  RenormalisedHOParameters(const unsigned int M_lat_,
                           const double T_final_,
                           const double m0_,
                           const double mu2_,
                           const bool perturbative_) :
    m0(m0_), mu2(mu2_), T_final(T_final_), M_lat(M_lat_),
    a_lat(T_final_/M_lat_), perturbative(perturbative_) {}
  /** @brief Renormalised coarse level mass \f$m_0^{(c)}\f$*/
  double m0_coarse() {
    if (perturbative) {
      return m0*(1.-0.5*a_lat*a_lat*mu2);
    } else {
      return m0/(1.+0.5*a_lat*a_lat*mu2);
    }
  }
  /** @brief Renormalised coarse level oscillator parameter
   * \f$\left(\mu^{(c)}\right)^2\f$ */
  double mu2_coarse() {
    return mu2*(1.+0.25*a_lat*a_lat*mu2);
  }
  
private:
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief Mass \f$m_0\f$ */
  const double m0;
  /** @brief Harmonic oscillator potential parameter \f$\mu^2\f$*/
  const double mu2;
  /** @brief Lattice spacing */
  const double a_lat;
  /** @brief Use perturbative expansion? */
  const bool perturbative;
};

#endif // RENORMALISATION_HH
