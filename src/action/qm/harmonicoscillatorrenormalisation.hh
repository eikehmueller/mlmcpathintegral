#ifndef HARMONICOSCILLATORRENORMALISATION_HH
#define HARMONICOSCILLATORRENORMALISATION_HH HARMONICOSCILLATORRENORMALISATION_HH
#include <memory>
#include "lattice/lattice1d.hh"
#include "action/renormalisation.hh"

/** @file harmonicoscillatorrenormalisation.hh
 * @brief Header file for harmonic oscillator action renormalisation class
 */

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
     * @param[in] lattice_ Underlying lattice
     * @param[in] m0_ Mass \f$m_0\f$
     * @param[in] mu2_ Harmonic oscillator potential parameter \f$\mu^2\f$
     * @param[in] renormalisation_ Type of renormalisation to use
     *              (0: none, 1: perturbative, 2: exact)
     */
    RenormalisedHOParameters(const std::shared_ptr<Lattice1D> lattice_,
                             const double m0_,
                             const double mu2_,
                             const RenormalisationType renormalisation_) :
        RenormalisedParameters(renormalisation_),
        lattice(lattice_), m0(m0_), mu2(mu2_) {}

    /** @brief Renormalised coarse level mass \f$m_0^{(c)}\f$*/
    double m0_coarse() {
        double a_lat = lattice->geta_lat();
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
        double a_lat = lattice->geta_lat();
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
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice1D> lattice;
};

#endif // HARMONICOSCILLATORRENORMALISATION_HH
