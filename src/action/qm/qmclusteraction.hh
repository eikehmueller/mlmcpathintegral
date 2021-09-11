#ifndef QMCLUSTERACTION_HH
#define QMCLUSTERACTION_HH QMCLUSTERACTION_HH

#include "lattice/lattice1d.hh"
#include "action/qm/qmaction.hh"

/** @file qmclusteraction.hh
 * @brief Header file for abstract quantum mechanical cluster action class
 */

/** @class QMClusterAction
 *
 * @brief Base class for QM cluster action
 *
 * Extends the action class by methods which are required to implement the
 * cluster algorithm for one-dimensional quantum problem with periodic
 * boundary conditions. This assumes that the energy can be written as the
 * sum over links
 * \f[
 *   E = \sum_{\ell} E_{\ell}(x^{(\ell)}_-,x^{(\ell)}_+)
 * \f]
 * where a link \f$\ell=(i,i+1)\f$ connects two neighbouring sites.
 * This, of course, implies that \f$x^{(\ell)}_-=x_i\f$ and
 * \f$x^{(\ell)}_+=x_{i+1}\f$.
 */
class QMClusterAction {
public:
    /** @brief Initialise class
     *
     * Create new instance of class.
     *
     * @param[in] generic_lattice_ Underlying generic lattice
     */
    QMClusterAction(const std::shared_ptr<Lattice> lattice_)
        : generic_lattice(lattice_) {}

    /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
     *
     * The probability to have a bond between sites \f$i\f$ and \f$i+1\f$
     * is given by \f$1-e^{\min(0,-S_{\ell})}\f$
     *
     * @param[in] x_m Value of \f$x^{(\ell)}_- = x_i\f$
     * @param[in] x_p Value of \f$x^{(\ell)}_+ = x_{i+1}\f$
     */
    virtual double S_ell(const double x_m, const double x_p) const = 0;

    /** @brief Set angle for the next step of the cluster algorithm
     */
    virtual void new_angle() const = 0;

    /** @brief Flip the current spin
     *
     * Return \f$hx\f$
     *
     * @param[in] x value of site \f$x\f$
     */
    virtual double flip(const double x) const = 0;

    /** @brief Initialise state */
    virtual void initialise_state(std::shared_ptr<SampleState> x_path) const = 0;
  
    /** @brief return size of samples */
    virtual unsigned int sample_size() const = 0;  

    const std::shared_ptr<Lattice> get_generic_lattice() const { return generic_lattice; }
  
protected:
    const std::shared_ptr<Lattice> generic_lattice;
      
};

#endif // QMCLUSTERACTION_HH
