#ifndef CLUSTERACTION_HH
#define CLUSTERACTION_HH CLUSTERACTION_HH

#include "common/samplestate.hh"
#include "lattice/lattice.hh"

/** @file clusteraction.hh
 * @brief Header file for abstract cluster action class
 */

/** @class ClusterAction
 *
 * @brief Base class for QMcluster action
 *
 * Extends the action class by methods which are required to implement the
 * cluster algorithm. This assumes that the energy can be written as the
 * sum over links
 * \f[
 *   E = \sum_{\ell} E_{\ell}(x^{(\ell)}_-,x^{(\ell)}_+)
 * \f]
 * where a link \f$\ell=(x_i,x_j)\f$ connects two neighbouring sites.
 * This, of course, implies that \f$x^{(\ell)}_-=x_i\f$ and
 * \f$x^{(\ell)}_+=x_{j}\f$.
 */
class ClusterAction {
public:
  /** @brief Initialise class
   *
   * Create new instance of class.
   *
   * @param[in] generic_lattice_ Underlying generic lattice
   */
  ClusterAction(const std::shared_ptr<Lattice> lattice_)
      : generic_lattice(lattice_) {}

  /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
   *
   * The probability to have a bond between sites \f$i\f$ and \f$j\f$
   * is given by \f$1-e^{\min(0,-S_{\ell})}\f$
   *
   * @param[in] x_path Sample state
   * @param[in] i index of first vertex
   * @param[in] i index of second vertex
   */
  virtual double S_ell(const std::shared_ptr<SampleState> x_path,
                       const unsigned int i, const unsigned int j) const = 0;

  /** @brief Set reflection direction for the next step of the cluster algorithm
   */
  virtual void new_reflection() const = 0;

  /** @brief Flip the spin at a given site
   *
   * @param[inout] x_path State to process
   * @param[in] ell Vertex at which to flip the spin
   */
  virtual void flip(std::shared_ptr<SampleState> x_path,
                    const unsigned int ell) const = 0;

  /** @brief Initialise state */
  virtual void initialise_state(std::shared_ptr<SampleState> x_path) const = 0;

  /** @brief return size of samples */
  virtual unsigned int sample_size() const = 0;

  const std::shared_ptr<Lattice> get_generic_lattice() const {
    return generic_lattice;
  }

protected:
  const std::shared_ptr<Lattice> generic_lattice;
};

#endif // CLUSTERACTION_HH
