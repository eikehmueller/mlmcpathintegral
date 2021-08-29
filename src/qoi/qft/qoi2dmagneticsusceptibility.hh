#ifndef QOI2DMAGNETICSUSCEPTIBILITY_HH
#define QOI2DMAGNETICSUSCEPTIBILITY_HH QOI2DMAGNETICSUSCEPTIBILITY_HH
#include <memory>
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/qft/nonlinearsigmaaction.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoi2dmagneticsusceptibility.hh
 * @brief Header file for average squared magnetisation of non-linear O(3) sigma model
 */

/** @class QoI2DMagneticSusceptibility
 *
 * @brief class for calculating the magnetic susceptibility of the non-linear O(3) sigma model
 *
 * The returned magnetic susceptibility (average squared magnetisation)  is defined as
 *
 * \f[
 *    \chi_t = 1/M*\langle \mu^2 \rangle
 * \f]
 *
 * where the 3-vector \f$\mu\f$ is the sum of all spins on the lattice and \f$M\f$ is the number
 * of lattice sites.
 *
 */

class QoI2DMagneticSusceptibility : public QoI {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Lattice
     * @param[in] rotated_ is the lattice action rotated?
     */
    QoI2DMagneticSusceptibility(const std::shared_ptr<Lattice2D> lattice_,
                       const bool rotated_) :
        lattice(lattice_),
        rotated(rotated_),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()) {}

    /** @brief Destructor */
    virtual ~QoI2DMagneticSusceptibility() {}

    /** @brief Evaluate on a state
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
     */
    const double virtual evaluate(const std::shared_ptr<SampleState> phi_state);

private:
    /** @brief Temporal extent of lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Is the lattice action rotated? */
    const bool rotated;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of lattice points in spatial direction */
    const unsigned int Mx_lat;
};

/** @class QoI2DMagneticSusceptibilityFactory
 *
 * @brief Factory for constructing the QoI for a particular action
 */
class QoI2DMagneticSusceptibilityFactory : public QoIFactory {
public:
    /** @brief Return QoI for a specific  action
     *
     * @param[in] action Action to use
     */
    virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) {
        std::shared_ptr<NonlinearSigmaAction> nonlinear_sigma_action;
        nonlinear_sigma_action = std::dynamic_pointer_cast<NonlinearSigmaAction>(action);
        std::shared_ptr<Lattice2D> lattice = nonlinear_sigma_action->get_lattice();
        bool rotated = nonlinear_sigma_action->is_rotated();
        return std::make_shared<QoI2DMagneticSusceptibility>(lattice,rotated);
    }
};

#endif // QOI2DMAGNETICSUSCEPTIBILITY_HH