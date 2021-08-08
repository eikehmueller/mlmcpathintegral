#ifndef QOI2DMAGNETISATION_HH
#define QOI2DMAGNETISATION_HH QOI2DMAGNETISATION_HH
#include <memory>
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "lattice/lattice2d.hh"
#include "action/qft/qftaction.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoi2dmagnetisation.hh
 * @brief Header file for average squared magnetisation of non-linear O(3) sigma model
 */

/** @class QoI2DMagnetisation
 *
 * @brief class for calculating the average squared magnetisation of the non-linear O(3) sigma model
 *
 * The returned magnetisation is defined as
 *
 * \f[
 *    V\chi_t = 1/M*\langle \mu^2 \rangle
 * \f]
 *
 * where the 3-vector \f$\mu\f$ is the sum of all spins on the lattice and \f$M\f$ is the number
 * of lattice sites.
 *
 */

class QoI2DMagnetisation : public QoI {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Lattice
     */
    QoI2DMagnetisation(const std::shared_ptr<Lattice2D> lattice_) :
        lattice(lattice_),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()) {}

    /** @brief Destructor */
    virtual ~QoI2DMagnetisation() {}

    /** @brief Evaluate on a state
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
     */
    const double virtual evaluate(const std::shared_ptr<SampleState> phi_state);

private:
    /** @brief Temporal extent of lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of lattice points in spatial direction */
    const unsigned int Mx_lat;
};

/** @class QoI2DMagnetisationFactory
 *
 * @brief Factory for constructing the QoI for a particular action
 */
class QoI2DMagnetisationFactory : public QoIFactory {
public:
    /** @brief Return QoI for a specific  action
     *
     * @param[in] action Action to use
     */
    virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) {
        std::shared_ptr<Lattice2D> lattice = std::dynamic_pointer_cast<QFTAction>(action)->get_lattice();
        return std::make_shared<QoI2DMagnetisation>(lattice);
    }
};

#endif // QOI2DMAGNETISATION_HH
