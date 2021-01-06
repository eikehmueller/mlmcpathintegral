#ifndef QOI2DSUSCEPTIBILITY_HH
#define QOI2DSUSCEPTIBILITY_HH QOI2DSUSCEPTIBILITY_HH
#include <memory>
#include <cmath>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoi2dsusceptibility.hh
 * @brief Header file for topological susceptibility in the 2D Schwinger model
 */

/** @class QoI2DSusceptibility
 *
 * @brief class for calculating the (dimensionless) topological susceptibility in the 2D Schwinger model
 *
 * The definition used here is
 *
 * \f[
 *    \chi_t = \frac{a^2}{LT} \langle Q^2\rangle = \frac{\langle Q^2 \rangle}{M_t M_x}
 * \f]
 *
 * where \f$Q\f$ is the topological charge.
 *
 */

class QoI2DSusceptibility : public QoI {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Lattice
     */
    QoI2DSusceptibility(const std::shared_ptr<Lattice2D> lattice_) :
        lattice(lattice_),
        four_pi2_P_inv(0.25/(M_PI*M_PI*lattice_->getNcells())),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()) {}

    /** @brief Destructor */
    virtual ~QoI2DSusceptibility() {}

    /** @brief Evaluate on a state
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
     */
    const double virtual evaluate(const std::shared_ptr<SampleState> phi_state);

private:
    /** @brief Temporal extent of lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Scaling factor \f$1/(4\pi^2 P)\f$ */
    const double four_pi2_P_inv;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of lattice points in spatial direction */
    const unsigned int Mx_lat;
};

#endif // QOI2DSUSCEPTIBILITY_HH
