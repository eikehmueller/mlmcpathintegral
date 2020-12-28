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
 * @brief class for calculating the topological susceptibility in the 2D Schwinger model
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
        four_pi2_inv(0.25/(M_PI*M_PI)),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()),
        vol_lat(lattice_->getT_lat()*lattice_->getL_lat()) {}

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
    /** @brief Scaling factor \f$1/(4\pi^2)\f$ */
    const double four_pi2_inv;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of lattice points in spatial direction */
    const unsigned int Mx_lat;
    /** @brief Lattice volume */
    const double vol_lat;
};

#endif // QOI2DSUSCEPTIBILITY_HH
