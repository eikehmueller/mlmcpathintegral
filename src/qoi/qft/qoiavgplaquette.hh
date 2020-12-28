#ifndef QOIAVGPLAQUETTE_HH
#define QOIAVGPLAQUETTE_HH QOIAVGPLAQUETTE_HH
#include <memory>
#include <cmath>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoiqvgplaquette.hh
 * @brief Header file for average plaquette in the 2D Schwinger model
 */

/** @class QoIAvgPlaquette
 *
 * @brief class for calculating the average value of the plaquette in the 2D Schwinger model
 *
 */

class QoIAvgPlaquette : public QoI {
public:
    /** @brief Create new instance
     *
     * @param[in] lattice_ Lattice
     */
    QoIAvgPlaquette(const std::shared_ptr<Lattice2D> lattice_) :
        lattice(lattice_),
        Mt_lat(lattice_->getMt_lat()),
        Mx_lat(lattice_->getMx_lat()) {}

    /** @brief Destructor */
    virtual ~QoIAvgPlaquette() {}

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

#endif // QOIAVGPLAQUETTE_HH
