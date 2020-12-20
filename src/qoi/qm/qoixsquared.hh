#ifndef QOIXSQUARED_HH
#define QOIXSQUARED_HH QOIXSQUARED_HH
#include <memory>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
#include "qoi/quantityofinterest.hh"

/** @file qoixsquared.hh
 * @brief Header file for QoI X^2
 */

/** @class QoIXsquared
 *
 * @brief Class for calculating \f$X_0^2\f$
 *
 */

class QoIXsquared : public QoI {
public:
    /** @brief Create new instance  */
    QoIXsquared(const std::shared_ptr<Lattice1D> lattice) :
        M_lat(lattice->getM_lat()) {}

    /** @brief Destructor */
    virtual ~QoIXsquared() {}

    /** @brief Evaluate on a path
     *
     * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
     */
    const double virtual evaluate(const std::shared_ptr<SampleState> x_path);
private:
    /** @brief Number of lattice points */
    const unsigned int M_lat;
};

#endif // QOIXSQUARED_HH
