#ifndef QUANTITYOFINTEREST_HH
#define QUANTITYOFINTEREST_HH QUANTITYOFINTEREST_HH
#include <memory>
#include "common/samplestate.hh"
#include "action/action.hh"

/** @file quantityofinterest.hh
 * @brief Header file for quantities of interest
 */

/** @class QoI
 *
 * @brief Abstract base class for state-dependent quantity of interest
 *
 */
class QoI {
public:
    /** @brief Create new instance
     *
     */
    QoI() {}

    /** @brief Evaluate on a state
     *
     * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
     */
    const double virtual evaluate(const std::shared_ptr<SampleState> phi_state) = 0;

};

/** Base class for QoI factory */
class QoIFactory {
public:
    /** @brief extract a qoi for a given action */
    virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) = 0;
};


#endif // QUANTITYOFINTEREST_HH
