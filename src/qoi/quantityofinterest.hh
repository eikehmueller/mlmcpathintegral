#ifndef QUANTITYOFINTEREST_HH
#define QUANTITYOFINTEREST_HH QUANTITYOFINTEREST_HH
#include <memory>
#include "common/samplestate.hh"

/** @file quantityofinterest.hh
 * @brief Header file for quantities of interest
 */

/** @class QoI
 *
 * @brief Abstract base class for path-dependent quantity of interest
 *
 */
class QoI {
public: 
  /** @brief Create new instance 
   *
   */
  QoI() {}

  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<SampleState> x_path) = 0;
  
};

#endif // QUANTITYOFINTEREST_HH
