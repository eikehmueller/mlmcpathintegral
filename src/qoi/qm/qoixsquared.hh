#ifndef QOIXSQUARED_HH
#define QOIXSQUARED_HH QOIXSQUARED_HH
#include <memory>
#include "fields/path.hh"
#include "common/auxilliary.hh"
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
  /** @brief Create new instance 
   *
   * @param[in] M_lat_ Number of time slices \f$M\f$
   */
  QoIXsquared() {}

  /** @brief Destructor */
  virtual ~QoIXsquared() {}
  
  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<Path> x_path);
};

#endif // QOIXSQUARED_HH
