#ifndef QUANTITYOFINTEREST_HH
#define QUANTITYOFINTEREST_HH QUANTITYOFINTEREST_HH
#include "path.hh"

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
   * @param[in] M_lat_ Number of time slices \f$M\f$
   */
  QoI(unsigned int M_lat_) : M_lat(M_lat_) {}

  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const Path* x_path) = 0;
  
protected:
  /** @brief Number of time slices */
  const unsigned int M_lat;
};

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
  QoIXsquared(unsigned int M_lat_) : QoI(M_lat_) {}

  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const Path* x_path);
};


#endif // QUANTITYOFINTEREST_HH
