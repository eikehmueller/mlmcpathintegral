#ifndef QUANTITYOFINTEREST_HH
#define QUANTITYOFINTEREST_HH QUANTITYOFINTEREST_HH
/** @class QoI
 *
 * @brief Abstract base class for path-dependent quantity of interest
 *
 */

class QoI {
public: 
  /** @brief Create new instance 
   *
   * @param M_lat Number of time slices \f$M\f$
   */
  QoI(unsigned int M_lat_) : M_lat(M_lat_) {}

  /** @brief Evaluate on a path
   *
   * @param x_path Path \f$\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const double* x_path) = 0;
  
protected:
  const unsigned int M_lat;
};

/** @class QoIXsuared
 *
 * @brief Class for calculating \f$X_0^2\f$
 *
 */

class QoIXsquared : public QoI {
public: 
  /** @brief Create new instance 
   *
   * @param M_lat Number of time slices \f$M\f$
   */
  QoIXsquared(unsigned int M_lat_) : QoI(M_lat_) {}

  /** @brief Evaluate on a path
   *
   * @param x_path Path \f$\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const double* x_path);
};


#endif // QUANTITYOFINTEREST_HH
