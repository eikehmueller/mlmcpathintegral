#ifndef QUANTITYOFINTEREST_HH
#define QUANTITYOFINTEREST_HH QUANTITYOFINTEREST_HH
#include <memory>
#include <cmath>
#include "path.hh"
#include "auxilliary.hh"

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
  const double virtual evaluate(const std::shared_ptr<Path> x_path) = 0;
  
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
  QoIXsquared() {}

  /** @brief Destructor */
  virtual ~QoIXsquared() {}
  
  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<Path> x_path);
};

/** @class QoISusceptibility
 *
 * @brief class for calculating the susceptibility for the QM rotor
 *
 * For a path \f$X\f$ the susceptibility is defined as
 * \f$\chi_t = Q[X]^2/T \f$ where
 * \f[
 *   Q[X] = \frac{1}{2\pi} \sum_{j=0}^{M_{lat}-1} (X_j - X_{j-1}) \mod [-\pi,\pi]
 * \f]
 * and the time is \f$T=a_{lat}M_{lat}\f$.
 */

class QoISusceptibility : public QoI {
public: 
  /** @brief Create new instance 
   *
   * @param[in] M_lat_ Number of time slices \f$M\f$
   */
  QoISusceptibility() : four_pi2_inv(0.25/(M_PI*M_PI)) {}

  /** @brief Destructor */
  virtual ~QoISusceptibility() {}

  /** @brief Evaluate on a path
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<Path> x_path);
  
private:
  /** @brief Constant \f$1/(4\pi^2)\f$ */
  const double four_pi2_inv;
};

#endif // QUANTITYOFINTEREST_HH
