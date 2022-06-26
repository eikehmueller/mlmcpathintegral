#ifndef QOI2DPHISQUARED_HH
#define QOI2DPHISQUARED_HH QOI2PHISQUARED_HH
#include "action/qft/qftaction.hh"
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice2d.hh"
#include "mpi/mpi_wrapper.hh"
#include "qoi/quantityofinterest.hh"
#include <cmath>
#include <memory>

/** @file qoi2dphisquared.hh
 * @brief Header file for averade squared value of 2d scalar field
 */

/** @class QoI2DPhiSquared
 *
 * @brief class for calculating the average squared value of a 2d scalar field
 *
 * The returned quantity is defined as
 *
 * \f[
 *    QoI = \frac{1}{M}\sum_{j=0}^{M-1} \phi_j^2
 * \f]
 *
 * where \f$M\f$ is the number of lattice sites.
 *
 */

class QoI2DPhiSquared : public QoI {
public:
  /** @brief Create new instance
   *
   * @param[in] lattice_ Lattice
   */
  QoI2DPhiSquared(const std::shared_ptr<Lattice2D> lattice_)
      : lattice(lattice_), M_lat(lattice_->getNvertices()) {}

  /** @brief Destructor */
  virtual ~QoI2DPhiSquared() {}

  /** @brief Evaluate on a state
   *
   * @param[in] phi_state State \f$\phi\f$ on which to evaluate the QoI
   */
  const double virtual evaluate(const std::shared_ptr<SampleState> phi_state);

private:
  /** @brief Temporal extent of lattice */
  const std::shared_ptr<Lattice2D> lattice;
  /** @brief Number of vertices */
  const unsigned int M_lat;
};

/** @class QoI2DPhiSquaredFactory
 *
 * @brief Factory for constructing the QoI for a particular action
 */
class QoI2DPhiSquaredFactory : public QoIFactory {
public:
  /** @brief Return QoI for a specific  action
   *
   * @param[in] action Action to use
   */
  virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) {
    std::shared_ptr<Lattice2D> lattice =
        std::dynamic_pointer_cast<QFTAction>(action)->get_lattice();
    return std::make_shared<QoI2DPhiSquared>(lattice);
  }
};

#endif // QOI2DPHISQUARED_HH
