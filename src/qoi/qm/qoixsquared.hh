#ifndef QOIXSQUARED_HH
#define QOIXSQUARED_HH QOIXSQUARED_HH
#include "action/qm/qmaction.hh"
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
#include "mpi/mpi_wrapper.hh"
#include "qoi/quantityofinterest.hh"
#include <memory>

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
  QoIXsquared(const std::shared_ptr<Lattice1D> lattice)
      : M_lat(lattice->getM_lat()) {}

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

/** @class QoIXsquaredFactory
 *
 * @brief Factory for constructing the QoI for a particular action
 */
class QoIXsquaredFactory : public QoIFactory {
public:
  /** @brief Return QoI for a specific  action
   *
   * @param[in] action Action to use
   */
  virtual std::shared_ptr<QoI> get(std::shared_ptr<Action> action) {
    std::shared_ptr<QMAction> qmaction =
        std::dynamic_pointer_cast<QMAction>(action);
    std::shared_ptr<Lattice1D> lattice = qmaction->get_lattice();
    return std::make_shared<QoIXsquared>(lattice);
  }
};

#endif // QOIXSQUARED_HH
