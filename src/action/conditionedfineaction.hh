#ifndef CONDITIONEDFINEACTION_HH
#define CONDITIONEDFINEACTION_HH CONDITIONEDFINEACTION_HH
#include <random>
#include <memory>
#include <algorithm>
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_random.hh"
#include "action/action.hh"

/** @file conditionedfineaction.hh
 * @brief Header file for conditioned fine action classes
 */

/** @class ConditionedFineAction
 *
 * @brief Base class for conditioned fine actions
 *
 * Derived classes of this type provide the following functionality:
 * 
 * - Given a state \f$\phi\f$ of length \f$M\f$, for which only the
 *   coarse points have been set, fill in all fine points by sampling from
 *   a suitable conditioned probability distribution
 *   \f[
 *     p(\phi_{fine}|\phi_{coarse}) = Z(\phi_{coarse})^{-1} \exp\left[-S^{cond}(\phi_{fine};\phi_{coarse}) \right].
 *   \f]
 *   Note that \f$Z(\phi_{coarse})\f$ is a normalisation constant which guarantees
 *   that this is indeed a probability density.
 * 
 * - Given a state \f$\phi\f$ for which the fine points have been set
 *   (for example with the above method), calculate the value of the
 *   conditioned action (including the normalisation constant) as
 *   \f[
 *     S^{cond}(\phi_{fine};\phi_{coarse}) + \log(Z(\phi_{coarse}))
 *   \f]
 */

class ConditionedFineAction {
public:
  /** @brief Fill in fine points
   * 
   * Given a phi \f$\phi\f$ for which the coarse points have been set,
   * fill in the fine points by sampling from the conditioned action.
   *
   * @param[inout] phi_state State \f$\phi\f$ to fill
   */
  virtual void fill_fine_points(std::shared_ptr<SampleState> phi_state) const = 0;

  /** @brief Evaluate conditioned action at fine points
   * 
   * Given a state \f$\phi\f$ for which all points have been set,
   * evaluate the conditioned action.
   *
   * @param[inout] phi_state State \f$\phi\f$ to fill
   */
  virtual double evaluate(const std::shared_ptr<SampleState> phi_state) const = 0;  
};

struct ConditionedFineActionFactory {
  /** @brief Abstract method for constructing a conditioned fine action
   * @param[in] action Coarse level action
   */
  virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) = 0;
};

#endif // CONDITIONEDFINEACTION_HH
