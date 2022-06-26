#ifndef ROTORCONDITIONEDFINEACTION_HH
#define ROTORCONDITIONEDFINEACTION_HH ROTORCONDITIONEDFINEACTION_HH
#include "action/action.hh"
#include "action/conditionedfineaction.hh"
#include "action/qm/rotoraction.hh"
#include "common/auxilliary.hh"
#include "common/samplestate.hh"
#include "distribution/expsin2distribution.hh"
#include "lattice/lattice1d.hh"
#include "mpi/mpi_random.hh"
#include <algorithm>
#include <memory>
#include <random>

/** @file rotorconditionedfineaction.hh
 * @brief Header file for rotor conditioned fine action class
 */

/** @class RotorConditionedFineAction
 *
 * @brief Conditioned fine action for the QM rotor
 *
 * The conditioned action is given by
 * \f[
 *    S^{cond}(X;X_-,X_+) =
 * 2W''(X_-,X_+)\sin^2\left(X-\overline{X}(X_-,X_+)\right) \f] where the
 * curvature \f$W''(X_-,X_+)\f$ and mean \f$\overline{X}(X_-,X_+)\f$ are given
 * by the getWcurvature() and getWminimum() methods of underlying action class
 * which is passed to the constructor.
 */
class RotorConditionedFineAction : public ConditionedFineAction {
public:
  /** @brief Constructor
   *
   * Construct new instance
   *
   * @param[in] action_ Underlying action class
   */
  RotorConditionedFineAction(const std::shared_ptr<RotorAction> action_)
      : action(action_) {
    engine.seed(11897197);
  }

  /** @brief Destructor */
  virtual ~RotorConditionedFineAction() {}

  /** @brief Fill in fine points
   *
   * Given a path \f$X\f$ for which the coarse points have been set, fill in
   * the fine points by sampling from the conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual void fill_fine_points(std::shared_ptr<SampleState> x_path) const;

  /** @brief Evaluate conditioned action at fine points
   *
   * Given a path \f$X\f$ for which all points have been set, evaluate the
   * conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual double evaluate(const std::shared_ptr<SampleState> x_path) const;

private:
  /** @brief Underlying action class */
  const std::shared_ptr<RotorAction> action;
  /** @brief Random number engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** Probability distribution */
  const ExpSin2Distribution exp_sin2_dist;
};

struct RotorConditionedFineActionFactory : public ConditionedFineActionFactory {
  /** @brief Destructor */
  virtual ~RotorConditionedFineActionFactory() {}

  /** @brief Method for constructing a conditioned fine action
   * @param[in] action Coarse level action
   */
  virtual std::shared_ptr<ConditionedFineAction>
  get(std::shared_ptr<Action> action) {
    return std::make_shared<RotorConditionedFineAction>(
        std::dynamic_pointer_cast<RotorAction>(action));
  }
};

#endif // ROTORCONDITIONEDFINEACTION_HH
