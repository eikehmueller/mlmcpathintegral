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
 * - Given a path \f$X\f$ of length \f$M\f$, for which only the coarse points
 *   \f$X_{2i}\f$ with \f$i=0,\dots,M/2-1\f$ have been set,
 *   fill in all fine points \f$X_{2i+1}\f$ by sampling from a suitable
 *   conditioned probability distribution
 *   \f[
 *     p(X_{2i+1}|X_{2i},X_{2i+2}) = Z(X_{2i,2i+1})^{-1} \exp\left[-S^{cond}(X_{2i+1;X_{2i},X_{2i+2}}) \right] \qquad\text{for all $i=0,\dots,M/2-1$}.
 *   \f]
 *   Note that \f$Z(2i,2i+2)\f$ is a normalisation constant which guarantees
 *   that this is indeed a probability density.
 * 
 * - Given a path \f$X\f$ of length \f$M\f$ for which the fine points
 *   \f$X_{2i+1}\f$ with \f$i=0,\dots,M/2-1\f$ have been set (for example with 
 *   the above method), calculate the value of the conditioned action
 *   (including the normalisation constant) as
 *   \f[
 *     \sum_{i=0}^{M/2-1} S^{cond}(X_{2i+1;X_{2i},X_{2i+2}}) + \log(Z(X_{2i},X_{2i+2}) 
 *   \f]
 */

class ConditionedFineAction {
public:
  /** @brief Fill in fine points
   * 
   * Given a path \f$X\f$ for which the coarse points have been set, fill in
   * the fine points by sampling from the conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual void fill_fine_points(std::shared_ptr<SampleState> x_path) const = 0;

  /** @brief Evaluate conditioned action at fine points
   * 
   * Given a path \f$X\f$ for which all points have been set, evaluate the
   * conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual double evaluate(const std::shared_ptr<SampleState> x_path) const = 0;  
};

struct ConditionedFineActionFactory {
  /** @brief Abstract method for constructing a conditioned fine action
   * @param[in] action Coarse level action
   */
  virtual std::shared_ptr<ConditionedFineAction> get(std::shared_ptr<Action> action) = 0;
};

#endif // CONDITIONEDFINEACTION_HH
