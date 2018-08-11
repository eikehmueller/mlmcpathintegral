#ifndef CONDITIONEDFINEACTION_HH
#define CONDITIONEDFINEACTION_HH CONDITIONEDFINEACTION_HH
#include <random>
#include "path.hh"
#include "action.hh"

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
 *     p(X_{2i+1}|X_{2i},X_{2i+2}) = Z(X_{2i,2i+1}) \exp\left[-S^{cond}_(X_{2i+1;X_{2i},X_{2i+2}}) \right] \qquad\text{for all $i=0,\dots,M/2-1$}.
 *   \f]
 *   Note that \f$Z(2i,2i+2)\f$ is a normalisation constant which guarantees
 *   that this is indeed a probability density.
 * 
 * - Given a path \f$X\f$ of length \f$M\f$ for which the fine points
 *   \f$X_{2i+1}\f$ with \f$i=0,\dots,M/2-1\f$ have been set (for example with 
 *   the above method), calculate the value of the conditioned action
 *   \f[
 *     \sum_{i=0}^{M/2-1} \log(Z(X_{2i},X_{2i+2})) - S^{cond}(X_{2i+1;X_{2i},X_{2i+2}})
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
  virtual void fill_fine_points(Path* x_path) const = 0;

  /** @brief Evaluate conditioned action at fine points
   * 
   * Given a path \f$X\f$ for which all points have been set, evaluate the
   * conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual double evaluate(const Path* x_path) const = 0;  

  /** @brief Set external engine (debug only, to ensure bitreproducibility) */
  virtual void set_ext_engine(std::mt19937_64* ext_engine_) const = 0;
};

/** @class GaussianConditionedFineAction
 *
 * @brief Gaussian conditioned fine action
 *
 * The conditioned action is given by
 * \f[
 *    S^{cond}(X;X_-,X_+) = \frac{1}{2}W''(X_-,X_+)\left(X-\overline{X}(X_-,X_+)\right)^2
 * \f]
 * where the curvature \f$W''(X_-,X_+)\f$ and mean \f$\overline{X}(X_-,X_+)\f$
 * are given by the getWcurvature() and getWminimum() methods of underlying
 * action class which is passed to the constructor.
 */
class GaussianConditionedFineAction : public ConditionedFineAction {
public:
  /** @brief Constructor
   * 
   * Construct new instance
   *
   * @param[in] action_ Underlying action class
   */
  GaussianConditionedFineAction (const Action& action_) : action(action_) {
    engine.seed(11897197);
  }

  /** @brief Fill in fine points
   * 
   * Given a path \f$X\f$ for which the coarse points have been set, fill in
   * the fine points by sampling from the conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual void fill_fine_points(Path* x_path) const;

  /** @brief Evaluate conditioned action at fine points
   * 
   * Given a path \f$X\f$ for which all points have been set, evaluate the
   * conditioned action.
   *
   * @param[inout] x_path Path \f$X\f$ to fill
   */
  virtual double evaluate(const Path* x_path) const;

  /** @brief Set external engine (debug only, to ensure bitreproducibility) */
  virtual void set_ext_engine(std::mt19937_64* ext_engine_) const {
    ext_engine = ext_engine_;
  }
  
private:
  /** @brief Underlying action class */
  const Action& action;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** External engine (for debugging) */
  mutable Engine* ext_engine;
  /** @brief Type of normal distribution */
  typedef std::normal_distribution<double> Normal;
  /** @brief Normal distribution for drawing from distribution */
  mutable Normal normal_dist;

};

#endif // CONDITIONEDFINEACTION_HH
