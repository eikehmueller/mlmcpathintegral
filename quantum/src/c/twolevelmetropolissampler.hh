#ifndef TWOLEVELMETROPOLISSAMPLER_HH
#define TWOLEVELMETROPOLISSAMPLER_HH TWOLEVELMETROPOLISSAMPLER_HH

#include "action.hh"
#include "sampler.hh"

/** @file twolevelmetropolissampler.hh
 * @brief Header file for two level Metropolis sampler
 */

/** @class TwoLevelMetropolisSampler 
 *
 * @brief Two level sampler
 *
 * Given a way of creating coarse samples \f$\theta_{\ell-1}\f$ from a 
 * sampler on level \f$\ell-1\f$, create fine samples \f$\theta_{\ell}\f$
 * using a two level Metropolis Hastings process.
 * For this, states on level \f$\ell\f$ is written as 
 * \f$\theta_\ell = [\theta_{\ell,C},\theta_{\ell,F}]\f$ where
 * \f$\theta_{\ell,C}\f$ are the states on the even timeslices and 
 * \f$\theta_{\ell,F}\f$ are the states on the odd timeslices.
 *
 * To construct a new trial state, take the even modes from the next level
 * \f$\ell\f$ sample, \f$\theta'_{\ell,C}=\theta_{\ell}^{n+1}\f$. The 
 * odd modes are drawn from the free action, conditioned on the two
 * neighbouring odd modes \f$x_-\f$ and \f$x_+\f$, i.e. 
 *
 * \f[
   p(x) \sim \exp\left[-\frac{m_0}{a} \left(x-\frac{x_-+x_+}{2}\right)^2\right].
   \f]
 *
 * Finally, accept and reject according to \f$\min\left\{1,\beta\right\}\f$ with
 *
 * \f[
      \beta = \frac{\pi^{\ell}(\theta'_\ell)\pi^{\ell-1}(\theta^n_{\ell,C})\pi_{free}^{\ell}(\theta^n_{\ell,F}|\theta^n_{\ell,C})}{\pi^{\ell}(\theta^n_\ell)\pi^{\ell-1}(\theta'_{\ell,C})\pi_{free}^{\ell}(\theta'_{\ell,F}|\theta'_{\ell,C})}
 * \f]
 *
 * This guarantees that the fine level samples have the correct distribution.
 */
class TwoLevelMetropolisSampler {
public:
  /** @brief Create a new instance
   *
   * @param[in] coarse_sampler_ Sampler on coarse level \f$\ell-1\f$
   * @param[in] coarse_action_ Action on coarse level \f$\ell-1\f$
   * @param[in] fine_action_ Action on fine level \f$\ell\f$
   */
  TwoLevelMetropolisSampler(Sampler& coarse_sampler_,
                            const Action& coarse_action_,
                            const Action& fine_action_) :
    coarse_sampler(coarse_sampler_),
    coarse_action(coarse_action_),
    fine_action(fine_action_) {
    assert(2*coarse_action.getM_lat()==fine_action.getM_lat());
    assert(coarse_sampler.getM_lat()==coarse_action.getM_lat());
    theta_coarse = new Path(coarse_action.getM_lat());
    theta_fine = new Path(fine_action.getM_lat());
    theta_fine_C = new Path(coarse_action.getM_lat());
    theta_prime = new Path(fine_action.getM_lat());
    engine.seed(89216491);
  }

  /** @brief Destroy intance
   * 
   * Deallocate all temporary memory
   */
  ~TwoLevelMetropolisSampler() {
    delete theta_coarse;
    delete theta_fine;
    delete theta_fine_C;
    delete theta_prime;
  }  

  /** @brief draw new fine-/coarse-level path pair
   *
   * @param[out] x_path Vector of pointers to fine- and coarse- path
   */
  virtual void draw(std::vector<Path*> x_path);
  
private:
  /** @brief Evaluate conditioned free action
   *
   * Let
   * \f[
     S_{free}[\theta] = \frac{m_0}{a}\sum_{j=0}^{M/2}(\theta_{2j}-\frac{\theta_{2j}+\theta_{2j+2}}{2})^2
     \f]
   * 
   * This can be used to calculate the conditioned free probability as
   * \f$\pi_{free}^{\ell}(\theta_F|\theta_C) = e^{-S_{free}[\theta]}\f$
   *
   * which is required in the accept-/reject-step
   *
   * @param[in] x_path Path \f$\theta\f$ on which to evaluate
   * \f$\pi_{free}^{\ell}(\theta)\f$
   */
  const double conditioned_free_action(const Path* x_path);
                            
protected:
  /** @brief Sampler on coarse level */
  Sampler& coarse_sampler;
  /** @brief Action on coarse level */
  const Action& coarse_action;
  /** @brief Action on fine level */
  const Action& fine_action;
  /** @brief Temporary state vector on coarse level \f$\theta^n_{\ell-1}\f$ */
  mutable Path* theta_coarse;
  /** @brief Temporary state vector on fine level \f$\theta^n_{\ell}\f$ */
  mutable Path* theta_fine;
  /** @brief Coarse part of state vector on fine level \f$\theta^n_{\ell,C}\f$ */
  mutable Path* theta_fine_C;
  /** @brief Trial state vector on fine level \f$\theta'_{\ell}\f$ */
  mutable Path* theta_prime;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of normal distribution */
  typedef std::normal_distribution<double> Normal;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Normal distribution for drawing from conditioned free distribution */
  mutable Normal normal_dist;
  /** @brief Normal distribution in [0,1] for accept/reject step */
  mutable Uniform uniform_dist;

};
#endif // TWOLEVELMETROPOLISSAMPLER_HH
