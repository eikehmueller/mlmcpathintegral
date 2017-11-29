#ifndef HARMONICOSCILLATORACTION_HH
#define HARMONICOSCILLATORACTION_HH HARMONICOSCILLATORACTION_HH
#include <cassert>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "path.hh"
#include "sampler.hh"
#include "action.hh"

/** @file harmonicoscillatoraction.hh
 * @brief Header file for harmonic oscillator action base class
 */


/** @class HarmonicOscillatorAction
 *
 * @brief Action class for harmonic oscillator 
 *
 * Action class for potential \f$V(x)=\frac{m_0}{2}\mu^2x^2\f$
 */
class HarmonicOscillatorAction : public Action, public Sampler {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] m0_ Mass of particle \f$m_0\f$
   * @param[in] mu2_ Frequency \f$\mu^2\f$
   */
  HarmonicOscillatorAction(const unsigned int M_lat_,
                           const double T_final_,
                           const double m0_,
                           const double mu2_)
    : Action(M_lat_,T_final_,m0_), Sampler(M_lat_), mu2(mu2_) {
    assert(mu2>0.0);
    build_covariance();
    engine.seed(124129017);
    y_tmp = new Path(M_lat_);
  }

  /** @brief Tidy up
   * 
   * Delete any temporary arrays
   */
  ~HarmonicOscillatorAction() {
    delete y_tmp;
  }
  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(const Path* x_path) const;

  /** @brief Draw sample from distribution
   *
   * Generate path which is distruted according to
   * \f$\pi(X) \sim e^{-S[X]}\f$
   *
   * @param[out] x_path: path to populate
   */
  const virtual void draw(std::vector<Path*> x_path);

  /** @brief Exact expression for expectation value of \f$X^2\f$
   *
   * Return the exact expression from Creutz and Freedman:
   *
   \f[
   \langle X^2 \rangle = \frac{1}{2m_0\mu\sqrt{1+\frac{a^2\mu^2}{4}}}\cdot \frac{1+R^M}{1-R^M}
   \f]
   * where
   \f[
   R = 1 + \frac{a^2\mu^2}{2}-a\mu\sqrt{1+\frac{a^2\mu^2}{4}}
   \f]
   */
  const double Xsquared_exact();

private:
  /** @brief Construct Cholesky factor of covariance matrix */
  void build_covariance();

private:
  /** @brief Oscillator frequency */
  const double mu2;
  /** @brief Eigen matrix type */
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Matrix;
  /** @brief Eigen vector type */
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> Vector;
  /** @brief Cholesky factor of the covariance matrix \f$\Sigma=L^TL\f$ */
  Matrix L_cov;
  /** @brief Type for Mersenne twister engine */
  typedef std::mt19937_64 Engine;
  /** @brief Random number engine */
  mutable Engine engine;
  /** @brief Type for normal distribution */
  typedef std::normal_distribution<double> Normal;
  /** @brief Normal distribution */
  mutable Normal normal_dist;
  /** @brief temporary array used for direct sampling */
  Path* y_tmp; 
};

#endif // HARMONICOSCILLATORACTION_HH
