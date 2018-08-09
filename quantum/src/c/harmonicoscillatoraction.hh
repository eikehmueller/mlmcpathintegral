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
    : Action(M_lat_,T_final_,m0_), Sampler(false), mu2(mu2_),
      Wcurvature((2./a_lat + a_lat*mu2)*m0),
      Wminimum_scaling(0.5/(1.+0.5*a_lat*a_lat*mu2)) {
    assert(mu2>0.0);
    build_covariance();
    engine.seed(124129017);
    y_tmp = new Path(M_lat_,T_final_);
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

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path \f$P\f$.
   *
   * Note that for this action we have
     \f[
         P_j = \frac{m_0}{a}\left(2X_j-X_{j-1}-X_{j+1}\right) + \frac{1}{2}am_0\mu^2 X_j
     \f]
   *
   * @param x_path Path \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const Path* x_path,
                     Path* p_path) const;

  /** @brief Initialise path 
   *
   * Set initial values of path to zero
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_path(Path* x_path) const {
    std::fill(x_path->data,x_path->data+M_lat,0.0);
  }  
  
  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned action
   *
   * For the harmonic oscillator potential the curvature of the modified
   * action (see Action::getWcurvature()) is 
   \f[
   W''_{x_-,x_+} = \frac{2m_0}{a}+am_0\mu^2
   \f]
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const {
    return Wcurvature;
  }

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * For the harmonic oscillator potential the minimum of the modified
   * action (see Action::getWminimum()) can be found at
   \f[
      x_0 = \left(1+\frac{1}{2}a^2\mu^2\right)^{-1}\overline{x}
   \f]
   * where \f$\overline{x}=\frac{x_++x_-}{2}\f$.
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const {
    return Wminimum_scaling*(x_m+x_p);
  }
    
  /** @brief Draw sample from distribution
   *
   * Generate path which is distruted according to
   * \f$\pi(X) \sim e^{-S[X]}\f$
   *
   * @param[out] x_path: path to populate
   */
  virtual void draw(std::vector<Path*> x_path);

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

  /** @brief Continuum limitof expectation value of \f$X^2\f$
   *
   * Return the continuum limit of \f$\langle X^2\rangle\f$ from
   * Creutz and Freedman:
   *
   \f[
   \lim_{a\rightarrow 0}\langle X^2 \rangle = \frac{1}{2m_0\mu}\cdot \frac{1-e^{-\mu T}}{1+e^{-\mu T}}
   \f]
   * where
   \f[
   R = 1 + \frac{a^2\mu^2}{2}-a\mu\sqrt{1+\frac{a^2\mu^2}{4}}
   \f]
   */
  const double Xsquared_exact_continuum();

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
  /** @brief Curvature of modified potential */
  const double Wcurvature;
  /** @brief Scaling factor for calculation of mimimum of modified potential */
  const double Wminimum_scaling;
};

#endif // HARMONICOSCILLATORACTION_HH
