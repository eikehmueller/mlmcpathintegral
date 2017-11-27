#ifndef HARMONICOSCILLATORACTION_HH
#define HARMONICOSCILLATORACTION_HH HARMONICOSCILLATORACTION_HH
#include <cassert>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "path.hh"
#include "sampler.hh"
#include "action.hh"

/** @class HarmonicOscillatorAction
 * @brief Action class for harmonic oscillator 
 *
 * Action class for potential \f$V(x)=\frac{m_0}{2}\mu^2x^2\f$
 */
class HarmonicOscillatorAction : public Action, public Sampler {
public:
  /* @brief Initialise class
   *
   * 
   * @param[in] M_lat Number of time slices \f$M\f$
   * @param[in] T_final Final time \f$T\f$
   * @param[in] m0: Mass of particle \f$m_0\f$
   * @param[in] mu2: Frequency \f$\mu^2\f$
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
   * Delete temporary arrays
   */
  ~HarmonicOscillatorAction() {
    delete y_tmp;
  }
  /* @brief Evaluate action for a specific path
   * 
   * Calculate \$S[X]\f$ for a specific path
   *
   * @param x_path: path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(const Path* x_path) const;

  /* @brief Draw sample from distribution
   *
   * Generate path which is distruted according to
   * \f$\pi(X) \sim e^{-S[X]}\f$
   *
   * @param x_path: path to populate
   */
  const virtual void draw(std::vector<Path*> x_path);

  /** @brief Construct Cholesky factor of covariance matrix */
  void build_covariance();
private:
  const double mu2;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> Vector;
  /** Cholesky factor of the covariance matrix */
  Matrix L_cov;
  /** Random number engine */
  typedef std::mt19937_64 Engine;
  mutable Engine engine;
  /** Normal distribution */
  typedef std::normal_distribution<double> Normal;
  mutable Normal normal_dist;
  Path* y_tmp; // temporary array used for direct sampling
};

#endif // HARMONICOSCILLATORACTION_HH
