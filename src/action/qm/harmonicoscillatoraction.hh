#ifndef HARMONICOSCILLATORACTION_HH
#define HARMONICOSCILLATORACTION_HH HARMONICOSCILLATORACTION_HH
#include "action/qm/harmonicoscillatorrenormalisation.hh"
#include "action/qm/qmaction.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
#include "mpi/mpi_random.hh"
#include "sampler/sampler.hh"
#include <Eigen/Dense>
#include <cassert>
#include <memory>
#include <random>
#include <vector>

/** @file harmonicoscillatoraction.hh
 * @brief Header file for harmonic oscillator action base class
 */

/** @class HarmonicOscillatorParameters
 *
 * @brief Class for storing parameters of harmonic oscillator action
 *
 * This stores the mass \f$m_0\f$ and curvature \f$\mu_2\f$ of the
 * harmonic oscillator action with potential \f$V(x)=\frac{m_0}{2}\mu^2x^2\f$
 */
class HarmonicOscillatorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  HarmonicOscillatorParameters()
      : Parameters("harmonicoscillator"), m0_(1.0), mu2_(1.0),
        renormalisation_(RenormalisationNone) {
    addKey("m0", Double, Positive);
    addKey("mu2", Double);
    addKey("renormalisation", String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
      mu2_ = getContents("mu2")->getDouble();
      std::string renormalisation_str =
          getContents("renormalisation")->getString();
      if (renormalisation_str == "none") {
        renormalisation_ = RenormalisationNone;
      } else if (renormalisation_str == "perturbative") {
        renormalisation_ = RenormalisationPerturbative;
      } else if (renormalisation_str == "nonperturbative") {
        renormalisation_ = RenormalisationNonperturbative;
      }
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return renormalisation */
  RenormalisationType renormalisation() const { return renormalisation_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Renormalisation */
  RenormalisationType renormalisation_;
};

/** @class HarmonicOscillatorAction
 *
 * @brief Action class for harmonic oscillator
 *
 * Action class for potential \f$V(x)=\frac{m_0}{2}\mu^2x^2\f$
 */
class HarmonicOscillatorAction : public QMAction, public Sampler {
public:
  /** @brief Initialise class
   *
   *
   * @param[in] lattice_ Underlying lattice
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] m0_ Mass of particle \f$m_0\f$
   * @param[in] mu2_ Frequency \f$\mu^2\f$
   */
  HarmonicOscillatorAction(const std::shared_ptr<Lattice1D> lattice_,
                           const RenormalisationType renormalisation_,
                           const double m0_, const double mu2_)
      : QMAction(lattice_, renormalisation_, m0_), Sampler(), mu2(mu2_),
        Wcurvature((2. / a_lat + a_lat * mu2) * m0),
        Wminimum_scaling(0.5 / (1. + 0.5 * a_lat * a_lat * mu2)) {
    build_covariance();
    engine.seed(124129017);
    y_tmp = std::make_shared<SampleState>(lattice_->getM_lat());
  }

  /** @brief Tidy up
   *
   * Delete any temporary arrays
   */
  ~HarmonicOscillatorAction() {}

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  std::shared_ptr<Action> virtual coarse_action() {
    RenormalisedHOParameters c_param(lattice, m0, mu2, renormalisation);
    std::shared_ptr<Action> new_action;
    new_action = std::make_shared<HarmonicOscillatorAction>(
        lattice->coarse_lattice(), renormalisation, c_param.m0_coarse(),
        c_param.mu2_coarse());
    return new_action;
  };

  /** @brief Evaluate action for a specific path
   *
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(
      const std::shared_ptr<SampleState> x_path) const;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path \f$P\f$.
   *
   * Note that for this action we have
     \f[
         P_j = \frac{m_0}{a}\left(2X_j-X_{j-1}-X_{j+1}\right) +
   \frac{1}{2}am_0\mu^2 X_j \f]
   *
   * @param x_path SampleState \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const std::shared_ptr<SampleState> x_path,
                     std::shared_ptr<SampleState> p_path) const;

  /** @brief Initialise path
   *
   * Set initial values of path to zero
   *
   * @param[out] x_path SampleState \f$X\f$ to be set
   */
  void virtual initialise_state(std::shared_ptr<SampleState> x_path) const {
    std::fill(x_path->data.data(), x_path->data.data() + x_path->data.size(),
              0.0);
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
  double virtual inline getWminimum(const double x_m, const double x_p) const {
    return Wminimum_scaling * (x_m + x_p);
  }

  /** @brief Draw sample from distribution
   *
   * Generate path which is distruted according to
   * \f$\pi(X) \sim e^{-S[X]}\f$
   *
   * @param[out] x_path: path to populate
   */
  virtual void draw(std::shared_ptr<SampleState> x_path);

  /** @brief Perturbative expression for expectation value of \f$X^2\f$
   *
   * Return the analytical expression from Creutz and Freedman:
   *
   \f[
   \langle X^2 \rangle = \frac{1}{2m_0\mu\sqrt{1+\frac{a^2\mu^2}{4}}}\cdot
   \frac{1+R^M}{1-R^M} \f]
   * where
   \f[
   R = 1 + \frac{a^2\mu^2}{2}-a\mu\sqrt{1+\frac{a^2\mu^2}{4}}
   \f]
   */
  const double Xsquared_analytical();

  /** @brief Continuum limitof expectation value of \f$X^2\f$
   *
   * Return the continuum limit of \f$\langle X^2\rangle\f$ from
   * Creutz and Freedman:
   *
   \f[
   \lim_{a\rightarrow 0}\langle X^2 \rangle = \frac{1}{2m_0\mu}\cdot
   \frac{1-e^{-\mu T}}{1+e^{-\mu T}} \f]
   * where
   \f[
   R = 1 + \frac{a^2\mu^2}{2}-a\mu\sqrt{1+\frac{a^2\mu^2}{4}}
   \f]
   */
  const double Xsquared_analytical_continuum();

  /** @brief Set current state to particular value
   *
   *  @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<SampleState> x_path){};

private:
  /** @brief Construct Cholesky factor of covariance matrix */
  void build_covariance();

private:
  /** @brief Oscillator frequency */
  const double mu2;
  /** @brief Eigen matrix type */
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  /** @brief Eigen vector type */
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  /** @brief Cholesky factor of the covariance matrix \f$\Sigma=L^TL\f$ */
  Matrix L_cov;
  /** @brief Type for Mersenne twister engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Random number engine */
  mutable Engine engine;
  /** @brief Type for normal distribution */
  typedef std::normal_distribution<double> Normal;
  /** @brief Normal distribution */
  mutable Normal normal_dist;
  /** @brief temporary array used for direct sampling */
  std::shared_ptr<SampleState> y_tmp;
  /** @brief Curvature of modified potential */
  const double Wcurvature;
  /** @brief Scaling factor for calculation of mimimum of modified potential */
  const double Wminimum_scaling;
};

class HarmonicOscillatorSamplerFactory : public SamplerFactory {
public:
  /** @brief Destructor */
  virtual ~HarmonicOscillatorSamplerFactory() {}

  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::dynamic_pointer_cast<Sampler>(action);
  }
};

#endif // HARMONICOSCILLATORACTION_HH
