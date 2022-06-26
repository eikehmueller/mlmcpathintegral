#ifndef ROTORACTION_HH
#define ROTORACTION_HH ROTORACTION_HH
#include "action/clusteraction.hh"
#include "action/qm/qmaction.hh"
#include "action/qm/rotorrenormalisation.hh"
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "config.h"
#include "distribution/expsin2distribution.hh"
#include "lattice/lattice1d.hh"
#include "mpi/mpi_random.hh"
#include "mpi/mpi_wrapper.hh"
#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

/** @file rotoraction.hh
 * @brief Header file for quantum mechanical rotor action class
 */

/** @class RotorParameters
 *
 * @brief Class for storing parameters of quantum rotor action
 *
 * This stores the mass \f$m_0\f$ of the quantum mechanical rotor.
 */
class RotorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  RotorParameters()
      : Parameters("rotor"), m0_(1.0), renormalisation_(RenormalisationNone) {
    addKey("m0", Double, Positive);
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
  /** @brief Return renormalisation */
  RenormalisationType renormalisation() const { return renormalisation_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Renormalisation */
  RenormalisationType renormalisation_;
};

/** @class RotorAction
 *
 * @brief Action class for the quantum mechanical rotor descibed in
 *        <a href="https://arxiv.org/abs/1503.05088">arXiv/1503.05088</a>
 *
 * Action class for free particle moving on a circle of radius \f$R\f$ with
 * moment of inertia (angular mass) \f$I=m_0R^2\f$. The discrete action is
 * given by
 * \f[
 *   S = \frac{I}{a}\sum_{i=0}^{M_{lat}-1} \left(1-\cos(x_i-x_{i-1})\right)
 * \f]
 * with periodic boundary conditions \f$x_{M_{lat}}=x_0\f$.
 *
 * This action can be use in the cluster algorithm. The group \f$G\f$ of
 * transformations is given by shifts \f$x\mapsto (x + z) \mod [-\pi,\pi)\f$
 * with \f$\overline{x}\in[-\pi,\\pi)\f$ and reflections \f$x\mapsto -x\f$.
 * The group is generated by the subgroups
 * \f$H_{\overline{x}}=\{h_1,h_2\}\f$ with \f$\overline{x}\in[-\pi,\pi)\f$. The
 * two subgroup elements are \f$h_1(x)=x\f$ and
 * \f$h_2(x)=\pi+2\overline{x}-x\f$.
 */
class RotorAction : public QMAction, public ClusterAction {
public:
  /** @brief Initialise class
   *
   *
   * @param[in] lattice_ Underlying lattice
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] m0_ Moment of inertia (angular mass) \f$I\f$
   */
  RotorAction(const std::shared_ptr<Lattice1D> lattice_,
              const RenormalisationType renormalisation_, const double m0_)
      : QMAction(lattice_, renormalisation_, m0_), ClusterAction(lattice_),
        uniform_dist(-M_PI, M_PI) {
    engine.seed(21172817);
  }

  /** @brief Tidy up
   *
   * Delete any temporary arrays
   */
  virtual ~RotorAction() {}

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  std::shared_ptr<Action> virtual coarse_action() {
    RenormalisedRotorParameters c_param(lattice, m0, renormalisation);
    std::shared_ptr<Action> new_action = std::make_shared<RotorAction>(
        lattice->coarse_lattice(), renormalisation, c_param.m0_coarse());
    return new_action;
  };

  /** @brief return size of samples */
  virtual unsigned int sample_size() const { return M_lat; }

  /** @brief Evaluate action for a specific path
   *
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(
      const std::shared_ptr<SampleState> x_path) const;

  /** @brief Draw local value of state from heat bath
   *
   * Update the local entry at position j of the path using a heat bath defined
   * by the neighbouring sites
   *
   *  @param[inout] phi_state Path to update
   *  @param[in] ell index of dof to update
   */
  virtual void heatbath_update(std::shared_ptr<SampleState> x_path,
                               const unsigned int ell);

  /** @brief Perform local overrelaxation update
   *
   * Update the local entry at position j of the state using overrelaxation
   *
   *  @param[inout] phi_state State to update
   *  @param[in] ell index of dof to update
   */
  virtual void overrelaxation_update(std::shared_ptr<SampleState> x_path,
                                     const unsigned int ell);

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path \f$P\f$.
   *
   * Note that for this action we have
     \f[
         P_j = \frac{I}{a}\left(\sin(x_j-x_{j-1})+\sin(x_j-x_{j+1})\right)
     \f]
   *
   * @param x_path Path \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const std::shared_ptr<SampleState> x_path,
                     std::shared_ptr<SampleState> p_path) const;

  /** @brief Initialise path
   *
   * Set initial values of path to random values in the range [-pi,pi)
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_state(std::shared_ptr<SampleState> x_path) const;

  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned action
   *
   * For the quantum mechanical rotor the curvature of the modified
   * action (see Action::getWcurvature()) is
   \f[
     W''_{x_-,x_+} = \frac{I}{a}\left|\cos\left(\frac{x_+-x_-}{2}\right)\right|
   \f]
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const {
    return 2.0 * m0 / a_lat * fabs(cos(0.5 * (x_p - x_m)));
  }

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * For the quantum mechanical rotor the minimum of the modified
   * action (see Action::getWminimum()) can be found at
   \f[
      \tan(x_0) = \frac{\sin(x_-)+\sin(x_+)}{\cos(x_-)+\cos(x_+)}
   \f]
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m, const double x_p) const {
    return atan2(sin(x_p) + sin(x_m), cos(x_p) + cos(x_m));
  }

  /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
   *
   * For this action we have
   * \f[
   *    S_{\ell} = -2\frac{I}{a}\cos(\phi_i-\alpha_r)\cos(\phi_{i+1}-\alpha_r)
   * \f]
   *
   * @param[in] x_path Sample state
   * @param[in] i index of first vertex
   * @param[in] j index of second vertex (j=i+1)
   */
  virtual double S_ell(const std::shared_ptr<SampleState> x_path,
                       const unsigned int i, const unsigned int j) const {
    double x_m = x_path->data[i];
    double x_p = x_path->data[j];
    return -2.0 * m0 / a_lat * cos(x_m - xbar) * cos(x_p - xbar);
  }

  /** @brief Pick angle for next step of cluster algorithm
   */
  virtual void new_reflection() const { xbar = uniform_dist(engine); }

  /** @brief Flip the spin at a given site
   *
   * @param[inout] x_path State to process
   * @param[in] ell Vertex at which to flip the spin
   */

  /** @brief Flip a site
   *
   * Return \f$hx\f$
   *
   * @param[in] x value of site \f$x\f$
   */
  virtual void flip(std::shared_ptr<SampleState> x_path,
                    const unsigned int ell) const {
    double x = x_path->data[ell];
    x_path->data[ell] = mod_2pi(M_PI + 2. * xbar - x);
  }

  /** @brief Exact analytical expression for topological susceptibility
   *
   * Return analytical value for finite values of \f$a\f$
   */
  double chit_exact() const;

  /** @brief Perturbative expression for topological susceptibility
   *
   * Return perturbative value (up to corrections of \f$O((a/I)^2\f$) for finite
   * values of \f$a\f$
   */
  double chit_perturbative() const;

  /** @brief Analytical expression for topological susceptibility in the
   * continuum limit
   *
   * Return analytical expression of \f$\chi_t\f$ in the continuum limit
   * \f$a\rightarrow0\f$.
   */
  double chit_continuum() const;

protected:
  /** @brief reflection angle for current subgroup \f$H_{\overline{x}}\f$ */
  mutable double xbar;
  /** @brief Random number engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Uniform distribution used for selecting subgroup (elements) */
  mutable Uniform uniform_dist;
  /** @brief distribution for drawing from heat bath */
  const ExpSin2Distribution exp_sin2_dist;
};

#endif // ROTORACTION_HH
