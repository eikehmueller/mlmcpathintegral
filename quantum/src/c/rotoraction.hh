#ifndef ROTORACTION_HH
#define ROTORACTION_HH ROTORACTION_HH
#include "config.h"
#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include "auxilliary.hh"
#include "path.hh"
#include "clusteraction.hh"
#include "parameters.hh"

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
  RotorParameters() :
    Parameters("rotor"),
    m0_(1.0),
    renormalisation_(RenormalisationNone) {
    addKey("m0",Double,Positive);
    addKey("renormalisation",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      m0_ = getContents("m0")->getDouble();
      std::string renormalisation_str = getContents("renormalisation")->getString();
      std::cout << "renormalisation_str " << renormalisation_str << std::endl;
      if (renormalisation_str == "none") {
        renormalisation_ = RenormalisationNone;
      } else if (renormalisation_str == "perturbative") {
        renormalisation_ = RenormalisationPerturbative;
      } else if (renormalisation_str == "exact") {
        renormalisation_ = RenormalisationExact;
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
class RotorAction : public ClusterAction {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] m0_ Moment of inertia (angular mass) \f$I\f$
   */
  RotorAction(const unsigned int M_lat_,
              const double T_final_,
              const RenormalisationType renormalisation_,
              const double m0_)
    : ClusterAction(M_lat_,T_final_,renormalisation_,m0_),
      uniform_dist(-M_PI,M_PI) {
    engine.seed(21172817);
  }

  /** @brief Tidy up
   * 
   * Delete any temporary arrays
   */
  ~RotorAction() {}

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  std::shared_ptr<Action> virtual coarse_action() {
    if (M_lat%2) {
      std::cerr << "ERROR: cannot coarsen action, number of lattice sites is odd." << std::endl;
      exit(1);
    }
    RenormalisedRotorParameters c_param(M_lat,T_final,m0,renormalisation);
    return std::make_shared<RotorAction>(M_lat/2,
                                         T_final,
                                         renormalisation,
                                         c_param.m0_coarse());
  };
  
  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path path \f$X\f$, has to be am array of length \f$M\f$
   */
  const double virtual evaluate(const std::shared_ptr<Path> x_path) const;

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
  void virtual force(const std::shared_ptr<Path> x_path,
                     std::shared_ptr<Path> p_path) const;

  /** @brief Initialise path 
   *
   * Set initial values of path to random values in the range [-pi,pi)
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_path(std::shared_ptr<Path> x_path) const;
  
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
    return 2.0*m0/a_lat*fabs(cos(0.5*(x_p-x_m)));
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
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const {
    return atan2(sin(x_p)+sin(x_m),cos(x_p)+cos(x_m));
  }

  /** @brief Change \f$S_{\ell}\f$ in energy used in bonding probabilities
   *
   * For this action we have
   * \f[
   *    S_{\ell} = -2\frac{I}{a}\cos(\phi_i-\alpha_r)\cos(\phi_{i+1}-\alpha_r)
   * \f]
   *
   * @param[in] x_m Value of \f$x^{(\ell)}_- = x_i\f$
   * @param[in] x_p Value of \f$x^{(\ell)}_+ = x_{i+1}\f$
   */
  virtual double S_ell(const double x_m, const double x_p) const {
    return -2.0*m0/a_lat*cos(x_m-xbar)*cos(x_p-xbar);
  }

  /** @brief Pick angle for next step of cluster algorithm
   */
  virtual void new_angle() const {
    xbar = uniform_dist(engine);
  }

  /** @brief Flip a site
   * 
   * Return \f$hx\f$
   *
   * @param[in] x value of site \f$x\f$
   */
  virtual double flip(const double x) const {
    return mod_2pi(M_PI+2.*xbar-x);
  }

protected:
  /** @brief reflection angle for current subgroup \f$H_{\overline{x}}\f$ */
  mutable double xbar;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Uniform distribution used for selecting subgroup (elements) */
  mutable Uniform uniform_dist;
};

#endif // ROTORACTION_HH
