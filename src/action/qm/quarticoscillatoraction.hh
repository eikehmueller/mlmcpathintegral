#ifndef QUARTICOSCILLATORACTION_HH
#define QUARTICOSCILLATORACTION_HH QUARTICOSCILLATORACTION_HH
#include <memory>
#include <vector>
#include "lattice/lattice1d.hh"
#include "fields/path.hh"
#include "action/action.hh"
#include "common/parameters.hh"
#include "mpi/mpi_wrapper.hh"

/** @file quarticoscillatoraction.hh
 * @brief Header file for quartic oscillator action class
 */

/** @class QuarticOscillatorParameters
 *
 * @brief Class for storing parameters of quartic oscillator action
 *
 * This stores the mass \f$m_0\f$ and parameters \f$\mu_2\f$, \f$\lambda\f$,
 * \f$x_0\f$ of the double well action with potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\frac{\lambda}{4}(x-x_0)^4\f$.
 */
class QuarticOscillatorParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  QuarticOscillatorParameters() :
    Parameters("quarticoscillator"),
    m0_(1.0),
    mu2_(1.0),
    lambda_(1.0),
    x0_(0.0) {
    addKey("m0",Double,Positive);
    addKey("mu2",Double);
    addKey("lambda",Double,NonNegative);
    addKey("x0",Double);
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
      lambda_ = getContents("lambda")->getDouble();
      x0_ = getContents("x0")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return parameter \f$\lambda\f$ */
  double lambda() const { return lambda_; }
  /** @brief Return parameter \f$x_0\f$ */
  double x0() const { return x0_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Parameter \f$\lambda\f$ */
  double lambda_;
  /** @brief Parameter \f$x_0\f$ */
  double x0_;
};

/** @class QuarticOscillatorAction
 *
 * @brief Action class for quartic oscillator 
 *
 * Action class for potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\frac{\lambda}{4}(x-x_0)^4\f$
 */
class QuarticOscillatorAction : public Action {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] lattice_ Underlying lattice
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] m0_ Mass of particle \f$m_0\f$
   * @param[in] mu2_ Frequency \f$\mu^2\f$
   * @param[in] lambda_ Coefficient of quartic term, \f$\lambda\f$
   * @param[in] x0_ Shift quartic term, \f$x_0\f$
   */
  QuarticOscillatorAction(const std::shared_ptr<Lattice1D> lattice_,
                          const RenormalisationType renormalisation_,
                          const double m0_,
                          const double mu2_,
                          const double lambda_,
                          const double x0_)
    : Action(lattice_,renormalisation_,m0_),
      M_lat(lattice_->getM_lat()),
      a_lat(lattice_->geta_lat()),
      mu2(mu2_), lambda(lambda_), x0(x0_) {
  }

  /** @brief Destructor */
  virtual ~QuarticOscillatorAction() {};
  
  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  std::shared_ptr<Action> virtual coarse_action() {
    std::shared_ptr<Action> new_action;
    new_action = std::make_shared<QuarticOscillatorAction>(lattice->coarse_lattice(),
                                                           renormalisation,
                                                           m0,
                                                           mu2,
                                                           lambda,
                                                           x0);
    return new_action;
  }
  
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
         P_j = \frac{m_0}{a}\left(2X_j-X_{j-1}-X_{j+1}\right) + am_0\mu^2 X_j + a\lambda (X_j-x_0)^3
     \f]
   *
   * @param x_path Path \f$X\f$ on which to evaluate the force
   * @param p_path Resulting force \f$P\f$ at every point
   */
  void virtual force(const std::shared_ptr<Path> x_path,
                     std::shared_ptr<Path> p_path) const;

    /** @brief Initialise path 
   *
   * Set initial values of path to zero.
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  void virtual initialise_path(std::shared_ptr<Path> x_path) const {
    x_path->fill([](){return 0.0;});
  }

  /** @brief Second derivative \f$W''_{\overline{x}}(x)\f$ of conditioned action
   *
   * For the harmonic oscillator potential the curvature of the modified
   * action (see Action::getWcurvature()) is 
   \f[
   W''_{x_-,x_+} = \frac{2m_0}{a}+am_0\mu^2 + 3a\lambda x^2
   \f]
   * where \f$x=\frac{x_-+x_+}{2}\f$.
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const {
    double x = 0.5*(x_m+x_p);
    return (2./a_lat + a_lat*mu2)*m0 + 3.*lambda*a_lat*(x-x0)*(x-x0);
  }

  /** @brief Find minimum of conditioned action \f$W_{\overline{x}}(x)\f$
   *
   * For the quartic oscillator potential the minimum \f$x_0\f$ of the modified
   * action (see Action::getWminimum()) can be found as the solution of
   \f[
      x_0\left(1+\frac{1}{2}a^2\mu^2\right) + \frac{\lambda a^2}{2m_0}  x_0^3 = \overline{x}
   \f]
   * where \f$x=\frac{x_-+x_+}{2}\f$. Here, we calculate an approximate value
   * by setting \f$x_0^{(0)} = \overline{x}\f$ and then iterating
   \f[
     x_0^{(n+1)} = \left(1+\frac{1}{2}a^2\mu^2\right)^{-1}\left(\overline{x} - \frac{\lambda a^2}{2m_0} \left(x_0^{(n)}\right)^3\right)
   \f]
   * for \f$n=0,\dots,2\f$
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const {
    const double xbar = 0.5*(x_m+x_p);
    const double rho = 1./(1.+0.5*a_lat*a_lat*mu2);
    double x = xbar;
    for (int i=0;i<4;++i) {
      double x_shifted = x-x0;
      x = rho*(xbar - 0.5*a_lat*a_lat*lambda/m0*x_shifted*x_shifted*x_shifted);
    }
    return x;
  }

private:
  /** @brief Number of lattice points */
  const unsigned int M_lat;
  /** @brief Lattice spacing */
  const double a_lat;
  /** @brief Oscillator frequency */
  const double mu2;
  /** @brief Coefficient of quartic term in potential */
  const double lambda;
  /** @brief Shift of quartic term in potential */
  const double x0;
};

#endif // QUARTICOSCILLATORACTION_HH
