#ifndef DOUBLEWELLACTION_HH
#define DOUBLEWELLACTION_HH DOUBLEWELLACTION_HH
#include <memory>
#include <vector>
#include "path.hh"
#include "action.hh"
#include "parameters.hh"
#include "mpi_wrapper.hh"

/** @file doublewellaction.hh
 * @brief Header file for double well action class
 */

/** @class DoubleWellParameters
 *
 * @brief Class for storing parameters of double well action
 *
 * This stores the mass \f$m_0\f$ and parameters \f$\mu_2\f$, \f$\lambda\f$
 * and \f$\sigma\f$ of the double well action with potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\lambda\exp\left(-\frac{x^2}{2\sigma^2}\right)\f$
 */
class DoubleWellParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  DoubleWellParameters() :
    Parameters("doublewell"),
    m0_(1.0),
    mu2_(1.0),
    lambda_(1.0),
    sigma_(1.0) {
    addKey("m0",Double,Positive);
    addKey("mu2",Double);
    addKey("lambda",Double);
    addKey("sigma",Double);
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
      sigma_ = getContents("sigma")->getDouble();
    }
    return readSuccess;
  }

  /** @brief Return unrenormalised mass \f$m_0\f$ */
  double m0() const { return m0_; }
  /** @brief Return parameter \f$\mu^2\f$ */
  double mu2() const { return mu2_; }
  /** @brief Return parameter \f$\lambda\f$ */
  double lambda() const { return lambda_; }
  /** @brief Return parameter \f$\sigma\f$ */
  double sigma() const { return sigma_; }

private:
  /** @brief Unrenormalised mass \f$m_0\f$ */
  double m0_;
  /** @brief Parameter \f$\mu^2\f$ */
  double mu2_;
  /** @brief Parameter \f$\lambda\f$ */
  double lambda_;
  /** @brief Parameter \f$\sigma\f$ */
  double sigma_;
};

/** @class DoubleWellAction
 *
 * @brief Action class for double well potential
 *
 * Action class for potential
 * \f$V(x)=\frac{m_0}{2}\mu^2x^2+\lambda\exp\left(-\frac{x^2}{2\sigma^2}\right)\f$
 */
class DoubleWellAction : public Action {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] M_lat_ Number of time slices \f$M\f$
   * @param[in] T_final_ Final time \f$T\f$
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] m0_ Mass of particle \f$m_0\f$
   * @param[in] mu2_ Frequency \f$\mu^2\f$
   * @param[in] lambda_ Coefficient of exponential term, \f$\lambda\f$
   * @param[in] sigma_ Width of exponential term, \f$\sigma\f$
   */
  DoubleWellAction(const unsigned int M_lat_,
                   const double T_final_,
                   const RenormalisationType renormalisation_,
                   const double m0_,
                   const double mu2_,
                   const double lambda_,
                   const double sigma_)
    : Action(M_lat_,T_final_,renormalisation_,m0_), mu2(mu2_), lambda(lambda_),
      sigma(sigma_),
      inv_sigma2(1./(sigma_*sigma_)) {
  }

  /** @brief destructor */
  virtual ~DoubleWellAction() {};
  
  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  std::shared_ptr<Action> virtual coarse_action() {
    if (M_lat%2) {
      mpi_parallel::cerr << "ERROR: cannot coarsen action, number of lattice sites is odd." << std::endl;
      exit(1);
    }
    return std::make_shared<DoubleWellAction>(M_lat/2,
                                              T_final,
                                              renormalisation,
                                              m0,
                                              mu2,
                                              lambda,
                                              sigma);
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
         P_j = \frac{m_0}{a}\left(2X_j-X_{j-1}-X_{j+1}\right) + am_0\mu^2 X_j - \frac{a\lambda X_j}{\sigma^2}\exp\left(-\frac{X_j^2}{2\sigma^2}\right)
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
    std::fill(x_path->data,x_path->data+M_lat,0.0);
  }

  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned action
   * at the minimum
   *
   * For the harmonic oscillator potential the curvature of the modified
   * action (see Action::getWcurvature()) at the minimum is 
   \f[
   W''_{x_-,x_+} = \frac{2m_0}{a}+am_0\mu^2 + \frac{a\lambda}{\sigma^2}\left(\frac{\overline{x}^2}{\sigma^2}-1\right)\exp\left(-\frac{\overline{x}^2}{2\sigma^2}\right)
   \f]
   * where \f$\overline{x}=\frac{x_-+x_+}{2}\f$.
   * 
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const {
    double x = getWminimum(x_m,x_p);
    return (2./a_lat + a_lat*mu2)*m0 + lambda*a_lat*inv_sigma2*(inv_sigma2*x*x-1.)*exp(-0.5*inv_sigma2*x*x);
  }

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * For the quartic oscillator potential the minimum \f$x_0\f$ of the modified
   * action (see Action::getWminimum()) can be found as the solution of
   \f[
      x_0\left(1+\frac{1}{2}a^2\mu^2\right) -\frac{\lambda a^2}{2m_0\sigma^2}\exp\left(-\frac{x^2}{2\sigma^2}\right) = \overline{x}
   \f]
   * Here, we calculate an approximate value by setting
   * \f$x_0^{(0)} = \overline{x}\f$ and iterating
   \f[
     x_0^{(n+1)} = \left(1+\frac{1}{2}a^2\mu^2\right)^{-1}\left(\overline{x} + \frac{\lambda a^2}{2m_0 \sigma^2} x_0^{(n)}\exp\left(-\frac{\left(x_0^{(n)}\right)^2}{2\sigma^2}\right)\right)
   \f]
   * for \f$n=0,\dots,2\f$ where \f$\overline{x}=\frac{x_++x_-}{2}\f$.
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const {
    double xbar = 0.5*(x_m+x_p);
    const double rho = 1./(1.+0.5*a_lat*a_lat*mu2);
    double x = xbar;
    for (int i=0;i<4;++i) {
      x = rho*(xbar + 0.5*a_lat*a_lat*lambda*inv_sigma2/m0*x*exp(-0.5*inv_sigma2*x*x));
    }
    return x;
  }

private:
  /** @brief Oscillator frequency */
  const double mu2;
  /** @brief Coefficient of quartic term in potential */
  const double lambda;
  /** @brief Parameter sigma in potential */
  const double sigma;
  /** @brief Inverse squared width of exponential term in potential */
  const double inv_sigma2;
};

#endif // DOUBLEWELLACTION_HH
