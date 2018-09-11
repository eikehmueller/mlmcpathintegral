#ifndef HMCSAMPLER_HH
#define HMCSAMPLER_HH HMCSAMPLER_HH
#include <memory>
#include <random>
#include <vector>
#include "path.hh"
#include "action.hh"
#include "parameters.hh"

/** @file hmcsampler.hh
 * @brief Header file for Hybrid Monte Carlo sampler base class
 */

/** @class HMCParameters
 *
 * @brief Class for storing parameters of Hybrid MC sampler
 */
class HMCParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  HMCParameters() :
    Parameters("hmc"),
    T_(1.0),
    dt_(0.1),
    n_burnin_(100) {
    addKey("T",Double,Positive);
    addKey("dt",Double,Positive);
    addKey("n_burnin",Integer,Positive);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      T_ = getContents("T")->getDouble();
      dt_ = getContents("dt")->getDouble();
      n_burnin_ = getContents("n_burnin")->getInt();
    }
    return readSuccess;
  }
  /** @brief Return length of integration interval \f$T\f$ */
  double T() const { return T_; }
  /** @brief Return integration timestep \f$dt\f$ */
  double dt() const { return dt_; }
  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
private:
  /** @brief Integration time interval \f$T\f$ */
  double T_;
  /** @brief Integration time step \f$dt\f$ */
  double dt_;
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
};

/** @class HMCSampler
 * @brief Hybrid Monte Carlo sampler
 *
 * Implementation of HMC sampling. Generate states by following
 * deterministic trajectories of the system
 * \f$H(P,X) = \frac{1}{2}P^2 + S(X)\f$ with a simple sympectic integrator
 * and the accept/reject the resulting state. The length of the deterministic
 * trajectories is \f$T_{HMC}\f$ and a timestep size of \f$\Delta t_{HMC}\f$ 
 * is used with a simple Symplectic-Euler method.
 */
class HMCSampler : public Sampler {
public:
  /** @brief Base class */
  typedef Sampler Base;
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] T_hmc_ Length \f$T_{HMC}\f$ of HMC trajectories 
   * @param[in] dt_hmc_ Step size \f$\Delta t_{HMC}\f$ of HMC trajectories 
   * @param[in] n_burnin_ Number of burnin steps
   */
  HMCSampler(const std::shared_ptr<Action> action_,
             const double T_hmc_,
             const double dt_hmc_,
             const unsigned int n_burnin_) :
    Base(true),
    action(action_),
    T_hmc(T_hmc_),
    dt_hmc(dt_hmc_),
    n_burnin(n_burnin_)
  {
    engine.seed(8923759233);
    const unsigned int M_lat = action->getM_lat();
    const double T_final = action->getT_final();
    // Create temporary workspace
    // current position
    x_path_cur = std::make_shared<Path>(M_lat,T_final);
    // current (conjugate) momentum
    p_path_cur = std::make_shared<Path>(M_lat,T_final);
    // Trial path
    x_path_trial = std::make_shared<Path>(M_lat,T_final);
    // Momentum change from force term
    dp_path = std::make_shared<Path>(M_lat,T_final);
    action->initialise_path(x_path_cur);
    // Burn in
    std::vector<std::shared_ptr<Path>> x_path_tmp;
    x_path_tmp.push_back(std::make_shared<Path>(M_lat,T_final));
    for (unsigned int i=0;i<n_burnin;++i) {
      draw(x_path_tmp);
    }
  }

  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~HMCSampler() {}

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::vector<std::shared_ptr<Path>> x_path);

protected:
  /** @brief Action to sample from */
  const std::shared_ptr<Action> action;
  /** @brief Length of deterministic trajectories */
  const double T_hmc;
  /** @brief Time step size of determininistic trajectories */
  const double dt_hmc;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
  /** @brief Current state (path) */
  mutable std::shared_ptr<Path> x_path_cur;
  /** @brief temporary for momenta */
  mutable std::shared_ptr<Path> p_path_cur;
  /** @brief Trial state (path) */
  mutable std::shared_ptr<Path> x_path_trial;
  /** @brief temporary increment for momenta */
  mutable std::shared_ptr<Path> dp_path;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of normal distribution */
  typedef std::normal_distribution<double> Normal;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Normal distribution used for momentum sampling */
  Normal normal_dist;
  /** @brief Uniform distribution used for Metropolis accept/reject step */
  Uniform uniform_dist;
};

#endif // HMCSAMPLER_HH
