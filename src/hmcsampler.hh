#ifndef HMCSAMPLER_HH
#define HMCSAMPLER_HH HMCSAMPLER_HH
#include <memory>
#include <random>
#include <vector>
#include "path.hh"
#include "action.hh"
#include "parameters.hh"
#include "sampler.hh"
#include "mpi_wrapper.hh"
#include "mpi_random.hh"

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
    nt_(100),
    dt_(0.1),
    n_rep_(1),
    n_burnin_(100) {
    addKey("nt",Integer,Positive);
    addKey("dt",Double,Positive);
    addKey("n_rep",Integer,Positive);
    addKey("n_burnin",Integer,Positive);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      nt_ = getContents("nt")->getInt();
      dt_ = getContents("dt")->getDouble();
      n_rep_ = getContents("n_rep")->getInt();
      n_burnin_ = getContents("n_burnin")->getInt();
    }
    return readSuccess;
  }
  /** @brief Return number of integration steps */
  unsigned int nt() const { return nt_; }
  /** @brief Return number of repetitions */
  unsigned int n_rep() const { return n_rep_; }
  /** @brief Return integration timestep \f$dt\f$ */
  double dt() const { return dt_; }
  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
private:
  /** @brief Number of integration steps */
  unsigned int nt_;
  /** @brief Integration time step \f$dt\f$ */
  double dt_;
  /** @brief Number of repetitions */
  unsigned int n_rep_;
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
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] nt_hmc_ Number of HMC steps 
   * @param[in] dt_hmc_ Step size \f$\Delta t_{HMC}\f$ of HMC trajectories 
   * @param[in] n_burnin_ Number of burnin steps
   */
  HMCSampler(const std::shared_ptr<Action> action_,
             const HMCParameters hmc_param_) :
    Sampler(),
    action(action_),
    nt_hmc(hmc_param_.nt()),
    dt_hmc(hmc_param_.dt()),
    n_rep(hmc_param_.n_rep()),
    n_burnin(hmc_param_.n_burnin())
  {
    engine.seed(8923759);
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
    std::shared_ptr<Path> x_path_tmp =
      std::make_shared<Path>(M_lat,T_final);
    for (unsigned int i=0;i<n_burnin;++i) {
      draw(x_path_tmp);
    }
    autotune_stepsize(0.8);
    // Reset counters
    reset_stats();
  }

  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~HMCSampler() {}

  /** @brief auto-tune step size
   *
   * Adjust HMC step to achieve a desired acceptance rate
   *
   * @param[in] p_accept_target target acceptance rate
   */
  void autotune_stepsize(const double p_accept_target);
  
  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<Path> x_path);

  /** @brief Set current state to particular value
   *
   * @param[in] x_path
   */
  virtual void set_state(std::shared_ptr<Path> x_path);

private:
  
  /** @brief Integrate a single HMC trajectory */
  bool single_step();
  
protected:
  /** @brief Action to sample from */
  const std::shared_ptr<Action> action;
  /** @brief Number of steps in trajectory */
  const unsigned int nt_hmc;
  /** @brief Time step size of determininistic trajectories */
  mutable double dt_hmc;
  /** @brief Number of repetitions */
  const unsigned int n_rep;
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
  typedef mpi_parallel::mt19937_64 Engine;
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
