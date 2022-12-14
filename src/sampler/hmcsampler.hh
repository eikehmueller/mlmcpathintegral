#ifndef HMCSAMPLER_HH
#define HMCSAMPLER_HH HMCSAMPLER_HH
#include "action/action.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_random.hh"
#include "mpi/mpi_wrapper.hh"
#include "sampler/sampler.hh"
#include <memory>
#include <random>
#include <vector>

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
  HMCParameters()
      : Parameters("hmc"), nt_(100), dt_(0.1), n_rep_(1), n_burnin_(100) {
    addKey("nt", Integer, Positive);
    addKey("dt", Double, Positive);
    addKey("n_rep", Integer, Positive);
    addKey("n_burnin", Integer, Positive);
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
   * @param[in] hmc_param_ HMC parameters
   */
  HMCSampler(const std::shared_ptr<Action> action_,
             const HMCParameters hmc_param_)
      : Sampler(), action(action_), nt_hmc(hmc_param_.nt()),
        dt_hmc(hmc_param_.dt()), n_rep(hmc_param_.n_rep()),
        n_burnin(hmc_param_.n_burnin()) {
    engine.seed(8923759);
    // Create temporary workspace
    // current position
    phi_state_cur = std::make_shared<SampleState>(action->sample_size());
    // current (conjugate) momentum
    p_state_cur = std::make_shared<SampleState>(action->sample_size());
    // Trial state
    phi_state_trial = std::make_shared<SampleState>(action->sample_size());
    // Momentum change from force term
    dp_state = std::make_shared<SampleState>(action->sample_size());
    action->initialise_state(phi_state_cur);
    // Burn in
    std::shared_ptr<SampleState> phi_state_tmp =
        std::make_shared<SampleState>(action->sample_size());
    for (unsigned int i = 0; i < n_burnin; ++i) {
      draw(phi_state_tmp);
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
   * returns a sample state \f$\phi\f$
   *
   * @param[out] phi_state Path \f$\phi\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<SampleState> phi_state);

  /** @brief Set current state to particular value
   *
   * @param[in] phi_state
   */
  virtual void set_state(std::shared_ptr<SampleState> phi_state);

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
  /** @brief Current state */
  mutable std::shared_ptr<SampleState> phi_state_cur;
  /** @brief temporary for momenta */
  mutable std::shared_ptr<SampleState> p_state_cur;
  /** @brief Trial state */
  mutable std::shared_ptr<SampleState> phi_state_trial;
  /** @brief temporary increment for momenta */
  mutable std::shared_ptr<SampleState> dp_state;
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

class HMCSamplerFactory : public SamplerFactory {
public:
  /** @brief Create new instance
   *
   * @param[in] param_hmc HMC parameters
   */
  HMCSamplerFactory(const HMCParameters param_hmc_) : param_hmc(param_hmc_) {}

  /** @brief Destructor */
  virtual ~HMCSamplerFactory() {}

  /** @brief Return sampler for a specific  action
   *
   * @param[in] action Action to sample from
   */
  virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
    return std::make_shared<HMCSampler>(action, param_hmc);
  }

private:
  /** HMC parameters */
  const HMCParameters param_hmc;
};

#endif // HMCSAMPLER_HH
