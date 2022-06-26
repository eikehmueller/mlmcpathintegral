#ifndef QUENCHEDSCHWINGERACTION_HH
#define QUENCHEDSCHWINGERACTION_HH QUENCEDSCHWINGERACTION_HH
#include "action/action.hh"
#include "action/qft/qftaction.hh"
#include "action/qft/quenchedschwingerrenormalisation.hh"
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "distribution/expcosdistribution.hh"
#include "lattice/lattice2d.hh"
#include "mpi/mpi_random.hh"
#include "mpi/mpi_wrapper.hh"
#include <cassert>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

/** @file quenchedschwingeraction.hh
 * @brief Header file for action of the quenched 2D Schwinger model
 */

/** @class SchwingerParameters
 *
 * @brief Class for storing parameters of Schwinger action
 *
 * This stores the (physical) coupling parameter \f$\beta\f$ of the 2D Schwinger
 * model
 */
class SchwingerParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  SchwingerParameters()
      : Parameters("schwinger"), beta_(1.0),
        renormalisation_(RenormalisationNone) {
    addKey("beta", Double, Positive);
    addKey("renormalisation", String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      beta_ = getContents("beta")->getDouble();
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

  /** @brief Return coupling parameter \f$\beta\f$ */
  double beta() const { return beta_; }
  /** @brief Return renormalisation */
  RenormalisationType renormalisation() const { return renormalisation_; }

private:
  /** @brief Coupling parameter \f$\beta\f$ */
  double beta_;
  /** @brief Renormalisation */
  RenormalisationType renormalisation_;
};

/** @class QuenchedSchwingerAction
 *
 * @brief Class for action of quenched 2D Schwinger model
 *
 * Allows calculation with the pure gauge action
 * \f[
 *      S[U]=\beta \sum_P \left[1-Re(U_P) \right]
 * \f]
 * where \f$P = U_{n,0} U_{n+\hat{0},1}
 * U^\dagger_{n+\hat{1},0}U^\dagger_{0,1}\f$ is the plaquette at the site (or
 * lattice point) \f$n=(i,j)\f$. \f$beta=\frac{1}{a_t a_x g^2}\f$ is the
 * dimensionless coupling constant.
 *
 * As the link at site \f$n\f$ in direction \f$\mu\f$ is defined as
 * \f[
 *       U_{n,\mu} = e^{i \theta_{n,\mu}}
 * \f]
 * the action can be expressed in terms of the real variables
 * \f$\theta_{n,\mu}\f$ as \f[ S[\theta]=\beta \sum_P
 * \left[1-\cos\left(\theta_{n,0}+\theta_{n+\hat{0},1}-\theta_{n+\hat{1},0}-\theta_{n,1}\right))
 * \right] \f]
 *
 *  Here if \f$n=(i,j)\f$ then \f$n+\hat{0}=(i+1,j)\f$, \f$n+\hat{1}=(i,j+1)\f$
 * etc.
 */

class QuenchedSchwingerAction : public QFTAction {
public:
  /** @brief Initialise class
   *
   * @param[in] lattice_ Underlying two-dimensional lattice
   * @param[in] lattice_ Underlying fine level two-dimensional lattice
   * @param[in] renormalisation_ Type of renormalisation
   * @param[in] beta_ non-dimensionalised coupling constant \f$\beta=1/(g^2 a_t
   * a_x)\f$
   */
  QuenchedSchwingerAction(const std::shared_ptr<Lattice2D> lattice_,
                          const std::shared_ptr<Lattice2D> fine_lattice_,
                          const RenormalisationType renormalisation_,
                          const double beta_)
      : QFTAction(lattice_, fine_lattice_, renormalisation_), beta(beta_),
        exp_cos_dist(beta) {
    engine.seed(2481317);
    CoarseningType coarsening_type = lattice->get_coarsening_type();
    if (not((coarsening_type == CoarsenBoth) or
            (coarsening_type == CoarsenTemporal) or
            (coarsening_type == CoarsenSpatial) or
            (coarsening_type == CoarsenAlternate))) {
      mpi_parallel::cerr
          << "ERROR: invalid coarsening for quenched Schwinger model."
          << std::endl;
      mpi_parallel::cerr
          << "Has to be 'both', 'temporal', 'spatial' or 'alternate'."
          << std::endl;
      mpi_exit(EXIT_FAILURE);
      throw std::runtime_error("...");
    }
  }

  /** @brief Return coupling constant \f$beta\f$ */
  double getbeta() const { return beta; }

  /** @brief return size of a sample, i.e. the number of links on the lattice */
  virtual unsigned int sample_size() const { return lattice->getNedges(); };

  /** @brief Cost of one action evaluation */
  virtual double evaluation_cost() const { return sample_size(); }

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the lattice hierarchy.
   */
  virtual std::shared_ptr<Action> coarse_action() {
    RenormalisedQuenchedSchwingerParameters c_param(lattice, beta,
                                                    renormalisation);
    // Construct coarse lattice
    std::shared_ptr<Lattice2D> coarse_lattice = lattice->get_coarse_lattice();
    if (coarse_lattice == nullptr) {
      mpi_parallel::cerr << "ERROR: cannot coarsen 2d lattice with M_{t,lat} = "
                         << lattice->getMt_lat();
      mpi_parallel::cerr << " , M_{x,lat} = " << lattice->getMx_lat() << "."
                         << std::endl;
      mpi_exit(EXIT_FAILURE);
      throw std::runtime_error("...");
    }
    // Construct coarse action
    std::shared_ptr<Action> new_action;
    new_action = std::make_shared<QuenchedSchwingerAction>(
        coarse_lattice, lattice, renormalisation, c_param.beta_coarse());
    return new_action;
  };

  /** @brief Evaluate action for a specific state
   *
   * Calculate \f$S[\phi]\f$ for a specific state
   *
   * @param[in] phi_state Sample state \f$\phi\f$
   */
  virtual const double
  evaluate(const std::shared_ptr<SampleState> phi_state) const;

  /** @brief Draw local value of state from heat bath
   *
   * Update the local entry at position \f$\ell\f$ of the state using a heat
   * bath defined by the neighbouring sites
   *
   *  @param[inout] phi_state State to update
   *  @param[in] ell index of dof to update
   */
  virtual void heatbath_update(std::shared_ptr<SampleState> phi_state,
                               const unsigned int ell);

  /** @brief Perform local overrelaxation update
   *
   * Update the local entry at position \f$\ell\f$ of the state using
   * overrelaxation
   *
   *  @param[inout] phi_state State to update
   *  @param[in] ell index of dof to update
   */
  virtual void overrelaxation_update(std::shared_ptr<SampleState> phi_state,
                                     const unsigned int ell);

  /** @brief Calculate force for HMC integrator for a specific state
   *
   * Calculate \f$P = \frac{\partial S[\phi]}{\partial \phi}\f$ for a specific
   * state and return the resulting force as a state.
   *
   * @param[in] phi_state State \f$\phi\f$ on which to evaluate the force
   * @param[out] p_state Resulting force \f$P\f$ at every point
   *
   */
  virtual void force(const std::shared_ptr<SampleState> phi_state,
                     std::shared_ptr<SampleState> p_state) const;

  /** @brief Copy coarse data points from sample on coarser level
   *
   * @param[in] phi_coarse Coarse sample to copy from
   * @param[in] phi_state Fine state to copy to (sample level as action)
   */
  virtual void copy_from_coarse(const std::shared_ptr<SampleState> phi_coarse,
                                std::shared_ptr<SampleState> phi_state);

  /** @brief Copy coarse data points from state on finer level
   *
   * @param[in] phi_fine Fine state to copy from
   * @param[in] phi_coarse Coarse state to copy to (same level as action)
   */
  virtual void copy_from_fine(const std::shared_ptr<SampleState> phi_fine,
                              std::shared_ptr<SampleState> phi_state);

  /** @brief Initialise state
   *
   * Set initial values of state, those values will be used to start the
   * sampling process
   *
   * @param[out] phi_state State \f$\phi\f$ to be set
   */
  virtual void initialise_state(std::shared_ptr<SampleState> phi_state) const;

  /** @brief Action information string
   *
   * return some information on this instance of the action
   */
  virtual std::string info_string() const;

private:
  /** @brief Compute staple angles \f$\theta_{n,\mu}^{(+)}\f$ and
   * \f$\theta_{n,\mu}^{(-)}\f$
   *
   * Let \f$\nu=1-mu\f$. Then the staple angles are defined as
   * \f[
   *    \theta_{n,\mu}^{(+)} = \theta_{n,\nu} + \theta_{n+\hat{\nu},\mu} -
   * \theta_{n+\hat{\mu},\nu} \f] and \f[ \theta_{n,\mu}^{(-)} =
   * \theta_{n-\hat{\nu}+\hat{\mu},\nu} + \theta_{n-\hat{\nu},\mu} -
   * \theta_{n-\hat{\nu},\nu} \f] Note that with those angles the action can be
   * written as a sum over all links \f[ S[U] = \frac{\beta}{4}\sum_{n}
   * \sum_{\mu=0,1} \left[2 - \cos(\theta_{n,\mu}-\theta_{n,\mu}^{(+)})
   *                                           -
   * \cos(\theta_{n,\mu}-\theta_{n,\mu}^{(-)})\right] \f]
   *
   *  @param[in] phi_state State vector, containing the angles
   * \f$\theta_{n,\mu}\f$
   *  @param[in] i Temporal index of coordinate \f$n\f$
   *  @param[in] j Temporal index of coordinate \f$n\f$
   *  @param[in] mu Direction \f$\mu\f$
   *  @param[out] theta_p Resulting value of \f$\theta_{n,\mu}^{(+)}\f$
   *  @param[out] theta_m Resulting value of \f$\theta_{n,\mu}^{(-)}\f$
   */
  void compute_staple_angles(std::shared_ptr<SampleState> phi_state,
                             const int i, const int j, const int mu,
                             double &theta_p, double &theta_m);

protected:
  /** @brief Dimensionless coupling constant \f$\beta=1/(a_t a_x g^2)\f$*/
  const double beta;
  /** @brief Random number engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief distribution for drawing from heat bath */
  const ExpCosDistribution exp_cos_dist;
};

#endif // QUENCHEDSCHWINGERACTION_HH
