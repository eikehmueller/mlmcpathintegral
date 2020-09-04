#ifndef QMACTION_HH
#define QMACTION_HH QMACTION_HH

#include <string>
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "lattice/lattice1d.hh"
#include "action/action.hh"
#include "action/renormalisation.hh"
#include "mpi/mpi_wrapper.hh"

/** @file qmaction.hh
 * @brief Header file for quantum mechanical action base class
 */

/** @brief Enum for different actions
 *  - 0: Harmonic Oscillator
 *  - 1: Quartic Oscillator
 *  - 2: Quantum mechanical rotor
*/
enum ActionType {
  ActionHarmonicOscillator = 0,
  ActionQuarticOscillator = 1,
  ActionRotor = 2
};

/** @class QMParameters
 *
 * @brief Class for storing general quantum mechanics parameter
 *
 */
class QMParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  QMParameters() :
    Parameters("quantummechanics"),
    action_(ActionHarmonicOscillator) {
    addKey("action",String);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      std::string action_str = getContents("action")->getString();
      if (action_str=="harmonicoscillator") {
        action_ = ActionHarmonicOscillator;
      } else if (action_str=="quarticoscillator") {
        action_ = ActionQuarticOscillator;
      } else if (action_str=="rotor") {
        action_ = ActionRotor;
      } else {
      }
    }
    return readSuccess;
  }

  /** @brief Return the action type */
  ActionType action() const { return action_; }

private:
  /** @brief Type of action */
  ActionType action_;
};


/** @class QMAction
 *
 * @brief Base class for quantum mechanical actions
 *
 * Allows calculation of action
 *   \f[
        S[X]=\sum_{j=0}^{M-1}\left(\frac{m_0}{2}\frac{(X_{j+1}-X_j)}{a^2}+V(x)\right)
      \f]
 * for one-dimensional quantum problem with periodic boundary conditions.
 */
class QMAction : public Action {
public:
  /** @brief Initialise class
   *
   * 
   * @param[in] lattice_ Underlying lattice
   * @param[in] renormalisation_ Type of renormalisation
  
   * @param[in] m0_ Mass of particle \f$m_0\f$
   */
  QMAction(const std::shared_ptr<Lattice1D> lattice_,
           const RenormalisationType renormalisation_,
           const double m0_)
    : Action(lattice_,renormalisation_),
      M_lat(lattice_->getM_lat()),
      a_lat(lattice_->geta_lat()),
      m0(m0_) {
    assert(m0>0.0);
  }
  
  /** @brief return size of samples */
  virtual unsigned int sample_size() { return M_lat; }

  /** @brief Return underlying lattice */
  std::shared_ptr<Lattice1D> get_lattice() const { return lattice; }

  /** @brief Return mass \f$m_0\f$ */
  double getm0() const { return m0;}

  /** @brief Cost of one action evaluation */
  virtual double evaluation_cost() const { return lattice->getM_lat(); }

  /** @brief Construct coarsened version of action
   *
   * This returns a coarsened version of the action on the next level
   * of the multigrid hierarchy.
   */
  virtual std::shared_ptr<Action> coarse_action() {
    mpi_parallel::cerr << "ERROR: cannot coarsen action" << std::endl;
    mpi_exit(EXIT_FAILURE);
    throw std::runtime_error("...");
  };
    
  /** @brief Evaluate action for a specific path
   * 
   * Calculate \f$S[X]\f$ for a specific path
   *
   * @param[in] x_path Path, has to be an array of length \f$M\f$
   */
  virtual const double evaluate(const std::shared_ptr<SampleState> x_path) const = 0;

  /** @brief Calculate force for HMC integrator for a specific path
   *
   * Calculate \f$P = \frac{\partial S[X]}{\partial X}\f$ for a specific
   * path and return the resulting force as a path.
   *
   * @param[in] x_path Path \f$X\f$ on which to evaluate the force
   * @param[out] p_path Resulting force \f$P\f$ at every point
   *
   */
  virtual void force(const std::shared_ptr<SampleState> x_path,
                     std::shared_ptr<SampleState> p_path) const = 0;

  /** @brief Initialise path 
   *
   * Set initial values of path, those values will be used to start the 
   * sampling process
   *
   * @param[out] x_path Path \f$X\f$ to be set
   */
  virtual void initialise_path(std::shared_ptr<SampleState> x_path) const = 0;
  
  /** @brief Second derivative \f$W''_{x_-,x_+}(x)\f$ of conditioned
   * action at its minimum.
   *
   * The second derivative (=curvature) of the conditioned action
   * \f$W_{x_-,x_+}(x)=S(\dots,x_-,x,x_+,\dots)\f$ at its minimum. In the
   * special case of a Lagrangian of the form \f$\frac{m_0}{2}\dot{x}^2+V(x)\f$
   * this becomes
   \f[ 
     W_{\overline{x}}(x)=\frac{m_0}{2a}\left((x-x_+)^2+(x-x_-)^2\right)+aV(x)
   \f]
   * where \f$\overline{x}=\frac{x_++x_-}{2}\f$.
   * This quantity is required for sampling of the fine lattice sites.
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWcurvature(const double x_m,
                                      const double x_p) const = 0;

  /** @brief Find minimum of conditioned action \f$W_{x_-,x_+}(x)\f$
   *
   * Given \f$x_-\f$ and \f$x_+\f$, find the minimum \f$x_0\f$ 
   * of the conditioned action \f$W_{\overline{x}}(x)\f$
   *
   * @param[in] x_m Value of \f$x_-\f$
   * @param[in] x_p Value of \f$x_+\f$
   */
  double virtual inline getWminimum(const double x_m,
                                    const double x_p) const = 0;
    
  /** @brief Copy coarse data points from sample on coarser level
   *
   * @param[in] x_coarse Coarse sample to copy from
   * @param[in] x_path Fine path to copy to (same level as action)
   */
  virtual void copy_from_coarse(const std::shared_ptr<SampleState> x_coarse,
                                std::shared_ptr<SampleState> x_path);
  
  /** @brief Copy coarse data points from path on finer level
   *
   * @param[in] x_fine Fine path to copy from
   * @param[in] x_path Coarse path to copy to (same level as action)
   */
  virtual void copy_from_fine(const std::shared_ptr<SampleState> x_fine,
                              std::shared_ptr<SampleState> x_path);

protected:
  /** @brief Number of lattice points */
  const unsigned int M_lat;
  /** @brief Lattice spacing */
  const double a_lat;
  /** @brief Particle mass */
  const double m0;
};

#endif // QMACTION_HH
