#ifndef CLUSTERSAMPLER_HH
#define CLUSTERSAMPLER_HH CLUSTERSAMPLER_HH
#include <random>
#include <vector>
#include <memory>
#include "path.hh"
#include "clusteraction.hh"
#include "parameters.hh"
#include "sampler.hh"
#include "mpi_random.hh"

/** @file clustersampler.hh
 * @brief Header file for sampler based on cluster algorithm
 */

/** @class ClusterParameters
 *
 * @brief Class for storing parameters of Cluster algorithm
*/
class ClusterParameters : public Parameters {
public:
  /** @brief Construct a new instance */
  ClusterParameters() :
    Parameters("clusteralgorithm"),
    n_burnin_(100),
    n_updates_(10) {
    addKey("n_burnin",Integer,Positive);
    addKey("n_updates",Integer,Positive);
  }

  /** @brief Read parameters from file
   *
   * @param[in] filename Name of file to read
   */
  int readFile(const std::string filename) {

    int readSuccess = Parameters::readFile(filename);
    if (!readSuccess) {
      n_burnin_ = getContents("n_burnin")->getInt();
      n_updates_ = getContents("n_updates")->getInt();
    }
    return readSuccess;
  }
  /** @brief Return number of burnin samples */
  unsigned int n_burnin() const { return n_burnin_; }
  /** @brief Return number of updates between steps */
  unsigned int n_updates() const { return n_updates_; }
private:
  /** @brief Number of burnin samples */
  unsigned int n_burnin_;
  /** @brief Number of cluster updates between steps */
  unsigned int n_updates_;
};

/** @class ClusterSampler
 * @brief Cluster algorithm sampler
 *
 * Generates samples by using the cluster algorithm described in
 * Wolff, U., 1989. "Collective Monte Carlo updating for spin systems".
 * Physical Review Letters, 62(4), p.361. (see also
 * <a href="https://arxiv.org/abs/hep-lat/9704009">arXiv/hep-lat/9704009</a>). 
 */
class ClusterSampler : public Sampler {
public:
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] n_burnin_ Number of burnin steps
   * @param[in] n_updates_ Number of cluster updates between steps
   */
  ClusterSampler(const std::shared_ptr<ClusterAction> action_,
                 const unsigned int n_burnin_,
                 const unsigned int n_updates_) :
    Sampler(),
    action(action_),
    n_burnin(n_burnin_),
    n_updates(n_updates_),
    uniform_dist(0.0,1.0),
    uniform_int_dist(0,action_->getM_lat()-1)
  {
    engine.seed(2141517);
    const unsigned int M_lat = action->getM_lat();
    const double T_final = action->getT_final();
    // Create temporary workspace
    x_path_cur = std::make_shared<Path>(M_lat,T_final);
    action->initialise_path(x_path_cur);
    // Burn in
    std::shared_ptr<Path> x_path_tmp = std::make_shared<Path>(M_lat,T_final);
    for (unsigned int i=0;i<n_burnin;++i) {
      draw(x_path_tmp);
    }
  }

  /** @brief Destroy instance
   *
   * Deallocate memory
   */
  virtual ~ClusterSampler() {}

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::shared_ptr<Path> x_path);

private:
  
/** @brief Process next link
 *
 * Depending on direction, the neighbour is either site \f$i+1\f$ or 
 * \f$i-1\f$. Connect the sites with the probability
 * \f$1-e^{\min(0,-S_{\ell})}\f$ and return a tuple containing the bond
 * value and the index of the neighbouring site.
 *
 * @param[in] Lattice site 
 * @param[in] direction Direction to neighbour (has to be +1 or -1)
 */
  std::pair<bool,int> process_link(const int i,
                                   const int direction);
  
protected:
  /** @brief Action to sample from */
  const std::shared_ptr<ClusterAction> action;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
  /** @brief Number of updates per step */
  const unsigned int n_updates;
  /** @brief Current state (path) */
  mutable std::shared_ptr<Path> x_path_cur;
  /** @brief Random number engine */
  typedef mpi_parallel::mt19937_64 Engine;
  /** @brief Type of Mersenne twister engine */
  mutable Engine engine;
  /** @brief Type of uniform distribution */
  typedef std::uniform_real_distribution<double> Uniform;
  /** @brief Uniform distribution used for setting bonds  */
  mutable Uniform uniform_dist;
  /** @brief Type for integer uniform distribution */
  typedef std::uniform_int_distribution<unsigned int> UniformInt;
  /** @brief Uniform int distribution for picking first site */
  mutable UniformInt uniform_int_dist;
};

#endif // CLUSTERSAMPLER_HH
