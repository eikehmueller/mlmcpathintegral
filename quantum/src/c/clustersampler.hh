#ifndef CLUSTERSAMPLER_HH
#define CLUSTERSAMPLER_HH CLUSTERSAMPLER_HH
#include "path.hh"
#include "clusteraction.hh"
#include <random>
#include <vector>

/** @file clustersampler.hh
 * @brief Header file for sampler based on cluster algorithm
 */

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
  /** @brief Base class */
  typedef Sampler Base;
  /** @brief Create new instance
   *
   * @param[in] action_ Action to sample from
   * @param[in] n_burnin_ Number of burnin steps
   */
  ClusterSampler(const ClusterAction& action_,
                 const unsigned int n_burnin_) :
    Base(action_.getM_lat()),
    action(action_),
    n_burnin(n_burnin_),
    uniform_dist(0.0,1.0),
    uniform_int_dist(0,action_.getM_lat())
  {
    engine.seed(2141517);
    const unsigned int M_lat = action.getM_lat();
    const double T_final = action.getT_final();
    // Create temporary workspace
    x_path_cur = std::make_shared<Path>(M_lat,T_final);
    action.initialise_path(x_path_cur);
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
  virtual ~ClusterSampler() {}

  /** @brief Draw a sample 
   *
   * returns a sample path \f$X\f$
   *
   * @param[out] x_path Path \f$X\f$ drawn from distribution
   */
  virtual void draw(std::vector<std::shared_ptr<Path>> x_path);

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
  const ClusterAction& action;
  /** @brief Number of burn-in steps */
  const unsigned int n_burnin;
  /** @brief Current state (path) */
  mutable std::shared_ptr<Path> x_path_cur;
  /** @brief Random number engine */
  typedef std::mt19937_64 Engine;
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
