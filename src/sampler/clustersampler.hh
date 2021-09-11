#ifndef CLUSTERSAMPLER_HH
#define CLUSTERSAMPLER_HH CLUSTERSAMPLER_HH
#include <random>
#include <vector>
#include <memory>
#include "common/parameters.hh"
#include "common/timer.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_random.hh"
#include "lattice/lattice.hh"
#include "action/action.hh"
#include "action/qm/qmclusteraction.hh"
#include "sampler/sampler.hh"

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
    unsigned int n_burnin() const {
        return n_burnin_;
    }
    /** @brief Return number of updates between steps */
    unsigned int n_updates() const {
        return n_updates_;
    }
private:
    /** @brief Number of burnin samples */
    unsigned int n_burnin_;
    /** @brief Number of cluster updates between steps */
    unsigned int n_updates_;
};

/** @class QMClusterSampler
 * @brief Cluster algorithm sampler for quantum mechanical systems
 *
 * Generates samples by using the cluster algorithm described in
 * Wolff, U., 1989. "Collective Monte Carlo updating for spin systems".
 * Physical Review Letters, 62(4), p.361. (see also
 * <a href="https://arxiv.org/abs/hep-lat/9704009">arXiv/hep-lat/9704009</a>).
 */
class QMClusterSampler : public Sampler {
public:
    /** @brief Create new instance
     *
     * @param[in] action_ Action to sample from
     * @param[in] cluster_param Parameters of cluster sample
     */
    QMClusterSampler(const std::shared_ptr<QMClusterAction> action_,
                     const ClusterParameters cluster_param);

    /** @brief Destroy instance
     *
     * Deallocate memory
     */
    virtual ~QMClusterSampler() {}

    /** @brief Draw a sample
     *
     * returns a sample path \f$X\f$
     *
     * @param[out] x_path Path \f$X\f$ drawn from distribution
     */
    virtual void draw(std::shared_ptr<SampleState> x_path);

    /** Return cost per sample */
    virtual double cost_per_sample() {
        return cost_per_sample_;
    }

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
    /** @brief Set current state to particular value
     *
     * @param[in] x_path
     */
    virtual void set_state(std::shared_ptr<SampleState> x_path) {};

protected:
    /** @brief Action to sample from */
    const std::shared_ptr<QMClusterAction> action;
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice> lattice;
    /** @brief Number of vertices on lattice */
    const unsigned int n_vertices;
    /** @brief Number of burn-in steps */
    const unsigned int n_burnin;
    /** @brief Number of updates per step */
    const unsigned int n_updates;
    /** @brief Current state (path) */
    mutable std::shared_ptr<SampleState> x_path_cur;
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
    /** @brief cost per sample */
    double cost_per_sample_;
};

class QMClusterSamplerFactory : public SamplerFactory {
public:
    /** @brief Create new instance
     *
     * @param[in] param_cluster Custer sampler parameters
     */
    QMClusterSamplerFactory(const ClusterParameters param_cluster_) :
        param_cluster(param_cluster_) {}

    /** @brief Destructor */
    virtual ~QMClusterSamplerFactory() {}

    /** @brief Return sampler for a specific  action
     *
     * @param[in] action Action to sample from
     */
    virtual std::shared_ptr<Sampler> get(std::shared_ptr<Action> action) {
        return std::make_shared<QMClusterSampler>(std::dynamic_pointer_cast<QMClusterAction>(action),
                                                  param_cluster);
    }
private:
    /** Cluster sampler parameters */
    const ClusterParameters param_cluster;
};

#endif // CLUSTERSAMPLER_HH
