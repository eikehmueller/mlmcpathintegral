#ifndef MPI_RANDOM
#define MPI_RANDOM MPI_RANDOM

#include <set>
#include <random>
#include <algorithm>
#include <type_traits>
#include "mpi/mpi_wrapper.hh"

/** @file mpi_random.hh
 * @brief Header file for parallel random number generation
 */

/** @class parallel_mt19937_64
 *
 * @brief Parallel Mersenne Twister engine
 *
 * Random number engine. An independent std::mt19937_64 engine is run by each
 * process, but the seeds differ. The initial seeds are obtained by generating
 * random numbers with a std::midstd_rand engine and scattering them to all
 * processors.
 * Note that the class is derived from the std::mt19937_64 class, and has the
 * same interface as the base class, i.e. random numbers are generated by
 * calling the operator() method.
 */
class parallel_mt19937_64 : public std::mt19937_64 {
public:
    /** @brief Constructor */
    parallel_mt19937_64(result_type value = default_seed) { seed(value); };
    
    /** @brief Seed random number generator
     *
     * @param[in] value Seed to be used on master process.
     */
    void seed(result_type value = default_seed);

    /** @brief Return next random number */
    result_type operator()() {
        return base_engine();
    }

    /** @brief Advances the internal state a number of times
     *
     * @param[in] z Number of random numbers to discard
     */
    void discard( unsigned long long z ) {
        base_engine.discard(z);
    }


protected:
    std::mt19937_64 base_engine;
};

namespace mpi_parallel {
#ifdef USE_MPI
typedef parallel_mt19937_64 mt19937_64;
#else
typedef std::mt19937_64 mt19937_64;
#endif
}

#endif // MPI_RANDOM
