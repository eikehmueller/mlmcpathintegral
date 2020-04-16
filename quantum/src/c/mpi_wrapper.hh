#ifndef MPI_WRAPPER_HH
#define MPI_WRAPPER_HH MPI_WRAPPER_HH

#include <string>
#include <iostream>
#include <ostream>
#include <vector>
#include <cmath>
#include <stdlib.h>

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

/** @brief Initialise MPI */
void mpi_init();

/** @brief Finalise MPI */
void mpi_finalize();

/** @brief Return MPI rank */
int mpi_comm_rank();

/** @brief Return total number of MPI ranks */
int mpi_comm_size();

/** @brief Return true on MPI master */
bool mpi_master();

/** @brief Parallel sum for scalar valued quantities
 * 
 * @param x quantity to sum across all processes
 */ 
double mpi_allreduce_sum(const double x);

/** @brief Parallel average for scalar valued quantities
 * 
 * @param x quantity to average across all processes
 */ 
double mpi_allreduce_avg(const double x);

/** @brief Parallel average for vector valued quantities
 * 
 * @param x quantity to average across all processes
 */
std::vector<double> mpi_allreduce_avg(const std::vector<double> x);

/** @brief Parallel minimum for scalar valued quantities
 * 
 * @param x quantity to take minimum of across all processes
 */ 
int mpi_allreduce_min(const int x);

/** @brief Parallel sum for scalar valued quantities
 * 
 * @param x quantity to sum across all processes
 */ 
int mpi_allreduce_sum(const int x);

/** @brief Parallel sum for scalar valued quantities
 * 
 * @param x quantity to sum across all processes
 */ 
unsigned int mpi_allreduce_sum(const unsigned int x);

/** @brief Parallel logical and for scalar valued quantities
 * 
 * @param x quantity to take logical and of across all processes
 */ 
bool mpi_allreduce_and(const bool x);

/** @brief Broadcast double value to all processes 
 *
 * @param x value to broadcast
 * @param root processor to broadcast from
 */
void mpi_bcast(double& x,const int root=0);

/** @brief Broadcast integer value to all processes 
 *
 * @param x value to broadcast
 * @param root processor to broadcast from
 */
void mpi_bcast(int& x,const int root=0);

/** @brief Broadcast bool value to all processes 
 *
 * @param x value to broadcast
 * @param root processor to broadcast from
 */
void mpi_bcast(bool& x,const int root=0);

/** @brief Broadcast string value to all processes 
 *
 * @param x value to broadcast
 * @param root processor to broadcast from
 */
void mpi_bcast(std::string& x,const int root=0);

/** @brief Scatter values to all processes 
 * 
 * @param[in] data List to scatter
 * @param[out] x Resulting values of all processes
 * 
 * This takes the list [x[0],x[1],x[2],...,x[nproc-1]] and scatters it.
 */
void mpi_scatter(unsigned int* data, unsigned int& x);

/** @brief Call mpi_finalize and exit
 *
 * @param[in] exit_code Exit code
 */
void mpi_exit(const int exit_code);

/** @brief Barrier */
void mpi_barrier();

/** @brief distribute a number evenly between processors
 *
 * Split the number n as evenly as possible between processors, i.e.
 * return n_p on processor p, such that n_0 + n_1 + ... + n_{nproc-1} = n.
 *
 * @param[in] n Number to distribute
 */
unsigned int distribute_n(const unsigned int n);

static bool mpi_initialised=false;

/* The following code implements output stream that can be used in parallel 
 * such that the output is only printed to std::cout/std::cerr on the 
 * master process.
 */
namespace mpi_parallel {

  /** @class MPIMasterStream 
   * 
   * @brief Class for wrapping output stream to print on MPI master only
   *
   * Based on the ideas for the LoggedStream class described at
   * https://stackoverflow.com/questions/772355/how-to-inherit-from-stdostream
   */
  class MPIMasterStream {
  public:
    /** @brief Constructor 
     * 
     * @param[in] out_ Stream to wrap
     */
    MPIMasterStream(std::ostream& out_) : is_dummy(true) {
      mpi_init();
      if (not mpi_master()) {
        out = new(std::ostream)(nullptr);
        out->setstate(std::ios_base::badbit);
      } else {
        is_dummy = false;
        out = &out_;
      }
    }

    /** @brief Destructor, free memory if it has been allocated */
    ~MPIMasterStream() {
      if (is_dummy) {
        delete(out);
      }
    }
    
    template <typename T>
    /** @brief Output only on MPI master
     *
     * @param[in] v Object to pass to << operator
     */
    const MPIMasterStream& operator<<(const T& v) const {
      *out << v;
      return *this;
    }

    /** @brief Output only on MPI master
     *
     * @param[in] F Functor to use
     */
    MPIMasterStream const& operator<<(std::ostream& (*F)(std::ostream&)) const
    {
      F(*out);
      return *this;
    }
    
  protected:
    std::ostream* out; /** @brief Wrapped output stream */
    bool is_dummy; /** @brief Set to true if not master process */ 
  };

  /** @brief wrapped std::cout object */
  static MPIMasterStream cout(std::cout);
  /** @brief wrapped std::cerr object */
  static MPIMasterStream cerr(std::cerr);
}

#endif // MPI_HH
