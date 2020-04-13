#ifndef MPI_WRAPPER_HH
#define MPI_WRAPPER_HH MPI_WRAPPER_HH

#define USE_MPI 1

#include <string>
#include <iostream>
#include <ostream>

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
double mpi_allreduce_sum_scalar(const double x);

/** @brief Parallel sum for scalar valued quantities
 * 
 * @param x quantity to sum across all processes
 */ 
int mpi_allreduce_sum_scalar(const int x);

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
    MPIMasterStream(std::ostream& out_) : out(out_) {}
    
    template <typename T>
    /** @brief Output only on MPI master
     *
     * @param[in] v Object to pass to << operator
     */
    const MPIMasterStream& operator<<(const T& v) const {
      if (mpi_master()) {
        out << v;
      }
      return *this;
    }

    /** @brief Output only on MPI master
     *
     * @param[in] F Functor to use
     */
    MPIMasterStream const& operator<<(std::ostream& (*F)(std::ostream&)) const
    {
      if (mpi_master()) {
        F(out);
      }
      return *this;
    }
    
  protected:
    std::ostream& out; /** @brief Wrapped output stream */
  };

  /** @brief wrapped std::cout object */
  static MPIMasterStream cout(std::cout);
  /** @brief wrapped std::cerr object */
  static MPIMasterStream cerr(std::cerr);
}

#endif // MPI_HH
