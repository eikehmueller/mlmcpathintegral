#ifndef PATH_HH
#define PATH_HH PATH_HH
#include <fstream>
#include <algorithm>
#include <memory>
#include <cassert>
#include "lattice/lattice1d.hh"

/** @file path.hh
 * @brief Header file for path class
 */

/** @class Path
 *
 * @brief Path on \f$M\f$ timeslices
 *
 * Ligthweight wrapper around vector of length \f$M\f$
 */
class Path {
public:
  /** @brief Create new instance
   *
   * Allocate memory
   * 
   * @param[in] lattice_ Underlying lattice
   */
  Path(const std::shared_ptr<Lattice1D> lattice_) :
    lattice(lattice_),
    M_lat(lattice_->getM_lat()) {
    data = new double[M_lat];
    std::fill(data,data+M_lat,0.0);
  }

  /** @brief Destroy instance
   * 
   * Deallocate memory
   */
  ~Path() { delete[] data; }

  /** @brief Save path to disk
   *
   * Save path to disk by writing a file which contains the position at
   * different times.
   * First line contains \f$M\f$ and \f$T\f$, the second line contains the
   * positions separated by spaces.
   *
   * @param[in] filename Name of file to save to
   */
  void save_to_disk(const std::string filename);

  /** @brief Axpy operation
   *
   * Add \f$\alpha x\f$ to path, i.e computer \f$y\mapsto y+\alpha x\f$
   *
   * @param[in] alpha Scaling factor \f$\alpha\f$
   * @param[in] other Other path \f$x\f$
   */
  void axpy(const double alpha, const std::shared_ptr<Path> other) {
    for(unsigned int j=0;j<M_lat;++j) {
      data[j] += alpha*other->data[j];
    }
  }
  
  /** @brief Squared Euclidean norm
   */
  double norm2() const {
    double nrm2 = 0.0;
    for(unsigned int j=0;j<M_lat;++j) {
      double tmp=data[j];
      nrm2 += tmp*tmp;
    }
    return nrm2;
  }
  
  /** @brief Fill with numbers drawn from lambda function
   *
   * @param[in] f Function to draw from
   */
  template <class Lambda>
  void fill(Lambda&& f) {
    for(unsigned int j=0;j<M_lat;++j) {
      data[j] = f();
    }
  }
  
  /** @brief Copy data from other path
   * @param[in] other Path to copy from
   */
  void copy(const std::shared_ptr<Path> other);

  /** @brief Copy coarse data points from path on coarser level
    *
   * @param[in] other Coarse path to copy from
   */
  void copy_from_coarse(const std::shared_ptr<Path> other);
  
  /** @brief Copy coarse data points from path on finer level
    *
   * @param[in] other Fine path to copy from
   */
  void copy_from_fine(const std::shared_ptr<Path> other);
  
  /** @brief Return underlying lattice */
  std::shared_ptr<Lattice1D> get_lattice() { return lattice; }

  /** @brief Data array */
  double* data;

private:
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Underlying lattice */
  const std::shared_ptr<Lattice1D> lattice;
};

#endif // PATH_HH
