#ifndef PATH_HH
#define PATH_HH PATH_HH
#include <fstream>
#include <algorithm>
#include <memory>
#include <cassert>

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
   * @param[in] M_lat_ Number of lattice sites
   * @param[in] T_final_ Final time
   */
  Path(const unsigned int M_lat_,
       const double T_final_) : M_lat(M_lat_), 
                                T_final(T_final_) {
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
  
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief Data array */
  double* data;
};

#endif // PATH_HH
