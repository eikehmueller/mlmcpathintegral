#ifndef SAMPLESTATE_HH
#define SAMPLESTATE_HH SAMPLESTATE_HH
#include <fstream>
#include <algorithm>
#include <memory>
#include <cassert>
#include <Eigen/Dense>

/** @file samplestate.hh
 * @brief Header file for sample state class
 */

/** @class SampleState
 *
 * @brief state that can be used in the sampler- and MC algorithms
 *
 * Ligthweight wrapper around Eigen-vector
 */
class SampleState {
public:
  /** @brief Create new instance
   *
   * Allocate memory
   * 
   * @param[in] M_ number of unknowns
   */
  SampleState(const unsigned int M_) : M(M_), data(M_) {
    for (unsigned int j=0;j<M;++j) data[j] = 0.0;
  }

  /** @brief Destroy instance
   */
  ~SampleState() {}

  /** @brief Save path to disk
   *
   * Save path to disk by writing a file which contains the position at
   * different times.
   * First line contains \f$M\f$, the second line contains the
   * positions separated by spaces.
   *
   * @param[in] filename Name of file to save to
   */
  void save_to_disk(const std::string filename);
    
  /** @brief Data array */
  Eigen::VectorXd data;

private:
  /** @brief Size of vector */
  const unsigned int M;
};

#endif // SAMPLEPATH_HH
