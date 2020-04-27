#ifndef PATH_HH
#define PATH_HH PATH_HH
#include <fstream>
#include <algorithm>

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

  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Final time */
  const double T_final;
  /** @brief Data array */
  double* data;
};

#endif // PATH_HH
