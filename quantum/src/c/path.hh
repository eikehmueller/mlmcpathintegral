#ifndef PATH_HH
#define PATH_HH PATH_HH

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
   */
  Path(const unsigned int M_lat_) : M_lat(M_lat_) {
    data = new double[M_lat];
  }

  /** @brief Destroy instance
   * 
   * Deallocate memory
   */
  ~Path() { delete[] data; }
  /** @brief Number of time slices */
  const unsigned int M_lat;
  /** @brief Data array */
  double* data;
};

#endif // PATH_HH