#ifndef LATTICE_HH
#define LATTICE_HH LATTICE_HH
#include "mpi/mpi_wrapper.hh"
#include <cassert>
#include <memory>
#include <stdexcept>

/** @file lattice.hh
 * @brief Header file for abstract lattice base class
 */

/** @class Lattice
 *
 * @brief Abstract lattice base class
 *
 */
class Lattice {
public:
  /** @brief Initialise class
   *
   * @param[in] coarsening_level_ Coarsening level (0=finest)
   * @param[in] dimension_ Lattice dimension
   */
  Lattice(const int coarsening_level_ = 0, const int dimension_ = -1)
      : coarsening_level(coarsening_level_), dimension(dimension_) {}

  /** @brief Return coarsening level of action */
  const int get_coarsening_level() const { return coarsening_level; }

  /** @brief Return the number of vertices */
  virtual const unsigned int getNvertices() const = 0;

  /** @brief return list of neighbour vertices
   *
   * Each entry is a list which contains the linear indices of the direct
   * neighbour vertices.
   */
  const std::vector<std::vector<unsigned int>> &get_neighbour_vertices() {
    return neighbour_vertices;
  }

  /** Lattice dimension */
  const int dimension;

protected:
  /** @brief coarsening level (0=finest level) */
  mutable int coarsening_level;
  /** @brief list of direct neighbour vertices */
  std::vector<std::vector<unsigned int>> neighbour_vertices;
};

#endif // LATTICE_HH
