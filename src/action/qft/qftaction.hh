#ifndef QFTACTION_HH
#define QFTACTION_HH QFTACTION_HH
#include <memory>
#include <cassert>
#include <random>
#include <vector>
#include <iostream>
#include "common/auxilliary.hh"
#include "common/parameters.hh"
#include "common/samplestate.hh"
#include "mpi/mpi_wrapper.hh"
#include "mpi/mpi_random.hh"
#include "lattice/lattice2d.hh"
#include "action/action.hh"

/** @file qftaction.hh
 * @brief Header file for base class of QFT actions
 */

/** @class QFTAction
 *
 * @brief Base class for action 2d QFT actions
 *
 */

class QFTAction : public Action {
public:
    /** @brief Initialise class
     *
     * @param[in] lattice_ Underlying two-dimensional lattice
     * @param[in] coarsening_type_ Type of lattice coarsening
     * @param[in] renormalisation_ Type of renormalisation
     */
    QFTAction(const std::shared_ptr<Lattice2D> lattice_,
              const CoarseningType coarsening_type_,
              const RenormalisationType renormalisation_)
        : Action(renormalisation_),
          lattice(lattice_),
          fine_lattice(lattice->fine_lattice()),
          coarsening_type(coarsening_type_),
          coarse_lattice(lattice->coarse_lattice(((coarsening_type_==CoarsenBoth) or
                                                  (coarsening_type_==CoarsenTemporal))?2:1,
                                                 ((coarsening_type_==CoarsenBoth) or
                                                  (coarsening_type_==CoarsenSpatial))?2:1,
                                                 false)),
          Mt_lat(lattice->getMt_lat()),
          Mx_lat(lattice->getMx_lat()) { }
    
    /** @brief Get underlying lattice */
    std::shared_ptr<Lattice2D> get_lattice() {
        return lattice;
    }
    
    /** @brief Get coarsening level
     *
     * This will return the coarsening level of the underlying lattice */
    virtual int get_coarsening_level() const {
        return lattice->get_coarsening_level();
    }
    
    /** @brief Return coarsening type */
    CoarseningType get_coarsening_type() const {
        return coarsening_type;
    }

    /** @brief Check whether action supports number of coarsening steps
     *
     * @param[in] n_level Number of additional coarsening steps (can be zero)
     */
    virtual void check_coarsening_is_permitted(const unsigned int n_level);

protected:
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Underlying refined lattice */
    const std::shared_ptr<Lattice2D> fine_lattice;
    /** @brief Underlying coarsened lattice */
    const std::shared_ptr<Lattice2D> coarse_lattice;
    /** @brief Coarsening type for lattice (can be both directions, temporal-only or spatial-only) */
    const CoarseningType coarsening_type;
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of points in spatial direction */
    const unsigned int Mx_lat;

};

#endif // QFTACTION_HH
