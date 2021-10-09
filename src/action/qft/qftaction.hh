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

/** @brief Enum for different actions
 *  - 0: Quenched Schwinger model
 *  - 1: Nonlinear sigma model
 *  - 2: Gaussian Free Field (GFF)
 */
enum QFTActionType {
    ActionQuenchedSchwinger = 0,
    ActionNonlinearSigma = 1,
    ActionGFF = 2
};

/** @class QFTParameters
 *
 * @brief Class for storing general quantum field theory parameters
 *
 */
class QFTParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    QFTParameters() :
        Parameters("quantumfieldtheory"),
        action_(ActionQuenchedSchwinger) {
        addKey("action",String);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            std::string action_str = getContents("action")->getString();
            if (action_str=="quenchedschwinger") {
                action_ = ActionQuenchedSchwinger;
            } else if (action_str=="nonlinearsigma") {
                action_ = ActionNonlinearSigma;
            } else if (action_str=="gff") {
                action_ = ActionGFF;
            } else {
            }
        }
        return readSuccess;
    }

    /** @brief Return the action type */
    QFTActionType action() const {
        return action_;
    }

private:
    /** @brief Type of action */
    QFTActionType action_;
};

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
     * @param[in] fine_lattice_ Two-dimensional lattice on next finer level
     * @param[in] renormalisation_ Type of renormalisation
     */
    QFTAction(const std::shared_ptr<Lattice2D> lattice_,
              const std::shared_ptr<Lattice2D> fine_lattice_,
              const RenormalisationType renormalisation_)
        : Action(renormalisation_),
          lattice(lattice_),
          fine_lattice(fine_lattice_),
          coarse_lattice(lattice->get_coarse_lattice()) { }
    
    /** @brief Get underlying lattice */
    std::shared_ptr<Lattice2D> get_lattice() {
        return lattice;
    }
    
    /** @brief Get underlying coarse lattice */
    std::shared_ptr<Lattice2D> get_coarse_lattice() {
        return coarse_lattice;
    }
    
    /** @brief Get coarsening level
     *
     * This will return the coarsening level of the underlying lattice */
    virtual int get_coarsening_level() const {
        return lattice->get_coarsening_level();
    }
        
    /** @brief Action information string
     *
     * return some information on this instance of the action
     */
    virtual std::string info_string() const;

protected:
    /** @brief Underlying lattice */
    const std::shared_ptr<Lattice2D> lattice;
    /** @brief Underlying refined lattice */
    const std::shared_ptr<Lattice2D> fine_lattice;
    /** @brief Underlying coarsened lattice */
    const std::shared_ptr<Lattice2D> coarse_lattice;
};

#endif // QFTACTION_HH
