#ifndef LATTICE1D_HH
#define LATTICE1D_HH LATTICE1D_HH
#include <memory>
#include <stdexcept>
#include <cassert>
#include "mpi/mpi_wrapper.hh"
#include "common/parameters.hh"
#include "lattice/lattice.hh"

/** @file lattice1d.hh
 * @brief Header file for one-dimensional lattice class
 */

/** @class Lattice1DParameters
 *
 * @brief Class for storing parameters of lattice
 *
 * This stores the number \f$M_{lat}\f$ of lattice sites and the
 * final time \f$T_{final}\f$
 */
class Lattice1DParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    Lattice1DParameters() :
        Parameters("lattice"),
        M_lat_(8),
        T_final_(1.0) {
        addKey("M_lat",Integer,Positive);
        addKey("T_final",Double,Positive);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            M_lat_ = getContents("M_lat")->getInt();
            T_final_ = getContents("T_final")->getDouble();
        }
        return readSuccess;
    }

    /** @brief Return number of lattice sites */
    unsigned int M_lat() const {
        return M_lat_;
    }
    /** @brief Return final time */
    double T_final() const {
        return T_final_;
    }
        
private:
    /** @brief Number of lattice sites */
    unsigned int M_lat_;
    /** @brief Final time */
    double T_final_;
};

/** @class Lattice1D
 *
 * @brief Class for 1d lattice
 *
 */
class Lattice1D : public Lattice {
public:
    /** @brief Initialise class
     *
     *
     * @param[in] M_lat_ Number of time slices \f$M\f$
     * @param[in] T_final_ Final time \f$T\f$
     */
    Lattice1D(const unsigned int M_lat_,
              const double T_final_,
              const int coarsening_level_=0);
    
    /** @brief Return number of timeslices \f$M\f$ */
    unsigned int getM_lat() const {
        return M_lat;
    }

    /** @return final time \f$T\f$ */
    double getT_final() const {
        return T_final;
    }

    /** @brief Return lattice spacing \f$a\f$ */
    double geta_lat() const {
        return a_lat;
    }

    virtual std::shared_ptr<Lattice1D> coarse_lattice() {
        if (M_lat%2) {
            mpi_parallel::cerr << "ERROR: cannot coarsen 1d lattice with M = " << M_lat << " points." << std::endl;
            mpi_exit(EXIT_FAILURE);
            throw std::runtime_error("...");
        }
        return std::make_shared<Lattice1D>(M_lat/2,T_final,coarsening_level+1);
    };
    
    /** @brief Return the number of vertices */
    virtual const unsigned int getNvertices() const { return M_lat; }
    
protected:
    /** @brief Number of time slices */
    const unsigned int M_lat;
    /** @brief Final time */
    const double T_final;
    /** @brief lattice spacing */
    const double a_lat;
};

#endif // LATTICE1D_HH
