#ifndef LATTICE2D_HH
#define LATTICE2D_HH LATTICE2D_HH
#include <memory>
#include <stdexcept>
#include <cassert>
#include "mpi/mpi_wrapper.hh"
#include "common/parameters.hh"
#include "lattice/lattice.hh"

/** @file lattice2d.hh
 * @brief Header file for two-dimensional lattice class
 */

/** @class Lattice2DParameters
 *
 * @brief Class for storing parameters of a two-dimensional lattice
 *
 * This stores the number \f$M_{t,lat}\f$, \f$M_{x,lat}\f$ of lattice sites in the
 * temporal and spatial direction, as well as the correponding extents of the lattice
 * \f$T_{lat}\f$ and \f$L_{lat}.\f$
 */
class Lattice2DParameters : public Parameters {
public:
    /** @brief Construct a new instance */
    Lattice2DParameters() :
        Parameters("lattice"),
        Mt_lat_(8),
        Mx_lat_(8),
        T_lat_(1.0),
        L_lat_(1.0) {
        addKey("Mt_lat",Integer,Positive);
        addKey("Mx_lat",Integer,Positive);
        addKey("T_lat",Double,Positive);
        addKey("L_lat",Double,Positive);
    }

    /** @brief Read parameters from file
     *
     * @param[in] filename Name of file to read
     */
    int readFile(const std::string filename) {

        int readSuccess = Parameters::readFile(filename);
        if (!readSuccess) {
            Mt_lat_ = getContents("Mt_lat")->getInt();
            Mx_lat_ = getContents("Mx_lat")->getInt();
            T_lat_ = getContents("T_lat")->getDouble();
            L_lat_ = getContents("L_lat")->getDouble();
        }
        return readSuccess;
    }

    /** @brief Return number of lattice points in temporal direction  */
    unsigned int Mt_lat() const {
        return Mt_lat_;
    }

    /** @brief Return number of lattice points in spatial direction  */
    unsigned int Mx_lat() const {
        return Mx_lat_;
    }

    /** @brief Return temporal extent of lattice */
    double T_lat() const {
        return T_lat_;
    }

    /** @brief Return spatial extent of lattice */
    double L_lat() const {
        return L_lat_;
    }

    
private:
    /** @brief Number of lattice points in temporal direction */
    unsigned int Mt_lat_;
    /** @brief Number of lattice points in spatial direction */
    unsigned int Mx_lat_;
    /** @brief Temporal extent of lattice */
    double T_lat_;
    /** @brief Spatial extent of lattice */
    double L_lat_;
};

/** @class Lattice2D
 *
 * @brief Class for 2d lattice
 *
 */
class Lattice2D : public Lattice {
public:
    /** @brief Initialise class
     *
     *
     * @param[in] Mt_lat Number of time slices
     * @param[in] Mx_lat Number of points in spatial direction
     * @param[in] T_lat Temporal extent of lattice
     * @param[in] L_lat Spatial extent of lattice
     * @param[in] coarsening_level_ Coarsening level
     */
    Lattice2D(const unsigned int Mt_lat_,
              const unsigned int Mx_lat_,
              const double T_lat_,
              const double L_lat_,
              const int coarsening_level_=0) :
        Lattice(coarsening_level_),
        Mt_lat(Mt_lat_),
        Mx_lat(Mx_lat_),
        T_lat(T_lat_),
        L_lat(L_lat_),
        at_lat(T_lat_/Mt_lat_),
        ax_lat(L_lat_/Mx_lat_) {}

    /** @brief Return number of timeslices \f$M_{t,lat}\f$ */
    unsigned int getMt_lat() const {
        return Mt_lat;
    }

    /** @brief Return number of points in spatial direction \f$M_{x,lat}\f$ */
    unsigned int getMx_lat() const {
        return Mx_lat;
    }

    /** @return Return temporal extent of lattice  \f$T_{lat}\f$ */
    double getT_lat() const {
        return T_lat;
    }

    /** @return Return spatial extent of lattice  \f$L_{lat}\f$ */
    double getL_lat() const {
        return L_lat;
    }
    
    /** @brief Return temporal lattice spacing \f$a_t\f$ */
    double getat_lat() const {
        return at_lat;
    }

    /** @brief Return spatial lattice spacing \f$a_x\f$ */
    double getax_lat() const {
        return ax_lat;
    }
    
    /** @brief Return number of edges (assuming periodic boundaries) */
    unsigned int getNedges() const {
        return 2*Mt_lat*Mx_lat;
    }

    /** @brief Return number of vertices (assuming periodic boundaries) */
    unsigned int getNvertices() const {
        return Mt_lat*Mx_lat;
    }

    /** @brief Return number of cells (assuming periodic boundaries) */
    unsigned int getNcells() const {
        return Mt_lat*Mx_lat;
    }
    
    /** @brief Convert cartesian lattice index of link to linear index
     *
     * Given a lattice site \f$n=(i,j)\f$ and a direction \f$\mu\f$, work out the
     * corresponding linear index \f$\ell\f$ of the link starting at \f$n\f$
     * and pointing in direction \f$\mu\f$. The links are arranged such that
     * \f[
     *   \ell = 2M_{t,lat}j + 2i + \mu
     * \f]
     *
     * ...             ...            ...            ...           ...
     *  |              |              |              |             |
     *  9           11           13           15          (9)
     *  |              |              |              |             |
     *  o--- 8 ---o--- 10 --o--- 12 --o-- 14--o
     *  |              |              |              |             |
     *  1            3             5             7           (1)
     *  |              |              |              |             |
     *  o--- 0 ---o--- 2 ---o--- 4 ---o--- 6 ---o
     *
     * temporal direction (index i) = horizontal
     * spatial direction (index j) = vertical
     * 
     * NOTE THAT j + Mx_lat AND i + Mx_lat HAVE TO BE POSITIVE. This is to ensure
     * unpredictable behaviour when computing a % b for negative values of a,
     *
     * @param[in] i Position index in temporal direction
     * @param[in] j Position index in spatial direction
     * @param[in] mu Direction \f$\mu\f$, must be 0 (temporal) or 1 (spatial)
     */
    inline unsigned int link_cart2lin(const int i, const int j, const int mu) const {
        return 2*Mt_lat*((j+Mx_lat)%Mx_lat) + 2*((i+Mt_lat)%Mt_lat)+mu;
    }

    /** @brief Convert linear index to lattice index
     *
     * Given \f$\ell = (2 j+\mu) M_{t,lat} + i\f$, work out cartesian index \f$(i,j)\f$
     * of site and direction \f$\mu\f$
     
     * @param[in] ell Linear index \f$\ell\f$
     * @param[out] i Position index in temporal direction
     * @param[out] j Position index in spatial direction
     * @param[out] mu Direction \f$\mu\f$
     */
    inline void link_lin2cart(const unsigned int ell, int& i, int& j, int& mu) const {
        j = ell / (2*Mt_lat); // j = ell / (2*M_{t,lat})
        unsigned int r = ell - (2*Mt_lat)*j;
        i = r >> 1; // i = r/2
        mu = r & 1; // mu = r % 2
    }

    /** @brief Construct coarsened lattice
     *
     * Returns lattice with half the lattice spacing
     */
    virtual std::shared_ptr<Lattice2D> fine_lattice() {
        return std::make_shared<Lattice2D>(2*Mt_lat,2*Mx_lat,T_lat,L_lat,coarsening_level-1);
    };

    /** @brief Construct coarsened lattice
     *
     * @param[in] exit_on_failure Abort with an error message if construction is not possible?
     *
     * Returns lattice with twice the lattice spacing
     */
    virtual std::shared_ptr<Lattice2D> coarse_lattice(const bool exit_on_failure=true) {
        if ( (Mt_lat%2) or (Mx_lat%2) ) {
            if (exit_on_failure ) {
                mpi_parallel::cerr << "ERROR: cannot coarsen 2d lattice with M_{t,lat} = " << Mt_lat;
                mpi_parallel::cerr << " , M_{x,lat} = " << Mx_lat << std::endl;
                mpi_exit(EXIT_FAILURE);
                throw std::runtime_error("...");
            } else {
                return nullptr;
            }
        } else {
            return std::make_shared<Lattice2D>(Mt_lat/2,Mx_lat/2,T_lat,L_lat,coarsening_level+1);
        }
    };

protected:
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of points in spatial direction */
    const unsigned int Mx_lat;
    /** @brief Temporal extent of lattice */
    const double T_lat;
    /** @brief Spatial extent of lattice */
    const double L_lat;
    /** @brief lattice spacing in temporal direction */
    const double at_lat;
    /** @brief lattice spacing in spatial direction */
    const double ax_lat;
};

#endif // LATTICE2D_HH
