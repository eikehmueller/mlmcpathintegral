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

/** type of coarsening */
enum CoarseningType {
    CoarsenUnspecified = -1,   // Unspecified coarsening
    CoarsenBoth = 0,           // Coarsen in both directions
    CoarsenTemporal = 1,       // Coarsen in temporal direction only
    CoarsenSpatial = 2,        // Coarsen in spatial direction only
    CoarsenAlternate = 3,      // Alternate coarsening in both directions
    CoarsenRotate = 4          // Coarsen in both directions while
                               // rotating the lattice
};

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
        coarsening_type_(CoarsenUnspecified) {
        addKey("Mt_lat",Integer,Positive);
        addKey("Mx_lat",Integer,Positive);
        addKey("coarsening",String);
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
            std::string coarsening_str = getContents("coarsening")->getString();
            if (coarsening_str == "both") {
                coarsening_type_ = CoarsenBoth;
            } else if (coarsening_str == "temporal") {
                coarsening_type_ = CoarsenTemporal;
            } else if (coarsening_str == "spatial") {
                coarsening_type_ = CoarsenSpatial;
            } else if (coarsening_str == "alternate") {
                coarsening_type_ = CoarsenAlternate;
            } else if (coarsening_str == "rotate") {
                coarsening_type_ = CoarsenRotate;
            }
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

    /** @brief Return coarsening type */
    CoarseningType coarsening_type() const {
        return coarsening_type_;
    }
    
private:
    /** @brief Number of lattice points in temporal direction */
    unsigned int Mt_lat_;
    /** @brief Number of lattice points in spatial direction */
    unsigned int Mx_lat_;
    /** @brief Coarsening type */
    CoarseningType coarsening_type_;
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
     * @param[in] Coarsening type
     * @param[in] coarsening_level_ Coarsening level     
     */
    Lattice2D(const unsigned int Mt_lat_,
              const unsigned int Mx_lat_,
              const CoarseningType coarsening_type_,
              const int coarsening_level_=0) :
        Lattice(coarsening_level_),
        Mt_lat(Mt_lat_),
        Mx_lat(Mx_lat_),
        coarsening_type(coarsening_type_),
        rotated(false),
        last_coarsening(-1) {}

    /** @brief Return number of timeslices \f$M_{t,lat}\f$ */
    unsigned int getMt_lat() const {
        return Mt_lat;
    }

    /** @brief Return number of points in spatial direction \f$M_{x,lat}\f$ */
    unsigned int getMx_lat() const {
        return Mx_lat;
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
    
    /** @brief Convert cartesian lattice index of vertex to linear index
     *
     * Given a lattice site \f$n=(i,j)\f$, work out the corresponding linear index \f$\ell\f$.
     * The links are arranged such that
     * \f[
     *   \ell = M_{t,lat}j + i
     * \f]
     *
     *    ...   ...   ...   ...  
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  4-----5-----6-----7----(8)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  0-----1-----2-----3----(0)
     *
     * temporal direction (index i) = horizontal
     * spatial direction (index j) = vertical
     *
     *
     * @param[in] i Position index in temporal direction
     * @param[in] j Position index in spatial direction
     */
    inline unsigned int vertex_cart2lin(const int i, const int j) const {
        return Mt_lat*((j+Mx_lat)%Mx_lat) + ((i+Mt_lat)%Mt_lat);
    }

    /** @brief Convert linear index of vertex to lattice index
     *
     * Given \f$\ell = j M_{t,lat} + i\f$, work out cartesian index \f$(i,j)\f$.
     
     * @param[in] ell Linear index \f$\ell\f$
     * @param[out] i Position index in temporal direction
     * @param[out] j Position index in spatial direction
     */
    inline void vertex_lin2cart(const unsigned int ell, int& i, int& j) const {
        j = ell / Mt_lat;
        i = ell - Mt_lat*j;
    }
    
    /** @brief Convert Cartesian rotated lattice index to linear index.
     * 
     * This assumes that (i,j) describes a point on the rotated lattice, i.e.
     * i+j is a multiple of two.
     * 
     * @param i Position index in temporal direction
     * @param j Position index in spatial direction
     * 
     * @return linear index oCoarsenf vertex
     */
    inline unsigned int diag_vertex_cart2lin(const int i, const int j) const {
        assert((i+j)%2==0);        
        /* Rotated indices are only allowed if the lattice contains an
         * even number of cells in both directions
         */
        assert(Mx_lat%2==0);
        assert(Mt_lat%2==0);        
        int Mt_lat_half = Mt_lat/2;
        int Mx_lat_half = Mx_lat/2;
        int i_shifted = ((i+Mt_lat)-(i&1))/2;
        int j_shifted = ((j+Mx_lat)-(j&1))/2;
        // The first Mt_lat*Mx_lat/4 entries are reserved for points which
        // have both even i and j indices.
        int offset = (Mt_lat*Mx_lat/4)*(i&1);
        return Mt_lat_half*(j_shifted%Mx_lat_half)+i_shifted%Mt_lat_half+offset;
    }
    
    /** @brief Convert linear index of vertex to rotated lattice index
     *
     * Given the linear indexCoarsen \f$\ell\f$, work out cartesian index \f$(i,j)\f$.
     
     * @param[in] ell Linear index \f$\ell\f$
     * @param[out] i Position index in temporal direction
     * @param[out] j Position index in spatial direction
     */
    inline void diag_vertex_lin2cart(const unsigned int ell, int& i, int& j) const {
        /* Rotated indices are only allowed if the lattice contains an
         * even number of cells in both directions
         */
        int Mt_lat_half = Mt_lat/2;
        int Mx_lat_half = Mx_lat/2;
        int parity = ell / (Mt_lat*Mx_lat/4);
        unsigned int ell_half = ell - (Mt_lat*Mx_lat/4)*parity;
        int j_half = ell_half / Mt_lat_half;
        j = 2*j_half+parity;
        i = 2*(ell_half - Mt_lat_half*j_half)+parity;
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
     *      ...       ...       ...       ...      
     *  !         !         !         !         !
     *  9        11        13        15        (9)
     *  !         !         !         !         !
     *  o--- 8 ---o--- 10 --o--- 12 --o--- 14 --o
     *  !         !         !         !         !
     *  1         3         5         7        (1)
     *  !         !         !         !         !
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

    /** @brief Convert linear index of link to lattice index
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
     * @param[in] rho_coarsen_t coarsening factor in temporal direction
     * @param[in] rho_coarsen_x coarsening factor in spatial direction
     * @param[in] exit_on_failure Abort with an error message if construction is not possible?
     *
     * Returns lattice with twice the lattice spacing
     */
    virtual std::shared_ptr<Lattice2D> coarse_lattice(const bool exit_on_failure) {
        int rho_coarsen_t;   // temporal coarsening factor
        int rho_coarsen_x;   // spatial coarsening factor
        bool coarse_rotated; // Is the coarse lattice rotated?
        int this_coarsening; // Coarsening used to obtain this lattice (if coarsening in
                             // alternate directions
        switch (coarsening_type) {
            case CoarsenBoth:
                rho_coarsen_t = 2;
                rho_coarsen_x = 2;
            break;
            case CoarsenTemporal:
                rho_coarsen_t = 2;
                rho_coarsen_x = 1;
            break;
            case CoarsenSpatial:
                rho_coarsen_t = 1;
                rho_coarsen_x = 2;
            break;
            case CoarsenAlternate:
                if (not (last_coarsening==1)) {
                    rho_coarsen_t = 2;
                    rho_coarsen_x = 1;
                    this_coarsening = 1;
                } else {
                    rho_coarsen_t = 1;
                    rho_coarsen_x = 2;
                    this_coarsening = 0;
                }
            break;
            case CoarsenRotate:
            if (rotated) {
                rho_coarsen_t = 2;
                rho_coarsen_x = 2;
                coarse_rotated = false;
            } else {
                rho_coarsen_t = 1;
                rho_coarsen_x = 1;
                coarse_rotated = true;
            }
            break;
        }
        unsigned int Mt_lat_coarse = Mt_lat;
        unsigned int Mx_lat_coarse = Mx_lat;
        if ( rho_coarsen_t > 1 ) {
            if (Mt_lat%rho_coarsen_t) {
                if (exit_on_failure ) {
                    mpi_parallel::cerr << "ERROR: cannot coarsen 2d lattice with M_{t,lat} = " << Mt_lat;
                    mpi_parallel::cerr << " , M_{x,lat} = " << Mx_lat << " in temporal direction." << std::endl;
                    mpi_exit(EXIT_FAILURE);
                    throw std::runtime_error("...");
                } else {
                    return nullptr;
                }
            }
            Mt_lat_coarse = Mt_lat/rho_coarsen_t;
        }
        if ( rho_coarsen_x > 1 ) {
            if (Mx_lat%rho_coarsen_x) {
                if (exit_on_failure ) {
                    mpi_parallel::cerr << "ERROR: cannot coarsen 2d lattice with M_{t,lat} = " << Mt_lat;
                    mpi_parallel::cerr << " , M_{x,lat} = " << Mx_lat << " in spatial direction." << std::endl;
                    mpi_exit(EXIT_FAILURE);
                    throw std::runtime_error("...");
                } else {
                    return nullptr;
                }
            }
            Mx_lat_coarse = Mx_lat/rho_coarsen_x;
        }
        std::shared_ptr<Lattice2D> coarse_lattice;
        coarse_lattice = std::make_shared<Lattice2D>(Mt_lat_coarse,
                                                     Mx_lat_coarse,
                                                     coarsening_type,
                                                     coarsening_level+1);
        // Rotate coarse lattice if necessary
        if (coarsening_type == CoarsenRotate) {
            coarse_lattice->set_rotated(coarse_rotated);
        }
        // Set last coarsening of coarse lattice if necessary
        if (coarsening_type == CoarsenAlternate) {
            coarse_lattice->set_last_coarsening(this_coarsening);
        }
        return coarse_lattice;
    };
    
    /** @brief Is this lattice a rotated lattice? */
    const bool is_rotated() const { return rotated; }
     
    /** @brief return coarsening type */   
    const CoarseningType get_coarsening_type() const { return coarsening_type; }
    
protected:
    /** @brief Set rotation
     * 
     * @param[in] new_rotated new rotation direction
     */
     void set_rotated(const bool new_rotated) {
         if ( (new_rotated == true) and (not (coarsening_type == CoarsenRotate) ) ) {
            mpi_parallel::cerr << "ERROR: can only rotate lattice if coarsening type is \"rotate\"" << std::endl;
            mpi_exit(EXIT_FAILURE);
            throw std::runtime_error("...");
         }
         rotated = new_rotated;
     }
     
     /** @brief Set direction of coarsening used to obtain this lattice (if 
      *         alternate coarsening is used)
      * 
      * @param[in] last_coarsening_ coarsening direction (0=temporal, 1=spatial_
      */
      void set_last_coarsening(const int last_coarsening_) {
          last_coarsening = last_coarsening_;
      }

protected:
    /** @brief Number of time slices */
    const unsigned int Mt_lat;
    /** @brief Number of points in spatial direction */
    const unsigned int Mx_lat;
    /** @brief coarsening type of lattice */
    const CoarseningType coarsening_type;
    /** @brief is this lattice a rotated lattice? */
    bool rotated;
    /** @brief coarsening direction used to obtain this lattice */
    int last_coarsening;
};

#endif // LATTICE2D_HH
