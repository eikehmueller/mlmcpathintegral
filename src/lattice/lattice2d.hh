#ifndef LATTICE2D_HH
#define LATTICE2D_HH LATTICE2D_HH
#include <memory>
#include <stdexcept>
#include <cassert>
#include <vector>
#include <array>
#include <algorithm>
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
 * temporal and spatial direction, as well as the type of coarsening to be used.\f$
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
 * 
 *
 */
class Lattice2D : public Lattice {
public:
    /** @brief Initialise class
     *
     * The lattice can either be unrotated or rotated, as shown in the following
     * figures:
     * 
     *   +---+---+---+---+      +       +       +
     *   !   !   !   !   !        \   /   \   /
     *   +---+---+---+---+          +       +
     *   !   !   !   !   !        /   \   /   \ 
     *   +---+---+---+---!      +       +       +
     *   !   !   !   !   !        \   /   \   /
     *   +---+---+---+---+          +       +
     *   !   !   !   !   !        /   \   /   \
     *   +---+---+---+---+      +       +       +
     *     unrotated                 rotated
     * 
     * Rotated lattices are represented by only using some of the vertices in an
     * unrotated lattice. More specifically, only vertices (i,j) with even i+j are
     * taken into account.
     *
     * @param[in] Mt_lat_ Number of time slices
     * @param[in] Mx_lat_ Number of points in spatial direction
     * @param[in] coarsening_type_ (see enum CoarseningType )
     * @param[in] coarsening_level_ Coarsening level (0=finest lattice)
     */
    Lattice2D(const unsigned int Mt_lat_,
              const unsigned int Mx_lat_,
              const CoarseningType coarsening_type_,
              const int coarsening_level_=0);
    
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
        if (rotated) {
            return Mt_lat*Mx_lat;        
        } else {
            return 2*Mt_lat*Mx_lat;
        }
    }

    /** @brief Return number of vertices (assuming periodic boundaries) */
    virtual const unsigned int getNvertices() const {
        if (rotated) {
            return Mt_lat*Mx_lat/2;
        } else {
            return Mt_lat*Mx_lat;
        }
    }

    /** @brief Return number of cells (assuming periodic boundaries) */
    unsigned int getNcells() const {
        if (rotated) {
            return Mt_lat*Mx_lat/2;
        } else {
            return Mt_lat*Mx_lat;            
        }
    }
    
    /** @brief Convert cartesian lattice index of vertex to linear index
     *
     * Given a lattice site \f$n=(i,j)\f$, work out the corresponding linear index \f$\ell\f$.
     * 
     * For an unrotated lattice the vertices are arranged such that
     * 
     * \f[
     *   \ell = M_{t,lat}j + i
     * \f]
     * 
     *  This is illustrated in the following figure:
     *
     * (0)---(1)---(2)---(3)---(0)
     *  !     !     !     !     !
     *  !     !     !     !     !
     * 12----13----14----15---(12)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  8-----9----10----11----(8)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  4-----5-----6-----7----(4)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  0-----1-----2-----3----(0)
     *
     * temporal direction (index i) = horizontal
     * spatial direction (index j) = vertical
     * 
     * For a rotated lattice the arrangement is this:
     * 
     * (0)----+----(1)----+----(0)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  +-----6-----+-----7-----+
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  2-----+-----3-----+----(8)
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  +-----4-----+-----5-----+
     *  !     !     !     !     !
     *  !     !     !     !     !
     *  0-----+-----1-----+----(0)
     * 
     * In other words, first the indices are distributed over the lattice where both
     * i and j are even, followed by the vertices where both i and j are odd.
     * 
     * The lattice can be coarsened in different ways:
     * 
     *   1. in both directions simultaneouly
     *   2.(a) in the temporal direction only
     *   2.(b) in the spatial direction only
     *   3. in the temporal and spatial direction in subsequent coarsenings
     *   4. by rotating the lattice by 45 degrees
     * 
     * The lattice also stores the following information:
     *   1. a list of all vertices that are present on the next-coarser lattice
     *   2. a list of all vertices that are NOT present on the next-coarser lattice
     *   3. a list of all neighbouring vertices for a given vertex
     *   4. a map from the coarse-only indices in 1. to the vertices on the 
     *      next-coarser lattice
     *
     * @param[in] i Position index in temporal direction
     * @param[in] j Position index in spatial direction
     */
    inline unsigned int vertex_cart2lin(const int i, const int j) const {
        if (rotated) {
            assert((i+j)%2==0);        
            int Mt_lat_half = Mt_lat/2;
            int Mx_lat_half = Mx_lat/2;
            int i_shifted = ((i+Mt_lat)-(i&1))/2;
            int j_shifted = ((j+Mx_lat)-(j&1))/2;
            // The first Mt_lat*Mx_lat/4 entries are reserved for points which
            // have both even i and j indices.
            int offset = (Mt_lat*Mx_lat/4)*(i&1);
            return Mt_lat_half*(j_shifted%Mx_lat_half)+i_shifted%Mt_lat_half+offset;
        } else {
            return Mt_lat*((j+Mx_lat)%Mx_lat) + ((i+Mt_lat)%Mt_lat);
        }
    }

    /** @brief Convert linear index of vertex to lattice index
     *
     * Given \f$\ell = j M_{t,lat} + i\f$, work out cartesian index \f$(i,j)\f$.
     
     * @param[in] ell Linear index \f$\ell\f$
     * @param[out] i Position index in temporal direction
     * @param[out] j Position index in spatial direction
     */
    inline void vertex_lin2cart(const unsigned int ell, int& i, int& j) const {
        if (rotated) {
            int Mt_lat_half = Mt_lat/2;
            int Mx_lat_half = Mx_lat/2;
            int parity = ell / (Mt_lat*Mx_lat/4);
            unsigned int ell_half = ell - (Mt_lat*Mx_lat/4)*parity;
            int j_half = ell_half / Mt_lat_half;
            j = 2*j_half+parity;
            i = 2*(ell_half - Mt_lat_half*j_half)+parity;            
        } else {
            j = ell / Mt_lat;
            i = ell - Mt_lat*j;            
        }
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
        // Currently links can only be handled on non-rotated lattices
        assert (not rotated);
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
        // Currently links can only be handled on non-rotated lattices
        assert (not rotated);
        j = ell / (2*Mt_lat); // j = ell / (2*M_{t,lat})
        unsigned int r = ell - (2*Mt_lat)*j;
        i = r >> 1; // i = r/2
        mu = r & 1; // mu = r % 2
    }

    /** @brief Return coarse lattice
     *
     * Returns coarsened version of lattice
     */
    virtual std::shared_ptr<Lattice2D> get_coarse_lattice() {
        return coarse_lattice;
    }
    
    /** @brief Is this lattice a rotated lattice? */
    const bool is_rotated() const { return rotated; }
     
    /** @brief return coarsening type */   
    const CoarseningType get_coarsening_type() const { return coarsening_type; }
    
    /** @brief return list of coarse vertices
     * 
     * This list contains the linear indices of all vertices that correspond to a
     * vertex on the next-coarser lattice.
     */
    const std::vector<unsigned int>& get_coarse_vertices() {
        return coarse_vertices;
    }
    
    /** @brief return list of fine-only dofs
     * 
     * This list contains the linear indices of all vertices that DO NOT have
     * correspond to a vertex on the next-coarser lattice.

     */
    const std::vector<unsigned int>& get_fineonly_vertices() {
        return fineonly_vertices;
    }
    
    /** @brief return fine-to-coarse map 
     * 
     * For each vertex on the current lattice that corresponds to a coarse-lattice
     * site
     * */
    const std::map<unsigned int, unsigned int>& get_fine2coarse_map() {
        return fine2coarse_map;
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
    /** @brief coarse lattice */
    std::shared_ptr<Lattice2D> coarse_lattice;
    /** @brief list of vertices that also exist on next-coarser lattice */
    std::vector<unsigned int> coarse_vertices;
    /** @brief list of vertices that DO NOT exist on next-coarser lattice */
    std::vector<unsigned int> fineonly_vertices;
    /** @brief map of coarse vertices to the corrresponding vertex on coarse lattice */
    std::map<unsigned int, unsigned int> fine2coarse_map;
};

#endif // LATTICE2D_HH
