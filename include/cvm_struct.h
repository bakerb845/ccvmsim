#include <limits.h>
#ifndef _cvm_struct_h__
#define _cvm_struct_h__ 1

enum mesh_element_type
{
    HEX8 = 1,        /*!< Linear hexahedron element type */
    HEX27 = 2,       /*!< Quadratic hexahedron element type */ 
    TET4 = 3,        /*!< Linear tetrahedron */
    TET8 = 4         /*!< Quadratic tetrahedron */
};

enum coarsen_type
{
    COARSEN_NONE = 1,   /*!< No coarsening */
    COARSEN_DOUBLE = 2, /*!< Doubling layer */
    COARSEN_TRIPLE = 3  /*!< Tripling layer */
};

enum mesh_bc_enum
{
    NO_BC = 0,        /*!< This is not a boundary condition */

    SOUTH_BDRY  = 1,  /*!< This is a south boundary */ 
    EAST_BDRY   = 2,  /*!< This is an east boundary */
    NORTH_BDRY  = 3,  /*!< This is a north boundary */
    WEST_BDRY   = 4,  /*!< This is a west boundary */
    BOTTOM_BDRY = 5,  /*!< This is a bottom boundary */
    TOP_BDRY    = 6,  /*!< This is a top boundary */ 

    SOUTH_FREE  = 21, /*!< Zero traction Neumann condition on face 1 */
    EAST_FREE   = 22, /*!< Zero traction Neumann condition on face 2 */
    NORTH_FREE  = 23, /*!< Zero traction Neumann condition on face 3 */
    WEST_FREE   = 24, /*!< Zero traction Neumann condition on face 4 */
    BOTTOM_FREE = 25, /*!< Zero traction Neumann condition on face 5 */
    TOP_FREE    = 26, /*!< Zero traction Neumann condition on face 6 */

    SOUTH_CE    = 31, /*!< Non-reflecting Neumann condition on face 1 */
    EAST_CE     = 32, /*!< Non-reflecting Neumann condition on face 2 */
    NORTH_CE    = 33, /*!< Non-reflecting Neumann condition on face 3 */
    WEST_CE     = 34, /*!< Non-reflecting Neumann condition on face 4 */
    BOTTOM_CE   = 35, /*!< Non-reflecting Neumann condition on face 5 */
    TOP_CE      = 36, /*!< Non-reflecting Neumann condition on face 6 */

    SOUTH_INT   = 41, /*!< Internal boundary on face 1 */
    EAST_INT    = 42, /*!< Internal boundary on face 2 */
    NORTH_INT   = 43, /*!< Internal boundary on face 3 */
    WEST_INT    = 44, /*!< Internal boundary on face 4 */
    BOTTOM_INT  = 45, /*!< Internal boundary on face 5 */
    TOP_INT     = 46, /*!< Internal boundary on face 6 */

};

struct mesh_element_struct
{
    double *vp;                  /*!< Compressional velocity of ia'th anchor
                                      node (m/s) */ 
    double *vs;                  /*!< Shear velocity of ia'th anchor node
                                       (m/s) */
    double *dens;                /*!< Density of ia'th anchor node (kg/m**3) */
    double *Qp;                  /*!< Compressional velocity quality factor at
                                      ia'th node */
    double *Qs;                  /*!< Shear velocity quality factor at ia'th
                                      node */
    double *x;                   /*!< x position of ia'th anchor node */
    double *y;                   /*!< y position of ia'th anchor node */
    double *z;                   /*!< z position of ia'th anchor node */
    int *ien_face;               /*!< Maps from ia'th node on the iface'th
                                      element face node to global anchor
                                      node (leading dimension 4 for hexes
                                      and 3 for tets) */
    int *ien;                    /*!< Maps from ia'th element node to global
                                      anchor node */
    int *bc;                     /*!< Boundary condition for iface'th face */
    int *neighbor;               /*!< Maps from iface'th element face to 
                                      neighbor element number */
    int nface;                   /*!< Number of faces on element */
    int ngnod;                   /*!< Number of anchor nodes on element */
    int ngnod_face;              /*!< Number of anchor nodes per element face */
    enum mesh_element_type type; /*!< Element type (probably HEX8) */
};

struct mesh_struct
{
    struct mesh_element_struct *element; /*!< Element structure */
    double *xlocs;        /*!< X anchor node locations (m) [nnpg] */
    double *ylocs;        /*!< Y anchor node locations (m) [nnpg] */
    double *zlocs;        /*!< Z anchor node locations (m) [nnpg] */
    double *vp;           /*!< Compressional velocity at anchor 
                               nodes (m/s) [nnpg] */
    double *vs;           /*!< Shear velocity at anchor nodes (m/s) [nnpg] */
    double *dens;         /*!< Density at anchor nodes (kg/m**3) [nnpg] */
    double *Qp;           /*!< Compressional velocity quality factor at
                               ia'th node [nnpg] */
    double *Qs;           /*!< Shear velocity quality factor at ia'th
                               node [nnpg] */
    int nnpg;             /*!< Number of anchor nodes in mesh */
    int nelem;            /*!< Number of elements in mesh */
    bool lptr_only;       /*!< If true element contains integer pointers only
                               If false element contains integer pointers
                               and material information */
    bool lhomog;          /*!< True then the mesh is homogeneous */
    char pad[6];
};


struct cvm_parms_struct
{
    char cvm_moddir[PATH_MAX];   /*!< CVM directory */
    char cvm_outputdir[PATH_MAX];  /*!< CVM output directory */
    char mesh_outputdir[PATH_MAX]; /*!< Mesh output directory */
    char nll_outputdir[PATH_MAX];  /*!< NLL output directory */
    char topofl[PATH_MAX];       /*!< topo30 topography file */
    double *dxl_cvm;             /*!< CVM grid spacing (m) in x [nlay_cvm] */
    double *dyl_cvm;             /*!< CVM grid spacing (m) in y [nlay_cvm] */
    double *dzl_cvm;             /*!< CVM grid spacing (m) in z [nlay_cvm] */
    double *zcoarsen;            /*!< Interfaces at which to coarsen mesh
                                      (m). [ncoarsen]  */
    double *z0_cvm;              /*!< CVM layer bottom interfaces (m) 
                                      [nlay_cvm] */
    int *nxl_cvm;                /*!< CVM x grid points in each layer
                                      [nlay_cvm] */
    int *nyl_cvm;                /*!< CVM y grid points in each layer
                                      [nlay_cvm] */
    int *nzl_cvm;                /*!< CVM z grid points in each layer
                                      [nlay_cvm] */
    enum coarsen_type *coarsen;  /*!< Refinement type for layer i
                                            to layer i + 1 [ncoarsen]*/
    double dx_fem;               /*!< Mesh grid spacing in x in top layer (m) */
    double dy_fem;               /*!< Mesh grid spacing in y in top layer (m) */
    double dz_fem;               /*!< Mesh grid spacing in z in top layer (m) */
    double common_mult;          /*!< Common multiple grid spacing (m) for
                                      all layers */
    double utm_x0_cvm;           /*!< CVM lower left corner (m) - zone 10 */
    double utm_y0_cvm;           /*!< CVM lower left corner (m) - zone 10 */
    double utm_x1_cvm;           /*!< CVM upper right corner (m) - zone 10 */
    double utm_y1_cvm;           /*!< CVM upper right corner (m) - zone 10 */
    double utm_x0;               /*!< Simulation lower left corner (m) */
    double utm_y0;               /*!< Simulation lower left corner (m) */ 
    double utm_x1;               /*!< Simulation upper right corner (m) */
    double utm_y1;               /*!< Simulation upper right corner (m) */
    double lat0_cvm;             /*!< CVM lower left corner (degrees) */
    double lon0_cvm;             /*!< CVM lower left corner (degrees) */
    double lat1_cvm;             /*!< CVM upper right corner (degrees) */
    double lon1_cvm;             /*!< CVM upper right corner (degrees) */
    double lat0;                 /*!< Lower left bounding box for simulation
                                      region (degrees) */
    double lon0;                 /*!< Lower left bounding box for simulation
                                      region (degrees) - positive east */
    double lat1;                 /*!< Upper right bounding box for simulation 
                                      region (degrees) */
    double lon1;                 /*!< Upper right bounding box for simulation 
                                      region (degrees) */
    double zmin;                 /*!< Minimum model depth (m) */
    double zmax;                 /*!< Maximum model depth (m) */
    double vp_min;               /*!< Minimum allowable P velocity (m/s) */
    double vp_max;               /*!< Maximum allowable P velocity (m/s) */
    double vs_min;               /*!< Minimum allowable S velocity (m/s) */
    double vs_max;               /*!< Maximum allowable S velocity (m/s) */
    double dens_min;             /*!< Minimum allowable density (m/s) */
    double dens_max;             /*!< Maximum allowable density (m/s) */
    double vpvs_min;             /*!< Minimum Vp/Vs ratio allowed */
    double vpvs_max;             /*!< Maximum Vp/Vs ratio allowed */
    double vpvs_ratio;           /*!< Fixed Vp/Vs ratio */
    double ztopo_min;            /*!< Depth (km) below flat free surface at which
                                      which the mesh deformation will no longer
                                      be linearly interpolated to mimic
                                      topography */
    int ncoarsen;                /*!< Number of coarsening interfaces */
    int utmzone_cvm;             /*!< CVM utm zone - should be 10 */
    int nlay_cvm;                /*!< Number of layers in CVM */
    int utm_zone;                /*!< UTM zone for unpacking topo30 - 
                                      should be 10 */
    bool ltopo;                  /*!< If true then set topography from topo30 */
    bool lthresh_vp;             /*!< If true then threshold vp */
    bool lthresh_vs;             /*!< If true then threshold vs */
    bool lthresh_vpvs;           /*!< If true then threshold vp/vs ratios */
    bool lthresh_dens;           /*!< If true then threshold density */
    bool setvp_from_vs;          /*!< If true then compute Vp from Vs using
                                      vpvs_ratio */
    bool setvs_from_vp;          /*!< If true then compute Vs from Vp using
                                      vpvs_ratoi */ 
    char pad[1];
};

struct cvm_model_struct
{
    double *vp;    /*!< Compressional (m/s) velocity at nodes in CVM.
                        The packing is fastest in x increasing east,
                        intermediate in y increasing north, and
                        slowest in z increasing down from the free
                        surface [npts] */
    double *vs;    /*!< Shear velocity (m/s) at nodes in CVM [npts] */
    double *dens;  /*!< Density (kg/m**3) at nodes in CVM [npts] */
    double *Qp;    /*!< Compressional velocity quality factor [npts] */
    double *Qs;    /*!< Shear velocity quality factor [npts] */
    double *xlocs; /*!< x node locations in CVM [npts] */
    double *ylocs; /*!< y node locations in CVM [npts] */
    double *zlocs; /*!< z node locations in CVM [npts] */
    double x0;     /*!< x origin (upper southwest corner) */
    double y0;     /*!< y origin (upper southwest corner) */
    double z0;     /*!< z origin (upper southwest corner) */
    double z1;     /*!< z max depth */ 
    double dx;     /*!< Grid spacing in x (m) */
    double dy;     /*!< Grid spacing in y (m) */
    double dz;     /*!< Grid spacing in z (m) */
    int nx;        /*!< Number of x grid points */
    int ny;        /*!< Number of y grid points */
    int nz;        /*!< Number of z grid points */
    int npts;      /*!< Number of points in CVM */
};
#endif
