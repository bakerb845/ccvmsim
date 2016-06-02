#include <stdio.h>
#include <stdlib.h>
#include <iscl/log/log.h>
#include "cvm.h"

/*!
 * @brief Determines number of elements in a regular mesh
 *
 * @param[in] nx       number of grid points in x (> 1)
 * @param[in] ny       number of grid points in y (> 1)
 * @param[in] nz       number of grid points in z (> 1)
 *
 * @param[out] nelem   number of elements in hexahedral mesh
 *
 * @result number of elements in regular hex mesh
 *
 * @author Ben Baker, ISTI
 *
 */
int regmesh_getNumberOfElements(int nx, int ny, int nz, int *nelem)
{
    const char *fcnm = "regmesh_getNumberOfElements\0";
    *nelem = 0;
    if (nx < 2){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, nx);
        return -1;
    }
    if (ny < 2){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, ny);
        return -1;
    }
    if (nz < 2){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, nz);
        return -1;
    }
    *nelem = (nx - 1)*(ny - 1)*(nz - 1);
    return 0;
}
//============================================================================//
/*!
 * @brief Makes the connectivity for a regular hex mesh
 *
 * @param[in] nx        number of x grid points (> 1)
 * @param[in] ny        number of y grid points (> 1)
 * @param[in] nz        number of z grid points (> 1)
 * 
 * @param[out] ien      maps from ia'th anchor node on ielem'th element to
 *                      global anchor node [(nx-1)*(ny-1)*(nz-1)*8]
 *
 * @author Ben Baker, ISTI
 *
 */
int regmesh_makeHexIEN(int nx, int ny, int nz, int *ien)
{ 
    const char *fcnm = "regmesh_makeHexIEN\0";
    int ielem, indx, ix, iy, iz, nnpg_x, nnpg_xy;
    const int ngnod = 8;
    if (nx < 2 || ny < 2 || nz < 2){
        printf("%s: Invalid number of points (nx,ny,nz)=(%d,%d,%d)\n",
               fcnm, nx, ny, nz);
        return -1;
    }
    // Define some constants for regmesh 
    nnpg_x  = nx; 
    nnpg_xy = nx*ny;
    // Loop on mesh 
    ielem = 0;
    for (iz=0; iz<nz-1; iz++){
        for (iy=0; iy<ny-1; iy++){
            for (ix=0; ix<nx-1; ix++){
                // Global anchor nodes
                indx = ngnod*ielem;
                ien[indx+0] = nnpg_xy*iz       + nnpg_x*iy       + ix; 
                ien[indx+1] = nnpg_xy*iz       + nnpg_x*iy       + ix + 1;
                ien[indx+2] = nnpg_xy*iz       + nnpg_x*(iy + 1) + ix + 1;
                ien[indx+3] = nnpg_xy*iz       + nnpg_x*(iy + 1) + ix; 
                ien[indx+4] = nnpg_xy*(iz + 1) + nnpg_x*iy       + ix; 
                ien[indx+5] = nnpg_xy*(iz + 1) + nnpg_x*iy       + ix + 1;
                ien[indx+6] = nnpg_xy*(iz + 1) + nnpg_x*(iy + 1) + ix + 1;
                ien[indx+7] = nnpg_xy*(iz + 1) + nnpg_x*(iy + 1) + ix;
                ielem = ielem + 1;
            } // Loop on x elements
        } // Loop on y elements
    } // Loop on z elements
    return 0;
}
//============================================================================//
/*!
 * @brief Generates pointers for a regular hexahedral mesh
 *
 * @reference http://prod.sandia.gov/techlib/access-control.cgi/1992/922137.pdf
 *
 * @param[in] nx         number of grid points in x (> 1)
 * @param[in] ny         number of grid points in y (> 1)
 * @param[in] nz         number of grid points in z (> 1)
 * @param[in] nelem      number of elements in mesh ( (nx-1)*(ny-1)*(nz-1) )
 *
 * @param[out] element   on successful exit determines the element connectivity,
 *                       boundary conditions, and neighbors for all elements 
 *                       [nelem]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int regmesh_makeHexMeshPointers(int nx, int ny, int nz,
                                int nelem, struct mesh_element_struct *element)
{
    const char *fcnm = "regmesh_makeHexMeshPointers\0";
    int ielem, ierr, ix, iy, iz,
        nelem_ref, nelem_x, nelem_xy, nnpg_x, nnpg_xy;
    ierr = regmesh_getNumberOfElements(nx, ny, nz, &nelem_ref); 
    if (ierr != 0){
        log_errorF("%s: Error in size estimation\n", fcnm);
        return -1;
    }
    if (nelem_ref != nelem){
        if (nelem_ref > nelem){
            log_errorF("%s: Invalid space in element\n", fcnm);
            return -1;
        }else{
            log_warnF("%s: Inconsistent element sizing\n", fcnm);
        }
    }
    // Define some constants for regmesh 
    nnpg_x  = nx;
    nnpg_xy = nx*ny;
    nelem_x = (nx - 1);
    nelem_xy = (nx - 1)*(nz - 1);
    // Loop on mesh 
    ielem = 0;
    for (iz=0; iz<nz-1; iz++){
        for (iy=0; iy<ny-1; iy++){
            for (ix=0; ix<nx-1; ix++){
                // Element type 
                element[ielem].type = HEX8;
                // Global anchor nodes
                element[ielem].ien[0] = nnpg_xy*iz       + nnpg_x*iy       + ix;
                element[ielem].ien[1] = nnpg_xy*iz       + nnpg_x*iy       + ix + 1;
                element[ielem].ien[2] = nnpg_xy*iz       + nnpg_x*(iy + 1) + ix + 1;
                element[ielem].ien[3] = nnpg_xy*iz       + nnpg_x*(iy + 1) + ix;
                element[ielem].ien[4] = nnpg_xy*(iz + 1) + nnpg_x*iy       + ix;
                element[ielem].ien[5] = nnpg_xy*(iz + 1) + nnpg_x*iy       + ix + 1;
                element[ielem].ien[6] = nnpg_xy*(iz + 1) + nnpg_x*(iy + 1) + ix + 1;
                element[ielem].ien[7] = nnpg_xy*(iz + 1) + nnpg_x*(iy + 1) + ix;
                // Faces (table 2); face 1
                element[ielem].ien_face[0]  = element[ielem].ien[0];
                element[ielem].ien_face[1]  = element[ielem].ien[1];
                element[ielem].ien_face[2]  = element[ielem].ien[5];
                element[ielem].ien_face[3]  = element[ielem].ien[4];
                element[ielem].neighbor[0] = ielem - nelem_x;
                element[ielem].bc[0] = NO_BC;
                if (iy == 0){element[ielem].neighbor[0] = 0;}
                if (iy == 0){element[ielem].bc[0] = SOUTH_BDRY;}
                // face 2
                element[ielem].ien_face[4]  = element[ielem].ien[1];
                element[ielem].ien_face[5]  = element[ielem].ien[2];
                element[ielem].ien_face[6]  = element[ielem].ien[6];
                element[ielem].ien_face[7]  = element[ielem].ien[5];
                element[ielem].neighbor[1] = ielem + 1;
                element[ielem].bc[1] = NO_BC;
                if (ix == nx - 2){element[ielem].neighbor[1] = 0;}
                if (ix == nx - 2){element[ielem].bc[1] = EAST_BDRY;}
                // face 3
                element[ielem].ien_face[8]  = element[ielem].ien[2];
                element[ielem].ien_face[9]  = element[ielem].ien[3];
                element[ielem].ien_face[10] = element[ielem].ien[7];
                element[ielem].ien_face[11] = element[ielem].ien[6];
                element[ielem].neighbor[2] = ielem + nelem_x;
                element[ielem].bc[2] = NO_BC;
                if (iy == ny - 2){element[ielem].neighbor[2] = 0;}
                if (iy == ny - 2){element[ielem].bc[2] = NORTH_BDRY;}
                // face 4
                element[ielem].ien_face[12] = element[ielem].ien[0];
                element[ielem].ien_face[13] = element[ielem].ien[4];
                element[ielem].ien_face[14] = element[ielem].ien[7];
                element[ielem].ien_face[15] = element[ielem].ien[3];
                element[ielem].neighbor[3] = ielem - 1;
                element[ielem].bc[3] = NO_BC;
                if (ix == 0){element[ielem].neighbor[3] = 0;}
                if (ix == 0){element[ielem].bc[3] = WEST_BDRY;}
                // face 5
                element[ielem].ien_face[16] = element[ielem].ien[0];
                element[ielem].ien_face[17] = element[ielem].ien[3];
                element[ielem].ien_face[18] = element[ielem].ien[2];
                element[ielem].ien_face[19] = element[ielem].ien[1];
                element[ielem].neighbor[4] = ielem - nelem_xy;
                element[ielem].bc[4] = NO_BC;
                if (iz == 0){element[ielem].neighbor[4] = 0;}
                if (iz == 0){element[ielem].bc[4] = BOTTOM_BDRY;}
                // face 6
                element[ielem].ien_face[20] = element[ielem].ien[4];
                element[ielem].ien_face[21] = element[ielem].ien[5];
                element[ielem].ien_face[22] = element[ielem].ien[6];
                element[ielem].ien_face[23] = element[ielem].ien[7];
                element[ielem].neighbor[5] = ielem + nelem_xy;
                element[ielem].bc[5] = NO_BC; 
                if (iz == nz - 2){element[ielem].neighbor[5] = 0;}
                if (iz == nz - 2){element[ielem].bc[5] = TOP_BDRY;}
                // update element counter
                ielem = ielem + 1;
            } // Loop on x
        } // Loop on y
    } // Loop on z
    // Fidelity check
    if (ielem != nelem){
        log_errorF("%s: Failed to initialize elements\n", fcnm);
        ierr = 1;
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief For debugging purposes this makes a constant material Earth model
 *
 * @param[in] vp           constant p-velocity (km/s)
 * @param[in] vs           constant s-velocity (km/s)
 * @param[in] dens         constant density (kg/m**3)
 * @param[in] Qp           P wave quality factor
 * @param[in] Qs           S wave quality factor
 *
 * @param[inout] element   on input contains the number of anchor nodes
 *                         on each element
 *                         on output now includes the compressional
 *                         and shear velocities and the density at 
 *                         each element's on the anchor nodes
 *
 * @author Ben Baker, ISTI
 *
 */
void regmesh_constantMaterialModel(double vp, double vs, double dens,
                                   double Qp, double Qs,
                                   int nelem,
                                   struct mesh_element_struct *element)
{
    int ia, ielem, ngnod;
    for (ielem=0; ielem<nelem; ielem++){
        ngnod = element[ielem].ngnod;
        element[ielem].vp   = (double *)calloc(ngnod, sizeof(double));
        element[ielem].vs   = (double *)calloc(ngnod, sizeof(double));
        element[ielem].dens = (double *)calloc(ngnod, sizeof(double));
        element[ielem].Qp = (double *)calloc(ngnod, sizeof(double));
        element[ielem].Qs = (double *)calloc(ngnod, sizeof(double));
        for (ia=0; ia<ngnod; ia++){
            element[ielem].vp[ia] = vp;
            element[ielem].vs[ia] = vs;
            element[ielem].dens[ia] = dens;
            element[ielem].Qp[ia] = Qp;
            element[ielem].Qs[ia] = Qs;
        }
    }
    return;
}
//============================================================================//
/*!
 * @brief Copies a regular model ordered fastest in x, intermediate in y,
 *        and slowest in z to the element structure
 *
 * @param[in] lflip     If True then flip the z axis
 * @param[in] nx        number of x grid points
 * @param[in] ny        number of y grid points
 * @param[in] nz        number of z grid points
 * @param[in] nelem     number of elements in mesh
 * @param[in] vp        compressional velocity to copy to element struct
 *                      [nx*ny*nz]
 * @param[in] vs        shear velocity to copy to element structure
 *                      [nx*ny*nz] 
 * @param[in] dens      density to copy to element structure [nx*ny*nz]
 *
 * @param[out] vpw      vp velocities with appropriate orientation [nx*ny*nz] 
 * @param[out] vsw      vs velocities with appropriate orientation [nx*ny*nz]
 * @param[out] densw    densities with appropriate orientation [nx*ny*nz]
 * @param[out] Qpw      Qp with appropriate orientation [nx*ny*nz]
 * @param[out] Qsw      Qs with appropriate orientation [nx*ny*nz]
 *
 * @author Ben Baker, ISTI
 * 
 */
void __regmesh_copyRegularModel(bool lflip, int nx, int ny, int nz,
                                const double *__restrict__ vp,
                                const double *__restrict__ vs,
                                const double *__restrict__ dens,
                                const double *__restrict__ Qp,
                                const double *__restrict__ Qs,
                                double *__restrict__ vpw,
                                double *__restrict__ vsw,
                                double *__restrict__ densw,
                                double *__restrict__ Qpw,
                                double *__restrict__ Qsw)
{
    int indx, inpg, ixy, iz, nnpg, nxy;
    if (lflip){
        nxy = nx*ny;
        #pragma omp parallel for \
         firstprivate (nxy), \
         private(indx, inpg), \
         shared(dens, densw, Qp, Qpw, Qs, Qsw, vp, vpw, vs, vsw)
        for (iz=0; iz<nz; iz++){
            indx = (nz-1-iz)*nxy;
            inpg = iz*nxy;
            #pragma omp simd
            for (ixy=0; ixy<nxy; ixy++){
                vpw[inpg] =   vp[indx];
                vsw[inpg] =   vs[indx];
                densw[inpg] = dens[indx];
                Qpw[inpg] = Qp[indx];
                Qsw[inpg] = Qs[indx];
                indx = indx + 1;
                inpg = inpg + 1;
            }
        }
    }else{
        nnpg = nx*ny*nz; 
        #pragma omp simd
        for (inpg = 0; inpg<nnpg; inpg++){
            vpw[inpg] =   vp[inpg];
            vsw[inpg] =   vs[inpg];
            densw[inpg] = dens[inpg];
            Qpw[inpg] = Qp[inpg];
            Qsw[inpg] = Qs[inpg];
        }
    }
    return;
}
//============================================================================//
/*!
 * @brief Copies a regular model ordered fastest in x, intermediate in y,
 *        and slowest in z to the element structure
 *
 * @param[in] lflip       If True then flip the z axis
 * @param[in] nx          number of x grid points
 * @param[in] ny          number of y grid points
 * @param[in] nz          number of z grid points
 * @param[in] nelem       number of elements in mesh
 * @param[in] nnpg        number of anchor nodes in mesh (nx*ny*nz)
 * @param[in] vp          compressional velocity to copy to element struct [nnpg]
 * @param[in] vs          shear velocity to copy to element structure [nnpg] 
 * @param[in] dens        density to copy to element structure [nnpg]
 * @param[in] Qp          P velocity quality factor to copy to element struct
 *                        [nnpg]
 * @param[in] Qs          S velocity quality factor to copy to element struct
 *                        [nnpg] 
 *
 * @param[inout] element  on input holds the IEN array for each element
 *                        and space for the material arrays
 *                        on output holds the materials on each element
 *
 * @author Ben Baker, ISTI
 */
int regmesh_copyRegularModel(bool lflip, int nx, int ny, int nz,
                             int nelem, int nnpg,
                             double *__restrict__ vp,
                             double *__restrict__ vs,
                             double *__restrict__ dens,
                             double *__restrict__ Qp,
                             double *__restrict__ Qs,
                             struct mesh_element_struct *element)
{
    const char *fcnm = "regmesh_copyMaterial\0";
    double *vpw, *vsw, *Qpw, *Qsw, *densw;
    int ierr;
    if (nnpg != nx*ny*nz){
        printf("%s: Error size inconsistent\n", fcnm);
        return -1;
    }
    if (nelem < 1){
        printf("%s: Error no elements\n", fcnm);
        return -1;
    }
    if (nnpg < 1){
        printf("%s: Error no anchor ndoes\n", fcnm);
        return -1;
    }
    // Pack the flipped model 
    if (lflip){
        vpw  = (double *)calloc(nnpg, sizeof(double));
        vsw  = (double *)calloc(nnpg, sizeof(double));
        densw = (double *)calloc(nnpg, sizeof(double));
        Qpw = (double *)calloc(nnpg, sizeof(double));
        Qsw = (double *)calloc(nnpg, sizeof(double));
        __regmesh_copyRegularModel(lflip, nx, ny, nz, 
                                   vp, vs, dens, Qp, Qs,
                                   vpw, vsw, densw, Qpw, Qsw);
        // Fill the model
        ierr = mesh_element__setAnchorNodeProperties(nelem, nnpg,
                                                     vpw, vsw, densw, Qpw, Qsw,
                                                     element);
        // Free space
        free(vpw);
        free(vsw);
        free(densw); 
        free(Qpw);
        free(Qsw);
        vpw   = NULL;
        vsw   = NULL;
        densw = NULL;
        Qpw = NULL;
        Qsw = NULL;
    }else{
        // Fill the model
        ierr = mesh_element__setAnchorNodeProperties(nelem, nnpg,
                                                     vp, vs, dens, Qp, Qs,
                                                     element);
    }
    if (ierr != 0){
        printf("%s: Error setting anchor node properties!\n", fcnm);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Determines number of anchor nodes in a regular mesh
 *
 * @param[in] nx       number of grid points in x (> 0)
 * @param[in] ny       number of grid points in y (> 0)
 * @param[in] nz       number of grid points in z (> 0)
 *
 * @param[out] nelem   number of elements in hexahedral mesh
 *
 * @result number of elements in regular hex mesh
 *
 * @author Ben Baker, ISTI
 *
 */
int regmesh_getNumberOfAnchorNodes(int nx, int ny, int nz, int *nnpg)
{
    const char *fcnm = "regmesh_getNumberOfAnchorNodes\0";
    if (nx < 1){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, nx);
        return -1;
    }
    if (ny < 1){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, ny);
        return -1;
    }
    if (nz < 1){
        log_errorF("%s: Invalid number of x grid points %d\n", fcnm, nz);
        return -1;
    }
    *nnpg = nx*ny*nz;
    return 0;
}
//============================================================================//
/*!
 * @brief Function for computing the regular spaced nodes packed
 *        with x fastest, y intermediate, and z slowest
 *
 * @param[in] nx          number of grid points in x (> 0)
 * @param[in] ny          number of grid points in y (> 0)
 * @param[in] nz          number of grid points in z (> 0)
 * @param[in] dx          grid spacing in x
 * @param[in] dy          grid spacing in y
 * @param[in] dz          grid spacing in z
 * @param[in] x0          x origin
 * @param[in] y0          y origin
 * @param[in] z0          z origin
 *
 * @param[out] xlocs      x anchor node locations [nx*ny*nz]
 * @param[out] ylocs      y anchor node locations [nx*ny*nz]
 * @param[out] zlocs      z anchor node locations [nx*ny*nz]
 *
 * @author Ben Baker, ISTI
 *
 */
void __regmesh_makeRegularNodes(int nx, int ny, int nz,
                                double dx, double dy, double dz,
                                double x0, double y0, double z0,
                                double *__restrict__ xlocs,
                                double *__restrict__ ylocs,
                                double *__restrict__ zlocs)
{
    int inpg, ix, iy, iz, nxy;
/*
    #pragma omp parallel for collapse(3) private(inpg, ix, iy, iz) \
     shared (xlocs, ylocs, zlocs)
*/
    nxy = nx*ny;
    #pragma omp simd collapse(3)
    for (iz=0; iz<nz; iz++){
        for (iy=0; iy<ny; iy++){
            for (ix=0; ix<nx; ix++){
                inpg = iz*nxy + iy*nx + ix;
                xlocs[inpg] = x0 + (double) ix*dx;
                ylocs[inpg] = y0 + (double) iy*dy;
                zlocs[inpg] = z0 + (double) iz*dz;
            }
        }
    }
    return;
}
//============================================================================//
/*!
 * @brief Determines anchor node locations in a regular mesh 
 *
 * @param[in] nx          number of grid points in x (> 0)
 * @param[in] ny          number of grid points in y (> 0)
 * @param[in] nz          number of grid points in z (> 0)
 * @param[in] dx          grid spacing in x
 * @param[in] dy          grid spacing in y
 * @param[in] dz          grid spacing in z
 * @param[in] x0          x origin
 * @param[in] y0          y origin
 * @param[in] z0          z origin
 * @param[in] nelem       number of elements in mesh
 *
 * @param[inout] element  on input holds the number of elements and
 *                        space for the (x,y,z) locations [nelem]
 *                        on output contains the positions of each 
 *                        anchor node on all elements 
 *
 * @result number of elements in regular hex mesh
 *
 * @author Ben Baker, ISTI
 *
 */
int regmesh_makeRegularNodes(int nx, int ny, int nz,
                             double dx, double dy, double dz,
                             double x0, double y0, double z0,
                             int nelem,
                             struct mesh_element_struct *element)
{
    const char *fcnm = "regmesh_makeRegularNodes\0";
    double *xlocs, *ylocs, *zlocs;
    int ierr, nnpg;
    nnpg = nx*ny*nz;
    if (nnpg < 1){
        printf("%s: Error no elements in mesh\n", fcnm);
        return -1;
    }
    xlocs = (double *)calloc(nnpg, sizeof(double));
    ylocs = (double *)calloc(nnpg, sizeof(double));
    zlocs = (double *)calloc(nnpg, sizeof(double));
    // Set the points
    __regmesh_makeRegularNodes(nx, ny, nz, 
                               dx, dy, dz, 
                               x0, y0, z0, 
                               xlocs, ylocs, zlocs);
    // Set the anchor node locations
    ierr = mesh_element__setAnchorNodeLocations(nelem, nnpg,
                                                xlocs, ylocs, zlocs,
                                                element); 
    if (ierr != 0){
        printf("%s: Error setting (x,y,z) locations on element\n", fcnm);
    }
    // Free space
    free(xlocs);
    free(ylocs);
    free(zlocs);
    xlocs = NULL;
    ylocs = NULL;
    zlocs = NULL;
    return ierr;
}

