#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
//#include <mpi.h>
#include <cblas.h>
#include <math.h>
#include <iscl/log/log.h>
#include "cvm.h"
#include "iscl/interpolate/interpolate.h"

/*!
 * @brief Utility function for generating SPECFEM3D meshes directly from the
 *        Cascadia Communicate Velocity Model.  The flow of control is as
 *        follows:
 *          (1) Define a subset of the 
 *          (2) 
 *
 * @author Ben Baker, ISTI
 *
 */ 
int main()
{
    const char *fcnm = "cvm_mesher\0";
    char *projnm = "cvm\0";
    char ini_file[PATH_MAX];
    struct cvm_parms_struct parms;
    struct cvm_model_struct *cvm_model;
    struct mesh_struct mesh;
    double x0, y0, z0;
    int ierr, lay, nx, ny, nz;
    // Initialize
    ierr = 0;
    cvm_model = NULL;
    memset(&parms, 0, sizeof(struct cvm_parms_struct));
    memset(ini_file, 0, sizeof(ini_file));
    memset(&mesh, 0, sizeof(mesh));
    // Read the ini file
    log_infoF("%s: Reading ini file...\n", fcnm);
    ierr = cvm_readini(projnm, &parms);
    if (ierr != 0)
    {
        log_errorF("%s: Error reading ini file\n", fcnm);
        goto ERROR;
    }
    // Load the model files into memory
    log_infoF("%s: Loading the CVM...\n", fcnm);
    cvm_model = (struct cvm_model_struct *)
                calloc(parms.nlay_cvm, sizeof(struct cvm_model_struct));
    for (lay=0; lay<parms.nlay_cvm; lay++)
    {
        log_infoF("%s: Loading layer: %d\n", fcnm, lay+1);
        ierr = cvmio_readLayer(lay+1, parms, &cvm_model[lay]);
        if (ierr != 0)
        {
            log_errorF("%s: There was an error reading the CVM\n", fcnm);
            goto ERROR;
        }
        log_infoF("\n");
    }
    // Extend the model so there are no gaps between layers
    ierr = cvm_extendBase(parms.nlay_cvm, cvm_model);
    // Dump the community velocity model 
/*
    log_infoF("%s: Archiving extracted model...\n", fcnm);
    ierr = cvmio_write__h5(parms.cvm_outputdir, projnm,
                           parms.nlay_cvm, cvm_model);
*/
    for (lay=0; lay<parms.nlay_cvm; lay++)
    {
        // Generate the mesh for this layer and dump it

    }
    // Build the mesh
    log_infoF("%s: Generating the mesh...\n", fcnm);
    if (parms.ncoarsen == 0)
    {
        /* TODO: make a regular mesh driver function in regmesh.c */
        log_infoF("%s: Generating regular mesh...\n", fcnm);
        ierr = cvm_estimateNumberOfGridPoints(parms.nlay_cvm, cvm_model,
                                              parms.dx_fem,
                                              parms.dy_fem,
                                              parms.dz_fem,
                                              &nx, &ny, &nz);
        if (ierr != 0)
        {
            log_errorF("%s: Error estimating regular mesh size\n", fcnm);
            goto ERROR;
        }
        log_infoF("%s: Grid points in regular mesh: (%d,%d,%d)\n",
                  fcnm, nx, ny, nz);
        // Make a regular mesh; get the number of elements
        ierr = regmesh_getNumberOfElements(nx, ny, nz, &mesh.nelem);
        mesh.element = (struct mesh_element_struct *)
                       calloc(mesh.nelem, sizeof(struct mesh_element_struct));
        // Because the mesh could be big don't use the redundant storage
        // method which puts material pointers on each element 
        mesh.lptr_only = true;
        // Homogeneous hex mesh
        mesh.lhomog = true;
        // Set the pointers
        log_infoF("%s: Setting regular mesh pointers...\n", fcnm);
        ierr = mesh_element_memory__allocateIntegerPointers(HEX8, mesh.nelem,
                                                            mesh.element);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to set memory for pointers\n", fcnm);
            goto ERROR;
        }
        ierr = regmesh_makeHexMeshPointers(nx, ny, nz, mesh.nelem,
                                           mesh.element);
        if (ierr != 0)
        {
            log_errorF("%s: Error setting pointers!\n", fcnm);
            goto ERROR;
        }
        // Set the anchor node locations
        log_infoF("%s: Setting the regular anchor node locations...\n", fcnm);
        ierr = regmesh_getNumberOfAnchorNodes(nx, ny, nz, &mesh.nnpg);
        if (ierr != 0 || mesh.nnpg < 1)
        {
            if (mesh.nnpg < 1){ierr = 1;}
            log_infoF("%s: Failed to get number of anchor nodes\n", fcnm);
            goto ERROR;
        }
        x0 = cvm_model[0].x0;
        y0 = cvm_model[0].y0;
        z0 = 0.0;
        if (mesh.lptr_only)
        {
            mesh.xlocs = (double *)calloc(mesh.nnpg, sizeof(double));
            mesh.ylocs = (double *)calloc(mesh.nnpg, sizeof(double));
            mesh.zlocs = (double *)calloc(mesh.nnpg, sizeof(double));
            __regmesh_makeRegularNodes(nx, ny, nz,
                                       parms.dx_fem, parms.dy_fem, parms.dz_fem,
                                       x0, y0, z0,
                                       mesh.xlocs,
                                       mesh.ylocs,
                                       mesh.zlocs);
        }
        else
        {
            ierr = regmesh_makeRegularNodes(nx, ny, nz,
                                            parms.dx_fem,
                                            parms.dy_fem,
                                            parms.dz_fem,
                                            x0, y0, z0,
                                            mesh.nelem, mesh.element);
        }
        if (ierr != 0)
        {
            log_errorF("%s: Error creating regular nodes\n", fcnm);
            goto ERROR;
        }
    
    }
    else // Mesh which coarsens with depth
    {
        mesh = layeredMesh_driver(parms,
                                  cvm_model,
                                  &ierr);
        if (ierr != 0)
        {
            log_errorF("%s: Error creating layered mesh!\n", fcnm);
            goto ERROR;
        }
    }
/*
    // Set the pointers 
    log_infoF("%s: Setting mesh pointers...\n", fcnm);
    ierr = mesh_element_memory__allocateIntegerPointers(HEX8, mesh.nelem,
                                                        mesh.element);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to set memory for pointers\n", fcnm);
        goto ERROR;
    }
    ierr = regmesh_makeHexMeshPointers(nx, ny, nz, mesh.nelem, mesh.element);
    if (ierr != 0)
    {
        log_errorF("%s: Error setting pointers!\n", fcnm);
        goto ERROR;
    }
    // Set the anchor node locations
    log_infoF("%s: Setting the anchor node locations...\n", fcnm);
    ierr = regmesh_getNumberOfAnchorNodes(nx, ny, nz, &mesh.nnpg);
    if (ierr != 0 || mesh.nnpg < 1)
    {
        if (mesh.nnpg < 1){ierr = 1;}
        log_infoF("%s: Failed to get number of anchor nodes\n", fcnm);
        goto ERROR;
    }
    x0 = cvm_model[0].x0;
    y0 = cvm_model[0].y0;
    z0 = 0.0;
    if (mesh.lptr_only)
    {
        mesh.xlocs = (double *)calloc(mesh.nnpg, sizeof(double));
        mesh.ylocs = (double *)calloc(mesh.nnpg, sizeof(double));
        mesh.zlocs = (double *)calloc(mesh.nnpg, sizeof(double));
        __regmesh_makeRegularNodes(nx, ny, nz,
                                   parms.dx_fem, parms.dy_fem, parms.dz_fem,
                                   x0, y0, z0,
                                   mesh.xlocs,
                                   mesh.ylocs,
                                   mesh.zlocs);
    }
    else
    {
        ierr = regmesh_makeRegularNodes(nx, ny, nz,
                                   parms.dx_fem, parms.dy_fem, parms.dz_fem,
                                   x0, y0, z0,
                                   mesh.nelem, mesh.element);
        if (ierr != 0){
            log_errorF("%s: Error creating regular nodes\n", fcnm);
            goto ERROR;
        }
    }
*/
    // Set the material properties (put this before the topography deformation)
    log_infoF("%s: Setting materials...\n", fcnm);
    ierr = cvm_cvm2meshMaterials(parms.nlay_cvm, cvm_model, &mesh);
    if (ierr != 0){
        log_errorF("%s: Error setting material properties\n", fcnm);
        goto ERROR;
    }
    // Reverse the orientation of the mesh
    ierr = cvm_cvm2freeSurfaceToZero(&mesh);
    if (ierr != 0){
        log_errorF("%s: Error reversing mesh\n", fcnm);
        goto ERROR;
    }
    // Dump a NLL grid
/*
    log_infoF("%s: Writing the NLL mesh...\n", fcnm);
printf("%f %f %f %f\n", cvm_model[0].x0, cvm_model[0].y0, x0, y0);
    ierr = meshio_write__NLLGrid(parms.nll_outputdir, projnm,
                                nx - 1, ny - 1, nz - 1,
                                x0, y0,
                                parms.dx_fem, parms.dy_fem, parms.dz_fem,
                                parms.utm_zone, 0.0,
                                mesh);
return 0;
*/
    // Add in topography?
    if (parms.ltopo)
    {
        log_infoF("%s: Mapping topography...\n", fcnm);

        ierr = topo30_deformMesh(parms.topofl, parms.utm_zone,
                                 parms.ztopo_min, 
                                 &mesh); 
        if (ierr != 0)
        {
            log_errorF("%s: Error deforming mesh\n", fcnm);
        }
    }
    log_infoF("%s: Writing h5 mesh...\n", fcnm);
    ierr = meshio_write__h5(parms.mesh_outputdir, projnm, mesh);
    log_infoF("%s: Writing specfem mesh...\n", fcnm);
    ierr = meshio_write__specfem3d(parms.mesh_outputdir, projnm, mesh);
    if (ierr != 0){
        log_errorF("%s: Error writing SPECFEM3D files!\n", fcnm);
        goto ERROR;
    }
/*
struct mesh_element_struct *element = NULL;
 ierr = regmesh_getNumberOfElements(nx, ny, nz, &nelem);
 element = (struct mesh_element_struct *)
           calloc(nelem, sizeof(struct mesh_element_struct));
 ierr = mesh_element__allocate(HEX8, nelem, element);
 ierr = regmesh_makeHexMeshPointers(nx, ny, nz, nelem, element);
 ierr = regmesh_getNumberOfAnchorNodes(nx, ny, nz, &mesh.nnpg);
 ierr = regmesh_makeRegularNodes(nx, ny, nz,
                                 dx, dy, dz,
                                 x0, y0, z0,
                                 nelem, element);
 
int meshio_write__h5(char *meshdir, char *projnm,
                     int nelem, struct mesh_element_struct *element);
 meshio_write__h5(parms.mesh_outputdir, projnm,
                       nelem, element);
 mesh_element__free(nelem, element);
*/

    // If topography is present then load the topo30 model
/*
    if (parms.ltopo){
        log_infoF("%s: Mapping topography...\n", fcnm);
    }
*/

ERROR:; // Break ahead for error
    if (ierr != 0){log_errorF("%s: Error occurred\n", fcnm);}
    // Free memory
    mesh_element_memory__free(mesh.nelem, mesh.element);
    mesh_memory__free(&mesh);
    cvm_memory_free__model(parms.nlay_cvm, cvm_model);
    cvm_memory_free__parms(&parms);
    return 0;
}

/*!
 * @brief Extends the CVM so that there are no `spaces' in between layers
 *
 * @param[in] nlay     number of layers in model
 *
 * @param[inout] cvm_model   on input holds the cvm_model read from disk
 *                           on exit holds the cvm_model s.t. the bases of
 *                           first and second layers are extended to the 
 *                           top of the second and third layers respectively
 *                           [nlay]
 *
 * @author Ben Baker, ISTI
 *
 */
int cvm_extendBase(int nlay, struct cvm_model_struct *cvm_model)
{
    const char *fcnm = "cvm_extendBase\0";
    double *xtop, *ytop, *xbot, *ybot, *ref;
    int ierr, ilay, indx, ix, iy, jndx, nx_bot, nx_top, nxy_bot, nxy_top,
        nxyz_top, ny_bot, ny_top;
    // While unlikely we avoid spacial aliasing by interpolating the fine 
    // grid on the coarse grid
    for (ilay=0; ilay<nlay-1; ilay++){
        // Get the top
        nx_top = cvm_model[ilay].nx;
        ny_top = cvm_model[ilay].ny;
        nxy_top = nx_top*ny_top;
        nxyz_top = cvm_model[ilay].nz*nxy_top;
        cvm_model[ilay].npts = cvm_model[ilay].npts + nxy_top;
        cvm_model[ilay].nz = cvm_model[ilay].nz + 1; // add layer
        cvm_model[ilay].z1 = cvm_model[ilay+1].z0;
        xtop = (double *)calloc(nx_top, sizeof(double));
        ytop = (double *)calloc(ny_top, sizeof(double));
        for (ix=0; ix<nx_top; ix++){
            xtop[ix] = cvm_model[ilay].x0 + cvm_model[ilay].dx*(double) ix;
        }
        for (iy=0; iy<ny_top; iy++){
            ytop[iy] = cvm_model[ilay].y0 + cvm_model[ilay].dy*(double) iy;
        }
        // Get the bottom
        nx_bot = cvm_model[ilay+1].nx;
        ny_bot = cvm_model[ilay+1].ny;
        nxy_bot = nx_bot*ny_bot;
        xbot = (double *)calloc(nx_bot, sizeof(double));
        ybot = (double *)calloc(ny_bot, sizeof(double));
        for (ix=0; ix<nx_bot; ix++){
            xbot[ix] = cvm_model[ilay+1].x0 + cvm_model[ilay+1].dx*(double) ix;
        }
        for (iy=0; iy<ny_bot; iy++){
            ybot[iy] = cvm_model[ilay+1].y0 + cvm_model[ilay+1].dy*(double) iy;
        }
        // Interpolate vp, vs, density
        ref = (double *)calloc(nxyz_top, sizeof(double));
        // Resize vp
        cblas_dcopy(nxyz_top, cvm_model[ilay].vp, 1, ref, 1);
        free(cvm_model[ilay].vp);
        cvm_model[ilay].vp = (double *)
                             calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].vp, 1);
        // Resize vs
        cblas_dcopy(nxyz_top, cvm_model[ilay].vs, 1, ref, 1); 
        free(cvm_model[ilay].vs);
        cvm_model[ilay].vs = (double *)
                             calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].vs, 1);
        // Resize density 
        cblas_dcopy(nxyz_top, cvm_model[ilay].dens, 1, ref, 1);
        free(cvm_model[ilay].dens);
        cvm_model[ilay].dens = (double *)
                               calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].dens, 1);
        // Resize Qp 
        cblas_dcopy(nxyz_top, cvm_model[ilay].Qp, 1, ref, 1); 
        free(cvm_model[ilay].Qp);
        cvm_model[ilay].Qp = (double *)
                             calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].Qp, 1); 
        // Resize Qs
        cblas_dcopy(nxyz_top, cvm_model[ilay].Qs, 1, ref, 1); 
        free(cvm_model[ilay].Qs);
        cvm_model[ilay].Qs = (double *)
                             calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].Qs, 1); 
        // Resize x nodes
        cblas_dcopy(nxyz_top, cvm_model[ilay].xlocs, 1, ref, 1); 
        free(cvm_model[ilay].xlocs);
        cvm_model[ilay].xlocs = (double *)
                                calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].xlocs, 1);
        // Resize y nodes
        cblas_dcopy(nxyz_top, cvm_model[ilay].ylocs, 1, ref, 1);
        free(cvm_model[ilay].ylocs);
        cvm_model[ilay].ylocs = (double *)
                                calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].ylocs, 1);
        // Resize z nodes
        cblas_dcopy(nxyz_top, cvm_model[ilay].zlocs, 1, ref, 1);
        free(cvm_model[ilay].zlocs);
        cvm_model[ilay].zlocs = (double *)
                                calloc(cvm_model[ilay].npts, sizeof(double));
        cblas_dcopy(nxyz_top, ref, 1, cvm_model[ilay].zlocs, 1);
        // Interpolate vp
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].vp[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].vp[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0){
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate vs
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].vs[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].vs[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0){ 
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate density
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].dens[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].dens[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0){ 
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate Qp
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].Qp[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].Qp[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0){
            log_errorF("%s: Error interpolating Qp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate Qs
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].Qs[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].Qs[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0){
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Fix nodes
        indx = 0;
        jndx = nxyz_top;
        for (iy=0; iy<ny_top; iy++){
            for (ix=0; ix<nx_top; ix++){
                cvm_model[ilay].xlocs[jndx] = cvm_model[ilay].xlocs[indx];
                cvm_model[ilay].ylocs[jndx] = cvm_model[ilay].ylocs[indx];
                cvm_model[ilay].zlocs[jndx] = cvm_model[ilay+1].z0; 
                indx = indx + 1;
                jndx = jndx + 1;
            }
        }
        // Free space
        free(ref);
        free(xbot);
        free(ybot);
        free(xtop);
        free(ytop);
        ref = NULL;
        xbot = NULL;
        ybot = NULL;
        xtop = NULL;
        ytop = NULL;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Estimates number of grid points for making finite element mesh 
 *
 * @param[in] nlay        number of layers in CVM (3)
 * @param[in] cvm_model   holds the layers in the Cascadia CVM [nlay]
 * @param[in] dx_fem      grid spacing (m) in x for finite element mesh
 * @param[in] dy_fem      grid spacing (m) in y for finite element mesh
 * @param[in] dz_fem      grid spacing (m) in z for finite element mesh 
 *
 * @param[out] nx         number of x grid points in finite element mesh 
 * @param[out] ny         number of y grid points in finite element mesh 
 * @param[out] nz         number of z grid points in finite element mesh
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 */
int cvm_estimateNumberOfGridPoints(int nlay, struct cvm_model_struct *cvm_model,
                                   double dx_fem, double dy_fem, double dz_fem,
                                   int *nx, int *ny, int *nz)
{
    const char *fcnm = "cvm_estimateNumberOfGridPoints\0";
    double xoff, yoff, zoff;
    int ierr, lay, nx0, ny0;
    //------------------------------------------------------------------------//
    ierr = 0;
    lay = 0;
    xoff = (double) (cvm_model[lay].nx - 1)*cvm_model[lay].dx;
    yoff = (double) (cvm_model[lay].ny - 1)*cvm_model[lay].dy;
    *nx = (int) (xoff/dx_fem + 0.5) + 1;
    if ((double) (*nx - 1)*dx_fem > xoff){*nx = *nx - 1;}
    *ny = (int) (yoff/dy_fem + 0.5) + 1;
    if ((double) (*ny - 1)*dy_fem > yoff){*ny = *ny - 1;}
    nx0 = *nx;
    ny0 = *ny;
    for (lay=1; lay<nlay; lay++)
    {
        xoff = (double) (cvm_model[lay].nx - 1)*cvm_model[lay].dx;
        yoff = (double) (cvm_model[lay].ny - 1)*cvm_model[lay].dy;
        *nx = (int) (xoff/dx_fem + 0.5) + 1;
        if ((double) (*nx - 1)*dx_fem > xoff){*nx = *nx - 1;}
        *ny = (int) (yoff/dy_fem + 0.5) + 1;
        if ((double) (*ny - 1)*dy_fem > yoff){*ny = *ny - 1;}
        if (*nx != nx0)
        {
            log_errorF("%s: There is an inconsistency in x\n", fcnm);
            ierr = 1;
            break;
        }
        if (*ny != ny0)
        {
            log_errorF("%s: There is an inconsistency in y\n", fcnm);
            ierr = 1;
            break;
        }
    }
    // Tally up the z grid points
    zoff = cvm_model[nlay-1].z1 - cvm_model[0].z0;
    *nz = (int) (zoff/dz_fem + 0.5) + 1;
    if ((double) (*nz - 1)*dz_fem > zoff){*nz = *nz - 1;}
    return ierr; 
}
//============================================================================//
/*!
 * @brief Modifies the mesh so that the free surface is at z=0.
 *
 * @param[inout] mesh      on input holds the mesh z locations where the free
 *                         surface is located at max(mesh->zlocs)
 *                         on output holds the mesh z locations are at 0
 *
 * @result 0 indicates success
 *
 */
int cvm_cvm2freeSurfaceToZero(struct mesh_struct *mesh)
{
    const char *fcnm = "cvm_cvm2freeSurfaceToZero\0";
    double zmax_mesh, zmin_mesh;
    int inpg;
    // Get the max/min of the mesh
    zmax_mesh = mesh->zlocs[0];
    zmin_mesh = mesh->zlocs[0];
    for (inpg=1; inpg<mesh->nnpg; inpg++)
    {
        zmax_mesh = fmax(mesh->zlocs[inpg], zmax_mesh);
        zmin_mesh = fmin(mesh->zlocs[inpg], zmin_mesh);
    }
    log_infoF("%s: Min and max points in mesh: %f %f\n",
              fcnm, zmin_mesh, zmax_mesh);
    for (inpg=0; inpg<mesh->nnpg; inpg++)
    {
        mesh->zlocs[inpg] = mesh->zlocs[inpg] - zmax_mesh;
    }
    return 0;
}
//============================================================================//
/*!
 *
 * @brief Interpolates the CVM model properties onto the mesh
 *
 */
int cvm_cvm2meshMaterials(int nlay, struct cvm_model_struct *cvm_model,
                          struct mesh_struct *mesh)
{
    const char *fcnm = "cvm_cvm2meshMaterials\0";
    double *dens, *Qp, *Qs, *vp, *vq, *vs, *x, *xq, *y, *yq, *z, *zq,
           zmax, zmin, zmax_mesh, zmin_mesh;
    int *idest, ierr, inpg, ixyzq, ix, iy, iz, lay, nxyzq, nx, nxyz, ny, nz;
    bool *linit;
    const double eps = 1.e-6; // tolerance at micrometer level
    bool lflip = true; // FEM is +up but CVM is +down
    bool lqptr_grid = false; // Interpolation points won't be gridded
    //------------------------------------------------------------------------//
    ierr = 0;
    if (mesh->nnpg < 1){
        log_errorF("%s: Error no anchor nodes in mesh!\n", fcnm);
        return -1;
    }
    // Set space
    linit = (bool *)calloc(mesh->nnpg, sizeof(bool));
    mesh->vp = (double *)calloc(mesh->nnpg, sizeof(double));
    mesh->vs = (double *)calloc(mesh->nnpg, sizeof(double));
    mesh->dens = (double *)calloc(mesh->nnpg, sizeof(double));
    mesh->Qp = (double *)calloc(mesh->nnpg, sizeof(double));
    mesh->Qs = (double *)calloc(mesh->nnpg, sizeof(double));
    // Get the max/min of the mesh
    zmax_mesh = mesh->zlocs[0];
    zmin_mesh = mesh->zlocs[0];
    for (inpg=1; inpg<mesh->nnpg; inpg++){
        zmax_mesh = fmax(mesh->zlocs[inpg], zmax_mesh);
        zmin_mesh = fmin(mesh->zlocs[inpg], zmin_mesh);
    }
    log_infoF("%s: Min and max points in mesh: %f %f\n",
              fcnm, zmin_mesh, zmax_mesh);
    // Loop on layers
    for (lay=0; lay<nlay; lay++){
        // Initialize
        x = NULL;
        y = NULL;
        z = NULL;
        xq = NULL;
        yq = NULL;
        xq = NULL;
        zq = NULL;
        vp = NULL;
        vs = NULL;
        dens = NULL;
        idest = NULL;
        nx = cvm_model[lay].nx;
        ny = cvm_model[lay].ny;
        nz = cvm_model[lay].nz;
        nxyz = nx*ny*nz;
        x = (double *)calloc(nx, sizeof(double));
        y = (double *)calloc(ny, sizeof(double));
        z = (double *)calloc(nz, sizeof(double));
        // Get the x, y, z locations
        #pragma omp simd
        for (ix=0; ix<nx; ix++){
            x[ix] = cvm_model[lay].x0 + (double) ix*cvm_model[lay].dx;
        }
        #pragma omp simd
        for (iy=0; iy<ny; iy++){
            y[iy] = cvm_model[lay].y0 + (double) iy*cvm_model[lay].dy;
        }
//printf("%d %f %f\n\n", nz, cvm_model[lay].z0, cvm_model[lay].dz);
        for (iz=0; iz<nz; iz++){
            z[iz] = cvm_model[lay].z0 + (double) (nz-1-iz)*cvm_model[lay].dz;
            if (iz == 0){
                z[iz] = cvm_model[lay].z1;
            }
//printf("---------%f %f\n",z[iz], (zmax_mesh - zmin_mesh) - z[iz]);
            z[iz] = (zmax_mesh - zmin_mesh) - z[iz];
        }
        // Set the material properties
        vp = (double *)calloc(nxyz, sizeof(double));
        vs = (double *)calloc(nxyz, sizeof(double));
        dens = (double *)calloc(nxyz, sizeof(double));
        Qp = (double *)calloc(nxyz, sizeof(double));
        Qs = (double *)calloc(nxyz, sizeof(double));
        __regmesh_copyRegularModel(lflip, nx, ny, nz,
                                   cvm_model[lay].vp,
                                   cvm_model[lay].vs,
                                   cvm_model[lay].dens,
                                   cvm_model[lay].Qp,
                                   cvm_model[lay].Qs,
                                   vp, vs, dens, Qp, Qs);
        // Get the (x, y, z) anchor nodes in this layer
        zmin = z[0] - eps;
        zmax = z[nz-1] + eps;
        nxyzq = 0;
        for (inpg=0; inpg<mesh->nnpg; inpg++)
        {
            if (linit[inpg]){continue;}
            if (zmin <= mesh->zlocs[inpg] && mesh->zlocs[inpg] <= zmax)
            {
                nxyzq = nxyzq + 1;
            }
        }
        // There are points to interpolate so interpolate them
        if (nxyzq > 0)
        {
            log_infoF("%s: Interpolating %d points in layer %d...\n",
                      fcnm, nxyzq, lay+1);
            idest = (int *)calloc(nxyzq, sizeof(int));
            xq = (double *)calloc(nxyzq, sizeof(double));
            yq = (double *)calloc(nxyzq, sizeof(double));
            zq = (double *)calloc(nxyzq, sizeof(double));
            vq = (double *)calloc(nxyzq, sizeof(double));
            ixyzq = 0;
            for (inpg=0; inpg<mesh->nnpg; inpg++){
                if (linit[inpg]){continue;}
                if (zmin <= mesh->zlocs[inpg] && mesh->zlocs[inpg] <= zmax){
                    xq[ixyzq] = mesh->xlocs[inpg];
                    yq[ixyzq] = mesh->ylocs[inpg];
                    zq[ixyzq] = mesh->zlocs[inpg];
                    idest[ixyzq] = inpg;
                    ixyzq = ixyzq + 1;
                }
            }
            // Interpolate vp
            ierr = __interpolate_interp3d(nx, x, ny, y, nz, z,
                                          nxyz, vp,
                                          nxyzq, xq, nxyzq, yq, nxyzq, zq,
                                          nxyzq, vq,
                                          TRILINEAR, lqptr_grid); 
            if (ierr != 0){
                log_errorF("%s: Error interpolating vp!\n", fcnm);
                return -1;
            }
            // Copy
            for (ixyzq=0; ixyzq<nxyzq; ixyzq++){
                inpg = idest[ixyzq];
                mesh->vp[inpg] = vq[ixyzq];
                linit[inpg] = true; 
            }
            // Interpolate vs
            ierr = __interpolate_interp3d(nx, x, ny, y, nz, z,
                                          nxyz, vs,
                                          nxyzq, xq, nxyzq, yq, nxyzq, zq,
                                          nxyzq, vq,
                                          TRILINEAR, lqptr_grid);
            if (ierr != 0){ 
                log_errorF("%s: Error interpolating vs!\n", fcnm);
                return -1; 
            }
            // Copy
            for (ixyzq=0; ixyzq<nxyzq; ixyzq++){
                inpg = idest[ixyzq];
                mesh->vs[inpg] = vq[ixyzq];
                linit[inpg] = true; 
            }
            // Interpolate density
            ierr = __interpolate_interp3d(nx, x, ny, y, nz, z,
                                          nxyz, dens,
                                          nxyzq, xq, nxyzq, yq, nxyzq, zq,
                                          nxyzq, vq,
                                          TRILINEAR, lqptr_grid);
            if (ierr != 0){ 
                log_errorF("%s: Error interpolating density!\n", fcnm);
                return -1; 
            }
            // Copy
            for (ixyzq=0; ixyzq<nxyzq; ixyzq++){
                inpg = idest[ixyzq];
                mesh->dens[inpg] = vq[ixyzq];
                linit[inpg] = true; 
            }
            // Interpolate Qp
            ierr = __interpolate_interp3d(nx, x, ny, y, nz, z,
                                          nxyz, Qp,
                                          nxyzq, xq, nxyzq, yq, nxyzq, zq, 
                                          nxyzq, vq, 
                                          TRILINEAR, lqptr_grid);
            if (ierr != 0){ 
                log_errorF("%s: Error interpolating Qp!\n", fcnm);
                return -1; 
            }
            // Copy
            for (ixyzq=0; ixyzq<nxyzq; ixyzq++){
                inpg = idest[ixyzq];
                mesh->Qp[inpg] = vq[ixyzq];
                linit[inpg] = true; 
            }
             // Interpolate Qs
            ierr = __interpolate_interp3d(nx, x, ny, y, nz, z,
                                          nxyz, Qs,
                                          nxyzq, xq, nxyzq, yq, nxyzq, zq, 
                                          nxyzq, vq, 
                                          TRILINEAR, lqptr_grid);
            if (ierr != 0){ 
                log_errorF("%s: Error interpolating Qs!\n", fcnm);
                return -1; 
            }
            // Copy
            for (ixyzq=0; ixyzq<nxyzq; ixyzq++){
                inpg = idest[ixyzq];
                mesh->Qs[inpg] = vq[ixyzq];
                linit[inpg] = true; 
            }
            // Free space
            free(vq);
            vq = NULL;

        }
        // Free space
        if (vp != NULL){free(vp);}
        if (vs != NULL){free(vs);}
        if (dens != NULL){free(dens);}
        if (Qp != NULL){free(Qp);}
        if (Qs != NULL){free(Qs);}
        if (x != NULL){free(x);}
        if (y != NULL){free(y);}
        if (z != NULL){free(z);}
        if (xq != NULL){free(xq);}
        if (yq != NULL){free(yq);}
        if (zq != NULL){free(zq);}
        if (idest != NULL){free(idest);}
        vp = NULL;
        vs = NULL;
        dens = NULL;
        Qp = NULL;
        Qs = NULL;
        x = NULL;
        y = NULL;
        z = NULL;
        xq = NULL;
        yq = NULL;
        xq = NULL;
        idest = NULL;
    }
    // Check for un-initialized nodes
    if (ierr == 0){
        for (inpg=0; inpg<mesh->nnpg; inpg++)
        {
            if (!linit[inpg])
            {
                ierr = ierr + 1;
            }
        }
        if (ierr > 0)
        {
            log_infoF("%s: Failed to initialize %d points\n", fcnm, ierr);
            ierr = 1;
        }
    }
    // Free memory
    free(linit);
    linit = NULL;
    return ierr;
}
