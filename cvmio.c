#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <iscl/os/os.h>
#include <iscl/log/log.h>
#include "cvm.h"

#define VP_DEFAULT 2300.0
#define VS_DEFAULT 1333.0
#define VMAX_CVM 10000.0
#define VMIN_CVM 1.0
#define NAN_CVM 1.e20
#define IS_BIG_ENDIAN (*(uint16_t *)"\0\xff" < 0x100)

float __cvm_io_qcFloat_vp(float valin); //, float val0);
float __cvm_io_qcFloat_vs(float valin); //, float val0);
void cvmio_write__xdmf(char *dirnm, char *projnm,
                       int nlay, int *nelem, int *nnpg);

/*!
 * @brief Loads a layer into into memory.  This function will apply all 
 *        thresholding to the model.
 *
 * @param[in] lay            layer [1, 2, or 3]
 * @param[in] parms          parameter structure
 *
 * @param[out] cvm_model     holds the regularly sampled velocity model
 *                           for this layer bounding the computational
 *                           domain 
 *                           
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int cvmio_readLayer(int lay,
                    struct cvm_parms_struct parms,
                    struct cvm_model_struct *cvm_model)
{
    char *fcnm = "cvmio_readLayer\0";
    char root[PATH_MAX], vpfl[PATH_MAX], vsfl[PATH_MAX], app[64];
    char *model_vp, *model_vs;
    double common_mult, lat0, lat1, lon0, lon1, utmx, utmy, vpvs, vsvp,
           z, zbeg, zlay_thickness;
    float pval, pval0, sval, sval0;
    int indx_cvm, indx_mod, ix, ixl, ixmax, ixmin, iy, iyl, iymax, iymin,
        iz, izmax, izmin, jx, jy, jz, nbytes, nx, nx_mod, ny,
        ny_mod, nz, nz_mod, nxyz;
    int iway = 1; // UTM -> (lat,lon)
    size_t lenos;
    bool lswap;
    bool suppress = false;
    // Null out the model
    memset(cvm_model, 0, sizeof(struct cvm_model_struct));
    // Set the model file names
    memset(root, 0, sizeof(root));
    memset(vpfl, 0, sizeof(vpfl));
    memset(vsfl, 0, sizeof(vsfl));
    strcpy(root, parms.cvm_moddir); 
    lenos = strlen(root);
    if (lenos > 0)
    {
        if (root[lenos-1] != '/'){strcat(root, "/\0");}
    }
    else
    {
        strcpy(root, "./\0");
    }
    // Vp file name
    strcpy(vpfl, root);
    memset(app, 0, sizeof(app));
    sprintf(app, "Vpl%dc.bin", lay);
    strcat(vpfl, app);
    // Vs file name
    strcpy(vsfl, root);
    memset(app, 0, sizeof(app));
    sprintf(app, "Vsl%dc.bin", lay);
    strcat(vsfl, app);
    // Detect my endianness
    lswap = false; // Assume this is a little endian machine
    if (IS_BIG_ENDIAN){lswap = true;}
    // Get the sizes 
    nx = parms.nxl_cvm[lay-1];
    ny = parms.nyl_cvm[lay-1];
    nz = parms.nzl_cvm[lay-1];
    nxyz = nx*ny*nz;
    // Read the Vp model
    model_vp = cvmio_readBinaryFile(vpfl, &nbytes);
    if (nbytes < 1)
    {
        printf("%s: Error loading file %s\n", fcnm, vpfl);
        return -1;
    }
    if (nbytes != 4*nxyz)
    {
        printf("%s: Failed to estimate byte size for Vp %d %d\n",
               fcnm, nbytes, 4*nxyz);
        return -1; 
    }
    // Read the Vs model
    model_vs = cvmio_readBinaryFile(vsfl, &nbytes);
    if (nbytes < 1)
    {
        printf("%s: Error loading file %s\n", fcnm, vsfl);
        return -1;
    }
    if (nbytes != 4*nxyz)
    {
        printf("%s: Failed to estimate byte size for Vs %d %d\n",
               fcnm, nbytes, 4*nxyz);
        return -1;
    }
    // Get the min location remembering we want the corners for all
    // layers aligned
    ixmin = 0;
    ixmax = parms.nxl_cvm[lay-1] - 1;
    iymin = 0;
    iymax = parms.nyl_cvm[lay-1] - 1;
    izmin = 0;
    izmax = parms.nzl_cvm[lay-1] - 1;
    common_mult = parms.common_mult;
    // Get the lower left corner
    ixl = (int) (common_mult/parms.dxl_cvm[lay-1] + 0.5);
    iyl = (int) (common_mult/parms.dyl_cvm[lay-1] + 0.5);
    for (ixmin=0; ixmin<parms.nxl_cvm[lay-1]; ixmin=ixmin+ixl)
    {
        utmx = parms.utm_x0_cvm + (double) (ixmin + ixl)*parms.dxl_cvm[lay-1];
        if (utmx > parms.utm_x0){break;}
    }
    for (iymin=0; iymin<parms.nyl_cvm[lay-1]; iymin=iymin+iyl)
    {
        utmy = parms.utm_y0_cvm + (double) (iymin + iyl)*parms.dyl_cvm[lay-1];
        if (utmy > parms.utm_y0){break;}
    }
    for (ixmax=ixmin; ixmax<parms.nxl_cvm[lay-1]; ixmax=ixmax+ixl)
    {
        utmx = parms.utm_x0_cvm + (double) ixmax*parms.dxl_cvm[lay-1];
        if (utmx > parms.utm_x1){break;}
    }
    for (iymax=iymin; iymax<parms.nyl_cvm[lay-1]; iymax=iymax+iyl)
    {
        utmy = parms.utm_y0_cvm + (double) iymax*parms.dyl_cvm[lay-1];
        if (utmy > parms.utm_y1){break;}
    }
    if (ixmin > parms.nxl_cvm[lay-1] || ixmax > parms.nxl_cvm[lay-1])
    {
        printf("%s: Error computing ixmin/ixmax %d %d\n",
        fcnm, iymin, iymax);
        return -1; 
    }
    if (iymin > parms.nyl_cvm[lay-1] || iymax > parms.nyl_cvm[lay-1])
    {
        printf("%s: Error computing iymin/iymax %d %d\n",
               fcnm, iymin, iymax);
        return -1;
    }
    // Figure out the max depth (we need the shallowest layers)
    if (parms.zmin > 0.0)
    {
        printf("%s: I'm skipping the zmin issue for now\n", fcnm);
    }
    zlay_thickness = (double) (parms.nzl_cvm[lay-1] - 1)*parms.dzl_cvm[lay-1];
    if (parms.zmax < parms.z0_cvm[lay-1] + zlay_thickness)
    {
        zbeg = parms.z0_cvm[lay-1];
        for (izmax=izmin; izmax<parms.nzl_cvm[lay-1]; izmax++)
        {
            z = zbeg + (double) izmax*parms.dzl_cvm[lay-1];
            if (z > parms.zmax)
            {
                break;
            }
        }
    }
    // For the sake of the user tell them where we are extracting the model
    utmx = parms.utm_x0_cvm + (double) ixmin*parms.dxl_cvm[lay-1];
    utmy = parms.utm_y0_cvm + (double) iymin*parms.dyl_cvm[lay-1];
    cvm_model->x0 = utmx;
    cvm_model->y0 = utmy;
    cvm_model->z0 = parms.z0_cvm[lay-1];
    utm_geo(&lon0, &lat0, &utmx, &utmy, &parms.utmzone_cvm, &iway, &suppress);
    if (lon0 < 0.0){lon0 = lon0 + 360.0;} 
    utmx = parms.utm_x0_cvm + (double) ixmax*parms.dxl_cvm[lay-1];
    utmy = parms.utm_y0_cvm + (double) iymax*parms.dyl_cvm[lay-1];
    utm_geo(&lon1, &lat1, &utmx, &utmy, &parms.utmzone_cvm, &iway, &suppress);
    if (lon1 < 0.0){lon1 = lon1 + 360.0;}
    printf("%s: Lower left corner of bounding box: (%f,%f)\n",
              fcnm, lat0, lon0);
    printf("%s: Upper right corner of bounding box: (%f,%f)\n",
              fcnm, lat1, lon1);
    //printf("%d %d %d %d %d %d\n", ixl, iyl, ixmin, iymin, ixmax, iymax);
    // Initializations for Art's quality control
    pval0 = VP_DEFAULT;
    sval0 = VS_DEFAULT;
    // Set the space
    nx_mod = ixmax - ixmin + 1;
    ny_mod = iymax - iymin + 1;
    nz_mod = izmax - izmin + 1;
    cvm_model->npts = nx_mod*ny_mod*nz_mod;
    // 2 x 2 x 2 is bad
    if (cvm_model->npts < 8)
    {
        printf("%s: Error insufficient number of points in model\n", fcnm);
        return -1;
    }
    cvm_model->nx = nx_mod;
    cvm_model->ny = ny_mod;
    cvm_model->nz = nz_mod;
    cvm_model->dx = parms.dxl_cvm[lay-1];
    cvm_model->dy = parms.dyl_cvm[lay-1];
    cvm_model->dz = parms.dzl_cvm[lay-1];
    cvm_model->z1 = cvm_model->z0
                  + (double) (izmax - izmin + 1 - 1)*cvm_model->dz;
    cvm_model->vp   = (double *)
                      calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->vs   = (double *)
                      calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->dens = (double *)
                      calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->Qp   = (double *)
                      calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->Qs   = (double *)
                      calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->xlocs = (double *)
                       calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->ylocs = (double *)
                       calloc((size_t) cvm_model->npts, sizeof(double));
    cvm_model->zlocs = (double *)
                       calloc((size_t) cvm_model->npts, sizeof(double));
    printf("%s: There will be %d points in extracted model\n",
           fcnm, cvm_model->npts);
    // Unpack the CVM binary file. It appears to be packed from deep to shallow.
    // We want the vp and vs models packed with the orientation that z increases
    // down so indx_cvm is reversed like as is done in Art Frankel's original
    // Fortran code
    printf("%s: Unpacking models...\n", fcnm);
    indx_mod = 0;
    jz = 0;
    for (iz=izmin; iz<=izmax; iz++)
    {
        jy = 0;
        for (iy=iymin; iy<=iymax; iy++)
        {
            jx = 0;
            for (ix=ixmin; ix<=ixmax; ix++)
            {
                indx_cvm = nx*ny*(nz - 1 - jz) + nx*iy + ix;
                pval = unpack_float(&model_vp[4*indx_cvm], lswap);
                sval = unpack_float(&model_vs[4*indx_cvm], lswap);
                // Quality control
                pval = __cvm_io_qcFloat_vp(pval); //, pval0);
                sval = __cvm_io_qcFloat_vs(sval); //, sval0);
                // Update
                pval0 = pval;
                sval0 = sval;
                // Fill material vectors
                cvm_model->vp[indx_mod] = (double) pval0;
                cvm_model->vs[indx_mod] = (double) sval0;
                cvm_model->xlocs[indx_mod] = cvm_model->x0
                                           + cvm_model->dx*(double) jx;
                cvm_model->ylocs[indx_mod] = cvm_model->y0
                                           + cvm_model->dy*(double) jy; 
                cvm_model->zlocs[indx_mod] = cvm_model->z0
                                           + cvm_model->dz*(double) jz;
                indx_mod = indx_mod + 1;
                jx = jx + 1;
            } // Loop on x
            jy = jy + 1;
        } // Loop on y
        jz = jz + 1;
    } // Loop on z
    if (indx_mod != cvm_model->npts)
    {
        printf("%s: Failed to initialize points\n", fcnm);
        return -1;
    } 
    // Apply clipping
    if (parms.lthresh_vp)
    {
        printf("%s: Clipping Vp to range [%f, %f]\n",
               fcnm, parms.vp_min, parms.vp_max);
        #pragma omp simd
        for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
        {
            cvm_model->vp[indx_mod] = fmax(cvm_model->vp[indx_mod],
                                           parms.vp_min);
            cvm_model->vp[indx_mod] = fmin(cvm_model->vp[indx_mod],
                                           parms.vp_max);
        }
    }
    if (parms.lthresh_vs)
    {
        printf("%s: Clipping Vs to range [%f, %f]\n",
               fcnm, parms.vs_min, parms.vs_max);
        #pragma omp simd
        for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
        {
            cvm_model->vs[indx_mod] = fmax(cvm_model->vs[indx_mod],
                                           parms.vs_min);
            cvm_model->vs[indx_mod] = fmin(cvm_model->vs[indx_mod],
                                           parms.vs_max);
        }
    }
    if (parms.setvp_from_vs || parms.setvs_from_vp)
    {
        if (parms.setvp_from_vs){
            printf("%s: Setting Vp = %f*Vs\n", fcnm, parms.vpvs_ratio);
            vpvs = parms.vpvs_ratio;
            #pragma omp simd
            for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
            {
                cvm_model->vp[indx_mod] = cvm_model->vs[indx_mod]*vpvs;
            }
        }
        else
        {
            printf("%s: Setting Vs = Vp/%f\n", fcnm, parms.vpvs_ratio);
            vsvp = 1.0/parms.vpvs_ratio;
            #pragma omp simd
            for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
            {
                cvm_model->vs[indx_mod] = cvm_model->vp[indx_mod]*vsvp;
            }
        } 
    // Could be thresholding
    }
    else
    {
        if (parms.lthresh_vpvs)
        {
            printf("%s: Clipping Vp/Vs to range [%f, %f]\n",
                      fcnm, parms.vpvs_min, parms.vpvs_max);
            // Thresholded on Vs so set Vp from that
            if (parms.lthresh_vs && !parms.lthresh_vp)
            {
                printf("%s: Clipping sets Vp from Vs\n", fcnm);
                for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
                {
                     vpvs = cvm_model->vp[indx_mod]/cvm_model->vs[indx_mod];
                     if (vpvs < parms.vpvs_min)
                     {
                         cvm_model->vp[indx_mod] = cvm_model->vs[indx_mod]
                                                  *parms.vpvs_min;
                     }
                     if (vpvs > parms.vpvs_max)
                     {
                         cvm_model->vp[indx_mod] = cvm_model->vs[indx_mod]
                                                  *parms.vpvs_max;
                     }
                }
            }
            else
            {
                printf("%s: Clipping sets Vs from Vp\n", fcnm);
                for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++){
                    vpvs = cvm_model->vp[indx_mod]/cvm_model->vs[indx_mod];
                    if (vpvs < parms.vpvs_min)
                    {
                        cvm_model->vs[indx_mod] =  cvm_model->vp[indx_mod]
                                                  /parms.vpvs_min;
                    }
                    if (vpvs > parms.vpvs_max)
                    {
                        cvm_model->vs[indx_mod] =  cvm_model->vp[indx_mod]
                                                  /parms.vpvs_max;
                    }
                }
            }
        }
    }
    // Compute the density
    printf("%s: Computing density...\n", fcnm);
    for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
    {
        cvm_model->dens[indx_mod] = density_Brocher(cvm_model->vp[indx_mod]);
    }
    if (parms.lthresh_dens)
    {
       printf("%s: Thresholding density...\n", fcnm);
       #pragma omp simd
       for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
       {
           cvm_model->dens[indx_mod] = fmax(cvm_model->dens[indx_mod],
                                            parms.dens_min);
           cvm_model->dens[indx_mod] = fmin(cvm_model->dens[indx_mod],
                                            parms.dens_max);
       }
    }
    // Compute the quality factor
    printf("%s: Computing the quality factor...\n", fcnm);
    for (indx_mod=0; indx_mod<cvm_model->npts; indx_mod++)
    {
        qualityFactor_Frankel(cvm_model->vs[indx_mod],
                              &cvm_model->Qp[indx_mod],
                              &cvm_model->Qs[indx_mod]); 
    }
    // Free space
    free(model_vp);
    free(model_vs);
    model_vp = NULL;
    model_vs = NULL;
    return 0;
}
//============================================================================//
/*!
 * @brief Writes the layer CVM to disk.  This is for data preservation
 *        reasons.
 *
 * @param[in] dirnm       name of directory to write H5 model
 * @param[in] projnm      name of project
 * @param[in] nlay        number of layers in CVM (3)
 * @param[in] cvm_model   holds the full CVM model for all layers [nlay]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int cvmio_write__h5(char *dirnm, char *projnm,
                    int nlay, struct cvm_model_struct *cvm_model)
{
    const char *fcnm = "cvmio_write__h5\0";
    char arcfile[PATH_MAX], projnm_work[512], layer[64], cwork[128];
    double *xlocs, *ylocs, *zlocs, *vp, *vs, *Qp, *Qs, *dens,
           dx, dy, dz, x0, y0, z0;
    float *xlocs4, *ylocs4, *zlocs4, *vp4, *vs4, *Qp4, *Qs4, *dens4;
    int *ien, nelemSave[3], nnpgSave[3],
        ierr, ilay, inpg, ncon, nelem, nnpg, nx, ny, nz;
    hid_t file_id, group_id;
    herr_t status;
    bool lflip = false; // Try to pack file as it was read
    // Set the HDF5 archive name
    memset(projnm_work, 0, sizeof(projnm_work));
    strcpy(projnm_work, projnm);
    strcat(projnm_work, "_cvmloc");
    meshio_write__setFilename(dirnm, projnm_work, ".h5\0", arcfile);
    // Open the h5 file 
    file_id = H5Fcreate(arcfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
 
    for (ilay=0; ilay<nlay; ilay++){
        printf("%s: Archiving layer: %d\n", fcnm, ilay+1);
        // Create the group
        memset(layer, 0, sizeof(layer));
        sprintf(layer, "Layer_%d", ilay+1);
        group_id = H5Gcreate2(file_id, layer, H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT);
        // Extract layer descriptor
        nx = cvm_model[ilay].nx;
        ny = cvm_model[ilay].ny;
        nz = cvm_model[ilay].nz;
        dx = cvm_model[ilay].dx;
        dy = cvm_model[ilay].dy;
        dz = cvm_model[ilay].dz;
        x0 = cvm_model[ilay].x0;
        y0 = cvm_model[ilay].y0;
        z0 = cvm_model[ilay].z0;
        printf("%s: Layer statistics: \n", fcnm);
        printf("                 Grid nx=%d, ny=%d, nz=%d\n", nx, ny, nz);
        printf("                 Grid spacing dx=%f, dy=%f, dz=%f\n", dx, dy, dz);
        printf("                 Origin x0=%f, y0=%f, z0=%f\n", x0, y0, z0);
        // Compute the number of elements
        ierr = regmesh_getNumberOfElements(nx, ny, nz, &nelem);
        if (ierr != 0){
            printf("%s: Error getting number of elements\n", fcnm);
            return -1;
        }
        nelemSave[ilay] = nelem;
        // Allocate space for connectivity
        ncon = 8*nelem;
        ien = (int *)calloc((size_t) ncon, sizeof(int));
        ierr = regmesh_makeHexIEN(nx, ny, nz, ien);
        if (ierr != 0){
            printf("%s: Error making hex mesh!\n", fcnm);
            return -1;
        }
        // Write the connectivity
        memset(cwork, 0, sizeof(cwork)); 
        sprintf(cwork, "%s/Connectivity", layer);
        h5_write_array__int(cwork, file_id, ncon, ien);
        free(ien);
        ien = NULL;
        // Generate the (x, y, z) locations 
        ierr = regmesh_getNumberOfAnchorNodes(nx, ny, nz, &nnpg);
        nnpgSave[ilay] = nnpg;
        xlocs = (double *)calloc((size_t) nnpg, sizeof(double));
        ylocs = (double *)calloc((size_t) nnpg, sizeof(double));
        zlocs = (double *)calloc((size_t) nnpg, sizeof(double));
        /*
        __regmesh_makeRegularNodes(nx, ny, nz,
                                   dx, dy, dz,
                                   x0, y0, z0,
                                   xlocs, ylocs, zlocs);
        */
        for (inpg=0; inpg<nnpg; inpg++){
            xlocs[inpg] = cvm_model[ilay].xlocs[inpg];
            ylocs[inpg] = cvm_model[ilay].ylocs[inpg];
            zlocs[inpg] = cvm_model[ilay].zlocs[inpg];
        }
        xlocs4 = (float *)calloc((size_t) nnpg, sizeof(float));
        ylocs4 = (float *)calloc((size_t) nnpg, sizeof(float));
        zlocs4 = (float *)calloc((size_t) nnpg, sizeof(float));
        for (inpg=0; inpg<nnpg; inpg++){
            xlocs4[inpg] = (float) xlocs[inpg];
            ylocs4[inpg] = (float) ylocs[inpg];
            zlocs4[inpg] = (float) zlocs[inpg];
        }
        // Write the anchor nodes
        memset(cwork, 0, sizeof(cwork)); 
        sprintf(cwork, "%s/XCoordinates", layer);
        h5_write_array__float(cwork, file_id, nnpg, xlocs4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/YCoordinates", layer);
        h5_write_array__float(cwork, file_id, nnpg, ylocs4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/ZCoordinates", layer);
        h5_write_array__float(cwork, file_id, nnpg, zlocs4);
        
        free(xlocs);
        free(ylocs);
        free(zlocs);
        free(xlocs4);
        free(ylocs4);
        free(zlocs4);
        xlocs = NULL;
        ylocs = NULL;
        zlocs = NULL;
        xlocs4 = NULL;
        ylocs4 = NULL;
        zlocs4 = NULL;

        // Set the material properties
        vp = (double *)calloc((size_t) nnpg, sizeof(double));
        vs = (double *)calloc((size_t) nnpg, sizeof(double));
        dens = (double *)calloc((size_t) nnpg, sizeof(double));
        Qp = (double *)calloc((size_t) nnpg, sizeof(double));
        Qs = (double *)calloc((size_t) nnpg, sizeof(double));
        __regmesh_copyRegularModel(lflip, nx, ny, nz,
                                   cvm_model[ilay].vp,
                                   cvm_model[ilay].vs,
                                   cvm_model[ilay].dens,
                                   cvm_model[ilay].Qp,
                                   cvm_model[ilay].Qs,
                                   vp, vs, dens, Qp, Qs);
        vp4 = (float *)calloc((size_t) nnpg, sizeof(float));
        vs4 = (float *)calloc((size_t) nnpg, sizeof(float));
        dens4 = (float *)calloc((size_t) nnpg, sizeof(float));
        Qp4 = (float *)calloc((size_t) nnpg, sizeof(float));
        Qs4 = (float *)calloc((size_t) nnpg, sizeof(float));
        for (inpg=0; inpg<nnpg; inpg++){
            vp4[inpg]   = (float) vp[inpg];
            vs4[inpg]   = (float) vs[inpg];
            dens4[inpg] = (float) dens[inpg];
            Qp4[inpg] = (float) Qp[inpg];
            Qs4[inpg] = (float) Qs[inpg];
        }

        memset(cwork, 0, sizeof(cwork)); 
        sprintf(cwork, "%s/CompressionalVelocity", layer);
        h5_write_array__float(cwork, file_id, nnpg, vp4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/ShearVelocity", layer);
        h5_write_array__float(cwork, file_id, nnpg, vs4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/Density", layer);
        h5_write_array__float(cwork, file_id, nnpg, dens4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/Qp", layer);
        h5_write_array__float(cwork, file_id, nnpg, Qp4);

        memset(cwork, 0, sizeof(cwork));
        sprintf(cwork, "%s/Qs", layer);
        h5_write_array__float(cwork, file_id, nnpg, Qs4);

        free(vp);
        free(vs);
        free(dens);
        free(Qp);
        free(Qs);
        free(vp4);
        free(vs4);
        free(dens4);
        free(Qp4);
        free(Qs4);
        vp = NULL;
        vs = NULL;
        dens = NULL;
        Qp = NULL;
        Qs = NULL;
        vp4 = NULL;
        vs4 = NULL;
        dens4 = NULL;
        Qp4 = NULL;
        Qs4 = NULL;
        // Close the group
        status = H5Gclose(group_id);
    }
    printf("%s: Closing file...\n", fcnm);
    status = H5Fclose(file_id);
    if (status == 0){
        printf("%s: Writing Xdmf file...\n", fcnm);
        cvmio_write__xdmf(dirnm, projnm, nlay, nelemSave, nnpgSave);
    }
    return status;
}
//============================================================================//
/*!
 * @brief Writes the Xdmf file corresponding to the CVM archive
 */
void cvmio_write__xdmf(char *dirnm, char *projnm,
                       int nlay, int *nelem, int *nnpg)
{
    FILE *xfl;
    char h5fl[PATH_MAX], xdfl[PATH_MAX], projnm_work[512], layer[64];
    int ilay;
    // Set the HDF5 archive name
    memset(projnm_work, 0, sizeof(projnm_work));
    strcpy(projnm_work, projnm);
    strcat(projnm_work, "_cvmloc");
    meshio_write__setFilename(dirnm, projnm_work, ".xdmf\0", xdfl);
    meshio_write__setFilename("./\0", projnm_work, ".h5\0", h5fl);
    // Open the file for writing
    xfl = fopen(xdfl, "w");
    fprintf(xfl, "<?xml version=\"1.0\" ?>\n");
    fprintf(xfl, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xfl, "<Xdmf Version=\"2.0\">\n");
    fprintf(xfl, "  <Domain>\n");
    for (ilay=0; ilay<nlay; ilay++){
        memset(layer, 0, sizeof(layer));
        sprintf(layer, "Layer_%d", ilay+1);
        fprintf(xfl, "    <Grid Collection=\"whole\" Name=\"CVM Archive Layer %d\">\n",
                ilay+1);
        fprintf(xfl, "      <Topology TopologyType=\"Hexahedron\" NumberOfElements=\"%d\">\n",
                nelem[ilay]);
        fprintf(xfl, "         <DataItem Dimensions=\"%d %d\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
                nelem[ilay], 8);
        fprintf(xfl, "           %s:/%s/Connectivity\n", h5fl, layer);
        fprintf(xfl, "         </DataItem>\n");
        fprintf(xfl, "      </Topology>\n");
        fprintf(xfl, "      <Geometry GeometryType=\"X_Y_Z\">\n");
        fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/XCoordinates\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/YCoordinates\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "        <DataItem ItemType=\"Function\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" Function=\"-1*$0\"> \n",
                nnpg[ilay]);
        fprintf(xfl, "          <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "             %s:/%s/ZCoordinates\n", h5fl, layer);
        fprintf(xfl, "          </DataItem>\n");
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Geometry>\n");
        fprintf(xfl, "      <Attribute Name=\"Compressional Velocity (m/s)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/CompressionalVelocity\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "      <Attribute Name=\"Shear Velocity (m/s)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/ShearVelocity\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "      <Attribute Name=\"Density (kg/m**3)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/Density\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "      <Attribute Name=\"Qp\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/Qp\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "      <Attribute Name=\"Qs\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "           %s:/%s/Qs\n", h5fl, layer);
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "      <Attribute Name=\"VpVs Ratio\" AttributeType=\"Scalar\" Center=\"Node\">\n");
        fprintf(xfl, "        <DataItem ItemType=\"Function\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Function=\"$0 / $1\">\n",
                nnpg[ilay]);

        fprintf(xfl, "          <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl, "             %s:/%s/CompressionalVelocity\n", h5fl, layer);
        fprintf(xfl, "          </DataItem>\n");
        fprintf(xfl, "          <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
                nnpg[ilay]);
        fprintf(xfl,   "           %s:/%s/ShearVelocity\n", h5fl, layer);
        fprintf(xfl, "          </DataItem>\n");
        fprintf(xfl, "        </DataItem>\n");
        fprintf(xfl, "      </Attribute>\n");
        fprintf(xfl, "    </Grid>\n");
    }
    fprintf(xfl, "  </Domain>\n");
    fprintf(xfl, "</Xdmf>\n");
    fclose(xfl);
    return;
} 
//============================================================================//
/*!
 * @brief This is a routine for quality controlling Vp  floats.  The values were
 *        lifted from Art Frankel's Fortran routine.  If the input value
 *        is NAN_CVM then we hardwire to a slow velocity.  If the number is 
 *        less than 1 or greater than VMAX_CVM we set the value to the 
 *        previous grid point. Otherwise, we return valin
 *
 * @param[in] valin     Vp value to quality control
 *
 * @result the quality controlled p velocity
 * 
 */
float __cvm_io_qcFloat_vp(float valin)
{
     float val = valin;
     // Fix NaN's
     if ((double) val >= NAN_CVM - 1.e-4){val = VP_DEFAULT;} //2300.0;}
     if ((double) val < VMIN_CVM || (double) val > VMAX_CVM)
     {
         val = VP_DEFAULT;
     }
     return val;
}
//============================================================================//
/*!
 * @brief This is a routine for quality controlling Vs floats.  The values were
 *        lifted from Art Frankel's Fortran routine.  If the input value
 *        is NAN_CVM then we hardwire to a slow velocity.  If the number is 
 *        less than 1 or greater than VMAX_CVM we set the value to the 
 *        previous grid point. Otherwise, we return valin
 *
 * @param[in] valin     Vs value to quality control
 *
 * @result the quality controlled s velocity 
 * 
 */
float __cvm_io_qcFloat_vs(float valin) //, float val0)
{
     float val = valin;
     // Fix NaN's
     if ((double) val >= NAN_CVM - 1.e-4){val = VS_DEFAULT;} //1333.0;}
     if ((double) val < VMIN_CVM || (double) val > VMAX_CVM)
     {
         val = VS_DEFAULT;
     }
     return val;
}
//============================================================================//
/*!
 * Reads a binary file into memory 
 * 
 * @param[in] filename     name of binary file to read into memory
 *
 * @param[out] nbytes      file size in bytes bytes; 0 if error
 * 
 * @result binary file read as a char* of sizes nbytes for unpacking 
 *         by something like utilsFilesUnpackFloat
 *
 * @author Ben Baker, ISTI
 *
 */
char *cvmio_readBinaryFile(char *filename, int *nbytes)
{
    const char *fcnm = "cvmio_readBinaryFile\0";
    FILE *bfile;
    char *buffer;
    int nread;
    //------------------------------------------------------------------------//
    //
    // Get size
    *nbytes = 0;
    *nbytes = cvmio_getBinaryFileSize(filename);
    if (*nbytes < 1)
    {
        printf("%s: Error reading binary file\n",fcnm);
        return NULL;
    }
    // Set space
    buffer = (char *)calloc((size_t) (*nbytes), sizeof(char));
    if (buffer == NULL)
    {
        *nbytes =-1;
        printf("%s: Error loading file\n", fcnm);
        return NULL;
    }
    // Load file 
    bfile = fopen(filename, "rb");
    nread = (int) (fread(buffer,sizeof(buffer), (size_t) *nbytes, bfile));
    if (nread < 1){printf("%s: Possible no items were read\n", fcnm);}
    fclose(bfile);
    return buffer;
}
//============================================================================//
/*! 
 * @brief Returns the size of a binary file in bytes
 *
 * @param[in]  filename    name of file of which to determine size
 * 
 * @result size of file in bytes
 *
 * @author Ben Baker, ISTI
 *
 */
int cvmio_getBinaryFileSize(char *filename)
{
    const char *fcnm = "cvmio_getBinaryFileSize\0";
    FILE *f;
    int size;
    if (!os_path_isfile(filename)){
        printf("%s: File %s does not exist!\n", fcnm, filename);
        return 0;
    }
    f = fopen(filename, "rb");
    if (f == NULL) return -1;
    fseek(f, 0, SEEK_END);
    size = (int) (ftell(f));
    fclose(f);
    return size;
}
