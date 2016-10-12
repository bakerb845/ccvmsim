#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <cblas.h>
#include <math.h>
#include <iscl/log/log.h>
#include "cvm.h"
#include "iscl/interpolate/interpolate.h"
/*!
 * @brief Extends the CVM so that there are no `spaces' in between layers
 *
 * @param[in] nlay     number of layers in model
 *
 * @param[in,out] cvm_model  on input holds the cvm_model read from disk
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
    for (ilay=0; ilay<nlay-1; ilay++)
    {
        // Get the top
        nx_top = cvm_model[ilay].nx;
        ny_top = cvm_model[ilay].ny;
        nxy_top = nx_top*ny_top;
        nxyz_top = cvm_model[ilay].nz*nxy_top;
        cvm_model[ilay].npts = cvm_model[ilay].npts + nxy_top;
        cvm_model[ilay].nz = cvm_model[ilay].nz + 1; // add layer
        cvm_model[ilay].z1 = cvm_model[ilay+1].z0;
        xtop = (double *)calloc((size_t) nx_top, sizeof(double));
        ytop = (double *)calloc((size_t) ny_top, sizeof(double));
        for (ix=0; ix<nx_top; ix++)
        {
            xtop[ix] = cvm_model[ilay].x0 + cvm_model[ilay].dx*(double) ix;
        }
        for (iy=0; iy<ny_top; iy++)
        {
            ytop[iy] = cvm_model[ilay].y0 + cvm_model[ilay].dy*(double) iy;
        }
        // Get the bottom
        nx_bot = cvm_model[ilay+1].nx;
        ny_bot = cvm_model[ilay+1].ny;
        nxy_bot = nx_bot*ny_bot;
        xbot = (double *)calloc((size_t) nx_bot, sizeof(double));
        ybot = (double *)calloc((size_t) ny_bot, sizeof(double));
        for (ix=0; ix<nx_bot; ix++)
        {
            xbot[ix] = cvm_model[ilay+1].x0 + cvm_model[ilay+1].dx*(double) ix;
        }
        for (iy=0; iy<ny_bot; iy++)
        {
            ybot[iy] = cvm_model[ilay+1].y0 + cvm_model[ilay+1].dy*(double) iy;
        }
        // Interpolate vp, vs, density
        ref = (double *)calloc((size_t) nxyz_top, sizeof(double));
        // Resize vp
        memcpy(ref, cvm_model[ilay].vp, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].vp);
        cvm_model[ilay].vp = (double *)
                             calloc((size_t) cvm_model[ilay].npts,
                                    sizeof(double));
        memcpy(cvm_model[ilay].vp, ref, (size_t) nxyz_top*sizeof(double));
        // Resize vs
        memcpy(ref, cvm_model[ilay].vs, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].vs);
        cvm_model[ilay].vs = (double *)
                             calloc((size_t) cvm_model[ilay].npts,
                                    sizeof(double));
        memcpy(cvm_model[ilay].vs, ref, (size_t) nxyz_top*sizeof(double));
        // Resize density 
        memcpy(ref, cvm_model[ilay].dens, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].dens);
        cvm_model[ilay].dens = (double *)
                               calloc((size_t) cvm_model[ilay].npts,
                                      sizeof(double));
        memcpy(cvm_model[ilay].dens, ref, (size_t) nxyz_top*sizeof(double));
        // Resize Qp 
        memcpy(ref, cvm_model[ilay].Qp, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].Qp);
        cvm_model[ilay].Qp = (double *)
                             calloc((size_t) cvm_model[ilay].npts,
                                    sizeof(double));
        memcpy(cvm_model[ilay].Qp, ref, (size_t) nxyz_top*sizeof(double));
        // Resize Qs
        memcpy(ref, cvm_model[ilay].Qs, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].Qs);
        cvm_model[ilay].Qs = (double *)
                             calloc((size_t) cvm_model[ilay].npts,
                                    sizeof(double));
        memcpy(cvm_model[ilay].Qs, ref, (size_t) nxyz_top*sizeof(double));
        // Resize x nodes
        memcpy(ref, cvm_model[ilay].xlocs, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].xlocs);
        cvm_model[ilay].xlocs = (double *)
                                calloc((size_t) cvm_model[ilay].npts,
                                       sizeof(double));
        memcpy(cvm_model[ilay].xlocs, ref, (size_t) nxyz_top*sizeof(double));
        // Resize y nodes
        memcpy(ref, cvm_model[ilay].ylocs, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].ylocs);
        cvm_model[ilay].ylocs = (double *)
                                calloc((size_t) cvm_model[ilay].npts,
                                       sizeof(double));
        memcpy(cvm_model[ilay].ylocs, ref, (size_t) nxyz_top*sizeof(double));
        // Resize z nodes
        memcpy(ref, cvm_model[ilay].zlocs, (size_t) nxyz_top*sizeof(double));
        free(cvm_model[ilay].zlocs);
        cvm_model[ilay].zlocs = (double *)
                                calloc((size_t) cvm_model[ilay].npts,
                                       sizeof(double));
        memcpy(cvm_model[ilay].zlocs, ref, (size_t) nxyz_top*sizeof(double));
        free(ref);
        // Interpolate vp
        printf("%s: Interpolating vp in layer %d...\n", fcnm, ilay+1);
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].vp[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].vp[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0)
        {
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate vs
        printf("%s: Interpolating vs in layer %d...\n", fcnm, ilay+1);
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].vs[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].vs[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0)
        {
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate density
        printf("%s: Interpolating density in layer %d...\n", fcnm, ilay+1);
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].dens[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].dens[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0)
        {
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate Qp
        printf("%s: Interpolating Qp in layer %d...\n", fcnm, ilay+1);
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].Qp[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].Qp[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0)
        {
            log_errorF("%s: Error interpolating Qp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Interpolate Qs
        printf("%s: Interpolating Qs in layer %d...\n", fcnm, ilay+1);
        ierr = __interpolate_interp2d(nx_bot, xbot,
                                      ny_bot, ybot,
                                      nxy_bot, &cvm_model[ilay+1].Qs[0],
                                      nx_top, xtop,
                                      ny_top, ytop,
                                      nxy_top, &cvm_model[ilay].Qs[nxyz_top],
                                      BILINEAR, true);
        if (ierr != 0)
        {
            log_errorF("%s: Error interpolating vp in layer: %d\n",
                       fcnm, ilay+1);
        }
        // Fix nodes
        indx = 0;
        jndx = nxyz_top;
        for (iy=0; iy<ny_top; iy++)
        {
            for (ix=0; ix<nx_top; ix++)
            {
                indx = iy*nx_top + ix;
                jndx = nxyz_top + iy*nx_top + ix;
                cvm_model[ilay].xlocs[jndx] = cvm_model[ilay].xlocs[indx];
                cvm_model[ilay].ylocs[jndx] = cvm_model[ilay].ylocs[indx];
                cvm_model[ilay].zlocs[jndx] = cvm_model[ilay+1].z0; 
                //indx = indx + 1;
                //jndx = jndx + 1;
            }
        }
        // Free space
NEXT:;
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
