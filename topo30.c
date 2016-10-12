#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <hdf5.h>
#include <gsl/gsl_interp.h>
#include "cvm.h"
#include "mesh.h"
#include "h5_cinter.h"
#include "iscl/array/array.h"
#include "iscl/interpolate/interpolate.h"
#include "iscl/log/log.h"
#include "iscl/os/os.h"
#include "iscl/sorting/sorting.h"

struct topoGrid_struct
{
    double *topo;    /*!< Eleveation above sea-level (m) at each (lat,lon).
                          The i'th longitude j'th latitude are accessed by
                          j*nlon + i [ntopo] */
    double *lats;    /*!< Latitudes (degrees, [-90,90]) stored in increasing
                          order [nlat] */
    double *lons;    /*!< Longitudes (degrees, [0,360]) stored in increasing
                          order [nlon] */
    int nlat;        /*!< Number of latitudes */
    int nlon;        /*!< Number of longitudes */
    int ntopo;       /*!< Number of points in topo */
};

static double __tform(double a, double b,
                      double c, double d,
                      double x, int *ierr);

/*!
 * @brief Reads the topography contained in the box define by a southwest
 *        corner [xlat0in, xlon0in] and northeast corner [xlat1in, xlon1in] 
 *
 * @param[in] fname       netCDF topo30 file
 * @param[in] xlat0in     latitude of southwest corner of box
 *                        in degrees [-90,90]
 * @param[in] xlon0in     longitude of southwest corner of box
 *                        in degrees [-180,360]
 * @param[in] xlat1in     latitude of northwest corner of box
 *                        in degrees [-90,90]
 * @param[in] xlon1in     longitude of northwest corner of box
 *                        in degrees [-180,360]
 *
 * @param[out]            topography structure
 *
 * @result 0 indicates success
 * 
 */
int topo30_read__netCDF(const char *fname,
                        double xlat0in, double xlon0in,
                        double xlat1in, double xlon1in,
                        struct topoGrid_struct *topo)
{
    const char *fcnm = "topo30_read__netCDF\0";
    double *lats, *lons, xlat0, xlon0, xlat1, xlon1; 
    int ierr, ilat, ilat0, ilat1, ilon, ilon0, ilon1, indx,
        nlat, nlat_read, nlon, nlon_read;
    int16_t *topo_data;
    hid_t dataset_id, dataspace_id, file_id, memspace_id;
    hsize_t block[2], count[2], count_out[2], dimsm[2],
            stride[2], offset[2], offset_out[2];
    herr_t status;
    //------------------------------------------------------------------------//
    //
    // Initialize
    ierr = 0;
    status = 0;
    topo_data = NULL;
    lats = NULL;
    lons = NULL;
    memset(topo, 0, sizeof(struct topoGrid_struct));
    // Require the file exists
    if (!os_path_isfile(fname))
    {
        log_errorF("%s: Error file %s does not exist\n", fcnm, fname);
        return -1;
    }
    // Copy inputs and switch longitude to [0,360] incase user forgot 
    xlon0 = xlon0in;
    xlat0 = xlat0in;
    xlon1 = xlon1in;
    xlat1 = xlat1in;
    if (xlon0 < 0.0){xlon0 = xlon0 + 360.0;}
    if (xlon1 < 0.0){xlon1 = xlon1 + 360.0;}
    // Check the inputs
    if (xlon0 < 0.0 || xlon0 > 360.0)
    {
        log_errorF("%s: xlon0 invalid %f\n", fcnm, xlon0);
        return -1;
    }
    if (xlat0 < -90.0 || xlat0 > 90.0)
    {
        log_errorF("%s: xlat0 invalid %f\n", fcnm, xlat0); 
        return -1;
    }
    if (xlon1 < 0.0 || xlon1 > 360.0)
    {   
        log_errorF("%s: xlon1 invalid %f\n", fcnm, xlon1);
        return -1;
    }   
    if (xlat1 < -90.0 || xlat1 > 90.0)
    {   
        log_errorF("%s: xlat1 invalid %f\n", fcnm, xlat1); 
        return -1;
    }
    if (xlat0 > xlat1)
    {
        log_errorF("%s: Error xlat0 > xlat1 %f %f\n", fcnm, xlat0, xlat1);
        return -1;
    }
    if (xlon0 > xlon1)
    {
        log_errorF("%s: Error xlon0 > xlon1 %f %f\n", fcnm, xlon0, xlon1);
        return -1;
    }
    // Open the netCDF file for reading
    file_id = h5_open_rdonly(fname);
    // Read the sizes of the
    nlon = h5_get_array_size("/lon\0", file_id);
    nlat = h5_get_array_size("/lat\0", file_id);
    if (nlat < 1 || nlon < 1)
    {
        if (nlat < 1){log_errorF("%s: Error no lat points!\n", fcnm);}
        if (nlon < 1){log_errorF("%s: Error no lon points!\n", fcnm);}
        H5Fclose(file_id);
        return -1;
    }
    // Read the latitudes / longitudes
    lats = (double *)calloc(nlat, sizeof(double));
    lons = (double *)calloc(nlon, sizeof(double));
    ierr = h5_read_array__double("lat\0", file_id, nlat, lats);
    if (ierr != 0)
    {
        log_errorF("%s: Error reading latitude\n", fcnm);
        goto ERROR;
    }
    ierr = h5_read_array__double("lon\0", file_id, nlon, lons);
    if (ierr != 0)
    {   
        log_errorF("%s: Error reading longitudes\n", fcnm);
        goto ERROR;
    }
    // Locate the indices bounding this box: xv[indx] <= x < xv[indx+1]
    ilat0 = gsl_interp_bsearch(lats, xlat0, 0, nlat);
    ilat1 = gsl_interp_bsearch(lats, xlat1, 0, nlat);
    ilon0 = gsl_interp_bsearch(lons, xlon0, 0, nlon);
    ilon1 = gsl_interp_bsearch(lons, xlon1, 0, nlon);
    if (xlat1 > lats[ilat1] && ilat1 < nlat - 1){ilat1 = ilat1 + 1;}
    if (xlon1 > lons[ilon1] && ilon1 < nlon - 1){ilon1 = ilon1 + 1;}
    // Set the locations for the hyperslab
    nlat_read = ilat1 - ilat0 + 1; 
    nlon_read = ilon1 - ilon0 + 1;
    offset[0] = ilat0;
    offset[1] = ilon0;
    count[0] = nlat_read;
    count[1] = nlon_read; 
    dimsm[0] = nlat_read;
    dimsm[1] = nlon_read;
    stride[0] = 1;
    stride[1] = 1;
    block[0] = 1;
    block[1] = 1;
    topo_data = (int16_t *)calloc(nlat_read*nlon_read, sizeof(int16_t));
    // Open the dataset for reading
    dataset_id = H5Dopen2(file_id, "/z\0", H5P_DEFAULT);
    dataspace_id = H5Dget_space(dataset_id);
    // Define hyperslab in dataset
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, stride,
                                 count, NULL);
    // Define the memory space
    memspace_id = H5Screate_simple(2, dimsm, NULL);   
    // Define the hyperslab
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0]  = nlat_read;
    count_out[1]  = nlon_read;
    status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset_out, NULL,
                                 count_out, block);
    // Read data from hyperslab in the file into the hyperslab in 
    status = H5Dread(dataset_id, H5T_NATIVE_INT16, memspace_id, dataspace_id,
                     H5P_DEFAULT, topo_data);
/*
    status += H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset,
                                  stride, count, block);
    // Read it
getchar();
printf("go\n");
    status += H5Dread(dataset_id, H5T_NATIVE_INT16, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, topo_data);
printf("back\n");
*/
    if (status != 0)
    {
        log_errorF("%s: Error reading topo data\n", fcnm);
        goto ERROR;
    }
    // Close it
    status = H5Sclose(memspace_id);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    // Copy the topography
    topo->nlat = nlat_read;
    topo->nlon = nlon_read;
    topo->ntopo = topo->nlat*topo->nlon;
    topo->lats = (double *)calloc(topo->nlat, sizeof(double));
    topo->lons = (double *)calloc(topo->nlon, sizeof(double)); 
    cblas_dcopy(topo->nlat, &lats[ilat0], 1, topo->lats, 1);
    cblas_dcopy(topo->nlon, &lons[ilon0], 1, topo->lons, 1);
    topo->topo = (double *)calloc(topo->ntopo, sizeof(double));
    #pragma omp simd collapse(2)
    for (ilat=0; ilat<nlat_read; ilat++)
    {
        for (ilon=0; ilon<nlon_read; ilon++)
        {
            indx = ilat*nlon_read + ilon;
            topo->topo[indx] = (double) topo_data[indx];
        }
    }
ERROR:;
    ierr = h5_close(file_id);
    if (topo_data != NULL){free(topo_data);}
    return ierr;
}
//============================================================================//
/*!
 * @brief Deforms a mesh to the topography in the topo30 file
 *
 * @author Ben Baker (ISTI)
 * 
 * @bugs the bisection hunt doesn't work with an irregular mesh and will
 *       usually default to a linear search
 *
 */
int topo30_deformMesh(const char *topofl, int utm_zone,
                      double ztopo_min,
                      struct mesh_struct *mesh )
{
    const char *fcnm = "mesh_deform2topo30\0";
    struct topoGrid_struct topo;
    double *lat_int, *lon_int, *topo_int, *x, *xl, *xlocs, *y, *ylocs,
           xlat0, xlon0, xlat1, xlon1, xutm_max,
           xutm_min, yutm_min, yutm_max, zbase,
           zmax, zmin, zint_base;
    int *iptrx, *iperm, *mesh2zind,
        i, i1, i2, ierr, iloc, ilocx, ilocy, indx, inpg, nloc,
        nlinear_search, nlocx, npts, npsort, ns;
    bool *lmask, lsearch, lxdif, lydif;
    const double extra_ll = 0.1; // additional latitude/longitude to add
    const double tol = 1.e-2;    // centimeter accuracy
    //------------------------------------------------------------------------//
    if (mesh == NULL)
    {
        log_errorF("%s: Error mesh is NULL\n", fcnm);
        return -1;
    }
    if (mesh->xlocs == NULL || mesh->ylocs == NULL || mesh->zlocs == NULL)
    {
        if (mesh->xlocs == NULL)
        {
            log_errorF("%s: Error xlocs must be set\n", fcnm);
        }
        if (mesh->ylocs == NULL)
        {
            log_errorF("%s: Error ylocs must be set\n", fcnm);
        }
        if (mesh->zlocs == NULL)
        {
            log_errorF("%s: Error zlocs must be set\n", fcnm);
        }
        return -1;
    }
    if (mesh->nnpg < 1)
    {
        log_errorF("%s: Error no points in mesh!\n", fcnm);
        return -1;
    }
    // Get the minimum and maximum corners of the mesh 
    xutm_min = array_min__double(mesh->nnpg, mesh->xlocs);
    yutm_min = array_min__double(mesh->nnpg, mesh->ylocs);
    xutm_max = array_max__double(mesh->nnpg, mesh->xlocs);
    yutm_max = array_max__double(mesh->nnpg, mesh->ylocs);
    zmin = array_min__double(mesh->nnpg, mesh->zlocs);
    zmax = array_max__double(mesh->nnpg, mesh->zlocs);
    // Convert to lats/lons
    utm_geo_utm2ll(&xutm_min, &yutm_min, &utm_zone, &xlat0, &xlon0);
    utm_geo_utm2ll(&xutm_max, &yutm_max, &utm_zone, &xlat1, &xlon1);
    if (xlon0 < 0.0){xlon0 = xlon0 + 360.0;}
    if (xlon1 < 0.0){xlon1 = xlon1 + 360.0;}
    // Try to fix any stretching weirdness so that it doesn't effect
    // the subsequent interpolation
    xlat0 = xlat0 - extra_ll;
    xlon0 = xlon0 - extra_ll;
    xlat1 = xlat1 + extra_ll;
    xlon1 = xlon1 + extra_ll;
    // Get the topography
    ierr = topo30_read__netCDF(topofl,
                               xlat0, xlon0,
                               xlat1, xlon1,
                               &topo);
    if (ierr != 0)
    {
        log_errorF("%s: Error getting topography\n");
        return -1;
    }
    // Set a logical mask for all points above the 
    zint_base = fmax(zmin, zmax - ztopo_min*1.e3);
    zbase = DBL_MAX;
    lmask = (bool *)calloc(mesh->nnpg, sizeof(bool));
    npts = 0;
    for (inpg=0; inpg<mesh->nnpg; inpg++)
    {
        lmask[inpg] = true;
        if (mesh->zlocs[inpg] >= zint_base)
        {
           zbase = fmin(zbase, mesh->zlocs[inpg]);
           lmask[inpg] = false;
           npts = npts + 1;
        }
    }
    if (npts == 0)
    {
        log_errorF("%s: Error no points to deform\n", fcnm);
        return -1;
    }
    // Extract the points to deform
    xlocs = (double *)calloc(npts+1, sizeof(double));
    ylocs = (double *)calloc(npts+1, sizeof(double));
    i = 0;
    for (inpg=0; inpg<mesh->nnpg; inpg++)
    {
        if (!lmask[inpg])
        {
            xlocs[i] = mesh->xlocs[inpg];
            ylocs[i] = mesh->ylocs[inpg];
            i = i + 1;
        }
    }
    // Create a set of unique x, y points
    iperm = sorting_argsort__double(npts, xlocs, ASCENDING, &ierr);
    if (ierr != 0)
    {
        log_errorF("%s: Error sorting xlocs\n", fcnm);
        return -1;
    }
    ierr = 0;
    ierr += __sorting_applyPermutation__double(npts, iperm, xlocs, xlocs);
    ierr += __sorting_applyPermutation__double(npts, iperm, ylocs, ylocs);
    if (ierr != 0)
    {
        log_errorF("%s: Error permuting xlocs and/or ylocs\n", fcnm);
        return -1;
    }
    free(iperm);
    // Step through list and sort y's
    xlocs[npts] = xlocs[0] - 100.0; // dirty trick to make next section work
    ylocs[npts] = ylocs[0] - 100.0; // dirty trick to make next section work
    i1 = 0;
    for (i=0; i<npts; i++)
    {
        if (fabs(xlocs[i+1] - xlocs[i1]) > tol)
        {
            i2 = i;
            npsort = i2 - i1 + 1;
            if (npsort > 0)
            {
                ierr = __sorting_sort__double(npsort, &ylocs[i1], ASCENDING);
                if (ierr != 0)
                {
                    log_errorF("%s: Error sorting ylocs\n", fcnm);
                    return -1;
                }
            }
            i1 = i + 1;
        }
    }
    // Count the number of unique points
    nloc = 0;
    nlocx = 0;
    for (i=0; i<npts; i++)
    {
        lxdif = false;
        lydif = false;
        if (fabs(xlocs[i+1] - xlocs[i]) > tol){lxdif = true;}
        if (fabs(ylocs[i+1] - ylocs[i]) > tol){lydif = true;}
        if (lxdif || lydif)
        {
            nloc = nloc + 1;
            if (lxdif){nlocx = nlocx + 1;}
        }
    }
    if (nloc < 1 || nlocx < 1)
    {
        if (nloc < 1){log_errorF("%s: nloc is zero\n", fcnm);}
        if (nlocx < 1){log_errorF("%s: nlocx is zero\n", fcnm);}
        return -1;
    }
    // And save them
    x = (double *)calloc(nloc, sizeof(double));
    y = (double *)calloc(nloc, sizeof(double));
    xl = (double *)calloc(nlocx, sizeof(double));
    iptrx = (int *)calloc(nlocx+1, sizeof(int));
    iloc = 0;
    ilocx = 0;
    for (i=0; i<npts; i++)
    {
        lxdif = false;
        lydif = false;
        if (fabs(xlocs[i+1] - xlocs[i]) > tol){lxdif = true;}
        if (fabs(ylocs[i+1] - ylocs[i]) > tol){lydif = true;}
        if (lxdif || lydif)
        {
            x[iloc] = xlocs[i];
            y[iloc] = ylocs[i];
            iloc = iloc + 1;
            if (lxdif)
            {
                iptrx[ilocx+1] = iloc;
                xl[ilocx] = xlocs[i];
                ilocx = ilocx + 1;
            }
        }
    }
    free(xlocs);
    free(ylocs);
    // Convert the (x,y) pairs to (lon,lat)
    lat_int  = (double *)calloc(nloc, sizeof(double));
    lon_int  = (double *)calloc(nloc, sizeof(double));
    for (i=0; i<nloc; i++)
    {
        utm_geo_utm2ll(&x[i], &y[i], &utm_zone, &lat_int[i], &lon_int[i]);
        if (lon_int[i] < 0.0){lon_int[i] = lon_int[i] + 360.0;}
    }
    // Now interpolate each point
    topo_int = interpolate_interp2d_gsl(topo.nlon, topo.lons,
                                        topo.nlat, topo.lats,
                                        topo.ntopo, topo.topo,
                                        nloc, lon_int,
                                        nloc, lat_int,
                                        BILINEAR,
                                        false,
                                        &nloc, &ierr);
    if (ierr != 0)
    {
        log_errorF("%s: Error interpolating topography\n", fcnm);
        return -1;
    }
    // Map each point in the mesh to a point in the interpolated points
    nlinear_search = 0;
    mesh2zind = (int *)calloc(mesh->nnpg, sizeof(int));
    for (inpg=0; inpg<mesh->nnpg; inpg++)
    {
        mesh2zind[inpg] =-1;
        if (lmask[inpg]){continue;}
        // Hunt for the x point
        lsearch = false;
        ilocx = gsl_interp_bsearch(xl, mesh->xlocs[inpg], 0, nlocx);
        iloc = iptrx[ilocx];
        indx =-1;
        if (fabs(xl[ilocx] - mesh->xlocs[inpg]) < tol)
        {
            i1 = iptrx[ilocx];
            i2 = iptrx[ilocx+1];
            ns = i2 - i1;
            ilocy = i1 + gsl_interp_bsearch(&y[i1], mesh->ylocs[inpg], 0, ns);
            if (fabs(x[ilocy] - mesh->xlocs[inpg]) < tol &&
                fabs(y[ilocy] - mesh->ylocs[inpg]) < tol)
            {
                indx = ilocy;
            }
            else
            {
                //log_warnF("%s: Linearly search in x/y\n", fcnm);
                nlinear_search = nlinear_search + 1;
                lsearch = true;
            }
        }
        else // Linearly search for it
        {
            //log_warnF("%s: Linearly searching in x\n", fcnm);
            nlinear_search = nlinear_search + 1;
            lsearch = true;
        } 
        // My clever binary search scheme failed - look the hard way 
        if (lsearch)
        {
            for (i=0; i<nloc; i++)
            {
                if (fabs(x[i] - mesh->xlocs[inpg]) < tol && 
                    fabs(y[i] - mesh->ylocs[inpg]) < tol)
                {
                    indx = i;
                    break;
                }
            }
            if (indx < 0)
            {
                log_errorF("%s: Couldn't find %f %f\n", fcnm,
                           mesh->xlocs[inpg], mesh->ylocs[inpg]);
                return -1;
            }
        }
        mesh2zind[inpg] = indx; 
    }
    if (nlinear_search > 0)
    {
        log_infoF("%s: Number of linear searches: %d\n", fcnm);
    }
    // Now interpolate each z value from the interface to the free surface
    for (inpg=0; inpg<mesh->nnpg; inpg++)
    {
        if (lmask[inpg]){continue;}
        if (mesh2zind[inpg] < 0)
        {
            log_errorF("%s: This is a surprising place to be\n", fcnm);
            return -1;
        }
        indx = mesh2zind[inpg];
        mesh->zlocs[inpg] = __tform(zbase, zmax,
                                    zbase, zmax + topo_int[indx],
                                    mesh->zlocs[inpg], &ierr);
/*
        mesh->zlocs[inpg] = __tform(zint_base, zmax_deformed,
                                    zbase, zmax,
                                    zmax + topo_int[indx], &ierr);
*/
//printf("%f\n", mesh->zlocs[inpg]);
//getchar();
    } 
    //zmax = 
    free(iptrx);
    free(lmask);
    free(mesh2zind);
    free(xl);
    free(lat_int);
    free(lon_int);
    free(topo_int);
    free(x);
    free(y);
    free(topo.topo);
    free(topo.lats);
    free(topo.lons);
    return 0;
}

static double __tform(double a, double b,
                      double c, double d,
                      double x, int *ierr)
{
    const char *fcnm = "__tform\0";
    double c1, c2, det, xi;
    *ierr = 0;
    if (a == b)
    {
        log_errorF("%s: Determinant undefind\n", fcnm);
        *ierr = 1;
        return 0.0;
    }
    det = 1.0/(b - a);
    c1 = det*(b*c - a*d);
    c2 = det*(d - c);
    xi = c1 + x*c2;
    return xi;
}
