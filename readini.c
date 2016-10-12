#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wpadded"
#endif
#include <iniparser.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#include <iscl/os/os.h>
#include <iscl/geodetic/geodetic.h>
#include "cvm.h"
#include "cvm_constants.h"

/*!
 * @brief Reads the CVM mesher initialization file
 *
 * @param[in] projnm       project name
 *
 * @param[out] parms       CVM mesher parameters
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int cvm_readini(char *projnm, struct cvm_parms_struct *parms)
{
    const char *fcnm = "cvm_readini\0";
    char ini_file[PATH_MAX], scoarsen[64];
    dictionary *ini;
    const char *s;
    double lat1, lon1, utmx, utmy, xoff_max, yoff_max;
    int ic, ictype, k;
    int iway = 1; // UTM -> (lat,lon)
    bool suppress = false; // don't suprress the UTM computation
    //------------------------------------------------------------------------//
    //
    // Null out then set constants 
    memset(parms, 0, sizeof(struct cvm_parms_struct));
    strcpy(ini_file, projnm);
    strcat(ini_file, ".ini\0");
    /* TODO make an CVM model init function */
    parms->nlay_cvm = nlay_cvm;
    parms->utm_x0_cvm = utm_x0_cvm;
    parms->utm_y0_cvm = utm_y0_cvm;
    parms->utmzone_cvm = utmzone_cvm;
    parms->lat0_cvm = lat0_cvm;
    parms->lon0_cvm = lon0_cvm;
    parms->common_mult = common_mult;
    parms->dxl_cvm = (double *)calloc((size_t) parms->nlay_cvm, sizeof(double));
    parms->dyl_cvm = (double *)calloc((size_t) parms->nlay_cvm, sizeof(double));
    parms->dzl_cvm = (double *)calloc((size_t) parms->nlay_cvm, sizeof(double));
    parms->z0_cvm  = (double *)calloc((size_t) parms->nlay_cvm, sizeof(double));
    parms->nxl_cvm = (int *)calloc((size_t) parms->nlay_cvm, sizeof(int));
    parms->nyl_cvm = (int *)calloc((size_t) parms->nlay_cvm, sizeof(int));
    parms->nzl_cvm = (int *)calloc((size_t) parms->nlay_cvm, sizeof(int));
    xoff_max = (double) (nxl_cvm[0]-1)*dxl_cvm[0];
    yoff_max = (double) (nyl_cvm[0]-1)*dyl_cvm[0];
    for (k=0; k<nlay_cvm; k++)
    {
        parms->dxl_cvm[k] = dxl_cvm[k];
        parms->dyl_cvm[k] = dyl_cvm[k];
        parms->dzl_cvm[k] = dzl_cvm[k];
        parms->z0_cvm[k]  = z0_cvm[k];
        parms->nxl_cvm[k] = nxl_cvm[k];
        parms->nyl_cvm[k] = nyl_cvm[k];
        parms->nzl_cvm[k] = nzl_cvm[k];
        // min prevents an error with interpolation out of grid since some 
        // layers go further than others
        xoff_max = fmin(xoff_max, (double) (nxl_cvm[k]-1)*dxl_cvm[k] );
        yoff_max = fmin(yoff_max, (double) (nyl_cvm[k]-1)*dyl_cvm[k] );
    }
    // Compute usable upper right corner w.r.t. interpolation
    utmx = utm_x0_cvm + xoff_max;
    utmy = utm_y0_cvm + yoff_max;
    utm_geo(&parms->lon1_cvm, &parms->lat1_cvm, &utmx, &utmy,
            &parms->utmzone_cvm, &iway, &suppress); 
    if (parms->lon1_cvm < 0.0){parms->lon1_cvm = parms->lon1_cvm + 360.0;}
    // Compute the min/max lat/lons of the CVM by computing four corners of 
    // Compute maximum lat and minimum lon
    utmx = utm_x0_cvm + 0.0;
    utmy = utm_y0_cvm + yoff_max;
    utm_geo(&lon1, &lat1,  &utmx, &utmy, &parms->utmzone_cvm, &iway, &suppress);
    if (lon1 < 0.0){lon1 = lon1 + 360.0;}
    parms->lat1_cvm = fmin(parms->lat1_cvm, lat1);
    parms->lon0_cvm = fmax(parms->lon0_cvm, lon1); 
    // Compute minimum lat and maximum longitude
    utmx = utm_x0_cvm + xoff_max;
    utmy = utm_y0_cvm;
    utm_geo(&lon1, &lat1,  &utmx, &utmy, &parms->utmzone_cvm, &iway, &suppress);
    if (lon1 < 0.0){lon1 = lon1 + 360.0;}
    parms->lat0_cvm = fmax(parms->lat0_cvm, lat1);
    parms->lon1_cvm = fmax(parms->lon1_cvm, lon1);
    // have ini parser load the ini file
    ini = iniparser_load(ini_file);
    if (ini == NULL)
    {
        printf("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    s = iniparser_getstring(ini, "cvm_h5repack:cvm_moddir\0", "./\0");
    strcpy(parms->cvm_moddir, s);
    if (!os_path_isdir(parms->cvm_moddir))
    {
        printf("%s: CVM directory doesn't exist\n", fcnm);
        return -1;
    }
    // get the mesh directory
    strcpy(parms->cvm_outputdir, "./\0");
    s = iniparser_getstring(ini, "cvm_h5repack:cvm_outputdir", "./\0");
    if (s != NULL)
    {
        if (strlen(s) > 0)
        {
            memset(parms->cvm_outputdir, 0, sizeof(parms->cvm_outputdir));
            strcpy(parms->cvm_outputdir, s); 
        }
    }   
    if (os_makedirs(parms->cvm_outputdir) != 0)
    {
        printf("%s: Failed to make cvm output directory\n", fcnm);
        return -1; 
    }
    // get the nll directory
    strcpy(parms->nll_outputdir, "./\0");
    s = iniparser_getstring(ini, "cvm_h5repack:nll_outputdir", "./\0");
    if (s != NULL)
    {
        if (strlen(s) > 0)
        {
            memset(parms->nll_outputdir, 0, sizeof(parms->nll_outputdir));
            strcpy(parms->nll_outputdir, s); 
        }
    }
    if (os_makedirs(parms->nll_outputdir) != 0)
    {
        printf("%s: Failed to make nll output directory\n", fcnm);
        return -1; 
    }
    // get lat/lon/max depth
    parms->lat0 = iniparser_getdouble(ini, "cvm_h5repack:lat0\0",
                                      parms->lat0_cvm);
    parms->lon0 = iniparser_getdouble(ini, "cvm_h5repack:lon0\0",
                                      parms->lon0_cvm);
    parms->lat1 = iniparser_getdouble(ini, "cvm_h5repack:lat1\0",
                                      parms->lat1_cvm);
    parms->lon1 = iniparser_getdouble(ini, "cvm_h5repack:lon1\0",
                                      parms->lon1_cvm);
    parms->zmin = iniparser_getdouble(ini, "cvm_h5repack:zmin\0", 0.0);
    parms->zmax = iniparser_getdouble(ini, "cvm_h5repack:zmax\0",
                                      z0_cvm[3]*1.e-3);
    if (parms->lon0 < 0.0){parms->lon0 = parms->lon0 + 360.0;}
    if (parms->lon1 < 0.0){parms->lon1 = parms->lon1 + 360.0;}
    if (parms->lat0 < parms->lat0_cvm)
    {
        printf("%s: Changing minimum latitude to %f\n",
                  fcnm, parms->lat0_cvm);
        parms->lat0 = parms->lat0_cvm;
    }
    if (parms->lon0 < parms->lon0_cvm)
    {
        printf("%s: Changing minimum longitude to %f\n",
                  fcnm, parms->lon0_cvm);
        parms->lon0 = parms->lon0_cvm;
    }
    if (parms->lat1 > parms->lat1_cvm)
    {
        printf("%s: Changing maximum latitude to %f\n",
                  fcnm, parms->lat1_cvm);
        parms->lat1 = parms->lat1_cvm;
    }
    if (parms->lon1 > parms->lon1_cvm)
    {
        printf("%s: Changing maximum longitude to %f\n",
                  fcnm , parms->lon1_cvm);
        parms->lon1 = parms->lon1_cvm;
    } 
    parms->zmin = parms->zmin*1.e3; // convert to m
    parms->zmax = parms->zmax*1.e3; // convert to m
    if (parms->zmin != 0.0)
    {
        printf("%s: Error haven't done zmin > 0.0 yet\n", fcnm);
        return -1;
    }
    if (parms->zmax < parms->zmin)
    {
        printf("%s: Error zmax can't be shallower than zmin\n", fcnm);
        return -1;
    }
    // Compute the UTMS of the corresponding bounding box
    iway = 0;
    utm_geo(&parms->lon0, &parms->lat0,
            &parms->utm_x0, &parms->utm_y0,
            &parms->utmzone_cvm, &iway, &suppress);
    utm_geo(&parms->lon1, &parms->lat1,
            &parms->utm_x1, &parms->utm_y1,
            &parms->utmzone_cvm, &iway, &suppress);
    if (parms->lat0 > parms->lat1)
    {
        printf("%s: Error lat0 > lat1\n", fcnm);
        return -1;
    }
    // Get the thresholding
    parms->vp_min = 1.e20;
    parms->vp_max = 0.0;
    parms->lthresh_vp
         = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vp\0", false);
    if (parms->lthresh_vp)
    {
        parms->vp_min = iniparser_getdouble(ini,
                                            "cvm_h5repack:vp_min\0", 1.e20);
        parms->vp_max = iniparser_getdouble(ini,
                                            "cvm_h5repack:vp_max\0", 0.0);
    }
    parms->vs_min = 1.e20;
    parms->vs_min = 0.0;
    parms->lthresh_vs
        = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vs\0", false);
    if (parms->lthresh_vs)
    {
        parms->vs_min = iniparser_getdouble(ini, 
                                            "cvm_h5repack:vs_min\0", 1.e20);
        parms->vs_max = iniparser_getdouble(ini, 
                                            "cvm_h5repack:vs_max\0", 0.0);
    }
    parms->dens_min = 1.e20;
    parms->dens_max = 0.0;
    parms->lthresh_dens
        = iniparser_getboolean(ini, "cvm_h5repack:lthresh_dens\0", false);
    if (parms->lthresh_dens)
    {
        parms->dens_min = iniparser_getdouble(ini, 
                                              "cvm_h5repack:dens_min\0", 1.e20);
        parms->dens_max = iniparser_getdouble(ini, 
                                              "cvm_h5repack:dens_max\0", 0.0);
    }
    parms->vpvs_min = 1.71;
    parms->vpvs_max = 1.71;
    parms->lthresh_vpvs
       = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vpvs\0", false);
    if (parms->lthresh_vpvs)
    {
        parms->vpvs_min = iniparser_getdouble(ini,
                                              "cvm_h5repack:vpvs_min\0", 1.71);
        parms->vpvs_max = iniparser_getdouble(ini,
                                              "cvm_h5repack:vpvs_max\0", 1.71);
    }
    // set the Vp from the Vs
    parms->vpvs_ratio = iniparser_getdouble(ini, "cvm_h5repack:vpvs_ratio\0",
                                            1.71);
    if (parms->vpvs_ratio <= 1.0)
    {
        printf("%s: Invalid vpvs ratio %f\n", fcnm, parms->vpvs_ratio);
        return -1;
    }
    parms->setvp_from_vs = iniparser_getboolean(ini,
                                                "cvm_h5repack:setvp_from_vs\0",
                                                false);
    parms->setvs_from_vp = iniparser_getboolean(ini,
                                                "cvm_h5repack:setvs_from_vp\0",
                                                false);
    //------------------------------------------------------------------------//
    //                                cvm_mesh                                //
    //------------------------------------------------------------------------//
    strcpy(parms->mesh_outputdir, "./\0");
    s = iniparser_getstring(ini, "cvm_mesh:mesh_outputdir\0", "./\0");
    if (s != NULL)
    {
        if (strlen(s) > 0){
            memset(parms->mesh_outputdir, 0, sizeof(parms->mesh_outputdir));
            strcpy(parms->mesh_outputdir, s);
        }
    }
    if (os_makedirs(parms->mesh_outputdir) != 0)
    {
        printf("%s: Failed to make mesh output directory\n", fcnm);
        return -1;
    }
    // coarsening?
    parms->ncoarsen = iniparser_getint(ini, "cvm_mesh:ncoarsen\0", 0);
    if (parms->ncoarsen > 0)
    {
        parms->zcoarsen = (double *)
                          calloc((size_t) parms->ncoarsen, sizeof(double));
        parms->coarsen = (enum coarsen_type *)
                         calloc((size_t) parms->ncoarsen,
                                sizeof(enum coarsen_type));
        for (ic=0; ic<parms->ncoarsen; ic++)
        {
            memset(scoarsen, 0, sizeof(scoarsen));
            sprintf(scoarsen, "cvm_mesh:coarsen_%d", ic+1);
            s = iniparser_getstring(ini, scoarsen, NULL);
            if (s == NULL)
            {
                printf("%s: Error reading coarsening layer: %d\n",
                           fcnm, ic+1);
            }
            sscanf(s, "%lf %d", &parms->zcoarsen[ic], &ictype);
            parms->coarsen[ic] = (enum coarsen_type) ictype; 
        }
    } 
    parms->dx_fem = iniparser_getdouble(ini, "cvm_mesh:dx_fem\0", 2500.0);
    if (parms->dx_fem <= 0.0)
    {
        printf("%s: Error dz grid spacing must be positive \n", fcnm);
        return -1;
    }
    parms->dy_fem = iniparser_getdouble(ini, "cvm_mesh:dy_fem\0", 2500.0);
    if (parms->dy_fem <= 0.0)
    {
        printf("%s: Error dz grid spacing must be positive \n", fcnm);
        return -1;
    }
    parms->dz_fem = iniparser_getdouble(ini, "cvm_mesh:dz_fem\0", 2500.0);
    if (parms->dz_fem <= 0.0)
    {
        printf("%s: Error dz grid spacing must be positive \n", fcnm);
        return -1;
    }
    // Deform FEM mesh to topography?
    parms->ltopo
        = iniparser_getboolean(ini, "cvm_mesh:setTopography\0", false);
    if (parms->ltopo)
    {
        s = iniparser_getstring(ini, "cvm_mesh:topo_file\0", NULL);
        if (s == NULL)
        {
            printf("%s: Error topography file not specified\n", fcnm);
            return -1;
        }
        strcpy(parms->topofl, s);
        if (!os_path_isfile(parms->topofl)) 
        {
            printf("%s: Error topography file doesn't exist %s\n",
                       fcnm, parms->topofl);
            return -1;
        }
    }
    parms->utm_zone = iniparser_getint(ini, "cvm_mesh:utm_zone\0",
                                       parms->utmzone_cvm);
    if (parms->utm_zone < 1 || parms->utm_zone > 60)
    {
        printf("%s: Invalid UTM zone\n", fcnm);
        return -1;
    }
    if (parms->utm_zone != parms->utmzone_cvm)
    {
        printf("%s: Warning topography UTM zone inconsistent with CVM\n",
                  fcnm);
    }
    parms->ztopo_min = iniparser_getdouble(ini, "cvm_mesh:ztopo_min\0",
                                           (parms->zmin + parms->zmax)/2.0);
    // have ini parser free the ini dictionary
    iniparser_freedict(ini);
    return 0;
}
