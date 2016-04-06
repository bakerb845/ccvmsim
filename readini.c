#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iniparser.h>
#include <iscl/os/os.h>
#include <iscl/geodetic/geodetic.h>
#include <iscl/log/log.h>
#include "cvm.h"
#include "cvm_constants.h"

/*!
 * @brief Reads the CVM mesher initialization file
 *
 * @param[in] ini_file     name of ini file
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
    char ini_file[PATH_MAX];
    dictionary *ini;
    const char *s;
    double lat1, lon1, utmx, utmy, xoff_max, yoff_max;
    int k;
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
    parms->dxl_cvm = (double *)calloc(parms->nlay_cvm, sizeof(double));
    parms->dyl_cvm = (double *)calloc(parms->nlay_cvm, sizeof(double));
    parms->dzl_cvm = (double *)calloc(parms->nlay_cvm, sizeof(double));
    parms->z0_cvm  = (double *)calloc(parms->nlay_cvm, sizeof(double));
    parms->nxl_cvm = (int *)calloc(parms->nlay_cvm, sizeof(int));
    parms->nyl_cvm = (int *)calloc(parms->nlay_cvm, sizeof(int));
    parms->nzl_cvm = (int *)calloc(parms->nlay_cvm, sizeof(int));
    xoff_max = (double) (nxl_cvm[0]-1)*dxl_cvm[0];
    yoff_max = (double) (nyl_cvm[0]-1)*dyl_cvm[0];
    for (k=0; k<nlay_cvm; k++){
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
    if (ini == NULL){
        log_errorF("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    s = iniparser_getstring(ini, "cvm_h5repack:cvm_moddir\0", "./\0");
    strcpy(parms->cvm_moddir, s);
    if (!os_path_isdir(parms->cvm_moddir)){
        log_errorF("%s: CVM directory doesn't exist\n", fcnm);
        return -1;
    }
    // get the mesh directory
    s = iniparser_getstring(ini, "cvm_h5repack:mesh_outputdir", "./\0");
    strcpy(parms->mesh_outputdir, "./\0");
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
                                      z0_cvm[2]*1.e-3);
    if (parms->lon0 < 0.0){parms->lon0 = parms->lon0 + 360.0;}
    if (parms->lon1 < 0.0){parms->lon1 = parms->lon1 + 360.0;}
    if (parms->lat0 < parms->lat0_cvm){
        log_warnF("%s: Changing minimum latitude to %f\n",
                  fcnm, parms->lat0_cvm);
        parms->lat0 = parms->lat0_cvm;
    }
    if (parms->lon0 < parms->lon0_cvm){
        log_warnF("%s: Changing minimum longitude to %f\n",
                  fcnm, parms->lon0_cvm);
        parms->lon0 = parms->lon0_cvm;
    }
    if (parms->lat1 > parms->lat1_cvm){
        log_warnF("%s: Changing maximum latitude to %f\n",
                  fcnm, parms->lat1_cvm);
        parms->lat1 = parms->lat1_cvm;
    }
    if (parms->lon1 > parms->lon1_cvm){
        log_warnF("%s: Changing maximum longitude to %f\n",
                  fcnm , parms->lon1_cvm);
        parms->lon1 = parms->lon1_cvm;
    } 
    parms->zmin = parms->zmin*1.e3; // convert to m
    parms->zmax = parms->zmax*1.e3; // convert to m
    if (parms->zmin != 0.0){
        log_errorF("%s: Error haven't done zmin > 0.0 yet\n", fcnm);
        return -1;
    }
    if (parms->zmax < parms->zmin){
        log_errorF("%s: Error zmax can't be shallower than zmin\n", fcnm);
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
    if (parms->lat0 > parms->lat1){
        log_errorF("%s: Error lat0 > lat1\n", fcnm);
        return -1;
    }
    // Get the thresholding
    parms->vp_min = 1.e20;
    parms->vp_max = 0.0;
    parms->lthresh_vp
         = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vp\0", false);
    if (parms->lthresh_vp){
        parms->vp_min = iniparser_getdouble(ini,
                                            "cvm_h5repack:vp_min\0", 1.e20);
        parms->vp_max = iniparser_getdouble(ini,
                                            "cvm_h5repack:vp_max\0", 0.0);
    }
    parms->vs_min = 1.e20;
    parms->vs_min = 0.0;
    parms->lthresh_vs
        = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vs\0", false);
    if (parms->lthresh_vs){
        parms->vs_min = iniparser_getdouble(ini, 
                                            "cvm_h5repack:vs_min\0", 1.e20);
        parms->vs_max = iniparser_getdouble(ini, 
                                            "cvm_h5repack:vs_max\0", 0.0);
    }
    parms->dens_min = 1.e20;
    parms->dens_max = 0.0;
    parms->lthresh_dens
        = iniparser_getboolean(ini, "cvm_h5repack:lthresh_dens\0", false);
    if (parms->lthresh_dens){
        parms->dens_min = iniparser_getdouble(ini, 
                                              "cvm_h5repack:dens_min\0", 1.e20);
        parms->dens_max = iniparser_getdouble(ini, 
                                              "cvm_h5repack:dens_max\0", 0.0);
    }
    parms->vpvs_min = 1.71;
    parms->vpvs_max = 1.71;
    parms->lthresh_vpvs
       = iniparser_getboolean(ini, "cvm_h5repack:lthresh_vpvs\0", false);
    if (parms->lthresh_vpvs){
        parms->vpvs_min = iniparser_getdouble(ini,
                                              "cvm_h5repack:vpvs_min\0", 1.71);
        parms->vpvs_max = iniparser_getdouble(ini,
                                              "cvm_h5repack:vpvs_max\0", 1.71);
    }
    // set the Vp from the Vs
    parms->vpvs_ratio = iniparser_getdouble(ini, "cvm_h5repack:vpvs_ratio\0",
                                            1.71);
    if (parms->vpvs_ratio <= 1.0){
        log_errorF("%s: Invalid vpvs ratio %f\n", fcnm, parms->vpvs_ratio);
        return -1;
    }
    parms->setvp_from_vs = iniparser_getboolean(ini,
                                                "cvm_h5repack:setvp_from_vs\0",
                                                false);
    parms->setvs_from_vp = iniparser_getboolean(ini,
                                                "cvm_h5repack:setvs_from_vp\0",
                                                false);
    // have ini parser free the ini dictionary
    iniparser_freedict(ini);
    return 0;
}
