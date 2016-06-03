#include <stdbool.h>
#include <omp.h>
#include "cvm_struct.h"
#include "topo30.h"
#include "mesh.h"
#ifndef __CVM_H__
#define __CVM_H__
#ifdef __cplusplus
extern "C"
{
#endif

/* Unpack value from binary file */ 
float unpack_float(char *s, bool lswap);
double unpack_double(char *s, bool lswap);
/* Sets the free surface to 0 (m) while keeping the orientation of
   mesh positive up */
int cvm_cvm2freeSurfaceToZero(struct mesh_struct *mesh);
/* Extends bases in CVM so that there are no gaps between layers */
int cvm_extendBase(int nlay, struct cvm_model_struct *cvm_model);
/* Utility function for getting number of grid points in regular mesh */
int cvm_estimateNumberOfGridPoints(int nlay, struct cvm_model_struct *cvm_model,                                   double dx_fem, double dy_fem, double dz_fem,
                                   int *nx, int *ny, int *nz);
/* Put CVM material properties onto mesh material properties */
int cvm_cvm2meshMaterials(int nlay, struct cvm_model_struct *cvm_model,
                          struct mesh_struct *mesh);
/* Read CVM binary model */
int cvmio_readLayer(int lay,
                    struct cvm_parms_struct parms,
                    struct cvm_model_struct *cvm_model);
char *cvmio_readBinaryFile(char *filename, int *nbytes);
int cvmio_getBinaryFileSize(char *filename);
/* Write CVM model */
int cvmio_write__h5(char *dirnm, char *projnm,
                    int nlay, struct cvm_model_struct *cvm_model);
/* Utility routines for allocating and freeing memory */
int cvm_memory_allocate__element(enum mesh_element_type type,
                                 int nelem,
                                 struct mesh_element_struct *element);
void cvm_memory_free__model(int nlay, struct cvm_model_struct *cvm_model);
void cvm_memory_free__parms(struct cvm_parms_struct *parms);
/* Read the ini mesh file */
int cvm_readini(char *projnm, struct cvm_parms_struct *parms);
/* Estimate density from Vp */
#pragma omp declare simd
double density_Brocher(double Vp_ms);
double density_DarcyMcphee(double Vp_ms);
/* Estimate Qp qnd Qs from Vs */
#pragma omp declare simd
void qualityFactor_Frankel(double vs, double *Qp, double *Qs);
/* Convert to and from UTMs in zone 10 */
void utm_geo_utm2ll(double *rx4, double *ry4, int *UTM_PROJECTION_ZONE,
                    double *rlat4, double *rlon4);
void utm_geo_ll2utm(double *rlat4, double *rlon4, int *UTM_PROJECTION_ZONE,
                    double *rx4, double *ry4);
void utm_geo(double *rlon4, double *rlat4,
             double *rx4, double *ry4,
             int *UTM_PROJECTION_ZONE,
             int *iway, bool *SUPPRESS_UTM_PROJECTION);

#ifdef __cplusplus
}
#endif
#endif /* ifndef __CVM_H__ */
