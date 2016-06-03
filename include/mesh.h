#include <stdbool.h>
#include "cvm_struct.h"
#ifndef __MESH_H__
#define __MESH_H__
#ifdef __cplusplus
extern "C" 
{
#endif

/* Vertical coarsening mesh prototypes */
int layeredMesh_driver(struct cvm_parms_struct parms,
                       struct cvm_model_struct *cvm_model,
                       struct mesh_struct *mesh);

/* Regular mesh prototypes */
int regmesh_getNumberOfElements(int nx, int ny, int nz, int *nelem);
int regmesh_makeHexMeshPointers(int nx, int ny, int nz, 
                                int nelem, struct mesh_element_struct *element);
int regmesh_getNumberOfAnchorNodes(int nx, int ny, int nz, int *nnpg);
int mesh_element__getNumberOfBoundaryElements(int nelem,
                                           enum mesh_bc_enum side,
                                           struct mesh_element_struct *element);
int mesh_element__getIENBoundary(int nelem, bool lhomog,
                                 enum mesh_bc_enum side,
                                 struct mesh_element_struct *element,
                                 int nelem_bdry, int *ien_bdry_ptr,
                                 int *bdry2glob_elem,
                                 int lenien_bdry, int *ien_bdry);
int regmesh_makeRegularNodes(int nx, int ny, int nz,
                             double dx, double dy, double dz,
                             double x0, double y0, double z0,
                             int nelem,
                             struct mesh_element_struct *element);
void __regmesh_makeRegularNodes(int nx, int ny, int nz, 
                                double dx, double dy, double dz, 
                                double x0, double y0, double z0, 
                                double *__restrict__ xlocs,
                                double *__restrict__ ylocs,
                                double *__restrict__ zlocs);
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
                                double *__restrict__ Qsw);
int regmesh_copyRegularModel(bool lflip, int nx, int ny, int nz,
                             int nelem, int nnpg,
                             double *__restrict__ vp,
                             double *__restrict__ vs,
                             double *__restrict__ dens,
                             double *__restrict__ Qp,
                             double *__restrict__ Qs,
                             struct mesh_element_struct *element);
int regmesh_makeHexIEN(int nx, int ny, int nz, int *ien);

int mesh_element__getAnchorNodeLocations(int nelem, int nnpg,
                                         struct mesh_element_struct *element,
                                         double *__restrict__ xlocs,
                                         double *__restrict__ ylocs,
                                         double *__restrict__ zlocs);
int mesh_element__getAnchorNodeProperties(int nelem, int nnpg,
                                          struct mesh_element_struct *element,
                                          double *__restrict__ vp, 
                                          double *__restrict__ vs, 
                                          double *__restrict__ dens,
                                          double *__restrict__ Qp,
                                          double *__restrict__ Qs);
int mesh_element__getNumberOfAnchorNodes(bool cnum, int nelem,
                                         struct mesh_element_struct *element);
int mesh_element__setAnchorNodeLocations(int nelem, int nnpg,
                                         double *__restrict__ xlocs,
                                         double *__restrict__ ylocs,
                                         double *__restrict__ zlocs,
                                         struct mesh_element_struct *element);
int mesh_element__setAnchorNodeProperties(int nelem, int nnpg,
                                         double *__restrict__ vp, 
                                         double *__restrict__ vs, 
                                         double *__restrict__ dens,
                                         double *__restrict__ Qp,
                                         double *__restrict__ Qs,
                                         struct mesh_element_struct *element);
int mesh_element_type2numAnchorNodesPerFace(enum mesh_element_type type);
int mesh_element_type2numFaces(enum mesh_element_type type);
int mesh_element_type2numAnchorNodes(enum mesh_element_type type);

bool mesh_element__ishomog(int nelem, struct mesh_element_struct *element);
int mesh_element_memory__allocateMaterialPointers(
    enum mesh_element_type type,
    int nelem,
    struct mesh_element_struct *element);
int mesh_element_memory__allocateIntegerPointers(
    enum mesh_element_type type,
    int nelem,
    struct mesh_element_struct *element);
void mesh_element_memory__free(int nelem,
                               struct mesh_element_struct *element);
void mesh_memory__free(struct mesh_struct *mesh);
int mesh_element__getIENSize(int nelem, bool lhomog,
                             struct mesh_element_struct *element);
int mesh_element__getIEN(int nelem, int len_ien, bool lhomog,
                         struct mesh_element_struct *element,
                         int *ien_ptr, int *ien);
int mesh_element__getBoundaryIENSize(int nelem,
                                     enum mesh_bc_enum side,
                                     struct mesh_element_struct *element);
void mesh_template__hexTriple(struct mesh_element_struct *tmplate);
int meshio_write__NLLGrid(char *dirnm, char *projnm,
                          int nelemx, int nelemy, int nelemz,
                          double x0, double y0,
                          double dx, double dy, double dz, 
                          int utm_zone, double dep0,
                          struct mesh_struct mesh);
int meshio_write__h5(char *meshdir, char *projnm,
                     struct mesh_struct mesh);
int meshio_write__xdmf(char *meshdir, char *projnm,
                       int nelem, int ngnod, int nnpg);
int meshio_write__setFilename(char *meshdir, char *projnm, char *app,
                              char fname[PATH_MAX]);
int meshio_write__specfem3d(char *meshdir, char *projnm,
                            struct mesh_struct mesh);

#ifdef __cplusplus
}
#endif
#endif /* __MESH_H__ */
