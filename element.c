#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include "cvm.h"
#include "iscl/sorting/sorting.h"

int __mesh_element__setBoundaryConditionList(enum mesh_bc_enum side,
                                             int *bc_list)
{
    const char *fcnm = "__mesh_element__setBoundaryConditionList\0";
    bc_list[0] =-12345;
    bc_list[1] =-12345;
    bc_list[2] =-12345;
    if (side == SOUTH_BDRY)
    {
        bc_list[0] = SOUTH_BDRY;
        bc_list[1] = SOUTH_FREE;
        bc_list[2] = SOUTH_CE;
    }
    else if (side == EAST_BDRY)
    {
        bc_list[0] = EAST_BDRY;
        bc_list[1] = EAST_FREE;
        bc_list[2] = EAST_CE;
    }
    else if (side == NORTH_BDRY)
    {
        bc_list[0] = NORTH_BDRY;
        bc_list[1] = NORTH_FREE;
        bc_list[2] = NORTH_CE;
    }
    else if (side == WEST_BDRY)
    {
        bc_list[0] = WEST_BDRY;
        bc_list[1] = WEST_FREE;
        bc_list[2] = WEST_CE;
    }
    else if (side == BOTTOM_BDRY)
    {
        bc_list[0] = BOTTOM_BDRY;
        bc_list[1] = BOTTOM_FREE;
        bc_list[2] = BOTTOM_CE;
    }
    else if (side == TOP_BDRY)
    {
        bc_list[0] = TOP_BDRY;
        bc_list[1] = TOP_FREE;
        bc_list[2] = TOP_CE;
    }
    else
    {
        printf("%s: I do not know which side to get!\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Gets the size of the IEN array describing the side'th mesh side
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] side      desired side for which to extract length of IEN
 * @param[in] element   holds the boundary conditions for each element in mesh
 *                      [nelem]
 *
 * @result length of IENbdry structure for this side
 *
 */
int mesh_element__getBoundaryIENSize(int nelem,
                                     enum mesh_bc_enum side,
                                     struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element__getBoundaryIENSize\0";
    const int nbc = 3;
    int bc_list[3], ibc, iface, ielem, ierr, len_IENside;
    len_IENside = 0;
    ierr = __mesh_element__setBoundaryConditionList(side, bc_list);
    if (ierr != 0){
        printf("%s: Error setting boundary list\n", fcnm);
        return len_IENside;
    }
    // Loop on mesh
    #pragma omp parallel for firstprivate(nbc, nelem), \
     private(ielem, iface, ibc), shared(element), \
     reduction(+:len_IENside)
    for (ielem=0; ielem<nelem; ielem++){
        for (iface=0; iface<element[ielem].nface; iface++){
            if (element[ielem].bc[iface] == NO_BC){continue;}
            for (ibc=0; ibc<nbc; ibc++){
                if (bc_list[ibc] == element[ielem].bc[iface]){
                    len_IENside = len_IENside + element[ielem].ngnod_face;
                    break;
                }
            } // Loop on available boundary conditions
        } // Loop on faces on element
    } // Loop on elements
    return len_IENside;
}
//============================================================================//
/*!
 * @brief Gets the number of boundary elements on the side'th side
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] side      desired side for which to extract number of elements
 * @param[in] element   holds the boundary conditions for each element in mesh
 *                      [nelem]
 *
 * @result number of elements on the side'th side
 *
 */
int mesh_element__getNumberOfBoundaryElements(int nelem,
                                            enum mesh_bc_enum side,
                                            struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element__getNumberOfBoundaryElements\0";
    const int nbc = 3;
    int bc_list[3], ibc, iface, ielem, ierr, nelem_bdry;
    nelem_bdry = 0;
    ierr = __mesh_element__setBoundaryConditionList(side, bc_list);
    if (ierr != 0){ 
        printf("%s: Error setting boundary list\n", fcnm);
        return nelem_bdry;
    }   
    // Loop on mesh
    #pragma omp parallel for firstprivate(nbc, nelem), \
     private(ielem, iface, ibc), shared(element), \
     reduction(+:nelem_bdry)
    for (ielem=0; ielem<nelem; ielem++){
        for (iface=0; iface<element[ielem].nface; iface++){
            if (element[ielem].bc[iface] == NO_BC){continue;}
            for (ibc=0; ibc<nbc; ibc++){
                if (bc_list[ibc] == element[ielem].bc[iface]){
                    nelem_bdry = nelem_bdry + 1;
                    break;
                }
            } // Loop on available boundary conditions
        } // Loop on faces on element
    } // Loop on elements
    return nelem_bdry;
}
//============================================================================//
int mesh_element__getIENBoundary(int nelem, bool lhomog,
                                 enum mesh_bc_enum side,
                                 struct mesh_element_struct *element,
                                 int nelem_bdry, int *ien_bdry_ptr,
                                 int *bdry2glob_elem,
                                 int lenien_bdry, int *ien_bdry)
{
    const char *fcnm = "mesh_element__getIENBoundary\0";
    const int nbc = 3;
    int bc_list[3], ia, ibc, ielem, ielem_bdry, ierr,
        iface, indx_bdry, indx_face, ngnod_face;
    ierr = __mesh_element__setBoundaryConditionList(side, bc_list);
    if (ierr != 0){
        printf("%s: Error setting boundary list!\n", fcnm);
        return -1;
    }
    ien_bdry_ptr[0] = 0;
    if (lhomog){
        // Set the boundary element ien pointer
        ngnod_face = element[0].ngnod_face;
        #pragma omp simd
        for (ielem_bdry=1; ielem_bdry<nelem_bdry+1; ielem_bdry++)
        {
            ien_bdry_ptr[ielem_bdry] = ngnod_face*ielem_bdry;
        }
        // Extract the anchor nodes on the boundary
        ielem_bdry = 0;
        for (ielem=0; ielem<nelem; ielem++){
            for (iface=0; iface<element[ielem].nface; iface++){
                if (element[ielem].bc[iface] == NO_BC){continue;}
                indx_face = iface*element[ielem].ngnod_face;
                for (ibc=0; ibc<nbc; ibc++){
                    if (bc_list[ibc] == element[ielem].bc[iface]){
                        indx_bdry = ien_bdry_ptr[ielem_bdry];
                        for (ia=0; ia<element[ielem].ngnod_face; ia++){
                            ien_bdry[indx_bdry]
                               = element[ielem].ien_face[indx_face];
                            indx_face = indx_face + 1;
                            indx_bdry = indx_bdry + 1;
                        }
                        bdry2glob_elem[ielem_bdry] = ielem;
                        ielem_bdry = ielem_bdry + 1;
                        break;
                    } // End check on boundary condition match
                } // Loop on available boundary conditions
            } // Loop on faces on element
        } // Loop on elements
    }else{
printf("not done!\n");
return -1;
        ielem_bdry = 0;
        for (ielem=0; ielem<nelem; ielem++){
            if (element[ielem].bc[iface] == NO_BC){continue;}
            for (iface=0; iface<element[ielem].nface; iface++){
                for (ibc=0; ibc<nbc; ibc++){
                    if (bc_list[ibc] == element[ielem].bc[iface]){
                        bdry2glob_elem[ielem_bdry] = ielem;
                        ielem_bdry = ielem_bdry + 1;
                        break;
                    }
                }
            }
        } // Loop on elements
    }
    return 0;
}
//============================================================================//
int mesh_setBoundaryMesh(struct mesh_struct mesh,
                         enum mesh_bc_enum side,
                         struct mesh_struct *bdry)
{
    const char *fcnm = "mesh_setBoundaryMesh\0";
    memset(bdry, 0, sizeof(struct mesh_struct));
    // Get the size of the boundary mesh
    bdry->nelem = mesh_element__getNumberOfBoundaryElements(mesh.nelem,
                                                            side, mesh.element);
    if (bdry->nelem < 1){
        printf("%s: No elements on this side\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the length of the IEN array
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] lhomog    If True then the mesh has only one element type
 *                      If False then the mesh has mixed element types
 *
 * @result the length of the IEN array
 *
 */
int mesh_element__getIENSize(int nelem, bool lhomog,
                             struct mesh_element_struct *element)
{
    int ielem, len_ien;
    if (nelem < 1){return 0;}
    if (lhomog)
    {
        len_ien = nelem*element[0].ngnod;
    }
    else
    {
        len_ien = 0;
        for (ielem=0; ielem<nelem; ielem++)
        {
            len_ien = len_ien + element[ielem].ngnod;
        }
    }
    return len_ien;
}
//============================================================================//
/*!
 * @brief Builds the element node to global node map
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] len_ien   length of the ien array
 * @param[in] lhomog    If true then all elements in mesh are the same
 *                      If false then multiple element types exist in mesh
 * @param[in] element   holds the IEN array for each element and 
 *                      number of anchor nodes per element [nelem]
 *
 * @param[out] ien      maps from ia'th node on ielem'th element to 
 *                      global anchor node number [len_ien]
 *
 */
int mesh_element__getIEN(int nelem, int len_ien, bool lhomog,
                         struct mesh_element_struct *element,
                         int *ien_ptr, int *ien)
{
    const char *fcnm = "mesh_element__getIEN\0";
    int ia, ielem, indx, ngnod;
    indx = 0;
    ien_ptr[0] = 0;
    if (lhomog){
        ngnod = element[0].ngnod;
        #pragma omp simd
        for (ielem=0; ielem<nelem+1; ielem++){
            ien_ptr[ielem] = ngnod*ielem;
        } 
        #pragma omp parallel for collapse(2) private(ia, ielem, indx), \
         firstprivate(ngnod), shared(element, ien)
        for (ielem=0; ielem<nelem; ielem++){
            for (ia=0; ia<ngnod; ia++){
                indx = ielem*ngnod + ia;
                ien[indx] = element[ielem].ien[ia];
            }   
        }
    }else{
        for (ielem=0; ielem<nelem; ielem++){
            ngnod = element[ielem].ngnod;
            for (ia=0; ia<ngnod; ia++){
                ien[indx] = element[ielem].ien[ia];
                indx = indx + 1;
            }
            ien_ptr[ielem+1] = indx;
        }
    }
    if (ien_ptr[nelem] != len_ien){
        printf("%s: ien_ptr %d and len_ien %d inconsistent\n",
               fcnm, ien_ptr[nelem], len_ien);
        return -1;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the number of anchor nodes in the mesh from the
 *        element struct
 *
 * @param[in] cnum      If true then the element structure is C numbered
 *                      If false the the element structure is Fortran numbered
 * @param[in] nelem     number of elements in mesh
 * @param[in] element   element structure [nelem]
 *
 * @result number of anchor nodes in mesh
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__getNumberOfAnchorNodes(bool cnum, int nelem,
                                         struct mesh_element_struct *element)
{
    int ia, ielem, ngnod, nnpg;
    nnpg = 0;
    #pragma omp parallel for private(ia, ielem, ngnod), \
     shared(element), reduction(max:nnpg) 
    for (ielem=0; ielem<nelem; ielem++)
    {
        ngnod = element[ielem].ngnod;
        for (ia=0; ia<ngnod; ia++)
        {
            nnpg = fmax(nnpg, element[ielem].ien[ia]);
        }
    }
    if (cnum){nnpg = nnpg + 1;}
    return nnpg;
}
//============================================================================//
/*!
 * @brief Returns the (x,y,z) anchor node locations from the element
 *        structure
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] nnpg      number of anchor nodes in mesh
 * @param[in] element   element structure from which to extract (x,y,z)
 *                      triples onto the mesh anchor node location
 **                     arrays [nelem]
 *
 * @param[out] xlocs    x anchor node locations [nnpg]
 * @param[out] ylocs    y anchor node locations [nnpg]
 * @param[out] zlocs    z anchor node locations [nnpg]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__getAnchorNodeLocations(int nelem, int nnpg,
                                         struct mesh_element_struct *element,
                                         double *__restrict__ xlocs,
                                         double *__restrict__ ylocs,
                                         double *__restrict__ zlocs)
{
    const char *fcnm = "mesh_element__getAnchorNodeLocations\0";
    int ia, ielem, inpg, ngnod;
    if (nnpg < 1){
        printf("%s: Error no anchor nodes in mesh\n", fcnm);
        return -1;
    }
    #pragma omp simd
    for (inpg=0; inpg<nnpg; inpg++){
        xlocs[inpg] = 0.0;
        ylocs[inpg] = 0.0;
        zlocs[inpg] = 0.0;
    }
    #pragma omp parallel for private(ia, ielem, inpg, ngnod), \
     shared(element, xlocs, ylocs, zlocs)
    for (ielem=0; ielem<nelem; ielem++){
        ngnod = element[ielem].ngnod;
        for (ia=0; ia<ngnod; ia++){
            inpg = element[ielem].ien[ia];
            xlocs[inpg] = element[ielem].x[ia];
            ylocs[inpg] = element[ielem].y[ia];
            zlocs[inpg] = element[ielem].z[ia]; 
        } // Loop on anchor nodes
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @brief Extracts the desired boundary surface
 *
 * @param[in] nelem            number of elements in mesh
 * @param[in] bdry             desired boundary
 * @param[in] element          holds the element faces' boundary conditions
 *                             [nelem]
 * @param[in] nelem_bdry       number of elements on this boundary
 * @param[out] ien_bdry_ptr    maps from ielem_bdry'th element to start 
 *                             index of ien_bdry [nelem_bdry+1]
 *
 * @param[in] len_ien_bdry     length of ien_bdry
 * @param[out] ien_bdry        holds the global anchor node number of the 
 *                             ia'th anchor node on the ielem_bdry'th 
 *                             boundary element
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__getBoundarySurface(int nelem, enum mesh_bc_enum bdry,
                                     struct mesh_element_struct *element,
                                     int nelem_bdry, int *ien_bdry_ptr,
                                     int len_ien_bdry, int *ien_bdry)
{
    int ia, ielem, ielem_bdry, iface, indx, jndx;
    ielem_bdry = 0;
    ien_bdry_ptr[0] = 0;
    // Loop on elements
    for (ielem=0; ielem<nelem; ielem++){
        // Loop on element faces
        for (iface=0; iface<element[ielem].nface; iface++){
            // Element face matches the boundary condition
            if (element[ielem].bc[iface] == bdry){
                // Update pointer
                ien_bdry_ptr[ielem_bdry+1] = ien_bdry_ptr[ielem_bdry]
                                           + element[ielem].ngnod_face;
                indx = ien_bdry_ptr[ielem_bdry];
                jndx = iface*element[ielem].ngnod_face;
                ielem_bdry = ielem_bdry + 1;
                // Copy the anchor nodes
                for (ia=0; ia<element[ielem].ngnod_face; ia++){
                    ien_bdry[indx] = element[ielem].ien_face[jndx];
                    indx = indx + 1;
                    jndx = jndx + 1;
                }
            } // End check on boundary condition match
        } // Loop on element faces
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the array sizes for describing the boundary element
 *        surfaces
 *
 * @param[in] nelem           number of elements in mesh
 * @param[in] bdry            desired boundary 
 * @param[in] element         holds the element faces' boundary conditions
 *                            [nelem]
 *
 * @param[out] nelem_bdry     number of elements on boundary
 * @param[out] len_ien_bdry   length of the ien_bdry array
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 */
int mesh_element__getIENBoundarySize(int nelem, enum mesh_bc_enum bdry,
                                     struct mesh_element_struct *element,
                                     int *nelem_bdry, int *len_ien_bdry)
{
    int ielem, iface;
    *nelem_bdry = 0;
    *len_ien_bdry = 0;
    // Loop on elements
    for (ielem=0; ielem<nelem; ielem++)
    {
        // Loop on element faces
        for (iface=0; iface<element[ielem].nface; iface++)
        {
            // Element face matches the boundary condition
            if (element[ielem].bc[iface] == bdry)
            {
                *nelem_bdry = *nelem_bdry + 1;
                *len_ien_bdry = *len_ien_bdry + element[ielem].ngnod_face;
            }
        } // Loop on element faces
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @brief Generates a mapping from anchor nodes to elements to which they
 *        are connected
 *
 * @param[in] nelem               number of elements in mesh
 * @param[in] nnpg                number of anchor nodes in mesh
 * @param[in] element             holds the number of anchor nodes and
 *                                anchor node to global anchor node (ien)
 *                                map on each element [nelem]
 *
 * @param[out] node2element_ptr   maps from the inpg'th node to the
 *                                start index of node2element [nnpg+1]
 * @param[out] ierr               0 indicates success
 *
 * @result node2element[i1:i2] are the elements attached to the 
 *         inpg'th node
 *
 * @author Ben Baker, ISTI
 *
 */
int *mesh_element__getNode2ElementMap(int nelem, int nnpg,
                                      struct mesh_element_struct *element,
                                      int *node2element_ptr,
                                      int *ierr)
{
    const char *fcnm = "mesh_element__getNode2ElementMap\0";
    int *iwork, *isort, *node2element,
        i, i1, i2, ia, ielem, inpg, j, maxcons, nwork;
    // Initialize 
    *ierr = 0;
    iwork = NULL;
    isort = NULL; 
    node2element = NULL;
    if (nnpg < 1 || nelem < 1)
    {
        if (nelem < 1){printf("%s: Error no elements in mesh\n", fcnm);}
        if (nnpg < 1){printf("%s: Error no anchor nodes in mesh\n", fcnm);}
        *ierr = 1;
        goto ERROR;
    }
    iwork = (int *)calloc(nnpg, sizeof(int));
    // Compute the workspace size
    maxcons = 0;
    for (ielem=0; ielem<nelem; ielem++)
    {
        for (ia=0; ia<element[ielem].ngnod; ia++)
        {
            inpg = element[ielem].ien[ia];
            iwork[inpg] = iwork[inpg] + 1; 
            maxcons = fmax(iwork[inpg], maxcons);
        }
    }
    free(iwork);
    iwork = NULL;
    // Set the connectivity workspace
    iwork = (int *)calloc(nnpg*maxcons, sizeof(int));
    for (i=0; i<nnpg*maxcons; i++)
    {
        iwork[i] =-1;
    }
    // Fill the workspace with the node to element connectivity
    nwork = 0;
    for (ielem=0; ielem<nelem; ielem++)
    {
        for (ia=0; ia<element[ielem].ngnod; ia++)
        {
            inpg = element[ielem].ien[ia];
            i1 = maxcons*inpg;
            i2 = maxcons*(inpg + 1); 
            for (i=i1; i<i2; i++)
            {
                if (iwork[i] == ielem){break;} // Already have it
                // It's new so add it
                if (iwork[i] ==-1)
                {
                    iwork[i] = ielem;
                    nwork = nwork + 1;
                    break;
                }
            }
        }
    }
    // Now tally up the connectivity and store it
    isort = (int *)calloc(maxcons, sizeof(int));
    node2element = (int *)calloc(nwork, sizeof(int));
    node2element_ptr[0] = 0;
    for (inpg=0; inpg<nnpg; inpg++)
    {
        i1 = maxcons*inpg;
        i2 = maxcons*(inpg + 1);  
        j = 0;
        for (i=i1; i<i2; i++)
        {
            if (iwork[i] ==-1){break;}
            isort[j] = iwork[i];
            j = j + 1;
        }
        if (j == 0)
        {
            printf("%s: This is an internal error; j can't be zero\n", fcnm);
        }
        // Sort it for my sanity
        *ierr = __sorting_sort__int(j, isort, ASCENDING);
        if (*ierr != 0)
        {
            printf("%s: Error sorting %d points\n", fcnm, j);
            goto ERROR;
        }
        // Copy it
        node2element_ptr[inpg+1] = node2element_ptr[inpg] + j;
        i1 = node2element_ptr[inpg];
        i2 = node2element_ptr[inpg+1];
        j = 0;
        for (i=i1; i<i2; i++)
        {
            node2element[i1+j] = isort[j];
            j = j + 1;
        }
    }
ERROR:;
    if (iwork != NULL){free(iwork);}
    if (isort != NULL){free(isort);}
    if (*ierr != 0)
    {
        free(node2element);
        node2element = NULL;
    }
    return node2element;
}
//============================================================================//
/*!
 * @brief Returns the material properties (vp, vs, dens, ...) at the
 *        anchor nodes from the element structure
 *
 * @param[in] nelem     number of elements in mesh
 * @param[in] nnpg      number of anchor nodes in mesh
 * @param[in] element   element structure from which to extract material
 *                      properties onto the mesh anchor nodes material
 *                      arrays [nelem]
 *
 * @param[out] vp       compressional velocity at anchor nodes [nnpg]
 * @param[out] vs       shear velocity at anchor nodes [nnpg]
 * @param[out] dens     density at anchor nodes [nnpg]
 * @param[out] Qp       P quality factor at anchor nodes [nnpg]
 * @param[out] Qs       S quality factor at anchor nodes [nnpg]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__getAnchorNodeProperties(int nelem, int nnpg,
                                          struct mesh_element_struct *element,
                                          double *__restrict__ vp,
                                          double *__restrict__ vs,
                                          double *__restrict__ dens,
                                          double *__restrict__ Qp,
                                          double *__restrict__ Qs)
{
    const char *fcnm = "mesh_element__getAnchorNodeProperties\0";
    int ia, ielem, inpg, ngnod;
    if (nnpg < 1)
    {
        printf("%s: Error no anchor nodes in mesh\n", fcnm);
        return -1; 
    }
    #pragma omp simd
    for (inpg=0; inpg<nnpg; inpg++)
    {
        vp[inpg] = 0.0;
        vs[inpg] = 0.0;
        dens[inpg] = 0.0;
    }
    #pragma omp parallel for private(ia, ielem, inpg, ngnod), \
     shared(element, vp, vs, dens)
    for (ielem=0; ielem<nelem; ielem++)
    {
        ngnod = element[ielem].ngnod;
        for (ia=0; ia<ngnod; ia++)
        {
            inpg = element[ielem].ien[ia];
            vp[inpg]   = element[ielem].vp[ia];
            vs[inpg]   = element[ielem].vs[ia];
            dens[inpg] = element[ielem].dens[ia];
            Qp[inpg] = element[ielem].Qp[ia];
            Qs[inpg] = element[ielem].Qs[ia];
        } // Loop on anchor nodes
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the element anchor node locations
 *
 * @param[in] nelem       number of elements in mesh
 * @param[in] nnpg        number of anchor nodes in mesh
 * @param[in] xlocs       x anchor node locations [nnpg]
 * @param[in] ylocs       y anchor node locations [nnpg]
 * @param[in] zlocs       z anchor node locations [nnpg]
 *
 * @param[inout] element  on entry contains the element element structure
 *                        with number of anchor nodes per element and 
 *                        array space pre-allocated for x, y, and z 
 *                        on exit the x, y, and z arrays for each element
 *                        
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__setAnchorNodeLocations(int nelem, int nnpg,
                                         double *__restrict__ xlocs,
                                         double *__restrict__ ylocs,
                                         double *__restrict__ zlocs,
                                         struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element__setAnchorNodeLocations\0";
    int ia, ielem, inpg, ngnod;
    if (nelem < 1)
    {
        printf("%s: There are no elements in the mesh\n", fcnm);
        return -1;
    }
    if (nnpg < 1)
    {
        printf("%s: There are no anchor nodes in the mesh\n", fcnm);
        return -1;
    }
    #pragma omp parallel for private(ia, ielem, inpg, ngnod), \
     shared(element, xlocs, ylocs, zlocs)
    for (ielem=0; ielem<nelem; ielem++)
    {
        ngnod = element[ielem].ngnod;
        for (ia=0; ia<ngnod; ia++)
        {
            inpg = element[ielem].ien[ia];
            element[ielem].x[ia] = xlocs[inpg];
            element[ielem].y[ia] = ylocs[inpg];
            element[ielem].z[ia] = zlocs[inpg];
        } // Loop on anchor nodes
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @param[in] nelem       number of elements in mesh 
 * @param[in] vp          compressional velocity at anchor nodes [nnpg]
 * @param[in] vs          shear velocity at anchor nodes [nnpg]
 * @param[in] dens        density at anchor nodes [nnpg]
 * @param[in] Qp          P quality factor at anchor nodes [nnpg]
 * @param[in] Qs          S quality factor at anchor nodes [nnpg]
 *
 * @param[inout] element  on entry contains the element element structure
 *                        with number of anchor nodes per element and 
 *                        array space pre-allocated for vp, vs, and dens
 *                        on exit the vp, vs, dens, Qp, and Qs, arrays for
 *                        each element
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element__setAnchorNodeProperties(int nelem, int nnpg,
                                         double *__restrict__ vp,
                                         double *__restrict__ vs,
                                         double *__restrict__ dens,
                                         double *__restrict__ Qp,
                                         double *__restrict__ Qs,
                                         struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element__setAnchorNodeProperties\0";
    int ia, ielem, inpg, ngnod;
    if (nelem < 1)
    {
        printf("%s: There are no elements in the mesh\n", fcnm);
        return -1;
    }   
    if (nnpg < 1)
    {
        printf("%s: There are no anchor nodes in the mesh\n", fcnm);
        return -1;
    }   
    #pragma omp parallel for private(ia, ielem, inpg, ngnod), \
     shared(element, vp, vs, dens, Qp, Qs)
    for (ielem=0; ielem<nelem; ielem++)
    {
        ngnod = element[ielem].ngnod;
        for (ia=0; ia<ngnod; ia++)
        {
            inpg = element[ielem].ien[ia];
            element[ielem].vp[ia]   = vp[inpg];
            element[ielem].vs[ia]   = vs[inpg];
            element[ielem].dens[ia] = dens[inpg];
            element[ielem].Qp[ia] = Qp[inpg];
            element[ielem].Qs[ia] = Qs[inpg];
        } // Loop on anchor nodes
    } // Loop on elements
    return 0;
}
//============================================================================//
/*!
 * @brief Releses the memory on the element structure
 *
 * @param[in] nelem       number of elements in mesh (> 1)
 *
 * @param[inout] element  freed array of element structures [nelem]
 *
 * @author Ben Baker, ISTI
 *
 */
void mesh_element_memory__free(int nelem,
                               struct mesh_element_struct *element)
{
    int ielem;
    if (element == NULL){return;}
    for (ielem=0; ielem<nelem; ielem++)
    {
        if (element[ielem].vp       != NULL){free(element[ielem].vp);}
        if (element[ielem].vs       != NULL){free(element[ielem].vs);}
        if (element[ielem].dens     != NULL){free(element[ielem].dens);}
        if (element[ielem].Qp       != NULL){free(element[ielem].Qp);}
        if (element[ielem].Qs       != NULL){free(element[ielem].Qs);}
        if (element[ielem].x        != NULL){free(element[ielem].x);}
        if (element[ielem].y        != NULL){free(element[ielem].y);}
        if (element[ielem].z        != NULL){free(element[ielem].z);}
        if (element[ielem].ien_face != NULL){free(element[ielem].ien_face);}
        if (element[ielem].ien      != NULL){free(element[ielem].ien);}
        if (element[ielem].bc       != NULL){free(element[ielem].bc);}
        if (element[ielem].neighbor != NULL){free(element[ielem].neighbor);}
        memset(&element[ielem], 0, sizeof(struct mesh_element_struct));
    }
    return;
}
//============================================================================//
/*!
 * @brief Frees the memory on the mesh structure
 *
 * @param[inout] mesh      on input contains the mesh structure
 *                         on output all pointers on the mesh structure
 *                         have been freed
 *
 * @author Ben Baker, ISTI
 *
 */
void mesh_memory__free(struct mesh_struct *mesh)
{
    if (mesh->element != NULL)
    {
        mesh_element_memory__free(mesh->nelem, mesh->element);
        free(mesh->element);
        mesh->element = NULL;
    }
    if (mesh->xlocs != NULL){free(mesh->xlocs);}
    if (mesh->ylocs != NULL){free(mesh->ylocs);}
    if (mesh->zlocs != NULL){free(mesh->zlocs);}
    memset(mesh, 0, sizeof(struct mesh_struct));
    return;
}
//============================================================================//
/*!
 * @brief Determines if the mesh has mixed elements
 *
 * @param[in] nelem      number of elements in mesh
 * @param[in] element    holds the element type of each element [nelem]
 *
 * @result If True then the mesh has multiple elemen types
 *         Otherwise the mesh is of uniform element type
 *
 * @author Ben Baker, ISTI
 *
 */
bool mesh_element__ishomog(int nelem, struct mesh_element_struct *element)
{
    enum mesh_element_type type0; 
    int ielem;
    bool lhomog;
    lhomog = true;
    if (nelem < 1){return lhomog;}
    type0 = element[0].type;
    for (ielem=1; ielem<nelem; ielem++)
    {
        if (element[ielem].type != type0)
        {
            lhomog = false;
            break;
        } 
    }
    return lhomog;
}
//============================================================================//
/*!
 * @brief Returns the number of anchor nodes associated with the element type
 *
 * @param[in] type     element type
 * @param[in] lis3d    if true then this is a 3D mesh.
 *                     otherwise this is a 2D mesh
 *
 * @result number of anchor nodes on element
 *
 */
int mesh_element_type2numAnchorNodes(enum mesh_element_type type, bool lis3d)
{
    const char *fcnm = "mesh_element__type2numAnchorNodes\0";
    int ngnod;
    if (lis3d)
    {
        ngnod = 8; // hexes
        if (type == TET4)
        {
            ngnod = 4;
        }
        else if (type == HEX27)
        {
            ngnod = 27;
        }
        else if (type == TET8)
        {
            ngnod = 11; 
        }
        else
        {
            if (type != HEX8)
            {
                printf("%s: Defaulting to HEX8\n", fcnm);
            }
        }
    }
    else
    {
        printf("%s: Error 2d not yet done\n", fcnm);
        return -1;
    }
    return ngnod;
}
//============================================================================//
/*!
 * @brief Returns the number of faces for the element type
 *
 * @param[in] type     element type
 *
 * @result number of faces on the element
 *
 */
int mesh_element_type2numFaces(enum mesh_element_type type)
{
    const char *fcnm = "mesh_element_type2numFaces\0";
    int nface;
    nface = 6; // hexes have 6 faces
    if (type == TET4)
    {
        nface = 4;
    }
    else if (type == HEX27)
    {
        nface = 6;
    }
    else if (type == TET8)
    {
        nface = 4;
    }
    else
    {
        if (type != HEX8)
        {
            printf("%s: Defaulting to HEX8\n", fcnm);
        }
    }
    return nface;
}
//============================================================================//
/*!
 * @brief Returns the number of anchor nodes per element face for this
 *        element type
 *
 * @param[in] type     element type
 *
 * @result number of anchor nodes per face 
 *
 */
int mesh_element_type2numAnchorNodesPerFace(enum mesh_element_type type)
{
    const char *fcnm = "mesh_element_type2numAnchorNodesPerFace\0";
    int ngnod_face;
    ngnod_face = 4; // each face has 4 nodes
    if (type == TET4)
    {
        ngnod_face = 3;
    }
    else if (type == HEX27)
    {
        ngnod_face = 9;
    }
    else if (type == TET8)
    {
        ngnod_face = 6;
    }
    else
    {
        if (type != HEX8)
        {
            printf("%s: Defaulting to HEX8\n", fcnm);
        }
    }
    return ngnod_face;
}

//============================================================================//
/*!
 * @brief Allocates integer pointers for each element in the mesh 
 *
 * @param[in] type         defines the element type (use HEX8)
 * @param[in] nelem        number of elements in mesh (length of element)
 *
 * @param[inout] element   on input this is the nelem length element structure
 *                         on output each integer pointer has been allocated
 *                         [nelem]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element_memory__allocateIntegerPointers(
    enum mesh_element_type type,
    int nelem,
    struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element_memory__allocateIntegerPointers\0";
    int ia, ielem, ngnod, ngnod_face, nface, npface;
    ngnod = mesh_element_type2numAnchorNodes(type, true);
    nface = mesh_element_type2numFaces(type); 
    ngnod_face = mesh_element_type2numAnchorNodesPerFace(type);
    if (type != HEX8)
    {
        printf("%s: Error only considered hex8's\n", fcnm);
        return -1;
    }
    npface = ngnod_face*nface;
    for (ielem=0; ielem<nelem; ielem++)
    {
        element[ielem].ien_face = (int *)calloc(npface, sizeof(int));
        element[ielem].ien      = (int *)calloc(ngnod,  sizeof(int));
        element[ielem].bc       = (int *)calloc(nface,  sizeof(int));
        element[ielem].neighbor = (int *)calloc(nface,  sizeof(int));
        for (ia=0; ia<nface; ia++){element[ielem].neighbor[ia] =-1;}
        element[ielem].nface = nface;
        element[ielem].ngnod = ngnod;
        element[ielem].ngnod_face = ngnod_face;
        element[ielem].type = type;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Allocates material pointers for each element in the mesh.  This
 *        is probably okay for small meshes but for large meshes this is 
 *        a serious memory sink and one would be better suited putting the
 *        material and position vectors in a single array.
 *
 * @param[in] type         defines the element type (use HEX8)
 * @param[in] nelem        number of elements in mesh (length of element)
 *
 * @param[inout] element   on input this is the nelem length element structure
 *                         on output each material pointer has been allocated
 *                         [nelem]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int mesh_element_memory__allocateMaterialPointers(
    enum mesh_element_type type,
    int nelem,
    struct mesh_element_struct *element)
{
    const char *fcnm = "mesh_element_memory__allocateMaterialPointers\0";
    int ielem, ngnod;
    ngnod = mesh_element_type2numAnchorNodes(type, true);
    if (type != HEX8)
    {
        printf("%s: Error only considered hex8's\n", fcnm);
        return -1;
    }
    for (ielem=0; ielem<nelem; ielem++)
    {
        element[ielem].vp   = (double *)calloc(ngnod, sizeof(double));
        element[ielem].vs   = (double *)calloc(ngnod, sizeof(double));
        element[ielem].dens = (double *)calloc(ngnod, sizeof(double));
        element[ielem].Qp   = (double *)calloc(ngnod, sizeof(double));
        element[ielem].Qs   = (double *)calloc(ngnod, sizeof(double)); 
        element[ielem].x    = (double *)calloc(ngnod, sizeof(double));
        element[ielem].y    = (double *)calloc(ngnod, sizeof(double));
        element[ielem].z    = (double *)calloc(ngnod, sizeof(double));
    }
    return 0;
}
