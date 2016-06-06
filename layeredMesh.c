#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "cvm.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/sorting/sorting.h"
#include "visit_writer.h"

static double __tform(double a, double b,
                      double c, double d,
                      double x, int *ierr);

struct mesh_struct layeredMesh_healMesh(int nmeshes,
                                        struct mesh_struct *meshes,
                                        int *ierr)
{
    const char *fcnm = "layeredMesh_healMesh\0";
    struct mesh_struct mesh;
    double *x, *y, *z, xl, yl, zl;
    int *iperm, itemp[9], jtemp[9], i, ia, i1, i2, ielem, iface, imesh, inpg,
        indx, jelem, jface, nelem, nlocy, nlocz, nnpg, npsort, nwork;
    bool lfound, lmatch, lxdif, lydif, lzdif;
    const int faceMask[24] = {0, 1, 5, 4,
                              1, 2, 6, 5,
                              2, 3, 7, 6,
                              0, 4, 7, 3,
                              0, 3, 2, 1,
                              4, 5, 6, 7};
    const double tol = 1.e-2; // centimeter accuracy
    //------------------------------------------------------------------------//
    *ierr = 0;
    x = NULL;
    y = NULL;
    z = NULL;
    mesh.xlocs = NULL;
    mesh.ylocs = NULL;
    mesh.zlocs = NULL;
    iperm = NULL;
    memset(&mesh, 0, sizeof(struct mesh_struct));
    if (nmeshes < 1)
    {
        printf("%s: Error no meshes!\n", fcnm);
        *ierr = 1;
        goto ERROR;
    } 
    // Tally up the max number of anchor nodes
    nwork = 0;
    nelem = 0;
    for (imesh=0; imesh<nmeshes; imesh++)
    {
        nelem = nelem + meshes[imesh].nelem;
        for (ielem=0; ielem<meshes[imesh].nelem; ielem++)
        {
            nwork = nwork + meshes[imesh].element[ielem].ngnod;
        }
    }
    if (nwork < 1 || nelem < 1)
    {
        if (nwork < 1){printf("%s: No points in mesh!\n", fcnm);}
        if (nelem < 1){printf("%s: No elements in mesh!\n", fcnm);}
        *ierr = 1;
        goto ERROR;
    }
    mesh.nelem = nelem;
    mesh.element = (struct mesh_element_struct *)
                   calloc(mesh.nelem, sizeof(struct mesh_element_struct));
    x = (double *)calloc(nwork+1, sizeof(double));
    y = (double *)calloc(nwork+1, sizeof(double));
    z = (double *)calloc(nwork+1, sizeof(double));
    // Extract the anchor nodes
    i = 0;
    for (imesh=0; imesh<nmeshes; imesh++)
    {
        for (ielem=0; ielem<meshes[imesh].nelem; ielem++)
        {
            for (ia=0; ia<meshes[imesh].element[ielem].ngnod; ia++)
            {
                x[i] = meshes[imesh].element[ielem].x[ia];
                y[i] = meshes[imesh].element[ielem].y[ia];
                z[i] = meshes[imesh].element[ielem].z[ia];
                i = i + 1;
            }
        }
    }
    if (i != nwork)
    {
        printf("%s: Error lost count in nwork %d %d\n", fcnm, i, nwork);
        *ierr = 1;
        goto ERROR;
    }
    // Sort the anchor nodes with increasing z
    iperm = sorting_argsort__double(nwork, z, ASCENDING, ierr);
    if (*ierr != 0)
    {
        printf("%s: Error sorting z\n", fcnm);
        goto ERROR;
    } 
    *ierr  = sorting_applyPermutation__double(nwork, iperm, x, x);
    *ierr += sorting_applyPermutation__double(nwork, iperm, y, y);
    *ierr += sorting_applyPermutation__double(nwork, iperm, z, z);
    if (*ierr != 0)
    {
        printf("%s: Error applying permutation\n", fcnm);
        goto ERROR;
    }
    memory_free__int(&iperm);
    // Fix the last point so that we don't get a match
    x[nwork] =-DBL_MAX;
    y[nwork] =-DBL_MAX;
    z[nwork] =-DBL_MAX; 
    // Sort the anchor nodes with increasing y
    i1 = 0;
    for (i=0; i<nwork; i++)
    {
        if (fabs(z[i+1] - z[i]) > tol)
        {
            i2 = i;
            npsort = i2 - i1 + 1;
            if (npsort > 0)
            {
                iperm = sorting_argsort__double(npsort, &y[i1], ASCENDING, ierr);
                if (*ierr != 0)
                {
                    printf("%s: Error sorting y\n", fcnm);
                    goto ERROR;
                }
                *ierr  = sorting_applyPermutation__double(npsort, iperm,
                                                          &x[i1], &x[i1]);
                *ierr += sorting_applyPermutation__double(npsort, iperm,
                                                          &y[i1], &y[i1]);
                if (*ierr != 0)
                {
                    printf("%s: Error applying permutation 2\n", fcnm);
                    goto ERROR;
                }
                memory_free__int(&iperm);
            }
            i1 = i + 1;
        }
    }
    // Sort the anchor nodes with increasing x
    i1 = 0;
    for (i=0; i<nwork; i++)
    {
        if (fabs(z[i+1] - z[i]) > tol || fabs(y[i+1] - y[i]) > tol)
        {
            i2 = i;
            npsort = i2 - i1 + 1;
            if (npsort > 0)
            {
                *ierr = __sorting_sort__double(npsort, &x[i1], ASCENDING);
                if (*ierr != 0)
                {
                    printf("%s: Error sorting x\n", fcnm);
                    goto ERROR;
                }
            }
            i1 = i + 1;
        }
    }
    // Count the number of unique points and make a z pointer to speed up search
    nnpg = 0;
    nlocy = 0;
    nlocz = 0;
    for (i=0; i<nwork; i++)
    {
        lxdif = false;
        lydif = false;
        lzdif = false;
        if (fabs(x[i+1] - x[i]) > tol){lxdif = true;}
        if (fabs(y[i+1] - y[i]) > tol){lydif = true;}
        if (fabs(z[i+1] - z[i]) > tol){lzdif = true;}
        if (lxdif || lydif || lzdif)
        {
            nnpg = nnpg + 1;
            if (lydif){nlocy = nlocy + 1;}
            if (lzdif){nlocz = nlocz + 1;}
        }
    }
    if (nnpg < 1 || nlocy < 1 || nlocz < 1)
    {
        if (nnpg < 1){printf("%s: nnpg is zero\n", fcnm);}
        if (nlocy < 1){printf("%s: nlocy is zero\n", fcnm);}
        if (nlocz < 1){printf("%s: nlocz is zero\n", fcnm);}
        *ierr = 1;
        goto ERROR;
    }
    printf("%s: Found %d anchor nodes\n", fcnm, nnpg);
    // Extract them
    mesh.nnpg = nnpg;
    mesh.xlocs = (double *)calloc(nnpg, sizeof(double));
    mesh.ylocs = (double *)calloc(nnpg, sizeof(double));
    mesh.zlocs = (double *)calloc(nnpg, sizeof(double));
    inpg = 0;
    for (i=0; i<nwork; i++)
    {
        lxdif = false;
        lydif = false;
        lzdif = false;
        if (fabs(x[i+1] - x[i]) > tol){lxdif = true;}
        if (fabs(y[i+1] - y[i]) > tol){lydif = true;}
        if (fabs(z[i+1] - z[i]) > tol){lzdif = true;}
        if (lxdif || lydif || lzdif)
        {
            mesh.xlocs[inpg] = x[i];
            mesh.ylocs[inpg] = y[i];
            mesh.zlocs[inpg] = z[i];
            inpg = inpg + 1;
        }
    }
    // Because the mesh is small put locations on pointers in element
    mesh.lptr_only = true;
    // Only hexes
    mesh.lhomog = true;
    // Set the memory for the pointers on the elements
    *ierr = mesh_element_memory__allocateIntegerPointers(HEX8,
                                                         mesh.nelem,
                                                         mesh.element);
    if (*ierr != 0)
    {
        printf("%s: Error setting integer points\n", fcnm);
        goto ERROR;
    }
    // Now hunt for unique nodes
    printf("%s: Hunting for nodes...\n", fcnm);
printf("%s TODO use at least a zptr\n", fcnm);
    nelem = 0;
    for (imesh=0; imesh<nmeshes; imesh++)
    {
        for (ielem=0; ielem<meshes[imesh].nelem; ielem++)
        {
            for (ia=0; ia<meshes[imesh].element[ielem].ngnod; ia++)
            {
                xl = meshes[imesh].element[ielem].x[ia];
                yl = meshes[imesh].element[ielem].y[ia];
                zl = meshes[imesh].element[ielem].z[ia];
                lfound = false;
                for (inpg=0; inpg<mesh.nnpg; inpg++)
                {
                    if (fabs(mesh.xlocs[inpg] - xl) < tol &&
                        fabs(mesh.ylocs[inpg] - yl) < tol &&
                        fabs(mesh.zlocs[inpg] - zl) < tol)
                    {
                        lfound = true;
                        mesh.element[nelem].ien[ia] = inpg; 
                        break;
                    }
                }
                if (!lfound)
                {
                    printf("%s: Error failed to find point\n", fcnm);
                    goto ERROR;
                }
            } // Loop on anchor nodes
            nelem = nelem + 1;
        } // Loop on elements in mesh
    } // Loop on on meshes
    // This is going to hurt but find the neighbors
    printf("%s: Searching for neighbors...\n", fcnm);
    for (ielem=0; ielem<mesh.nelem; ielem++)
    {
        // don't have the right faceMask
        if (mesh.element[ielem].type != HEX8)
        {
            printf("%s: Need tet facemask\n", fcnm);
            *ierr = 1;
            goto ERROR;
        }
        for (iface=0; iface<mesh.element[ielem].nface; iface++)
        {
            // Already done
            if (mesh.element[ielem].neighbor[iface] >-1){continue;}
            // Get the anchor nodes on this face
            for (ia=0; ia<mesh.element[ielem].ngnod_face; ia++)
            {
                indx = iface*4 + ia;
                itemp[ia] = mesh.element[ielem].ien[faceMask[indx]];
            }
            *ierr = __sorting_sort__int(mesh.element[ielem].ngnod_face,
                                        itemp, ASCENDING);
            if (*ierr != 0){goto ERROR;}
            // Loop on other elements in mesh
            for (jelem=0; jelem<mesh.nelem; jelem++)
            {
                if (ielem == jelem){continue;} // Not attached to myself
                // This would indicate an unconforming mesh
                if (mesh.element[ielem].type != mesh.element[jelem].type)
                {
                    printf("%s: Mesh must be all hexes are tets\n", fcnm);
                    *ierr = 1;
                    goto ERROR;
                }
                for (jface=0; jface<mesh.element[jelem].nface; jface++)
                {
                    for (ia=0; ia<mesh.element[jelem].ngnod_face; ia++)
                    {
                        indx = iface*4 + ia;
                        jtemp[ia] = mesh.element[jelem].ien[faceMask[indx]];
                    }
                    *ierr = __sorting_sort__int(mesh.element[jelem].ngnod_face,
                                                jtemp, ASCENDING);
                    if (*ierr != 0){goto ERROR;}
                    lmatch = true;
                    for (ia=0; ia<mesh.element[jelem].ngnod_face; ia++)
                    {
                        if (itemp[ia] != jtemp[ia])
                        {
                            lmatch = false;
                            break; 
                        }
                    } // Loop on itemp == jtemp
                    if (lmatch)
                    {
                        mesh.element[ielem].neighbor[iface] = jelem;
                        //mesh.element[jelem].neighbor[jface] = ielem;
                    }
                } // Loop on faces on candidate neighbor
            } // Loop on other elements in mesh
        } // Loop on faces on home element 
    } // Loop on faces on home element
ERROR:;
    if (iperm != NULL){free(iperm);}
    if (x != NULL){free(x);}
    if (y != NULL){free(y);}
    if (z != NULL){free(z);}
    return mesh;
}
//============================================================================//
/*!
 * @brief Generates a template vertical column
 *
 * @param[in] nc      number of coarsening layers
 * @param[in] dx      grid spacing (m) in first layer in x
 * @param[in] dy      grid spacing (m) in first layer in y
 * @param[in] dz      grid spacing (m) in first layer in z
 * @param[in] zmax    max modeling depth (m)
 * @param[in] zmin    min modeling depth (m)
 *
 */
int layeredMesh_makeVerticalColumn(int nc,
                                   double dx, double dy, double dz,
                                   double zmax, double zmin,
                                   enum coarsen_type *coarsen, 
                                   double *zcoarsen )
{
    const char *fcnm = "layeredMesh_makeVerticalColumn\0";
    struct mesh_struct *mesh, meshColumn;
    struct mesh_element_struct *tmplate;
    double *zint, dxl, dyl, dzl, x, x1, y, y1, x2, y2,
           z, z0l, z1l, z1, z2;
    int *nzint, ia, ielemx, ielemy, ielemxy,
        ierr, iface, il, ioff, ipad, it, jl,
        lastElem, nelemx, nelemxy, nelemy, nelemz, nfact,
        nl, nx, ny, nz, nzl;
    const double x0 = 0.0;
    const double y0 = 0.0;
    const int nzmax = 100000000;
    const int nelem_template = 13;
    //------------------------------------------------------------------------//
    //
    // Set the mesh structure for each layer
    dzl = dz;
    z0l = zmin;
    nl = nc + 1;
    nz = 0;
    dxl = dx;
    dyl = dy;
    tmplate = (struct mesh_element_struct *)
              calloc(nelem_template, sizeof(struct mesh_element_struct));
    mesh = (struct mesh_struct *)calloc(nl, sizeof(struct mesh_struct));
    zint = (double *)calloc(nl+1, sizeof(double));
    nzint = (int *)calloc(nl, sizeof(int)); 
    // Set the vertical tripling template
    mesh_template__hexTriple(tmplate);
    // Figure out the widths of each layer and number of z points
    for (il=0; il<nl; il++)
    {
        nfact = 1;
        z1l = zmax;
        if (il < nl - 1)
        {
            nfact = coarsen[il];
            z1l = zcoarsen[il];
        }
        // Tally up the z grid points in this layer
        if (il < nl - 1)
        {
            for (nzl=0; nzl<nzmax; nzl++)
            {
                // safer to push the interface deeper
                if (z0l + (double) (nzl - 1)*dzl >= z1l)
                {
                    break;
                }
            }
            if (nzl >= nzmax)
            {
                printf("%s: Failed to count nz in layer %d\n", fcnm, il);
                return -1;
            }
        }
        else // Figure it out from the model base
        {
            nzl = (int) ((zmax - z0l)/dzl + 0.5) + 1;
            if ((double) (nzl - 1)*dzl > zmax - z0l){nzl = nzl - 1;}
        }
        // Make sure we finish with at least one point 
        if (nzl < 1)
        {
            printf("%s: Error no points in layer %d\n", fcnm, il+1);
            return -1;
        }
        z1l = z0l + (double) (nzl - 1)*dzl;
        zint[il+1] = z1l;
        nzint[il] = nzl;
        nz = nz + nzl;
        // Update for next layer
        z0l = z1l;
        dzl = (double) nfact*dzl;
        dxl = dxl*(double) nfact;
        dyl = dyl*(double) nfact;
    }
    // Now mesh to each interface starting deep and moving shallow.  
    // At the top of each layer apply refinement template, keep
    // track of the number of elements in each element, and anchor
    // node locations.  The subsequent mesh healing step will correct
    // all my errant bookkeeping by renumbering based on unique anchor nodes.
    nelemx = 1;
    nelemy = 1;
    jl = 0;
    for (il=nl-1; il>=0; il--)
    {
        z1l = fmax(zint[nl] - zint[il], 0.0);
        z0l = zint[nl] - zint[il+1];
        nz = nzint[il];
        nfact = 1;
        if (il > 0){nfact = coarsen[il-1];} 
        nx = nelemx + 1;
        ny = nelemy + 1;
        nelemz = nz - 1;
        if (il < nl - 1){printf("\n");}
        printf("%s: Meshing layer: %d\n", fcnm, il+1);
        printf("%s: This layer goes from %f to %f\n", fcnm, z0l, z1l);
        printf("%s: Layer thickness is %f\n", fcnm, z1l - z0l);
        printf("%s: Grid spacing is %f %f %f\n",
               fcnm, dxl, dyl, dzl);
        printf("%s: Number of elements is %d %d %d\n",
               fcnm, nelemx, nelemy, nelemz);
        // Get number of elements
        ierr = regmesh_getNumberOfElements(nx, ny, nz, &mesh[jl].nelem);
        if (ierr != 0)
        {
            printf("%s: Error counting number of elements\n", fcnm);
            return -1;
        }
        // Set space leaving space for template
        ipad = 0;
        nelemxy = 0;
        if (il > 0)
        {
            nelemxy = nelemx*nelemy;
            // top element is replaced so don't double count
            ipad = nelemxy*(nelem_template - 1);
        }
        mesh[jl].element = (struct mesh_element_struct *)
                           calloc(mesh[jl].nelem + ipad,
                                  sizeof(struct mesh_element_struct));
        // Because the mesh is small put locations on pointers in element
        mesh[jl].lptr_only = false;
        // Homogeneous hex mesh in this layer
        mesh[jl].lhomog = true;
        // Set the memory for the pointers on the elements
        ierr = mesh_element_memory__allocateIntegerPointers(HEX8,
                                                          mesh[jl].nelem + ipad,
                                                          mesh[jl].element);
        if (ierr != 0)
        {
            printf("%s: Failed to set memory for integer pointers\n", fcnm);
            return -1;
        }
        ierr = mesh_element_memory__allocateMaterialPointers(
                                        HEX8,
                                        mesh[jl].nelem + ipad,
                                        mesh[jl].element); 
        if (ierr != 0)
        {
            printf("%s: Error setting space for material pointers\n", fcnm);
            return -1;
        }
        // Set the integer pointers
        ierr = regmesh_makeHexMeshPointers(nx, ny, nz, mesh[jl].nelem,
                                           mesh[jl].element);
        if (ierr != 0)
        {
            printf("%s: Error setting integer pointers!\n", fcnm);
            return -1;
        }
        // Set the anchor node locations
        ierr = regmesh_makeRegularNodes(nx, ny, nz,
                                        dxl,
                                        dyl,
                                        dzl,
                                        x0, y0, z0l,
                                        mesh[jl].nelem, mesh[jl].element);
        if (ierr != 0)
        {
            printf("%s: Failed setting anchor node locations\n", fcnm);
            return -1;
        }
        // Apply the template to the top element in this layer
        if (ipad > 0)
        {
            // Get the nnpg offset.  Note, the regular mesher works z
            // positive up so take the second to last element in the mesh
            if (mesh[jl].nelem == 1)
            {
                lastElem = 0;
            }
            else
            {
                lastElem = (nelemz - 1)*nelemxy;
            }
            // Figure out the min/max z locations in the top element
            z1 = array_min__double(mesh[jl].element[lastElem].ngnod,
                                   mesh[jl].element[lastElem].z);
            z2 = array_max__double(mesh[jl].element[lastElem].ngnod,
                                   mesh[jl].element[lastElem].z);
            // Copy the templates and fix the anchor nodes and element numbers
            // for the top elements in the (x, y) plane
            ielemx =-1;
            ielemy = 0;
            ioff = lastElem;
            for (ielemxy=0; ielemxy<nelemxy; ielemxy++)
            {
                ielemx = ielemx + 1;
                if (ielemx == nelemx)
                {
                    ielemy = ielemy + 1;
                    ielemx = 0; 
                }
                // Figure out min and max (x, y) locations for this element
                x1 = array_min__double(mesh[jl].element[ielemxy].ngnod,
                                       mesh[jl].element[ielemxy].x);
                x2 = array_max__double(mesh[jl].element[ielemxy].ngnod,
                                       mesh[jl].element[ielemxy].x);
                y1 = array_min__double(mesh[jl].element[ielemxy].ngnod,
                                       mesh[jl].element[ielemxy].y);
                y2 = array_max__double(mesh[jl].element[ielemxy].ngnod,
                                       mesh[jl].element[ielemxy].y);

                for (it=0; it<nelem_template; it++)
                {
                    // Basic information
                    mesh[jl].element[ioff+it].nface = tmplate[it].nface;
                    mesh[jl].element[ioff+it].ngnod = tmplate[it].ngnod;
                    mesh[jl].element[ioff+it].ngnod_face
                        = tmplate[it].ngnod_face;
                    mesh[jl].element[ioff+it].type = tmplate[it].type;
                    // Anchor node numbers
                    for (ia=0; ia<tmplate[it].ngnod; ia++)
                    {
                        mesh[jl].element[ioff+it].ien[ia]
                           = tmplate[it].ien[ia];
                    }
                    // Boundary conditions and neighbors
                    for (iface=0; iface<tmplate[it].nface; iface++)
                    {
                        mesh[jl].element[ioff+it].bc[iface]
                           = tmplate[it].bc[iface]; 
                        mesh[jl].element[ioff+it].neighbor[iface]
                           = tmplate[it].neighbor[iface];
                    }
                    // Fix the anchor node locations by mapping from the 
                    // unit box to the box defined by the mesh
                    for (ia=0; ia<tmplate[it].ngnod; ia++)
                    {
                        x = __tform(-1.0, 1.0,
                                    x1, x2,
                                    tmplate[it].x[ia], &ierr);
                        y = __tform(-1.0, 1.0,
                                    y1, y2,
                                    tmplate[it].y[ia], &ierr);
                        z = __tform(-1.0, 1.0,
                                    z1, z2,
                                    tmplate[it].z[ia], &ierr);
                        mesh[jl].element[ioff+it].x[ia] = x;
                        mesh[jl].element[ioff+it].y[ia] = y;
                        mesh[jl].element[ioff+it].z[ia] = z;
                    }
                }
                ioff = ioff + nelem_template;
            } // loop on elements in (x, y) plane
        }
        else // Not refining top element in this layer
        {
            if (il != 0)
            {
                printf("%s: Warning not refining layer %d\n", fcnm, jl+1);
            }
        }
        // Fix the pointers in the previous layer
        mesh[jl].nelem = mesh[jl].nelem + ipad;
        printf("%s: Number of elements in this layer: %d\n",
               fcnm, mesh[jl].nelem);
        // Update the number of elements and grid spacing 
        nelemx = nelemx*nfact;
        nelemy = nelemy*nfact;
        dxl = dxl/(double) nfact;
        dyl = dyl/(double) nfact;
        dzl = dzl/(double) nfact;
        jl = jl + 1;
    }
    // Heal the mesh
    printf("%s: Healing the columnar mesh...\n", fcnm);
    meshColumn = layeredMesh_healMesh(nl,
                                      mesh,
                                      &ierr);
    // Generate connectivity
    float *pts = (float *)calloc(3*meshColumn.nnpg, sizeof(float));
    int *zonetypes = (int *)calloc(meshColumn.nelem, sizeof(int));
    float *elems = (float *)calloc(meshColumn.nelem, sizeof(float));
    float *nodes = (float *)calloc(meshColumn.nnpg, sizeof(float));
    int *connectivity = NULL;
    int ielem, inpg, ncon, nnpg;
    for (ielem=0; ielem<meshColumn.nelem; ielem++)
    {
        zonetypes[ielem] = VISIT_HEXAHEDRON;
        elems[ielem] = (float) ielem + 1; 
    }
    for (inpg=0; inpg<meshColumn.nnpg; inpg++)
    {
        nodes[inpg] = (float) inpg + 1;
        pts[3*inpg+0] = meshColumn.xlocs[inpg];
        pts[3*inpg+1] = meshColumn.ylocs[inpg];
        pts[3*inpg+2] = meshColumn.zlocs[inpg];
    }
    ncon = meshio_mesh2connectivity__getSize(meshColumn.lhomog,
                                             meshColumn.nelem,
                                             meshColumn.element);
    connectivity = (int *)calloc(ncon, sizeof(int));
    ierr = meshio_mesh2connectivity(true, meshColumn.lhomog, meshColumn.nelem,
                                    meshColumn.element,
                                    ncon, connectivity);
    int nvars = 2;
    int vardims[2] = {1, 1};
    int centering[2] = {0, 1};
    const char *varnames[2] = {"ElementNumbers\0", "NodeNumbers\0"};
    float *vars[] = {elems, nodes};
    write_unstructured_mesh("column.vtk\0", 1, meshColumn.nnpg,
                            pts, meshColumn.nelem,
                            zonetypes, connectivity, nvars, vardims, centering,
                            varnames, vars);
    free(pts);
    free(zonetypes);
    free(elems);
    free(nodes);
    free(connectivity);
    // Free space
    printf("%s: Freeing memory...\n", fcnm);
    mesh_element_memory__free(nelem_template, tmplate);
    free(tmplate);
    for (il=0; il<nl; il++)
    {
        mesh_memory__free(&mesh[il]);
    }
    free(mesh);
    free(zint);
    free(nzint);
    return 0;
}
//============================================================================//
/*!
 * @brief Generates the layered mesh with corasening in the vertical
 *        direction
 */
int layeredMesh_driver(struct cvm_parms_struct parms,
                       struct cvm_model_struct *cvm_model,
                       struct mesh_struct *mesh)
{
    const char *fcnm = "layeredMesh_driver\0";
    double rlat, rlon, x1, y1;
    int ic, ierr, ix, iy, ncom_denom, nx, nx0, ny, ny0, nz;
    memset(mesh, 0, sizeof(struct mesh_struct));
    // First fix the number of x and y grid points in the top layer
    ncom_denom = 1;
    for (ic=0; ic<parms.ncoarsen; ic++)
    {
        ncom_denom = ncom_denom*parms.coarsen[ic];
    }
    if (ncom_denom < 1)
    {
        printf("%s: Invalid common denominator %d\n", fcnm, ncom_denom);
        return -1;
    }
    // Estimate nx and ny grid points in shallowest regular mesh 
    ierr = cvm_estimateNumberOfGridPoints(parms.nlay_cvm, cvm_model,
                                          parms.dx_fem,
                                          parms.dy_fem,
                                          parms.dz_fem,
                                          &nx, &ny, &nz);
    if (ierr != 0)
    {
        printf("%s: Error estimating nx/ny grid points in top layer\n", fcnm);
        return -1;
    }
    // Now walk the nx and ny grid points so they are divisible by ncom_denom
    nx0 = nx;
    for (ix=0; ix<nx; ix++)
    {
        if (fmod(nx-ix, ncom_denom) == 0)
        {
            nx = nx - ix;
            break;
        }
    }
    ny0 = ny;
    for (iy=0; iy<ny; iy++)
    {
        if (fmod(ny-iy, ncom_denom) == 0)
        {   
            ny = ny - iy; 
            break;
        }   
    } 
    if (nx0 != nx || ny0 != ny)
    {
        x1 = parms.utm_x0 + (double) (nx - 1)*parms.dx_fem;
        y1 = parms.utm_y0 + (double) (ny - 1)*parms.dy_fem;
        utm_geo_utm2ll(&x1, &y1, &parms.utm_zone, &rlat, &rlon);
        if (rlon < 0.0){rlon = rlon + 360.0;}
        printf("%s: Moving (x1,y1) to (%f, %f)\n", fcnm, rlat, rlon);
    }
    // Make the vertical column mesh
     layeredMesh_makeVerticalColumn(parms.ncoarsen, 
                                    parms.dx_fem, parms.dy_fem, parms.dz_fem, 
                                    parms.zmax, parms.zmin,
                                    parms.coarsen, 
                                    parms.zcoarsen );
printf("%f %f\n", parms.zmax, parms.zmin);
printf("%d %d\n", nx, ny);
    //nx = parmsutm_x0
    //parms->utm_x0 + parms->dx_fem;
    
    // Walk the y grid points back so the ensuing interpolation works

    return 0;
}
//============================================================================//

static double __tform(double a, double b,
                      double c, double d,
                      double x, int *ierr)
{
    const char *error = "__tform: Determinant undefined\n\0";
    double c1, c2, det, xi; 
    *ierr = 0;
    if (a == b)
    {   
        perror(error);
        *ierr = 1;
        return 0.0;
    }   
    det = 1.0/(b - a); 
    c1 = det*(b*c - a*d);
    c2 = det*(d - c); 
    xi = c1 + x*c2;
    return xi; 
}

/*!
 * @brief Maps a point on the unit box into the isoparametric element
 *        for a 8 node hexahedral element.  The hexahedron is meshed
 *        using the convetion of T.J.R. Hughes - The Finite Element
 *        Method pg 136 (2001).  
 *
 * @param[in] xi      xi point (analogous to x) on unit box [-1,1]
 * @param[in] eta     eta point (analogous to y) on unit box [-1,1]
 * @param[in] zeta    zeta point (analogous to z) on unit box [-1,1]
 * @param[in] x       anchor node locations on unit box ordered in a
 *                    form that is consistent with Hughes 2001 [8].
 *
 * @result interpolated point from unit box onto isoparametric 
 *
 */
double element_isoparametric_hex8__evaluateShapeFunction(
    const double xi,
    const double eta,
    const double zeta,
    const double *__restrict__ x)
{
    double t1, t2, t3, t4, t5, t6, t7, t8, xint;
    const double  one = 1.0;
    const double eighth = 0.125;
    t1 = eighth*(one - xi)*(one - eta)*(one - zeta);
    t2 = eighth*(one + xi)*(one - eta)*(one - zeta);
    t3 = eighth*(one + xi)*(one + eta)*(one - zeta);
    t4 = eighth*(one - xi)*(one + eta)*(one - zeta);
    t5 = eighth*(one - xi)*(one - eta)*(one + zeta);
    t6 = eighth*(one + xi)*(one - eta)*(one + zeta);
    t7 = eighth*(one + xi)*(one + eta)*(one + zeta);
    t8 = eighth*(one - xi)*(one + eta)*(one + zeta);
    xint = t1*x[0] + t2*x[1] + t3*x[2] + t4*x[3]
         + t5*x[4] + t6*x[5] + t7*x[6] + t8*x[7];
    return xint;
}
