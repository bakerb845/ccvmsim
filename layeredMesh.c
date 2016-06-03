#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "cvm.h"

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
    struct mesh_struct *mesh;
    struct mesh_element_struct tmplate[13];
    double *zint, dxl, dyl, dzl, z0l, z1l, zdep;
    int *nzint, ia, ic, ierr, il, indx, ipad, it, jl,
        nelem_last, nelemx, nelemy, nelemz, nfact, nnpg, nl, nx, ny, nz, nzl;
    const int nzmax = 100000000;
    //------------------------------------------------------------------------//
    //
    // Set the mesh structure for each layer
    dzl = dz;
    z0l = zmin;
    nl = nc + 1;
    nz = 0;
    dxl = dx;
    dyl = dy;
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
    // Now mesh to each interface starting deep and moving shallow
    nelemx = 1;
    nelemy = 1;
    jl = 0;
    for (il=nl-1; il>=0; il--)
    {
        z0l = fmax(zint[il] - zint[nl+1], 0.0);
        z1l = zint[il+1] - zint[nl+1];
        nz = nzint[il];
        nfact = 1;
        if (il > 0){nfact = coarsen[il-1];} 
        nx = nelemx + 1;
        ny = nelemy + 1;
        nelemz = nz - 1;
        printf("%s: Meshing layer: %d:\n", fcnm, il+1);
        printf("%s: This layer goes from %f to %f\n", fcnm, z0l, z1l);
        printf("%s: Grid spacing is %f %f %f\n",
              fcnm, dxl, dyl, dzl);
        printf("%s: Number of elements is %d %d %d\n",
               fcnm, nelemx, nelemy, nelemz);
        printf("\n");
        // Get number of elements
        ierr = regmesh_getNumberOfElements(nx, ny, nz, &mesh[jl].nelem);
        if (ierr != 0)
        {
            printf("%s: Error counting number of elements\n", fcnm);
            return -1;
        }
        // Set space leaving space for template
        ipad = 0;
        if (il > 0){ipad = 13 - 1;}
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
            printf("%s: Failed to set memory for pointers\n", fcnm);
            return -1;
        }
        ierr = regmesh_makeHexMeshPointers(nx, ny, nz, mesh[jl].nelem,
                                           mesh[jl].element);
        if (ierr != 0)
        {
            printf("%s: Error setting pointers!\n", fcnm);
            return -1;
        }
        // Apply the template
        if (ipad > 0)
        {
            // Get the nnpg offset.  Note, the regular mesher works z
            // positive up so take the second to last element in the mesh
            if (mesh[jl].nelem == 1)
            {
                nnpg = 0;
                nelem_last = 0; 
            }
            else
            {
                nnpg = mesh_element__getNumberOfAnchorNodes(true,
                                                            mesh[jl].nelem-1,
                                                            mesh[jl].element);
                nelem_last = mesh[jl].nelem - 1;
            }
            for (it=0; it<13; it++)
            {
                indx = mesh[jl].nelem + it;
                mesh[jl].element[nelem_last+it].nface = tmplate[it].nface;
                mesh[jl].element[nelem_last+it].ngnod = tmplate[it].ngnod;
                mesh[jl].element[nelem_last+it].ngnod_face
                    = tmplate[it].ngnod_face;
                mesh[jl].element[nelem_last+it].type = tmplate[it].type;
                for (ia=0; ia<mesh[jl].element[nelem_last].ngnod; ia++)
                {
                    mesh[jl].element[nelem_last+it].ien[ia]
                       = tmplate[it].ien[ia] + nnpg;
                }
                for (ia=0; ia<mesh[jl].element[nelem_last].ngnod; ia++)
                {
                    mesh[jl].element[nelem_last+it].bc[ia]
                       = tmplate[it].bc[ia]; 
                }
                //mesh[jl].element[indx];
printf("%d\n", nnpg); //mesh[jl].nnpg);
                
            }
        }
         
        // Set the regular mesh
        nelemx = nelemx*nfact;
        nelemy = nelemy*nfact;
        dxl = dxl/(double) nfact;
        dyl = dyl/(double) nfact;
        jl = jl + 1;
    }
return 0;
    for (ic=0; ic<nc; ic++)
    {
        zdep = zcoarsen[ic];
        ierr = regmesh_getNumberOfElements(nx, ny, nzl, &mesh[ic].nelem);
    }
    mesh_element_memory__free(13, tmplate);
    // Free the mesh
    for (il=0; il<nl; il++)
    {
        mesh_memory__free(&mesh[il]);
    }
    free(zint);
    free(nzint);
    free(mesh);
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
