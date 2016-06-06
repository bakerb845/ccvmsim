#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include <omp.h>
#include <hdf5.h>
#include "cvm.h"
#include "h5_cinter.h"
#include "iscl/os/os.h"

int meshio_write__setFilenameNLL(char *dirnm, char *projnm, char *file_type,
                                 bool isbuf, char fname[PATH_MAX]);

//============================================================================//
/*!
 * @brief Comptues the size of the connectivity array for writing to HDF5
 *        and contextualizing with XDMF for viewing with Paraview
 *
 * @param[in] lhomog     If True then all element types are the same
 *                       If False then the mesh contains different element
 *                       types
 * @param[in] nelem      number of elements in mesh
 * @param[in] element    element structure holding number of anchor nodes
 *                       per element
 *
 * @result size of connnectivity array
 *
 * @author Ben Baker, ISTI
 *
 */
int meshio_mesh2connectivity__getSize(bool lhomog, int nelem,
                                      struct mesh_element_struct *element)
{
    int iadd, ielem, ncon, ngnod;
    iadd = 0;
    if (!lhomog){iadd = 1;}
    ncon = 0;
    #pragma omp parallel for private(ielem, ngnod), \
     firstprivate(iadd), reduction(+:ncon)
    for (ielem=0; ielem<nelem; ielem++)
    {
        ngnod = element[ielem].ngnod;
        ncon = ncon + ngnod + iadd;
    }
    return ncon;
}
//============================================================================//
/*!
 * @brief Computes the connectivity array for writing to HDF5 and 
 *         contextualizing with XDMF for viewing with Paraview
 *
 * @param[in] cnum            If True -> the elements use a C numbering
 *                            (i.e. they start counting at index 0)
 *                            If False -> the elements use a Fortran numbering
 *                            (i.e. they start counting at index 1)
 * @param[in] lhomog          If True -> the mesh is homogeneous
 *                            If False -> the mesh contains multiple element
 *                            types
 * @param[in] ncon            length of the connectivity array
 *
 * @param[out] connectivity   holds the element connectivity for writing
 *                            to HDF5 [ncon]
 *
 * @result 0 indicates success
 *
 */ 
int meshio_mesh2connectivity(bool cnum, bool lhomog, int nelem,
                             struct mesh_element_struct *element,
                             int ncon, int *connectivity)
{
    const char *fcnm = "meshio_mesh2connectivity\0";
    int *ien, *ien_ptr, ia, ielem, ierr, indx, ishift, len_ien, ngnod;
    if (nelem < 1)
    {
        printf("%s: No elements in mesh\n", fcnm);
        return -1;
    }
    ishift = 0;
    if (!cnum){ishift = 1;}
    // Get the size of the IEN array
    len_ien = mesh_element__getIENSize(nelem, lhomog, element);
    ien_ptr = (int *)calloc(nelem + 1, sizeof(int));
    // Homogeneous meshes are pretty easy to write
    if (lhomog)
    {
        ngnod = element[0].ngnod;
        ierr = mesh_element__getIEN(nelem, len_ien, lhomog,
                                    element, ien_ptr, connectivity);
    // Mixed element meshes require more care
    }
    else
    {
        ien = (int *)calloc(len_ien, sizeof(int));
        ierr = mesh_element__getIEN(nelem, len_ien, lhomog,
                                    element, ien_ptr, ien);
        indx = 0;
        for (ielem=0; ielem<nelem; ielem++)
        {
            ngnod = ien_ptr[ielem+1] - ien_ptr[ielem];
            if (element[ielem].type == HEX8)
            {
                connectivity[indx] = 9;
            }
            else if(element[ielem].type == TET4)
            {
                connectivity[indx] = 6;
            }
            else
            {
                printf("%s: This isn't really done yet\n", fcnm);
                connectivity[indx] = 3;
                indx = indx + 1;
                connectivity[indx] = ngnod;
            }
            indx = ien_ptr[ielem];
            for (ia=0; ia<ngnod; ia++)
            {
                connectivity[indx] = ien[indx] - ishift;
                indx = indx + 1;
            }
        }
        free(ien);
        ien = NULL;
    }
    // Free memory
    free(ien_ptr);
    ien_ptr = NULL;
    return ierr;
}
//============================================================================//
/*!
 * @brief Writes a NonLinLoc format grid file
 *
 * @param[in] dirnm         directory to write NLL model
 * @param[in] projnm        project name
 * @param[in] isp           if True then write the P velocity model.
 *                          if false then write the S velocity model
 *                          corresponding to a P or S wave model
 * @param[in] nelemx        number of elements in x
 * @param[in] nelemy        number of elements in y
 * @param[in] nelemz        number of elements in z
 * @param[in] x0            x (UTM easting) origin (m)
 * @param[in] y0            y (UTM northing) origin (m)
 * @param[in] dx            grid spacing in x (m)
 * @param[in] dy            grid spacing in y (m)
 * @param[in] dz            grid spacing in z (m)
 * @param[in] utm_zone      UTM zone in which model is defined (should be 10
 *                          for Cascadia CVM)
 * @param[in] dep0          depth origin (m).  this should be 0 because
 *                          the NLL model starts at the free surface and
 *                          extends down 
 * @param[in] iengv_ptr     points from the ielem'th element in the 
 *                          regular mesh to the start index of iengv
 *                          [nelem+1]
 * @param[in] iengv         points from the ia'th anchor node on the 
 *                          ielem'th element to the global anchor node index
 * @param[in] xmod          the Vp (m/s) or Vs (m/s) model corresponding to isp 
 *
 * @result 0 indicate success
 *
 * @author Ben Baker, ISTI
 *
 */
int meshio_write__NLLGrid(char *dirnm, char *projnm,
                          int nelemx, int nelemy, int nelemz,
                          double x0, double y0,
                          double dx, double dy, double dz, 
                          int utm_zone, double dep0,
                          struct mesh_struct mesh)
{
    const char *fcnm = "regmesh_write__NLLGrid\0";
    FILE *fpio;
    char fname[PATH_MAX], chr_type[256], file_type[64];
    const char *grid_type = "SLOW_LEN\0";
    void *buffer = NULL;
    double *xmod, lat0, lon0, vsum, vavg;
    float *slow_len = NULL;
    int ia1, ia2, ia3, ia4, ia5, ia6, ia7, ia8,
        ielem, ielemx, ielemy, ielemz, ierr, imod, jndx;
    bool isp;
    const double ngnod8i = 1.0/8.0;
    // Basic checks
    if (dx != dy || dx != dz){
        printf("%s: Error dx must equal dy must equal dz %f %f %f\n",
               fcnm, dx, dy, dz);
        return -1;
    }
    if (nelemx < 1 || nelemy < 1 || nelemz < 1){
        printf("%s: Error no points to write %d %d %d\n",
               fcnm, nelemx, nelemy, nelemz);
        return -1;
    }
    // Set space
    slow_len = (float *)calloc(nelemx*nelemy*nelemz, sizeof(float));
    buffer = (void *)calloc(nelemx*nelemy*nelemz, sizeof(float));
    // Get the lat0, lon0
    utm_geo_utm2ll(&x0, &y0, &utm_zone, &lat0, &lon0);
    if (lon0 > 180.0){lon0 = lon0 - 360.0;}
    // Loop on vp and vs model
    for (imod=0; imod<2; imod++){
        isp = true;
        if (imod == 1){isp = false;}
        // Set the filename
        memset(file_type, 0, sizeof(file_type));
        if (isp){
            strcpy(file_type, "P\0");
        }else{
            strcpy(file_type, "S\0");
        }
        ierr = meshio_write__setFilenameNLL(dirnm, projnm, file_type,
                                            true, fname);
        if (ierr != 0){
            printf("%s: Error setting buffer file name\n", fcnm);
            return -1;
        }
        if (isp){
            xmod = mesh.vp;
        }else{
            xmod = mesh.vs;
        }
        // Pack the model noting the reversed z conventions
        jndx = 0;
        for (ielemx=0; ielemx<nelemx; ielemx++){
            for (ielemy=0; ielemy<nelemy; ielemy++){
                // The model is packed bottom up - NLL is the opposite system
                for (ielemz=nelemz-1; ielemz>=0; ielemz--){
                    ielem = ielemz*nelemx*nelemy + ielemy*nelemx + ielemx;
                    ia1 = mesh.element[ielem].ien[0];
                    ia2 = mesh.element[ielem].ien[1];
                    ia3 = mesh.element[ielem].ien[2];
                    ia4 = mesh.element[ielem].ien[3];
                    ia5 = mesh.element[ielem].ien[4];
                    ia6 = mesh.element[ielem].ien[5];
                    ia7 = mesh.element[ielem].ien[6];
                    ia8 = mesh.element[ielem].ien[7];
                    vsum = xmod[ia1] + xmod[ia2] + xmod[ia3] + xmod[ia4]
                         + xmod[ia5] + xmod[ia6] + xmod[ia7] + xmod[ia8];
                    vavg = vsum*ngnod8i;
                    // note dx is units meters and vavg is units m/s
                    slow_len[jndx] = (float) (dx/vavg);
                    jndx = jndx + 1;
                }
            }
        }
        memcpy(buffer, slow_len, nelemx*nelemy*nelemz*sizeof(float));
        // Write the model
        fpio = fopen(fname, "w");
        ierr = fwrite((char *)buffer, nelemx*nelemy*nelemz*sizeof(float),
                      1, fpio);
        if (ierr != 1){
            printf("%s: Error writing model\n", fcnm);
            return -1;
        }else{
            ierr = 0;
        }
        fclose(fpio);
        // Set the type for this model (lat0, lon0, no rotation)
        memset(chr_type, 0, sizeof(chr_type));
        sprintf(chr_type, "TRANS SIMPLE %8.3f %8.3f %8.3f", lat0, lon0, 0.0);
        // Set the header filename
        ierr = meshio_write__setFilenameNLL(dirnm, projnm, file_type,
                                            false, fname);
        if (ierr != 0){
            printf("%s: Error setting header file name\n", fcnm);
            return -1;
        }
        fpio = fopen(fname, "w");
        fprintf(fpio, "%d %d %d  %lf %lf %lf  %lf %lf %lf %s\n",
                nelemx, nelemy, nelemz, x0*1.e-3, y0*1.e-3, dep0*1.e-3,
                dx*1.e-3, dy*1.e-3, dz*1.e-3, grid_type);
        fprintf(fpio, " FLOAT\n");
        fprintf(fpio, "%s\n", chr_type);
        fclose(fpio);
        xmod = NULL;
    } // Loop on vp/vs model
    free(slow_len);
    free(buffer);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the filename for outputting grid files for NonLinLoc
 */
int meshio_write__setFilenameNLL(char *dirnm, char *projnm, char *file_type,
                                 bool isbuf, char fname[PATH_MAX])
{
    const char *fcnm = "meshio_write_filenameNLL\0";
    char file_root[256];
    memset(file_root, 0, sizeof(file_root));
    if (dirnm != NULL){
        if (strlen(dirnm) > 0){
            strcpy(file_root, dirnm);
        }else{
            strcpy(file_root, ".");
        }
    }else{
        strcpy(file_root, ".");
    }
    if (file_root[strlen(file_root)-1] != '/'){
        strcat(file_root, "/");
    }
    if (!os_path_isdir(file_root)){
        printf("%s: Error output directory %s doesn't exist\n",
               fcnm, file_root);
        return -1;
    }
    strcat(file_root, projnm);
    if (file_type != NULL){
        if (isbuf){
            sprintf(fname, "%s.%s.mod.buf", file_root, file_type);
        }else{
            sprintf(fname, "%s.%s.mod.hdr", file_root, file_type);
        }
    }else{
        if (isbuf){
            sprintf(fname, "%s.mod.buf", file_root);
        }else{
            sprintf(fname, "%s.mod.hdr", file_root);
        }
    }
    return 0;
}
//============================================================================//

int meshio_write__setFilename(char *meshdir, char *projnm, char *app,
                              char fname[PATH_MAX])
{
    const char *fcnm = "meshio_write__setFilename\0";
    int ierr;
    memset(fname, 0, PATH_MAX*sizeof(char));
    if (meshdir != NULL){
        if (strlen(meshdir) == 0){
            strcpy(fname, ".\0");
        }else{
            strcpy(fname, meshdir);
            if (fname[strlen(fname)-1] == '/'){fname[strlen(fname)-1] = '\0';}
        }
    }else{
        strcpy(fname, ".\0");
    }
    strcat(fname, "/");
    // Make the directory
    if (!os_path_isdir(fname)){
        ierr = os_makedirs(fname);
        if (ierr != 0){
            printf("%s: Error making directory!\n", fcnm);
            return -1;
        }
    }
    // Append the exodus identifier
    if (strlen(projnm) > 0){strcat(fname, projnm);}
    if (strlen(projnm) > 0 && app[0] != '.'){strcat(fname, ".");}
    strcat(fname, app);
    return 0;
}
//============================================================================//
/*!
 * @brief Writes the mesh in an HDF5 file format with an accompanying
 *        XDMF file for viewing
 *
 * @param[in] meshdir    name of directory where the mesh is to reside
 * @param[in] projnm     name of project
 * @param[in] mesh       holds the mesh pointers and elements to write
 *
 * @result 0 indicates success
 *
 */
int meshio_write__h5(char *meshdir, char *projnm,
                     struct mesh_struct mesh)
{
    const char *fcnm = "meshio_write__h5\0";
    enum mesh_element_type type;
    const char *chex8 = "HEX8\0";
    char exofl[PATH_MAX];
    double *coords, *xmat;
    int *connectivity, ierr, ncon, ngnod;
    bool cnum = true;  // We use C based numbering
    const int nmat = 5;
    hsize_t dims[2];
    hid_t file_id, dataset_id, dataspace_id; 
    herr_t status;
    //------------------------------------------------------------------------//
    ierr = 0;
    status = 0;
    coords = NULL;
    connectivity = NULL;
    xmat = NULL;
    // Verify there's something in the mesh
    if (mesh.nelem < 1 || mesh.nnpg < 1 || mesh.element == NULL){
        if (mesh.nelem < 1){printf("%s: No elements!\n", fcnm);}
        if (mesh.nnpg < 1){printf("%s: No anchor nodes!\n", fcnm);}
        if (mesh.element == NULL){printf("%s: Element is null!\n", fcnm);}
        return -1;
    }
    // Set the mesh file
    ierr = meshio_write__setFilename(meshdir, projnm, ".h5\0", exofl);
    if (ierr != 0){
        printf("%s: Failed to set mesh filename!\n", fcnm);
        return -1;
    }
    // Make the hdf5 file
    file_id = H5Fcreate(exofl, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //-------------------------------coordinates------------------------------//
    if (!mesh.lptr_only){
        coords = (double *)calloc(3*mesh.nnpg, sizeof(double));
        ierr = mesh_element__getAnchorNodeLocations(mesh.nelem, mesh.nnpg,
                                                    mesh.element,
                                                    &coords[0],
                                                    &coords[mesh.nnpg],
                                                    &coords[2*mesh.nnpg]);
        if (ierr != 0){
            printf("%s: Error getting nodal locations!\n", fcnm);
            goto ERROR;
        }
        h5_write_array__double("XCoordinates\0", file_id, mesh.nnpg,
                               &coords[0]);
        h5_write_array__double("YCoordinates\0", file_id, mesh.nnpg,
                               &coords[mesh.nnpg]);
        h5_write_array__double("ZCoordinates\0", file_id, mesh.nnpg,
                               &coords[2*mesh.nnpg]);
        free(coords);
        coords = NULL;
    }else{
        h5_write_array__double("XCoordinates\0", file_id, mesh.nnpg,
                               mesh.xlocs);
        h5_write_array__double("YCoordinates\0", file_id, mesh.nnpg,
                               mesh.ylocs);
        h5_write_array__double("ZCoordinates\0", file_id, mesh.nnpg,
                               mesh.zlocs);
    }
    //----------------------------material properties-------------------------//
    if (!mesh.lptr_only){
        xmat = (double *)calloc(nmat*mesh.nnpg, sizeof(double));
        ierr = mesh_element__getAnchorNodeProperties(mesh.nelem, mesh.nnpg,
                                                     mesh.element,
                                                     &xmat[0],           // vp
                                                     &xmat[mesh.nnpg],   // vs
                                                     &xmat[2*mesh.nnpg], // dens
                                                     &xmat[3*mesh.nnpg], // Qp
                                                     &xmat[4*mesh.nnpg]);// Qs 
        if (ierr != 0){
            printf("%s: Error getting material properties!\n", fcnm);
            goto ERROR;
        }
        h5_write_array__double("CompressionalVelocity\0", file_id, mesh.nnpg,
                               &xmat[0]); 
        h5_write_array__double("ShearVelocity\0", file_id, mesh.nnpg,
                               &xmat[mesh.nnpg]);
        h5_write_array__double("Density\0", file_id, mesh.nnpg,
                               &xmat[2*mesh.nnpg]);
        h5_write_array__double("Qp\0", file_id, mesh.nnpg,
                               &xmat[3*mesh.nnpg]);
        h5_write_array__double("Qs\0", file_id, mesh.nnpg,
                               &xmat[4*mesh.nnpg]);
        free(xmat);
        xmat = NULL;
    }else{
        h5_write_array__double("CompressionalVelocity\0", file_id, mesh.nnpg,
                               mesh.vp);
        h5_write_array__double("ShearVelocity\0", file_id, mesh.nnpg,
                               mesh.vs);
        h5_write_array__double("Density\0", file_id, mesh.nnpg,
                               mesh.dens);
        h5_write_array__double("Qp\0", file_id, mesh.nnpg,
                               mesh.Qp);
        h5_write_array__double("Qs\0", file_id, mesh.nnpg,
                               mesh.Qs);
    }
    //-------------------------------connectivity-----------------------------//
    // Get the connectivity
    if (!mesh.lhomog){
        printf("%s: Cannot handle mixed element types!\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    ncon = meshio_mesh2connectivity__getSize(mesh.lhomog, mesh.nelem,
                                             mesh.element);
    if (ncon < 1){
        printf("%s: Error no points in connectivity\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    connectivity = (int *)calloc(ncon, sizeof(int));
    ierr = meshio_mesh2connectivity(cnum, mesh.lhomog, mesh.nelem,
                                    mesh.element,
                                    ncon, connectivity);
    if (ierr != 0){
        printf("%s: Error computing connectivity\n", fcnm);
        goto ERROR;
    }
    type = mesh.element[0].type;
    if (type != HEX8 && type != TET4){
        printf("%s: Warning not sure about this\n", fcnm);
    }
    ngnod = mesh.element[0].ngnod;
    dims[0] = mesh.nelem;
    dims[1] = ngnod;
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, "Connectivity\0", H5T_NATIVE_INT,
                            dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status += H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, connectivity);
    status += h5_write_attribute__char("elem_type\0", dataset_id,
                                       1, &chex8);
    status += H5Sclose(dataspace_id);
    status += H5Dclose(dataset_id);
    if (status != 0){
        printf("%s: Error writing connectivity!\n", fcnm);
        goto ERROR;
    }
    // Free space
ERROR:;
    if (coords != NULL){
        free(coords);
        coords = NULL;
    }
    if (xmat != NULL){
        free(xmat);
        xmat = NULL;
    }
    // Close the hdf5 file
    status = H5Fclose(file_id);
    // Write the Xdmf file
    if (ierr == 0){
        meshio_write__xdmf(meshdir, projnm, mesh.nelem, ngnod, mesh.nnpg);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Writes the mesh for consumption by SPECFEM3D
 *
 * @param[in] meshdir     name of directory to write mesh
 * @param[in] projnm      name of project
 *
 * @param[in] mesh        holds the mesh structure
 *
 * @result 0 indicates success
 *
 */
int meshio_write__specfem3d(char *meshdir, char *projnm,
                            struct mesh_struct mesh)
{
    const char *fcnm = "meshio_write__specfem3d\0";
    FILE *ofl;
    char app[64];
    enum mesh_bc_enum side;
    double *xlocs, *ylocs, *zlocs, vp_avg, vs_avg, Qp_avg, Qs_avg, dn_avg;
    char meshfl[PATH_MAX];
    int *bdry2glob_elem, *ien_bdry, *ien_bdry_ptr,
        ia, ielem, ielem_bdry, ierr, inpg, iside, lenien_bdry, nelem_bdry;
    const int elastic = 2;    // this is an elastic simulation
    const int anisotropy = 0; // no anisotropic models
    const double Q_kappa = 9999.0; // turn off p-wave attenuation
    const double Q_mu = 9999.0;    // turn off s-wave attenuation
    //------------------------------------------------------------------------//
    //
    // Initialize
    ierr = 0;
    xlocs = NULL;
    ylocs = NULL;
    zlocs = NULL;
    // Verify there's something in the mesh
    if (mesh.nelem < 1 || mesh.nnpg < 1 || mesh.element == NULL){
        if (mesh.nelem < 1){printf("%s: No elements!\n", fcnm);}
        if (mesh.nnpg < 1){printf("%s: No anchor nodes!\n", fcnm);}
        if (mesh.element == NULL){printf("%s: Element is null!\n", fcnm);}
        return -1;
    }
    //---------------------------write the mesh file--------------------------//
    ierr = meshio_write__setFilename(meshdir, "\0", // projnm
                                     "mesh_file\0", meshfl);
    if (ierr != 0){
        printf("%s: Error setting mesh file!\n", fcnm);
        return -1;
    }
    ofl = fopen(meshfl, "w");
    fprintf(ofl, "%d\n", mesh.nelem); 
    for (ielem=0; ielem<mesh.nelem; ielem++){
        if (mesh.element[ielem].type != HEX8 &&
            mesh.element[ielem].type != HEX27)
        {
            printf("%s: Error SPECFEM3D only accepts hex meshes\n", fcnm);
            ierr = 1;
            break;
        }
        fprintf(ofl, "%d", ielem + 1); // C->Fortran
        for (ia=0; ia<mesh.element[ielem].ngnod; ia++){
            fprintf(ofl, " %d", mesh.element[ielem].ien[ia] + 1);// C->Fortran
        }
        fprintf(ofl, "\n");
    }
    fclose(ofl);
    //-----------------------write the nodes coords file----------------------//
    ierr = meshio_write__setFilename(meshdir, "\0", //projnm,
                                     "nodes_coords_file\0", meshfl);
    if (ierr != 0){
        printf("%s: Error setting node coordinates file!\n", fcnm);
        return -1;
    }
    ofl = fopen(meshfl, "w");
    fprintf(ofl, "%d\n", mesh.nnpg);
    if (mesh.lptr_only){
        for (inpg=0; inpg<mesh.nnpg; inpg++){
            fprintf(ofl, "%12d  %18.6f %18.6f %18.6f\n",
                    inpg + 1, // c -> fortran
                    mesh.xlocs[inpg], mesh.ylocs[inpg], mesh.zlocs[inpg]);
        }
    }else{
        xlocs = (double *)calloc(mesh.nnpg, sizeof(double));
        ylocs = (double *)calloc(mesh.nnpg, sizeof(double));
        zlocs = (double *)calloc(mesh.nnpg, sizeof(double));
        ierr = mesh_element__getAnchorNodeLocations(mesh.nelem, mesh.nnpg,
                                                    mesh.element,
                                                    xlocs, ylocs, zlocs);
        
        for (inpg=0; inpg<mesh.nnpg; inpg++){
            fprintf(ofl, "%12d  %18.6f %18.6f %18.6f\n",
                    inpg + 1, // c -> fortran
                    xlocs[inpg], ylocs[inpg], zlocs[inpg]);
        }
        free(xlocs);
        free(ylocs);
        free(zlocs);
        xlocs = NULL;
        ylocs = NULL;
        zlocs = NULL;
    }
    fclose(ofl);
    //---------------------------write the nummaterial file-------------------//
    ierr = meshio_write__setFilename(meshdir, "\0", //projnm,
                                     "nummaterial_velocity_file\0", meshfl);
    if (ierr != 0){
        printf("%s: Error setting nummaterial file!\n", fcnm);
        return -1;
    }
    ofl = fopen(meshfl, "w");
    if (mesh.lptr_only){
        for (ielem=0; ielem<mesh.nelem; ielem++){
            vp_avg = 0.0;
            vs_avg = 0.0;
            dn_avg = 0.0;
            Qp_avg = 0.0;
            Qs_avg = 0.0;
            for (ia=0; ia<mesh.element[ielem].ngnod; ia++){
                inpg = mesh.element[ielem].ien[ia];
                vp_avg = vp_avg + mesh.vp[inpg];
                vs_avg = vs_avg + mesh.vs[inpg];
                dn_avg = dn_avg + mesh.dens[inpg];
                Qp_avg = Qp_avg + mesh.Qp[inpg];
                Qs_avg = Qs_avg + mesh.Qs[inpg]; 
            }
            vp_avg = vp_avg/(double) mesh.element[ielem].ngnod;
            vs_avg = vs_avg/(double) mesh.element[ielem].ngnod;
            dn_avg = dn_avg/(double) mesh.element[ielem].ngnod;
            Qp_avg = Qp_avg/(double) mesh.element[ielem].ngnod;
            Qs_avg = Qs_avg/(double) mesh.element[ielem].ngnod;
            if (Qp_avg <= 0.0){Qp_avg = Q_kappa;}
            if (Qs_avg <= 0.0){Qs_avg = Q_mu;}
            fprintf(ofl, "%d %12d %18.6f %18.6f %18.6f %18.6f %18.6f %d\n",
                    elastic, ielem+1,
                    dn_avg, vp_avg, vs_avg, Qp_avg, Qs_avg,
                    anisotropy);
        }
    }else{
        for (ielem=0; ielem<mesh.nelem; ielem++){
            vp_avg = 0.0;
            vs_avg = 0.0;
            dn_avg = 0.0;
            Qp_avg = 0.0;
            Qs_avg = 0.0;
            for (ia=0; ia<mesh.element[ielem].ngnod; ia++){
                vp_avg = vp_avg + mesh.element[ielem].vp[ia];
                vs_avg = vs_avg + mesh.element[ielem].vs[ia];
                dn_avg = dn_avg + mesh.element[ielem].dens[ia];
                Qp_avg = Qp_avg + mesh.element[ielem].Qp[ia];
                Qs_avg = Qs_avg + mesh.element[ielem].Qs[ia];
            }
            vp_avg = vp_avg/(double) mesh.element[ielem].ngnod;
            vs_avg = vs_avg/(double) mesh.element[ielem].ngnod;
            dn_avg = dn_avg/(double) mesh.element[ielem].ngnod;
            Qp_avg = Qp_avg/(double) mesh.element[ielem].ngnod;
            Qs_avg = Qs_avg/(double) mesh.element[ielem].ngnod;
            if (Qp_avg <= 0.0){Qp_avg = Q_kappa;}
            if (Qs_avg <= 0.0){Qs_avg = Q_mu;}
            fprintf(ofl, "%1d %12d %18.6f %18.6f %18.6f %18.6f %18.6f %d\n",
                    elastic, ielem+1,
                    dn_avg, vp_avg, vs_avg, Qp_avg, Qs_avg,
                    anisotropy);
        }
    }
    fclose(ofl);
    //-------------------------------materials file---------------------------//
    ierr = meshio_write__setFilename(meshdir, "\0", //projnm,
                                     "materials_file\0", meshfl);
    if (ierr != 0){
        printf("%s: Error setting materials file!\n", fcnm);
        return -1;
    }
    ofl = fopen(meshfl, "w");
    for (ielem=0; ielem<mesh.nelem; ielem++){
        fprintf(ofl, "%d %d\n", ielem+1, ielem+1);
    }
    fclose(ofl);
    //----------------------------boundary conditions-------------------------//
    for (iside=0; iside<6; iside++){
        memset(app, 0, sizeof(app));
        if (iside == 0)
        {
            strcpy(app, "absorbing_surface_file_ymin\0");
            side = SOUTH_BDRY;
        }
        else if (iside == 1)
        {
            strcpy(app, "absorbing_surface_file_xmin\0");
            side = WEST_BDRY;
        }
        else if (iside == 2)
        {
            strcpy(app, "absorbing_surface_file_ymax\0");
            side = NORTH_BDRY;
        }
        else if (iside == 3)
        {
            strcpy(app, "absorbing_surface_file_xmax\0");
            side = EAST_BDRY;
        }
        else if (iside == 4){
            strcpy(app, "absorbing_surface_file_bottom\0");
            side = BOTTOM_BDRY;
        }
        else if (iside == 5)
        {
            strcpy(app, "free_or_absorbing_surface_file_zmax\0");
            side = TOP_BDRY;
        }else{
            printf("%s: Unknown side\n", fcnm);
            continue;
        }
        nelem_bdry = mesh_element__getNumberOfBoundaryElements(mesh.nelem,
                                                               side, 
                                                               mesh.element);
        if (nelem_bdry == 0){continue;}
        if (!mesh.lhomog){
            lenien_bdry = mesh_element__getBoundaryIENSize(mesh.nelem,
                                                           side,
                                                           mesh.element);
        }else{
            lenien_bdry = nelem_bdry*mesh.element[0].ngnod_face;
        }
        if (lenien_bdry == 0){
            printf("%s: Error getting ien size!\n", fcnm);
            continue;
        }
        bdry2glob_elem = (int *)calloc(nelem_bdry, sizeof(int)); 
        ien_bdry = (int *)calloc(lenien_bdry, sizeof(int));
        ien_bdry_ptr = (int *)calloc(nelem_bdry + 1, sizeof(int));
        ierr = mesh_element__getIENBoundary(mesh.nelem, mesh.lhomog,
                                            side,
                                            mesh.element,
                                            nelem_bdry, ien_bdry_ptr,
                                            bdry2glob_elem,
                                            lenien_bdry, ien_bdry);
        if (ierr != 0){
            printf("%s: Error setting IENboundary\n", fcnm);
            goto ERROR;
        }
        // Open file for writing
        ierr = meshio_write__setFilename(meshdir, "\0", //projnm,
                                         app, meshfl);
        if (ierr != 0){ 
            printf("%s: Error setting boundary file!\n", fcnm);
            goto ERROR;
        }   
        ofl = fopen(meshfl, "w");
        fprintf(ofl, "%d\n", nelem_bdry);
        for (ielem_bdry=0; ielem_bdry<nelem_bdry; ielem_bdry++){
            fprintf(ofl, "%d", bdry2glob_elem[ielem_bdry] + 1); //c -> fortran
            for (ia=ien_bdry_ptr[ielem_bdry];
                 ia<ien_bdry_ptr[ielem_bdry+1]; ia++){
                fprintf(ofl, " %d", ien_bdry[ia] + 1); //c -> fortran
            }
            fprintf(ofl, "\n");
        }
        fclose(ofl);
ERROR:;
        free(bdry2glob_elem);
        free(ien_bdry_ptr);
        free(ien_bdry);
        bdry2glob_elem = NULL;
        ien_bdry_ptr = NULL;
        ien_bdry = NULL;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Writes the Xdmf file accompanying the H5 mesh dump for viewing with
 *        Paraview
 *
 * @param[in] meshdir    directory containing the mesh
 * @param[in] projnm     name of project
 * @param[in] nelem      number of elements in mesh
 * @param[in] ngnod      number of anchor nodes on each element
 * @param[in] nnpg       number of anchor nodes in mesh
 *
 * @author Ben Baker, ISTI
 */
int meshio_write__xdmf(char *meshdir, char *projnm,
                       int nelem, int ngnod, int nnpg)
{
    FILE *xfl;
    const char *fcnm = "meshio_write__xdmf\0";
    char xdmfl[PATH_MAX], h5fl[PATH_MAX], shape[64];
    int ierr;
    // Set the mesh file
    ierr = meshio_write__setFilename(meshdir, projnm, ".xdmf\0", xdmfl);
    ierr = meshio_write__setFilename("./\0", projnm, ".h5\0", h5fl);
    if (ierr != 0){ 
        printf("%s: Failed to set mesh filename!\n", fcnm);
        return -1; 
    } 
    memset(shape, 0, sizeof(shape));
    if (ngnod == 8){
        strcpy(shape, "Hexahedron\0"); 
    }else if (ngnod == 4){
        strcpy(shape, "Tetrahedron\0");
    }
    // Open the file
    xfl = fopen(xdmfl, "w");
    fprintf(xfl, "<?xml version=\"1.0\" ?>\n");
    fprintf(xfl, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xfl, "<Xdmf Version=\"2.0\">\n");
    fprintf(xfl, "  <Domain>\n");
    fprintf(xfl, "    <Grid Collection=\"whole\" Name=\"Mesh\">\n");
    fprintf(xfl, "      <Topology TopologyType=\"%s\" NumberOfElements=\"%d\">\n",
            shape, nelem);
    fprintf(xfl, "        <DataItem Dimensions=\"%d %d\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
            nelem, ngnod);
    fprintf(xfl, "           %s:Connectivity\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Topology>\n");
    fprintf(xfl, "      <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg); 
    fprintf(xfl, "           %s:XCoordinates\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg); 
    fprintf(xfl, "           %s:YCoordinates\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg); 
    fprintf(xfl, "           %s:ZCoordinates\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Geometry>\n");
    fprintf(xfl, "      <Attribute Name=\"Compressional Velocity (m/s)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "           %s:CompressionalVelocity\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "      <Attribute Name=\"Shear Velocity (m/s)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "           %s:ShearVelocity\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "      <Attribute Name=\"Density (kg/m**3)\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "           %s:Density\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "      <Attribute Name=\"Qp\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "           %s:Qp\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "      <Attribute Name=\"Qs\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "           %s:Qs\n", h5fl);
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "      <Attribute Name=\"VpVs Ratio\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xfl, "        <DataItem ItemType=\"Function\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\" Function=\"$0 / $1\">\n",
            nnpg);

    fprintf(xfl, "          <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl, "             %s:CompressionalVelocity\n", h5fl);
    fprintf(xfl, "          </DataItem>\n");
    fprintf(xfl, "          <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
            nnpg);
    fprintf(xfl,   "           %s:ShearVelocity\n", h5fl);
    fprintf(xfl, "          </DataItem>\n");
    fprintf(xfl, "        </DataItem>\n");
    fprintf(xfl, "      </Attribute>\n");
    fprintf(xfl, "    </Grid>\n");
    fprintf(xfl, "  </Domain>\n");
    fprintf(xfl, "</Xdmf>\n");
    fclose(xfl);
    return 0;
}
