#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "mesh.h"

/*!
 * @brief Triples a hexahedron element in the vertical direction with 
 *        the top of the element being fined and the bottom of the
 *        element being coarse
 *
 * @param[out] tmplate    returns the 13 element tripling template
 *                        defined on the unit box [-1,1] x [-1,1] x [-1,1]
 *                        with the most refined elements at z = 1
 *                        and the least refined element at z =-1 [13]
 *
 * @author Ben Baker, ISTI
 *
 * @reference http://www.robertschneiders.de/papers/vki.pdf
 */
void mesh_template__hexTriple(struct mesh_element_struct *tmplate)
{
     const double xt[32] = {   1.000000000,  -1.000000000,   1.000000000,
                              -1.000000000,   1.000000000,   1.000000000,
                              -1.000000000,  -1.000000000,  -1.000000000,
                              -1.000000000,   1.000000000,   1.000000000,
                              -0.333333333,   0.333333333,   0.333333333,
                              -0.333333333,   0.333333333,  -0.333333333,
                               0.333333333,  -0.333333333,  -0.333333333,
                               0.333333333,   1.000000000,   1.000000000,
                              -1.000000000,  -1.000000000,  -0.333333333,
                               0.333333333,  -0.333333333,   0.333333333,
                              -0.333333333,   0.333333333};
     const double yt[32] = {   1.000000000,   1.000000000,  -1.000000000,
                              -1.000000000,   1.000000000,  -1.000000000,
                              -1.000000000,   1.000000000,   0.333333333,
                              -0.333333333,  -0.333333333,   0.333333333,
                              -1.000000000,  -1.000000000,   1.000000000,
                               1.000000000,  -0.333333333,  -0.333333333,
                               0.333333333,   0.333333333,   1.000000000,
                               1.000000000,  -0.333333333,   0.333333333,
                              -0.333333333,   0.333333333,   0.333333333,
                               0.333333333,  -0.333333333,  -0.333333333,
                              -1.000000000,  -1.000000000};
     const double zt[32] = {  -1.000000000,  -1.000000000,  -1.000000000,
                              -1.000000000,   1.000000000,   1.000000000,
                               1.000000000,   1.000000000,   1.000000000,
                               1.000000000,   1.000000000,   1.000000000,
                               1.000000000,   1.000000000,   1.000000000,
                               1.000000000,   1.000000000,   1.000000000,
                               1.000000000,   1.000000000,   0.000000000,
                               0.000000000,   0.000000000,   0.000000000,
                               0.000000000,   0.000000000,   0.500000000,
                               0.500000000,   0.500000000,   0.500000000,
                               0.000000000,   0.000000000};
     const int ien_face[312] = {24, 22, 23, 25, 
                                22,  2,  0, 23, 
                                 2,  3,  1,  0, 
                                24, 25,  1,  3, 
                                24,  3,  2, 22, 
                                25, 23,  0,  1, 
                                 7,  8, 19, 15, 
                                 8, 25, 26, 19, 
                                25,  1, 20, 26, 
                                 7, 15, 20,  1, 
                                 7,  1, 25,  8, 
                                15, 19, 26, 20, 
                                15, 19, 18, 14, 
                                19, 26, 27, 18, 
                                26, 20, 21, 27, 
                                15, 14, 21, 20, 
                                15, 20, 26, 19, 
                                14, 18, 27, 21, 
                                14, 18, 11,  4, 
                                18, 27, 23, 11, 
                                27, 21,  0, 23, 
                                14,  4,  0, 21, 
                                14, 21, 27, 18, 
                                 4, 11, 23,  0, 
                                20, 26, 27, 21, 
                                26, 25, 23, 27, 
                                25,  1,  0, 23, 
                                20, 21,  0,  1, 
                                20,  1, 25, 26, 
                                21, 27, 23,  0, 
                                 8,  9, 17, 19, 
                                 9, 24, 28, 17, 
                                24, 25, 26, 28, 
                                 8, 19, 26, 25, 
                                 8, 25, 24,  9, 
                                19, 17, 28, 26, 
                                19, 17, 16, 18, 
                                17, 28, 29, 16, 
                                28, 26, 27, 29, 
                                19, 18, 27, 26, 
                                19, 26, 28, 17, 
                                18, 16, 29, 27, 
                                18, 16, 10, 11, 
                                16, 29, 22, 10, 
                                29, 27, 23, 22, 
                                18, 11, 23, 27, 
                                18, 27, 29, 16, 
                                11, 10, 22, 23, 
                                26, 28, 29, 27, 
                                28, 24, 22, 29, 
                                24, 25, 23, 22, 
                                26, 27, 23, 25, 
                                26, 25, 24, 28, 
                                27, 29, 22, 23, 
                                 9,  6, 12, 17, 
                                 6,  3, 30, 12, 
                                 3, 24, 28, 30, 
                                 9, 17, 28, 24, 
                                 9, 24,  3,  6, 
                                17, 12, 30, 28, 
                                17, 12, 13, 16, 
                                12, 30, 31, 13, 
                                30, 28, 29, 31, 
                                17, 16, 29, 28, 
                                17, 28, 30, 12, 
                                16, 13, 31, 29, 
                                16, 13,  5, 10, 
                                13, 31,  2,  5, 
                                31, 29, 22,  2, 
                                16, 10, 22, 29, 
                                16, 29, 31, 13, 
                                10,  5,  2, 22, 
                                28, 30, 31, 29, 
                                30,  3,  2, 31, 
                                 3, 24, 22,  2, 
                                28, 29, 22, 24, 
                                28, 24,  3, 30, 
                                29, 31,  2, 22};
     const int ien[104] = {24, 22,  2,  3, 25, 23,  0,  1, 
                            7,  8, 25,  1, 15, 19, 26, 20, 
                           15, 19, 26, 20, 14, 18, 27, 21, 
                           14, 18, 27, 21,  4, 11, 23,  0, 
                           20, 26, 25,  1, 21, 27, 23,  0, 
                            8,  9, 24, 25, 19, 17, 28, 26, 
                           19, 17, 28, 26, 18, 16, 29, 27, 
                           18, 16, 29, 27, 11, 10, 22, 23, 
                           26, 28, 24, 25, 27, 29, 22, 23, 
                            9,  6,  3, 24, 17, 12, 30, 28, 
                           17, 12, 30, 28, 16, 13, 31, 29, 
                           16, 13, 31, 29, 10,  5,  2, 22, 
                           28, 30,  3, 24, 29, 31,  2, 22};
     const int bc[78] = { 0,  2,  5,  4,  0,  0, 
                          6,  0,  0,  3,  4,  0, 
                          6,  0,  0,  3,  0,  0, 
                          6,  0,  0,  3,  0,  2, 
                          0,  0,  0,  3,  0,  0, 
                          6,  0,  0,  0,  4,  0, 
                          6,  0,  0,  0,  0,  0, 
                          6,  0,  0,  0,  0,  2, 
                          0,  0,  0,  0,  0,  0, 
                          6,  1,  0,  0,  4,  0, 
                          6,  1,  0,  0,  0,  0, 
                          6,  1,  0,  0,  0,  2, 
                          0,  1,  0,  0,  0,  0};
     const int neighbor[78] = { 8, -1, -1, -1, 12,  4, 
                               -1,  5,  4, -1, -1,  2, 
                               -1,  6,  4, -1,  1,  3, 
                               -1,  7,  4, -1,  2, -1, 
                                2,  8,  0, -1,  1,  3, 
                               -1,  9,  8,  1, -1,  6, 
                               -1, 10,  8,  2,  5,  7, 
                               -1, 11,  8,  3,  6, -1, 
                                6, 12,  0,  4,  5,  7, 
                               -1, -1, 12,  5, -1, 10, 
                               -1, -1, 12,  6,  9, 11, 
                               -1, -1, 12,  7, 10, -1, 
                               10, -1,  0,  8,  9, 11};
     const int ngnod = 8;
     const int nface = 6;
     const int ngnod_face = 4;
     const int elem_type = HEX8;
     const int nelem_template = 13;
     int ia, ielem, iface, indx, inpg, jndx;
     memset(tmplate, 0, nelem_template*sizeof(tmplate));
     for (ielem=0; ielem<nelem_template; ielem++)
     {
         tmplate[ielem].x = (double *)calloc(ngnod, sizeof(double));
         tmplate[ielem].y = (double *)calloc(ngnod, sizeof(double));
         tmplate[ielem].z = (double *)calloc(ngnod, sizeof(double));
         tmplate[ielem].ien_face = (int *)calloc(ngnod_face*nface,
                                                 sizeof(int));
         tmplate[ielem].ien = (int *)calloc(ngnod, sizeof(int));
         tmplate[ielem].bc = (int *)calloc(nface, sizeof(int));
         tmplate[ielem].neighbor = (int *)calloc(nface, sizeof(int));
         for (ia=0; ia<8; ia++)
         {
             inpg = ien[8*ielem+ia];
             tmplate[ielem].ien[ia] = inpg;
             tmplate[ielem].x[ia] = xt[inpg];
             tmplate[ielem].y[ia] = yt[inpg];
             tmplate[ielem].z[ia] = zt[inpg];
         }
         for (iface=0; iface<nface; iface++)
         {
             tmplate[ielem].bc[iface] = bc[6*ielem+iface];
             tmplate[ielem].neighbor[iface] = neighbor[6*ielem+iface];
             for (ia=0; ia<ngnod_face; ia++)
             {
                 indx = iface*ngnod_face + ia;
                 jndx = ielem*ngnod_face*nface + iface*ngnod_face + ia;
                 tmplate[ielem].ien_face[indx] = ien_face[jndx];
             }
         }
         tmplate[ielem].nface = nface;
         tmplate[ielem].ngnod = ngnod;
         tmplate[ielem].ngnod_face = ngnod_face;
         tmplate[ielem].type = elem_type;
     }
     return;
}
/*!
 * @brief Decomposes a hex into five tetrahedra
 *
 * @author Ben Baker, ISTI
 *
 * @reference https://ieeexplore.ieee.org/ieee_pilot/articles/06/ttg2009061587/figures.html
 *
 */
/*
void mesh_template__hex2tet5()//struct element_struct hex, struct element_struct *tet)
{
    const int face_mask[12] = {0, 1, 3,    // Face 1
                               1, 2, 3,    // Face 2
                               0, 3, 2,    // Face 3 
                               0, 2, 1};   // Face 4
    const int tet_mask[20] = {7, 5, 4, 0,  // Tet 1
                              7, 6, 5, 2,  // Tet 2
                              0, 7, 2, 3,  // Tet 3
                              0, 7, 5, 2,  // Tet 4
                              0, 1, 2, 5}; // Tet 5
    const int nelem_tet = 5;
    int ielem;
    for (ielem=0; ielem<nelem_tet; ielem++){

    }
    return;
}
*/
//============================================================================//
/*!
 * @brief Computes the volume of a 4 node tetrahedron
 *
 * @param[in] x      x anchor node locations of tetrahedron [4] 
 * @param[in] y      y anchor node locations of tetrahedron [4]
 * @param[in] z      z anchor node locations of tetrahedron [4]
 *
 * @result the volume of the 4 noded tetrahedron
 *
 * @reference http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
 *
 * @author Ben Baker, ISTI
 *
 */
#pragma omp declare simd
double mesh_volume__tet4(double x[4], double y[4], double z[4])
{
    double J, volume, x21, x32, x43, y12, y23, y34, z12, z23, z34;
    const double sixth = 1.0/6.0;
    // x vectors
    x21 = x[1] - x[0];
    x32 = x[2] - x[1];
    x43 = x[3] - x[2];
    // y vectors
    y12 = y[0] - y[1];
    y23 = y[1] - y[2];
    y34 = y[2] - y[3];
    // z vectors
    z12 = z[0] - z[1];
    z23 = z[1] - z[2];
    z34 = z[2] - z[3];
    // Compute jacobian  
    J = x21*(y23*z34 - y34*z23)
      + x32*(y34*z12 - y12*z34)
      + x43*(y12*z23 - y23*z12);
    // V = J/6
    volume = sixth*J;
    return volume;
}
//============================================================================//
/*!
 * @brief Computes the volume of an irregular hexahedron by decomposing
 *        it into 5 tetrahedron
 *
 * @param[in] x       x anchor node locations of hexahedron [8]
 * @param[in] y       y anchor node locations of hexahedron [8]
 * @param[in] z       z anchor node locations of hexahedron [8]
 *
 * @result volume of the 8 noded hexahedron 
 *
 * @author Ben Baker, ISTI
 *
 */
double mesh_volume__hex8(double x[8], double y[8], double z[8])
{
    double x4[4], y4[4], z4[4], volume;
    int ia, ielem, indx, jndx;
    const int tet_mask[20] = {7, 5, 4, 0,  // Tet 1
                              7, 6, 5, 2,  // Tet 2
                              0, 7, 2, 3,  // Tet 3
                              0, 7, 5, 2,  // Tet 4
                              0, 1, 2, 5}; // Tet 5
    const int nelem_tet = 5;
    volume = 0.0;
    #pragma omp simd reduction(+:volume)
    for (ielem=0; ielem<nelem_tet; ielem++){
        jndx = 4*ielem;
        for (ia=0; ia<4; ia++){
            indx = tet_mask[jndx];
            x4[ia] = x[indx];
            y4[ia] = y[indx];
            z4[ia] = z[indx];
            jndx = jndx + 1;
        }
        volume = volume + mesh_volume__tet4(x4, y4, z4);
    }
    return volume;
}

/*
int main()
{
    double x[8], y[8], z[8], J;
    x[0] = 0.0; y[0] = 0.0; z[0] = 0.0;
    x[1] = 4.0; y[1] = 0.0; z[1] = 0.0;
    x[2] = 4.0; y[2] = 2.0; z[2] = 0.0;
    x[3] = 0.0; y[3] = 2.0; z[3] = 0.0;
    x[4] = 0.0; y[4] = 0.0; z[4] = 3.0;
    x[5] = 4.0; y[5] = 0.0; z[5] = 3.0;
    x[6] = 4.0; y[6] = 2.0; z[6] = 3.0;
    x[7] = 0.0; y[7] = 2.0; z[7] = 3.0;
int i;
for (i=0; i<1000000; i++){
    J = mesh_volume__hex8(x, y, z);
}
printf("%f\n",J);
    return 0;
}
*/
