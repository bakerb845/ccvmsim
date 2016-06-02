#ifndef __CVM_CONSTANTS_H__
#define __CVM_CONSTANTS_H__
// These are the constants from Art; the model is composed of 3 distinct
// layers:  The first, shallowest layer is from [0,1.2] km, the second 
// intermediate layer is from [1.5,9.9] km, and deepest layer is from 
// [10.8,59.4] km
const int nlay_cvm = 3;
const double dxl_cvm[3] = {200.0, 300.0, 900.0}; /*!< Grid spacing (m) in x in
                                                      each layer */ 
const double dyl_cvm[3] = {200.0, 300.0, 900.0}; /*!< Grid spacing (m) in y in
                                                      each layer */
const double dzl_cvm[3] = {100.0, 300.0, 900.0}; /*!< Grid spacing (m) in z in
                                                      each layer */
const double z0_cvm[4]  = {0.0, 1500.0, 10800.0, 59400.0}; /*!< Depth (m) of each
                                                               layer top */
const int nxl_cvm[3] = {3271, 2181,  727};  /*!< Number of x grid points in
                                                 each layer */
const double common_mult = 1800.0;    /*!< Common multiple grid spacing
                                           for all layers in x and y */
//const int nxl_cvm_use[3] = {3267, 2178, 726}; /*!< Number of usable points so the interpolation works */
const int nyl_cvm[3] = {5367, 3578, 1193};  /*!< Number of y grid points in
                                                 each layer */
//const int nyl_cvm_use[3] = {};
const int nzl_cvm[3] = {  13,   29,   55};  /*!< Number of z grid points in
                                                 each layer */
const int utmzone_cvm = 10;    /*!< CVM is defined in UTM zone 10 */
// Note these constants are from the model definition - they will be 
// recomputed by the program as there exists a maximum and minimum usable
// allowable geometry as the model is specified on a rectangular grid which
// is too large - look at Figure 1 of Velocity and Density Incorporating
// the Cascadia Subduction Zone for 3D Earthquake Ground Motion Simulations
// for an understanding the maxes and mins computed by the program will look
// like
const double utm_x0_cvm =-10800.0;   /*!< 40.2 N - lower left corner (m) */
const double utm_y0_cvm = 4467300.0; /*!< -129 W - lower left corner (m) */
const double lat0_cvm = 40.2;        /*!< CVM lower left corner (degrees) */
const double lon0_cvm = 231.0;       /*!< CVM lower left corner (degrees) */
const double lat0_cvm_usable = 40.344076;  /*!< Min usable lat (deg) */
const double lon0_cvm_usable = 231.000000; /*<! Min usable lon (deg) */
const double lat1_cvm_usable = 49.796694;  /*<! Max usable latitude (deg) */
const double lon1_cvm_usable = 238.989666; /*<! Max usable longitude (deg) */
//const double lat0_cvm =-
#endif /* #ifndef __CVM_CONSTANTS_H__ */
