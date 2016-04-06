#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "cvm.h"
/*!
 * @brief Borcher's empirical relationship between compressional velocity
 *        and density
 *
 * @reference Brocher, 2005.  
 *            Empirical relations between elastic wavesppeds and 
 *            density in the Earth's crust; BSSA v 95, no 6 p 2081-2092
 *
 * @param[in] Vp_ms  compressional velocity (m/s)
 *
 * @result corresponding density (kg/m**3)
 *
 * @author Ben Baker, ISTI
 *
 */
#pragma omp declare simd
double density_Brocher(double Vp_ms)
{
    double rho, Vp, Vp2, Vp3, Vp4, Vp5;
    Vp = Vp_ms*1.e-3;  // Convert to km/s
    Vp2 = Vp*Vp;
    Vp3 = Vp*Vp2;
    Vp4 = Vp*Vp3;
    Vp5 = Vp*Vp4; 
    rho = 1.6612*Vp - 0.4721*Vp2 + 0.0671*Vp3 - 0.0043*Vp4 + 0.000106*Vp5;
    rho = rho*1.e3; // Convert to kg/m**3
    return rho;
}
//============================================================================//
/*!
 * @brief Darcy-McPhees Empircal relationship between compressional velocity
 *        and density 
 *
 * @reference Roecker, Thurber, and McPhee (2004).  Joint inversion of 
 *            arrival time data from Parkfield: New constraints on 
 *            structure and hypocenter locations near SAFOD drill site.
 *            GRL, 31, L12S04.
 *
 * @param[in] Vp_ms  compressional velocity (m/s)
 *
 * @result corresponding density (kg/m**3)
 *
 * @author Steve Roecker, RPI, and Ben Baker, ISTI
 *
 */
double density_DarcyMcphee(double Vp_ms)
{
    const double dcoef[9] = {-21893558144627.e-7,  30290041149786.e-7,
                             -18300316791235.e-7,  63062641966165.e-8,
                             -13556725156168.e-8,  18616741015779.e-9,
                             -15948394116079.e-10, 77924797412933.e-12,
                             -16626306058716.e-13};
    double rho, pgrcm3, vkmps;
    int ij;
    vkmps = Vp_ms*1.e-3; // Convert to km/s
    if (vkmps >= 5.93){
        pgrcm3 = 0.9893 + 0.2891*vkmps;
    }else if (vkmps >= 5.5 && vkmps < 5.93){
        pgrcm3 = 0.0;
        for (ij=1; ij<=8; ij++){
            pgrcm3 = (pgrcm3 + dcoef[9-ij])*vkmps;
        }
        pgrcm3 = pgrcm3 + dcoef[0];
    }else{
        pgrcm3 = 1.7407*pow(vkmps, 0.25);
    } 
    rho = pgrcm3*1000.;
    return rho;
}
