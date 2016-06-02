#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "cvm.h"
/*!
 * @brief Computes the quality factor ala Art Frankel
 *
 * @param[in] vs      shear velocity (m/s)
 *
 * @param[out] Qp     P velocity quality factor
 * @param[out] Qs     S velocity quality factor
 *
 */
#pragma omp declare simd
void qualityFactor_Frankel(double vs, double *Qp, double *Qs)
{
    *Qs = 0.15*vs;
    if (vs <= 1000.0){*Qs = 0.16428*vs - 14.2857;}
    *Qp = 2.0**Qs;
    return;
}
