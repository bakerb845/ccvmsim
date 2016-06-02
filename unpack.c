#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "cvm.h"

/*!
 * @brief Converts a char from a binary file into a floating value
 *
 * @param[in] s      character string to unpack
 * @param[in] lswap  =True -> then swap bytes
 *                   =False -> then do not swap bytes (default)
 *
 * @result floating value corresponding to s from binary file
 *
 * @author Ben Baker, ISTI
 *
 */
float unpack_float(char *s, bool lswap)
{
    char s4[4];
    float f4;
    if (lswap)
    {
        s4[0] = s[3];
        s4[1] = s[2];
        s4[2] = s[1];
        s4[3] = s[0];
    }
    else
    {
        s4[0] = s[0];
        s4[1] = s[1];
        s4[2] = s[2];
        s4[3] = s[3];
    }
    memcpy(&f4, s4, 4);
    return f4;
}
//============================================================================//
/*!
 * @brief Converts a char from a binary file into a double value
 *
 * @param[in] s      character string to unpack
 * @param[in] lswap  =True -> then swap bytes
 *                   =False -> then do not swap bytes (default)
 *
 * @result double value corresponding to s from binary file
 *
 * @author Ben Baker, ISTI
 *
 */
double unpack_double(char *s, bool lswap)
{
    char s8[8];
    double f8; 
    if (lswap)
    {
        s8[0] = s[7];
        s8[1] = s[6];
        s8[2] = s[5];
        s8[3] = s[4];
        s8[4] = s[3];
        s8[5] = s[2];
        s8[6] = s[1];
        s8[7] = s[0];
    }
    else
    {
        s8[0] = s[0];
        s8[1] = s[1];
        s8[2] = s[2];
        s8[3] = s[3];
        s8[4] = s[4];
        s8[5] = s[5];
        s8[6] = s[6];
        s8[7] = s[7];
    }   
    memcpy(&f8, s8, 8); 
    return f8;
}

