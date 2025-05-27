#ifndef LST_H
#define LST_H

#include <iostream>
#include <cmath>
#include "J0.h"

/*
    This function calculates the local sidereal time.
    lst - local sidereal time (degrees)

    y   - year
    m   - month
    d   - day
    ut  - Universal Time (hours)
    EL  - east longitude (degrees)
    j0  - Julian day number at 0 hr UT
    j   - number of centuries since J2000
    g0  - Greenwich sidereal time (degrees) at 0 hr UT
    gst - Greenwich sidereal time (degrees) at the specified UT

    User h-function required: J0
    User subfunction required: zeroTo360
*/
inline double zeroTo360(double x)
{
    /*
        This subfunction reduces an angle to the range 0 - 360 degrees.
    
        x - The angle (degrees) to be reduced
        y - The reduced value 
    */
    if (x >= 360.0) 
    {
        x -= std::floor(x / 360.0) * 360.0;
    } 
    else if (x < 0.0) 
    {
        x -= (std::floor(x / 360.0) - 1.0) * 360.0;
    }
    return x;
}

inline double LST(int year, int month, int day, double ut, double EL) 
{
    //...Equation 5.48:
    double j0 = J0(year, month, day);

    //...Equation 5.49:
    double j = (j0 - 2451545.0) / 36525.0;

    //...Equation 5.50:
    double g0 = 100.4606184 + 36000.77004 * j + 0.000387933 * j * j - 2.583e-8 * j * j * j;

    //...Reduce g0 so it lies in the range 0 - 360 degrees
    g0 = zeroTo360(g0);

    //...Equation 5.51:
    double gst = g0 + 360.98564724 * ut / 24.0;

    //...Equation 5.52:
    double lst = gst + EL;

    //...Reduce lst to the range 0 - 360 degrees:
    lst = lst - 360.0 * std::floor(lst / 360.0);

    return lst;
}

#endif // LST_H
