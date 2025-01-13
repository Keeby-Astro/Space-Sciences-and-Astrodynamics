#ifndef J0_H
#define J0_H

#include <iostream>
#include <cmath>

/*
    This function computes the Julian day number at 0 UT for any year
    between 1900 and 2100 using Equation 5.48.

    j0      - Julian day at 0 hr UT (Universal Time)
    year    - range: 1901 - 2099
    month   - range: 1 - 12
    day     - range: 1 - 31

    User h-functions required: none
*/
inline double J0(int year, int month, int day) 
{
    double j0 = 367.0 * year - std::floor(7 * (year + std::floor((month + 9) / 12.0)) / 4.0)
           + std::floor(275 * month / 9.0) + day + 1721013.5;
    return j0;
}

#endif // J0_H
