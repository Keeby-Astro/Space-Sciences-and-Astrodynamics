#ifndef ATAN2D_0_360_H
#define ATAN2D_0_360_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#define M_PI 3.14159265358979323846

/*
    This function calculates the arc tangent of y/x in degrees
    and places the result in the range [0, 360].

    t - angle in degrees
*/
inline double atan2d_0_360(double y, double x) 
{
    if (x == 0) 
    {
        if (y == 0) 
        {
            return 0.0;
        } 
        else if (y > 0) 
        {
            return 90.0;
        } 
        else 
        {
            return 270.0;
        }
    }
    else if (x > 0) 
     {
        if (y >= 0) 
        {
            return std::atan2(y, x) * 180.0 / M_PI;
        }
        else 
        {
            return std::atan2(y, x) * 180.0 / M_PI + 360.0;
        }
    }
    else 
    {
        if (y == 0) 
        {
            return 180.0;
        } else 
        {
            return std::atan2(y, x) * 180.0 / M_PI + 180.0;
        }
    }
}

#endif // ATAN2D_0_360_H
