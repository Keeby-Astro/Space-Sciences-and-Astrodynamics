#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <cmath>
#include <algorithm>
/*
    ATMOSPHERE calculates density for altitudes from sea level
    through 1000 km using exponential interpolation.
*/
inline double atmosphere(double z) 
{
    //...Geometric altitudes (km):
    const double h[] = 
    {
        0, 25, 30, 40, 50, 60, 70,
        80, 90, 100, 110, 120, 130, 140,
        150, 180, 200, 250, 300, 350, 400,
        450, 500, 600, 700, 800, 900, 1000
    };

    //...Corresponding densities (kg/m^3) from USSA76:
    const double r[] = 
    {
        1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5,
        1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9,
        2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12,
        1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15
    };

    //...Scale heights (km):
    const double H[] = 
    {
        7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661,
        5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322,
        21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019,
        60.980, 65.654, 76.377, 100.587, 147.203, 208.020
    };

    //...Handle altitudes outside of the range:
    z = std::clamp(z, 0.0, 1000.0);

    //...Determine the interpolation interval:
    int i = 0;
    for (int j = 0; j < 27; ++j) 
    {
        if (z >= h[j] && z < h[j + 1]) 
        {
            i = j;
            break;
        }
    }
    if (z == 1000) 
    {
        i = 26;
    }

    //...Exponential interpolation:
    return r[i] * std::exp(-(z - h[i]) / H[i]);
}

#endif // ATMOSPHERE_H
