// ALGORITHM 10.3: DETERMINE WHETHER OR NOT A SATELLITE IS IN EARTH'S SHADOW

#ifndef LOS_H
#define LOS_H

#include <cmath>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

/*
    This function uses the ECI position vectors of the satellite (r_sat)
    and the sun (r_sun) to determine whether the earth is in the line of
    sight between the two.
    
    User h-functions required: None
*/
inline int los(const double r_sat[3], const double r_sun[3]) 
{
    const double RE = 6378.0; // Earth's radius (km)
    double rsat = std::sqrt(r_sat[0]*r_sat[0] + r_sat[1]*r_sat[1] + r_sat[2]*r_sat[2]);
    double rsun = std::sqrt(r_sun[0]*r_sun[0] + r_sun[1]*r_sun[1] + r_sun[2]*r_sun[2]);

    double dot_product = r_sat[0]*r_sun[0] + r_sat[1]*r_sun[1] + r_sat[2]*r_sun[2];

    //...Angle between sun and satellite position vectors:
    double theta = std::acos(dot_product / (rsat * rsun)) * 180.0 / M_PI;

    //...Angle between the satellite position vector and the radial to the point
    // of tangency with the earth of a line from the satellite:
    double theta_sat = std::acos(RE / rsat) * 180.0 / M_PI;

    //...Angle between the sun position vector and the radial to the point
    // of tangency with the earth of a line from the sun:    
    double theta_sun = std::acos(RE / rsun) * 180.0 / M_PI;

    //...Determine whether a line from the sun to the satellite
    // intersects the earth:
    if (theta_sat + theta_sun <= theta) 
    {
        return 0; //yes
    } 
    else 
    {
        return 1; //no
    }
}

#endif // LOS_H
