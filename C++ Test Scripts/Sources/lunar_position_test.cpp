#include <iostream>
#include "lunar_position.h"

/*
    Outputs the geocentric equatorial position vector of the moon
    given the Julian day.
    
    User h-functions required: None
*/
int main() 
{
    std::array<double, 3> r_moon_1 = {0.0, 0.0, 0.0};
    double jd_1 = 2451545.0; // Julian Date for J2000 Epoch
    lunar_position(jd_1, r_moon_1);
    double jd_2 = 2459200.5; // Julian Date for September 22, 2020
    std::array<double, 3> r_moon_2 = {0.0, 0.0, 0.0};
    lunar_position(jd_2, r_moon_2);
    double jd_3 = 2440587.5; // Julian Date for January 1, 1970
    std::array<double, 3> r_moon_3 = {0.0, 0.0, 0.0};
    lunar_position(jd_3, r_moon_3);

    // Display results
    std::cout << "Julian Date = " << jd_1 << "\n";
    std::cout << "Geocentric Position Vector (r_moon):\n";
    std::cout << "  x = " << r_moon_1[0] << " km\n";
    std::cout << "  y = " << r_moon_1[1] << " km\n";
    std::cout << "  z = " << r_moon_1[2] << " km\n";

    std::cout << "Julian Date = " << jd_2 << "\n";
    std::cout << "Geocentric Position Vector (r_moon):\n";
    std::cout << "  x = " << r_moon_2[0] << " km\n";
    std::cout << "  y = " << r_moon_2[1] << " km\n";
    std::cout << "  z = " << r_moon_2[2] << " km\n";

    std::cout << "Julian Date = " << jd_3 << "\n";
    std::cout << "Geocentric Position Vector (r_moon):\n";
    std::cout << "  x = " << r_moon_3[0] << " km\n";
    std::cout << "  y = " << r_moon_3[1] << " km\n";
    std::cout << "  z = " << r_moon_3[2] << " km\n";
}