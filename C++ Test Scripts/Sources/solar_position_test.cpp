#include <iostream>
#include "solar_position.h"

/*
    Outputs the geocentric equatorial position vector
    of the sun, given the julian date.
    
    User h-functions required: None
*/
int main() 
{
    // Julian date for 2023-12-21 00:00:00 UTC (Winter Solstice)
    double jd_test = 2460271.5;

    double lamda = 0.0;
    double eps = 0.0;
    std::array<double, 3> r_S = {0.0, 0.0, 0.0};

    // Call the solar_position function
    solar_position(jd_test, lamda, eps, r_S);

    // Display the results
    std::cout << "Julian Date = " << jd_test << "\n";
    std::cout << "Apparent Ecliptic Longitude (lamda): " << lamda << " degrees\n";
    std::cout << "Obliquity of the Ecliptic (eps): " << eps << " degrees\n";
    std::cout << "Geocentric Position Vector (r_S):\n";
    std::cout << "  x = " << r_S[0] << " km\n";
    std::cout << "  y = " << r_S[1] << " km\n";
    std::cout << "  z = " << r_S[2] << " km\n";

    return 0;
}