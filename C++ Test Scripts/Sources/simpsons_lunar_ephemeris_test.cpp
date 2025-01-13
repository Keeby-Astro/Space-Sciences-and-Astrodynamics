#include "simpsons_lunar_ephemeris.h"
#include <iostream>
#include <iomanip>

/*
    Outputs the state vector of the moon at a given time
    relative to the earth's geocentric equatorial frame using a curve fit
    to JPL's DE200 (1982) ephemeris model.
    
    jd      - julian date (days)
    pos     - position vector (km)
    vel     - velocity vector (km/s)
    a       - matrix of amplitudes (km)
    b       - matrix of frequencies (rad/century)
    c       - matrix of phase angles (rad)
    t       - time in centuries since J2000
    tfac    - no. of seconds in a Julian century (36525 days)

    User h-functions required: simpsons_lunar_ephemeris
*/
int main() {
    // Define a sample Julian date (e.g., January 1, 2024, at 12:00 TT)
    double jd_sample = 2460271.0; // Julian date

    std::array<double, 3> pos;
    std::array<double, 3> vel;

    // Call the simpsons_lunar_ephemeris function
    simpsons_lunar_ephemeris(jd_sample, pos, vel);

    // Display the Julian date
    std::cout << "Julian date:\n";
    std::cout << "JD: " << std::fixed << std::setprecision(6) << jd_sample << "\n";

    // Display the position vector
    std::cout << "\nPosition vector (km):\n";
    std::cout << "X: " << pos[0] << " km\n";
    std::cout << "Y: " << pos[1] << " km\n";
    std::cout << "Z: " << pos[2] << " km\n";

    // Display the velocity vector
    std::cout << "\nVelocity vector (km/s):\n";
    std::cout << "VX: " << vel[0] << " km/s\n";
    std::cout << "VY: " << vel[1] << " km/s\n";
    std::cout << "VZ: " << vel[2] << " km/s\n";

    return 0;
}
