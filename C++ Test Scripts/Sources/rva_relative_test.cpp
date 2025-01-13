#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include "sv_from_coe.h"
#include "rva_relative.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

/*
    This program uses the data to calculate the position,
    velocity and acceleration of an orbiting chaser B relative to an
    orbiting target A.

    mu               - gravitational parameter (km^3/s^2)
    deg              - conversion factor from degrees to radians

                    Spacecraft A & B:
    h_A, h_B         -   angular momentum (km^2/s)
    e_A, E_B         -   eccentricity
    i_A, i_B         -   inclination (radians)
    RAAN_A, RAAN_B   -   right ascension of the ascending node (radians)
    omega_A, omega_B -   argument of perigee (radians)
    theta_A, theta_A -   true anomaly (radians)

    rA, vA           - inertial position (km) and velocity (km/s) of A
    rB, vB           - inertial position (km) and velocity (km/s) of B
    r                - position (km) of B relative to A in A's
                        co-moving frame
    v                - velocity (km/s) of B relative to A in A's
                        co-moving frame
    a                - acceleration (km/s^2) of B relative to A in A's
                        co-moving frame

    User h-function required:   sv_from_coe, rva_relative
*/
int main() 
{
    const double mu = 398600.4418;
    const double deg = M_PI / 180.0;

    //...Input data:

    //   Spacecraft A
    double h_A = 52059.0;
    double e_A = 0.025724;
    double i_A = 60.0 * deg;
    double RAAN_A = 40.0 * deg;
    double omega_A = 30.0 * deg;
    double theta_A = 40.0 * deg;

    //   Spacecraft B
    double h_B = 52362.0;
    double e_B = 0.0072696;
    double i_B = 50.0 * deg;
    double RAAN_B = 40.0 * deg;
    double omega_B = 120.0 * deg;
    double theta_B = 40.0 * deg;

    //...End input data

    //...Compute the initial state vectors of A and B using Algorithm 4.5:
    auto [rA, vA] = sv_from_coe({h_A, e_A, RAAN_A, i_A, omega_A, theta_A}, mu);
    auto [rB, vB] = sv_from_coe({h_B, e_B, RAAN_B, i_B, omega_B, theta_B}, mu);

    //...Compute relative position, velocity and acceleration using
    //   Algorithm 7.1:
    auto [r_rel, v_rel, a_rel] = rva_relative(rA, vA, rB, vB, mu);

    //...Output
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n-----------------------------------------------\n";
    std::cout << "\nOrbital elements of spacecraft A:\n";
    std::cout << "  angular momentum = " << h_A << " (km^2/s)\n";
    std::cout << "  eccentricity = " << e_A << "\n";
    std::cout << "  inclination = " << i_A / deg << " (deg)\n";
    std::cout << "  RAAN = " << RAAN_A / deg << " (deg)\n";
    std::cout << "  argument of perigee = " << omega_A / deg << " (deg)\n";
    std::cout << "  true anomaly = " << theta_A / deg << " (deg)\n";

    std::cout << "\nState vector of spacecraft A:\n";
    std::cout << "  r = [" << rA[0] << ", " << rA[1] << ", " << rA[2] << "] km\n";
    std::cout << "      (magnitude = " << std::sqrt(rA[0] * rA[0] + rA[1] * rA[1] + rA[2] * rA[2]) << " km)\n";
    std::cout << "  v = [" << vA[0] << ", " << vA[1] << ", " << vA[2] << "] km/s\n";
    std::cout << "      (magnitude = " << std::sqrt(vA[0] * vA[0] + vA[1] * vA[1] + vA[2] * vA[2]) << " km/s)\n";

    std::cout << "\nOrbital elements of spacecraft B:\n";
    std::cout << "  angular momentum = " << h_B << " (km^2/s)\n";
    std::cout << "  eccentricity = " << e_B << "\n";
    std::cout << "  inclination = " << i_B / deg << " (deg)\n";
    std::cout << "  RAAN = " << RAAN_B / deg << " (deg)\n";
    std::cout << "  argument of perigee = " << omega_B / deg << " (deg)\n";
    std::cout << "  true anomaly = " << theta_B / deg << " (deg)\n";

    std::cout << "\nState vector of spacecraft B:\n";
    std::cout << "  r = [" << rB[0] << ", " << rB[1] << ", " << rB[2] << "] km\n";
    std::cout << "      (magnitude = " << std::sqrt(rB[0] * rB[0] + rB[1] * rB[1] + rB[2] * rB[2]) << " km)\n";
    std::cout << "  v = [" << vB[0] << ", " << vB[1] << ", " << vB[2] << "] km/s\n";
    std::cout << "      (magnitude = " << std::sqrt(vB[0] * vB[0] + vB[1] * vB[1] + vB[2] * vB[2]) << " km/s)\n";

    std::cout << "\nIn the co-moving fram attached to A:\n";
    std::cout << "  Position vector of B relative to A: [" << r_rel[0] << ", " << r_rel[1] << ", " << r_rel[2] << "] km\n";
    std::cout << "     (magnitude = " << std::sqrt(r_rel[0] * r_rel[0] + r_rel[1] * r_rel[1] + r_rel[2] * r_rel[2]) << " km)\n";
    std::cout << "  Velocity vector of B relative to A: [" << v_rel[0] << ", " << v_rel[1] << ", " << v_rel[2] << "] km/s\n";
    std::cout << "     (magnitude = " << std::sqrt(v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]) << " km/s)\n";
    std::cout << "  Acceleration vector of B relative to A: [" << a_rel[0] << ", " << a_rel[1] << ", " << a_rel[2] << "] km/s^2\n";
    std::cout << "     (magnitude = " << std::sqrt(a_rel[0] * a_rel[0] + a_rel[1] * a_rel[1] + a_rel[2] * a_rel[2]) << " km/s^2)\n";

    std::cout << "\n-----------------------------------------------\n";

    return 0;
}
