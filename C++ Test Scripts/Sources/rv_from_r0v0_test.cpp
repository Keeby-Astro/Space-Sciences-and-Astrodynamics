#include <iostream>
#include <vector>
#include "rv_from_r0v0.h"

/*
    This program computes the state vector (R,V) from the initial
    state vector (R0,V0) and the elapsed time using the data.

    mu - gravitational parameter (km^3/s^2)
    R0 - the initial position vector (km)
    V0 - the initial velocity vector (km/s)
    R  - the final position vector (km)
    V  - the final velocity vector (km/s)
    t  - elapsed time (s)

    User h-functions required: rv_from_r0v0
*/
int main() 
{
    const double mu = 398600.4418;

    //...Data declaration:
    std::vector<double> R0 = {7000.0, -12124.0, 0.0};
    std::vector<double> V0 = {2.6679, 4.6210, 0.0};
    double t = 3600.0;
    //...

    std::vector<double> R(3), V(3);

    //...Algorithm 3.4:
    rv_from_r0v0(R0, V0, t, R, V, mu);

    //...Echo the input data and output the results to the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Initial position vector (km):\n";
    std::cout << "r0 = (" << R0[0] << ", " << R0[1] << ", " << R0[2] << ")\n";

    std::cout << "\nInitial velocity vector (km/s):\n";
    std::cout << "v0 = (" << V0[0] << ", " << V0[1] << ", " << V0[2] << ")\n";

    std::cout << "\nElapsed time = " << t << " s\n";

    std::cout << "\nFinal position vector (km):\n";
    std::cout << "r = (" << R[0] << ", " << R[1] << ", " << R[2] << ")\n";

    std::cout << "\nFinal velocity vector (km/s):\n";
    std::cout << "v = (" << V[0] << ", " << V[1] << ", " << V[2] << ")\n";
    std::cout << "-----------------------------------------------------\n";

    return 0;
}
