#include "f_and_g.h"
#include <iostream>
#include <iomanip>

/*
    This function tests the f_and_g function by providing sample inputs
    and verifying the outputs. 

    User h-functions required: f_and_g, stumpC, stumpS
*/
int main() {

    // Define the global variable mu (gravitational parameter)
    double mu = 398600.4418; // km^3/s^2, standard for Earth

    //...Inputs
    double ro = 7000.0;       // km, radial position at time to
    double t = 500.0;         // s, time elapsed
    double x = 1.5;           // universal anomaly
    double a = 1.0 / 10000.0; // 1/km, reciprocal of the semimajor axis

    // Outputs
    double f = 0.0, g = 0.0;

    // Call the f_and_g function
    f_and_g(x, t, ro, a, f, g, mu);

    //...Results
    std::cout << "Results:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Computed f: " << f << "\n";
    std::cout << "Computed g: " << g << "\n";

    return 0;
}
