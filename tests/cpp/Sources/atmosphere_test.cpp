#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "atmosphere.h"

/*
    Tests the atmosphere function by plotting the density
    as a function of altitude from sea level through 1000 km.
*/
int main() {
    //...Geometric altitudes (km):
    const int num_points = 500;
    const double z_min = 0.0;
    const double z_max = 1000.0;
    std::vector<double> z(num_points);
    std::vector<double> density(num_points);

    for (int i = 0; i < num_points; ++i) {
        z[i] = z_min + i * (z_max - z_min) / (num_points - 1);
    }

    //...Calculate densities:
    for (int i = 0; i < num_points; ++i) {
        density[i] = atmosphere(z[i]);
    }

    //...Output the results:
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "Altitude (km)\tDensity (kg/m^3)\n";
    std::cout << "--------------------------------\n";
    for (int i = 0; i < num_points; ++i) {
        std::cout << z[i] << "\t" << density[i] << "\n";
    }

    return 0;
}
