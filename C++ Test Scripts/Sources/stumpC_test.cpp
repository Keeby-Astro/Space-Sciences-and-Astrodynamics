#include <iostream>
#include <iomanip>
#include <cassert>
#include "stumpC.h"

/*
    This function tests the stumpC function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User h-functions required: stumpC
*/
void stumpC_test() 
{
    // Define a range of test values for z
    double z_values[] = {-10, -1, -0.1, 0, 0.1, 1, 10};
    size_t n = sizeof(z_values) / sizeof(z_values[0]);

    // Display the results
    std::cout << "z\t\tC(z)\n";
    std::cout << "-------------------\n";

    // Loop through test values and compute C(z)
    for (size_t i = 0; i < n; ++i) 
    {
        double z = z_values[i];
        double c_values = stumpC(z);
        std::cout << std::fixed << std::setprecision(2) << z << "\t\t"
                  << std::setprecision(6) << c_values << "\n";
    }

    // Verify specific cases with known outcomes
    assert(std::abs(stumpC(0) - 0.5) < 1e-6 && "Test failed for z = 0");
    assert(std::abs(stumpC(1) - (1.0 - std::cos(std::sqrt(1))) / 1.0) < 1e-6 && "Test failed for z > 0");
    assert(std::abs(stumpC(-1) - (std::cosh(std::sqrt(1)) - 1.0) / 1.0) < 1e-6 && "Test failed for z < 0");

    std::cout << "\nAll tests passed successfully.\n";
}

int main() 
{
    stumpC_test();
    return 0;
}
