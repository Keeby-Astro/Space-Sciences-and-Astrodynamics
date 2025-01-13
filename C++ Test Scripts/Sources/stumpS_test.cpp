#include <iostream>
#include <iomanip>
#include <cassert>
#include "stumpS.h"

/*
    This function tests the stumpS function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User h-functions required: stumpS
*/
void stumpS_test() 
{
    // Define a range of test values for z
    double z_values[] = {-10, -1, -0.1, 0, 0.1, 1, 10};
    size_t n = sizeof(z_values) / sizeof(z_values[0]);

    // Display the results
    std::cout << "z\t\tS(z)\n";
    std::cout << "-------------------\n";

    // Loop through test values and compute S(z)
    for (size_t i = 0; i < n; ++i) 
    {
        double z = z_values[i];
        double s_values = stumpS(z);
        std::cout << std::fixed << std::setprecision(2) << z << "\t\t"
                  << std::setprecision(6) << s_values << "\n";
    }

    // Verify specific cases with known outcomes
    assert(std::abs(stumpS(0) - 1.0 / 6.0) < 1e-6 && "Test failed for z = 0");
    assert(std::abs(stumpS(1) - (std::sqrt(1) - std::sin(std::sqrt(1))) / std::pow(std::sqrt(1), 3)) < 1e-6 && "Test failed for z > 0");
    assert(std::abs(stumpS(-1) - (std::sinh(std::sqrt(1)) - std::sqrt(1)) / std::pow(std::sqrt(1), 3)) < 1e-6 && "Test failed for z < 0");

    std::cout << "\nAll tests passed successfully.\n";
}

int main() 
{
    stumpS_test();
    return 0;
}
