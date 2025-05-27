#include <iostream>
#include <vector>
#include "atan2d_0_360.h"

/*
    This script tests the atan2d_0_360 function for various inputs.

    t - angle in degrees

    User h-function required: atan2d_0_360
*/
int main() 
{
    //...Inputs
    std::vector<std::pair<double, double>> inputs = 
    {
        {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
        {1, 1}, {-1, -1}, {1, -1}, {-1, 1},
        {std::sqrt(3), 1}, {-std::sqrt(3), -1},
        {std::sqrt(3), -1}, {-std::sqrt(3), 1}
    };

    int case_num = 1;
    for (const auto& [y, x] : inputs) 
    {
        double result = atan2d_0_360(y, x);
        std::cout << "atan2d_0_360(" << y << ", " << x << ") = " << result << "\n";
    }

    return 0;
}
