#include "los.h"
#include <iostream>

/*
    Uses the ECI position vectors of the satellite (r_sat) and the sun (r_sun)
    and evaluates whether the earth is in the line of sight between the two.

    User h-functions required: None
*/
void los_test() 
{
    // Scenario 1: Satellite in sunlight
    double r_sat1[3] = {7000.0, 0.0, 0.0};  // Satellite position vector (ECI, km)
    double r_sun1[3] = {1.496e8, 0.0, 0.0}; // Sun position vector (ECI, km)
    int light_switch1 = los(r_sat1, r_sun1);

    if (light_switch1 == 1) 
    {
        std::cout << "Satellite is in sunlight\n";
    } 
    else 
    {
        std::cout << "Satellite is not in sunlight\n";
    }

    // Scenario 2: Satellite in Earth's shadow
    double r_sat2[3] = {7000.0, 0.0, 0.0};   // Satellite position vector (ECI, km)
    double r_sun2[3] = {-1.496e8, 0.0, 0.0}; // Sun position vector (ECI, km)
    int light_switch2 = los(r_sat2, r_sun2);

    if (light_switch2 == 0) 
    {
        std::cout << "Satellite is in Earth's shadow\n";
    } 
    else 
    {
        std::cout << "Satellite is not in Earth's shadow\n";
    }

    // Scenario 3: Satellite on Earth's surface
    const double RE = 6378.0;               // Earth's radius (km)
    double r_sat3[3] = {RE, 0.0, 0.0};      // Satellite on Earth's surface
    double r_sun3[3] = {1.496e8, 0.0, 0.0}; // Sun position vector (ECI, km)
    int light_switch3 = los(r_sat3, r_sun3);

    if (light_switch3 == 1) 
    {
        std::cout << "Satellite on Earth's surface is in sunlight\n";
    } 
    else 
    {
        std::cout << "Satellite on Earth's surface is not in sunlight\n";
    }
}

int main() 
{
    los_test();
    return 0;
}