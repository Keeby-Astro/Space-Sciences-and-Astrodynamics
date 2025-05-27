#include <iostream>
#include <string>
#include "month_planet_names.h"

/*
    Outputs the name of the month and the planet corresponding, respectively,
    to the numbers "month_id" and "planet_id".

    month_id  - the month number (1 - 12)
    planet_id - the planet number (1 - 9)

    User h-functions required: month_planet_names
*/
int main() {
    std::cout << "Month and Planet Names Test" << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Test month and planet IDs
    for (int month_id = 1; month_id <= 12; ++month_id) {
        for (int planet_id = 1; planet_id <= 9; ++planet_id) {
            std::string month, planet;
            month_planet_names(month_id, planet_id, month, planet);

            // Output
            std::cout << "Month ID: " << month_id << ", Month: " << month
                      << " | Planet ID: " << planet_id << ", Planet: " << planet << std::endl;
        }
    }
    return 0;
}
