#ifndef MONTH_PLANET_NAMES_H
#define MONTH_PLANET_NAMES_H

#include <string>

/*
    This function returns the name of the month and the planet
    corresponding, respectively, to the numbers "month_id" and
    "planet_id".

    months    - an array containing the names of the 12 months
    planets   - an array containing the names of the 9 planets
    month_id  - the month number (1 - 12)
    planet_id - the planet number (1 - 9)

    User h-functions required: none
*/

void month_planet_names(int month_id, int planet_id, std::string &month, std::string &planet) {
    static const char *months[] = 
    {
        "January  ", 
        "February ", 
        "March    ", 
        "April    ", 
        "May      ", 
        "June     ",
        "July     ", 
        "August   ", 
        "September", 
        "October  ", 
        "November ", 
        "December "
    };

    static const char *planets[] = 
    {
        "Mercury", 
        "Venus  ", 
        "Earth  ", 
        "Mars   ", 
        "Jupiter",
        "Saturn ", 
        "Uranus ", 
        "Neptune", 
        "Pluto  "
    };

    if (month_id >= 1 && month_id <= 12) 
    {
        month = months[month_id - 1];
    } 
    else 
    {
        month = "Invalid Month";
    }

    if (planet_id >= 1 && planet_id <= 9) 
    {
        planet = planets[planet_id - 1];
    } 
    else 
    {
        planet = "Invalid Planet";
    }
}

#endif // MONTH_PLANET_NAMES_H