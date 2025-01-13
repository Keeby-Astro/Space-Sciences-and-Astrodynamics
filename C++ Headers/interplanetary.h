// ALGORITHM 8.2: CALCULATION OF THE SPACECRAFT TRAJECTORY
// FROM PLANET 1 TO PLANET 2

#ifndef INTERPLANETARY_H
#define INTERPLANETARY_H

#include <vector>
#include <tuple>
#include <cmath>
#include "planet_elements_and_sv.h"
#include "lambert.h"

/*
    This function determines the spacecraft trajectory from the sphere
    of influence of planet 1 to that of planet 2 using Algorithm 8.2

    mu          - gravitational parameter of the sun (km^3/s^2)
    dum         - a dummy vector not required in this procedure

    planet_id   - planet identifier:
                  1 = Mercury
                  2 = Venus
                  3 = Earth
                  4 = Mars
                  5 = Jupiter
                  6 = Saturn
                  7 = Uranus
                  8 = Neptune
                  9 = Pluto

    year        - range: 1901 - 2099
    month       - range: 1 - 12
    day         - range: 1 - 31
    hour        - range: 0 - 23
    minute      - range: 0 - 60
    second      - range: 0 - 60

    jd1, jd2    - Julian day numbers at departure and arrival
    tof         - time of flight from planet 1 to planet 2 (s)
    Rp1, Vp1    - state vector of planet 1 at departure (km, km/s)
    Rp2, Vp2    - state vector of planet 2 at arrival (km, km/s)
    R1, V1      - heliocentric state vector of spacecraft at
                  departure (km, km/s)
    R2, V2      - heliocentric state vector of spacecraft at
                  arrival (km, km/s)

    depart      - [planet_id, year, month, day, hour, minute, second]
                  at departure
    arrive      - [planet_id, year, month, day, hour, minute, second]
                  at arrival

    planet1     - [Rp1, Vp1, jd1]
    planet2     - [Rp2, Vp2, jd2]
    trajectory  - [V1, V2]

    User h-functions required: planet_elements_and_sv, lambert
*/
inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
interplanetary(const std::vector<double>& depart,
               const std::vector<double>& arrive, double mu)
{
    int planet_id = static_cast<int>(depart[0]);
    int year      = static_cast<int>(depart[1]);
    int month     = static_cast<int>(depart[2]);
    int day       = static_cast<int>(depart[3]);
    int hour      = static_cast<int>(depart[4]);
    int minute    = static_cast<int>(depart[5]);
    int second    = static_cast<int>(depart[6]);

    //...Use Algorithm 8.1 to obtain planet 1's state vector (don't
    //...need its orbital elements ["dum"]):
    auto [dum1, Rp1, Vp1, jd1] = planet_elements_and_sv(
        planet_id, year, month, day, hour, minute, second, mu);

    planet_id = static_cast<int>(arrive[0]);
    year      = static_cast<int>(arrive[1]);
    month     = static_cast<int>(arrive[2]);
    day       = static_cast<int>(arrive[3]);
    hour      = static_cast<int>(arrive[4]);
    minute    = static_cast<int>(arrive[5]);
    second    = static_cast<int>(arrive[6]);

    //...Likewise use Algorithm 8.1 to obtain planet 2's state vector:
    auto [dum2, Rp2, Vp2, jd2] = planet_elements_and_sv(
        planet_id, year, month, day, hour, minute, second, mu);

    double tof = (jd2 - jd1) * 24.0 * 3600.0;

    //...Patched conic assumption:
    std::vector<double> R1 = Rp1;
    std::vector<double> R2 = Rp2;

    //...Use Algorithm 5.2 to find the spacecraft's velocity at
    // departure and arrival, assuming a prograde trajectory:
    auto [V1, V2] = lambert(R1, R2, tof, "pro", mu);

    std::vector<double> planet1(7), planet2(7), trajectory(6);

    for (int i = 0; i < 3; ++i) 
    {
        planet1[i] = Rp1[i];
        planet1[i+3] = Vp1[i];
        planet2[i] = Rp2[i];
        planet2[i+3] = Vp2[i];
        trajectory[i] = V1[i];
        trajectory[i+3] = V2[i];
    }
    planet1[6] = jd1;
    planet2[6] = jd2;

    return {planet1, planet2, trajectory};
}

#endif // INTERPLANETARY_H
