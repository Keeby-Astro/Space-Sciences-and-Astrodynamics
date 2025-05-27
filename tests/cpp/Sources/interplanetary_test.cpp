#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <cmath>
#include <string>
#include "interplanetary.h"
#include "coe_from_sv.h"
#include "month_planet_names.h"

static const double mu  = 1.327124e11; 
static const double PI  = 3.14159265358979;
static const double deg = PI / 180.0;

/*
    This program uses Algorithm 8.2 to solve interplanetary trajectory.

    mu           - gravitational parameter of the sun (km^3/s^2)
    deg          - conversion factor between degrees and radians
    pi           - 3.1415926...

    planet_id    - planet identifier:
                    1 = Mercury
                    2 = Venus
                    3 = Earth
                    4 = Mars
                    5 = Jupiter
                    6 = Saturn
                    7 = Uranus
                    8 = Neptune
                    9 = Pluto

    year         - range: 1901 - 2099
    month        - range: 1 - 12
    day          - range: 1 - 31
    hour         - range: 0 - 23
    minute       - range: 0 - 60
    second       - range: 0 - 60 

    depart       - [planet_id, year, month, day, hour, minute, second]
                    at departure
    arrive       - [planet_id, year, month, day, hour, minute, second]
                    at arrival

    planet1      - [Rp1, Vp1, jd1]
    planet2      - [Rp2, Vp2, jd2]
    trajectory   - [V1, V2]

    coe          - orbital elements [h e RA incl w TA]
                    where
                    h    = angular momentum (km^2/s)
                    e    = eccentricity
                    RA   = right ascension of the ascending
                            node (rad)
                    incl = inclination of the orbit (rad)
                    w    = argument of perigee (rad)
                    TA   = true anomaly (rad)
                    a    = semimajor axis (km)

    jd1, jd2     - Julian day numbers at departure and arrival
    tof          - time of flight from planet 1 to planet 2 (days)

    Rp1, Vp1     - state vector of planet 1 at departure (km, km/s)
    Rp2, Vp2     - state vector of planet 2 at arrival (km, km/s)
    R1, V1       - heliocentric state vector of spacecraft at
                    departure (km, km/s)
    R2, V2       - heliocentric state vector of spacecraft at
                    arrival (km, km/s)

    vinf1, vinf2 - hyperbolic excess velocities at departure
                    and arrival (km/s)

    User h-functions required: interplanetary, coe_from_sv,
                               month_planet_names
*/
int main()
{
    //...Data declaration:

    //...Departure
    std::vector<double> depart = 
    {
        3,     // Earth (planet_id)
        1996,  // year
        11,    // month
        7,     // day
        0,     // hour
        0,     // minute
        0      // second
    };
    std::string month_name_depart, planet_name_depart;
    int planet_id_depart = static_cast<int>(depart[0]);
    int month_depart = static_cast<int>(depart[2]);
    month_planet_names(month_depart, planet_id_depart, month_name_depart, planet_name_depart);

    //...Arrival
    std::vector<double> arrive = 
    {
        4,     // Mars (planet_id)
        1997,  // year
        9,     // month
        12,    // day
        0,     // hour
        0,     // minute
        0      // second
    };
    std::string month_name_arrive, planet_name_arrive;
    int planet_id_arrive = static_cast<int>(arrive[0]);
    int month_arrive = static_cast<int>(arrive[2]);
    month_planet_names(month_arrive, planet_id_arrive, month_name_arrive, planet_name_arrive);
    //...

    //...Algorithm 8.2:
    auto [planet1, planet2, trajectory] = interplanetary(depart, arrive, mu);

    std::vector<double> R1(3), Vp1(3), R2(3), Vp2(3), V1(3), V2(3);
    for (int i = 0; i < 3; ++i) 
    {
        R1[i]  = planet1[i];
        Vp1[i] = planet1[i+3];
        R2[i]  = planet2[i];
        Vp2[i] = planet2[i+3];
        V1[i]  = trajectory[i];
        V2[i]  = trajectory[i+3];
    }
    double jd1 = planet1[6];
    double jd2 = planet2[6];

    double tof = jd2 - jd1;

    //...Use Algorithm 4.2 to find the orbital elements of the
    //   spacecraft trajectory based on [Rp1, V1]...
    auto coe  = coe_from_sv(R1, V1, mu);
    //   ... and [R2, V2]
    auto coe2 = coe_from_sv(R2, V2, mu);

    //...Equations 8.94 and 8.95:
    std::vector<double> vinf1(3), vinf2(3);
    for (int i = 0; i < 3; ++i) 
    {
        vinf1[i] = V1[i] - Vp1[i];
        vinf2[i] = V2[i] - Vp2[i];
    }

    //...Echo the input data and output the solution to
    //   the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "\nDeparture:\n";
    std::cout << "  Planet: " << planet_name_depart << "\n";
    std::cout << "  Year  : " << depart[1] << "\n";
    std::cout << "  Month : " << month_name_depart << "\n";
    std::cout << "  Day   : " << depart[3] << "\n";
    std::cout << "  Hour  : " << depart[4] << "\n";
    std::cout << "  Minute: " << depart[5] << "\n";
    std::cout << "  Second: " << depart[6] << "\n\n";
    std::cout << "  Julian day: " << jd1 << "\n\n";

    std::cout << "  Planet position vector (km)    = ["
              << R1[0] << "  " << R1[1] << "  " << R1[2] << "]\n";
    double magR1 = std::sqrt(R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2]);
    std::cout << "  Magnitude                      = " << magR1 << "\n\n";

    std::cout << "  Planet velocity (km/s)         = ["
              << Vp1[0] << "  " << Vp1[1] << "  " << Vp1[2] << "]\n";
    double magVp1 = std::sqrt(Vp1[0]*Vp1[0] + Vp1[1]*Vp1[1] + Vp1[2]*Vp1[2]);
    std::cout << "  Magnitude                      = " << magVp1 << "\n\n";

    std::cout << "  Spacecraft velocity (km/s)     = ["
              << V1[0] << "  " << V1[1] << "  " << V1[2] << "]\n";
    double magV1 = std::sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2]);
    std::cout << "  Magnitude                      = " << magV1 << "\n\n";

    std::cout << "  v-infinity at departure (km/s) = ["
              << vinf1[0] << "  " << vinf1[1] << "  " << vinf1[2] << "]\n";
    double magVinf1 = std::sqrt(vinf1[0]*vinf1[0] + vinf1[1]*vinf1[1] + vinf1[2]*vinf1[2]);
    std::cout << "  Magnitude                      = " << magVinf1 << "\n\n";

    std::cout << "  Time of flight = " << tof << " days\n\n";

    std::cout << "Arrival:\n";
    std::cout << "  Planet: " << planet_name_arrive << "\n";
    std::cout << "  Year  : " << arrive[1] << "\n";
    std::cout << "  Month : " << month_name_arrive << "\n";
    std::cout << "  Day   : " << arrive[3] << "\n";
    std::cout << "  Hour  : " << arrive[4] << "\n";
    std::cout << "  Minute: " << arrive[5] << "\n";
    std::cout << "  Second: " << arrive[6] << "\n\n";
    std::cout << "  Julian day: " << jd2 << "\n\n";

    std::cout << "  Planet position vector (km)    = ["
              << R2[0] << "  " << R2[1] << "  " << R2[2] << "]\n";
    double magR2 = std::sqrt(R2[0]*R2[0] + R2[1]*R2[1] + R2[2]*R2[2]);
    std::cout << "  Magnitude                      = " << magR2 << "\n\n";

    std::cout << "  Planet velocity (km/s)         = ["
              << Vp2[0] << "  " << Vp2[1] << "  " << Vp2[2] << "]\n";
    double magVp2 = std::sqrt(Vp2[0]*Vp2[0] + Vp2[1]*Vp2[1] + Vp2[2]*Vp2[2]);
    std::cout << "  Magnitude                      = " << magVp2 << "\n\n";

    std::cout << "  Spacecraft Velocity (km/s)     = ["
              << V2[0] << "  " << V2[1] << "  " << V2[2] << "]\n";
    double magV2 = std::sqrt(V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2]);
    std::cout << "  Magnitude                      = " << magV2 << "\n\n";

    std::cout << "  v-infinity at arrival (km/s)   = ["
              << vinf2[0] << "  " << vinf2[1] << "  " << vinf2[2] << "]\n";
    double magVinf2 = std::sqrt(vinf2[0]*vinf2[0] + vinf2[1]*vinf2[1] + vinf2[2]*vinf2[2]);
    std::cout << "  Magnitude                      = " << magVinf2 << "\n\n";

    std::cout << "Orbital elements of flight trajectory:\n";
    std::cout << "  Angular momentum (km^2/s)                   = " << coe[0] << "\n";
    std::cout << "  Eccentricity                                = " << coe[1] << "\n";
    std::cout << "  Right ascension of ascending node (deg)     = " << (coe[2] / deg) << "\n";
    std::cout << "  Inclination (deg)                           = " << (coe[3] / deg) << "\n";
    std::cout << "  Argument of perihelion (deg)                = " << (coe[4] / deg) << "\n";
    std::cout << "  True anomaly at departure (deg)             = " << (coe[5] / deg) << "\n";
    std::cout << "  True anomaly at arrival (deg)               = " << (coe2[5] / deg) << "\n";
    std::cout << "  Semimajor axis (km)                         = " << coe[6] << "\n";
    // If the orbit is an ellipse, output the period:
    if (coe[1] < 1.0) {
        double period_seconds = 2.0 * PI * std::sqrt(std::pow(coe[6], 3) / mu);
        double period_days    = period_seconds / (24.0 * 3600.0);
        std::cout << "  Period (days)                               = " << period_days << "\n";
    }
    std::cout << "-----------------------------------------------------\n";

    return 0;
}