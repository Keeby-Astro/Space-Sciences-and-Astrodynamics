#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "planet_elements_and_sv.h"
#include "month_planet_names.h"

/*
    This program uses Algorithm 8.1 to compute the orbital elements
    and state vector of the earth at the date and time specified.
    To obtain the same results for Mars, set planet_id = 4.
    
    mu        - gravitational parameter of the sun (km^3/s^2)
    deg       - conversion factor between degrees and radians
    pi        - 3.1415926...
    
    coe       - vector of heliocentric orbital elements
                [h  e  RA  incl  w  TA  a  w_hat  L  M  E],
                where
                 h     = angular momentum                    (km^2/s)
                 e     = eccentricity
                 RA    = right ascension                     (deg)
                 incl  = inclination                         (deg)
                 w     = argument of perihelion              (deg)
                 TA    = true anomaly                        (deg)
                 a     = semimajor axis                      (km)
                 w_hat = longitude of perihelion ( = RA + w) (deg)
                 L     = mean longitude ( = w_hat + M)       (deg)
                 M     = mean anomaly                        (deg)
                 E     = eccentric anomaly                   (deg)
    
    r         - heliocentric position vector (km)
    v         - heliocentric velocity vector (km/s) 
    
    planet_id - planet identifier:
                 1 = Mercury
                 2 = Venus
                 3 = Earth
                 4 = Mars
                 5 = Jupiter
                 6 = Saturn
                 7 = Uranus
                 8 = Neptune
                 9 = Pluto
    
    year      - range: 1901 - 2099
    month     - range: 1 - 12
    day       - range: 1 - 31
    hour      - range: 0 - 23
    minute    - range: 0 - 60
    second    - range: 0 - 60 
    
    User h-functions required: planet_elements_and_sv, month_planet_names
*/
int main() 
{
    const double mu = 1.327124e11;
    
    //...Input data
    int planet_id = 3;
    int year = 2003;
    int month = 8;
    int day = 27;
    int hour = 12;
    int minute = 0;
    int second = 0;
    //...
    
    //...Algorithm 8.1:
    auto result = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second, mu);
    auto& [coe, r, v, jd] = result;
    
    //...Convert the planet_id and month numbers into names for output:
    std::string month_name, planet_name;
    month_planet_names(month, planet_id, month_name, planet_name);
    
    //...Echo the input data and output the solution to
    //   the command window:
    std::cout << "-----------------------------------------------------\n";
    std::cout << "\n Input data:\n";
    std::cout << "\n   Planet: " << planet_name;
    std::cout << "\n   Year  : " << year;
    std::cout << "\n   Month : " << month_name;
    std::cout << "\n   Day   : " << day;
    std::cout << "\n   Hour  : " << hour;
    std::cout << "\n   Minute: " << minute;
    std::cout << "\n   Second: " << second;
    std::cout << "\n\n Julian day: " << std::setprecision(12) << jd;
    
    std::cout << "\n\n Orbital elements:\n";
    std::cout << "\n  Angular momentum (km^2/s)                   = " << coe[0];
    std::cout << "\n  Eccentricity                                = " << coe[1];
    std::cout << "\n  Right ascension of the ascending node (deg) = " << coe[2];
    std::cout << "\n  Inclination to the ecliptic (deg)           = " << coe[3];
    std::cout << "\n  Argument of perihelion (deg)                = " << coe[4];
    std::cout << "\n  True anomaly (deg)                          = " << coe[5];
    std::cout << "\n  Semimajor axis (km)                         = " << coe[6];
    std::cout << "\n\n  Longitude of perihelion (deg)               = " << coe[7];
    std::cout << "\n  Mean longitude (deg)                        = " << coe[8];
    std::cout << "\n  Mean anomaly (deg)                          = " << coe[9];
    std::cout << "\n  Eccentric anomaly (deg)                     = " << coe[10];
    
    std::cout << "\n\n State vector:\n";
    std::cout << "\n  Position vector (km) = [" << r[0] << ", " << r[1] << ", " << r[2] << "]";
    std::cout << "\n  Magnitude            = " << std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    std::cout << "\n\n  Velocity (km/s)      = [" << v[0] << ", " << v[1] << ", " << v[2] << "]";
    std::cout << "\n  Magnitude            = " << std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    std::cout << "\n-----------------------------------------------------\n";
    
    return 0;
}