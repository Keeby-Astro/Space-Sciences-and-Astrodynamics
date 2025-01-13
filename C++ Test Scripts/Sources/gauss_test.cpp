#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "gauss.h"
#include "coe_from_sv.h"

/*
    This program uses Algorithms 5.5 and 5.6 (Gauss's method) to compute
    the state vector from the data provided.

    deg          - factor for converting between degrees and radians
    pi           - 3.1415926...
    mu           - gravitational parameter (km^3/s^2)
    Re           - earth's radius (km)
    f            - earth's flattening factor
    H            - elevation of observation site (km)
    phi          - latitude of site (deg)
    t            - vector of observation times t1, t2, t3 (s)
    ra           - vector of topocentric equatorial right ascensions
                    at t1, t2, t3 (deg)
    dec          - vector of topocentric equatorial right declinations
                    at t1, t2, t3 (deg)
    theta        - vector of local sidereal times for t1, t2, t3 (deg)
    R            - matrix of site position vectors at t1, t2, t3 (km)
    rho          - matrix of direction cosine vectors at t1, t2, t3
    fac1, fac2   - common factors
    r_old, v_old - the state vector without iterative improvement (km, km/s)
    r, v         - the state vector with iterative improvement (km, km/s)
    coe          - vector of orbital elements for r, v: 
                    [h, e, RA, incl, w, TA, a]
                    where h    = angular momentum (km^2/s)
                        e    = eccentricity
                        incl = inclination (rad)
                        w    = argument of perigee (rad)
                        TA   = true anomaly (rad)
                        a    = semimajor axis (km)
    coe_old     - vector of orbital elements for r_old, v_old       

    User h-functions required: gauss, coe_from_sv
*/
int main()
{
    double mu = 398600.4418;

    const double pi  = 3.141592653589793;
    const double deg = pi / 180.0;
    const double Re  = 6378.14;
    const double f   = 1.0 / 298.26;

    //...Data declaration:
    double H    = 1.0;       
    double phi  = 40.0 * deg; 
    std::vector<double> t = {0.0, 118.104, 237.577};
    std::vector<double> ra = {43.5365*deg, 54.4196*deg, 64.3178*deg};
    std::vector<double> dec = {-8.78334*deg, -12.0739*deg, -15.1054*deg};
    std::vector<double> theta = {44.5065*deg, 45.0*deg, 45.4992*deg};
    //...

    std::vector<std::vector<double>> R(3, std::vector<double>(3, 0.0));
    std::vector<std::vector<double>> rho(3, std::vector<double>(3, 0.0));

    //...Equations 5.64, 5.76 and 5.79:
    double fac1 = Re / std::sqrt(1.0 - (2.0*f - f*f)*std::pow(std::sin(phi), 2));
    double fac2 = (Re*(1.0 - f)*(1.0 - f) / 
                   std::sqrt(1.0 - (2.0*f - f*f)*std::pow(std::sin(phi), 2)) + H)
                  * std::sin(phi);
    for (int i = 0; i < 3; ++i)
    {
        R[i][0] = (fac1 + H) * std::cos(phi) * std::cos(theta[i]);
        R[i][1] = (fac1 + H) * std::cos(phi) * std::sin(theta[i]);
        R[i][2] = fac2;

        rho[i][0] = std::cos(dec[i]) * std::cos(ra[i]);
        rho[i][1] = std::cos(dec[i]) * std::sin(ra[i]);
        rho[i][2] = std::sin(dec[i]);
    }

    //...Algorithms 5.5 and 5.6:
    std::vector<double> r(3), v(3), r_old(3), v_old(3);
    gauss(rho[0], rho[1], rho[2],
          R[0],   R[1],   R[2],
          t[0],   t[1],   t[2],
          r, v, r_old, v_old, mu);

    //...Algorithm 4.2 for the initial estimate of the state vector
    //   and for the iteratively improved one:
    std::vector<double> coe_old = coe_from_sv(r_old, v_old, mu);
    std::vector<double> coe     = coe_from_sv(r,     v,     mu);

    //...Echo the input data and output the solution to
    //   the console:
    std::cout << "-----------------------------------------------------\n";
    std::cout << " Radius of earth (km)               = " << Re << "\n";
    std::cout << " Flattening factor                  = " << f << "\n";
    std::cout << " Gravitational parameter (km^3/s^2) = " << mu << "\n\n";
    
    std::cout << " Input data:\n";
    std::cout << "\n Latitude (deg)                = " << (phi / deg);
    std::cout << "\n Altitude above sea level (km) = " << H << "\n\n";
    
    std::cout << " Observations:\n";
    std::cout << "               Right                           Local\n";
    std::cout << "   Time (s)   Ascension (deg)   Declination (deg)   Sidereal time (deg)\n";
    for (int i = 0; i < 3; ++i) 
    {
        std::cout << " " << std::setw(9) << t[i]
                  << " " << std::setw(11) << (ra[i] / deg)
                  << " " << std::setw(19) << (dec[i] / deg)
                  << " " << std::setw(20) << (theta[i] / deg) << "\n";
    }
    
    std::cout << "\n Solution:\n";
    std::cout << "\n Without iterative improvement...\n\n";
    std::cout << " r (km)                          = ["
              << r_old[0] << ", " << r_old[1] << ", " << r_old[2] << "]\n";
    std::cout << " v (km/s)                        = ["
              << v_old[0] << ", " << v_old[1] << ", " << v_old[2] << "]\n\n";
    
    double h_old   = coe_old[0];
    double e_old   = coe_old[1];
    double RA_old  = coe_old[2];
    double incl_old= coe_old[3];
    double w_old   = coe_old[4];
    double TA_old  = coe_old[5];
    double a_old   = coe_old[6];
    
    std::cout << "   Angular momentum (km^2/s)     = " << h_old << "\n";
    std::cout << "   Eccentricity                  = " << e_old << "\n";
    std::cout << "   RA of ascending node (deg)    = " << (RA_old / deg) << "\n";
    std::cout << "   Inclination (deg)             = " << (incl_old / deg) << "\n";
    std::cout << "   Argument of perigee (deg)     = " << (w_old / deg) << "\n";
    std::cout << "   True anomaly (deg)            = " << (TA_old / deg) << "\n";
    std::cout << "   Semimajor axis (km)           = " << a_old << "\n";
    std::cout << "   Periapse radius (km)          = "
              << (h_old*h_old/mu/(1.0 + e_old)) << "\n";
    // If the orbit is an ellipse, output the period:
    if (e_old < 1.0)
    {
        double T = 2.0*pi / std::sqrt(mu) * std::pow(a_old, 1.5);
        std::cout << "   Period:\n";
        std::cout << "     Seconds                     = " << T << "\n";
        std::cout << "     Minutes                     = " << (T/60.0) << "\n";
        std::cout << "     Hours                       = " << (T/3600.0) << "\n";
        std::cout << "     Days                        = " << (T/86400.0) << "\n";
    }
    
    std::cout << "\n With iterative improvement...\n\n";
    std::cout << " r (km)                          = ["
              << r[0] << ", " << r[1] << ", " << r[2] << "]\n";
    std::cout << " v (km/s)                        = ["
              << v[0] << ", " << v[1] << ", " << v[2] << "]\n\n";

    double h_   = coe[0];
    double e_   = coe[1];
    double RA_  = coe[2];
    double incl_= coe[3];
    double w_   = coe[4];
    double TA_  = coe[5];
    double a_   = coe[6];

    std::cout << "   Angular momentum (km^2/s)     = " << h_ << "\n";
    std::cout << "   Eccentricity                  = " << e_ << "\n";
    std::cout << "   RA of ascending node (deg)    = " << (RA_ / deg) << "\n";
    std::cout << "   Inclination (deg)             = " << (incl_ / deg) << "\n";
    std::cout << "   Argument of perigee (deg)     = " << (w_ / deg) << "\n";
    std::cout << "   True anomaly (deg)            = " << (TA_ / deg) << "\n";
    std::cout << "   Semimajor axis (km)           = " << a_ << "\n";
    std::cout << "   Periapse radius (km)          = "
              << (h_*h_/mu/(1.0 + e_)) << "\n";
    // If the orbit is an ellipse, output the period:
    if (e_ < 1.0)
    {
        double T = 2.0*pi / std::sqrt(mu) * std::pow(a_, 1.5);
        std::cout << "   Period:\n";
        std::cout << "     Seconds                     = " << T << "\n";
        std::cout << "     Minutes                     = " << (T/60.0) << "\n";
        std::cout << "     Hours                       = " << (T/3600.0) << "\n";
        std::cout << "     Days                        = " << (T/86400.0) << "\n";
    }

    std::cout << "-----------------------------------------------------\n";

    return 0;
}
