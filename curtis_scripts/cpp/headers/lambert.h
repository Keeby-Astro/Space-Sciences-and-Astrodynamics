#ifndef LAMBERT_H
#define LAMBERT_H

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <stdexcept>
#include "stumpC.h"
#include "stumpS.h"

/*
    This function solves Lambert's problem.

    mu      - gravitational parameter (km^3/s^2)
    R1, R2  - initial and final position vectors (km)
    r1, r2  - magnitudes of R1 and R2
    t       - the time of flight from R1 to R2 (a constant) (s)
    V1, V2  - initial and final velocity vectors (km/s)
    c12     - cross product of R1 into R2
    theta   - angle between R1 and R2
    string  - 'pro' if the orbit is prograde
              'retro' if the orbit is retrograde
    A       - a constant given by Equation 5.35
    z       - alpha*x^2, where alpha is the reciprocal of the
              semimajor axis and x is the universal anomaly
    y(z)    - a function of z given by Equation 5.38
    F(z,t)  - a function of the variable z and constant t,
            - given by Equation 5.40
    dFdz(z) - the derivative of F(z,t), given by Equation 5.43
    ratio   - F/dFdz
    tol     - tolerance on precision of convergence
    nmax    - maximum number of iterations of Newton's procedure
    f, g    - Lagrange coefficients
    gdot    - time derivative of g
              C(z), S(z) - Stumpff functions
    dum     - a dummy variable

    User h-functions required: stumpC and stumpS
*/
std::pair<std::vector<double>, std::vector<double>> lambert(const std::vector<double>& R1, 
const std::vector<double>& R2, double t, const std::string& orbit_type, double mu) 
{
    if (R1.size() != 3 || R2.size() != 3) 
    {
        throw std::invalid_argument("Input vectors R1 and R2 must each have exactly three elements.");
    }

    //...Magnitudes of R1 and R2:
    double r1 = std::sqrt(R1[0]*R1[0] + R1[1]*R1[1] + R1[2]*R1[2]);
    double r2 = std::sqrt(R2[0]*R2[0] + R2[1]*R2[1] + R2[2]*R2[2]);

    std::vector<double> c12 = 
    {
        R1[1]*R2[2] - R1[2]*R2[1],
        R1[2]*R2[0] - R1[0]*R2[2],
        R1[0]*R2[1] - R1[1]*R2[0]
    };
    double dot_product = R1[0]*R2[0] + R1[1]*R2[1] + R1[2]*R2[2];
    double theta = std::acos(dot_product / (r1 * r2));

    //...Determine whether the orbit is prograde or retrograde:
    std::string orbit = orbit_type;
    if (orbit != "pro" && orbit != "retro") 
    {
        orbit = "pro";
        std::cerr << "** Prograde trajectory assumed. **\n";
    }

    if (orbit == "pro" && c12[2] <= 0.0) 
    {
        theta = 2.0 * M_PI - theta;
    } 
    else if (orbit == "retro" && c12[2] >= 0.0) 
    {
        theta = 2.0 * M_PI - theta;
    }

    //...Equation 5.35:
    double A = std::sin(theta) * std::sqrt(r1 * r2 / (1.0 - std::cos(theta)));

    //...Stumpff functions:
    auto C = [](double z) { return stumpC(z); };
    auto S = [](double z) { return stumpS(z); };

    //...Equation 5.38:
    auto y = [&](double z) 
    {
        double Cz = C(z);
        if (Cz <= 0) return std::numeric_limits<double>::infinity();
        return r1 + r2 + A * ((z * S(z)) - 1.0) / std::sqrt(Cz);
    };

    //...Determine approximately where F(z,t) changes sign, and
    //...use that value of z as the starting value for Equation 5.45:
    auto F = [&](double z, double t_) 
    {
        double yz = y(z);
        double Cz = C(z);
        if (Cz <= 0 || yz <= 0) return std::numeric_limits<double>::infinity();
        return std::pow(yz / Cz, 1.5) * S(z) + A * std::sqrt(yz) - std::sqrt(mu) * t_;
    };

    //...Equation 5.43:
    auto dFdz = [&](double z) 
    {
        if (std::fabs(z) < 1.e-12) z = 1.e-12;
        double Cz = C(z);
        double Sz = S(z);
        double yz = y(z);
        return std::pow(yz / Cz, 1.5) * ((1.0 / (2.0 * z)) * (Cz - (3.0 * Sz / (2.0 * Cz))) + 3.0 * std::pow(Sz, 2) / (4.0 * Cz)) + (A / 8.0) * (3.0 * (Sz / Cz) * std::sqrt(yz) + A * std::sqrt(Cz / yz));
    };

    double z = -100.0;
    //...Equation 5.40:
    while (F(z, t) < 0.0) 
    {
        z += 0.1;
        if (z > 1.e6) 
        {
            std::cerr << "** F(z,t) did not become positive. Possibly no single-rev solution. **\n";
            break;
        }
    }

    //...Set an error tolerance and a limit on the number of iterations:
    const double tol = 1.0e-8;
    const int nmax = 5000;
    int n = 0;
    double ratio = 1.0;

    z = std::fabs(z);
    while (std::fabs(ratio) > tol && n < nmax) 
    {
        n++;
        ratio = F(z, t) / dFdz(z);
        z -= ratio;
    }

    //...Report if the maximum number of iterations is exceeded:
    if (n >= nmax) 
    {
        std::cerr << "** Number of iterations (" << nmax << ") exceeded **\n";
    }

    double yz = y(z);

    //...Equation 5.46a:
    double f = 1.0 - (yz / r1);

    //...Equation 5.46b:
    double g = A * std::sqrt(yz / mu);

    //...Equation 5.46d:
    double gdot = 1.0 - (yz / r2);

    std::vector<double> V1(3), V2(3);
    for (int i = 0; i < 3; ++i) 
    {
        //...Equation 5.28:
        V1[i] = (1.0 / g) * (R2[i] - f * R1[i]);

        //...Equation 5.29:
        V2[i] = (1.0 / g) * (gdot * R2[i] - R1[i]);
    }

    return {V1, V2};
}

#endif // LAMBERT_H
