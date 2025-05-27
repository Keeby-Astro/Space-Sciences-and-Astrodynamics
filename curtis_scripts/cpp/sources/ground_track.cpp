#include <iostream>
#include <vector>
#include <cmath>
#include "sv_from_coe.h"
#include "kepler_E.h"
#include "ra_and_dec_from_r.h"

/*
    This program plots the ground track of an earth satellite
    for which the orbital elements are specified

    mu          - gravitational parameter (km^3/s^2)
    deg         - factor that converts degrees to radians
    J2          - second zonal harmonic
    Re          - earth s radius (km)
    we          - earth s angular velocity (rad/s)
    rP          - perigee of orbit (km)
    rA          - apogee of orbit (km)
    TA, TAo     - true anomaly, initial true anomaly of satellite (rad)
    RA, RAo     - right ascension, initial right ascension of the node (rad)
    incl        - orbit inclination (rad)
    wp, wpo     - argument of perigee, initial argument of perigee (rad)
    n_periods   - number of periods for which ground track is to be plotted
    a           - semimajor axis of orbit (km)
    T           - period of orbit (s)
    e           - eccentricity of orbit
    h           - angular momentum of orbit (km^2/s)
    E, Eo       - eccentric anomaly, initial eccentric anomaly (rad)
    M, Mo       - mean anomaly, initial mean anomaly (rad)
    to, tf      - initial and final times for the ground track (s)
    fac         - common factor in Equations 4.53 and 4.53
    RAdot       - rate of regression of the node (rad/s)
    wpdot       - rate of advance of perigee (rad/s)
    times       - times at which ground track is plotted (s)
    ra          - vector of right ascensions of the spacecraft (deg)
    dec         - vector of declinations of the spacecraft (deg)
    TA          - true anomaly (rad)
    r           - perifocal position vector of satellite (km)
    R           - geocentric equatorial position vector (km)
    R1          - DCM for rotation about z through RA
    R2          - DCM for rotation about x through incl
    R3          - DCM for rotation about z through wp
    QxX         - DCM for rotation from perifocal to geocentric equatorial
    Q           - DCM for rotation from geocentric equatorial
                  into earth-fixed frame
    r_rel       - position vector in earth-fixed frame (km)
    alpha       - satellite right ascension (deg)
    delta       - satellite declination (deg)
    n_curves    - number of curves comprising the ground track plot
    RA          - cell array containing the right ascensions for each of
                  the curves comprising the ground track plot
    Dec         - cell array containing the declinations for each of
                  the curves comprising the ground track plot

    User h-functions required: sv_from_coe, kepler_E, ra_and_dec_from_r
*/

//...Constants
constexpr double PI = 3.141592653589793;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double MU = 398600.4418;
constexpr double J2 = 0.00108263;
constexpr double RE = 6378.14;
constexpr double WE = (2 * PI + 2 * PI / 365.26) / (24 * 3600.0);

std::vector<double> ra, dec;
std::vector<std::vector<double>> RA, Dec;
size_t n_curves;

void find_ra_and_dec(double to, double tf, double T, double h, double e, double Wdot, double wpdot,
                     double Wo, double wpo, double incl);
void form_separate_curves();
void plot_ground_track();
void print_orbital_data(double h, double e, double a, double rP, double rA, double T, double incl, double TAo,
                        double to, double Wo, double Wdot, double wpo, double wpdot, const std::vector<double>& r0, const std::vector<double>& v0);

int main() 
{
    //...Data declaration:
    double rP = 6700.0;
    double rA = 10000.0;
    double TAo = 230.0 * DEG_TO_RAD;
    double Wo = 270.0 * DEG_TO_RAD;
    double incl = 60.0 * DEG_TO_RAD;
    double wpo = 45.0 * DEG_TO_RAD;
    double n_periods = 3.25;
    //...End data declaration

    //...Compute the initial time (since perigee) and
    // the rates of node regression and perigee advance
    double a = (rA + rP) / 2.0;
    double T = 2.0 * PI / std::sqrt(MU) * std::pow(a, 1.5);
    double e = (rA - rP) / (rA + rP);
    double h = std::sqrt(MU * a * (1 - e * e));

    double Eo = 2.0 * std::atan(std::tan(TAo / 2.0) * std::sqrt((1 - e) / (1 + e)));
    double Mo = Eo - e * std::sin(Eo);
    double to = Mo * (T / (2.0 * PI));
    double tf = to + n_periods * T;

    double fac = -1.5 * std::sqrt(MU) * J2 * std::pow(RE, 2) / std::pow((1 - e * e), 2) / std::pow(a, 3.5);
    double Wdot = fac * std::cos(incl);
    double wpdot = fac * (2.5 * std::sin(incl) * std::sin(incl) - 2);

    find_ra_and_dec(to, tf, T, h, e, Wdot, wpdot, Wo, wpo, incl);

    form_separate_curves();

    plot_ground_track();

    std::vector<double> coe = {h, e, Wo, incl, wpo, TAo};
    auto result = sv_from_coe(coe, MU);
    std::vector<double> r0 = result.first;
    std::vector<double> v0 = result.second;

    // Print orbital data
    print_orbital_data(h, e, a, rP, rA, T, incl, TAo, to, Wo, Wdot, wpo, wpdot, r0, v0);

    return 0;
}

void find_ra_and_dec(double to, double tf, double T, double h, double e, double Wdot, double wpdot,
                     double Wo, double wpo, double incl) 
{
    /*
    Propagates the orbit over the specified time interval, transforming
    the position vector into the earth-fixed frame and, from that,
    computing the right ascension and declination histories
    */
    size_t num_points = 1000;
    double dt = (tf - to) / (num_points - 1);

    ra.clear();
    dec.clear();

    for (size_t i = 0; i < num_points; ++i) 
    {
        double t = to + i * dt;
        double M = 2.0 * PI / T * t;
        double E = kepler_E(e, M);
        double TA = 2.0 * std::atan(std::tan(E / 2.0) * std::sqrt((1 + e) / (1 - e)));
        std::vector<double> r = {h * h / MU / (1 + e * std::cos(TA)) * std::cos(TA),
                                 h * h / MU / (1 + e * std::cos(TA)) * std::sin(TA),
                                 0.0};

        double W = Wo + Wdot * t;
        double wp = wpo + wpdot * t;

        std::vector<std::vector<double>> R1 = {{std::cos(W), std::sin(W), 0},
                                               {-std::sin(W), std::cos(W), 0},
                                               {0, 0, 1}};
        std::vector<std::vector<double>> R2 = {{1, 0, 0},
                                               {0, std::cos(incl), std::sin(incl)},
                                               {0, -std::sin(incl), std::cos(incl)}};
        std::vector<std::vector<double>> R3 = {{std::cos(wp), std::sin(wp), 0},
                                               {-std::sin(wp), std::cos(wp), 0},
                                               {0, 0, 1}};

        std::vector<double> R(3, 0.0);
        for (size_t row = 0; row < 3; ++row) 
        {
            for (size_t col = 0; col < 3; ++col) 
            {
                R[row] += R1[row][col] * R2[col][0] * r[0] + R2[col][1] * r[1];
            }
        }

        double theta = WE * (t - to);
        std::vector<std::vector<double>> Q = {{std::cos(theta), std::sin(theta), 0},
                                              {-std::sin(theta), std::cos(theta), 0},
                                              {0, 0, 1}};

        std::vector<double> r_rel(3, 0.0);
        for (size_t row = 0; row < 3; ++row) 
        {
            for (size_t col = 0; col < 3; ++col) 
            {
                r_rel[row] += Q[row][col] * R[col];
            }
        }

        auto [alpha, delta] = ra_and_dec_from_r(r_rel);
        ra.push_back(alpha);
        dec.push_back(delta);
    }
}

void form_separate_curves() 
{
    /*
    Breaks the ground track up into separate curves which start
    and terminate at right ascensions in the range [0,360 deg].
    */
    double tol = 100.0;
    n_curves = 1;
    double ra_prev = ra[0];
    std::vector<double> current_ra, current_dec;

    for (size_t i = 0; i < ra.size(); ++i) {
        if (std::abs(ra[i] - ra_prev) > tol) {
            RA.push_back(current_ra);
            Dec.push_back(current_dec);
            current_ra.clear();
            current_dec.clear();
            n_curves++;
        }
        current_ra.push_back(ra[i]);
        current_dec.push_back(dec[i]);
        ra_prev = ra[i];
    }

    if (!current_ra.empty()) {
        RA.push_back(current_ra);
        Dec.push_back(current_dec);
    }
}

void print_orbital_data(double h, double e, double a, double rP, double rA, double T, double incl, double TAo,
                        double to, double Wo, double Wdot, double wpo, double wpdot, 
                        const std::vector<double>& r0, const std::vector<double>& v0) 
{
    std::cout << "\n----------------------------------------------\n";
    std::cout << "Angular momentum     = " << h << " km^2/s\n";
    std::cout << "Eccentricity         = " << e << "\n";
    std::cout << "Semimajor axis       = " << a << " km\n";
    std::cout << "Perigee radius       = " << rP << " km\n";
    std::cout << "Apogee radius        = " << rA << " km\n";
    std::cout << "Period               = " << T / 3600.0 << " hours\n";
    std::cout << "Inclination          = " << incl / DEG_TO_RAD << " deg\n";
    std::cout << "Initial true anomaly = " << TAo / DEG_TO_RAD << " deg\n";
    std::cout << "Time since perigee   = " << to / 3600.0 << " hours\n";
    std::cout << "Initial RA           = " << Wo / DEG_TO_RAD << " deg\n";
    std::cout << "RA_dot               = " << Wdot / DEG_TO_RAD * T << " deg/period\n";
    std::cout << "Initial wp           = " << wpo / DEG_TO_RAD << " deg\n";
    std::cout << "wp_dot               = " << wpdot / DEG_TO_RAD * T << " deg/period\n";
    std::cout << "\nr0 = [" << r0[0] << " " << r0[1] << " " << r0[2] << "]\n";
    std::cout << "magnitude = " << std::sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]) << " km\n";
    std::cout << "v0 = [" << v0[0] << " " << v0[1] << " " << v0[2] << "]\n";
    std::cout << "magnitude = " << std::sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]) << " km/s\n";

    std::cout << "----------------------------------------------\n";
}

void plot_ground_track() {
    std::cout << "\n----------------------------------------------\n";
    std::cout << "Ground Track Data:\n";
    for (size_t i = 0; i < RA.size(); ++i) {
        std::cout << "\nCurve " << i + 1 << ":\n";
        for (size_t j = 0; j < RA[i].size(); ++j) {
            std::cout << "RA: " << RA[i][j] << " deg, Dec: " << Dec[i][j] << " deg\n";
        }
    }
    std::cout << "----------------------------------------------\n";
}