#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <functional>
#include "rkf45.h"
#include "ra_and_dec_from_r.h"
#include "simpsons_lunar_ephemeris.h"

/*
    This program presents the graphical solution of the motion of a
    spacecraft in the gravity fields of both the earth and the moon for
    the initial data provided in the input declaration below.
    Header file rkf45.h Runge-Kutta solver is used.

    deg                         - conversion factor, degrees to radians
    days                        - conversion factor, days to seconds
    Re, Rm                      - radii of earth and moon, respectively (km)
    m_e, m_m                    - masses of earth and moon, respectively (kg)
    mu_e, mu_m                  - gravitational parameters of earth and moon,
                                  respectively (km^3/s^2)
    D                           - semimajor axis of moon's orbit (km)
    I_, J_, K_                  - unit vectors of ECI frame
    RS                          - radius of moon's sphere of influence (km)
                                  year, month, hour, minute
    second                      - Date and time of spacecraft's lunar arrival
    t0                          - initial time on the trajectory (s)
    z0                          - initial altitude of the trajectory (km)
    alpha0, dec0                - initial right ascension and declination of
                                  spacecraft (deg)
    gamma0                      - initial flight path angle (deg)
    fac                         - ratio of spaccraft's initial speed to the
                                  escape speed.
    ttt                         - predicted time to perilune (s)
    tf                          - time at end of trajectory (s)
    jd0                         - julian date of lunar arrival
    rm0, vm0                    - state vector of the moon at jd0 (km, km/s)
    RA, Dec                     - right ascension and declination of the moon
                                  at jd0 (deg)
    hmoon_, hmoon               - moon's angular momentum vector and magnitude
                                  at jd0 (km^2/s)
    inclmoon                    - inclination of moon's orbit earth's
                                  equatorial plane (deg)
    r0                          - initial radius from earth's center to
                                  probe (km)
    r0_                         - initial ECI position vector of probe (km)
    vesc                        - escape speed at r0 (km/s)
    v0                          - initial ECI speed of probe (km/s)
    w0_                         - unit vector normal to plane of translunar
                                  orbit at time t0
    ur_                         - radial unit vector to probe at time t0
    uperp_                      - transverse unit vector at time t0
    vr                          - initial radial speed of probe (km/s)
    vperp                       - initial transverse speed of probe (km/s)
    v0_                         - initial velocity vector of probe (km/s)
    uv0_                        - initial tangential unit vector
    y0                          - initial state vector of the probe (km, km/s)
    t                           - vector containing the times from t0 to tf at
                                  which the state vector is evaluated (s)
    y                           - a matrix whose 6 columns contain the inertial
                                  position and velocity components evaluated
                                  at the times t(:) (km, km/s)
    X, Y, Z                     - the probe's inertial position vector history
    vX, vY, VZ                  - the probe's inertial velocity history
    x, y, z                     - the probe's position vector history in the
                                  Moon-fixed frame
    Xm, Ym, Zm                  - the Moon's inertial position vector history
    vXm, vYm, vZm               - the Moon's inertial velocity vector history
    ti                          - the ith time of the set [t0,tf] (s)
    r_                          - probe's inertial position vector at time ti
                                  (km)
    r                           - magnitude of r_ (km)
    jd                          - julian date of corresponding to ti (days)
    rm_, vm_                    - the moon's state vector at time ti (km,km/s)
    x_, y_, z_                  - vectors along the axes of the rotating
                                  moon-fixed at time ti (km)
    i_, j_, k_                  - unit vectors of the moon-fixed rotating frame
                                  at time ti
    Q                           - DCM of transformation from ECI to moon-fixed
                                  frame at time ti
    rx_                         - probe's inertial position vector in moon-
                                  fixed coordinates at time ti (km)
    rmx_                        - Moon's inertial position vector in moon-
                                  fixed coordinates at time ti (km)
    dist_                       - position vector of probe relative to the moon
                                  at time ti (km)
    dist                        - magnitude of dist_ (km)
    dist_min                    - perilune of trajectory (km)
    rmTLI_                      - Moon's position vector at TLI
    RATLI, DecTLI               - Moon's right ascension and declination at
                                  TKI (deg)
    v_atdmin_                   - Probe's velocity vector at perilune (km/s)
    rm_perilume, vm_perilune    - Moon's state vector when the probe is at
                                  perilune (km, km/s)
    rel_speed                   - Speed of probe relative to the Moon at
                                  perilune (km/s)
    RA_at_perilune              - Moon's RA at perilune arrival (deg)
    Dec_at_perilune             - Moon's Dec at perilune arrival (deg)
    target_error                - Distance between Moon's actual position at
                                  perilune arrival and its position after the
                                  predicted flight time, ttt (km).
    rms_                        - position vector of moon relative to
                                  spacecraft (km)
    rms                         - magnitude of rms_ (km)
    aearth_                     - acceleration of spacecraft due to
                                  earth (km/s^2)
    amoon_                      - acceleration of spacecraft due to
                                  moon (km/s^2)
    atot_                       - aearth_ + amoon_ (km/s^2)
    binormal_                   - unit vector normal to the osculating plane
    incl                        - angle between inertial Z axis and the
                                  binormal (deg)
    rend_                       - Position vector of end point of trajectory
                                  (km)
    alt_end                     - Altitude of end point of trajectory (km)
    ra_end, dec_end             - Right ascension and declination of end point
                                  of trajectory (km)

    User h-functions required: none
    User subfunctions required: rates
*/
//...general data
static constexpr double PI      = 3.14159265358979323846;
static constexpr double DEG2RAD = PI/180.0;
static constexpr double RAD2DEG = 180.0/PI;
static constexpr double DAY2SEC = 24.0*3600.0;
static constexpr double Re      = 6378.0;   
static constexpr double Rm      = 1737.0;     
static constexpr double mu_e    = 398600.4; 
static constexpr double mu_m    = 4902.8;    
double jd0;
double ttt; 

// Vector algebra utilities
using Vec3      = std::array<double,3>;
using state_type = std::vector<double>;

double norm(const Vec3 &v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
Vec3 operator+(const Vec3 &a, const Vec3 &b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}
Vec3 operator-(const Vec3 &a, const Vec3 &b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
Vec3 operator*(double s, const Vec3 &v) {
    return {s*v[0], s*v[1], s*v[2]};
}
Vec3 operator/(const Vec3 &v, double s) {
    return {v[0]/s, v[1]/s, v[2]/s};
}
Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

//...State vector of moon at target date:
double julian_date(int year, int month, int day,
                   int hour, int minute, double second)
{
    if(month <= 2) {
        year  -= 1;
        month += 12;
    }
    int A = year/100;
    int B = 2 - A + A/4;
    double day_frac = (hour + minute/60.0 + second/3600.0)/24.0;
    double jd = std::floor(365.25*(year + 4716))
              + std::floor(30.6001*(month + 1))
              + day + B - 1524.5 + day_frac;
    return jd;
}

//------------------------------------------------------------------------------
// The RHS of the restricted three-body problem
std::vector<double>
ode_fun(double t, const std::vector<double> &y)
{
    // Unpack state
    Vec3 r = { y[0], y[1], y[2] };
    Vec3 v = { y[3], y[4], y[5] };
    double r_norm = norm(r);

    // Current Julian date
    double jd = jd0 - (ttt - t)/DAY2SEC;

    // Moon state
    Vec3 rm, vm;
    simpsons_lunar_ephemeris(jd, rm, vm);

    double rm_norm  = norm(rm);
    Vec3   rms_vec  = rm - r;
    double rms_norm = norm(rms_vec);

    // Accelerations
    Vec3 a_e = (-mu_e/(r_norm*r_norm*r_norm)) * r;
    Vec3 a_m = mu_m * ( rms_vec/(rms_norm*rms_norm*rms_norm)
                      - rm/(rm_norm*rm_norm*rm_norm) );
    Vec3 a   = a_e + a_m;

    // Pack derivative
    std::vector<double> dydt(6);
    dydt[0] = v[0];
    dydt[1] = v[1];
    dydt[2] = v[2];
    dydt[3] = a[0];
    dydt[4] = a[1];
    dydt[5] = a[2];
    return dydt;
}

//------------------------------------------------------------------------------
// Main program
int main()
{
    //...Data declaration
    const std::string Title = "Lunar Trajectory";
    //   Date and time of lunar arrival:
    int    year   = 2020;
    int    month  = 5;
    int    day    = 4;
    int    hour   = 12;
    int    minute = 0;
    double second = 0.0;
    double t0     = 0.0;              
    double z0     = 320.0;     
    double alpha0 = 90.0 * DEG2RAD;
    double dec0   = 15.0 * DEG2RAD;
    double gamma0 = 40.0 * DEG2RAD;
    double fac    = 0.9924;              // Fraction of Vesc
    ttt           = 3 * DAY2SEC;
    double tf     = ttt + 2.667 * DAY2SEC;
    //...End data declaration

    //...State vector of moon at target date:
    jd0 = julian_date(year, month, day, hour, minute, second);
    Vec3 rm0, vm0;
    simpsons_lunar_ephemeris(jd0, rm0, vm0);

    auto [RA, Dec]   = ra_and_dec_from_r({rm0[0], rm0[1], rm0[2]});
    double distance  = norm(rm0);
    Vec3   hmoon_vec = cross(rm0, vm0);
    double hmoon_mag = norm(hmoon_vec);
    double inclmoon  = std::acos(hmoon_vec[2]/hmoon_mag)*RAD2DEG;

    //...Initial position vector of probe:
    double r0_mag = Re + z0;
    Vec3   I_     = {1,0,0}, J_ = {0,1,0}, K_ = cross(I_, J_);
    Vec3 r0_vec   = r0_mag * Vec3{std::cos(alpha0)*std::cos(dec0),
                                  std::sin(alpha0)*std::cos(dec0),
                                  std::sin(dec0)};
    double vesc  = std::sqrt(2*mu_e/r0_mag);
    double v0    = fac * vesc;
    Vec3 w0_vec  = cross(r0_vec, rm0);
    w0_vec       = w0_vec / norm(w0_vec);
    Vec3 ur      = r0_vec / norm(r0_vec);
    Vec3 uperp   = cross(w0_vec, ur);
    uperp        = uperp / norm(uperp);
    double vr    = v0 * std::sin(gamma0);
    double vperp = v0 * std::cos(gamma0);
    Vec3 v0_vec  = vr*ur + vperp*uperp;

    state_type y0(6);
    y0[0]=r0_vec[0]; y0[1]=r0_vec[1]; y0[2]=r0_vec[2];
    y0[3]=v0_vec[0]; y0[4]=v0_vec[1]; y0[5]=v0_vec[2];

    auto sol = rkf45(
        std::function<std::vector<double>(double,const std::vector<double>&)>(
            ode_fun), std::make_pair(t0, tf), y0, 1e-10);

    //...Compute the Moon's trajectory from an ephemeris, find perilune of the
    //   probe's trajectory, and project the probe's trajectory onto the axes
    //   of the Moon-fixed rotating frame:
    double dist_min = 1e30;
    size_t imin     = 0;

    for(size_t i=0; i<sol.tout.size(); ++i) {
        double t = sol.tout[i];
        Vec3   r = { sol.yout[i][0], sol.yout[i][1], sol.yout[i][2] };

        double jd_i = jd0 - (ttt - t)/DAY2SEC;
        Vec3   rm, vm;
        simpsons_lunar_ephemeris(jd_i, rm, vm);

        double d = norm(rm - r);
        if(d < dist_min) {
            dist_min = d;
            imin     = i;
        }
    }

    //...State vector and celestial position of moon when probe is at perilune:
    auto &y_per = sol.yout[imin];
    Vec3 r_per = { y_per[0], y_per[1], y_per[2] };
    Vec3 v_per = { y_per[3], y_per[4], y_per[5] };
    double t_per = sol.tout[imin];

    double jd_per = jd0 - (ttt - t_per)/DAY2SEC;
    Vec3   rm_per, vm_per;
    simpsons_lunar_ephemeris(jd_per, rm_per, vm_per);
    auto [RA_per, Dec_per] = ra_and_dec_from_r({rm_per[0], rm_per[1], rm_per[2]});

    Vec3 v_rel = { v_per[0]-vm_per[0], v_per[1]-vm_per[1], v_per[2]-vm_per[2] };

    //..Speed of probe relative to Moon at perilune:
    double rel_speed = norm(v_rel);

    //...End point of trajectory:
    auto &y_end = sol.yout.back();
    Vec3 rend = { y_end[0], y_end[1], y_end[2] };
    double alt_end = norm(rend) - Re;
    auto [RA_end, Dec_end] = ra_and_dec_from_r({rend[0], rend[1], rend[2]});

    //...Output:
    std::cout<<std::fixed<<std::setprecision(6)
             << "\n\n" << Title << "\n\n"
             << "Date and time of arrival at Moon: "
             << month<<"/"<<day<<"/"<<year<<" "
             <<hour<<":"<<minute<<":"<<int(second)<<"\n\n"

             << "Moon's position:\n"
             << "  Distance               = " << distance << " km\n"
             << "  Right Ascension        = " << RA       << " deg\n"
             << "  Declination            = " << Dec      << " deg\n"
             << "  Orbital inclination    = " << inclmoon << " deg\n\n"

             << "Probe at Earth departure (t = " << t0 << " s):\n"
             << "  Altitude               = " << z0       << " km\n"
             << "  RA                     = " << alpha0*RAD2DEG << " deg\n"
             << "  Dec                    = " << dec0*RAD2DEG   << " deg\n"
             << "  Flight-path angle      = " << gamma0*RAD2DEG << " deg\n"
             << "  Speed                  = " << v0        << " km/s\n"
             << "  Escape speed           = " << vesc      << " km/s\n"
             << "  v/vesc                 = " << v0/vesc   << "\n\n"

             << "Probe at perilune (i = " << imin << "):\n"
             << "  Altitude above Moon    = " << (dist_min - Rm) << " km\n"
             << "  Speed                  = " << norm(v_per)    << " km/s\n"
             << "  Relative speed         = " << rel_speed      << " km/s\n"
             << "  Moon RA at perilune    = " << RA_per << " deg\n"
             << "  Moon Dec at perilune   = " << Dec_per<< " deg\n"
             << "  Time from TLI to perilune = "<< t_per/3600.0 <<" hours\n\n"

             << "Final state at t = " << tf << " s:\n"
             << "  Earth altitude         = " << alt_end << " km\n"
             << "  RA                     = " << RA_end  << " deg\n"
             << "  Dec                    = " << Dec_end << " deg\n\n";

    return 0;
}
