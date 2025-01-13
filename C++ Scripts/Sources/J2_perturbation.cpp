#include <iostream>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846

/*
    This function use MATLAB's ode45 to numerically integrate Equations 10.89
    (the Gauss planetary equations) to determine the J2 perturbation of 
    the orbital elements.

    User h-functions required:  None
    User subfunctions required: rates     
*/
//...Constants:
const double mu = 398600;     // Gravitational parameter (km^3/s^2)
const double RE = 6378;       // Earth's radius (km)
const double J2 = 1082.63e-6; // Earth's J2

struct J2_Perturbation 
{
    double h, e, RA, i, w, TA;
};

std::vector<double> linspace(double start, double end, int num) 
{
    std::vector<double> linspaced;
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) 
    {
        linspaced.push_back(start + delta * i);
    }
    return linspaced;
}

J2_Perturbation rates(double t, const J2_Perturbation& coe) 
{
    /*
        This function calculates the time rates of the orbital elements
        from Gauss's variational equations (Equations 12.89).
    */
    J2_Perturbation dfdt;

    //...The orbital elements at time t:
    double h = coe.h;
    double e = coe.e;
    double RA = coe.RA;
    double i = coe.i;
    double w = coe.w;
    double TA = coe.TA;

    double r = h * h / mu / (1 + e * cos(TA));  // The radius
    double u = w + TA;                          // Argument of latitude                   

    //...Orbital element rates at time t (Equations 12.89):
    dfdt.h = -3.0 / 2.0 * J2 * mu * RE * RE / 
             pow(r, 3) * sin(i) * sin(i) * sin(2 * u);

    dfdt.e = 3.0 / 2.0 * J2 * mu * RE * RE / h /
             pow(r, 3) * (h * h / mu / r * (sin(u) * sin(i) * sin(i) *
             (3 * sin(TA) * sin(u) - 2 * cos(TA) * cos(u)) - sin(TA)) -
             sin(i) * sin(i) * sin(2 * u) * (e + cos(TA)));

    dfdt.e = 3.0 / 2.0 * J2 * mu * RE * RE / h /
             pow(r, 3) * (h * h / mu / r * sin(TA) * (3 * sin(i) * sin(i) * sin(u) * sin(u) - 1)
             - sin(2 * u) * sin(i) * sin(i) * ((2 + e * cos(TA)) * cos(TA) + e));

    dfdt.TA = h / (r * r) + 3.0 / 2.0 * J2 * mu * RE * RE / e / h / 
             pow(r, 3) * (h * h / mu / r * cos(TA) * (3 * sin(i) * sin(i) * sin(u) * sin(u) - 1) 
             + sin(2 * u) * sin(i) * sin(i) * sin(TA) * (h * h / mu / r + 1));

    dfdt.RA = -3.0 * J2 * mu * RE * RE / h / pow(r, 3) * sin(u) * sin(u) * cos(i);

    dfdt.i = -3.0 / 4.0 * J2 * mu * RE * RE / h / pow(r, 3) * sin(2 * u) * sin(2 * i);

    dfdt.w = 3.0 / 2.0 * J2 * mu * RE * RE / e / h / pow(r, 3) *
             (-h * h / mu / r * cos(TA) * (3 * sin(i) * sin(i) * sin(u) *
             sin(u) - 1) - sin(2 * u) * sin(i) * sin(i) * sin(TA) * (2 + e * cos(TA))
             + 2 * e * cos(i) * cos(i) * sin(u) * sin(u));

    return dfdt;
}

void integrate(double t0, double tf, int nout, const J2_Perturbation& coe0) 
{
    std::vector<double> tspan = linspace(t0, tf, nout);
    std::vector<J2_Perturbation> y(nout);
    y[0] = coe0;

    double dt = (tf - t0) / (nout - 1);
    for (int i = 1; i < nout; ++i) {
        J2_Perturbation k1 = rates(tspan[i - 1], y[i - 1]);
        J2_Perturbation k2 = rates(tspan[i - 1] + dt / 2, {y[i - 1].h + k1.h * dt / 2, y[i - 1].e + k1.e * dt / 2, y[i - 1].RA + k1.RA * dt / 2, y[i - 1].i + k1.i * dt / 2, y[i - 1].w + k1.w * dt / 2, y[i - 1].TA + k1.TA * dt / 2});
        J2_Perturbation k3 = rates(tspan[i - 1] + dt / 2, {y[i - 1].h + k2.h * dt / 2, y[i - 1].e + k2.e * dt / 2, y[i - 1].RA + k2.RA * dt / 2, y[i - 1].i + k2.i * dt / 2, y[i - 1].w + k2.w * dt / 2, y[i - 1].TA + k2.TA * dt / 2});
        J2_Perturbation k4 = rates(tspan[i - 1] + dt, {y[i - 1].h + k3.h * dt, y[i - 1].e + k3.e * dt, y[i - 1].RA + k3.RA * dt, y[i - 1].i + k3.i * dt, y[i - 1].w + k3.w * dt, y[i - 1].TA + k3.TA * dt});

        //...Assign the time histories mnemonic variable names:
        y[i].h = y[i - 1].h + dt / 6 * (k1.h + 2 * k2.h + 2 * k3.h + k4.h);
        y[i].e = y[i - 1].e + dt / 6 * (k1.e + 2 * k2.e + 2 * k3.e + k4.e);
        y[i].RA = y[i - 1].RA + dt / 6 * (k1.RA + 2 * k2.RA + 2 * k3.RA + k4.RA);
        y[i].i = y[i - 1].i + dt / 6 * (k1.i + 2 * k2.i + 2 * k3.i + k4.i);
        y[i].w = y[i - 1].w + dt / 6 * (k1.w + 2 * k2.w + 2 * k3.w + k4.w);
        y[i].TA = y[i - 1].TA + dt / 6 * (k1.TA + 2 * k2.TA + 2 * k3.TA + k4.TA);

        //...Ouput the time histories of the osculating elements:
        std::cout << "t = " << tspan[i] << ", h = " << y[i].h << ", e = " << y[i].e << ", RA = " << y[i].RA << ", i = " << y[i].i << ", w = " << y[i].w << ", TA = " << y[i].TA << std::endl;
    }
}

int main() 
{
    //...Constants:
    double hours = 3600;                        // Hours to seconds
    double days = 24 * hours;                   // Days to seconds
    double deg = M_PI / 180;                    // Degrees to radians

    //...Initial orbital parameters (given):
    double rp0 = RE + 300;                      // Perigee radius (km)
    double ra0 = RE + 3062;                     // Apogee radius (km)
    double RA0 = 45 * deg;                      // Right ascension of the node (radians)
    double i0 = 28 * deg;                       // Inclination (radians)
    double w0 = 30 * deg;                       // Argument of perigee (radians)
    double TA0 = 40 * deg;                      // True anomaly (radians)

    //...Initial orbital parameters (inferred):
    double e0 = (ra0 - rp0) / (ra0 + rp0);      // Eccentricity
    double h0 = sqrt(rp0 * mu * (1 + e0));      // Angular momentum (km^2/s)
    double a0 = (rp0 + ra0) / 2;                // Semimajor axis (km)
    double T0 = 2 * M_PI/sqrt(mu)*pow(a0, 1.5); // Period (s)

    //...Store initial orbital elements (from above) in the vector coe0:
    J2_Perturbation coe0 = {h0, e0, RA0, i0, w0, TA0};

    //...Use ODE to integrate the Gauss variational equations (Equations
    //   12.89) from t0 to tf:
    double t0 = 0;
    double tf = 2 * days;
    int nout = 5000; // Number of solution points to output for plotting purposes
    integrate(t0, tf, nout, coe0);

    return 0;
}