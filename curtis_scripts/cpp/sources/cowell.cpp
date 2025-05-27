#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>  
#include "sv_from_coe.h"   
#include "atmosphere.h"
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
using namespace std;

/*
    This function uses MATLAB's ode45 to numerically
    integrate Equation 10.2 for atmospheric drag.
    
    User h-functions required:  sv_from_coe, atmosphere
    User subfunctions required: rates, terminate    
*/
static double alt = 0.0;

vector<double> rates(double t, const vector<double> &f);
bool terminate(double currentAlt);

vector<double> rk4Step(double t, double dt, const vector<double> &y)
{
    vector<double> k1 = rates(t, y);

    vector<double> yk2(y.size());
    for (size_t i = 0; i < y.size(); i++)
        yk2[i] = y[i] + 0.5 * dt * k1[i];

    vector<double> k2 = rates(t + 0.5 * dt, yk2);

    vector<double> yk3(y.size());
    for (size_t i = 0; i < y.size(); i++)
        yk3[i] = y[i] + 0.5 * dt * k2[i];

    vector<double> k3 = rates(t + 0.5 * dt, yk3);

    vector<double> yk4(y.size());
    for (size_t i = 0; i < y.size(); i++)
        yk4[i] = y[i] + dt * k3[i];

    vector<double> k4 = rates(t + dt, yk4);

    vector<double> yout(y.size());
    for (size_t i = 0; i < y.size(); i++)
    {
        yout[i] = y[i] + (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    return yout;
}

int main()
{
    //...Conversion factors:
    double hours = 3600.0;                          // Hours to seconds
    double days  = 24.0 * hours;                    // Days to seconds
    double deg   = M_PI / 180.0;                    // Degrees to radians

    //...Constants:
    double mu    = 398600.0;                        // Gravitational parameter (km^3/s^2)
    double RE    = 6378.0;                          // Earth's radius (km)
    vector<double> wE = {0.0, 0.0, 7.2921159e-5};   // Earth's angular velocity (rad/s)
     
    //...Satellite data:
    double CD    = 2.2;                             // Drag coefficient
    double m     = 100.0;                           // Mass (kg)
    double A     = M_PI/4.0 * pow(1.0, 2.0);        // Frontal area (m^2)

    //...Initial orbital parameters (given):
    double rp    = RE + 215.0;                      // Perigee radius (km)
    double ra    = RE + 939.0;                      // Apogee radius (km)   
    double RA    = 339.94*deg;                      // Right ascension of the node (radians)
    double i     = 65.1*deg;                        // Inclination (radians)
    double w     = 58.0*deg;                        // Argument of perigee (radians)
    double TA    = 332.0*deg;                       // True anomaly (radians)

    //...Initial orbital parameters (inferred):
    double e     = (ra - rp) / (ra + rp);           // Eccentricity
    double a     = (rp + ra) / 2.0;                 // Semi-major axis (km)
    double h     = sqrt(mu * a * (1.0 - e*e));      // Angular momentum (km^2/s)
    double T     = 2.0*M_PI / sqrt(mu)*pow(a, 1.5); // Period (s)

    //...Store initial orbital elements (from above) in the vector coe0:
    vector<double> coe0(6);
    coe0[0] = h;
    coe0[1] = e;
    coe0[2] = RA;
    coe0[3] = i;
    coe0[4] = w;
    coe0[5] = TA;

    //...Obtain the initial state vector from Algorithm 4.5 (sv_from_coe):
    auto result = sv_from_coe(coe0, mu);
    std::vector<double> R0 = result.first;  // R0 is the initial position vector
    std::vector<double> V0 = result.second; // V0 is the initial velocity vector

    // Magnitudes of R0 and V0
    double r0 = sqrt(R0[0]*R0[0] + R0[1]*R0[1] + R0[2]*R0[2]);
    double v0 = sqrt(V0[0]*V0[0] + V0[1]*V0[1] + V0[2]*V0[2]);

    //...Use ODE to integrate the equations of motion d/dt(R,V) = f(R,V) 
    //   from t0 to tf:
    double t0    = 0.0;           // Initial time
    double tf    = 120.0 * days;  // Final time
    int    nout  = 40000;         // Number of solution points to output    
    vector<double> tspan(nout);
    double dt = (tf - t0) / (nout - 1);
    for (int i = 0; i < nout; i++)
    {
        tspan[i] = t0 + i * dt; // Integration time interval
    }

    vector<double> y0(6); // Initial state vector
    y0[0] = R0[0]; y0[1] = R0[1]; y0[2] = R0[2];
    y0[3] = V0[0]; y0[4] = V0[1]; y0[5] = V0[2];

    vector< vector<double> > yHist(nout, vector<double>(6));
    vector<double> altHist(nout, 0.0);

    //   Set error tolerances, initial step size, and termination event:
    yHist[0]   = y0;
    altHist[0] = sqrt(y0[0]*y0[0] + y0[1]*y0[1] + y0[2]*y0[2]) - RE;

    bool terminateIntegration = false;
    for(int iStep = 0; iStep < nout-1; iStep++)
    {
        double currentTime = tspan[iStep];
        vector<double> currentY = yHist[iStep];

        if (terminate(altHist[iStep])) 
        {
            terminateIntegration = true;
            break;
        }

        vector<double> nextY = rk4Step(currentTime, dt, currentY);

        yHist[iStep+1] = nextY;
        double rMag = sqrt(nextY[0]*nextY[0] + nextY[1]*nextY[1] + nextY[2]*nextY[2]);
        altHist[iStep+1] = rMag - RE;
    }

    //...Extract the locally extreme altitudes:
    vector<int> isMax(nout, 0); // Logical array for maxima
    vector<int> isMin(nout, 0); // Logical array for minima

    for(int i = 1; i < nout-1; i++)
    {
        if (altHist[i] > altHist[i-1] && altHist[i] > altHist[i+1])
        {
            isMax[i] = 1;
        }
        if (altHist[i] < altHist[i-1] && altHist[i] < altHist[i+1])
        {
            isMin[i] = 1;
        }
    }

    vector<double> tMax, altMax;
    vector<double> tMin, altMin;
    for(int i = 1; i < nout-1; i++)
    {
        if(isMax[i] == 1)
        {
            tMax.push_back(tspan[i]);
            altMax.push_back(altHist[i]);
        }
        if(isMin[i] == 1)
        {
            tMin.push_back(tspan[i]);
            altMin.push_back(altHist[i]);
        }
    }

    cout << fixed << setprecision(6);
    cout << "\nApogee (max altitude):\n";
    for (size_t k = 0; k < tMax.size(); k++)
    {
        cout << "  Time (days) = " 
             << (tMax[k] / days) << "  Altitude (km) = " << altMax[k] << "\n";
    }

    cout << "\nPerigee (min altitude):\n";
    for (size_t k = 0; k < tMin.size(); k++)
    {
        cout << "  Time (days) = "
             << (tMin[k] / days) << "  Altitude (km) = " << altMin[k] << "\n";
    }
    return 0;
}

vector<double> rates(double /*t*/, const vector<double> &f)
{
    /*
        This function calculates the spacecraft acceleration from its
        position and velocity at time t.
    */
    static const double RE  = 6378.14;    
    static const double mu  = 398600.4418;      
    static const double CD  = 2.2;         
    static const double m   = 100.0;      
    static const double A   = M_PI/4.0 * 1.0*1.0; 
    static const vector<double> wE = {0.0, 0.0, 7.2921159e-5};

    vector<double> R(3), V(3);
    R[0] = f[0]; R[1] = f[1]; R[2] = f[2];
    V[0] = f[3]; V[1] = f[4]; V[2] = f[5];

    double r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);  
    alt = r - RE;                                        

    double rho = atmosphere(alt);

    vector<double> crossWE_R(3);
    crossWE_R[0] = wE[1]*R[2] - wE[2]*R[1];
    crossWE_R[1] = wE[2]*R[0] - wE[0]*R[2];
    crossWE_R[2] = wE[0]*R[1] - wE[1]*R[0];

    vector<double> Vrel(3);
    Vrel[0] = V[0] - crossWE_R[0];
    Vrel[1] = V[1] - crossWE_R[1];
    Vrel[2] = V[2] - crossWE_R[2];

    double vrel = sqrt(Vrel[0]*Vrel[0] + Vrel[1]*Vrel[1] + Vrel[2]*Vrel[2]);

    vector<double> uv(3);
    if (vrel > 1e-12)
    {
        uv[0] = Vrel[0] / vrel;
        uv[1] = Vrel[1] / vrel;
        uv[2] = Vrel[2] / vrel;
    }
    else
    {
        uv[0] = uv[1] = uv[2] = 0.0; 
    }

    double factor = -CD*A/m * 0.5 * rho * (1000.0 * vrel)*(1000.0 * vrel);

    vector<double> ap(3);
    ap[0] = (factor * uv[0]) / 1000.0;
    ap[1] = (factor * uv[1]) / 1000.0;
    ap[2] = (factor * uv[2]) / 1000.0;

    double a0factor = -mu / (r*r*r);
    vector<double> a0(3);
    a0[0] = a0factor * R[0];
    a0[1] = a0factor * R[1];
    a0[2] = a0factor * R[2];

    vector<double> a(3);
    a[0] = a0[0] + ap[0];
    a[1] = a0[1] + ap[1];
    a[2] = a0[2] + ap[2];

    vector<double> dfdt(6);
    dfdt[0] = V[0]; 
    dfdt[1] = V[1];
    dfdt[2] = V[2];
    dfdt[3] = a[0];
    dfdt[4] = a[1];
    dfdt[5] = a[2];

    return dfdt;
}

bool terminate(double currentAlt)
{
    /*
        This function specifies the event at which ode terminates.
    */
    if (currentAlt <= 100.0)
        return true;
    return false;
}