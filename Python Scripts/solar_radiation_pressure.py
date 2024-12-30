import sys
import os
import importlib.util

std_lib_dir = os.path.dirname(os.__file__)  
std_bisect_path = os.path.join(std_lib_dir, 'bisect.py')
spec = importlib.util.spec_from_file_location("bisect", std_bisect_path)
real_bisect = importlib.util.module_from_spec(spec)
spec.loader.exec_module(real_bisect)
sys.modules["bisect"] = real_bisect

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fftpack import dct, idct
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
import sv_from_coe
import los
import solar_position

def solar_radiation_pressure():
    '''
    This function solves the Gauss planetary equations for
    solar radiation pressure (Equations 10.106).
    
    User M-functions required:  sv_from_coe, los, solar_position
    User subfunctions required: rates
    The M-function rsmooth may be found in Garcia, D: Robust Smoothing of Gridded Data in One
    and Higher Dimensions with Missing Values, Computational Statistics and Data Analysis,
    5 Vol. 54, 1167-1178, 2010.
    '''
    global JD0

    #...Conversion factors:
    hours = 3600                   # Hours to seconds
    days = 24 * hours              # Days to seconds
    deg = np.pi / 180              # Degrees to radians

    #...Constants:
    mu = 398600.4418               # Gravitational parameter (km^3/s^2)
    RE = 6378.14                   # Earth's radius (km)
    c = 2.998e8                    # Speed of light (m/s)
    S = 1367                       # Solar constant (W/m^2)
    Psr = S / c                    # Solar pressure (Pa)

    #...Satellite data:
    CR = 2                         # Radiation pressure coefficient
    m = 100                        # Mass (kg)
    As = 200                       # Frontal area (m^2)

    #...Initial orbital parameters (given):
    a0 = 10085.44                  # Semimajor axis (km)
    e0 = 0.025422                  # Eccentricity
    incl0 = 88.3924 * deg          # Inclination (radians)
    RA0 = 45.38124 * deg           # Right ascension of the node (radians)
    TA0 = 343.4268 * deg           # True anomaly (radians)
    w0 = 227.493 * deg             # Argument of perigee (radians)

    #...Initial orbital parameters (inferred):
    h0 = np.sqrt(mu * a0 * (1 - e0**2))    # Angular momentum (km^2/s)
    T0 = 2 * np.pi / np.sqrt(mu) * a0**1.5 # Period (s)
    rp0 = h0**2 / mu / (1 + e0)            # Perigee radius (km)
    ra0 = h0**2 / mu / (1 - e0)            # Apogee radius (km)

    #...Store initial orbital elements in the vector coe0:
    coe0 = [h0, e0, RA0, incl0, w0, TA0]

    #...Use solve_ivp to integrate Equations 12.106, the Gauss planetary equations
    #   from t0 to tf:
    JD0 = 2438400.5                # Initial Julian date (6 January 1964 0 UT)
    t0 = 0                         # Initial time (s)
    tf = 3 * 365 * days            # Final time (s)
    nout = 4000                    # Number of solution points to output
    t_eval = np.linspace(t0, tf, nout)

    def rates(t, f):
        # Update the Julian Date at time t:
        JD = JD0 + t / days

        # Compute the apparent position vector of the sun:
        lamda, eps, r_sun = solar_position.solar_position(JD)

        # Convert to radians:
        lamda = lamda * deg
        eps = eps * deg

        # Extract orbital elements:
        h, e, RA, i, w, TA = f
        u = w + TA  # Argument of latitude

        # Compute the state vector:
        coe = [h, e, RA, i, w, TA]
        R, V = sv_from_coe.sv_from_coe(coe, mu)

        # Calculate the magnitude of the radius vector:
        r = np.linalg.norm(R)

        # Compute shadow function and solar radiation perturbation:
        nu = los.los(R, r_sun)
        pSR = nu * (S / c) * CR * As / m / 1000

        # Calculate trigonometric functions:
        sl = np.sin(lamda)
        cl = np.cos(lamda)
        se = np.sin(eps)
        ce = np.cos(eps)
        sW = np.sin(RA)
        cW = np.cos(RA)
        si = np.sin(i)
        ci = np.cos(i)
        su = np.sin(u)
        cu = np.cos(u)
        sT = np.sin(TA)
        cT = np.cos(TA)

        # Earth-sun unit vector components:
        ur = sl * ce * cW * ci * su + sl * ce * sW * cu - cl * sW * ci * su + cl * cW * cu + sl * se * si * su
        us = sl * ce * cW * ci * cu - sl * ce * sW * su - cl * sW * ci * cu - cl * cW * su + sl * se * si * cu
        uw = -sl * ce * cW * si + cl * sW * si + sl * se * ci

        # Rates of change:
        hdot = -pSR * r * us
        edot = -pSR * (h / mu * sT * ur + 1 / mu / h * ((h**2 + mu * r) * cT + mu * e * r) * us)
        TAdot = h / r**2 - pSR / e / h * (h**2 / mu * cT * ur - (r + h**2 / mu) * sT * us)
        RAdot = -pSR * r / h / si * su * uw
        idot = -pSR * r / h * cu * uw
        wdot = -pSR * (-1 / e / h * (h**2 / mu * cT * ur - (r + h**2 / mu) * sT * us) - r * su / h / si * ci * uw)

        return [hdot, edot, RAdot, idot, wdot, TAdot]

    sol = solve_ivp(rates, [t0, tf], coe0, t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-8)

    #...Extract the solution:
    t = sol.t
    y = sol.y
    h, e, RA, incl, w, TA = y
    a = h**2 / mu / (1 - e**2)

    #...Smooth the data:
    h = rsmooth(h)
    e = rsmooth(e)
    RA = rsmooth(RA)
    incl = rsmooth(incl)
    w = rsmooth(w)
    a = rsmooth(a)

    #...Plot the results:
    plt.figure(figsize=(10, 12))

    plt.subplot(3, 2, 1)
    plt.plot(t / days, h - h0)
    plt.title('Angular Momentum (km^2/s)')
    plt.xlabel('Days')

    plt.subplot(3, 2, 2)
    plt.plot(t / days, e - e0)
    plt.title('Eccentricity')
    plt.xlabel('Days')

    plt.subplot(3, 2, 3)
    plt.plot(t / days, a - a0)
    plt.title('Semimajor axis (km)')
    plt.xlabel('Days')

    plt.subplot(3, 2, 4)
    plt.plot(t / days, (RA - RA0) / deg)
    plt.title('Right Ascension (deg)')
    plt.xlabel('Days')

    plt.subplot(3, 2, 5)
    plt.plot(t / days, (incl - incl0) / deg)
    plt.title('Inclination (deg)')
    plt.xlabel('Days')

    plt.subplot(3, 2, 6)
    plt.plot(t / days, (w - w0) / deg)
    plt.title('Argument of Perigee (deg)')
    plt.xlabel('Days')

    plt.tight_layout()
    plt.show()

def rsmooth(y):
    '''
    Apply recursive smoothing to the data.
    '''
    y = np.asarray(y)
    n = len(y)
    Lambda = -2 + 2 * np.cos(np.arange(n) * np.pi / n)
    Lambda_sq = Lambda**2 
    W = np.ones(n)
    z = y.copy()

    for _ in range(6):
        tol = float('inf')
        while tol > 1e-5:
            DCTy = dct(W * (y - z) + z, norm='ortho')
            GCVscore = lambda p: np.sum((W * (y - idct(1 / (1 + 10**p * Lambda_sq) * DCTy, norm='ortho')))**2)
            result = minimize_scalar(GCVscore, bounds=(-15, 38), method='bounded')
            s = 10**result.x
            Gamma = 1 / (1 + s * Lambda_sq)
            new_z = idct(Gamma * DCTy, norm='ortho') 
            tol = np.linalg.norm(new_z - z) / max(np.linalg.norm(z), 1e-10)
            z = new_z
    return z

if __name__ == '__main__':
    solar_radiation_pressure()
