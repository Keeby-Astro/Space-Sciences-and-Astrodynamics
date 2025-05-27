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
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.fftpack import dct, idct
from scipy.optimize import minimize_scalar
import sv_from_coe
import lunar_position

def lunar_perturbation_test():
    '''
    This function uses MATLAB's ode45 to integrate
    Equations 10.84, the Gauss variational equations, for a lunar
    gravitational perturbation.
    
    User py-functions required:  sv_from_coe, lunar_position
    User subfunctions required: solveit rates
    '''
    global JD

    # Conversion factors:
    hours = 3600  # Hours to seconds
    days = 24 * hours  # Days to seconds
    deg = np.pi / 180  # Degrees to radians

    # Constants:
    mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
    mu3 = 4904.8695  # Moon's gravitational parameter (km^3/s^2)
    RE = 6378.14  # Earth's radius (km)

    # Initial data for each of the three given orbits:
    orbit_types = ['GEO', 'HEO', 'LEO']

    # GEO
    n = 1
    a0 = 42164  # Semimajor axis (km)
    e0 = 0.0001  # Eccentricity
    w0 = 0  # Argument of perigee (rad)
    RA0 = 0  # Right ascension (rad)
    i0 = 1 * deg  # Inclination (rad)
    TA0 = 0  # True anomaly (rad)
    JD0 = 2454283  # Julian Day
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n)

    # HEO
    n = 2
    a0 = 26553.4
    e0 = 0.741
    w0 = 270 * deg
    RA0 = 0
    i0 = 63.4 * deg
    TA0 = 0
    JD0 = 2454283
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n)

    # LEO
    n = 3
    a0 = 6678.136
    e0 = 0.01
    w0 = 0
    RA0 = 0
    i0 = 28.5 * deg
    TA0 = 0
    JD0 = 2454283
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n)

def solveit(a0, e0, w0, RA0, i0, TA0, JD0, n):
    '''
    Calculations and plots common to all of the orbits
    '''
    # Constants:
    mu = 398600.4418
    mu3 = 4904.8695 
    days = 86400

    # Initial orbital parameters:
    h0 = np.sqrt(mu * a0 * (1 - e0**2))  # Angular momentum (km^2/s)
    T0 = 2 * np.pi / np.sqrt(mu) * a0**1.5  # Period (s)

    # Store initial orbital elements in coe0:
    coe0 = np.array([h0, e0, RA0, i0, w0, TA0])

    # Integration time interval:
    t0 = 0
    tf = 60 * days
    nout = 400
    tspan = np.linspace(t0, tf, nout)

    # Solve ODE:
    sol = solve_ivp(rates, [t0, tf], coe0, t_eval=tspan, args=(mu, mu3, JD0), method='RK45', atol=1e-8, rtol=1e-8)

    # Extract results:
    t = sol.t
    y = sol.y.T
    RA = y[:, 2]
    i = y[:, 3]
    w = y[:, 4]

    # Smooth the data to eliminate short-period variations:
    RA = rsmooth(RA)
    i = rsmooth(i)
    w = rsmooth(w)

    # Plot results:
    t_plot = t[2:-2] 
    RA_plot = RA[2:-2]
    i_plot = i[2:-2]
    w_plot = w[2:-2]

    plt.figure(n)
    plt.subplot(1, 3, 1)
    plt.plot(t_plot / days, (RA_plot - RA0) / (np.pi / 180))
    plt.title('Right Ascension vs Time')
    plt.xlabel('Time (days)')
    plt.ylabel('RA (deg)')
    plt.tight_layout()

    plt.subplot(1, 3, 2)
    plt.plot(t_plot / days, (i_plot - i0) / (np.pi / 180))
    plt.title('Inclination vs Time')
    plt.xlabel('Time (days)')
    plt.ylabel('Inclination (deg)')
    plt.tight_layout()

    plt.subplot(1, 3, 3)
    plt.plot(t_plot / days, (w_plot - w0) / (np.pi / 180))
    plt.title('Argument of Perigee vs Time')
    plt.xlabel('Time (days)')
    plt.ylabel('Argument of Perigee (deg)')
    plt.tight_layout()

    plt.show()

def rates(t, f, mu, mu3, JD0):
    days = 86400

    # Extract orbital elements:
    h, e, RA, i, w, TA = f
    phi = w + TA

    # State vector:
    coe = [h, e, RA, i, w, TA]
    R, V = sv_from_coe.sv_from_coe(coe, mu)

    # Unit vectors in the RSW system:
    r = np.linalg.norm(R)
    ur = R / r
    H = np.cross(R, V)
    uh = H / np.linalg.norm(H)
    us = np.cross(uh, ur)

    # Update Julian day:
    JD = JD0 + t / days

    # Moon's position:
    R_m = lunar_position.lunar_position(JD)
    r_m = np.linalg.norm(R_m)

    R_rel = R_m - R
    r_rel = np.linalg.norm(R_rel)

    # Perturbations:
    q = np.dot(R, (2 * R_m - R)) / r_m**2
    F = (q**2 - 3 * q + 3) * q / (1 + (1 - q)**1.5)

    ap = mu3 / r_rel**3 * (F * R_m - R)
    apr = np.dot(ap, ur)
    aps = np.dot(ap, us)
    aph = np.dot(ap, uh)

    # Gauss variational equations:
    hdot = r * aps
    edot = h / mu * np.sin(TA) * apr + 1 / mu / h * ((h**2 + mu * r) * np.cos(TA) + mu * e * r) * aps
    RAdot = r / h / np.sin(i) * np.sin(phi) * aph
    idot = r / h * np.cos(phi) * aph
    wdot = -h * np.cos(TA) / mu / e * apr + (h**2 + mu * r) / mu / e / h * np.sin(TA) * aps - r * np.sin(phi) / h / np.tan(i) * aph
    TAdot = h / r**2 + 1 / e / h * (h**2 / mu * np.cos(TA) * apr - (r + h**2 / mu) * np.sin(TA) * aps)

    return [hdot, edot, RAdot, idot, wdot, TAdot]

def rsmooth(y):
    '''
    Apply recursive smoothing to the data.
    '''
    y = np.asarray(y)
    n = len(y)
    Lambda = -2 + 2 * np.cos(np.arange(n) * np.pi / n)
    W = np.ones(n)
    z = y.copy()

    for _ in range(6):
        tol = float('inf')
        while tol > 1e-5:
            DCTy = dct(W * (y - z) + z)
            GCVscore = lambda p: np.sum((W * (y - idct(1 / (1 + 10**p * Lambda**2) * DCTy)))**2)
            s = 10**minimize_scalar(GCVscore, bounds=(-15, 38), method='bounded').x
            Gamma = 1 / (1 + s * Lambda**2)
            new_z = idct(Gamma * DCTy)
            tol = np.linalg.norm(new_z - z) / np.linalg.norm(z)
            z = new_z
    return z

if __name__ == "__main__":
    lunar_perturbation_test()