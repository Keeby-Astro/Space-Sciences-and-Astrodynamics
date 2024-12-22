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
from sv_from_coe import sv_from_coe
from solar_position import solar_position

def solar_perturbation():
    '''
    This function uses Python's solve_ivp to integrate
    the Gauss variational equations for a solar gravitational perturbation.
    
    User modules required: sv_from_coe, solar_position
    '''
    global JD  # Julian day

    # Conversion factors:
    hours = 3600  # Hours to seconds
    days = 24 * hours  # Days to seconds
    deg = np.pi / 180  # Degrees to radians

    # Constants
    mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
    mu3 = 132.712e9  # Sun's gravitational parameter (km^3/s^2)
    RE = 6378.14  # Earth's radius (km)

    # Initial data for each of the three given orbits:
    orbits = [
        {"name": "GEO", "a0": 42164, "e0": 0.0001, "i0": 1 * deg, "w0": 0, "RA0": 0, "TA0": 0},
        {"name": "HEO", "a0": 26553.4, "e0": 0.741, "i0": 63.4 * deg, "w0": 270, "RA0": 0, "TA0": 0},
        {"name": "LEO", "a0": 6678.136, "e0": 0.01, "i0": 28.5 * deg, "w0": 0, "RA0": 0, "TA0": 0}
    ]

    # Julian Day for all orbits
    JD0 = 2454283

    # Solve for each orbit
    for n, orbit in enumerate(orbits):
        solveit(orbit, JD0, mu, mu3, days, deg)

def solveit(orbit, JD0, mu, mu3, days, deg):
    '''
    Solve and plot the orbital perturbations for the given orbit.
    '''
    a0, e0, i0, w0, RA0, TA0 = orbit["a0"], orbit["e0"], orbit["i0"], orbit["w0"], orbit["RA0"], orbit["TA0"]

    # Initial orbital parameters:
    h0 = np.sqrt(mu * a0 * (1 - e0 ** 2))  # Angular momentum (km^2/s)
    coe0 = [h0, e0, RA0, i0, w0, TA0]  # Initial orbital elements

    # Time span for integration
    t0 = 0
    tf = 720 * days
    nout = 400
    tspan = np.linspace(t0, tf, nout)

    # Solve the ODE
    sol = solve_ivp(
        rates,
        [t0, tf],
        coe0,
        t_eval=tspan,
        args=(mu, mu3, JD0, days),
        rtol=1e-8,
        atol=1e-8
    )

    # Extract solution components
    t = sol.t
    y = sol.y
    RA, i, w = y[2], y[3], y[4]

    # Smooth data (placeholder, implement as needed)
    RA = rsmooth(RA)
    i = rsmooth(i)
    w = rsmooth(w)

    # Plot results
    plt.figure()

    plt.subplot(1, 3, 1)
    plt.plot(t / days, (RA - RA0) / deg)
    plt.title("Right Ascension vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("RA (deg)")
    plt.tight_layout()

    plt.subplot(1, 3, 2)
    plt.plot(t / days, (i - i0) / deg)
    plt.title("Inclination vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("Inclination (deg)")
    plt.tight_layout()

    plt.subplot(1, 3, 3)
    plt.plot(t / days, (w - w0) / deg)
    plt.title("Argument of Perigee vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("Argument of Perigee (deg)")
    plt.tight_layout()

    plt.show()

def rates(t, f, mu, mu3, JD0, days):
    '''
    Compute the rates of change of the orbital elements.
    '''
    h, e, RA, i, w, TA = f
    phi = w + TA  # Argument of latitude

    # State vector from orbital elements
    coe = [h, e, RA, i, w, TA]
    R, V = sv_from_coe(coe, mu)

    # Update Julian Day
    JD = JD0 + t / days

    # Position vector of the sun
    lamda, eps, R_S = solar_position(JD)
    r_S = np.linalg.norm(R_S)

    # Perturbation components
    R_rel = R_S - R
    r_rel = np.linalg.norm(R_rel)

    q = np.dot(R, (2 * R_S - R)) / r_S ** 2
    F = (q ** 2 - 3 * q + 3) * q / (1 + (1 - q) ** 1.5)
    ap = mu3 / r_rel ** 3 * (F * R_S - R)

    # Variational equations
    apr = np.dot(ap, R / np.linalg.norm(R))
    aps = np.dot(ap, np.cross(np.cross(R, V), R) / np.linalg.norm(np.cross(np.cross(R, V), R)))
    aph = np.dot(ap, np.cross(R, V) / np.linalg.norm(np.cross(R, V)))

    hdot = np.linalg.norm(R) * aps
    edot = (h / mu * np.sin(TA) * apr + (1 / mu / h) * ((h ** 2 + mu * np.linalg.norm(R)) * np.cos(TA) + mu * e * np.linalg.norm(R)) * aps)
    RAdot = (np.linalg.norm(R) / h / np.sin(i) * np.sin(phi) * aph)
    idot = (np.linalg.norm(R) / h * np.cos(phi) * aph)
    wdot = (-h * np.cos(TA) / mu / e * apr + (h ** 2 + mu * np.linalg.norm(R)) / mu / e / h * np.sin(TA) * aps - np.linalg.norm(R) * np.sin(phi) / h / np.tan(i) * aph)
    TAdot = (h / np.linalg.norm(R) ** 2 + 1 / e / h * (h ** 2 / mu * np.cos(TA) * apr - (np.linalg.norm(R) + h ** 2 / mu) * np.sin(TA) * aps))

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

if __name__ == '__main__':
    solar_perturbation()