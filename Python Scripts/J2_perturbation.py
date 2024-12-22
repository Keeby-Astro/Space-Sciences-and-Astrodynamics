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
import matplotlib.pyplot as plt

def J2_perturbation():
    '''
    This function use MATLAB's ode45 to numerically integrate Equations 10.89
    (the Gauss planetary equations) to determine the J2 perturbation of 
    the orbital elements.
    
    User py-functions required:  None
    User subfunctions required: rates    
    '''
    # Conversion factors:
    hours = 3600  # Hours to seconds
    days = 24 * hours  # Days to seconds
    deg = np.pi / 180  # Degrees to radians

    # Constants:
    mu = 398600.4418  # Gravitational parameter (km^3/s^2)
    RE = 6378.14  # Earth's radius (km)
    J2 = 1082.63e-6  # Earth's J2

    # Initial orbital parameters (given):
    rp0 = RE + 300  # Perigee radius (km)
    ra0 = RE + 3062  # Apogee radius (km)
    RA0 = 45 * deg  # Right ascension of the node (radians)
    i0 = 28 * deg  # Inclination (radians)
    w0 = 30 * deg  # Argument of perigee (radians)
    TA0 = 40 * deg  # True anomaly (radians)

    # Initial orbital parameters (inferred):
    e0 = (ra0 - rp0) / (ra0 + rp0)  # Eccentricity
    h0 = np.sqrt(rp0 * mu * (1 + e0))  # Angular momentum (km^2/s)
    a0 = (rp0 + ra0) / 2  # Semi-major axis (km)
    T0 = 2 * np.pi / np.sqrt(mu) * a0 ** 1.5  # Period (s)

    # Store initial orbital elements in the vector coe0:
    coe0 = [h0, e0, RA0, i0, w0, TA0]

    # Use solve_ivp to integrate the Gauss variational equations:
    t0 = 0
    tf = 2 * days
    nout = 5000  # Number of solution points to output for plotting purposes
    tspan = np.linspace(t0, tf, nout)

    def rates(t, f):
        '''
        This function calculates the time rates of the orbital elements
        from Gauss's variational equations (Equations 12.89).
        '''

        # The orbital elements at time t:
        h, e, RA, i, w, TA = f

        r = h ** 2 / mu / (1 + e * np.cos(TA))  # The radius
        u = w + TA  # Argument of latitude

        # Orbital element rates:
        hdot = -3 / 2 * J2 * mu * RE ** 2 / r ** 3 * np.sin(i) ** 2 * np.sin(2 * u)

        edot = 3 / 2 * J2 * mu * RE ** 2 / h / r ** 3 * (
            h ** 2 / mu / r * np.sin(TA) * (3 * np.sin(i) ** 2 * np.sin(u) ** 2 - 1)
            - np.sin(2 * u) * np.sin(i) ** 2 * ((2 + e * np.cos(TA)) * np.cos(TA) + e)
        )

        TAdot = h / r ** 2 + 3 / 2 * J2 * mu * RE ** 2 / e / h / r ** 3 * (
            h ** 2 / mu / r * np.cos(TA) * (3 * np.sin(i) ** 2 * np.sin(u) ** 2 - 1)
            + np.sin(2 * u) * np.sin(i) ** 2 * np.sin(TA) * (h ** 2 / mu / r + 1)
        )

        RAdot = -3 * J2 * mu * RE ** 2 / h / r ** 3 * np.sin(u) ** 2 * np.cos(i)

        idot = -3 / 4 * J2 * mu * RE ** 2 / h / r ** 3 * np.sin(2 * u) * np.sin(2 * i)

        wdot = 3 / 2 * J2 * mu * RE ** 2 / e / h / r ** 3 * (
            -h ** 2 / mu / r * np.cos(TA) * (3 * np.sin(i) ** 2 * np.sin(u) ** 2 - 1)
            - np.sin(2 * u) * np.sin(i) ** 2 * np.sin(TA) * (2 + e * np.cos(TA))
            + 2 * e * np.cos(i) ** 2 * np.sin(u) ** 2
        )

        return [hdot, edot, RAdot, idot, wdot, TAdot]

    sol = solve_ivp(rates, [t0, tf], coe0, t_eval=tspan, rtol=1e-8, atol=1e-8)

    # Assign the time histories mnemonic variable names:
    t = sol.t
    h, e, RA, i, w, TA = sol.y

    # Plot the time histories of the osculating elements:
    plt.figure(figsize=(10, 15))

    plt.subplot(5, 1, 1)
    plt.plot(t / hours, (RA - RA0) / deg)
    plt.title('Right Ascension (degrees)')
    plt.xlabel('Hours')
    plt.grid()

    plt.subplot(5, 1, 2)
    plt.plot(t / hours, (w - w0) / deg)
    plt.title('Argument of Perigee (degrees)')
    plt.xlabel('Hours')
    plt.grid()

    plt.subplot(5, 1, 3)
    plt.plot(t / hours, h - h0)
    plt.title('Angular Momentum (km^2/s)')
    plt.xlabel('Hours')
    plt.grid()

    plt.subplot(5, 1, 4)
    plt.plot(t / hours, e - e0)
    plt.title('Eccentricity')
    plt.xlabel('Hours')
    plt.grid()

    plt.subplot(5, 1, 5)
    plt.plot(t / hours, (i - i0) / deg)
    plt.title('Inclination (degrees)')
    plt.xlabel('Hours')
    plt.grid()

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    J2_perturbation()