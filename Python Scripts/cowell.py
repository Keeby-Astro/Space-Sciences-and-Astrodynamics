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
import atmosphere
import sv_from_coe

def cowell():
    '''
    This function uses MATLAB's ode45 to numerically
    integrate Equation 10.2 for atmospheric drag.
    
    User py-functions required: sv_from_coe, atmosphere
    User subfunctions required: rates, terminate                           
    '''
    # Conversion factors:
    hours = 3600                            # Hours to seconds
    days = 24 * hours                       # Days to seconds
    deg = np.pi / 180                       # Degrees to radians

    # Constants:
    mu = 398600.4418                        # Gravitational parameter (km^3/s^2)
    RE = 6378.14                            # Earth's radius (km)
    wE = np.array([0, 0, 7.2921159e-5])     # Earth's angular velocity (rad/s)

    # Satellite data:
    CD = 2.2                                # Drag coefficient
    m = 100                                 # Mass (kg)
    A = np.pi / 4 * (1 ** 2)                # Frontal area (m^2)

    # Initial orbital parameters (given):
    rp = RE + 215                           # Perigee radius (km)
    ra = RE + 939                           # Apogee radius (km)
    RA = 339.94 * deg                       # Right ascension of the node (radians)
    i = 65.1 * deg                          # Inclination (radians)
    w = 58 * deg                            # Argument of perigee (radians)
    TA = 332 * deg                          # True anomaly (radians)

    # Initial orbital parameters (inferred):
    e = (ra - rp) / (ra + rp)               # Eccentricity
    a = (rp + ra) / 2                       # Semimajor axis (km)
    h = np.sqrt(mu * a * (1 - e ** 2))      # Angular momentum (km^2/s)
    T = 2 * np.pi / np.sqrt(mu) * a ** 1.5  # Period (s)

    # Store initial orbital elements in the vector coe0:
    coe0 = [h, e, RA, i, w, TA]

    # Obtain the initial state vector from sv_from_coe:
    R0, V0 = sv_from_coe.sv_from_coe(coe0, mu)  # R0: Initial position vector, V0: Initial velocity vector
    r0 = np.linalg.norm(R0)                     # Magnitude of R0
    v0 = np.linalg.norm(V0)                     # Magnitude of V0

    # Use solve_ivp to integrate the equations of motion d/dt(R,V) = f(R,V)
    # from t0 to tf:
    t0 = 0
    tf = 120 * days
    y0 = np.hstack((R0, V0))            # Initial state vector
    nout = 40000                        # Number of solution points to output
    tspan = np.linspace(t0, tf, nout)   # Integration time interval

    # Termination event function
    def terminate(t, y):
        '''
        This function specifies the event at which ode45 terminates.
        '''
        R = y[:3]
        r = np.linalg.norm(R)
        alt = r - RE
        return alt - 100

    terminate.terminal = True
    terminate.direction = -1

    def rates(t, f):
        '''
        This function calculates the spacecraft acceleration from its
        position and velocity at time t.
        '''
        R = f[:3]                                            # Position vector (km/s)
        r = np.linalg.norm(R)                                # Distance from Earth's center (km)
        alt = r - RE                                         # Altitude (km)
        rho = atmosphere.atmosphere(alt)                     # Air density from US Standard Model (kg/m^3)
        V = f[3:]                                            # Velocity vector (km/s)
        Vrel = V - np.cross(wE, R)                           # Velocity relative to the atmosphere (km/s)
        vrel = np.linalg.norm(Vrel)                          # Speed relative to the atmosphere (km/s)
        uv = Vrel / vrel                                     # Relative velocity unit vector
        ap = -CD * A / m * rho * (1000 * vrel) ** 2 / 2 * uv # Drag acceleration (m/s^2)
        a0 = -mu * R / r ** 3                                # Gravitational acceleration (km/s^2)
        a = a0 + ap / 1000                                   # Total acceleration (km/s^2)
        return np.hstack((V, a))                             # Velocity and acceleration returned to solve_ivp

    solution = solve_ivp(rates, [t0, tf], y0, method='RK45', t_eval=tspan,
                         rtol=1e-8, atol=1e-8, events=terminate)

    t = solution.t
    y = solution.y.T

    # Extract the locally extreme altitudes:
    altitude = np.linalg.norm(y[:, :3], axis=1) - RE  # Altitude at each time
    max_altitude_indices = np.where((altitude[1:-1] > altitude[:-2]) & (altitude[1:-1] > altitude[2:]))[0] + 1
    min_altitude_indices = np.where((altitude[1:-1] < altitude[:-2]) & (altitude[1:-1] < altitude[2:]))[0] + 1

    apogee = np.column_stack((t[max_altitude_indices], altitude[max_altitude_indices]))
    perigee = np.column_stack((t[min_altitude_indices], altitude[min_altitude_indices]))

    # Plot perigee and apogee history:
    plt.figure()
    plt.plot(apogee[:, 0] / days, apogee[:, 1], 'b', linewidth=2, label='Apogee')
    plt.plot(perigee[:, 0] / days, perigee[:, 1], 'r', linewidth=2, label='Perigee')
    plt.grid(True, which='both')
    plt.xlabel('Time (days)')
    plt.ylabel('Altitude (km)')
    plt.ylim([0, 1000])
    plt.legend()
    plt.show()

if __name__ == "__main__":
    cowell()
