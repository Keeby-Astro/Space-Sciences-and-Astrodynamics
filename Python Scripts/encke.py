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
from sv_from_coe import sv_from_coe
from rv_from_r0v0 import rv_from_r0v0
from coe_from_sv import coe_from_sv
import matplotlib.pyplot as plt

def encke():
    '''
    This function use Encke's method together with MATLAB's ode45
    to integrate Equation 10.2 for a J2 gravitational perturbation 
    given by Equation 10.30.
    
    User py-functions required: sv_from_coe, coe_from_sv, rv_from_r0v0
    User subfunction required: rates 
    '''
    #...Conversion factors:
    hours = 3600                  # Hours to seconds
    days = 24 * hours             # Days to seconds
    deg = np.pi / 180             # Degrees to radians

    #...Constants:
    global mu
    mu = 398600.4418               # Gravitational parameter (km^3/s^2)
    RE = 6378.14                   # Earth's radius (km)
    J2 = 1082.63e-6

    #...Initial orbital parameters (given):
    zp0 = 300                     # Perigee altitude (km)
    za0 = 3062                    # Apogee altitude (km)
    RA0 = 45 * deg                # Right ascension of the node (radians)
    i0  = 28 * deg                # Inclination (radians)
    w0  = 30 * deg                # Argument of perigee (radians)
    TA0 = 40 * deg                # True anomaly (radians)

    #...Initial orbital parameters (inferred):
    rp0 = RE + zp0                          # Perigee radius (km)
    ra0 = RE + za0                          # Apogee radius (km)
    e0  = (ra0 - rp0) / (ra0 + rp0)         # Eccentricity
    a0  = (ra0 + rp0) / 2                   # Semimajor axis (km)
    h0  = np.sqrt(rp0 * mu * (1 + e0))      # Angular momentum (km^2/s)
    T0  = 2 * np.pi / np.sqrt(mu) * a0**1.5 # Period (s)

    t0 = 0
    tf = 2 * days                 # Initial and final time (s)

    # Store the initial orbital elements in the array coe0:
    coe0 = [h0, e0, RA0, i0, w0, TA0]

    #...Obtain the initial state vector from sv_from_coe:
    R0, V0 = sv_from_coe(coe0, mu)  # R0 is the initial position vector
                                    # V0 is the initial velocity vector
    r0 = np.linalg.norm(R0)         # Magnitude of R0
    v0 = np.linalg.norm(V0)         # Magnitude of V0

    del_t = T0 / 100                # Time step for Encke procedure

    #...Begin the Encke integration:
    t = t0                         # Initialize the time scalar
    tsave = [t0]                   # Initialize the vector of solution times
    y = [np.concatenate((R0, V0))] # Initialize the state vector
    del_y0 = np.zeros(6)           # Initialize the state vector perturbation

    t += del_t                     # First time step

    # Loop over the time interval [t0, tf] with equal increments del_t:
    while t <= tf + del_t / 2:

        def rates(t, f):
            '''
            This function calculates the time rates of Encke's deviation in position
            del_r and velocity del_v.
            '''
            del_r = f[:3]          # Position deviation
            del_v = f[3:]          # Velocity deviation

            # Compute the state vector on the osculating orbit at time t
            Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, mu)

            # Calculate the components of the state vector on the perturbed orbit:
            Rpp = Rosc + del_r
            Vpp = Vosc + del_v
            rosc = np.linalg.norm(Rosc)
            rpp = np.linalg.norm(Rpp)

            # Compute the J2 perturbing acceleration:
            xx, yy, zz = Rpp

            fac = 3 / 2 * J2 * (mu / rpp**2) * (RE / rpp)**2
            ap = -fac * np.array([
                (1 - 5 * (zz / rpp)**2) * (xx / rpp),
                (1 - 5 * (zz / rpp)**2) * (yy / rpp),
                (3 - 5 * (zz / rpp)**2) * (zz / rpp)
            ])

            # Compute the total perturbing acceleration:
            F = 1 - (rosc / rpp)**3
            del_a = -mu / rosc**3 * (del_r - F * Rpp) + ap

            return np.concatenate((del_v, del_a))

        sol = solve_ivp(rates, [t0, t], del_y0, max_step=del_t)
        z = sol.y[:, -1]

        # Compute the osculating state vector at time t:
        Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0, mu)

        # Rectify:
        R0 = Rosc + z[:3]
        V0 = Vosc + z[3:]
        t0 = t

        # Prepare for next time step:
        tsave.append(t)
        y.append(np.concatenate((R0, V0)))
        t += del_t
        del_y0 = np.zeros(6)

    # Extract the orbital elements from the state vector:
    t = np.array(tsave)
    y = np.array(y)

    n_times = len(t)
    r, v, h, e, RA, i, w, TA = ([] for _ in range(8))

    #...At each solution time extract the orbital elements from the state
    #   vector using Algorithm 4.2:
    for j in range(n_times):
        R = y[j, :3]
        V = y[j, 3:]
        r.append(np.linalg.norm(R))
        v.append(np.linalg.norm(V))
        coe = coe_from_sv(R, V, mu)
        h.append(coe[0])
        e.append(coe[1])
        RA.append(coe[2])
        i.append(coe[3])
        w.append(coe[4])
        TA.append(coe[5])

    # Convert lists to numpy arrays for easier handling:
    r, v, h, e, RA, i, w, TA = map(np.array, [r, v, h, e, RA, i, w, TA])

    #...Plot selected osculating elements:
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t / 3600, (RA - RA0) / deg)
    plt.title('Variation of Right Ascension')
    plt.xlabel('hours')
    plt.ylabel('\u0394\u03a9 (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(t / 3600, (w - w0) / deg)
    plt.title('Variation of Argument of Perigee')
    plt.xlabel('hours')
    plt.ylabel('\u0394\u03c9 (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.figure(2)
    plt.subplot(3, 1, 1)
    plt.plot(t / 3600, h - h0)
    plt.title('Variation of Angular Momentum')
    plt.xlabel('hours')
    plt.ylabel('\u0394h (km^2/s)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(3, 1, 2)
    plt.plot(t / 3600, e - e0)
    plt.title('Variation of Eccentricity')
    plt.xlabel('hours')
    plt.ylabel('\u0394e')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.subplot(3, 1, 3)
    plt.plot(t / 3600, (i - i0) / deg)
    plt.title('Variation of Inclination')
    plt.xlabel('hours')
    plt.ylabel('\u0394i (deg)')
    plt.grid(True, which='both')
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    encke()
