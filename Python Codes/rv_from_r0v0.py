# Algorithm 3.4: Calculation of the state vector given the initial state vector
# and the time lapse Δt.
import numpy as np
import kepler_U
import f_and_g
import fDot_and_gDot

def rv_from_r0v0(R0, V0, t, mu):
    """
    This function computes the state vector (R,V) from the
    initial state vector (R0,V0) and the elapsed time.

    mu - gravitational parameter (km^3/s^2)
    R0 - initial position vector (km)
    V0 - initial velocity vector (km/s)
    t  - elapsed time (s)
    R  - final position vector (km)
    V  - final velocity vector (km/s)

    % User py-functions required: kepler_U, f_and_g, fDot_and_gDot
    """
    
    # Magnitudes of R0 and V0
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)

    # Initial radial velocity
    vr0 = np.dot(R0, V0) / r0

    # Reciprocal of the semimajor axis (from the energy equation)
    alpha = 2 / r0 - v0**2 / mu

    # Compute the universal anomaly
    x = kepler_U(t, r0, vr0, alpha)

    # Compute the f and g functions
    f, g = f_and_g(x, t, r0, alpha)

    # Compute the final position vector
    R = f * R0 + g * V0

    # Compute the magnitude of R
    r = np.linalg.norm(R)

    # Compute the derivatives of f and g
    fdot, gdot = fDot_and_gDot(x, r, r0, alpha)

    # Compute the final velocity
    V = fdot * R0 + gdot * V0

    return R, V
