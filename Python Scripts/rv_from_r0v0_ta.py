# ALGORITHM 2.3: CALCULATE THE STATE VECTOR FROM THE INITIAL
# STATE VECTOR AND THE CHANGE IN TRUE ANOMALY

import numpy as np
import f_and_g_ta
import fDot_and_gDot_ta

def rv_from_r0v0_ta(r0, v0, dt, mu):
    '''
    This function computes the state vector (r,v) from the
    initial state vector (r0,v0) and the change in true anomaly.

    mu - gravitational parameter (km^3/s^2)
    r0 - initial position vector (km)
    v0 - initial velocity vector (km/s)
    dt - change in true anomaly (degrees)
    r  - final position vector (km)
    v  - final velocity vector (km/s)

    User py-functions required: f_and_g_ta, fDot_and_gDot_ta
    '''
    #...Compute the f and g functions and their derivatives
    f, g       = f_and_g_ta.f_and_g_ta(r0, v0, dt, mu)
    fdot, gdot = fDot_and_gDot_ta.fDot_and_gDot_ta(r0, v0, dt, mu)

    #...Compute the final position and velocity vectors
    r = f * np.array(r0) + g * np.array(v0)
    v = fdot * np.array(r0) + gdot * np.array(v0)

    return r, v
