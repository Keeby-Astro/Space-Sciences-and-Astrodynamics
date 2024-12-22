import numpy as np
import stumpC
import stumpS

def fDot_and_gDot(x, r, ro, a, mu):
    '''
    This function calculates the time derivatives of the
    Lagrange f and g coefficients.

    Parameters:
    mu   - the gravitational parameter (km^3/s^2)
    a    - reciprocal of the semimajor axis (1/km)
    ro   - the radial position at time to (km)
    t    - the time elapsed since initial state vector (s)
    r    - the radial position after time t (km)
    x    - the universal anomaly after time t (km^0.5)
    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)

    User py-functions required: stumpC, stumpS
    '''
    z = a * x**2

    # Equation 3.69c:
    fdot = np.sqrt(mu) / (r * ro) * (z * stumpS.stumpS(z) - 1) * x

    # Equation 3.69d:
    gdot = 1 - x**2 / r * stumpC.stumpC(z)

    return fdot, gdot