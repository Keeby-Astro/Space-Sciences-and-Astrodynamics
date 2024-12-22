import stumpC
import stumpS

def f_and_g(x, t, ro, a, mu):
    '''
    This function calculates the Lagrange f and g coefficients.

    mu - the gravitational parameter (km^3/s^2)
    a - reciprocal of the semimajor axis (1/km)
    ro - the radial position at time to (km)
    t - the time elapsed since ro (s)
    x - the universal anomaly after time t (km^0.5)
    f - the Lagrange f coefficient (dimensionless)
    g - the Lagrange g coefficient (s)

    User py-functions required: stumpC, stumpS
    '''
    z = a * x**2

    # Equation 3.69a:
    f = 1 - (x**2 / ro) * stumpC.stumpC(z)

    # Equation 3.69b:
    g = t - (1 / mu**0.5) * x**3 * stumpS.stumpS(z)

    return f, g