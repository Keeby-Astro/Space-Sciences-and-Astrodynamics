# Algorithm 3.3: Solution of the universal Kepler’s equation using Newton’s
# method.
import numpy as np
import stumpC
import stumpS

def kepler_U(dt, ro, vro, a, mu):
    """
    This function uses Newton’s method to solve the universal
    Kepler equation for the universal anomaly.

    mu   - gravitational parameter (km^3/s^2)
    x    - the universal anomaly (km^0.5)
    dt   - time since x = 0 (s)
    ro   - radial position (km) when x = 0
    vro  - radial velocity (km/s) when x = 0
    a    - reciprocal of the semimajor axis (1/km)
    z    - auxiliary variable (z = a*x^2)
    C    - value of Stumpff function C(z)
    S    - value of Stumpff function S(z)
    n    - number of iterations for convergence
    nMax - maximum allowable number of iterations

    User py-functions required: stumpC, stumpS
    """
    # Set an error tolerance and a limit on the number of iterations:
    error = 1.e-8
    nMax = 1000

    # Starting value for x:
    x = np.sqrt(mu) * abs(a) * dt

    # Iterate on the equation until convergence occurs within the error tolerance:
    n = 0
    ratio = 1
    while abs(ratio) > error and n <= nMax:
        n += 1
        C = stumpC(a * x**2)
        S = stumpS(a * x**2)
        F = ro * vro / np.sqrt(mu) * x**2 * C + (1 - a * ro) * x**3 * S + ro * x - np.sqrt(mu) * dt
        dFdx = ro * vro / np.sqrt(mu) * x * (1 - a * x**2 * S) + (1 - a * ro) * x**2 * C + ro
        ratio = F / dFdx
        x -= ratio

    # Deliver a value for x, but report that nMax was reached:
    if n > nMax:
        print(f'\n **No. iterations of Kepler equation = {n}')
        print(f'\n F/dFdx = {F/dFdx}\n')

    return x
