# ALGORITHM 3.3: SOLUTION OF THE UNIVERSAL KEPLER'S EQUATION
# USING NEWTON'S METHOD

import numpy as np
import stumpC
import stumpS

def kepler_U(dt, ro, vro, a, mu):
    '''
    This function uses Newton's method to solve the universal
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
    '''
    #...Check the supplied orbital parameters:
    if not np.isfinite([dt, ro, vro, a, mu]).all() or ro <= 0 or mu <= 0:
        raise ValueError("Universal propagation requires finite inputs and positive ro and mu.")
    if dt == 0:
        return 0.0

    error = 1.e-12
    nMax = 200

    #...Bracket the universal anomaly in the direction of propagation:
    direction = np.sign(dt)
    target = np.sqrt(mu) * abs(dt)

    #...Universal Kepler equation and its derivative:
    def equation(x):
        C     = stumpC.stumpC(a * x**2)
        S     = stumpS.stumpS(a * x**2)
        F     = ro * vro / np.sqrt(mu) * x**2 * C + (1 - a * ro) * x**3 * S + ro * x - np.sqrt(mu) * dt
        dFdx  = ro * vro / np.sqrt(mu) * x * (1 - a * x**2 * S) + (1 - a * ro) * x**2 * C + ro
        return F, dFdx
    
    lo = 0
    if a < 0:
        scale = np.sqrt(-1 / a)
        hi = max(1, scale * np.log1p(target / (ro * scale)))
    else:
        hi = max(1, target * max(a, 1 / ro))
    for n in range(100):
        F, _ = equation(direction * hi)
        if not np.isfinite(F):
            raise RuntimeError("Universal anomaly bracket overflowed.")
        if direction * F >= 0:
            break
        hi = 2 * hi
    else:
        raise RuntimeError("Could not bracket the universal anomaly.")

    #...Iterate on Equation 3.65 with a safeguarded Newton step
    #...until the residual or bracket width satisfies the error tolerance:
    u = (lo + hi) / 2
    for n in range(nMax):
        x = direction * u
        F, dFdx = equation(x)
        if not np.isfinite([F, dFdx]).all():
            raise RuntimeError("Universal anomaly iteration is nonfinite.")
        value = direction * F
        if abs(value) <= error * max(target, 1):
            return x
        if value > 0:
            hi = u
        else:
            lo = u
        if dFdx > 0:
            trial = u - value / dFdx
        else:
            trial = np.nan
        if not lo < trial < hi:
            trial = (lo + hi) / 2
        if hi - lo <= error * max(1, abs(u)):
            x = direction * trial
            return x
        u = trial

    raise RuntimeError("Universal Kepler equation failed to converge.")