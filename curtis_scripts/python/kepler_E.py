# ALGORITHM 3.1: SOLUTION OF KEPLER'S EQUATION BY NEWTON'S METHOD

import numpy as np

def kepler_E(e, M):
    '''
    This function uses Newton's method to solve Kepler's
    equation E - e*sin(E) = M for the eccentric anomaly,
    given the eccentricity and the mean anomaly.

    E  - eccentric anomaly (radians)
    e  - eccentricity, passed from the calling program
    M  - mean anomaly (radians), passed from the calling program
    pi - 3.1415926...

    User py-functions required: none
    '''
    #...Check the supplied orbital parameters:
    if not np.isfinite([e, M]).all() or e < 0 or e >= 1:
        raise ValueError("Elliptic Kepler equation requires 0 <= e < 1 and finite M.")

    #...Set an error tolerance and a limit on the number of iterations:
    error = 1.e-13
    nMax = 100

    #...Select a starting value for E and bracket the root:
    m = (M + np.pi) % (2 * np.pi) - np.pi
    lo = m - e
    hi = m + e
    E = m + np.sign(m) * 0.85 * e

    #...Iterate on Equation 3.17 until E is determined to within
    #...the error tolerance. Bisect when a Newton step leaves the bracket:
    for n in range(nMax):
        F = E - e * np.sin(E) - m
        if F == 0:
            break
        if F > 0:
            hi = E
        else:
            lo = E
        trial = E - F / (1 - e * np.cos(E))
        if trial <= lo or trial >= hi:
            trial = (lo + hi) / 2
        ratio = trial - E
        E = trial
        if abs(ratio) <= error * max(1, abs(E)):
            break
    else:
        raise RuntimeError("Elliptic Kepler equation failed to converge.")

    E = E + (M - m)
    return E