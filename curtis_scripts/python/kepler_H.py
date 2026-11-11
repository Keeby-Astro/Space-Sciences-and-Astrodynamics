# ALGORITHM 3.2: SOLUTION OF KEPLER'S EQUATION FOR THE
# HYPERBOLA USING NEWTON'S METHOD

import numpy as np

def kepler_H(e, M):
    '''
    This function uses Newton's method to solve Kepler's equation
    for the hyperbola e*sinh(F) - F = M for the hyperbolic
    eccentric anomaly, given the eccentricity and the hyperbolic
    mean anomaly.

    F - hyperbolic eccentric anomaly (radians)
    e - eccentricity, passed from the calling program
    M - hyperbolic mean anomaly (radians), passed from the
    calling program

    User py-functions required: none
    '''
    # Check the supplied orbital parameters:
    if not np.isfinite([e, M]).all() or e <= 1:
        raise ValueError("Hyperbolic Kepler equation requires e > 1 and finite M.")
    if M == 0:
        return 0.0

    error = 1.e-13
    nMax = 100

    #...Starting value for F and an upper bound on the root:
    m = abs(M)
    lo = 0
    hi = np.arcsinh(m / e) + 1
    while e * np.sinh(hi) - hi < m:
        hi = 2 * hi
    F = np.arcsinh(m / e)

    #...Iterate on Equation 3.45, using bisection when a Newton step leaves
    #...the bracket. Limit the iteration count and test the error tolerance:
    for n in range(nMax):
        residual = e * np.sinh(F) - F - m
        if residual == 0:
            break
        if residual > 0:
            hi = F
        else:
            lo = F
        trial = F - residual / (e * np.cosh(F) - 1)
        if trial <= lo or trial >= hi:
            trial = (lo + hi) / 2
        ratio = trial - F
        F = trial
        if abs(ratio) <= error * max(1, abs(F)):
            break
    else:
        raise RuntimeError("Hyperbolic Kepler equation failed to converge.")

    F = np.sign(M) * F
    return F