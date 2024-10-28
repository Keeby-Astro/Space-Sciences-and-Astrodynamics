# Algorithm 3.1: Solution of Kepler’s equation by Newton’s method.
import numpy as np

def kepler_E(e, M):
    """
    This function uses Newton’s method to solve Kepler’s
    equation E - e*sin(E) = M for the eccentric anomaly,
    given the eccentricity and the mean anomaly.

    E  - eccentric anomaly (radians)
    e  - eccentricity, passed from the calling program
    M  - mean anomaly (radians), passed from the calling program
    pi - 3.1415926...

    User py-functions required: none
    """
    # Set an error tolerance:
    error = 1.e-8

    # Select a starting value for E:
    if M < np.pi:
        E = M + e / 2
    elif M < np.pi:
        E = M - e / 2
    else:
        E = M

    # Iterate on Equation 3.17 until E is determined to within
    # the error tolerance:
    ratio = 1
    while abs(ratio) > error:
        ratio = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E = E - ratio
        E = E % (2*np.pi)

    return E
