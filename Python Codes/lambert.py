# Algorithm 5.2: Solution of Lambert’s problem.

import numpy as np
import stumpC
import stumpS

def lambert(R1, R2, t, string='pro'):
    """
    This function solves Lambert's problem.

    mu - gravitational parameter (km^3/s^2)
    R1, R2 - initial and final position vectors (km)
    r1, r2 - magnitudes of R1 and R2
    t - the time of flight from R1 to R2 (a constant) (s)
    V1, V2 - initial and final velocity vectors (km/s)
    c12 - cross product of R1 into R2
    theta - angle between R1 and R2
    string - 'pro' if the orbit is prograde
    'retro' if the orbit is retrograde
    A - a constant given by Equation 5.35
    z - alpha*x^2, where alpha is the reciprocal of the
    semimajor axis and x is the universal anomaly
    y(z) - a function of z given by Equation 5.38
    F(z,t) - a function of the variable z and constant t,
    - given by Equation 5.40
    dFdz(z) - the derivative of F(z,t), given by Equation 5.43
    ratio - F/dFdz
    tol - tolerance on precision of convergence
    nmax - maximum number of iterations of Newton’s procedure
    f, g - Lagrange coefficients
    gdot - time derivative of g
    C(z), S(z) - Stumpff functions
    dum - a dummy variable

    User py-functions required: stumpC and stumpS
    """
    global mu

    # Magnitudes of R1 and R2
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)
    c12 = np.cross(R1, R2)
    theta = np.arccos(np.dot(R1, R2) / (r1 * r2))

    # Determine whether the orbit is prograde or retrograde
    if string not in ['retro', 'pro']:
        string = 'pro'
        print('\n ** Prograde trajectory assumed.\n')

    if string == 'pro':
        if c12[2] <= 0:
            theta = 2 * np.pi - theta
    elif string == 'retro':
        if c12[2] >= 0:
            theta = 2 * np.pi - theta

    # Equation 5.35
    A = np.sin(theta) * np.sqrt(r1 * r2 / (1 - np.cos(theta)))

    # Determine approximately where F(z,t) changes sign, and
    # use that value of z as the starting value for Equation 5.45
    z = -100
    while F(z, t, r1, r2, A) < 0:
        z += 0.1

    # Set an error tolerance and a limit on the number of iterations
    tol = 1.e-8
    nmax = 5000

    # Iterate on Equation 5.45 until z is determined to within the error tolerance
    ratio = 1
    n = 0
    while abs(ratio) > tol and n <= nmax:
        n += 1
        ratio = F(z, t, r1, r2, A) / dFdz(z, r1, r2, A)
        z -= ratio

    # Report if the maximum number of iterations is exceeded
    if n >= nmax:
        print('\n\n **Number of iterations exceeds %g \n\n ' % nmax)

    # Equation 5.46a
    f = 1 - y(z, r1, r2, A) / r1

    # Equation 5.46b
    g = A * np.sqrt(y(z, r1, r2, A) / mu)

    # Equation 5.46d
    gdot = 1 - y(z, r1, r2, A) / r2

    # Equation 5.28
    V1 = (R2 - f * R1) / g

    # Equation 5.29
    V2 = (gdot * R2 - R1) / g

    return V1, V2

# Subfunctions used in the main body:
def y(z, r1, r2, A):
    return r1 + r2 + A * (z * S(z) - 1) / np.sqrt(C(z))

def F(z, t, r1, r2, A):
    return (y(z, r1, r2, A) / C(z))**1.5 * S(z) + A * np.sqrt(y(z, r1, r2, A)) - np.sqrt(mu) * t

def dFdz(z, r1, r2, A):
    if z == 0:
        return np.sqrt(2) / 40 * y(0, r1, r2, A)**1.5 + A / 8 * (np.sqrt(y(0, r1, r2, A)) + A * np.sqrt(1 / 2 / y(0, r1, r2, A)))
    else:
        return (y(z, r1, r2, A) / C(z))**1.5 * (1 / 2 / z * (C(z) - 3 * S(z) / 2 / C(z)) + 3 * S(z)**2 / 4 / C(z)) + A / 8 * (3 * S(z) / C(z) * np.sqrt(y(z, r1, r2, A)) + A * np.sqrt(C(z) / y(z, r1, r2, A)))

def C(z):
    return stumpC(z)

def S(z):
    return stumpS(z)
