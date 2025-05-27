import numpy as np
import stumpC
import stumpS

def lambert(R1, R2, t, orbit_type, mu):
    '''
    This function solves Lambert's problem.

    mu      - gravitational parameter (km^3/s^2)
    R1, R2  - initial and final position vectors (km)
    r1, r2  - magnitudes of R1 and R2
    t       - the time of flight from R1 to R2 (a constant) (s)
    V1, V2  - initial and final velocity vectors (km/s)
    c12     - cross product of R1 into R2
    theta   - angle between R1 and R2
    string  - 'pro' if the orbit is prograde
              'retro' if the orbit is retrograde
    A       - a constant given by Equation 5.35
    z       - alpha*x^2, where alpha is the reciprocal of the
              semimajor axis and x is the universal anomaly
    y(z)    - a function of z given by Equation 5.38
    F(z,t)  - a function of the variable z and constant t,
            - given by Equation 5.40
    dFdz(z) - the derivative of F(z,t), given by Equation 5.43
    ratio   - F/dFdz
    tol     - tolerance on precision of convergence
    nmax    - maximum number of iterations of Newton's procedure
    f, g    - Lagrange coefficients
    gdot    - time derivative of g
              C(z), S(z) - Stumpff functions
    dum     - a dummy variable

    User py-functions required: stumpC and stumpS
    '''
    #...Magnitudes of R1 and R2:
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    c12 = np.cross(R1, R2)
    theta = np.arccos(np.dot(R1, R2) / (r1 * r2))

    #...Determine whether the orbit is prograde or retrograde:
    if orbit_type not in ["pro", "retro"]:
        orbit_type = "pro"
        print("** Prograde trajectory assumed. **")

    if orbit_type == "pro":
        if c12[2] <= 0.0:
            theta = 2.0 * np.pi - theta
    elif orbit_type == "retro":
        if c12[2] >= 0.0:
            theta = 2.0 * np.pi - theta

    # Equation 5.35
    A = np.sin(theta) * np.sqrt(r1 * r2 / (1.0 - np.cos(theta)))

    #...Stumpff functions:
    def C(z):
        '''Stumpff C function.'''
        return stumpC.stumpC(z)

    def S(z):
        '''Stumpff S function.'''
        return stumpS.stumpS(z)

    # Equation 5.38
    def y(z):
        if C(z) <= 0:
            return np.inf
        return r1 + r2 + A * ((z * S(z)) - 1.0) / np.sqrt(C(z))

    # Equation 5.40:
    def F(z, t_):
        if C(z) <= 0 or y(z) <= 0:
            return np.inf  # Return a large value to signal an invalid state
        return ((y(z) / C(z)) ** 1.5) * S(z) + A * np.sqrt(y(z)) - np.sqrt(mu) * t_

    # Equation 5.43:
    def dFdz(z):
        if abs(z) < 1.e-12:
            z = 1.e-12 
        term1 = ((y(z) / C(z)) ** 1.5) * (
            (1.0 / (2.0 * z)) * (C(z) - (3.0 * S(z) / (2.0 * C(z)))) 
            + 3.0 * (S(z) ** 2) / (4.0 * C(z))
        )
        term2 = (A / 8.0) * (
            3.0 * (S(z) / C(z)) * np.sqrt(y(z))
            + A * np.sqrt(C(z) / y(z))
        )

        return term1 + term2

    #...Determine approximately where F(z,t) changes sign, and
    #...use that value of z as the starting value for Equation 5.45:
    z = -100.0
    while F(z, t) < 0.0:
        z += 0.1
        # If z gets too big, break to avoid infinite loop
        if z > 1.e6:
            print("** F(z,t) did not become positive. Possibly no single-rev solution. **")
            break

    #...Set an error tolerance and a limit on the number of iterations:
    tol = 1.0e-8
    nmax = 5000

    #...Iterate on Equation 5.45 until z is determined to within the
    #...error tolerance:
    n = 0
    z = abs(z)
    ratio = 1.0
    while (abs(ratio) > tol) and (n < nmax):
        n += 1
        ratio = F(z, t) / dFdz(z)
        z = z - ratio

    if n >= nmax:
        print(f"** Number of iterations ({nmax}) exceeded **")

    # Equation 5.46a
    f    = 1.0 - (y(z) / r1)

    # Equation 5.46b
    g    = A * np.sqrt(y(z) / mu)

    # Equation 5.46d
    gdot = 1.0 - (y(z) / r2)

    # Equation 5.28
    V1 = (1.0 / g) * (R2 - f * R1)

    # Equation 5.29
    V2 = (1.0 / g) * (gdot * R2 - R1)

    return V1, V2