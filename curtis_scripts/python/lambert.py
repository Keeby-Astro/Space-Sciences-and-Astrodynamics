import numpy as np
import stumpC
import stumpS

def lambert(R1, R2, t, orbit_type, mu):
    '''
    This function solves the zero-revolution Lambert problem.

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
    tol     - relative tolerance on the time-of-flight residual
    nmax    - maximum number of bisection iterations
    f, g    - Lagrange coefficients
    gdot    - time derivative of g
              C(z), S(z) - Stumpff functions
    dum     - a dummy variable

    User py-functions required: stumpC and stumpS
    '''
    R1 = np.asarray(R1, dtype=float)
    R2 = np.asarray(R2, dtype=float)
    if R1.shape != (3,) or R2.shape != (3,) or not np.isfinite(np.r_[R1, R2, t, mu]).all() or mu <= 0 or t <= 0:
        raise ValueError("Lambert requires finite three-vectors and positive time and mu.")

    #...Magnitudes of R1 and R2:
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    c12 = np.cross(R1, R2)
    if r1 == 0 or r2 == 0 or np.linalg.norm(c12) <= 1.e-12 * r1 * r2:
        raise ValueError("Collinear or zero position vectors require a separate Lambert formulation.")
    theta = np.arccos(np.clip(np.dot(R1, R2) / (r1 * r2), -1, 1))

    #...Determine whether the orbit is prograde or retrograde:
    if orbit_type not in ["pro", "retro"]:
        raise ValueError("Orbit type must be 'pro' or 'retro'.")

    if orbit_type == "pro":
        if c12[2] <= 0.0:
            theta = 2.0 * np.pi - theta
    elif orbit_type == "retro":
        if c12[2] >= 0.0:
            theta = 2.0 * np.pi - theta

    #...Equation 5.35
    A = np.sin(theta) * np.sqrt(r1 * r2 / (1.0 - np.cos(theta)))

    def C(z):
        '''Stumpff C function.'''
        return stumpC.stumpC(z)

    def S(z):
        '''Stumpff S function.'''
        return stumpS.stumpS(z)

    #...Equation 5.38
    def y(z):
        return r1 + r2 + A * ((z * S(z)) - 1.0) / np.sqrt(C(z))

    #...Equation 5.40:
    def F(z, t_):
        yz = y(z)
        if yz <= 0:
            return -np.sqrt(mu) * t_
        return ((yz / C(z)) ** 1.5) * S(z) + A * np.sqrt(yz) - np.sqrt(mu) * t_

    #...Bracket the zero-revolution branch, z < 4*pi^2:
    lo = -1.0
    hi = 4 * np.pi**2 * (1 - 1.e-7)
    for n in range(100):
        value = F(lo, t)
        if not np.isfinite(value):
            raise RuntimeError("Lambert lower bracket overflowed.")
        if value <= 0:
            break
        lo = 2 * lo
    else:
        raise RuntimeError("Could not bracket Lambert's solution.")
    if not np.isfinite(F(hi, t)) or F(hi, t) <= 0:
        raise RuntimeError("Time of flight exceeds this solver's zero-revolution bracket.")

    tol = 1.e-11
    nmax = 200

    #...Bisect the bracket until the time-of-flight residual satisfies the tolerance:
    for n in range(nmax):
        z = (lo + hi) / 2
        value = F(z, t)
        if not np.isfinite(value):
            raise RuntimeError("Lambert iteration is nonfinite.")
        if abs(value) <= tol * np.sqrt(mu) * t:
            break
        if value > 0:
            hi = z
        else:
            lo = z
    else:
        raise RuntimeError("Lambert solver failed to converge.")

    yz = y(z)
    if yz <= 0:
        raise RuntimeError("Lambert solution lies outside the physical domain.")

    #...Equation 5.46a
    f    = 1.0 - (yz / r1)

    #...Equation 5.46b
    g    = A * np.sqrt(yz / mu)

    #...Equation 5.46d
    gdot = 1.0 - (yz / r2)

    #...Equation 5.28
    V1 = (1.0 / g) * (R2 - f * R1)

    #...Equation 5.29
    V2 = (1.0 / g) * (gdot * R2 - R1)

    return V1, V2