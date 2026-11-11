import numpy as np

def coe_from_sv(R, V, mu):
    '''
    This function computes the classical orbital elements (coe)
    from the state vector (R, V) using Algorithm 4.1.

    mu      - gravitational parameter (km^3/s^2)
    R       - position vector in the geocentric equatorial frame (km)
    V       - velocity vector in the geocentric equatorial frame (km/s)
    r, v    - the magnitudes of R and V
    vr      - radial velocity component (km/s)
    H       - the angular momentum vector (km^2/s)
    h       - the magnitude of H (km^2/s)
    incl    - inclination of the orbit (rad)
    N       - the node line vector (km^2/s)
    n       - the magnitude of N
    cp      - cross product of N and R
    RA      - right ascension of the ascending node (rad)
    E       - eccentricity vector
    e       - eccentricity (magnitude of E)
    eps     - a small number below which the eccentricity is considered
              to be zero
    w       - argument of perigee (rad)
    TA      - true anomaly (rad)
    a       - semimajor axis (km)
    pi      - 3.1415926...
    coe     - vector of orbital elements [h e RA incl w TA a]
    
    User py-functions required: None
    '''
    eps = 1.e-10
    R = np.asarray(R, dtype=float)
    V = np.asarray(V, dtype=float)
    if R.shape != (3,) or V.shape != (3,) or not np.isfinite(np.r_[R, V, mu]).all() or mu <= 0:
        raise ValueError("Finite three-vectors and a positive mu are required.")

    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    if r == 0:
        raise ValueError("Classical elements require a nonzero position vector.")
    vr = np.dot(R, V) / r

    H = np.cross(R, V)
    h = np.linalg.norm(H)
    if h == 0:
        raise ValueError("Classical elements require nonzero angular momentum.")

    #...Equation 4.7:
    incl = np.arctan2(np.hypot(H[0], H[1]), H[2])

    #...Equation 4.8:
    N = np.cross([0, 0, 1], H)
    n = np.linalg.norm(N)

    #...Equation 4.9:
    if n > eps * h:
        RA = np.arctan2(N[1], N[0]) % (2 * np.pi)
    else:
        RA = 0

    #...Equation 4.10:
    E = (1 / mu) * ((v**2 - mu / r) * R - r * vr * V)
    e = np.linalg.norm(E)

    #...Equation 4.12 (incorporating the case e = 0):
    if n > eps * h:
        if e > eps:
            w = np.arctan2(np.dot(np.cross(N, E), H) / h, np.dot(N, E)) % (2 * np.pi)
        else:
            w = 0
    else:
        if e > eps:
            w = np.arctan2(np.sign(H[2]) * E[1], E[0]) % (2 * np.pi)
        else:
            w = 0

    #...Equation 4.13a (incorporating the case e = 0):
    if e > eps:
        TA = np.arctan2(np.dot(np.cross(E, R), H) / h, np.dot(E, R)) % (2 * np.pi)
    else:
        if n > eps * h:
            cp = np.cross(N, R)
            TA = np.arctan2(np.dot(cp, H) / h, np.dot(N, R)) % (2 * np.pi)
        else:
            TA = np.arctan2(np.sign(H[2]) * R[1], R[0]) % (2 * np.pi)

    #...Equation 4.62 (a < 0 for a hyperbola):
    den = 1 - e**2
    if abs(den) <= 32 * np.finfo(float).eps:
        a = np.inf
    else:
        a = h**2 / mu / den
    coe = [h, e, RA, incl, w, TA, a]

    return coe