# ALGORITHM 5.1: GIBBS' METHOD OF PRELIMINARY ORBIT DETERMINATION

import numpy as np

def gibbs(R1, R2, R3, mu):
    '''
    This function uses the Gibbs method of orbit determination to
    compute the velocity corresponding to the second of three
    supplied position vectors.

    mu              - gravitational parameter (km^3/s^2)
    R1, R2, R3      - three coplanar geocentric position vectors (km)
    r1, r2, r3      - the magnitudes of R1, R2 and R3 (km)
    c12, c23, c31   - three independent cross products among
                      R1, R2 and R3
    N, D, S         - vectors formed from R1, R2 and R3 during
                      the Gibbs  procedure
    tol             - tolerance for determining if R1, R2 and R3
                      are coplanar
    ierr            - = 0 if R1, R2, R3 are found to be coplanar
                      = 1 otherwise
    V2              - the velocity corresponding to R2 (km/s)

    User py-functions required: none
    '''
    tol = 1e-4
    ierr = 0

    # Magnitudes of R1, R2 and R3:
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)
    r3 = np.linalg.norm(R3)

    # Cross products among R1, R2 and R3:
    c12 = np.cross(R1, R2)
    c23 = np.cross(R2, R3)
    c31 = np.cross(R3, R1)

    # Check that R1, R2 and R3 are coplanar; if not set error flag:
    if abs(np.dot(R1, c23) / r1 / np.linalg.norm(c23)) > tol:
        ierr = 1

    # Equation 5.13:
    N = r1 * c23 + r2 * c31 + r3 * c12

    # Equation 5.14:
    D = c12 + c23 + c31

    # Equation 5.21:
    S = R1 * (r2 - r3) + R2 * (r3 - r1) + R3 * (r1 - r2)

    # Equation 5.22:
    V2 = np.sqrt(mu / np.linalg.norm(N) / np.linalg.norm(D)) * (np.cross(D, R2) / r2 + S)

    return V2, ierr
