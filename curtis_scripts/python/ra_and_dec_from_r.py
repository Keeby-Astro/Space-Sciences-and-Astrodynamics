# ALGORITHM 4.1: OBTAIN THE RIGHT ASCENSION AND DECLINATION
# FROM THE POSITION VECTOR

import numpy as np

def ra_and_dec_from_r(r):
    '''
    This function calculates the right ascension and the
    declination from the geocentric equatorial position vector.

    r       - position vector
    l, m, n - direction cosines of r
    ra      - right ascension (degrees)
    dec     - declination (degrees)
    '''
    r = np.asarray(r, dtype=float)
    if r.shape != (3,) or not np.isfinite(r).all() or np.linalg.norm(r) == 0:
        raise ValueError("Position must be a finite, nonzero three-vector.")

    l = r[0] / np.linalg.norm(r)
    m = r[1] / np.linalg.norm(r)
    n = r[2] / np.linalg.norm(r)

    dec = np.degrees(np.arctan2(n, np.hypot(l, m)))

    if m >= 0:
        ra = np.degrees(np.arctan2(m, l))
    else:
        ra = 360 + np.degrees(np.arctan2(m, l))
    ra = ra % 360

    return ra, dec