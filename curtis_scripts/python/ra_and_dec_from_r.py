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
    l = r[0] / np.linalg.norm(r)
    m = r[1] / np.linalg.norm(r)
    n = r[2] / np.linalg.norm(r)

    dec = np.degrees(np.arcsin(n))

    if m > 0:
        ra = np.degrees(np.arccos(l / np.cos(np.radians(dec))))
    else:
        ra = 360 - np.degrees(np.arccos(l / np.cos(np.radians(dec))))

    return ra, dec