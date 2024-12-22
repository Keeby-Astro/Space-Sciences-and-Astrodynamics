# ALGORITHM 10.2: CALCULATE THE GEOCENTRIC POSITION OF THE SUN AT A GIVEN EPOCH

import numpy as np

def solar_position(jd):
    '''
    This function alculates the geocentric equatorial position vector
    of the sun, given the julian date.
    
    User py-functions required: None
    '''
    # Astronomical unit (km)
    AU = 149597870.691

    # Julian days since J2000
    n = jd - 2451545

    # Julian centuries since J2000
    cy = n / 36525

    # Mean anomaly (deg)
    M = 357.528 + 0.9856003 * n
    M = M % 360

    # Mean longitude (deg)
    L = 280.460 + 0.98564736 * n
    L = L % 360

    # Apparent ecliptic longitude (deg)
    lamda = L + 1.915 * np.sin(np.radians(M)) + 0.020 * np.sin(np.radians(2 * M))
    lamda = lamda % 360

    # Obliquity of the ecliptic (deg)
    eps = 23.439 - 0.0000004 * n

    # Unit vector from earth to sun
    u = np.array([np.cos(np.radians(lamda)), 
                  np.sin(np.radians(lamda)) * np.cos(np.radians(eps)), 
                  np.sin(np.radians(lamda)) * np.sin(np.radians(eps))])

    # Distance from earth to sun (km)
    rS = (1.00014 - 0.01671 * np.cos(np.radians(M)) - 0.000140 * np.cos(np.radians(2 * M))) * AU

    # Geocentric position vector (km)
    r_S = rS * u

    return lamda, eps, r_S