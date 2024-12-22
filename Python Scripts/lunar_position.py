# ALGORITHM 10.4: CALCULATE THE GEOCENTRIC POSITION OF THE
# MOON AT A GIVEN EPOCH

import numpy as np

def lunar_position(jd):
    '''
    Calculates the geocentric equatorial position vector of the moon
    given the Julian day.
    
    User py-functions required: None
    '''
    # Earth's radius (km)
    RE = 6378.14

    # Time in centuries since J2000
    T = (jd - 2451545) / 36525

    # Ecliptic longitude (deg)
    e_long = 218.32 + 481267.881 * T \
             + 6.29 * np.sin(np.radians(135.0 + 477198.87 * T)) \
             - 1.27 * np.sin(np.radians(259.3 - 413335.36 * T)) \
             + 0.66 * np.sin(np.radians(235.7 + 890534.22 * T)) \
             + 0.21 * np.sin(np.radians(269.9 + 954397.74 * T)) \
             - 0.19 * np.sin(np.radians(357.5 + 35999.05 * T)) \
             - 0.11 * np.sin(np.radians(186.5 + 966404.03 * T))
    e_long = np.mod(e_long, 360)

    # Ecliptic latitude (deg)
    e_lat = 5.13 * np.sin(np.radians(93.3 + 483202.02 * T)) \
            + 0.28 * np.sin(np.radians(228.2 + 960400.89 * T)) \
            - 0.28 * np.sin(np.radians(318.3 + 6003.15 * T)) \
            - 0.17 * np.sin(np.radians(217.6 - 407332.21 * T))
    e_lat = np.mod(e_lat, 360)

    # Horizontal parallax (deg)
    h_par = 0.9508 \
            + 0.0518 * np.cos(np.radians(135.0 + 477198.87 * T)) \
            + 0.0095 * np.cos(np.radians(259.3 - 413335.36 * T)) \
            + 0.0078 * np.cos(np.radians(235.7 + 890534.22 * T)) \
            + 0.0028 * np.cos(np.radians(269.9 + 954397.74 * T))
    h_par = np.mod(h_par, 360)

    # Angle between earth’s orbit and its equator (deg)
    obliquity = 23.439291 - 0.0130042 * T

    # Direction cosines of the moon’s geocentric equatorial position vector
    l = np.cos(np.radians(e_lat)) * np.cos(np.radians(e_long))
    m = np.cos(np.radians(obliquity)) * np.cos(np.radians(e_lat)) * np.sin(np.radians(e_long)) \
        - np.sin(np.radians(obliquity)) * np.sin(np.radians(e_lat))
    n = np.sin(np.radians(obliquity)) * np.cos(np.radians(e_lat)) * np.sin(np.radians(e_long)) \
        + np.cos(np.radians(obliquity)) * np.sin(np.radians(e_lat))

    # Earth-moon distance (km)
    dist = RE / np.sin(np.radians(h_par))

    # Moon’s geocentric equatorial position vector (km)
    r_moon = dist * np.array([l, m, n])

    return r_moon