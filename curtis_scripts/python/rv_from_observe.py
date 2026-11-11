# ALGORITHM 5.4: CALCULATION OF THE STATE VECTOR FROM
# MEASUREMENTS OF RANGE, ANGULAR POSITION, AND THEIR RATES

import numpy as np

def rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H,
                    f=1/298.26, Re=6378.14, wE=7.292115e-5):
    '''
    This function calculates the geocentric equatorial position and
    velocity vectors of an object from radar observations of range,
    azimuth, elevation angle and their rates.

    deg     - conversion factor between degrees and radians
    pi      - 3.1415926...
    Re      - equatorial radius of the earth (km)
    f       - earth's flattening factor
    wE      - angular velocity of the earth (rad/s)
    omega   - earth's angular velocity vector (rad/s) in the
              geocentric equatorial frame
    theta   - local sidereal time (degrees) of tracking site
    phi     - geodetic latitude (degrees) of site
    H       - elevation of site (km)
    R       - geocentric equatorial position vector (km) of tracking site
    Rdot    - inertial velocity (km/s) of site
    rho     - slant range of object (km)
    rhodot  - range rate (km/s)
    A       - azimuth (degrees) of object relative to observation site
    Adot    - time rate of change of azimuth (degrees/s)
    a       - elevation angle (degrees) of object relative to observation site
    adot    - time rate of change of elevation angle (degrees/s)
    East    - local east unit vector in the geocentric equatorial frame
    North   - local north unit vector in the geocentric equatorial frame
    Up      - local up unit vector in the geocentric equatorial frame
    Rho     - unit vector from site to object
    Rhodot  - time rate of change of Rho (1/s)
    r       - geocentric equatorial position vector of object (km)
    v       - geocentric equatorial velocity vector of object (km/s)

    User py-functions required: none
    '''
    deg = np.pi / 180
    omega = np.array([0, 0, wE])

    #...Convert angular quantities from degrees to radians:
    A = A * deg
    Adot = Adot * deg
    a = a * deg
    adot = adot * deg
    theta = theta * deg
    phi = phi * deg

    #...Equation 5.56:
    R = np.array([
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.cos(theta),
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.sin(theta),
        (Re * (1 - f)**2 / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.sin(phi)
    ])

    #...Equation 5.66:
    Rdot = np.cross(omega, R)

    #...Resolve and differentiate the sight line in the east, north and up basis
    East = np.array([-np.sin(theta), np.cos(theta), 0])
    North = np.array([-np.sin(phi)*np.cos(theta), -np.sin(phi)*np.sin(theta), np.cos(phi)])
    Up = np.array([np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta), np.sin(phi)])
    Rho = np.cos(a)*np.sin(A)*East + np.cos(a)*np.cos(A)*North + np.sin(a)*Up
    Rhodot = (np.cross(omega, Rho)
              + Adot*np.cos(a)*(np.cos(A)*East - np.sin(A)*North)
              + adot*(-np.sin(a)*np.sin(A)*East - np.sin(a)*np.cos(A)*North + np.cos(a)*Up))

    #...Equations 5.69 and 5.72:
    r = R + rho*Rho

    #...Equation 5.64:
    v = Rdot + rhodot*Rho + rho*Rhodot
    return r, v