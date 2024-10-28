# Algorithm 5.4: Calculation of the state vector from measurements of range,
# position, and their rates.
import numpy as np

def rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H):
    """
    This function calculates the geocentric equatorial position and
    velocity vectors of an object from radar observations of range,
    azimuth, elevation angle and their rates.

    deg     - conversion factor between degrees and radians
    pi      - 3.1415926...
    Re      - equatorial radius of the earth (km)
    f       - earth’s flattening factor
    wE      - angular velocity of the earth (rad/s)
    omega   - earth’s angular velocity vector (rad/s) in the
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
    dec     - topocentric equatorial declination of object (rad)
    decdot  - declination rate (rad/s)
    h       - hour angle of object (rad)
    RA      - topocentric equatorial right ascension of object (rad)
    RAdot   - right ascension rate (rad/s)
    Rho     - unit vector from site to object
    Rhodot  - time rate of change of Rho (1/s)
    r       - geocentric equatorial position vector of object (km)
    v       - geocentric equatorial velocity vector of object (km)

    User py-functions required: none
    """
    
    global f, Re, wE
    deg = np.pi / 180
    omega = np.array([0, 0, wE])

    # Convert angular quantities from degrees to radians:
    A = A * deg
    Adot = Adot * deg
    a = a * deg
    adot = adot * deg
    theta = theta * deg
    phi = phi * deg

    # Equation 5.56:
    R = np.array([
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.cos(theta),
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.sin(theta),
        (Re * (1 - f)**2 / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.sin(phi)
    ])

    # Equation 5.66:
    Rdot = np.cross(omega, R)

    # Equation 5.83a:
    dec = np.arcsin(np.cos(phi) * np.cos(A) * np.cos(a) + np.sin(phi) * np.sin(a))

    # Equation 5.83b:
    h = np.arccos((np.cos(phi) * np.sin(a) - np.sin(phi) * np.cos(A) * np.cos(a)) / np.cos(dec))
    if (A > 0) and (A < np.pi):
        h = 2 * np.pi - h

    # Equation 5.83c:
    RA = theta - h

    # Equations 5.57:
    Rho = np.array([np.cos(RA) * np.cos(dec), np.sin(RA) * np.cos(dec), np.sin(dec)])

    # Equation 5.63:
    r = R + rho * Rho

    # Equation 5.84:
    decdot = (-Adot * np.cos(phi) * np.sin(A) * np.cos(a) + adot * (np.sin(phi) * np.cos(a)
              - np.cos(phi) * np.cos(A) * np.sin(a))) / np.cos(dec)

    # Equation 5.85:
    RAdot = (wE
             + (Adot * np.cos(A) * np.cos(a) - adot * np.sin(A) * np.sin(a)
                + decdot * np.sin(A) * np.cos(a) * np.tan(dec))
             / (np.cos(phi) * np.sin(a) - np.sin(phi) * np.cos(A) * np.cos(a)))

    # Equations 5.69 and 5.72:
    Rhodot = np.array([
        -RAdot * np.sin(RA) * np.cos(dec) - decdot * np.cos(RA) * np.sin(dec),
        RAdot * np.cos(RA) * np.cos(dec) - decdot * np.sin(RA) * np.sin(dec),
        decdot * np.cos(dec)
    ])

    # Equation 5.64:
    v = Rdot + rhodot * Rho + rho * Rhodot

    return r, v
