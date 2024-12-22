import numpy as np
import gibbs
import coe_from_sv

def gibbs_test():
    '''
    This program uses the Gibbs method and orbital element computation
    to obtain orbital elements from three coplanar position vectors.

    - Gravitational parameter (mu): km^3/s^2
    - r1, r2, r3: Three coplanar geocentric position vectors (km)
    - v2: Velocity corresponding to r2 (km/s)
    - coe: Orbital elements [h, e, RA, incl, w, TA, a]
        h: Angular momentum (km^2/s)
        e: Eccentricity
        RA: Right ascension of the ascending node (rad)
        incl: Orbit inclination (rad)
        w: Argument of perigee (rad)
        TA: True anomaly (rad)
        a: Semimajor axis (km)
        T: Period of elliptical orbit (s)

    User py-functions required: gibbs, coe_from_sv
    '''
    deg = np.pi / 180
    mu = 398600.4418  # Gravitational parameter (km^3/s^2)

    # Input data
    r1 = np.array([-294.32, 4265.1, 5986.7])
    r2 = np.array([-1365.5, 3637.6, 6346.8])
    r3 = np.array([-2940.3, 2473.7, 6555.8])

    # Echo the input data
    print('-----------------------------------------------------')
    print('\n\n Input data:\n')
    print(f'\n Gravitational parameter (km^3/s^2) = {mu}')
    print(f'\n r1 (km) = {r1}')
    print(f'\n r2 (km) = {r2}')
    print(f'\n r3 (km) = {r3}')
    print('\n\n')

    # Algorithm 5.1: Gibbs method
    v2, ierr = gibbs.gibbs(r1, r2, r3, mu)

    # If the vectors r1, r2, r3 are not coplanar, abort
    if ierr == 1:
        print('\n These vectors are not coplanar.\n\n')
        return

    # Algorithm 4.2: Compute orbital elements
    coe = coe_from_sv.coe_from_sv(r2, v2, mu)

    h    = coe[0]
    e    = coe[1]
    RA   = coe[2]
    incl = coe[3]
    w    = coe[4]
    TA   = coe[5]
    a    = coe[6]

    # Output the results
    print(' Solution:')
    print('\n')
    print(f'\n v2 (km/s) = {v2}')
    print('\n\n Orbital elements:')
    print(f'\n Angular momentum (km^2/s)  = {h}')
    print(f'\n Eccentricity               = {e}')
    print(f'\n Inclination (deg)          = {incl / deg}')
    print(f'\n RA of ascending node (deg) = {RA / deg}')
    print(f'\n Argument of perigee (deg)  = {w / deg}')
    print(f'\n True anomaly (deg)         = {TA / deg}')
    print(f'\n Semimajor axis (km)        = {a}')

    # If the orbit is an ellipse, output the period
    if e < 1:
        T = 2 * np.pi / np.sqrt(mu) * a ** 1.5
        print(f'\n Period (s)                 = {T}')
    print('\n-----------------------------------------------------\n')

if __name__ == '__main__':
    gibbs_test()