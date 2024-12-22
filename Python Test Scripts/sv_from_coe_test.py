import numpy as np
import sv_from_coe

def sv_from_coe_test():
    '''
    This program uses Algorithm 4.5 to obtain the state vector from
    the orbital elements.

    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    coe - orbital elements [h e RA incl w TA a]
          where h = angular momentum (km^2/s)
          e       = eccentricity
          RA      = right ascension of the ascending node (rad)
          incl    = orbit inclination (rad)
          w       = argument of perigee (rad)
          TA      = true anomaly (rad)
          a       = semimajor axis (km)
    r   - position vector (km) in geocentric equatorial frame
    v   - velocity vector (km) in geocentric equatorial frame

    User py-function required: sv_from_coe
    '''
    deg = np.pi / 180
    mu = 398600.4418

    # Data declaration (angles in degrees):
    h = 80000
    e = 1.4
    RA = 40
    incl = 30
    w = 60
    TA = 30

    # Orbital elements (convert angles to radians)
    coe = [h, e, RA * deg, incl * deg, w * deg, TA * deg]

    # Algorithm 4.5 (requires angular elements be in radians):
    r, v = sv_from_coe.sv_from_coe(coe, mu)

    # Echo the input data and output the results to the console:
    print('-----------------------------------------------------')
    print(f'\n Gravitational parameter (km^3/s^2)  = {mu}')
    print(f'\n Angular momentum (km^2/s)           = {h}')
    print(f'\n Eccentricity                        = {e}')
    print(f'\n Right ascension (deg)               = {RA}')
    print(f'\n Argument of perigee (deg)           = {w}')
    print(f'\n True anomaly (deg)                  = {TA}')
    print('\n\n State vector:')
    print(f'\n   r (km)   = [{r[0]} {r[1]} {r[2]}]')
    print(f'\n   v (km/s) = [{v[0]} {v[1]} {v[2]}]')
    print('-----------------------------------------------------')

if __name__ == '__main__':
    sv_from_coe_test()