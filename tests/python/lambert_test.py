import numpy as np
import lambert
import coe_from_sv

def lambert_test():
    '''
    This program uses Algorithm 5.2 to solve Lambert's problem for the
    data provided.

    deg - factor for converting between degrees and radians
    pi - 3.1415926...
    mu - gravitational parameter (km^3/s^2)
    r1, r2 - initial and final position vectors (km)
    dt - time between r1 and r2 (s)
    string - = 'pro' if the orbit is prograde
    = 'retro if the orbit is retrograde
    v1, v2 - initial and final velocity vectors (km/s)
    coe - orbital elements [h e RA incl w TA a]
    where h = angular momentum (km^2/s)
    e = eccentricity
    RA = right ascension of the ascending node (rad)
    incl = orbit inclination (rad)
    w = argument of perigee (rad)
    TA = true anomaly (rad)
    a = semimajor axis (km)
    TA1 - Initial true anomaly (rad)
    TA2 - Final true anomaly (rad)
    T - period of an elliptic orbit (s)

    User py-functions required: lambert, coe_from_sv
    '''
    mu = 398600.4418

    deg = np.pi / 180.0

    #...Data declaration:
    r1 = np.array([5000.0, 10000.0, 2100.0])
    r2 = np.array([-14600.0, 2500.0, 7000.0])
    dt = 3600.0
    string = 'pro'

    # ...Algorithm 5.2:
    v1, v2 = lambert.lambert(r1, r2, dt, string, mu)

    # ...Algorithm 4.1 (using r1 and v1):
    coe = coe_from_sv.coe_from_sv(r1, v1, mu)
    # ...Save the initial true anomaly:
    TA1 = coe[5]

    # ...Algorithm 4.1 (using r2 and v2):
    coe = coe_from_sv.coe_from_sv(r2, v2, mu)
    # ...Save the final true anomaly:
    TA2 = coe[5]

    # ...Echo the input data and output the results:
    print('-----------------------------------------------------')
    print('\n\n Input data:\n')

    print(f'   Gravitational parameter (km^3/s^2) = {mu}')
    print(f'\n   r1 (km)                       = [{r1[0]} {r1[1]} {r1[2]}]')
    print(f'   r2 (km)                       = [{r2[0]} {r2[1]} {r2[2]}]')
    print(f'   Elapsed time (s)              = {dt}')

    print('\n\n Solution:\n')
    print(f'   v1 (km/s)                     = [{v1[0]} {v1[1]} {v1[2]}]')
    print(f'   v2 (km/s)                     = [{v2[0]} {v2[1]} {v2[2]}]')

    print('\n\n Orbital elements:')
    # coe = [h, e, RA, incl, w, TA, a]  (typical structure)
    print(f'   Angular momentum (km^2/s)     = {coe[0]}')
    print(f'   Eccentricity                  = {coe[1]}')
    print(f'   Inclination (deg)             = {coe[3] / deg}')
    print(f'   RA of ascending node (deg)    = {coe[2] / deg}')
    print(f'   Argument of perigee (deg)     = {coe[4] / deg}')
    print(f'   True anomaly initial (deg)    = {TA1 / deg}')
    print(f'   True anomaly final (deg)      = {TA2 / deg}')
    print(f'   Semimajor axis (km)           = {coe[6]}')

    # Periapse radius = h^2 / mu / (1+e)
    rp = (coe[0] ** 2) / mu / (1.0 + coe[1])
    print(f'   Periapse radius (km)          = {rp}')

    # ...If the orbit is an ellipse, output its period:
    if coe[1] < 1.0:
        # T = 2*pi/sqrt(mu)*coe(7)^1.5
        T = 2.0 * np.pi / np.sqrt(mu) * (coe[6] ** 1.5)

        print('\n   Period:')
        print(f'     Seconds                 = {T}')
        print(f'     Minutes                 = {T/60.0}')
        print(f'     Hours                   = {T/3600.0}')
        print(f'     Days                    = {T/86400.0}')

    print('\n-----------------------------------------------------\n')

if __name__ == '__main__':
    lambert_test()
