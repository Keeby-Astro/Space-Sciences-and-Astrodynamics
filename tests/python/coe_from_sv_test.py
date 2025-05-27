import numpy as np
import coe_from_sv

def coe_from_sv_test():
    '''
    This program uses Algorithm 4.2 to obtain the orbital
    elements from the state vector.
    
    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    r   - position vector (km) in the geocentric equatorial frame
    v   - velocity vector (km/s) in the geocentric equatorial frame
    coe - orbital elements [h e RA incl w TA a]
        where h    = angular momentum (km^2/s)
                e    = eccentricity
                RA   = right ascension of the ascending node (rad)
                incl = orbit inclination (rad)
                w    = argument of perigee (rad)
                TA   = true anomaly (rad)
                a    = semimajor axis (km)
    T   - Period of an elliptic orbit (s)

    User py-function required: coe_from_sv
    '''
    deg = np.pi/180
    mu = 398600.4418

    #...Data declaration
    r = np.array([-6045, -3490, 2500])
    v = np.array([-3.457, 6.618, 2.533])

    #...Algorithm 4.2:
    coe = coe_from_sv.coe_from_sv(r, v, mu)

    #...Echo the input data and output results to the console
    print('-----------------------------------------------------')
    print(f'\n Gravitational parameter (km^3/s^2) = {mu}')
    print('\n State vector:')
    print(f'\n r (km)                            = [{r[0]} {r[1]} {r[2]}]')
    print(f'\n v (km/s)                          = [{v[0]} {v[1]} {v[2]}]')
    print(' ')
    print(f'\nAngular momentum (km^2/s)          = {coe[0]}')
    print(f'Eccentricity                       = {coe[1]}')
    print(f'Right ascension (deg)              = {coe[2] / deg}')
    print(f'Inclination (deg)                  = {coe[3] / deg}')
    print(f'Argument of perigee (deg)          = {coe[4] / deg}')
    print(f'True anomaly (deg)                 = {coe[5] / deg}')
    print(f'Semimajor axis (km):               = {coe[6]}')

    #...If the orbit is an ellipse, output its period (Equation 2.73)
    if coe[1] < 1:
        T = 2 * np.pi / np.sqrt(mu) * coe[6] ** 1.5
        print('\nPeriod:')
        print(f'Seconds = {T}')
        print(f'Minutes = {T / 60}')
        print(f'Hours   = {T / 3600}')
        print(f'Days    = {T / (24 * 3600)}')
    print('-----------------------------------------------------\n')

if __name__ == "__main__":
    coe_from_sv_test()
