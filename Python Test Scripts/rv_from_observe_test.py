import numpy as np
import rv_from_observe
import coe_from_sv

def rv_from_observe_test():
    '''
    This program uses Algorithms 5.4 and 4.2 to obtain the orbital
    elements from the observational data.
    
    deg    - conversion factor between degrees and radians
    pi     - 3.1415926...
    mu     - gravitational parameter (km^3/s^2)

    Re     - equatorial radius of the earth (km)
    f      - earth's flattening factor
    wE     - angular velocity of the earth (rad/s)
    omega  - earth's angular velocity vector (rad/s) in the
             geocentric equatorial frame

    rho    - slant range of object (km)
    rhodot - range rate (km/s)
    A      - azimuth (deg) of object relative to observation site
    Adot   - time rate of change of azimuth (deg/s)
    a      - elevation angle (deg) of object relative to observation site
    adot   - time rate of change of elevation angle (degrees/s)

    theta  - local sidereal time (deg) of tracking site
    phi    - geodetic latitude (deg) of site
    H      - elevation of site (km)

    r      - geocentric equatorial position vector of object (km)
    v      - geocentric equatorial velocity vector of object (km)

    coe    - orbital elements [h e RA incl w TA a]
             where
                 h    = angular momentum (km^2/s)
                 e    = eccentricity
                 RA   = right ascension of the ascending node (rad)
                 incl = inclination of the orbit (rad)
                 w    = argument of perigee (rad)
                 TA   = true anomaly (rad)
                 a    = semimajor axis (km)
    rp     - perigee radius (km)
    T      - period of elliptical orbit (s)
    '''
    deg = np.pi / 180
    f   = 1 / 298.256421867
    Re  = 6378.13655
    wE  = 7.292115e-5
    mu  = 398600.4418

    #...Data declaration:
    rho    = 2551
    rhodot = 0
    A      = 90
    Adot   = 0.1130
    a      = 30
    adot   = 0.05651
    theta  = 300
    phi    = 60
    H      = 0
    #...

    #...Algorithm 5.4:
    r, v = rv_from_observe.rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H)

    #...Algorithm 4.2:
    coe = coe_from_sv.coe_from_sv(r, v, mu)

    h    = coe[0]
    e    = coe[1]
    RA   = coe[2]
    incl = coe[3]
    w    = coe[4]
    TA   = coe[5]
    a    = coe[6]

    #...Equation 2.40:
    rp = h**2 / mu / (1 + e)

    # Echo the input data and output the solution to the console:
    print('-----------------------------------------------------')
    print('\n Input data:')
    print(f'\n Slant range (km)              = {rho}')
    print(f' Slant range rate (km/s)       = {rhodot}')
    print(f' Azimuth (deg)                 = {A}')
    print(f' Azimuth rate (deg/s)          = {Adot}')
    print(f' Elevation (deg)               = {a}')
    print(f' Elevation rate (deg/s)        = {adot}')
    print(f' Local sidereal time (deg)     = {theta}')
    print(f' Latitude (deg)                = {phi}')
    print(f' Altitude above sea level (km) = {H}')
    print('\n')

    print(' Solution:')

    print('\n State vector:')
    print(f'\n r (km)                        = [{r[0]}, {r[1]}, {r[2]}]')
    print(f' v (km/s)                      = [{v[0]}, {v[1]}, {v[2]}]')

    print('\n Orbital elements:')
    print(f'\n   Angular momentum (km^2/s)   = {h}')
    print(f'   Eccentricity                = {e}')
    print(f'   Inclination (deg)           = {incl / deg}')
    print(f'   RA of ascending node (deg)  = {RA / deg}')
    print(f'   Argument of perigee (deg)   = {w / deg}')
    print(f'   True anomaly (deg)          = {TA / deg}')
    print(f'   Semimajor axis (km)         = {a}')
    print(f'   Perigee radius (km)         = {rp}')

    # If the orbit is an ellipse, output its period:
    if e < 1:
        T = 2 * np.pi / np.sqrt(mu) * a**1.5
        print(f'\n   Period:')
        print(f'     Seconds = {T}')
        print(f'     Minutes = {T / 60}')
        print(f'     Hours   = {T / 3600}')
        print(f'     Days    = {T / 24 / 3600}')

    print('-----------------------------------------------------')

if __name__ == "__main__":
    rv_from_observe_test()
