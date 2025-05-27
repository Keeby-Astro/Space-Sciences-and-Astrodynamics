import numpy as np
import gauss
import coe_from_sv

def gauss_test():
    '''
    This program uses Algorithms 5.5 and 5.6 (Gauss's method) to compute
    the state vector from the data provided.

    deg          - factor for converting between degrees and radians
    pi           - 3.1415926...
    mu           - gravitational parameter (km^3/s^2)
    Re           - earth's radius (km)
    f            - earth's flattening factor
    H            - elevation of observation site (km)
    phi          - latitude of site (deg)
    t            - vector of observation times t1, t2, t3 (s)
    ra           - vector of topocentric equatorial right ascensions
                    at t1, t2, t3 (deg)
    dec          - vector of topocentric equatorial right declinations
                    at t1, t2, t3 (deg)
    theta        - vector of local sidereal times for t1, t2, t3 (deg)
    R            - matrix of site position vectors at t1, t2, t3 (km)
    rho          - matrix of direction cosine vectors at t1, t2, t3
    fac1, fac2   - common factors
    r_old, v_old - the state vector without iterative improvement (km, km/s)
    r, v         - the state vector with iterative improvement (km, km/s)
    coe          - vector of orbital elements for r, v: 
                    [h, e, RA, incl, w, TA, a]
                    where h    = angular momentum (km^2/s)
                          e    = eccentricity
                          incl = inclination (rad)
                          w    = argument of perigee (rad)
                          TA   = true anomaly (rad)
                          a    = semimajor axis (km)
    coe_old     - vector of orbital elements for r_old, v_old  
    
    User py-functions required: gauss, coe_from_sv     
    '''
    deg = np.pi/180
    mu  = 398600.4418
    Re  = 6378.14
    f   = 1/298.26

    #...Data declaration:
    H     = 1
    phi   = 40 * deg
    t     = np.array([0, 118.104, 237.577])
    ra    = np.array([43.5365, 54.4196, 64.3178]) * deg
    dec   = np.array([-8.78334, -12.0739, -15.1054]) * deg
    theta = np.array([44.5065, 45.000, 45.4992]) * deg

    #...Equations 5.64, 5.76 and 5.79
    fac1 = Re / np.sqrt(1 - (2 * f - f * f) * np.sin(phi)**2)
    fac2 = (Re * (1 - f)**2 / np.sqrt(1 - (2 * f - f * f) * np.sin(phi)**2) + H) * np.sin(phi)

    R = np.zeros((3, 3))
    rho = np.zeros((3, 3))

    for i in range(3):
        R[i, 0] = (fac1 + H) * np.cos(phi) * np.cos(theta[i])
        R[i, 1] = (fac1 + H) * np.cos(phi) * np.sin(theta[i])
        R[i, 2] = fac2
        rho[i, 0] = np.cos(dec[i]) * np.cos(ra[i])
        rho[i, 1] = np.cos(dec[i]) * np.sin(ra[i])
        rho[i, 2] = np.sin(dec[i])

    #...Algorithms 5.5 and 5.6
    r, v, r_old, v_old = gauss.gauss(rho[0, :], rho[1, :], rho[2, :],
                                R[0, :], R[1, :], R[2, :],
                                t[0], t[1], t[2])

    #...Algorithm 4.2 for the initial estimate of the state vector
    #   and for the iteratively improved one:
    coe_old = coe_from_sv.coe_from_sv(r_old, v_old, mu)
    coe = coe_from_sv.coe_from_sv(r, v, mu)

    #...Echo the input data and output the solution
    print('-----------------------------------------------------')
    print(' Example 5.11: Orbit determination by the Gauss method\n')
    print(f' Radius of earth (km)               = {Re}')
    print(f' Flattening factor                  = {f}')
    print(f' Gravitational parameter (km^3/s^2) = {mu}\n')

    print(' Input data:')
    print(f' Latitude (deg)                = {phi / deg}')
    print(f' Altitude above sea level (km) = {H}\n')

    print(' Observations:')
    print('               Right                                     Local')
    print('   Time (s)   Ascension (deg)   Declination (deg)   Sidereal time (deg)')
    for i in range(3):
        print(f' {t[i]:9.4g} {ra[i] / deg:11.4f} {dec[i] / deg:19.4f} {theta[i] / deg:20.4f}')

    print('\n Solution:\n')

    print(' Without iterative improvement...\n')
    print(f' r (km)                          = [{r_old[0]}, {r_old[1]}, {r_old[2]}]')
    print(f' v (km/s)                        = [{v_old[0]}, {v_old[1]}, {v_old[2]}]\n')

    print(f'   Angular momentum (km^2/s)     = {coe_old[0]}')
    print(f'   Eccentricity                  = {coe_old[1]}')
    print(f'   RA of ascending node (deg)    = {coe_old[2] / deg}')
    print(f'   Inclination (deg)             = {coe_old[3] / deg}')
    print(f'   Argument of perigee (deg)     = {coe_old[4] / deg}')
    print(f'   True anomaly (deg)            = {coe_old[5] / deg}')
    print(f'   Semimajor axis (km)           = {coe_old[6]}')
    print(f'   Periapse radius (km)          = {coe_old[0]**2 / mu / (1 + coe_old[1])}')

    if coe_old[1] < 1:
        T = 2 * np.pi / np.sqrt(mu) * coe_old[6]**1.5
        print(f'   Period:')
        print(f'     Seconds                     = {T}')
        print(f'     Minutes                     = {T / 60}')
        print(f'     Hours                       = {T / 3600}')
        print(f'     Days                        = {T / 24 / 3600}')

    print('\n With iterative improvement...\n')
    print(f' r (km)                          = [{r[0]}, {r[1]}, {r[2]}]')
    print(f' v (km/s)                        = [{v[0]}, {v[1]}, {v[2]}]\n')

    print(f'   Angular momentum (km^2/s)     = {coe[0]}')
    print(f'   Eccentricity                  = {coe[1]}')
    print(f'   RA of ascending node (deg)    = {coe[2] / deg}')
    print(f'   Inclination (deg)             = {coe[3] / deg}')
    print(f'   Argument of perigee (deg)     = {coe[4] / deg}')
    print(f'   True anomaly (deg)            = {coe[5] / deg}')
    print(f'   Semimajor axis (km)           = {coe[6]}')
    print(f'   Periapse radius (km)          = {coe[0]**2 / mu / (1 + coe[1])}')

    if coe[1] < 1:
        T = 2 * np.pi / np.sqrt(mu) * coe[6]**1.5
        print(f'   Period:')
        print(f'     Seconds                     = {T}')
        print(f'     Minutes                     = {T / 60}')
        print(f'     Hours                       = {T / 3600}')
        print(f'     Days                        = {T / 24 / 3600}')

    print('-----------------------------------------------------\n')

if __name__ == '__main__':
    gauss_test()