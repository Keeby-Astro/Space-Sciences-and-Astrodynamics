# ALGORITHM 4.6: CALCULATE THE GROUND TRACK OF A SATELLITE FROM ITS ORBITAL ELEMENTS

import numpy as np
import matplotlib.pyplot as plt
import sv_from_coe
import kepler_E
import ra_and_dec_from_r

def ground_track():
    '''
    This program plots the ground track of an earth satellite
    for which the orbital elements are specified

    mu          - gravitational parameter (km^3/s^2)
    deg         - factor that converts degrees to radians
    J2          - second zonal harmonic
    Re          - earth s radius (km)
    we          - earth s angular velocity (rad/s)
    rP          - perigee of orbit (km)
    rA          - apogee of orbit (km)
    TA, TAo     - true anomaly, initial true anomaly of satellite (rad)
    RA, RAo     - right ascension, initial right ascension of the node (rad)
    incl        - orbit inclination (rad)
    wp, wpo     - argument of perigee, initial argument of perigee (rad)
    n_periods   - number of periods for which ground track is to be plotted
    a           - semimajor axis of orbit (km)
    T           - period of orbit (s)
    e           - eccentricity of orbit
    h           - angular momentum of orbit (km^2/s)
    E, Eo       - eccentric anomaly, initial eccentric anomaly (rad)
    M, Mo       - mean anomaly, initial mean anomaly (rad)
    to, tf      - initial and final times for the ground track (s)
    fac         - common factor in Equations 4.53 and 4.53
    RAdot       - rate of regression of the node (rad/s)
    wpdot       - rate of advance of perigee (rad/s)
    times       - times at which ground track is plotted (s)
    ra          - vector of right ascensions of the spacecraft (deg)
    dec         - vector of declinations of the spacecraft (deg)
    TA          - true anomaly (rad)
    r           - perifocal position vector of satellite (km)
    R           - geocentric equatorial position vector (km)
    R1          - DCM for rotation about z through RA
    R2          - DCM for rotation about x through incl
    R3          - DCM for rotation about z through wp
    QxX         - DCM for rotation from perifocal to geocentric equatorial
    Q           - DCM for rotation from geocentric equatorial
                  into earth-fixed frame
    r_rel       - position vector in earth-fixed frame (km)
    alpha       - satellite right ascension (deg)
    delta       - satellite declination (deg)
    n_curves    - number of curves comprising the ground track plot
    RA          - cell array containing the right ascensions for each of
                  the curves comprising the ground track plot
    Dec         - cell array containing the declinations for each of
                  the curves comprising the ground track plot

    User py-functions required: sv_from_coe, kepler_E, ra_and_dec_from_r
    '''
    # Constants
    deg = np.pi / 180
    mu = 398600.4418
    J2 = 0.00108263
    Re = 6378.14
    we = (2 * np.pi + 2 * np.pi / 365.26) / (24 * 3600)

    # Data declaration for Example 4.12:
    rP = 6700
    rA = 10000
    TAo = 230 * deg
    Wo = 270 * deg
    incl = 60 * deg
    wpo = 45 * deg
    n_periods = 3.25
    # End data declaration

    # Compute the initial time (since perigee) and
    # the rates of node regression and perigee advance
    a = (rA + rP) / 2
    T = 2 * np.pi / np.sqrt(mu) * a**(3 / 2)
    e = (rA - rP) / (rA + rP)
    h = np.sqrt(mu * a * (1 - e**2))
    Eo = 2 * np.arctan(np.tan(TAo / 2) * np.sqrt((1 - e) / (1 + e)))
    Mo = Eo - e * np.sin(Eo)
    to = Mo * (T / (2 * np.pi))
    tf = to + n_periods * T
    fac = -3 / 2 * np.sqrt(mu) * J2 * Re**2 / (1 - e**2)**2 / a**(7 / 2)
    Wdot = fac * np.cos(incl)
    wpdot = fac * (5 / 2 * np.sin(incl)**2 - 2)

    ra, dec = find_ra_and_dec(to, tf, T, h, e, mu, Wo, Wdot, incl, wpo, wpdot)
    RA, Dec, n_curves = form_separate_curves(ra, dec)
    plot_ground_track(RA, Dec, ra, dec)
    print_orbital_data(h, e, Wo, incl, wpo, TAo, T, to, rP, rA, Wdot, wpdot, mu, a)

# –––––––––––––––––––––––––––
def find_ra_and_dec(to, tf, T, h, e, mu, Wo, Wdot, incl, wpo, wpdot):
    '''
    Propagates the orbit over the specified time interval, transforming
    the position vector into the earth-fixed frame and computing
    the right ascension and declination histories.
    '''
    we = (2 * np.pi + 2 * np.pi / 365.26) / (24 * 3600)

    times = np.linspace(to, tf, 1000)
    ra = []
    dec = []
    theta = 0
    for t in times:
        M = 2 * np.pi / T * t
        E = kepler_E.kepler_E(e, M)
        TA = 2 * np.arctan(np.tan(E / 2) * np.sqrt((1 + e) / (1 - e)))
        r = h**2 / mu / (1 + e * np.cos(TA)) * np.array([np.cos(TA), np.sin(TA), 0])

        W = Wo + Wdot * t
        wp = wpo + wpdot * t

        R1 = np.array([[np.cos(W), np.sin(W), 0],
                       [-np.sin(W), np.cos(W), 0],
                       [0, 0, 1]])

        R2 = np.array([[1, 0, 0],
                       [0, np.cos(incl), np.sin(incl)],
                       [0, -np.sin(incl), np.cos(incl)]])

        R3 = np.array([[np.cos(wp), np.sin(wp), 0],
                       [-np.sin(wp), np.cos(wp), 0],
                       [0, 0, 1]])

        QxX = (R3 @ R2 @ R1).T
        R = QxX @ r

        theta = we * (t - to)
        Q = np.array([[np.cos(theta), np.sin(theta), 0],
                      [-np.sin(theta), np.cos(theta), 0],
                      [0, 0, 1]])
        r_rel = Q @ R

        alpha, delta = ra_and_dec_from_r.ra_and_dec_from_r(r_rel)
        ra.append(alpha)
        dec.append(delta)
    
    return np.array(ra), np.array(dec)

# –––––––––––––––––––––––––––
def form_separate_curves(ra, dec):
    '''
    Breaks the ground track up into separate curves starting and
    terminating at right ascensions in the range [0, 360 deg].
    '''
    tol = 100
    n_curves = 1
    curve_no = 1
    k = 0
    ra_prev = ra[0]
    RA = {}
    Dec = {}

    for i in range(len(ra)):
        if abs(ra[i] - ra_prev) > tol:
            curve_no += 1
            n_curves += 1
            k = 0
        if curve_no not in RA:
            RA[curve_no] = []
            Dec[curve_no] = []
        RA[curve_no].append(ra[i])
        Dec[curve_no].append(dec[i])
        ra_prev = ra[i]
    
    return RA, Dec, n_curves

# –––––––––––––––––––––––––––
def plot_ground_track(RA, Dec, ra, dec):
    plt.figure()
    plt.xlabel('East Longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    plt.axis('equal')
    plt.grid(True)

    for curve_no in RA:
        plt.plot(RA[curve_no], Dec[curve_no])

    plt.axis([0, 360, -90, 90])
    plt.text(ra[0], dec[0], 'o Start')
    plt.text(ra[-1], dec[-1], 'o Finish')
    plt.axhline(0, color='k')  # the equator
    plt.show()

# –––––––––––––––––––––––––
def print_orbital_data(h, e, Wo, incl, wpo, TAo, T, to, rP, rA, Wdot, wpdot, mu, a):
    '''
    Prints the orbital data of the satellite.
    '''
    coe = [h, e, Wo, incl, wpo, TAo]
    ro, vo = sv_from_coe.sv_from_coe(coe, mu)
    deg = np.pi / 180

    print('\n-----------------------------------------------\n')
    print(f' Angular momentum     = {h} km^2/s')
    print(f' Eccentricity         = {e}')
    print(f' Semimajor axis       = {a} km')
    print(f' Perigee radius       = {rP} km')
    print(f' Apogee radius        = {rA} km')
    print(f' Period               = {T/3600} hours')
    print(f' Inclination          = {incl/deg} deg')
    print(f' Initial true anomaly = {TAo/deg} deg')
    print(f' Time since perigee   = {to/3600} hours')
    print(f' Initial RA           = {Wo/deg} deg')
    print(f' RA_dot               = {Wdot/deg*T} deg/period')
    print(f' Initial wp           = {wpo/deg} deg')
    print(f' wp_dot               = {wpdot/deg*T} deg/period')
    print('\n r0 = [{:12g}, {:12g}, {:12g}] (km)'.format(*ro))
    print(f' magnitude = {np.linalg.norm(ro)} km\n')
    print('\n v0 = [{:12g}, {:12g}, {:12g}] (km/s)'.format(*vo))
    print(f' magnitude = {np.linalg.norm(vo)} km/s\n')
    print('-----------------------------------------------\n')

if __name__ == '__main__':
    ground_track()