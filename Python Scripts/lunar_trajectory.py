import sys
import os
import importlib.util
import datetime

std_lib_dir = os.path.dirname(os.__file__)  
std_bisect_path = os.path.join(std_lib_dir, 'bisect.py')
spec = importlib.util.spec_from_file_location("bisect", std_bisect_path)
real_bisect = importlib.util.module_from_spec(spec)
spec.loader.exec_module(real_bisect)
sys.modules["bisect"] = real_bisect

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import simpsons_lunar_ephemeris
import ra_and_dec_from_r

def lunar_trajectory():
    '''
    This program presents the graphical solution of the motion of a
    spacecraft in the gravity fields of both the earth and the moon for
    the initial data provided in the input declaration below.
    MATLAB's ode45 Runge-Kutta solver is used.

    deg                         - conversion factor, degrees to radians
    days                        - conversion factor, days to seconds
    Re, Rm                      - radii of earth and moon, respectively (km)
    m_e, m_m                    - masses of earth and moon, respectively (kg)
    mu_e, mu_m                  - gravitational parameters of earth and moon,
                                    respectively (km^3/s^2)
    D                           - semimajor axis of moon's orbit (km)
    I_, J_, K_                  - unit vectors of ECI frame
    RS                          - radius of moon's sphere of influence (km)
                                    year, month, hour, minute
    second                      - Date and time of spacecraft's lunar arrival
    t0                          - initial time on the trajectory (s)
    z0                          - initial altitude of the trajectory (km)
    alpha0, dec0                - initial right ascension and declination of
                                    spacecraft (deg)
    gamma0                      - initial flight path angle (deg)
    fac                         - ratio of spaccraft's initial speed to the
                                    escape speed.
    ttt                         - predicted time to perilune (s)
    tf                          - time at end of trajectory (s)
    jd0                         - julian date of lunar arrival
    rm0, vm0                    - state vector of the moon at jd0 (km, km/s)
    RA, Dec                     - right ascension and declination of the moon
                                    at jd0 (deg)
    hmoon_, hmoon               - moon's angular momentum vector and magnitude
                                    at jd0 (km^2/s)
    inclmoon                    - inclination of moon's orbit earth's
                                    equatorial plane (deg)
    r0                          - initial radius from earth's center to
                                    probe (km)
    r0_                         - initial ECI position vector of probe (km)
    vesc                        - escape speed at r0 (km/s)
    v0                          - initial ECI speed of probe (km/s)
    w0_                         - unit vector normal to plane of translunar
                                    orbit at time t0
    ur_                         - radial unit vector to probe at time t0
    uperp_                      - transverse unit vector at time t0
    vr                          - initial radial speed of probe (km/s)
    vperp                       - initial transverse speed of probe (km/s)
    v0_                         - initial velocity vector of probe (km/s)
    uv0_                        - initial tangential unit vector
    y0                          - initial state vector of the probe (km, km/s)
    t                           - vector containing the times from t0 to tf at
                                    which the state vector is evaluated (s)
    y                           - a matrix whose 6 columns contain the inertial
                                    position and velocity components evaluated
                                    at the times t(:) (km, km/s)
    X, Y, Z                     - the probe's inertial position vector history
    vX, vY, VZ                  - the probe's inertial velocity history
    x, y, z                     - the probe's position vector history in the
                                    Moon-fixed frame
    Xm, Ym, Zm                  - the Moon's inertial position vector history
    vXm, vYm, vZm               - the Moon's inertial velocity vector history
    ti                          - the ith time of the set [t0,tf] (s)
    r_                          - probe's inertial position vector at time ti
                                    (km)
    r                           - magnitude of r_ (km)
    jd                          - julian date of corresponding to ti (days)
    rm_, vm_                    - the moon's state vector at time ti (km,km/s)
    x_, y_, z_                  - vectors along the axes of the rotating
                                    moon-fixed at time ti (km)
    i_, j_, k_                  - unit vectors of the moon-fixed rotating frame
                                    at time ti
    Q                           - DCM of transformation from ECI to moon-fixed
                                    frame at time ti
    rx_                         - probe's inertial position vector in moon-
                                    fixed coordinates at time ti (km)
    rmx_                        - Moon's inertial position vector in moon-
                                    fixed coordinates at time ti (km)
    dist_                       - position vector of probe relative to the moon
                                    at time ti (km)
    dist                        - magnitude of dist_ (km)
    dist_min                    - perilune of trajectory (km)
    rmTLI_                      - Moon's position vector at TLI
    RATLI, DecTLI               - Moon's right ascension and declination at
                                    TKI (deg)
    v_atdmin_                   - Probe's velocity vector at perilune (km/s)
    rm_perilume, vm_perilune    - Moon's state vector when the probe is at
                                    perilune (km, km/s)
    rel_speed                   - Speed of probe relative to the Moon at
                                    perilune (km/s)
    RA_at_perilune              - Moon's RA at perilune arrival (deg)
    Dec_at_perilune             - Moon's Dec at perilune arrival (deg)
    target_error                - Distance between Moon's actual position at
                                    perilune arrival and its position after the
                                    predicted flight time, ttt (km).
    rms_                        - position vector of moon relative to
                                    spacecraft (km)
    rms                         - magnitude of rms_ (km)
    aearth_                     - acceleration of spacecraft due to
                                    earth (km/s^2)
    amoon_                      - acceleration of spacecraft due to
                                    moon (km/s^2)
    atot_                       - aearth_ + amoon_ (km/s^2)
    binormal_                   - unit vector normal to the osculating plane
    incl                        - angle between inertial Z axis and the
                                    binormal (deg)
    rend_                       - Position vector of end point of trajectory
                                    (km)
    alt_end                     - Altitude of end point of trajectory (km)
    ra_end, dec_end             - Right ascension and declination of end point
                                    of trajectory (km)

    User M-functions required: none
    User subfunctions required: rates, plotit_XYZ, plotit_xyz
    '''
    #...general data
    deg  = np.pi/180.0
    days = 24.0*3600.0
    Re   = 6378.0
    Rm   = 1737.0
    m_e  = 5974.e21
    m_m  = 73.48e21
    mu_e = 398600.4
    mu_m = 4902.8
    D    = 384400.0
    RS   = D*((m_m/m_e)**(2.0/5.0))
    #...

    #...Data declaration
    Title = 'Lunar Trajectory'
    #   Date and time of lunar arrival:
    year    = 2020
    month   = 5
    day     = 4
    hour    = 12
    minute  = 0
    second  = 0
    t0      = 0.0
    z0      = 320.0
    alpha0  = 90.0
    dec0    = 15.0
    gamma0  = 40.0
    fac     = 0.9924  # Fraction of Vesc
    ttt     = 3.0*days
    tf      = ttt + 2.667*days
    #...End data declaration

    #...State vector of moon at target date:
    def juliandate_py(yr, mo, d, hr, mn, sec):
        date_obj = datetime.datetime(yr, mo, d, hr, mn, sec)
        ordinal = date_obj.toordinal()
        frac_day = (hr + mn/60.0 + sec/3600.0)/24.0
        jd_approx = ordinal + 1721424.5 + (date_obj.hour +
                                           date_obj.minute/60.0 +
                                           date_obj.second/3600.0)/24.0
        return jd_approx
    
    jd0 = juliandate_py(year, month, day, hour, minute, second)

    rm0_, vm0_ = simpsons_lunar_ephemeris.simpsons_lunar_ephemeris(jd0)
    RA, Dec    = ra_and_dec_from_r.ra_and_dec_from_r(rm0_)
    distance   = np.linalg.norm(rm0_)
    hmoon_     = np.cross(rm0_, vm0_)
    hmoon      = np.linalg.norm(hmoon_)
    inclmoon   = np.degrees(np.arccos(hmoon_[2]/hmoon))

    #...Initial position vector of probe:
    I_ = np.array([1.0, 0.0, 0.0])
    J_ = np.array([0.0, 1.0, 0.0])
    K_ = np.cross(I_, J_)
    r0 = Re + z0
    # In radians:
    alpha0_rad = alpha0*deg
    dec0_rad   = dec0*deg

    r0_ = r0*(np.cos(alpha0_rad)*np.cos(dec0_rad)*I_ +
              np.sin(alpha0_rad)*np.cos(dec0_rad)*J_ +
              np.sin(dec0_rad)*K_)
    vesc = np.sqrt(2.0*mu_e/r0)
    v0   = fac*vesc
    w0_  = np.cross(r0_, rm0_)
    w0_  = w0_/np.linalg.norm(w0_)

    #...Initial velocity vector of probe:
    ur_    = r0_/np.linalg.norm(r0_)
    uperp_ = np.cross(w0_, ur_)
    uperp_ = uperp_/np.linalg.norm(uperp_)

    gamma0_rad = gamma0*deg
    vr    = v0*np.sin(gamma0_rad)
    vperp = v0*np.cos(gamma0_rad)
    v0_   = vr*ur_ + vperp*uperp_
    uv0_  = v0_/v0

    #...Initial state vector of the probe:
    y0 = np.array([r0_[0], r0_[1], r0_[2],
                   v0_[0], v0_[1], v0_[2]])

    def rates(t, y, jd0, days, mu_m, mu_e, ttt):
        '''
        This function evaluates the 3D acceleration of the spacecraft in a
        restricted 3-body system at time t from its position and velocity
        and the position of the moon at that time.

        t           - time (s)
        ttt         - flight time, TLI to target point (s)
        jd0         - Julian Date on arrival at target (days)
        jd          - Julian Date at time t (days)
        X, Y, Z     - Components of spacecraft's geocentric position vector (km)
        vX, vY, vZ  - Components of spacecraft's geocentric velocity vector (km/s)
        aX, aY, aZ  - Components of spacecraft's geocentric acceleration
                        vector (km/s^2)
        y           - column vector containing the geocentric position and
                        velocity components of the spacecraft at time t
        r_          - geocentric position vector [X Y Z] of the spacecraft
        rm_         - geocentric position vector of the moon
        rms_        - rm_ - r_, the position of the moon relative to the
                        spacecraft
        aearth_     - spacecraft acceleration vector due to earth's gravity
        amoon_      - spacecraft acceleration vector due to lunar gravity
        a_          - total spacecraft acceleration vector
        dydt        - column vector containing the geocentric velocity and
                        acceleration components of the spacecraft at time t
        '''
        X = y[0]
        Y = y[1]
        Z = y[2]
        vX = y[3]
        vY = y[4]
        vZ = y[5]

        r_ = np.array([X, Y, Z])
        r = np.linalg.norm(r_)

        # Time-to-JD conversion:
        jd = jd0 - (ttt - t)/days

        # Moon ephemeris at current jd:
        rm_, vm_ = simpsons_lunar_ephemeris.simpsons_lunar_ephemeris(jd)
        rm = np.linalg.norm(rm_)

        rms_ = rm_ - r_
        rms = np.linalg.norm(rms_)

        aearth_ = -mu_e*r_/r**3
        amoon_  = mu_m*(rms_/rms**3 - rm_/rm**3)
        a_      = aearth_ + amoon_

        aX = a_[0]
        aY = a_[1]
        aZ = a_[2]

        return [vX, vY, vZ, aX, aY, aZ]

    sol = solve_ivp(fun=lambda t, y: rates(t, y, jd0, days, mu_m, mu_e, ttt),
                    t_span=(t0, tf),
                    y0=y0,
                    method='RK45',
                    rtol=1.e-10,
                    atol=1.e-10,
                    dense_output=True)

    # Extract solution times and states
    t = sol.t
    y_sol = sol.y.T  # shape: (len(t), 6)

    #...Spacecraft trajectory in ECI frame
    X = y_sol[:, 0]
    Y_ = y_sol[:, 1]
    Z_ = y_sol[:, 2]
    vX = y_sol[:, 3]
    vY = y_sol[:, 4]
    vZ = y_sol[:, 5]

    #...Prepare arrays for the rotating Moon-fixed frame
    x_tra = []
    y_tra = []
    z_tra = []

    #...Moon trajectory in ECI frame
    Xm = []
    Ym = []
    Zm = []
    vXm = []
    vYm = []
    vZm = []
    xm_tra = []
    ym_tra = []
    zm_tra = []

    #...Compute the Moon's trajectory from an ephemeris, find perilune of the
    #   probe's trajectory, and project the probe's trajectory onto the axes
    #   of the Moon-fixed rotating frame:
    dist_min = 1.e30
    imin = 0
    ist_min_ = None

    for i in range(len(t)):
        ti = t[i]
        r_ = np.array([X[i], Y_[i], Z_[i]])

        # Time -> JD
        jd = jd0 - (ttt - ti)/days

        rm_i, vm_i = simpsons_lunar_ephemeris.simpsons_lunar_ephemeris(jd)
        Xm.append(rm_i[0])
        Ym.append(rm_i[1])
        Zm.append(rm_i[2])
        vXm.append(vm_i[0])
        vYm.append(vm_i[1])
        vZm.append(vm_i[2])

        # Construct moon-fixed axes:
        x_axis = rm_i
        z_axis = np.cross(x_axis, vm_i)
        y_axis = np.cross(z_axis, x_axis)
        i_unit = x_axis/np.linalg.norm(x_axis)
        j_unit = y_axis/np.linalg.norm(y_axis)
        k_unit = z_axis/np.linalg.norm(z_axis)

        # DCM from ECI to moon-fixed
        Q = np.array([i_unit, j_unit, k_unit])

        # Position of probe in moon-fixed:
        rx_ = Q.dot(r_)
        x_tra.append(rx_[0])
        y_tra.append(rx_[1])
        z_tra.append(rx_[2])

        # Position of moon in moon-fixed:
        rmx_ = Q.dot(rm_i)
        xm_tra.append(rmx_[0])
        ym_tra.append(rmx_[1])
        zm_tra.append(rmx_[2])

        # Find perilune
        dist_ = r_ - rm_i
        dist  = np.linalg.norm(dist_)
        if dist < dist_min:
            imin = i
            ist_min_ = dist_
            dist_min = dist

    #...Location of the Moon at TLI:
    rmTLI_ = np.array([Xm[0], Ym[0], Zm[0]])
    RATLI, DecTLI = ra_and_dec_from_r.ra_and_dec_from_r(rmTLI_)

    #...Spacecraft velocity at perilune:
    v_atdmin_ = np.array([vX[imin], vY[imin], vZ[imin]])

    #...State vector and celestial position of moon when probe is at perilune:
    rm_perilune_ = np.array([Xm[imin], Ym[imin], Zm[imin]])
    vm_perilune_ = np.array([vXm[imin], vYm[imin], vZm[imin]])
    RA_at_perilune, Dec_at_perilune = ra_and_dec_from_r.ra_and_dec_from_r(rm_perilune_)
    target_error = np.linalg.norm(rm_perilune_ - rm0_)

    #...Speed of probe relative to Moon at perilune:
    rel_speed = np.linalg.norm(v_atdmin_ - vm_perilune_)

    #...End point of trajectory:
    rend_   = np.array([X[-1], Y_[-1], Z_[-1]])
    alt_end = np.linalg.norm(rend_) - Re
    ra_end, dec_end = ra_and_dec_from_r.ra_and_dec_from_r(rend_)

    #...Find the history of the trajectory's binormal:
    time_hist = []
    rms = []
    incl = []
    for i in range(imin+1): 
        time_hist.append(t[i])
        r__   = np.array([X[i], Y_[i], Z_[i]])
        r_mag = np.linalg.norm(r__)
        v__   = np.array([vX[i], vY[i], vZ[i]])
        rm__  = np.array([Xm[i], Ym[i], Zm[i]])
        rm_mag = np.linalg.norm(rm__)

        rms_ = rm__ - r__
        rms_val = np.linalg.norm(rms_)
        rms.append(rms_val)

        aearth_ = -mu_e*r__/(r_mag**3)
        amoon_  = mu_m*(rms_/(rms_val**3) - rm__/(rm_mag**3))
        atot_   = aearth_ + amoon_

        binormal_ = np.cross(v__, atot_)
        binormal_ = binormal_/np.linalg.norm(binormal_)
        binormalz = binormal_[2]
        incl.append(np.degrees(np.arccos(binormalz)))

    #...Output:
    print('\n\n' + Title + '\n')
    print('Date and time of arrival at moon: ', end='')
    print(str(month)+'/'+str(day)+'/'+str(year)+' '+str(hour)+':'+str(minute)+':'+str(second))
    print("Moon's position: ")
    print(' Distance                        = {:11g} km'.format(distance))
    print(' Right Ascension                 = {:11g} deg'.format(RA))
    print(' Declination                     = {:11g} deg '.format(Dec))
    print(" Moon's orbital inclination      = {:11g} deg".format(inclmoon))

    print('\nThe probe at earth departure (t = {} sec):'.format(t0))
    print(' Altitude                        = {:11g} km'.format(z0))
    print(' Right ascension                 = {:11g} deg'.format(alpha0))
    print(' Declination                     = {:11g} deg'.format(dec0))
    print(' Flight path angle               = {:11g} deg'.format(gamma0))
    print(' Speed                           = {:11g} km/s'.format(v0))
    print(' Escape speed                    = {:11g} km/s'.format(vesc))
    print(' v/vesc                          = {:11g}'.format(v0/vesc))

    inc_tlo = np.degrees(np.arccos(w0_[2]))
    print(' Inclination of translunar orbit = {:11g} deg'.format(inc_tlo))

    print('\nThe moon when the probe is at TLI:')
    print(' Distance                        = {:11g} km'.format(np.linalg.norm(rmTLI_)))
    print(' Right ascension                 = {:11g} deg'.format(RATLI))
    print(' Declination                     = {:11g} deg'.format(DecTLI))

    print('\nThe moon when the probe is at perilune: ')
    print(' Distance                        = {:11g} km'.format(np.linalg.norm(rm_perilune_)))
    print(' Speed                           = {:11g} km/s'.format(np.linalg.norm(vm_perilune_)))
    print(' Right ascension                 = {:11g} deg'.format(RA_at_perilune))
    print(' Declination                     = {:11g} deg'.format(Dec_at_perilune))
    print(' Target error                    = {:11g} km'.format(target_error))

    print('\n\nThe probe at perilune:')
    print(' Altitude                        = {:11g} km'.format(dist_min - Rm))
    print(' Speed                           = {:11g} km/s'.format(np.linalg.norm(v_atdmin_)))
    print(' Relative speed                  = {:11g} km/s'.format(rel_speed))
    print(' Inclination of osculating plane = {:11g} deg'.format(incl[imin]))
    print(' Time from TLI to perilune       = {:11g} hours ({:g} days)'.format(
          abs(t[imin])/3600.0, abs(t[imin])/3600.0/24.0))

    print('\n\nTotal time of flight             = {:11g} days'.format(t[-1]/days))
    print('Time to target point             = {:11g} days'.format(ttt/days))
    print('Final earth altitude             = {:11g} km'.format(alt_end))
    print('Final right ascension            = {:11g} deg'.format(ra_end))
    print('Final declination                = {:11g} deg'.format(dec_end))

    #...Graphical output:
    #   Plot the trajectory relative to the inertial frame:
    def plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin):
        '''
        % --------------------------------------
        % function plotit_XYZ(X,Y,Z,Xm,Ym,Zm,imin)
        % --------------------------------------
        % global Re Rm
        %
        % figure ('Name','Trajectories of Spacecraft (red) and Moon (green)', ...
        %         'Color', [1 1 1])
        % [xx, yy, zz] = sphere(128);
        % hold on
        %
        % ... etc ...
        '''
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        fig.suptitle('Trajectories of Spacecraft (red) and Moon (green)')

        # Earth radius sphere:
        # We mimic [xx, yy, zz] = sphere(128) using numpy:
        u, v = np.mgrid[0:2*np.pi:128j, 0:np.pi:128j]
        xx = np.cos(u)*np.sin(v)
        yy = np.sin(u)*np.sin(v)
        zz = np.cos(v)

        # Axes lines:
        L = 20*Re
        ax.plot([0, L], [0, 0], [0, 0], color='k')
        ax.text(L, 0, 0, 'X', fontsize=12)
        ax.plot([0, 0], [0, L], [0, 0], color='k')
        ax.text(0, L, 0, 'Y', fontsize=12)
        ax.plot([0, 0], [0, 0], [0, L], color='k')
        ax.text(0, 0, L, 'Z', fontsize=12)

        # Earth:
        ax.plot_surface(Re*xx, Re*yy, Re*zz, color='b', alpha=0.5, linewidth=0)

        # Spacecraft at TLI
        ax.scatter(X[0], Y[0], Z[0], color='k', s=10, marker='o')
        # Spacecraft at closest approach
        ax.scatter(X[imin], Y[imin], Z[imin], color='k', s=8, marker='o')
        # Spacecraft at tf
        ax.scatter(X[-1], Y[-1], Z[-1], color='r', s=10, marker='o')

        # Moon at TLI
        ax.text(Xm[0], Ym[0], Zm[0], 'Moon at TLI')
        ax.plot_surface(Rm*xx + Xm[0], Rm*yy + Ym[0], Rm*zz + Zm[0],
                        color='g', alpha=0.99, linewidth=0)

        # Moon at closest approach
        ax.plot_surface(Rm*xx + Xm[imin], Rm*yy + Ym[imin], Rm*zz + Zm[imin],
                        color='g', alpha=0.99, linewidth=0)

        # Moon at end
        ax.plot_surface(Rm*xx + Xm[-1], Rm*yy + Ym[-1], Rm*zz + Zm[-1],
                        color='g', alpha=0.99, linewidth=0)

        # Spacecraft trajectory
        ax.plot(X, Y, Z, color='r', linewidth=1.5)
        # Moon trajectory
        ax.plot(Xm, Ym, Zm, color='g', linewidth=0.5)

        # Make axes equal, turn off
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_box_aspect((1,1,1))  # equal aspect
        plt.show()

    def plotit_xyz(xvals, yvals, zvals, xmvals, ymvals, zmvals, imin):
        '''
        % --------------------------------------
        % function plotit_xyz(x,y,z,xm,ym,zm,imin)
        % --------------------------------------
        % global Re Rm Rm0_ Q0
        % ... etc ...
        '''
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        fig.suptitle('Spacecraft trajectory in Moon-fixed rotating frame')

        # mimic [xx, yy, zz] = sphere(128)
        u, v = np.mgrid[0:2*np.pi:128j, 0:np.pi:128j]
        xx = np.cos(u)*np.sin(v)
        yy = np.sin(u)*np.sin(v)
        zz = np.cos(v)

        # Spacecraft trajectory
        ax.plot(xvals, yvals, zvals, color='r', linewidth=2.0)

        # Moon trajectory
        ax.plot(xmvals, ymvals, zmvals, color='g', linewidth=0.5)

        # Earth (roughly, though in the rotating frame it might not be at origin)
        ax.plot_surface(Re*xx, Re*yy, Re*zz, color='b', alpha=0.5, linewidth=0)

        # Axes lines in geocentric moon-fixed
        L1 = 63*Re
        L2 = 20*Re
        L3 = 29*Re
        ax.plot([0, L1], [0, 0], [0, 0], color='k')
        ax.text(L1, 0, 0, 'x', fontsize=12)
        ax.plot([0, 0], [0, L2], [0, 0], color='k')
        ax.text(0, L2, 0, 'y', fontsize=12)
        ax.plot([0, 0], [0, 0], [0, L3], color='k')
        ax.text(0, 0, L3, 'z', fontsize=12)

        # Spacecraft at TLI
        ax.scatter(xvals[0], yvals[0], zvals[0], color='k', s=10, marker='o')
        # Spacecraft at closest approach
        ax.scatter(xvals[imin], yvals[imin], zvals[imin], color='k', s=8, marker='o')
        # Spacecraft at tf
        ax.scatter(xvals[-1], yvals[-1], zvals[-1], color='r', s=10, marker='o')

        # Moon at TLI
        ax.text(xmvals[0], ymvals[0], zmvals[0], 'Moon at TLI')
        ax.plot_surface(Rm*xx + xmvals[0],
                        Rm*yy + ymvals[0],
                        Rm*zz + zmvals[0],
                        color='g', alpha=0.99, linewidth=0)

        # Moon at spacecraft closest approach
        ax.plot_surface(Rm*xx + xmvals[imin],
                        Rm*yy + ymvals[imin],
                        Rm*zz + zmvals[imin],
                        color='g', alpha=0.99, linewidth=0)

        # Moon at end
        ax.plot_surface(Rm*xx + xmvals[-1],
                        Rm*yy + ymvals[-1],
                        Rm*zz + zmvals[-1],
                        color='g', alpha=0.99, linewidth=0)

        ax.set_box_aspect((1,1,1))
        plt.show()

    #   Plot the trajectory relative to the inertial frame:
    plotit_XYZ(X, Y_, Z_, Xm, Ym, Zm, imin)

    #   Plot inclination of the osculating plane vs distance from the Moon
    fig_incl = plt.figure()
    ax_incl = fig_incl.add_subplot(111)
    ax_incl.plot(np.array(rms)/RS, incl)
    ax_incl.axhline(y=90, color='r', linestyle='-')
    ax_incl.set_title('Osculating Plane Inclination vs Distance from Moon')
    ax_incl.set_xlabel('$r_{ms}$/$R_s$')
    ax_incl.set_ylabel('Inclination (deg)')
    ax_incl.grid(True, which='both')
    plt.show()

    #   Plot the trajectory relative to the rotating Moon-fixed frame:
    plotit_xyz(x_tra, y_tra, z_tra, xm_tra, ym_tra, zm_tra, imin)

    #...End graphical output
    return

if __name__ == '__main__':
    lunar_trajectory()