import numpy as np
import rkf45
import coe_from_sv
import rv_from_r0v0_ta

def integrate_thrust():
    '''
    This function uses rkf45 to numerically integrate Equation 6.26 during
    the delta-v burn and then find the apogee of the post-burn orbit.

    mu      - gravitational parameter (km^3/s^2)
    RE      - earth radius (km)
    g0      - sea level acceleration of gravity (m/s^2)
    T       - rated thrust of rocket engine (kN)
    Isp     - specific impulse of rocket engine (s)
    m0      - initial spacecraft mass (kg)
    r0      - initial position vector (km)
    v0      - initial velocity vector (km/s)
    t0      - initial time (s)
    t_burn  - rocket motor burn time (s)
    y0      - column vector containing r0, v0 and m0
    t       - column vector of the times at which the solution is found (s)
    y       - a matrix whose elements are:
                columns 1, 2 and 3:
                    The solution for the x, y and z components of the
                    position vector r at the times t
                columns 4, 5 and 6:
                    The solution for the x, y and z components of the
                    velocity vector v at the times t
                column 7:
                    The spacecraft mass m at the times t
    r1      - position vector after the burn (km)
    v1      - velocity vector after the burn (km/s)
    m1      - mass after the burn (kg)
    coe     - orbital elements of the post-burn trajectory
            (h e RA incl w TA a)
    ra      - position vector vector at apogee (km)
    va      - velocity vector at apogee (km)
    rmax    - apogee radius (km)

    User py-functions required: rkf45, coe_from_sv, rv_from_r0v0_ta
    User py-subfunctions required: rates, output
    '''
    #...Constants
    deg = np.pi / 180
    mu = 398600.1148 # gravitational parameter (km^3/s^2)
    RE = 6378.14     # earth radius (km)
    g0 = 9.807       # sea level acceleration of gravity (m/s^2)

    #...Input data
    r0 = np.array([RE + 480, 0, 0])  # initial position vector (km)
    v0 = np.array([0, 7.7102, 0])    # initial velocity vector (km/s)
    t0 = 0                           # initial time (s)
    t_burn = 261.1127                # burn time (s)
    m0 = 2000                        # initial mass (kg)
    T = 10                           # thrust (kN)
    Isp = 300                        # specific impulse (s)

    # Initial state vector
    y0 = np.hstack((r0, v0, m0))
    t, y = rkf45.rkf45(rates, [t0, t_burn], y0, 1.e-16)

    # State vector and mass after the burn
    r1 = y[-1, 0:3]
    v1 = y[-1, 3:6]
    m1 = y[-1, 6]
    coe = coe_from_sv.coe_from_sv(r1, v1, mu)
    e = coe[1]    # eccentricity
    TA = coe[5]   # true anomaly (radians)
    a = coe[6]    # semimajor axis (km)

    # Find state vector at apogee of the post-burn trajectory
    if TA <= np.pi:
        dtheta = np.pi - TA
    else:
        dtheta = 3 * np.pi - TA
    ra, va = rv_from_r0v0_ta.rv_from_r0v0_ta(r1, v1, np.degrees(dtheta), mu)
    rmax = np.linalg.norm(ra)

    # Output results
    output(r0, v0, m0, T, t_burn, m1, r1, v1, e, a, ra, va)

# ––––––––––––––––––––––––
def rates(t, f):
    '''
    This function calculates the acceleration vector using Equation 6.26.

    t           - time (s)
    f           - column vector containing the position vector, velocity
                  vector and the mass at time t
    x, y, z     - components of the position vector (km)
    vx, vy, vz  - components of the velocity vector (km/s)
    m           - mass (kg)
    r           - magnitude of the the position vector (km)
    v           - magnitude of the velocity vector (km/s)
    ax, ay, az  - components of the acceleration vector (km/s^2)
    mdot        - rate of change of mass (kg/s)
    dfdt        - column vector containing the velocity and acceleration
                  components and the mass rate
    '''
    mu = 398600  # gravitational parameter (km^3/s^2)
    g0 = 9.807   # sea level acceleration of gravity (m/s^2)
    Isp = 300    # specific impulse (s)
    T = 10       # thrust (kN)
    
    x, y, z = f[0], f[1], f[2]
    vx, vy, vz = f[3], f[4], f[5]
    m = f[6]

    r = np.linalg.norm([x, y, z])
    v = np.linalg.norm([vx, vy, vz])
    
    # Acceleration components
    ax = -mu * x / r**3 + T / m * vx / v
    ay = -mu * y / r**3 + T / m * vy / v
    az = -mu * z / r**3 + T / m * vz / v
    mdot = -T * 1000 / g0 / Isp  # Rate of change of mass (kg/s)

    return np.array([vx, vy, vz, ax, ay, az, mdot])

def output(r0, v0, m0, T, t_burn, m1, r1, v1, e, a, ra, va):
    print('\n\n------------------------------------------------------')
    print('\nBefore ignition:')
    print(f'  Mass = {m0} kg')
    print('  State vector:')
    print(f'    r = [{r0[0]:10g}, {r0[1]:10g}, {r0[2]:10g}] (km)')
    print(f'      Radius = {np.linalg.norm(r0)}')
    print(f'    v = [{v0[0]:10g}, {v0[1]:10g}, {v0[2]:10g}] (km/s)')
    print(f'      Speed = {np.linalg.norm(v0)}\n')
    print(f'Thrust          = {T:12g} kN')
    print(f'Burn time       = {t_burn:12.6f} s')
    print(f'Mass after burn = {m1:12.6E} kg\n')
    print('\nEnd-of-burn-state vector:')
    print(f'    r = [{r1[0]:10g}, {r1[1]:10g}, {r1[2]:10g}] (km)')
    print(f'      Radius = {np.linalg.norm(r1)}')
    print(f'    v = [{v1[0]:10g}, {v1[1]:10g}, {v1[2]:10g}] (km/s)')
    print(f'      Speed = {np.linalg.norm(v1)}\n')
    print('\nPost-burn trajectory:')
    print(f'    Eccentricity = {e}')
    print(f'  Semimajor axis = {a} km')
    print(f'  Apogee state vector:')
    print(f'    r = [{ra[0]:17.10E}, {ra[1]:17.10E}, {ra[2]:17.10E}] (km)')
    print(f'      Radius = {np.linalg.norm(ra)}')
    print(f'    v = [{va[0]:17.10E}, {va[1]:17.10E}, {va[2]:17.10E}] (km/s)')
    print(f'      Speed = {np.linalg.norm(va)}')
    print('\n\n------------------------------------------------------\n\n')

if __name__ == '__main__':
    integrate_thrust()