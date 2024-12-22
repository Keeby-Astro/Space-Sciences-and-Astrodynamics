import numpy as np
import matplotlib.pyplot as plt
import rkf45
import sv_from_coe
import rv_from_r0v0

def relative_motion():
    '''
    This function plots the motion of chaser B relative to target A
    for the data.

    mu       - gravitational parameter (km^3/s^2)
    RE       - radius of the Earth (km)

             Target orbit at time t = 0:
    rp       - perigee radius (km)
    e        - eccentricity
    i        - inclination (rad)
    RA       - right ascension of the ascending node (rad)
    omega    - argument of perigee (rad)
    theta    - true anomaly (rad)
    ra       - apogee radius (km)
    h        - angular momentum (km^2/s)
    a        - semimajor axis (km)
    T        - period (s)
    n        - mean motion (rad/s)

    dr0, dv0 - initial relative position (km) and relative velocity (km/s)
               of B in the co-moving frame
    t0, tf   - initial and final times (s) for the numerical integration
    R0, V0   - initial position (km) and velocity (km/s) of A in the
               geocentric equatorial frame
    y0       - column vector containing r0, v0

    User py-functions required:  sv_from_coe, rkf45
    User subfunctions required: rates
    '''
    global mu
    mu = 398600.4418
    RE = 6378.14

    # Input data:
    # Prescribed initial orbital parameters of target A:
    rp = RE + 300
    e = 0.1
    i = 0
    RA = 0
    omega = 0
    theta = 0

    # Additional computed parameters:
    ra = rp * (1 + e) / (1 - e)
    h = np.sqrt(2 * mu * rp * ra / (ra + rp))
    a = (rp + ra) / 2
    T = 2 * np.pi / np.sqrt(mu) * a**1.5
    n = 2 * np.pi / T

    # Prescribed initial state vector of chaser B in the co-moving frame:
    dr0 = np.array([-1, 0, 0])
    dv0 = np.array([0, -2 * n * dr0[0], 0])
    t0 = 0
    tf = 5 * T

    # Calculate the target's initial state vector using sv_from_coe:
    R0, V0 = sv_from_coe.sv_from_coe([h, e, RA, i, omega, theta], mu)

    # Initial state vector of B's orbit relative to A:
    y0 = np.concatenate((dr0, dv0))

    # Integrate Equations using rkf45:
    t, y = rkf45.rkf45(rates, [t0, tf], y0)

    plotit(t, y)

    return

def rates(t, f):
    '''
    This function computes the components of f(t,y).

    t             - time
    f             - column vector containing the relative position and
                    velocity vectors of B at time t
    R, V          - updated state vector of A at time t
    X, Y, Z       - components of R
    VX, VY, VZ    - components of V
    R_            - magnitude of R
    RdotV         - dot product of R and V
    h             - magnitude of the specific angular momentum of A

    dx , dy , dz  - components of the relative position vector of B
    dvx, dvy, dvz - components of the relative velocity vector of B
    dax, day, daz - components of the relative acceleration vector of B
    dydt          - column vector containing the relative velocity
                    and acceleration components of B at time t

    User py-function required: rv_from_r0v0
    '''
    global R0, V0

    mu = 398600.4418
    RE = 6378.14

    # Input data:
    # Prescribed initial orbital parameters of target A:
    rp = RE + 300
    e = 0.1
    i = 0
    RA = 0
    omega = 0
    theta = 0

    # Additional computed parameters:
    ra = rp * (1 + e) / (1 - e)
    h = np.sqrt(2 * mu * rp * ra / (ra + rp))

    # Calculate the target's initial state vector using sv_from_coe:
    R0, V0 = sv_from_coe.sv_from_coe([h, e, RA, i, omega, theta], mu)

    # Update the state vector of the target orbit using rv_from_r0v0:
    R, V = rv_from_r0v0.rv_from_r0v0(R0, V0, t, mu)

    X, Y, Z = R
    VX, VY, VZ = V

    R_ = np.linalg.norm([X, Y, Z])
    RdotV = np.dot([X, Y, Z], [VX, VY, VZ])
    h = np.linalg.norm(np.cross([X, Y, Z], [VX, VY, VZ]))

    dx, dy, dz = f[0:3]
    dvx, dvy, dvz = f[3:6]

    dax = (2 * mu / R_**3 + h**2 / R_**4) * dx - 2 * RdotV / R_**4 * h * dy + 2 * h / R_**2 * dvy
    day = -(mu / R_**3 - h**2 / R_**4) * dy + 2 * RdotV / R_**4 * h * dx - 2 * h / R_**2 * dvx
    daz = -mu / R_**3 * dz

    dydt = np.array([dvx, dvy, dvz, dax, day, daz])

    return dydt

def plotit(t, y):
    '''
    Plot the trajectory of B relative to A.
    '''
    plt.figure()
    plt.plot(y[:, 1], y[:, 0])
    plt.axis('equal')
    plt.axis([0, 40, -5, 5])
    plt.xlabel('y (km)')
    plt.ylabel('x (km)')
    plt.grid(True)
    plt.box(True)

    # Label the start of B's trajectory relative to A:
    plt.text(y[0, 1], y[0, 0], 'o')
    plt.show()

if __name__ == '__main__':
    relative_motion()