# ALGORITHM 2.2: NUMERICAL SOLUTION OF THE TWO-BODY RELATIVE MOTION PROBLEM

import numpy as np
import rkf45
import matplotlib.pyplot as plt
from datetime import datetime

def orbit():
    '''
    This function computes the orbit of a spacecraft by using rkf45 to
    numerically integrate Equation 2.22.

    It also plots the orbit and computes the times at which the maximum
    and minimum radii occur and the speeds at those times.

    hours     - converts hours to seconds
    G         - universal gravitational constant (km^3/kg/s^2)
    m1        - planet mass (kg)
    m2        - spacecraft mass (kg)
    mu        - gravitational parameter (km^3/s^2)
    R         - planet radius (km)
    r0        - initial position vector (km)
    v0        - initial velocity vector (km/s)
    t0, tf    - initial and final times (s)
    y0        - column vector containing r0 and v0
    t         - array of times at which the solution is found
    y         - a matrix whose columns are:
                columns 1, 2, and 3:
                    The solution for the x, y, and z components of the
                    position vector r at the times in t
                columns 4, 5, and 6:
                    The solution for the x, y, and z components of the
                    velocity vector v at the times in t
    r         - magnitude of the position vector at the times in t
    imax      - component of r with the largest value
    rmax      - largest value of r
    imin      - component of r with the smallest value
    rmin      - smallest value of r
    v_at_rmax - speed where r = rmax
    v_at_rmin - speed where r = rmin

    User py-function required: rkf45
    User subfunctions required: rates, output
    '''

    hours = 3600
    G     = 6.6742e-20

    #...Input data:
    #   Earth:
    m1 = 5.974e24  # Mass of Earth (kg)
    R  = 6378.14   # Radius of Earth (km)
    m2 = 1000      # Mass of spacecraft (kg)

    r0 = np.array([8000, 0, 6000])  # Initial position vector (km)
    v0 = np.array([0, 7, 0])        # Initial velocity vector (km/s)

    t0 = 0
    tf = 4 * hours
    #...End input data

    #...Numerical integration:
    mu = G * (m1 + m2)  # Gravitational parameter (km^3/s^2)
    y0 = np.concatenate((r0, v0))
    # t, y = rkf45(rates, [t0, tf], y0)
    t, y = rkf45.rkf45(rates, [t0, tf], y0)
    #...Output the results:
    output(t, y, mu, R, r0, v0, tf, hours)

    return


def rates(_, f):
    '''
    This function calculates the acceleration vector using Equation 2.22.
    
    t          - time
    f          - column vector containing the position vector and the
                 velocity vector at time t
    x, y, z    - components of the position vector r
    r          - the magnitude of the position vector
    vx, vy, vz - components of the velocity vector v
    ax, ay, az - components of the acceleration vector a
    dydt       - column vector containing the velocity and acceleration
                 components
    '''
    G = 6.6742e-20
    m1 = 5.974e24
    m2 = 1000

    mu = G * (m1 + m2)
        
    x, y, z = f[0], f[1], f[2]
    vx, vy, vz = f[3], f[4], f[5]

    r = np.linalg.norm([x, y, z])

    ax = -mu * x / r**3
    ay = -mu * y / r**3
    az = -mu * z / r**3

    dydt = np.array([vx, vy, vz, ax, ay, az])
    return dydt

def output(t, y_data, _, R, r0, v0, tf, hours):
    '''
    # This function computes the maximum and minimum radii, the times they
    # occur, and the speed at those times. It prints those results to
    # the command window and plots the orbit.
    '''

    r = np.linalg.norm(y_data[:, :3], axis=1)

    rmax = np.max(r)
    imax = np.argmax(r)
    rmin = np.min(r)
    imin = np.argmin(r)

    v_at_rmax = np.linalg.norm(y_data[imax, 3:])
    v_at_rmin = np.linalg.norm(y_data[imin, 3:])

    # Output to command window
    print("\n\n--------------------------------------------------------------")
    print("\n Earth Orbit")
    print(f" {datetime.now()}")
    print(f"\n The initial position is [{r0[0]}, {r0[1]}, {r0[2]}] (km).")
    print(f" Magnitude = {np.linalg.norm(r0)} km")
    print(f"\n The initial velocity is [{v0[0]}, {v0[1]}, {v0[2]}] (km/s).")
    print(f" Magnitude = {np.linalg.norm(v0)} km/s")
    print(f"\n Initial time = 0 h.\n Final time = {tf/hours} h.")
    print(f"\n The minimum altitude is {rmin-R} km at time = {t[imin]/hours} h.")
    print(f" The speed at that point is {v_at_rmin} km/s.")
    print(f"\n The maximum altitude is {rmax-R} km at time = {t[imax]/hours} h.")
    print(f" The speed at that point is {v_at_rmax} km/s.")
    print("\n--------------------------------------------------------------\n")

    # Plot the results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Draw the planet
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = R * np.cos(u) * np.sin(v)
    y = R * np.sin(u) * np.sin(v)
    z = R * np.cos(v)
    ax.plot_surface(x, y, z, color='blue', edgecolor='none')

    # Draw the orbit
    ax.plot(y_data[:, 0], y_data[:, 1], y_data[:, 2], 'k')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('Orbit Plot')
    max_range = np.ptp(y_data[:, :3], axis=0).max() / 2.0
    mid = np.mean(y_data[:, :3], axis=0)
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)

    # Viewpoint
    ax.view_init(elev=20, azim=135)
    plt.show()

if __name__ == "__main__":
    orbit()