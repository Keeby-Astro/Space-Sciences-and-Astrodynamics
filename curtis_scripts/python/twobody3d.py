# ALGORITHM 2.1: NUMERICAL SOLUTION OF THE TWO-BODY PROBLEM
# RELATIVE TO AN INERTIAL FRAME

import numpy as np
import matplotlib.pyplot as plt
import rkf45

def twobody3d():
    '''
    This function solves the inertial two-body problem in three dimensions
    numerically using the RKF4(5) method.

    G               - universal gravitational constant (km^3/kg/s^2)
    m1,m2           - the masses of the two bodies (kg)
    m               - the total mass (kg)
    t0              - initial time (s)
    tf              - final time (s)
    R1_0,V1_0       - 3 by 1 column vectors containing the components of tbe
                      initial position (km) and velocity (km/s) of m1
    R2_0,V2_0       - 3 by 1 column vectors containing the components of the
                      initial position (km) and velocity (km/s) of m2
    y0              - 12 by 1 column vector containing the initial values
                      of the state vectors of the two bodies:
                      [R1_0; R2_0; V1_0; V2_0]
    t               - column vector of the times at which the solution is found
    X1,Y1,Z1        - column vectors containing the X,Y and Z coordinates (km)
                      of m1 at the times in t
    X2,Y2,Z2        - column vectors containing the X,Y and Z coordinates (km)
                      of m2 at the times in t
    VX1, VY1, VZ1   - column vectors containing the X,Y and Z components
                      of the velocity (km/s) of m1 at the times in t
    VX2, VY2, VZ2   - column vectors containing the X,Y and Z components
                      of the velocity (km/s) of m2 at the times in t
    y               - a matrix whose 12 columns are, respectively,
                      X1,Y1,Z1; X2,Y2,Z2; VX1,VY1,VZ1; VX2,VY2,VZ2
    XG,YG,ZG        - column vectors containing the X,Y and Z coordinates (km)
                      the center of mass at the times in t

    User py-functions required: rkf45, 
    User subfunctions required: rates, output
    '''
    #...Input data
    m1 = 1.e26
    m2 = 1.e26
    t0 = 0
    tf = 480

    R1_0 = np.array([0, 0, 0])
    R2_0 = np.array([3000, 0, 0])

    V1_0 = np.array([10, 20, 30])
    V2_0 = np.array([0, 40, 0])
    #...End input data

    y0 = np.hstack((R1_0, R2_0, V1_0, V2_0))

    #...Integrate the equations of motion
    t, y = rkf45.rkf45(rates, [t0, tf], y0)

    #...Output the results
    output(t, y, m1, m2)

# ––––––––––––––
def rates(_, y):
# ––––––––––––––
    '''
    This function calculates the accelerations in Equations 2.19.

    t       - time
    y       - column vector containing the position and velocity vectors
              of the system at time t
    R1, R2  - position vectors of m1 & m2
    V1, V2  - velocity vectors of m1 & m2
    r       - magnitude of the relative position vector
    A1, A2  - acceleration vectors of m1 & m2
    dydt    - column vector containing the velocity and acceleration
              vectors of the system at time t
    '''
    G = 6.67259e-20 
    m1 = 1.e26
    m2 = 1.e26

    R1 = y[0:3]
    R2 = y[3:6]

    V1 = y[6:9]
    V2 = y[9:12]

    r = np.linalg.norm(R2 - R1)

    A1 = G * m2 * (R2 - R1) / r**3
    A2 = G * m1 * (R1 - R2) / r**3

    dydt = np.hstack((V1, V2, A1, A2))

    return dydt

# –––––––––––––––––––––––
def output(_, y, m1, m2):
# –––––––––––––––––––––––
    '''
    This function calculates the trajectory of the center of mass and
    plots
    (a) the motion of m1, m2 and G relative to the inertial frame
    (b) the motion of m2 and G relative to m1
    (c) the motion of m1 and m2 relative to G

    User subfunction required: common_axis_settings
    '''
    #...Extract the particle trajectories:
    X1, Y1, Z1 = y[:, 0], y[:, 1], y[:, 2]
    X2, Y2, Z2 = y[:, 3], y[:, 4], y[:, 5]

    #...Locate the center of mass at each time step:
    XG = (m1 * X1 + m2 * X2) / (m1 + m2)
    YG = (m1 * Y1 + m2 * Y2) / (m1 + m2)
    ZG = (m1 * Z1 + m2 * Z2) / (m1 + m2)

    #...Plot the trajectories:
    plt.figure(1)
    plt.title('Figure 2.3: Motion relative to the inertial frame')
    plt.plot(X1, Y1, Z1, '-r', label="m1")
    plt.plot(X2, Y2, Z2, '-g', label="m2")
    plt.plot(XG, YG, ZG, '-b', label="Center of Mass")
    plt.legend()
    common_axis_settings()

    plt.figure(2)
    plt.title('Figure 2.4a: Motion of m2 and G relative to m1')
    plt.plot(X2 - X1, Y2 - Y1, Z2 - Z1, '-g', label="m2 relative to m1")
    plt.plot(XG - X1, YG - Y1, ZG - Z1, '-b', label="G relative to m1")
    plt.legend()
    common_axis_settings()

    plt.figure(3)
    plt.title('Figure 2.4b: Motion of m1 and m2 relative to G')
    plt.plot(X1 - XG, Y1 - YG, Z1 - ZG, '-r', label="m1 relative to G")
    plt.plot(X2 - XG, Y2 - YG, Z2 - ZG, '-g', label="m2 relative to G")
    plt.legend()
    common_axis_settings()
    plt.show()

# ––––––––––––––––––––––––-
def common_axis_settings():
# ––––––––––––––––––––––––-
    '''
    This function establishes axis properties common to the several plots.
    '''
    ax = plt.figure().add_subplot(projection='3d')
    ax.text(0, 0, 0, 'o')
    ax.set_box_aspect([1,1,1])
    ax.view_init(elev=24, azim=60)
    ax.grid(True)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

if __name__ == '__main__':
    twobody3d()
