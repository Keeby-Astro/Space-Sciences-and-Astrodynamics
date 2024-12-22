import numpy as np
import matplotlib.pyplot as plt

def threebody():
    '''
    This program presents the graphical solution of the motion of three
    bodies in the plane for data provided in the input definitions below.

    Python's ode45 Runge-Kutta solver is used.

    G                           - gravitational constant (km^3/kg/s^2)
    t0, tf                      - initial and final times (s)
    m1, m2, m3                  - masses of the three bodies (kg)
    m                           - total mass (kg)
    X1,Y1; X2,Y2; X3,Y3         - coordinates of the three masses (km)
    VX1,VY1; VX2,VY2; VX3,VY3   - velocity components of the three
                                    masses (km/s)
    XG, YG                      - coordinates of the center of mass (km)
    y0                          - column vector of the initial conditions
    t                           - column vector of times at which the solution
                                    was computed
    y                           - matrix, the columns of which contain the
                                    position and velocity components evaluated at
                                    the times t(:):
                                    y(:,1) , y(:, 2) = X1(:), Y1(:)
                                    y(:,3) , y(:, 4) = X2(:), Y2(:)
                                    y(:,5) , y(:, 6) = X3(:), Y3(:)
                                    y(:,7) , y(:, 8) = VX1(:), VY1(:)
                                    y(:,9) , y(:,10) = VX2(:), VY2(:)
                                    y(:,11), y(:,12) = VX3(:), VY3(:)
    
    User py-functions required: none
    User py-subfunctions required: rates, plotit
    '''
    G = 6.67259e-20

    # Input data
    m1 = 1.e29
    m2 = 1.e29
    m3 = 1.e29

    t0 = 0
    tf = 67000
    dt = 1  # Time step for Euler method

    X1, Y1 = 0, 0
    X2, Y2 = 300000, 0
    X3, Y3 = 2 * X2, 0

    VX1, VY1 = 0, 0
    VX2, VY2 = 250, 250
    VX3, VY3 = 0, 0

    m = m1 + m2 + m3
    y0 = [X1, Y1, X2, Y2, X3, Y3, VX1, VY1, VX2, VY2, VX3, VY3]

    # Solve the ODE using Euler method
    t = np.arange(t0, tf, dt)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0

    for i in range(1, len(t)):
        dydt = rates(t[i-1], y[i-1], G, m1, m2, m3)
        y[i] = y[i-1] + np.array(dydt) * dt

    X1, Y1 = y[:, 0], y[:, 1]
    X2, Y2 = y[:, 2], y[:, 3]
    X3, Y3 = y[:, 4], y[:, 5]

    # Locate the center of mass at each time step
    XG = (m1 * X1 + m2 * X2 + m3 * X3) / m
    YG = (m1 * Y1 + m2 * Y2 + m3 * Y3) / m

    # Coordinates of each particle relative to the center of mass
    X1G, Y1G = X1 - XG, Y1 - YG
    X2G, Y2G = X2 - XG, Y2 - YG
    X3G, Y3G = X3 - XG, Y3 - YG

    plotit(X1, Y1, X2, Y2, X3, Y3, XG, YG, X1G, Y1G, X2G, Y2G, X3G, Y3G)

def rates(t, y, G, m1, m2, m3):
    '''
    This function evaluates the acceleration of each member of a planar
    3-body system at time t from their positions and velocities
    at that time.

    t                           - time (s)
    y                           - column vector containing the position and
                                    velocity components of the three masses
                                    at time t
    R12                         - cube of the distance between m1 and m2 (km^3)
    R13                         - cube of the distance between m1 and m3 (km^3)
    R23                         - cube of the distance between m2 and m3 (km^3)
    AX1,AY1; AX2,AY2; AX3,AY3   - acceleration components of each mass (km/s^2)
    dydt                        - column vector containing the velocity and
                                    acceleration components of the three
                                    masses at time t
    '''
    X1, Y1, X2, Y2, X3, Y3 = y[0], y[1], y[2], y[3], y[4], y[5]
    VX1, VY1, VX2, VY2, VX3, VY3 = y[6], y[7], y[8], y[9], y[10], y[11]

    R12 = np.linalg.norm([X2 - X1, Y2 - Y1])**3
    R13 = np.linalg.norm([X3 - X1, Y3 - Y1])**3
    R23 = np.linalg.norm([X3 - X2, Y3 - Y2])**3

    AX1 = G * m2 * (X2 - X1) / R12 + G * m3 * (X3 - X1) / R13
    AY1 = G * m2 * (Y2 - Y1) / R12 + G * m3 * (Y3 - Y1) / R13
    AX2 = G * m1 * (X1 - X2) / R12 + G * m3 * (X3 - X2) / R23
    AY2 = G * m1 * (Y1 - Y2) / R12 + G * m3 * (Y3 - Y2) / R23
    AX3 = G * m1 * (X1 - X3) / R13 + G * m2 * (X2 - X3) / R23
    AY3 = G * m1 * (Y1 - Y3) / R13 + G * m2 * (Y2 - Y3) / R23

    return [VX1, VY1, VX2, VY2, VX3, VY3, AX1, AY1, AX2, AY2, AX3, AY3]

def plotit(X1, Y1, X2, Y2, X3, Y3, XG, YG, X1G, Y1G, X2G, Y2G, X3G, Y3G):
    plt.figure(1)
    plt.title('Figure 2.4: Motion relative to the inertial frame', fontweight='bold', fontsize=12)
    plt.plot(XG, YG, '--k', linewidth=0.25)
    plt.plot(X1, Y1, 'r', linewidth=0.5)
    plt.plot(X2, Y2, 'g', linewidth=0.75)
    plt.plot(X3, Y3, 'b', linewidth=1.00)
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.grid(True)
    plt.axis('equal')

    plt.figure(2)
    plt.title('Figure 2.5: Motion relative to the center of mass', fontweight='bold', fontsize=12)
    plt.plot(X1G, Y1G, 'r', linewidth=0.5)
    plt.plot(X2G, Y2G, '--g', linewidth=0.75)
    plt.plot(X3G, Y3G, 'b', linewidth=1.00)
    plt.xlabel('X(km)')
    plt.ylabel('Y(km)')
    plt.grid(True)
    plt.axis('equal')

    plt.show()

if __name__ == "__main__":
    threebody()