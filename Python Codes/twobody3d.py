# Algorithm 2.1: Numerical solution for the motion of two bodies relative to
# an inertial frame.
import numpy as np
import matplotlib.pyplot as plt
import rkf45

def twobody3d():
    """
    This function solves the inertial two-body problem in three dimensions
    numerically using the RKF4(5) method.

    Parameters and Variables:
        G, m1, m2    - gravitational constant and masses of two bodies
        t0, tf       - initial and final time
        R1_0, V1_0   - initial position and velocity of m1
        R2_0, V2_0   - initial position and velocity of m2
        y0           - initial state vector
        t, y         - time vector and solution matrix
        X1, Y1, Z1   - position components of m1 over time
        X2, Y2, Z2   - position components of m2 over time
        XG, YG, ZG   - center of mass coordinates over time

    User Functions Required: rkf45, rates, output, common_axis_settings
    """
    
    # Constants
    G = 6.67259e-20  # universal gravitational constant (km^3/kg/s^2)

    # Input data
    m1 = 1.e26
    m2 = 1.e26
    t0 = 0
    tf = 480

    R1_0 = np.array([0, 0, 0])
    R2_0 = np.array([3000, 0, 0])
    V1_0 = np.array([10, 20, 30])
    V2_0 = np.array([0, 40, 0])

    # Initial state vector
    y0 = np.hstack((R1_0, R2_0, V1_0, V2_0))

    # Integrate the equations of motion
    t, y = rkf45(rates, [t0, tf], y0)

    # Output the results
    output(t, y, m1, m2)

# ––––––––––––––––––––––––
def rates(t, y):
    """
    This function calculates the accelerations based on the gravitational
    interaction between two bodies.

    Parameters:
        t : float
            Time (not used explicitly in this function).
        y : np.ndarray
            State vector containing positions and velocities of the two bodies.

    Returns:
        dydt : np.ndarray
            Derivative of the state vector.
    """
    
    G = 6.67259e-20  # gravitational constant (km^3/kg/s^2)
    m1 = 1.e26       # mass of the first body
    m2 = 1.e26       # mass of the second body

    # Extract positions and velocities from y
    R1 = y[0:3]
    R2 = y[3:6]
    V1 = y[6:9]
    V2 = y[9:12]

    # Distance between the two bodies
    r = np.linalg.norm(R2 - R1)

    # Accelerations due to gravitational force
    A1 = G * m2 * (R2 - R1) / r**3
    A2 = G * m1 * (R1 - R2) / r**3

    return np.hstack((V1, V2, A1, A2))

# –––––––––––––
def output(t, y, m1, m2):
    """
    Extracts and plots the trajectories of m1 and m2, as well as
    the center of mass.

    Parameters:
        t : np.ndarray
            Time vector.
        y : np.ndarray
            Solution matrix with positions and velocities of m1 and m2.
        m1, m2 : float
            Masses of the two bodies.
    """
    
    # Extract trajectories
    X1, Y1, Z1 = y[:, 0], y[:, 1], y[:, 2]
    X2, Y2, Z2 = y[:, 3], y[:, 4], y[:, 5]

    # Calculate center of mass
    XG = (m1 * X1 + m2 * X2) / (m1 + m2)
    YG = (m1 * Y1 + m2 * Y2) / (m1 + m2)
    ZG = (m1 * Z1 + m2 * Z2) / (m1 + m2)

    # Plot trajectories
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

# ––––––––––––––––––––––––---
def common_axis_settings():
    """
    Sets common axis properties for the 3D plots.
    """
    
    plt.text(0, 0, 0, 'o')
    plt.axis('equal')
    plt.view_init(elev=24, azim=60)  # Adjust for a better 3D view
    plt.grid(True)
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.gca().set_zlabel('Z (km)')
