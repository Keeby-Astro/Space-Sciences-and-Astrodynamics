import numpy as np
import matplotlib.pyplot as plt
import rv_from_r0v0
import sv_from_coe
import rva_relative

def rva_relative_test_plot():
    '''
    This function produces a 3D plot of the motion of spacecraft B
    relative to A.

    User py-functions required: rv_from_r0v0, sv_from_coe, rva_relative
    '''

    # Gravitational parameter and Earth radius
    mu = 398600.4418
    RE = 6378.14

    # Conversion factor from degrees to radians
    deg = np.pi / 180

    # Input data
    # Initial orbital parameters (angular momentum, eccentricity,
    # inclination, RAAN, argument of perigee, and true anomaly)
    # Spacecraft A
    h_A = 52059
    e_A = 0.025724
    i_A = 60 * deg
    RAAN_A = 40 * deg
    omega_A = 30 * deg
    theta_A = 40 * deg

    # Spacecraft B
    h_B = 52362
    e_B = 0.0072696
    i_B = 50 * deg
    RAAN_B = 40 * deg
    omega_B = 120 * deg
    theta_B = 40 * deg

    vdir = [1, 1, 1]

    # Compute the initial state vectors of A and B
    rA0, vA0 = sv_from_coe.sv_from_coe([h_A, e_A, RAAN_A, i_A, omega_A, theta_A], mu)
    rB0, vB0 = sv_from_coe.sv_from_coe([h_B, e_B, RAAN_B, i_B, omega_B, theta_B], mu)

    h0 = np.cross(rA0, vA0)

    # Period of A
    TA = 2 * np.pi / mu**2 * (h_A / np.sqrt(1 - e_A**2))**3

    # Number of time steps per period of A's orbit
    n = 100

    # Time step as a fraction of A's period
    dt = TA / n

    # Number of periods of A's orbit for which the trajectory will be plotted
    n_Periods = 60

    # Initialize the time
    t = -dt

    # Initialize arrays to store results
    x = []
    y = []
    z = []
    r = []
    T = []

    # Generate the trajectory of B relative to A
    for count in range(n_Periods * n):
        # Update the time
        t += dt

        # Update the state vector of both orbits
        rA, vA = rv_from_r0v0.rv_from_r0v0(rA0, vA0, t, mu)
        rB, vB = rv_from_r0v0.rv_from_r0v0(rB0, vB0, t, mu)

        # Compute r_rel
        r_rel, v_rel, a_rel = rva_relative.rva_relative(rA, vA, rB, vB, mu)

        # Store the components of the relative position vector
        x.append(r_rel[0])
        y.append(r_rel[1])
        z.append(r_rel[2])
        r.append(np.linalg.norm(r_rel))
        T.append(t)

    # Plot the trajectory of B relative to A
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z)
    ax.set_box_aspect([8, 8, 1])  # Equal aspect ratio

    # Draw the co-moving x, y, and z axes
    ax.quiver(0, 0, 0, 2000, 0, 0, color='red', label='x-axis')
    ax.quiver(0, 0, 0, 0, 3500, 0, color='green', label='y-axis')
    ax.quiver(0, 0, 0, 0, 0, 2000, color='blue', label='z-axis')

    # Label the origin of the moving frame attached to A
    ax.text(0, 0, 0, 'A', color='black')

    # Label the start of B's relative trajectory
    ax.text(x[0], y[0], z[0], 'B', color='black')

    # Draw the initial position vector of B
    ax.quiver(0, 0, 0, 0.5*x[0], 0.5*y[0], 0.5*z[0], color='purple', label='Initial position')

    # Set plot labels and grid
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.grid(True)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    rva_relative_test_plot()