import numpy as np
import sv_from_coe
import rva_relative

def rva_relative_test():
    """
    This program uses the data of Example 7.1 to calculate the position,
    velocity and acceleration of an orbiting chaser B relative to an
    orbiting target A.

    mu               - gravitational parameter (km^3/s^2)
    deg              - conversion factor from degrees to radians

                     Spacecraft A & B:
    h_A, h_B         -   angular momentum (km^2/s)
    e_A, e_B         -   eccentricity
    i_A, i_B         -   inclination (radians)
    RAAN_A, RAAN_B   -   right ascension of the ascending node (radians)
    omega_A, omega_B -   argument of perigee (radians)
    theta_A, theta_B -   true anomaly (radians)

    rA, vA           - inertial position (km) and velocity (km/s) of A
    rB, vB           - inertial position (km) and velocity (km/s) of B
    r                - position (km) of B relative to A in A's
                        co-moving frame
    v                - velocity (km/s) of B relative to A in A's
                        co-moving frame
    a                - acceleration (km/s^2) of B relative to A in A's
                        co-moving frame

    User py-function required:   sv_from_coe, rva_relative
    """
    mu  = 398600.4418
    deg = np.pi/180

    # Input data:

    # Spacecraft A:
    h_A     = 52059
    e_A     = 0.025724
    i_A     = 60 * deg
    RAAN_A  = 40 * deg
    omega_A = 30 * deg
    theta_A = 40 * deg

    # Spacecraft B:
    h_B     = 52362
    e_B     = 0.0072696
    i_B     = 50 * deg
    RAAN_B  = 40 * deg
    omega_B = 120 * deg
    theta_B = 40 * deg

    # Compute the initial state vectors of A and B using Algorithm 4.5:
    rA, vA = sv_from_coe.sv_from_coe([h_A, e_A, RAAN_A, i_A, omega_A, theta_A], mu)
    rB, vB = sv_from_coe.sv_from_coe([h_B, e_B, RAAN_B, i_B, omega_B, theta_B], mu)

    # Compute relative position, velocity, and acceleration using Algorithm 7.1:
    r, v, a = rva_relative.rva_relative(rA, vA, rB, vB, mu)

    # Output
    print("\n\n--------------------------------------------------------\n\n")
    print("\nOrbital parameters of spacecraft A:")
    print(f"   angular momentum    = {h_A} (km^2/s)")
    print(f"   eccentricity        = {e_A}")
    print(f"   inclination         = {i_A / deg} (deg)")
    print(f"   RAAN                = {RAAN_A / deg} (deg)")
    print(f"   argument of perigee = {omega_A / deg} (deg)")
    print(f"   true anomaly        = {theta_A / deg} (deg)\n")

    print("State vector of spacecraft A:")
    print(f"   r = [{rA[0]}, {rA[1]}, {rA[2]}]")
    print(f"       (magnitude = {np.linalg.norm(rA)})")
    print(f"   v = [{vA[0]}, {vA[1]}, {vA[2]}]")
    print(f"       (magnitude = {np.linalg.norm(vA)})\n")

    print("Orbital parameters of spacecraft B:")
    print(f"   angular momentum    = {h_B} (km^2/s)")
    print(f"   eccentricity        = {e_B}")
    print(f"   inclination         = {i_B / deg} (deg)")
    print(f"   RAAN                = {RAAN_B / deg} (deg)")
    print(f"   argument of perigee = {omega_B / deg} (deg)")
    print(f"   true anomaly        = {theta_B / deg} (deg)\n")

    print("State vector of spacecraft B:")
    print(f"   r = [{rB[0]}, {rB[1]}, {rB[2]}]")
    print(f"       (magnitude = {np.linalg.norm(rB)})")
    print(f"   v = [{vB[0]}, {vB[1]}, {vB[2]}]")
    print(f"       (magnitude = {np.linalg.norm(vB)})\n")

    print("In the co-moving frame attached to A:")
    print(f"   Position of B relative to A = [{r[0]}, {r[1]}, {r[2]}]")
    print(f"      (magnitude = {np.linalg.norm(r)})\n")
    print(f"   Velocity of B relative to A = [{v[0]}, {v[1]}, {v[2]}]")
    print(f"      (magnitude = {np.linalg.norm(v)})\n")
    print(f"   Acceleration of B relative to A = [{a[0]}, {a[1]}, {a[2]}]")
    print(f"      (magnitude = {np.linalg.norm(a)})\n")
    print("\n\n--------------------------------------------------------\n\n")

if __name__ == "__main__":
    rva_relative_test()
