import fDot_and_gDot

def fDot_and_gDot_test():
    '''
    Outputs the time derivatives of the Lagrange f and g coefficients.

    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)

    User py-functions required: fDot_and_gDot
    '''
    mu = 398600.4418

    # Inputs
    x = 1.5      # Universal anomaly (km^0.5)
    r = 7000     # Radial position after time t (km)
    ro = 6878    # Radial position at time to (km)
    a = 1 / 8000 # Reciprocal of semimajor axis (1/km)

    fdot, gdot = fDot_and_gDot.fDot_and_gDot(x, r, ro, a, mu)

    # Display results
    print("Results:")
    print(f"f_dot = {fdot:.6e} 1/s")
    print(f"g_dot = {gdot:.6f} (dimensionless)")

if __name__ == '__main__':
    fDot_and_gDot_test()