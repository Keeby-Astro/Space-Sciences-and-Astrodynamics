import numpy as np
import stumpC

def stumpC_test():
    '''
    This function tests the stumpC function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User py-functions required: stumpC
    '''
    # Define a range of test values for z
    z_values = [-10, -1, -0.1, 0, 0.1, 1, 10]

    # Initialize a list to store results
    c_values = []

    # Loop through test values and compute C(z)
    for z in z_values:
        c_values.append(stumpC.stumpC(z))

    # Display the results
    print("z\t\tC(z)")
    print("-------------------")
    for z, c in zip(z_values, c_values):
        print(f"{z:.2f}\t\t{c:.6f}")

    # Verify specific cases with known outcomes
    assert abs(stumpC.stumpC(0) - 0.5) < 1e-6, 'Test failed for z = 0'
    assert abs(stumpC.stumpC(1) - (1 - np.cos(np.sqrt(1))) / 1) < 1e-6, 'Test failed for z > 0'
    assert abs(stumpC.stumpC(-1) - (np.cosh(np.sqrt(1)) - 1) / 1) < 1e-6, 'Test failed for z < 0'

if __name__ == '__main__':
    stumpC_test()