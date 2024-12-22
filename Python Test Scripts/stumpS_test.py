import numpy as np
import stumpS

def stumpS_test():
    '''
    This function tests the stumpS function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User py-functions required: stumpS
    '''
    # Define a range of test values for z
    z_values = [-10, -1, -0.1, 0, 0.1, 1, 10]

    # Initialize a list to store results
    s_values = []

    # Loop through test values and compute C(z)
    for z in z_values:
        s_values.append(stumpS.stumpS(z))

    # Display the results
    print("z\t\tS(z)")
    print("-------------------")
    for z, s in zip(z_values, s_values):
        print(f"{z:.2f}\t\t{s:.6f}")

    # Verify specific cases with known outcomes
    assert abs(stumpS.stumpS(0) - (1/6)) < 1e-6, 'Test failed for z = 0'
    assert abs(stumpS.stumpS(1) - ((np.sqrt(1)) - np.sin(np.sqrt(1))) / (np.sqrt(1)**3)) < 1e-6, 'Test failed for z > 0'
    assert abs(stumpS.stumpS(-1) - (np.sinh(np.sqrt(1)) - np.sqrt(1) / (np.sqrt(1)**3))) < 1e-6, 'Test failed for z < 0'

if __name__ == '__main__':
    stumpS_test()