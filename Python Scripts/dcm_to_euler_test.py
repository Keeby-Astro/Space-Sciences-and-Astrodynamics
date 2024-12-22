import numpy as np
import dcm_to_euler

def dcm_to_euler_test():
    '''
    This script tests the dcm_to_euler function using various input
    direction cosine matrices and verifies the output Euler angles.

    Q     - direction cosine matrix

    User py-function required: dcm_to_euler
    '''
    # Input Matrices
    input_matrices = [
        np.eye(3),                                                      # Identity matrix (no rotation)
        np.array([[0, -1,  0], [1,  0,  0], [0,  0,  1]]),              # 90-degree rotation about Z-axis
        np.array([[0,  0, -1], [0,  1,  0], [1,  0,  0]]),              # 90-degree rotation about X-axis
        np.array([[0,  1,  0], [-1,  0,  0], [0,  0,  1]]),             # -90-degree rotation about Z-axis
        np.array([[np.cos(np.radians(45)), -np.sin(np.radians(45)), 0], # 45-degree rotation about Z-axis
                  [np.sin(np.radians(45)),  np.cos(np.radians(45)), 0],
                  [0,                     0,                     1]])
    ]

    for i, Q in enumerate(input_matrices):
        # Compute Euler Angles
        alpha, beta, gamma = dcm_to_euler.dcm_to_euler(Q)

        # Results
        print('Input DCM:')
        print(Q)
        print(f'Output Euler Angles: alpha = {alpha:.2f}, beta = {beta:.2f}, gamma = {gamma:.2f} (degrees)')
        print()

if __name__ == "__main__":
    dcm_to_euler_test()