import numpy as np
import dcm_to_ypr

def dcm_to_ypr_test():
    '''
    This script tests the dcm_to_ypr function using various input
    direction cosine matrices and verifies the output yaw, pitch, and roll angles.

    Q     - direction cosine matrix

    User py-function required: dcm_to_ypr
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
        # Compute Yaw, Pitch, and Roll Angles
        yaw, pitch, roll = dcm_to_ypr.dcm_to_ypr(Q)

        # Results
        print('Input DCM:')
        print(Q)
        print(f'Output Yaw, Pitch, and Roll Angles: yaw = {yaw:.2f}, pitch = {pitch:.2f}, roll = {roll:.2f} (degrees)')
        print()

if __name__ == "__main__":
    dcm_to_ypr_test()