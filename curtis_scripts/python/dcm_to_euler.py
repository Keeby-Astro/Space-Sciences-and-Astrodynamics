# ALGORITHM 4.3: OBTAIN THE CLASSICAL EULER ANGLE SEQUENCE
# FROM A DIRECTION COSINE MATRIX

import numpy as np
import atan2d_0_360

def dcm_to_euler(Q):
    '''
    This function finds the angles of the classical Euler sequence
    R3(gamma)*R1(beta)*R3(alpha) from the direction cosine matrix.

    Q     - direction cosine matrix
    alpha - first angle of the sequence (deg)
    beta  - second angle of the sequence (deg)
    gamma - third angle of the sequence (deg)

    User py-function required: atan2d_0_360
    '''
    alpha = atan2d_0_360.atan2d_0_360(Q[2,0], -Q[2,1])
    beta  = np.degrees(np.arccos(np.clip(Q[2,2], -1, 1)))
    gamma = atan2d_0_360.atan2d_0_360(Q[0,2], Q[1,2])
    
    if np.hypot(Q[2,0], Q[2,1]) < 1.e-12:
        alpha = atan2d_0_360.atan2d_0_360(Q[0,1], Q[0,0])
        gamma = 0.0
    return alpha, beta, gamma