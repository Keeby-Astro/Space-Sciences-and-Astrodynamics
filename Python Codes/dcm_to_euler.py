# Algorithm 4.3: Obtain the classical Euler angle sequence from a DCM.
import numpy as np
import atan2d_0_360

def dcm_to_euler(Q):
    """
    This function finds the angles of the classical Euler sequence
    R3(gamma)*R1(beta)*R3(alpha) from the direction cosine matrix.

    Q     - direction cosine matrix
    alpha - first angle of the sequence (deg)
    beta  - second angle of the sequence (deg)
    gamma - third angle of the sequence (deg)

    User py-function required: atan2d_0_360
    """
    alpha = atan2d_0_360(Q[2,0], -Q[2,1])
    beta  = np.degrees(np.arccos(Q[2,2]))
    gamma = atan2d_0_360(Q[0,2], Q[1,2])
    
    return alpha, beta, gamma
