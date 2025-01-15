# ALGORITHM 4.4: OBTAIN THE YAW, PITCH, AND ROLL ANGLES FROM A
# DIRECTION COSINE MATRIX

import numpy as np
import atan2d_0_360

def dcm_to_ypr(Q):
    '''
    This function finds the angles of the yaw-pitch-roll sequence
    R1(gamma)*R2(beta)*R3(alpha) from the direction cosine matrix.

    Q     - direction cosine matrix
    yaw   - yaw angle (deg)
    pitch - pitch angle (deg)
    roll  - roll angle (deg)

    User py-function required: atan2d_0_360
    '''
    yaw = atan2d_0_360.atan2d_0_360(Q[0, 1], Q[0, 0])
    pitch = np.degrees(np.arcsin(-Q[0, 2]))
    roll = atan2d_0_360.atan2d_0_360(Q[1, 2], Q[2, 2])
    
    return yaw, pitch, roll