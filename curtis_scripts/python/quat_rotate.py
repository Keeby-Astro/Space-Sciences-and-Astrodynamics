import numpy as np
import quatmultiply
import quatinv

def quat_rotate(q, v):
    '''
    quat_rotate rotates a vector by a unit quaternion.
    r = quat_rotate(q, v) calculates the rotated vector r for a
        quaternion q and a vector v.
    q   is a 1-by-4 array whose norm must be 1. q[0] is the scalar part
        of the quaternion.
    v   is a 1-by-3 array.
    r   is a 1-by-3 array.

    The 3-vector v is made into a pure quaternion 4-vector V = [0, v]. r is
    produced by the quaternion product R = q*V*qinv. r = [R[1], R[2], R[3]].

    py functions used: quatmultiply, quatinv.
    '''
    qinv = quatinv.quatinv(q)
    v_quat = np.concatenate(([0], v))
    r_quat = quatmultiply.quatmultiply(quatmultiply.quatmultiply(q, v_quat), qinv)
    r = r_quat[1:4]
    
    return r