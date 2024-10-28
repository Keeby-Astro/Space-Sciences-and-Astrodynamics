# Algorithm 11.2: Calculate the quaternion from the direction cosine matrix.
import numpy as np

def q_from_dcm(Q):
    """
    This function calculates the quaternion from the direction
    cosine matrix.

    Q - direction cosine matrix
    q - quaternion (where q[3] is the scalar part)
    """
    K3 = np.array([
        [Q[0,0]-Q[1,1]-Q[2,2], Q[1,0]+Q[0,1], Q[2,0]+Q[0,2], Q[1,2]-Q[2,1]],
        [Q[1,0]+Q[0,1], Q[1,1]-Q[0,0]-Q[2,2], Q[2,1]+Q[1,2], Q[2,0]-Q[0,2]],
        [Q[2,0]+Q[0,2], Q[2,1]+Q[1,2], Q[2,2]-Q[0,0]-Q[1,1], Q[0,1]-Q[1,0]],
        [Q[1,2]-Q[2,1], Q[2,0]-Q[0,2], Q[0,1]-Q[1,0], Q[0,0]+Q[1,1]+Q[2,2]]
    ]) / 3

    eigval, eigvec = np.linalg.eig(K3)

    i = np.argmax(np.diag(eigval))

    q = eigvec[:, i]
    return q
