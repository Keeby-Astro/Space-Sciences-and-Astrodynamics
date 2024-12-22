import numpy as np

def quatinv(q):
    '''
    Compute the inverse of a quaternion or an array of quaternions,
    replicating MATLAB's quatinv behavior.

    q    - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    qinv - m-by-4 array of quaternion inverses.
    '''
    q = np.asarray(q)
    
    if q.ndim == 1:
        # Single quaternion case
        w, x, y, z = q
        norm_squared = w**2 + x**2 + y**2 + z**2
        qinv = np.array([w, -x, -y, -z]) / norm_squared
        return qinv
    elif q.ndim == 2:
        # Multiple quaternions case
        w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
        norm_squared = w**2 + x**2 + y**2 + z**2
        conjugate = np.column_stack((w, -x, -y, -z))
        qinv = conjugate / norm_squared[:, np.newaxis]
        return qinv
    else:
        raise ValueError("Input must be a 1-by-4 or m-by-4 array of quaternions.")