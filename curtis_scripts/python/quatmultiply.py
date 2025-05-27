import numpy as np

def quatmultiply(q1, q2):
    '''
    Perform quaternion multiplication, replicating MATLAB's quatmultiply behavior.

    q1     - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    q2     - m-by-4 array of quaternions, where each quaternion is [w, x, y, z].
    qmulti - m-by-4 array of quaternion products.
    '''
    q1 = np.asarray(q1)
    q2 = np.asarray(q2)
    
    if q1.ndim == 1:
        q1 = q1[np.newaxis, :]
    if q2.ndim == 1:
        q2 = q2[np.newaxis, :]
    
    # Handle broadcasting if the number of rows differs
    if q1.shape[0] != q2.shape[0]:
        if q1.shape[0] == 1:
            q1 = np.tile(q1, (q2.shape[0], 1))
        elif q2.shape[0] == 1:
            q2 = np.tile(q2, (q1.shape[0], 1))
        else:
            raise ValueError("Input quaternions must have matching dimensions or one must be a single quaternion.")
    
    # Extract components of q1 and q2
    w1, x1, y1, z1 = q1[:, 0], q1[:, 1], q1[:, 2], q1[:, 3]
    w2, x2, y2, z2 = q2[:, 0], q2[:, 1], q2[:, 2], q2[:, 3]

    # Compute the product
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2
    z = w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2

    qmulti = np.column_stack((w, x, y, z))

    if qmulti.shape[0] == 1:
        return qmulti[0]
    return qmulti