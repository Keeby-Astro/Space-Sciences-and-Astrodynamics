import numpy as np
import quat_rotate

def quat_rotate_test():
    '''
    Evaluates known quaternions and vectors to compute
    that the rotation is performed properly.

    q - Quaternion to rotate by 1 x 4 array
    v - Vector to rotate by 1 x 3 array
    r - Rotated vector by 1 x 3 array

    User functions required: quat_rotate
    '''
    # Quaternions
    q1 = np.array([np.sqrt(2)/2, 0, 0, np.sqrt(2)/2])  # Quaternion for 90-degree Z-rotation
    q2 = np.array([0, 0, 1, 0])                        # Quaternion for 180-degree Y-rotation
    q3 = np.array([1, 0, 0, 0])                        # Identity quaternion (no rotation)
    q4 = np.array([np.sqrt(2)/2, np.sqrt(2)/2, 0, 0])  # Quaternion for 90-degree X-rotation

    # Vectors to rotate
    v1 = np.array([1, 0, 0])
    v2 = np.array([1, 0, 0])
    v3 = np.array([1, 2, 3])
    v4 = np.array([0, 1, 0])

    # Rotated vectors
    r1 = quat_rotate.quat_rotate(q1, v1)
    r2 = quat_rotate.quat_rotate(q2, v2)
    r3 = quat_rotate.quat_rotate(q3, v3)
    r4 = quat_rotate.quat_rotate(q4, v4)

    # Results
    print("90-degree rotation around Z-axis")
    print(f"Input Vector: [{v1[0]:.6f}, {v1[1]:.6f}, {v1[2]:.6f}]")
    print(f"Rotated Vector: [{r1[0]:.6f}, {r1[1]:.6f}, {r1[2]:.6f}]\n")

    print("180-degree rotation around Y-axis")
    print(f"Input Vector: [{v2[0]:.6f}, {v2[1]:.6f}, {v2[2]:.6f}]")
    print(f"Rotated Vector: [{r2[0]:.6f}, {r2[1]:.6f}, {r2[2]:.6f}]\n")

    print("No rotation")
    print(f"Input Vector: [{v3[0]:.6f}, {v3[1]:.6f}, {v3[2]:.6f}]")
    print(f"Rotated Vector: [{r3[0]:.6f}, {r3[1]:.6f}, {r3[2]:.6f}]\n")

    print("90-degree rotation around X-axis")
    print(f"Input Vector: [{v4[0]:.6f}, {v4[1]:.6f}, {v4[2]:.6f}]")
    print(f"Rotated Vector: [{r4[0]:.6f}, {r4[1]:.6f}, {r4[2]:.6f}]\n")

if __name__ == "__main__":
    quat_rotate_test()