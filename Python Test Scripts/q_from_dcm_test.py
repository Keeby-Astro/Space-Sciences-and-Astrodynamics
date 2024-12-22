import q_from_dcm

def q_from_dcm_test():
    '''
    Verifies the calculation of the quaternion from the direction
    cosine matrix using known inputs.

    Q - direction cosine matrix
    q - quaternion (where q[3] is the scalar part)

    User py-functions required: q_from_dcm
    '''
    Q1 = [[1, 0, 0],
          [0, 1, 0],
          [0, 0, 1]]  # Identity matrix
    Q2 = [[1,  0,  0],
          [0, -1,  0],
          [0,  0, -1]]  # 180-degree rotation about x-axis
    Q3 = [[-1,  0,  0],
          [ 0,  1,  0],
          [ 0,  0, -1]]  # 180-degree rotation about y-axis
    Q4 = [[-1,  0,  0],
          [ 0, -1,  0],
          [ 0,  0,  1]]  # 180-degree rotation about z-axis

    # Calculate the quaternions using the function
    q1 = q_from_dcm.q_from_dcm(Q1)
    q2 = q_from_dcm.q_from_dcm(Q2)
    q3 = q_from_dcm.q_from_dcm(Q3)
    q4 = q_from_dcm.q_from_dcm(Q4)

    # Display the results
    print("Results for direction cosine matrix Q1:")
    print(f"Quaternion: [{q1[0]:.6f}, {q1[1]:.6f}, {q1[2]:.6f}, {q1[3]:.6f}]")

    print("Results for direction cosine matrix Q2:")
    print(f"Quaternion: [{q2[0]:.6f}, {q2[1]:.6f}, {q2[2]:.6f}, {q2[3]:.6f}]")

    print("Results for direction cosine matrix Q3:")
    print(f"Quaternion: [{q3[0]:.6f}, {q3[1]:.6f}, {q3[2]:.6f}, {q3[3]:.6f}]")

    print("Results for direction cosine matrix Q4:")
    print(f"Quaternion: [{q4[0]:.6f}, {q4[1]:.6f}, {q4[2]:.6f}, {q4[3]:.6f}]")

if __name__ == "__main__":
    q_from_dcm_test()