import dcm_from_q

def dcm_from_q_test():
    '''
    Verifies the direction cosine matrix from known quaternions.

    q - quaternion (where q[3] is the scalar part)
    Q - direction cosine matrix

    User py-functions required: dcm_from_q
    '''
    q1 = [1, 0, 0, 0]  # Identity quaternion
    q2 = [0, 1, 0, 0]  # 180-degree rotation about x-axis
    q3 = [0, 0, 1, 0]  # 180-degree rotation about y-axis
    q4 = [0, 0, 0, 1]  # 180-degree rotation about z-axis

    # Calculate the direction cosine matrices using the function
    Q1 = dcm_from_q.dcm_from_q(q1)
    Q2 = dcm_from_q.dcm_from_q(q2)
    Q3 = dcm_from_q.dcm_from_q(q3)
    Q4 = dcm_from_q.dcm_from_q(q4)

    # Display the results
    print(f"Results for quaternion q1: {q1}")
    print("Direction Cosine Matrix:")
    for row in Q1:
        print(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}")

    print(f"Results for quaternion q2: {q2}")
    print("Direction Cosine Matrix:")
    for row in Q2:
        print(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}")

    print(f"Results for quaternion q3: {q3}")
    print("Direction Cosine Matrix:")
    for row in Q3:
        print(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}")

    print(f"Results for quaternion q4: {q4}")
    print("Direction Cosine Matrix:")
    for row in Q4:
        print(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}")

if __name__ == "__main__":
    dcm_from_q_test()