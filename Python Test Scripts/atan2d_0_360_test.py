import numpy as np
import atan2d_0_360

def atan2d_0_360_test():
    '''
    This script tests the atan2d_0_360 function for various inputs.

    User py-function required: atan2d_0_360
    '''
    # Inputs
    inputs = [
        [0, 0], [1, 0], [0, 1], [-1, 0],
        [0, -1], [1, 1], [-1, -1], [1, -1],
        [-1, 1], [np.sqrt(3), 1], [-np.sqrt(3), -1],
        [np.sqrt(3), -1], [-np.sqrt(3), 1]
    ]

    for i, (y, x) in enumerate(inputs, start=1):
        result = atan2d_0_360.atan2d_0_360(y, x)

        print(f"Test case {i}: ({y}, {x})")
        print(f"Result = {result}")

if __name__ == "__main__":
    atan2d_0_360_test()