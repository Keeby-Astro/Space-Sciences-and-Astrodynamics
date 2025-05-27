import numpy as np

def atan2d_0_360(y, x):
    '''
    This function calculates the arc tangent of y/x in degrees
    and places the result in the range [0, 360].

    t - angle in degrees
    
    User py-functions required: none
    '''
    if x == 0:
        if y == 0:
            t = 0
        elif y > 0:
            t = 90
        else:
            t = 270
    elif x > 0:
        if y >= 0:
            t = np.degrees(np.arctan2(y, x))
        else:
            t = np.degrees(np.arctan2(y, x)) + 360
    elif x < 0:
        if y == 0:
            t = 180
        else:
            t = np.degrees(np.arctan2(y, x)) + 180

    return t
