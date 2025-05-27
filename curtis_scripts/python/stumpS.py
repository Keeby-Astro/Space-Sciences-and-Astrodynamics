import numpy as np

def stumpS(z):
    '''
    This function evaluates the Stumpff function S(z) according
    to Equation 3.52.

    z - input argument
    s - value of S(z)

    User py-functions required: none
    '''
    if z > 0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z))**3
    elif z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z))**3
    else:
        s = 1/6
    return s