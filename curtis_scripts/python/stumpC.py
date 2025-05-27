import numpy as np

def stumpC(z):
    '''
    This function evaluates the Stumpff function C(z) according
    to Equation 3.53.

    z - input argument
    c - value of C(z)

    User py-functions required: none
    '''  
    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        c = 1 / 2
    return c