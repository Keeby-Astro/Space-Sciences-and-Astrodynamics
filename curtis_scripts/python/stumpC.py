import numpy as np

def stumpC(z):
    '''
    This function evaluates the Stumpff function C(z) according
    to Equation 3.53.

    z - input argument
    c - value of C(z)

    User py-functions required: none
    '''  
    if abs(z) < 1.e-3:
        return 1/2 - z/24 + z**2/720 - z**3/40320 + z**4/3628800
    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        c = 1 / 2
    return c