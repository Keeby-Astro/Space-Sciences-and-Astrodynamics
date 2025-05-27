# ALGORITHM 2.4: FIND THE ROOT OF A FUNCTION USING THE
# BISECTION METHOD

import numpy as np

def bisect(fun, xl, xu):
    '''
    This function evaluates a root of a function using
    the bisection method.

    tol  - error to within which the root is computed (default is 1.e-6)
    n    - number of iterations
    xl   - low end of the interval containing the root
    xu   - upper end of the interval containing the root
    i    - loop index
    xm   - mid-point of the interval from xl to xu
    fun  - name of the function whose root is being found
    fxl  - value of fun at xl
    fxm  - value of fun at xm
    root - the computed root

    User py-functions required: none
    '''
    tol=1.e-6
    n = np.ceil(np.log(abs(xu - xl) / tol) / np.log(2))

    for i in range(int(n)):
        xm = (xl + xu) / 2
        fxl = fun(xl)
        fxm = fun(xm)
        if fxl * fxm > 0:
            xl = xm
        else:
            xu = xm

    root = xm
    return root
