# ALGORITHM 2.4: FIND THE ROOT OF A FUNCTION USING THE
# BISECTION METHOD

import math

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
    tol = 1.e-6
    if not all(math.isfinite(value) for value in [xl, xu]) or xl > xu:
        raise ValueError("Bounds must be finite and ordered.")
    fxl = fun(xl)
    fxu = fun(xu)
    if not all(math.isfinite(value) for value in [fxl, fxu]):
        raise ValueError("Function values must be finite.")
    if fxl == 0:
        return xl
    if fxu == 0:
        return xu
    if (fxl < 0) == (fxu < 0):
        raise ValueError("Bounds must bracket a root.")

    while xu - xl > 2 * tol:
        xm = xl + (xu - xl) / 2
        if xm == xl or xm == xu:
            raise RuntimeError("Bisection cannot refine the interval at floating-point resolution.")
        fxm = fun(xm)
        if not math.isfinite(fxm):
            raise ValueError("Function values must be finite.")
        if fxm == 0:
            return xm
        if (fxl < 0) == (fxm < 0):
            xl = xm
            fxl = fxm
        else:
            xu = xm

    root = xl + (xu - xl) / 2
    return root