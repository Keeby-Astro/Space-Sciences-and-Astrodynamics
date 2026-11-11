# ALGORITHM 1.2: NUMERICAL INTEGRATION BY HEUN'S PREDICTOR-CORRECTOR METHOD

import numpy as np

def heun(ode_function, tspan, y0, h):
    '''
    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    ode_function    - handle for the user M-function in which the derivatives
                      f are computed
    t               - time
    t0              - initial time
    tf              - final time
    tspan           - the vector [t0 tf] giving the time interval for the
                      solution
    h               - time step
    y0              - column vector of initial values of the vector y
    tout            - column vector of the times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the correponding time in tout
    feval           - a built-in MATLAB function which executes ’ode_function’
                      at the arguments t and y
    tol             - Maximum allowable relative error for determining
                      convergence of the corrector
    itermax         - maximum allowable number of iterations for corrector
                      convergence
    iter            - iteration number in the corrector convergence loop
    t1              - time at the beginning of a time step
    y1              - value of y at the beginning of a time step
    f1              - derivative of y at the beginning of a time step
    f2              - derivative of y at the end of a time step
    favg            - average of f1 and f2
    y2p             - predicted value of y at the end of a time step
    y2              - corrected value of y at the end of a time step
    err             - maximum relative error (for all components) between y2p
                      and y2 for given iteration
    eps             - unit roundoff error (the smallest number for which
                      1 + eps > 1). Used to avoid a zero denominator.

    User py-function required: ode_function
    '''
    tol = 1.e-6
    itermax = 100

    t0, tf = tspan
    if not np.isfinite([t0, tf]).all() or tf < t0:
        raise ValueError("Integration requires finite times with tf >= t0.")
    if np.asarray(y0).ndim != 1 or len(y0) == 0 or not np.isfinite(y0).all():
        raise ValueError("Initial state must be a finite, nonempty vector.")
    if not np.isfinite(h) or h <= 0:
        raise ValueError("Step size must be finite and positive.")
    t = t0
    y = np.array(y0, dtype=float)
    tout = [t]
    yout = [y.copy()]

    while t < tf:
        h = min(h, tf - t)
        if t + h == t:
            raise RuntimeError("Integration step cannot advance time.")
        t1 = t
        y1 = y
        f1 = np.asarray(ode_function(t1, y1), dtype=float)
        if f1.shape != y.shape or not np.isfinite(f1).all():
            raise ValueError("Derivative must match the state and be finite.")
        y2 = y1 + f1 * h
        t2 = t1 + h
        err = tol + 1
        iter = 0
        while err > tol and iter <= itermax:
            y2p = y2
            f2 = np.asarray(ode_function(t2, y2p), dtype=float)
            if f2.shape != y.shape or not np.isfinite(f2).all():
                raise ValueError("Derivative must match the state and be finite.")
            favg = (f1 + f2) / 2
            y2 = y1 + favg * h
            err = np.max(np.abs((y2 - y2p) / (y2 + np.finfo(float).eps)))
            iter += 1

        if iter > itermax:
            raise RuntimeError("Heun corrector failed to converge.")

        t += h
        y = y2
        tout.append(t) # adds t to the bottom of the column vector tout
        yout.append(y.copy()) # adds y' to the bottom of the matrix yout

    return np.array(tout), np.array(yout)