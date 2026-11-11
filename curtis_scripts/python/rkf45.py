# ALGORITHM 1.3: NUMERICAL INTEGRATION OF A SYSTEM OF
# FIRST-ORDER DIFFERENTIAL EQUATIONS BY THE RUNGE-KUTTA-FEHLBERG
# 4(5) METHOD WITH ADAPTIVE SIZE CONTROL

import numpy as np

def rkf45(ode_function, tspan, y0, tolerance=1.e-8):
    '''
    This function uses the Runge-Kutta-Fehlberg 4(5) algorithm to
    integrate a system of first-order differential equations
    dy/dt = f(t,y).

    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    t               - time
    a               - Fehlberg coefficients for locating the six solution
                      points (nodes) within each time interval.
    b               - Fehlberg coupling coefficients for computing the
                      derivatives at each interior point
    c4              - Fehlberg coefficients for the fourth-order solution
    c5              - Fehlberg coefficients for the fifth-order solution
    tol             - allowable truncation error
    ode_function    - handle for user M-function in which the derivatives f
                      are computed
    tspan           - the vector [t0 tf] giving the time interval for the
                      solution
    t0              - initial time
    tf              - final time
    y0              - column vector of initial values of the vector y
    tout            - column vector of times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the correponding time in tout
    h               - time step
    hmin            - minimum allowable time step
    ti              - time at the beginning of a time step
    yi              - values of y at the beginning of a time step
    t_inner         - time within a given time step
    y_inner         - values of y witin a given time step
    te              - trucation error for each y at a given time step
    te_allowed      - allowable truncation error
    te_max          - maximum absolute value of the components of te
    ymax            - maximum absolute value of the components of y
    tol             - relative tolerance
    delta           - fractional change in step size
    eps             - unit roundoff error (the smallest number for which
                      1 + eps > 1)
    eps(x)          - the smallest number such that x + eps(x) = x

    User py-function required: ode_function
    '''
    a = np.array([0, 1/4, 3/8, 12/13, 1, 1/2])

    b = np.array([
        [0,        0,         0,         0,        0],
        [1/4,      0,         0,         0,        0],
        [3/32,     9/32,      0,         0,        0],
        [1932/2197, -7200/2197, 7296/2197, 0,      0],
        [439/216, -8,        3680/513,  -845/4104, 0],
        [-8/27,    2,       -3544/2565, 1859/4104, -11/40]])
    
    c4 = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
    c5 = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])

    t0, tf = tspan
    if not np.isfinite([t0, tf]).all() or tf < t0:
        raise ValueError("Integration requires finite times with tf >= t0.")
    if np.asarray(y0).ndim != 1 or len(y0) == 0 or not np.isfinite(y0).all():
        raise ValueError("Initial state must be a finite, nonempty vector.")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Tolerance must be finite and positive.")
    t = t0
    y = np.array(y0, dtype=float)
    tout = [t]
    yout = [y.copy()]
    h = (tf - t0) / 100  # Assumed initial time step.

    while t < tf:
        h = min(h, tf - t)
        if h <= 0 or t + h == t:
            raise RuntimeError("Integration step cannot advance time.")
        hmin = 16 * abs(np.spacing(t))
        ti = t
        yi = y.copy()

        #...Evaluate the time derivatives at six points within the current 
        # interval:
        f = np.zeros((len(y), 6))
        for i in range(6):
            t_inner = ti + a[i] * h
            y_inner = yi + h * sum(b[i, j] * f[:, j] for j in range(i))
            fi = np.asarray(ode_function(t_inner, y_inner), dtype=float)
            if fi.shape != y.shape or not np.isfinite(fi).all():
                raise ValueError("Derivative must match the state and be finite.")
            f[:, i] = fi

        #...Compute the maximum truncation error:
        te = h * f @ (c4 - c5)  # Difference between 4th and 
                                # 5th order solutions
        te_max = max(abs(te))

        #...Compute the allowable truncation error:
        ymax = max(abs(y))
        te_allowed = tolerance * max(ymax, 1.0)

        #...Compute the fractional change in step size:
        delta = 4.0 if te_max == 0 else 0.9 * (te_allowed / te_max)**(1/5)

        #...If the truncation error is in bounds, update the solution:
        if te_max <= te_allowed:
            t += h
            y = yi + h * f @ c5
            tout.append(t)
            yout.append(y.copy())

        #...Update the time step
        h = min(delta * h, 4 * h)
        if t < tf and h < hmin:
            raise RuntimeError("Required step size is below floating-point resolution.")

    return np.array(tout), np.array(yout)