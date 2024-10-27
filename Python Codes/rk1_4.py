# Algorithm 1.1: Numerical integration of a system of first-order differential
# equations by choice of Runge-Kutta methods RK1, RK2, RK3. or RK4.
import numpy as np

def rk1_4(ode_function, tspan, y0, h, rk):
    """
    This function uses a selected Runge-Kutta procedure to integrate
    a system of first-order differential equations dy/dt = f(t,y).
    
    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    t               - time
    rk              - = 1 for RK1; = 2 for RK2; = 3 for RK3; = 4 for RK4
    n_stages        - the number of points within a time interval that
                      the derivatives are to be computed
    a               - coefficients for locating the solution points within
                      each time interval
    b               - coefficients for computing the derivatives at each
    interior point
    c               - coefficients for the computing solution at the end of
                      the time step
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
    ti              - time at the beginning of a time step
    yi              - values of y at the beginning of a time step
    t_inner         - time within a given time step
    y_inner         - values of y within a given time step

    User py-function required: ode_function
    """
    
    # Determine which of the four Runge-Kutta methods to use:
    if rk == 1:
        n_stages = 1
        a = np.array([0])
        b = np.array([[0]])
        c = np.array([1])
    elif rk == 2:
        n_stages = 2
        a = np.array([0, 1])
        b = np.array([[0, 0], [1, 0]])
        c = np.array([1/2, 1/2])
    elif rk == 3:
        n_stages = 3
        a = np.array([0, 1/2, 1])
        b = np.array([[0, 0, 0], [1/2, 0, 0], [-1, 2, 0]])
        c = np.array([1/6, 2/3, 1/6])
    elif rk == 4:
        n_stages = 4
        a = np.array([0, 1/2, 1/2, 1])
        b = np.array([[0, 0, 0, 0], [1/2, 0, 0, 0], [0, 1/2, 0, 0], [0, 0, 1, 0]])
        c = np.array([1/6, 1/3, 1/3, 1/6])
    else:
        raise ValueError("The parameter rk must have the value 1, 2, 3, or 4.")
    
    t0, tf = tspan
    t = t0
    y = np.array(y0)
    tout = [t]
    yout = [y.copy()]

    while t < tf:
        ti = t
        yi = y.copy()

        # Evaluate the time derivative(s) at the `n_stages` points within the current interval:
        f = np.zeros((len(y), n_stages))
        for i in range(n_stages):
            t_inner = ti + a[i] * h
            y_inner = yi + h * sum(b[i, j] * f[:, j] for j in range(i))
            f[:, i] = ode_function(t_inner, y_inner)

        h = min(h, tf - t)
        t += h
        y = yi + h * f @ c

        tout.append(t) # Add current time to the output time vector
        yout.append(y) # Add current y to the output solution matrix

    return np.array(tout), np.array(yout)
