import numpy as np
import matplotlib.pyplot as plt
import heun

def heun_test():
    '''
    This program uses Heun's method with two different time steps to solve
    for and plot the response of a damped single degree of freedom
    spring-mass system to a sinusoidal forcing function, represented by

    x'' + 2*z*wn*x' + wn^2*x = (Fo/m)*sin(w*t)

    The numerical integration is done in the external function 'heun',
    which uses the subfunction 'rates' herein to compute the derivatives.

    x     - displacement (m)
    '     - shorthand for d/dt
    t     - time (s)
    wn    - natural circular frequency (radians/s)
    z     - damping factor
    Fo    - amplitude of the sinusoidal forcing function (N)
    m     - mass (kg)
    w     - forcing frequency (radians/s)
    t0    - initial time (s)
    tf    - final time (s)
    h     - uniform time step (s)
    tspan - row vector containing t0 and tf
    x0    - value of x at t0 (m)
    Dx0   - value of dx/dt at t0 (m/s)
    f0    - column vector containing x0 and Dx0
    t     - column vector of times at which the solution was computed
    f     - a matrix whose columns are:
            column 1: solution for x at the times in t
            column 2: solution for x' at the times in t

    User py-functions required: heun
    User py-subfunctions required: rates
    '''
    # System properties
    m = 1
    z = 0.03
    wn = 1
    Fo = 1
    w = 0.4 * wn

    # Time range
    t0 = 0
    tf = 110
    tspan = [t0, tf]

    # Initial conditions
    x0 = 0
    Dx0 = 0
    f0 = np.array([x0, Dx0])

    # Calculate and plot the solution for h = 1.0
    h = 1.0
    t1, f1 = heun.heun(lambda t, f: rates(t, f, m, z, wn, Fo, w), tspan, f0, h)

    # Calculate and plot the solution for h = 0.1
    h = 0.1
    t2, f2 = heun.heun(lambda t, f: rates(t, f, m, z, wn, Fo, w), tspan, f0, h)

    output(t1, f1, t2, f2)

def rates(t, f, m, z, wn, Fo, w):
    x = f[0]
    Dx = f[1]
    D2x = (Fo / m) * np.sin(w * t) - 2 * z * wn * Dx - wn**2 * x
    return np.array([Dx, D2x])

def output(t1, f1, t2, f2):
    plt.plot(t1, f1[:, 0], '-r', linewidth=0.5, label='h = 1.0')
    plt.plot(t2, f2[:, 0], '-k', linewidth=1, label='h = 0.1')
    plt.xlabel('time, s')
    plt.ylabel('x, m')
    plt.grid()
    plt.axis([0, 110, -2, 2])
    plt.legend()
    plt.show()

if __name__ == "__main__":
    heun_test()