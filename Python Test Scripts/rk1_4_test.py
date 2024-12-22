import numpy as np
import matplotlib.pyplot as plt
import rk1_4

def rk1_4_test():
    '''
    This function uses the RK1 through RK4 methods with two
    different time steps each to solve for and plot the response
    of a damped single degree of freedom spring-mass system to
    a sinusoidal forcing function, represented by

    x'' + 2*z*wn*x' + wn^2*x = (Fo/m)*sin(w*t)
    
    The numerical integration is done by the external
    function 'rk1_4', which uses the subfunction 'rates'
    herein to compute the derivatives.

    This function also plots the exact solution for comparison.

    x           - displacement (m)
    '           - shorthand for d/dt
    t           - time (s)
    wn          - natural circular frequency (radians/s)
    z           - damping factor
    wd          - damped natural frequency
    Fo          - amplitude of the sinusoidal forcing function (N)
    m           - mass (kg)
    w           - forcing frequency (radians/s)
    t0          - initial time (s)
    tf          - final time (s)
    h           - uniform time step (s)
    tspan       - a row vector containing t0 and tf
    x0          - value of x at t0 (m)
    x_dot0      - value of dx/dt at t0 (m/s)
    f0          - column vector containing x0 and x_dot0
    rk          - = 1 for RK1; = 2 for RK2; = 3 for RK3; = 4 for RK4
    t           - solution times for the exact solution
    t1, ...,t4  - solution times for RK1,...,RK4 for smaller
    t11,...,t41 - solution times for RK1,...,RK4 for larger h
    f1, ...,f4  - solution vectors for RK1,...,RK4 for smaller h
    f11,...,f41 - solution vectors for RK1,...,RK4 for larger h

    User py-functions required: rk1_4
    User py-subfunctions required: rates
    '''
    # Input data
    m = 1
    z = 0.03
    wn = 1
    Fo = 1
    w = 0.4 * wn

    x0 = 0
    x_dot0 = 0
    f0 = np.array([x0, x_dot0])

    t0 = 0
    tf = 110
    tspan = [t0, tf]
    # End input data

    # Solve using RK1 through RK4, using the same and a larger
    # time step for each method
    rk = 1
    h = 0.01
    t1, f1 = rk1_4(rates, tspan, f0, h, rk)
    h = 0.1
    t11, f11 = rk1_4(rates, tspan, f0, h, rk)

    rk = 2
    h = 0.1
    t2, f2 = rk1_4(rates, tspan, f0, h, rk)
    h = 0.5
    t21, f21 = rk1_4(rates, tspan, f0, h, rk)

    rk = 3
    h = 0.5
    t3, f3 = rk1_4(rates, tspan, f0, h, rk)
    h = 1.0
    t31, f31 = rk1_4(rates, tspan, f0, h, rk)

    rk = 4
    h = 1.0
    t4, f4 = rk1_4(rates, tspan, f0, h, rk)
    h = 2.0
    t41, f41 = rk1_4(rates, tspan, f0, h, rk)

    output(t0, tf, x0, x_dot0, wn, z, w, Fo, m, t1, f1, t11, f11, t2, f2, t21, f21, t3, f3, t31, f31, t4, f4, t41, f41)

def rates(t, f):
    m = 1
    z = 0.03
    wn = 1
    Fo = 1
    w = 0.4 * wn
    
    x = f[0]
    Dx = f[1]
    D2x = Fo/m * np.sin(w * t) - 2 * z * wn * Dx - wn**2 * x
    return np.array([Dx, D2x])

def rk1_4(rates, tspan, f0, h, rk):
    t0, tf = tspan
    t = np.arange(t0, tf, h)
    f = np.zeros((len(t), len(f0)))
    f[0, :] = f0

    for i in range(1, len(t)):
        if rk == 1:
            k1 = rates(t[i-1], f[i-1, :])
            f[i, :] = f[i-1, :] + h * k1
        elif rk == 2:
            k1 = rates(t[i-1], f[i-1, :])
            k2 = rates(t[i-1] + h/2, f[i-1, :] + h/2 * k1)
            f[i, :] = f[i-1, :] + h * k2
        elif rk == 3:
            k1 = rates(t[i-1], f[i-1, :])
            k2 = rates(t[i-1] + h/2, f[i-1, :] + h/2 * k1)
            k3 = rates(t[i-1] + h, f[i-1, :] - h * k1 + 2 * h * k2)
            f[i, :] = f[i-1, :] + h/6 * (k1 + 4 * k2 + k3)
        elif rk == 4:
            k1 = rates(t[i-1], f[i-1, :])
            k2 = rates(t[i-1] + h/2, f[i-1, :] + h/2 * k1)
            k3 = rates(t[i-1] + h/2, f[i-1, :] + h/2 * k2)
            k4 = rates(t[i-1] + h, f[i-1, :] + h * k3)
            f[i, :] = f[i-1, :] + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return t, f

def output(t0, tf, x0, x_dot0, wn, z, w, Fo, m, t1, f1, t11, f11, t2, f2, t21, f21, t3, f3, t31, f31, t4, f4, t41, f41):
    wd = wn * np.sqrt(1 - z**2)
    den = (wn**2 - w**2)**2 + (2 * w * wn * z)**2
    C1 = (wn**2 - w**2) / den * Fo / m
    C2 = -2 * w * wn * z / den * Fo / m
    A = x0 * wn / wd + x_dot0 / wd + (w**2 + (2 * z**2 - 1) * wn**2) / den * w / wd * Fo / m
    B = x0 + 2 * w * wn * z / den * Fo / m

    t = np.linspace(t0, tf, 5000)
    x = (A * np.sin(wd * t) + B * np.cos(wd * t)) * np.exp(-wn * z * t) + C1 * np.sin(w * t) + C2 * np.cos(w * t)

    # Plot solutions
    plt.figure(figsize=(10, 15))

    # Exact
    plt.subplot(5, 1, 1)
    plt.plot(t / max(t), x / max(x), 'k', linewidth=1)
    plt.title('Exact')
    plt.grid(False)
    plt.axis('tight')

    # RK1
    plt.subplot(5, 1, 2)
    plt.plot(t1 / max(t1), f1[:, 0] / max(f1[:, 0]), '-r', linewidth=1)
    plt.plot(t11 / max(t11), f11[:, 0] / max(f11[:, 0]), '-k')
    plt.title('RK1')
    plt.legend(['h = 0.01', 'h = 0.1'])
    plt.grid(False)
    plt.axis('tight')

    # RK2
    plt.subplot(5, 1, 3)
    plt.plot(t2 / max(t2), f2[:, 0] / max(f2[:, 0]), '-r', linewidth=1)
    plt.plot(t21 / max(t21), f21[:, 0] / max(f21[:, 0]), '-k')
    plt.title('RK2')
    plt.legend(['h = 0.1', 'h = 0.5'])
    plt.grid(False)
    plt.axis('tight')

    # RK3
    plt.subplot(5, 1, 4)
    plt.plot(t3 / max(t3), f3[:, 0] / max(f3[:, 0]), '-r', linewidth=1)
    plt.plot(t31 / max(t31), f31[:, 0] / max(f31[:, 0]), '-k')
    plt.title('RK3')
    plt.legend(['h = 0.5', 'h = 1.0'])
    plt.grid(False)
    plt.axis('tight')

    # RK4
    plt.subplot(5, 1, 5)
    plt.plot(t4 / max(t4), f4[:, 0] / max(f4[:, 0]), '-r', linewidth=1)
    plt.plot(t41 / max(t41), f41[:, 0] / max(f41[:, 0]), '-k')
    plt.title('RK4')
    plt.legend(['h = 1.0', 'h = 2.0'])
    plt.grid(False)
    plt.axis('tight')

    plt.show()

if __name__ == "__main__":
    rk1_4_test()