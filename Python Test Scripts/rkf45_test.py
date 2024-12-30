import numpy as np
import matplotlib.pyplot as plt
import rkf45

def rkf45_test():
    '''
    This program uses RKF4(5) with adaptive step size control
    to solve the differential equation

    x'' + mu/x^2 = 0

    The numerical integration is done by the function 'rkf45' which uses
    the subfunction 'rates' herein to compute the derivatives.
    
    x     - displacement (km)
    '     - shorthand for d/dt
    t     - time (s)
    mu    - = go*RE^2 (km^3/s^2), where go is the sea level gravitational
            acceleration and RE is the radius of the earth
    x0    - initial value of x
    v0    = initial value of the velocity (x')
    y0    - column vector containing x0 and v0
    t0    - initial time
    tf    - final time
    tspan - a row vector with components t0 and tf
    t     - column vector of the times at which the solution is found
    f     - a matrix whose columns are:
            column 1: solution for x at the times in t
            column 2: solution for x' at the times in t

    User py-function required: rkf45
    User py-subfunction required: rates
    '''
    # Parameters
    mu = 398600
    minutes = 60  # Conversion from minutes to seconds

    # Initial conditions
    x0 = 6500
    v0 = 7.8
    y0 = np.array([x0, v0])
    t0 = 0
    tf = 70 * minutes
    tspan = [t0, tf]

    # Solve the differential equation using rkf45
    t, f = rkf45.rkf45(rates, tspan, y0)

    # Plot results
    plotit(t, f, minutes)

def rates(t, f):
    '''
    This function calculates first and second time derivatives of x
    governed by the equation of two-body rectilinear motion.

    x'' + mu/x^2 = 0

    Dx   - velocity x'
    D2x  - acceleration x''
    f    - column vector containing x and Dx at time t
    dfdt - column vector containing Dx and D2x at time t
    
    User py-functions required: none
    '''
    mu = 398600  # Earth's gravitational parameter in km^3/s^2
    x = f[0]
    Dx = f[1]
    D2x = -mu / x**2
    return np.array([Dx, D2x])

def plotit(t, f, minutes):
    # Position vs time
    plt.figure(figsize=(8, 10))

    plt.subplot(2, 1, 1)
    plt.plot(t / minutes, f[:, 0], '-ok')
    plt.xlabel('time, minutes')
    plt.ylabel('position, km')
    plt.grid(True)
    plt.axis([t[0] / minutes, t[-1] / minutes, 5000, 15000])

    # Velocity vs time
    plt.subplot(2, 1, 2)
    plt.plot(t / minutes, f[:, 1], '-ok')
    plt.xlabel('time, minutes')
    plt.ylabel('velocity, km/s')
    plt.grid(True)
    plt.axis([t[0] / minutes, t[-1] / minutes, -10, 10])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    rkf45_test()
