import numpy as np
import matplotlib.pyplot as plt
import rkf45

def rkf45_test_lunar():
    '''
    This program uses the Runge-Kutta-Fehlberg 4(5) method to solve the
    earth-moon restricted three-body problem (Equations 2.192a and 2.192b)
    for the trajectory of a spacecraft.

    The numerical integration is done in the external function 'rkf45',
    which uses the subfunction 'rates' herein to compute the derivatives.

    days      - converts days to seconds
    G         - universal gravitational constant (km^3/kg/s^2)
    rmoon     - radius of the moon (km)
    rearth    - radius of the earth (km)
    r12       - distance from center of earth to center of moon (km)
    m1,m2     - masses of the earth and of the moon, respectively (kg)
    M         - total mass of the restricted 3-body system (kg)
    mu        - gravitational parameter of earth-moon system (km^3/s^2)
    mu1,mu2   - gravitational parameters of the earth and of the moon,
                respectively (km^3/s^2)
    pi_1,pi_2 - ratios of the earth mass and the moon mass, respectively,
                to the total earth-moon mass
    W         - angular velocity of moon around the earth (rad/s)
    x1,x2     - x-coordinates of the earth and of the moon, respectively,
                relative to the earth-moon barycenter (km)
    d0        - initial altitude of spacecraft (km)
    phi       - polar azimuth coordinate (degrees) of the spacecraft
                measured positive counterclockwise from the earth-moon line
    v0        - initial speed of spacecraft relative to rotating earth-moon
                system (km/s)
    gamma     - initial flight path angle (degrees)
    r0        - initial radial distance of spacecraft from the earth (km)
    x,y       - x and y coordinates of spacecraft in rotating earth-moon
                system (km)
    vx,vy     - x and y components of spacecraft velocity relative to
                rotating earth-moon system (km/s)
    f0        - column vector containing the initial values of x, y, vx and vy
    t0,tf     - initial time and final times (s)
    t         - column vector of times at which the solution was computed
    f         - a matrix whose columns are:
                column 1: solution for x at the times in t
                column 2: solution for y at the times in t
                column 3: solution for vx at the times in t
                column 4: solution for vy at the times in t
    xf,yf     - x and y coordinates of spacecraft in rotating earth-moon
                system at tf
    vxf, vyf  - x and y components of spacecraft velocity relative to
                rotating earth-moon system at tf
    df        - distance from surface of the moon at tf
    vf        - relative speed at tf

    User py-functions required: rkf45
    User subfunctions required: rates, circle
    '''
    # Constants and input data:
    days   = 24 * 3600
    G      = 6.6742e-20
    rmoon  = 1737
    rearth = 6378
    r12    = 384400
    m1     = 5.974e24
    m2     = 7.348e22

    M    = m1 + m2
    pi_1 = m1 / M
    pi_2 = m2 / M

    mu1 = 398600
    mu2 = 4903.02
    mu  = mu1 + mu2

    W  = np.sqrt(mu / r12**3)
    x1 = -pi_2 * r12
    x2 = pi_1 * r12

    #...Input data:
    d0    = 200
    phi   = -90
    v0    = 10.9148
    gamma = 20
    t0    = 0
    tf    = 3.16689 * days

    r0 = rearth + d0
    x  = r0 * np.cos(np.radians(phi)) + x1
    y  = r0 * np.sin(np.radians(phi))

    vx = v0 * (np.sin(np.radians(gamma)) * np.cos(np.radians(phi)) - np.cos(np.radians(gamma)) * np.sin(np.radians(phi)))
    vy = v0 * (np.sin(np.radians(gamma)) * np.sin(np.radians(phi)) + np.cos(np.radians(gamma)) * np.cos(np.radians(phi)))
    f0 = [x, y, vx, vy]

    #...Compute the trajectory:
    t, f = rkf45.rkf45(rates, [t0, tf], f0)
    x    = f[:, 0]
    y    = f[:, 1]
    vx   = f[:, 2]
    vy   = f[:, 3]

    xf   = x[-1]
    yf   = y[-1]

    vxf  = vx[-1]
    vyf  = vy[-1]

    df   = np.linalg.norm([xf - x2, yf - 0]) - rmoon
    vf   = np.linalg.norm([vxf, vyf])

    #...Output the results:
    output(d0, phi, gamma, tf, days, df, vf, x, y, x1, x2, rearth, rmoon)

# --------------
def rates(t, f):
# --------------
    '''
    This subfunction calculates the components of the relative acceleration
    for the restricted 3-body problem, using Equations 2.192a and 2.192b

    ax,ay - x and y components of relative acceleration (km/s^2)
    r1    - spacecraft distance from the earth (km)
    r2    - spacecraft distance from the moon (km)
    f     - column vector containing x, y, vx and vy at time t
    fdt   - column vector containing vx, vy, ax and ay at time t

    All other variables are defined above.

    User py-functions required: none
    '''
    m1 = 5.974e24
    m2 = 7.348e22

    M    = m1 + m2
    pi_1 = m1 / M
    pi_2 = m2 / M

    mu1 = 398600
    mu2 = 4903.02
    mu  = mu1 + mu2

    W  = np.sqrt(mu / 384400**3)
    x1 = -pi_2 * 384400
    x2 = pi_1 * 384400
    
    r12 = 384400

    x  = f[0]
    y  = f[1]
    vx = f[2]
    vy = f[3]

    r1 = np.linalg.norm([x + pi_2 * r12, y])
    r2 = np.linalg.norm([x - pi_1 * r12, y])

    ax = 2 * W * vy + W**2 * x - mu1 * (x - x1) / r1**3 - mu2 * (x - x2) / r2**3
    ay = -2 * W * vx + W**2 * y - (mu1 / r1**3 + mu2 / r2**3) * y

    dfdt = np.array([vx, vy, ax, ay])

    return dfdt

def output(d0, phi, gamma, tf, days, df, vf, x, y, x1, x2, rearth, rmoon):
    '''
    This subfunction echoes the input data and prints the results to the
    console. It also plots the trajectory.

    User py-functions required: none
    User subfunction required: circle
    '''
    print('-----------------------------------------------------------')
    print('\n Example 2.18: Lunar trajectory using the restricted')
    print(' three-body equations.')
    print(f'\n Initial Earth altitude (km)         = {d0}')
    print(f' Initial angle between radial')
    print(f' and earth-moon line (degrees)       = {phi}')
    print(f' Initial flight path angle (degrees) = {gamma}')
    print(f' Flight time (days)                  = {tf / days}')
    print(f' Final distance from the moon (km)   = {df}')
    print(f' Final relative speed (km/s)         = {vf}')
    print('-----------------------------------------------------------\n')

    #...Plot the trajectory and place filled circles representing the earth
    #   and moon on the plot:
    plt.plot(x, y)
    # Set plot display parameters
    xmin = -20.e3
    xmax = 4.e5
    ymin = -20.e3
    ymax = 1.e5
    plt.axis([xmin, xmax, ymin, ymax])
    plt.axis('equal')
    plt.xlabel('x, km')
    plt.ylabel('y, km')
    plt.title('Trajectory of the spacecraft')

    #...Plot the earth (blue) and moon (green) to scale
    earth = circle(x1, 0, rearth)
    moon  = circle(x2, 0, rmoon)
    plt.fill(earth[:, 0], earth[:, 1], 'b')
    plt.fill(moon[:, 0], moon[:, 1], 'grey')
    plt.show()

# -------------------------
def circle(xc, yc, radius):
# -------------------------
    '''
    This subfunction calculates the coordinates of points spaced
    0.1 degree apart around the circumference of a circle

    x,y    - x and y coordinates of a point on the circumference
    xc,yc  - x and y coordinates of the center of the circle
    radius - radius of the circle
    xy     - an array containing the x coordinates in column 1 and the
             y coordinates in column 2

    User py-functions required: none
    '''
    theta = np.radians(np.arange(0, 360.1, 0.1))
    x     = xc + radius * np.cos(theta)
    y     = yc + radius * np.sin(theta)
    xy    = np.column_stack((x, y))

    return xy

if __name__ == "__main__":
    rkf45_test_lunar()
