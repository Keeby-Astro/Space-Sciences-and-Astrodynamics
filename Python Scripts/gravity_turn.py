import numpy as np
import matplotlib.pyplot as plt
from atmosisa import atmosisa  # Assuming atmosisa is a valid Python function
from rkf45 import rkf45  # Assuming rkf45 is a valid Python function

def gravity_turn():
    '''
    This program numerically integrates Equations 13.6 through
    13.8 for a gravity turn trajectory.
    
    User py-functions required: rkf45, atmosisa
    User subfunction required:  rates
    '''
    deg = np.pi / 180  # Convert degrees to radians
    g0 = 9.81  # Sea-level acceleration of gravity (m/s)
    Re = 6378.14e3  # Radius of the earth (m)
    hscale = 7.5e3  # Density scale height (m)
    rho0 = 1.225  # Sea level density of atmosphere (kg/m^3)

    diam = 196.85 / 12 * 0.3048  # Vehicle diameter (m)
    A = np.pi / 4 * (diam) ** 2  # Frontal area (m^2)
    CD = 0.5  # Drag coefficient (assumed constant)
    m0 = 149912 * 0.4536  # Lift-off mass (kg)
    n = 7  # Mass ratio
    T2W = 1.4  # Thrust to weight ratio
    Isp = 390  # Specific impulse (s)

    mfinal = m0 / n  # Burnout mass (kg)
    Thrust = T2W * m0 * g0  # Rocket thrust (N)
    m_dot = Thrust / Isp / g0  # Propellant mass flow rate (kg/s)
    mprop = m0 - mfinal  # Propellant mass (kg)
    tburn = mprop / m_dot  # Burn time (s)
    hturn = 130  # Height at which pitchover begins (m)

    t0 = 0  # Initial time for the numerical integration
    tf = tburn  # Final time for the numerical integration
    tspan = [t0, tf]  # Range of integration

    # Initial conditions:
    v0 = 0  # Initial velocity (m/s)
    gamma0 = 89.85 * deg  # Initial flight path angle (rad)
    x0 = 0  # Initial downrange distance (km)
    h0 = 0  # Initial altitude (km)
    vD0 = 0  # Initial velocity loss due to drag (m/s)
    vG0 = 0  # Initial velocity loss due to gravity (m/s)

    # Initial conditions vector:
    f0 = [v0, gamma0, x0, h0, vD0, vG0]

    # Call to Runge-Kutta numerical integrator 'rkf45'
    # rkf45 solves the system of equations df/dt = f(t):
    t, f = rkf45(rates, tspan, f0)

    # Solution f(t) returned on the time interval [t0, tf]:
    v = f[:, 0] * 1.e-3  # Velocity (km/s)
    gamma = f[:, 1] / deg  # Flight path angle (degrees)
    x = f[:, 2] * 1.e-3  # Downrange distance (km)
    h = f[:, 3] * 1.e-3  # Altitude (km)
    vD = -f[:, 4] * 1.e-3  # Velocity loss due to drag (km/s)
    vG = -f[:, 5] * 1.e-3  # Velocity loss due to gravity (km/s)

    # Dynamic pressure vs time:
    q = []
    M = []
    for i in range(len(t)):
        Rho = rho0 * np.exp(-h[i] * 1000 / hscale)  # Air density (kg/m^3)
        q.append(0.5 * Rho * (v[i] * 1.e3) ** 2)  # Dynamic pressure (Pa)
        _, a, _, _ = atmosisa(h[i] * 1000)  # Speed of sound (m/s)
        M.append(1000 * v[i] / a)  # Mach number

    # Maximum dynamic pressure and corresponding parameters:
    maxQ = max(q)  # qMax
    imax = q.index(maxQ)
    tQ = t[imax]  # Time
    vQ = v[imax]  # Speed
    hQ = h[imax]  # Altitude
    _, aQ, _, _ = atmosisa(h[imax] * 1000)  # Speed of sound at altitude
    MQ = 1000 * vQ / aQ

    # Output results:
    print("\n\n -----------------------------------")
    print(f"\n Initial flight path angle = {gamma0 / deg:.3f} deg ")
    print(f"\n Pitchover altitude        = {hturn:.3f} m   ")
    print(f"\n Burn time                 = {tburn:.3f} s   ")
    print(f"\n Maximum dynamic pressure  = {maxQ * 9.869e-6:.3f} atm ")
    print(f"    Time                   = {tQ / 60:.3f} min ")
    print(f"    Speed                  = {vQ:.3f} km/s")
    print(f"    Altitude               = {hQ:.3f} km  ")
    print(f"    Mach Number            = {MQ:.3f}     ")
    print("\n At burnout:")
    print(f"    Speed                  = {v[-1]:.3f} km/s")
    print(f"    Flight path angle      = {gamma[-1]:.3f} deg ")
    print(f"    Altitude               = {h[-1]:.3f} km  ")
    print(f"    Downrange distance     = {x[-1]:.3f} km  ")
    print(f"    Drag loss              = {vD[-1]:.3f} km/s")
    print(f"    Gravity loss           = {vG[-1]:.3f} km/s")
    print("\n\n -----------------------------------")

    # Plot Trajectory and Dynamic Pressure
    plt.figure(figsize=(10, 8))

    # Subplot 1: Altitude vs Downrange Distance
    plt.subplot(2, 1, 1)
    plt.plot(x, h, label='Trajectory')
    plt.title('(a) Altitude vs Downrange Distance')
    plt.xlabel('Downrange Distance (km)')
    plt.ylabel('Altitude (km)')
    plt.axis('equal')
    plt.grid(True)

    # Subplot 2: Dynamic Pressure vs Altitude
    plt.subplot(2, 1, 2)
    plt.plot(h, np.array(q) * 9.869e-6, label='Dynamic Pressure', color='r')
    plt.title('(b) Dynamic Pressure vs Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Dynamic Pressure (atm)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# Define rates function
def rates(t, y):
    '''
    Calculates the time derivatives dy/dt for the gravity turn trajectory.
    '''
    # Constants
    deg = np.pi / 180  # Convert degrees to radians
    g0 = 9.81  # Sea-level acceleration of gravity (m/s)
    Re = 6378.14e3  # Radius of the earth (m)
    hscale = 7.5e3  # Density scale height (m)
    rho0 = 1.225  # Sea level density of atmosphere (kg/m^3)

    diam = 196.85 / 12 * 0.3048  # Vehicle diameter (m)
    A = np.pi / 4 * (diam) ** 2  # Frontal area (m^2)
    CD = 0.5
    m0 = 149912 * 0.4536  # Lift-off mass (kg)
    n = 7  # Mass ratio
    T2W = 1.4  # Thrust to weight ratio
    Isp = 390  # Specific impulse (s)

    mfinal = m0 / n  # Burnout mass (kg)
    Thrust = T2W * m0 * g0  # Rocket thrust (N)
    m_dot = Thrust / Isp / g0  # Propellant mass flow rate (kg/s)
    mprop = m0 - mfinal  # Propellant mass (kg)
    tburn = mprop / m_dot  # Burn time (s)
    hturn = 130  # Height at which pitchover begins (m)

    dydt = np.zeros(6)

    v, gamma, x, h, vD, vG = y

    if t < tburn:
        m = m0 - m_dot * t  # Current vehicle mass
        T = Thrust  # Current thrust
    else:
        m = m0 - m_dot * tburn  # Current vehicle mass
        T = 0  # Current thrust

    g = g0 / (1 + h / Re) ** 2  # Gravitational variation with altitude h
    rho = rho0 * np.exp(-h / hscale)  # Exponential density variation with altitude
    D = 0.5 * rho * v ** 2 * A * CD  # Drag

    if h <= hturn:
        gamma_dot = 0
        v_dot = T / m - D / m - g
        x_dot = 0
        h_dot = v
        vG_dot = -g
    else:
        v_dot = T / m - D / m - g * np.sin(gamma)
        gamma_dot = -1 / v * (g - v ** 2 / (Re + h)) * np.cos(gamma)
        x_dot = Re / (Re + h) * v * np.cos(gamma)
        h_dot = v * np.sin(gamma)
        vG_dot = -g * np.sin(gamma)

    vD_dot = -D / m

    dydt[0] = v_dot
    dydt[1] = gamma_dot
    dydt[2] = x_dot
    dydt[3] = h_dot
    dydt[4] = vD_dot
    dydt[5] = vG_dot

    return dydt

if __name__ == "__main__":
    gravity_turn()
