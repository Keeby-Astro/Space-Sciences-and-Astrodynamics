import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from rsmooth import rsmooth
import sv_from_coe
import solar_position

def solar_perturbation():
    '''
    This function uses Python's solve_ivp to integrate
    the Gauss variational equations for a solar gravitational perturbation.
    
    User modules required: sv_from_coe, solar_position
    '''
    global JD  # Julian day

    #...Conversion factors:
    hours = 3600      # Hours to seconds
    days = 24 * hours # Days to seconds
    deg = np.pi / 180 # Degrees to radians

    #...Constants
    mu = 398600       # Earth's gravitational parameter (km^3/s^2)
    mu3 = 132.712e9   # Sun's gravitational parameter (km^3/s^2)
    RE = 6378         # Earth's radius (km)

    #...Initial data for each of the three given orbits:
    orbit_types = ['GEO', 'HEO', 'LEO']

    #...GEO
    n = 1
    a0 = 42164     # Semimajor axis (km)
    e0 = 0.0001    # Eccentricity
    w0 = 0         # Argument of perigee (rad)
    RA0 = 0        # Right ascension (rad)
    i0 = 1 * deg   # Inclination (rad)
    TA0 = 0        # True anomaly (rad)
    JD0 = 2454283  # Julian Day
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n, mu, mu3, days, deg)

    #...HEO
    n = 2
    a0 = 26553.4
    e0 = 0.741
    w0 = 270 * deg
    RA0 = 0
    i0 = 63.4 * deg
    TA0 = 0
    JD0 = 2454283
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n, mu, mu3, days, deg)

    #...LEO
    n = 3
    a0 = 6678.136
    e0 = 0.01
    w0 = 0
    RA0 = 0
    i0 = 28.5 * deg
    TA0 = 0
    JD0 = 2454283
    solveit(a0, e0, w0, RA0, i0, TA0, JD0, n, mu, mu3, days, deg)

#...Subfunctions:

def solveit(a0, e0, w0, RA0, i0, TA0, JD0, n, mu, mu3, days, deg):
    '''
    Calculations and plots common to all of the orbits
    '''

    #...Initial orbital parameters (calculated from the given data):
    h0 = np.sqrt(mu * a0 * (1 - e0 ** 2)) # Angular momentum (km^2/s)
    coe0 = [h0, e0, RA0, i0, w0, TA0]     # Initial orbital elements

    #...Use the ODE to integrate the Equations 12.84, the Gauss variational
    #   equations with lunar gravity as the perturbation, from t0 to tf:
    t0 = 0
    tf = 720 * days
    nout = 400
    tspan = np.linspace(t0, tf, nout)

    sol = solve_ivp(
        rates,
        [t0, tf],
        coe0,
        t_eval=tspan,
        args=(mu, mu3, JD0, days),
        rtol=1e-8,
        atol=1e-8
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    #...Time histories of the right ascension, inclination and argument of
    #   perigee:
    t = sol.t
    y = sol.y
    RA, i, w = y[2], y[3], y[4]

    #...Smooth the data to eliminate short period variations:
    RA = rsmooth(RA)
    i = rsmooth(i)
    w = rsmooth(w)

    plt.figure(n)

    plt.subplot(1, 3, 1)
    plt.plot(t / days, (RA - RA0) / deg)
    plt.title("Right Ascension vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("RA (deg)")
    plt.tight_layout()

    plt.subplot(1, 3, 2)
    plt.plot(t / days, (i - i0) / deg)
    plt.title("Inclination vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("Inclination (deg)")
    plt.tight_layout()

    plt.subplot(1, 3, 3)
    plt.plot(t / days, (w - w0) / deg)
    plt.title("Argument of Perigee vs Time")
    plt.xlabel("Time (days)")
    plt.ylabel("Argument of Perigee (deg)")
    plt.tight_layout()

    plt.show()

def rates(t, f, mu, mu3, JD0, days):
    '''
    Compute the rates of change of the orbital elements.
    '''
    #...The orbital elements at time t:
    h, e, RA, i, w, TA = f
    phi = w + TA  # Argument of latitude

    #...Obtain the state vector at time t from Algorithm 4.5:
    coe = [h, e, RA, i, w, TA]
    R, V = sv_from_coe.sv_from_coe(coe, mu)

    #...Update Julian Day
    JD = JD0 + t / days

    #...Find and normalize the position vector of the sun:
    lamda, eps, R_S = solar_position.solar_position(JD)
    r_S = np.linalg.norm(R_S)

    R_rel = R_S - R # Position vector of sun with respect to satellite
    r_rel = np.linalg.norm(R_rel)

    #...See Appendix F
    q = np.dot(R, (2 * R_S - R)) / r_S ** 2
    F = (q ** 2 - 3 * q + 3) * q / (1 + (1 - q) ** 1.5)

    #...Gravitationl perturbation of the sun (Equation 12.130):
    ap = mu3 / r_rel ** 3 * (F * R_S - R)

    #...Perturbation components in the rsw system:
    apr = np.dot(ap, R / np.linalg.norm(R))
    aps = np.dot(ap, np.cross(np.cross(R, V), R) / np.linalg.norm(np.cross(np.cross(R, V), R)))
    aph = np.dot(ap, np.cross(R, V) / np.linalg.norm(np.cross(R, V)))

    #...Gauss variational equations (Equations 12.84):
    hdot = np.linalg.norm(R) * aps
    edot = (h / mu * np.sin(TA) * apr + (1 / mu / h) * ((h ** 2 + mu * np.linalg.norm(R)) * np.cos(TA) + mu * e * np.linalg.norm(R)) * aps)
    RAdot = (np.linalg.norm(R) / h / np.sin(i) * np.sin(phi) * aph)
    idot = (np.linalg.norm(R) / h * np.cos(phi) * aph)
    wdot = (-h * np.cos(TA) / mu / e * apr + (h ** 2 + mu * np.linalg.norm(R)) / mu / e / h * np.sin(TA) * aps - np.linalg.norm(R) * np.sin(phi) / h / np.tan(i) * aph)
    TAdot = (h / np.linalg.norm(R) ** 2 + 1 / e / h * (h ** 2 / mu * np.cos(TA) * apr - (np.linalg.norm(R) + h ** 2 / mu) * np.sin(TA) * aps))

    #...Return rates to solve_ivp in the array dfdt:
    dfdt = np.array([hdot, edot, RAdot, idot, wdot, TAdot])

    return dfdt

if __name__ == '__main__':
    solar_perturbation()