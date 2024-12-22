import numpy as np

def f_and_g_ta(r0, v0, dt, mu):
    '''
    This function calculates the Lagrange f and g coefficients from the
    change in true anomaly since time t0.

    mu  - gravitational parameter (km^3/s^2)
    dt  - change in true anomaly (degrees)
    r0  - position vector at time t0 (km)
    v0  - velocity vector at time t0 (km/s)
    h   - angular momentum (km^2/s)
    vr0 - radial component of v0 (km/s)
    r   - radial position after the change in true anomaly
    f   - the Lagrange f coefficient (dimensionless)
    g   - the Lagrange g coefficient (s)

    User py-functions required: None
    '''
    h = np.linalg.norm(np.cross(r0, v0))
    vr0 = np.dot(v0, r0) / np.linalg.norm(r0)
    r0_norm = np.linalg.norm(r0)
    s = np.sin(np.radians(dt))
    c = np.cos(np.radians(dt))

    # Equation 2.152:
    r = h**2 / mu / (1 + (h**2 / mu / r0_norm - 1) * c - h * vr0 * s / mu)

    # Equations 2.158a & b:
    f = 1 - mu * r * (1 - c) / h**2
    g = r * r0_norm * s / h

    return f, g