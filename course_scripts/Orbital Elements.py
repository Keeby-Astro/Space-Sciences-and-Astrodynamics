# Orbital Elements and the state vectors
import numpy as np

mu = 398600.4418 # km^3/s^2

pi = np.pi
pi2 = 2*np.pi

r_vec = np.array([-6045, -3490, 2500])
print(f"r_vec = {r_vec} (km)")
v_vec = np.array([-3.457, 6.618, 2.533]) 
print(f"v_vec = {v_vec} (km/s)")
print("")

r = np.linalg.norm(r_vec) # magnitude of the position vector
print(f"r = {r} km")

v = np.linalg.norm(v_vec) # magnitude of the velocity vector
print(f"v = {v} km/s")

if r == 0:
    raise ValueError("Position must be nonzero.")

vr = np.dot(r_vec, v_vec)/r # radial velocity
print(f"vr = {vr} km/s")
print("")

h_vec = np.cross(r_vec, v_vec) # vector angular momentum
print(f"h_vec = {h_vec} (km^2/s)")

h = np.linalg.norm(h_vec) # magnitude of the angular momentum
print(f"h = {h} km^2/s")
print("")

if h == 0:
    raise ValueError("Classical orbital elements require nonzero angular momentum.")

i = np.arctan2(np.hypot(h_vec[0], h_vec[1]), h_vec[2]) # inclination
print(f"i = {np.rad2deg(i)} deg")
print("")

# Defining the direction k and the Nodal vector
k_dir = np.array([0.0, 0.0, 1.0])

Nodal_vec = np.cross(k_dir, h_vec)
print(f"Nodal_vec = {Nodal_vec} (km^2/s)")

Nodal = np.linalg.norm(Nodal_vec)
print(f"Nodal = {Nodal} km^2/s")
print("")

# Longitude of the ascending node

if Nodal > 1e-10*h:
    Omega = np.mod(np.arctan2(Nodal_vec[1], Nodal_vec[0]), pi2)
else:
    Omega = 0.0  # Equatorial reference direction

print(f"\u03A9 = {np.rad2deg(Omega)} deg")
print("")

# Eccentricity vector
e_vec = np.cross(v_vec, h_vec)/mu - r_vec/r
print(f"e_vec = {e_vec}")

# Eccentricity
e = np.sqrt(e_vec[0]**2 + e_vec[1]**2 + e_vec[2]**2)
print(f"e = {e}")
print("")

# Argument of the perigee
if e > 1e-10:
    if Nodal > 1e-10*h:
        w = np.mod(np.arctan2(np.dot(np.cross(Nodal_vec,e_vec),h_vec)/h, np.dot(Nodal_vec,e_vec)), pi2)
    else:
        w = np.mod(np.arctan2(np.sign(h_vec[2])*e_vec[1], e_vec[0]), pi2)
else:
    w = 0.0  # Circular reference direction
print(f"\u03C9 = {np.rad2deg(w)} deg")

# True anomaly, argument of latitude, or true longitude for a circular orbit
if e > 1e-10:
    theta = np.mod(np.arctan2(np.dot(np.cross(e_vec,r_vec),h_vec)/h, np.dot(e_vec,r_vec)), pi2)
elif Nodal > 1e-10*h:
    theta = np.mod(np.arctan2(np.dot(np.cross(Nodal_vec,r_vec),h_vec)/h, np.dot(Nodal_vec,r_vec)), pi2)
else:
    theta = np.mod(np.arctan2(np.sign(h_vec[2])*r_vec[1], r_vec[0]), pi2)
print(f"\u03B8 = {np.rad2deg(theta)} deg")

# Semi-major axis
den = 1 - e**2
a = np.inf if abs(den) <= 32*np.finfo(float).eps else h**2/(mu*den)
print(f"a = {a} km")

if den > 32*np.finfo(float).eps:
    # Elliptic eccentric anomaly and mean anomaly
    Ea = np.mod(2*np.arctan2(np.sqrt(1-e)*np.sin(theta/2), np.sqrt(1+e)*np.cos(theta/2)), pi2)
    Me = Ea - e*np.sin(Ea)
    print(f"Ea = {np.rad2deg(Ea)} deg")
    print(f"Me = {np.rad2deg(Me)} deg")
elif den < -32*np.finfo(float).eps:
    # Hyperbolic eccentric anomaly and mean anomaly
    Ha = np.arcsinh(np.sqrt(e**2-1)*np.sin(theta)/(1+e*np.cos(theta)))
    Mh = e*np.sinh(Ha)-Ha
    print(f"Ha = {Ha} rad")
    print(f"Mh = {Mh} rad")
else:
    # Barker parameter for parabolic time of flight
    D = np.tan(theta/2)
    print(f"Barker parameter = {D + D**3/3}")