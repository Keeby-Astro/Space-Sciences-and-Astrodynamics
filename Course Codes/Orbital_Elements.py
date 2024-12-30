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

vr = np.dot(r_vec, v_vec)/r # radial velocity
print(f"vr = {vr} km/s")
print("")

h_vec = np.cross(r_vec, v_vec) # vector angular momentum
print(f"h_vec = {h_vec} (km^2/s)")

h = np.linalg.norm(h_vec) # magnitude of the angular momentum
print(f"h = {h} km^2/s")
print("")

i = np.arccos(h_vec[2]/h) # inclination
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

if Nodal_vec[1] >= 0.0:
    Omega = np.arccos(Nodal_vec[0]/Nodal)
else:
    Omega = pi2 - np.arccos(Nodal_vec[0]/Nodal)
    
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
if e_vec[2] >= 0.0:
    w = np.arccos(np.dot(Nodal_vec,e_vec)/(Nodal*e))
else:
    w = pi2 - np.arccos(np.dot(Nodal_vec,e_vec)/(Nodal*e))
print(f"\u03C9 = {np.rad2deg(w)} deg")

# True anomaly
if np.dot(r_vec,v_vec) >= 0.0:
    theta = np.arccos(np.dot(e_vec,r_vec)/(e*r))
else:
    theta = pi2 - np.arccos(np.dot(e_vec,r_vec)/(e*r))
print(f"\u03B8 = {np.rad2deg(theta)} deg")

# Semi-major axis
a = h**2/(mu*(1 - e**2))
print(f"a = {a} km")

# Eccentric anomaly
Ea = np.arccos((e + np.cos(theta))/(1 + e*np.cos(theta)))
print(f"Ea = {np.rad2deg(Ea)} deg")

# Mean anomaly
Me = Ea - e*np.sin(Ea)
print(f"Me = {np.rad2deg(Me)} deg")
