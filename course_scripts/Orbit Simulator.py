import numpy as np
import timeit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import os
import csv

# Functions with orbital elements <-> state vector conversion
from functions import orb2xyz
from functions import xyz2orb

# Constants
mu = 398600.4418  # km^3/s^2
Re = 6378.14  # km

# -------------------### Spacecraft's initial conditions ###---------------------
a = Re + 50000.0  # Semi-major axis (km)
e = 0.7  # Eccentricity
i = np.deg2rad(5)  # Inclination w.r.t. Earth's equatorial plane (rad)
w = np.deg2rad(27)  # Argument of the perigee (rad)
OM = np.deg2rad(0.0)  # Longitude of the ascending node (rad)
theta = np.deg2rad(0.0)  # True anomaly (rad)

print("\n--------------### Spacecraft's initial conditions ###----------------")
print(f"a = {a:.5e} km       e = {e:.5e}      i = {np.rad2deg(i):.3f} deg")
print(f"w = {np.rad2deg(w):.2f} deg            OM = {np.rad2deg(OM):.2f} deg     theta = {np.rad2deg(theta):.2f} deg\n")

# Store the orbital elements in a list
oe = [a, e, i, w, OM, theta]

# Orbital elements to Cartesian coordinates (state vector)
init_sc = orb2xyz(mu, oe)

initial_state = init_sc

# -------------#### Event Functions to control the integration ####--------------
critical = 0

def event1(t, f):
    return np.linalg.norm(f[:3]) - Re

event1.direction = -1
event1.terminal = True

# -----------------#### Function with derivatives (model) ####-------------------
def Derivs(t, f):
    GM = mu
    x_sat = np.zeros(6)
    x_sat[0] = f[0]
    x_sat[1] = f[1]
    x_sat[2] = f[2]
    x_sat[3] = f[3]
    x_sat[4] = f[4]
    x_sat[5] = f[5]
    r = np.sqrt(x_sat[0] ** 2 + x_sat[1] ** 2 + x_sat[2] ** 2)

    # Spacecraft EOM
    dxdt = x_sat[3]
    dydt = x_sat[4]
    dzdt = x_sat[5]

    ddxdt = - GM * x_sat[0] / r ** 3
    ddydt = - GM * x_sat[1] / r ** 3
    ddzdt = - GM * x_sat[2] / r ** 3

    return [dxdt, dydt, dzdt, ddxdt, ddydt, ddzdt]

# -------------------------------------------------------------------------------
# Integration process
# inital time t0 = 0 seconds
tf = 20 * 24 * 3600  # final time in seconds
dt = 1.0  # time step
nt = int(tf / dt)  # number of steps
tspan = np.arange(0.0, tf + dt, dt)

start = timeit.default_timer()
print('\nRunning...\n')

solution = solve_ivp(Derivs, [0.0, tf], initial_state, events=[event1], method='DOP853', \
                     t_eval=tspan, rtol=1e-8, atol=1e-8)

if not solution.success:
    raise RuntimeError(solution.message)

state = solution.y
times = solution.t

# Check if state has data
if state.size == 0:
    print("No data in state variable. Check the solver's output.")
else:
    print("State variable contains data.")

x = state[0, :]
y = state[1, :]
z = state[2, :]

# Print the range of data for debugging
print(f"x range: {min(x)} to {max(x)}")
print(f"y range: {min(y)} to {max(y)}")
print(f"z range: {min(z)} to {max(z)}")

stop = timeit.default_timer()
runtime = stop - start
if runtime < 60.0:
    print(f"Runtime = {runtime:.2f} seconds.\n")
elif runtime >= 60.0 and runtime < 3600.0:
    print(f"Runtime = {runtime / 60.0:.2f} minutes.\n")
else:
    print(f"Runtime = {runtime / 3600.0:.2f} hours.\n")

################################################################################
# Data analysis
plt.plot(state[0, :], state[1, :])
# plt.plot(state[6, :], state[7, :])
plt.show()

# Save trajectories in file
file_name_1 = "course-state.dat"
file_name_2 = "course-orbital.csv"

isExist = os.path.exists(file_name_1)
if isExist:
    print(file_name_1, "already exists. Overwriting.")
    os.remove(file_name_1)
with open(file_name_1, "a") as file1:
    np.savetxt(file1, np.transpose(state))
file1.close()

# Spacecraft: state vector to orbital elements and CSV output
orb = []
with open(file_name_2, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Time (sec)', 'a (km)', 'e', 'i (deg)', 'w (deg)', 'OM (deg)', 'theta (deg)'])
    for j in range(len(times)):
        elements = xyz2orb(mu, state[0:3, j], state[3:6, j])
        orb.append(elements)
        csvwriter.writerow([times[j]] + elements)

orbit = np.array(orb)

# Creating a 3D plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Plotting the orbit
ax.plot(x, y, z, label='Spacecraft Orbit')

# Setting labels and title
ax.set_xlabel('X Position (km)')
ax.set_ylabel('Y Position (km)')
ax.set_zlabel('Z Position (km)')
ax.set_title('3D Plot of Spacecraft Orbit')

# Setting the axes limits (adjust these based on your data range)
ax.set_xlim([min(x), max(x)])
ax.set_ylim([min(y), max(y)])
ax.set_zlim([min(z), max(z)])

# Adding a grid and legend
ax.grid(True)
ax.legend()

# Show the plot
plt.show()