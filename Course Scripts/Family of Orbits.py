import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation

# Gravitational parameter
mu = 3.986004418e5 # km^3/s^2

# Equatorial radius of the Earth
R_e = 6378.14 # km

numDataPoints = 100
theta_0 = 0
orbit = np.zeros((numDataPoints + 1, 2))

def f(a, e):
    for i in range(0, numDataPoints + 1):
        theta = theta_0 + float(i)*2*np.pi/numDataPoints
    
        r = a*(1 - e**2)/(1 + e*np.cos(theta))
    
        x = r*np.cos(theta)
        y = r*np.sin(theta)
    
        orbit[i, 0] = x
        orbit[i, 1] = y

    return orbit

###### Animated plot

a = 40000.0

def animate_func(num):
    ax.clear()
    
    # orbit1
    e1 = 0
    orbit1 = f(a, e1)
    ax.plot(orbit1[:num+1,0], orbit1[:num+1,1], color='black')
    ax.scatter(orbit1[num, 0], orbit1[num, 1], marker='o', color='black')

    # orbit2
    e2 = 0.2
    orbit2 = f(a, e2)
    ax.plot(orbit2[:num+1,0], orbit2[:num+1,1], color='purple')
    ax.scatter(orbit2[num, 0], orbit2[num, 1], marker='o', color='purple')
    
    # orbit3
    e3 = 0.4
    orbit3 = f(a, e3)
    ax.plot(orbit3[:num+1,0], orbit3[:num+1,1], color='blue')
    ax.scatter(orbit3[num, 0], orbit3[num, 1], marker='o', color='blue')
    
    # orbit4
    e4 = 0.6
    orbit4 = f(a, e4)
    ax.plot(orbit4[:num+1,0], orbit4[:num+1,1], color='green')
    ax.scatter(orbit4[num, 0], orbit4[num, 1], marker='o', color='green')
    
    # orbit5
    e5 = 0.8
    orbit5 = f(a, e5)
    ax.plot(orbit5[:num+1,0], orbit5[:num+1,1], color='red')
    ax.scatter(orbit5[num, 0], orbit5[num, 1], marker='o', color='red')

    
    # Earth
    ax.scatter(0, 0, marker='o', color='blue', s = 1000)
    
    ax.set_xlabel('X, km', fontsize=12)
    ax.set_ylabel('Y, km', fontsize=12)

    ax.set_xlim(-80000,50000)
    ax.set_ylim(-42000,42000)
    
    ax.set_aspect('equal')

    ax.set_title('a = 40,000 km', fontsize=18)
    
    plt.tight_layout()

if __name__ == "__main__":
    # Plotting the Animation
    fig = plt.figure()
    ax = plt.axes()
    animation_1 = animation.FuncAnimation(fig, animate_func, interval=1, frames=numDataPoints)
    plt.show()
