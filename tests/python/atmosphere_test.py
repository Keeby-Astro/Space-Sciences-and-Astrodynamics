import numpy as np
import matplotlib.pyplot as plt
import atmosphere

def atmosphere_test():
    '''
    Tests the atmosphere function by plotting the density
    as a function of altitude from sea level through 1000 km.
    '''
    #...Geometric altitudes (km):
    z = np.linspace(0, 1000, 1000)

    #...Calculate densities:
    density = np.array([atmosphere.atmosphere(zi) for zi in z])

    #...Plot the results:
    plt.figure()
    plt.semilogy(z, density, linewidth=1.5)
    plt.grid(True)
    plt.title('Atmospheric Density vs. Altitude')
    plt.xlabel('Altitude (km)')
    plt.ylabel('Density (kg/m^3)')
    plt.legend(['Density'], loc='upper right')
    plt.show()

if __name__ == '__main__':
    atmosphere_test()
