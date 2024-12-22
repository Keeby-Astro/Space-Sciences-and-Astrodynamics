import numpy as np
import rv_from_r0v0

def rv_from_r0v0_test():
    '''
    This program computes the state vector (R,V) from the initial
    state vector (R0,V0) and the elapsed time using the data.
    
    mu - gravitational parameter (km^3/s^2)
    R0 - the initial position vector (km)
    V0 - the initial velocity vector (km/s)
    R  - the final position vector (km)
    V  - the final velocity vector (km/s)
    t  - elapsed time (s)
    
    User py-functions required: rv_from_r0v0
    '''
    #...Data declaration:
    mu = 398600.4418  
    R0 = np.array([7000, -12124, 0]) 
    V0 = np.array([2.6679, 4.6210, 0])
    t = 3600 

    #...Algorithm 3.4:
    R, V = rv_from_r0v0.rv_from_r0v0(R0, V0, t, mu)

    #...Echo the input data and output the results to the console:
    print('-----------------------------------------------------')
    print('\n Initial position vector (km):')
    print(f' r0 = ({R0[0]}, {R0[1]}, {R0[2]})\n')
    print('\n Initial velocity vector (km/s):')
    print(f' v0 = ({V0[0]}, {V0[1]}, {V0[2]})')
    print(f'\n\n Elapsed time = {t} s\n')
    print('\n Final position vector (km):')
    print(f' r = ({R[0]}, {R[1]}, {R[2]})\n')
    print('\n Final velocity vector (km/s):')
    print(f' v = ({V[0]}, {V[1]}, {V[2]})')
    print('-----------------------------------------------------')

if __name__ == "__main__":
    rv_from_r0v0_test()