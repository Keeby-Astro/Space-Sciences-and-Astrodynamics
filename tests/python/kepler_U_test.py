import kepler_U

def kepler_U_test():
    '''
    This program uses Algorithm 3.3 and the data
    to solve the universal Kepler's equation.
    
    mu  - gravitational parameter (km^3/s^2)
    x   - the universal anomaly (km^0.5)
    dt  - time since x = 0 (s)
    ro  - radial position when x = 0 (km)
    vro - radial velocity when x = 0 (km/s)
    a   - semimajor axis (km)

    User py-function required: kepler_U
    '''
    mu = 398600.4418

    # Data declaration:
    ro = 10000
    vro = 3.0752
    dt = 3600
    a = -19655

    # Pass the input data to the function kepler_U, which returns x
    # (Universal Kepler's requires the reciprocal of semimajor axis):
    x = kepler_U.kepler_U(dt, ro, vro, 1/a, mu)

    # Echo the input data and output the results to the console:
    print('-----------------------------------------------------')
    print(f'Initial radial coordinate (km) = {ro}')
    print(f'Initial radial velocity (km/s) = {vro}')
    print(f'Elapsed time (seconds)         = {dt}')
    print(f'Semimajor axis (km)            = {a}')
    print(f'Universal anomaly (km^0.5)     = {x}')
    print('-----------------------------------------------------')

if __name__ == '__main__':
    kepler_U_test()