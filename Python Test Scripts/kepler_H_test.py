import kepler_H

def kepler_H_test():
    '''
    This program uses Algorithm 3.2 and the data to solve
    Kepler's equation for the hyperbola.

    e - eccentricity
    M - hyperbolic mean anomaly (dimensionless)
    E - hyperbolic eccentric anomaly (dimensionless)

    User py-function required: kepler_H
    '''
    # ...Data declaration:
    e = 2.7696
    M = 40.69
    # ...

    # ...Pass the input data to the function kepler_E, which returns E:
    E = kepler_H.kepler_H(e, M)

    # ...Echo the input data and output to the command window:
    print('---------------------------------------------------')
    print(f'\n Eccentricity                 = {e}')
    print(f'\n Hyperbolic mean anomaly      = {M}')
    print(f'\n Hyperbolic eccentric anomaly = {E}')
    print('---------------------------------------------------')

if __name__ == '__main__':
    kepler_H_test()