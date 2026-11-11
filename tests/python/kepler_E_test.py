import sys
from pathlib import Path

# Locate companion scripts relative to this test file
_SCRIPTS_DIR = str(Path(__file__).resolve().parents[2] / "curtis_scripts" / "python")
if _SCRIPTS_DIR not in sys.path:
    sys.path.append(_SCRIPTS_DIR)

import kepler_E

def kepler_E_test():
    '''
    This program uses Algorithm 3.1 and the data to solve
    Kepler's equation.

    e - eccentricity
    M - mean anomaly (rad)
    E - eccentric anomaly (rad)

    User function required: kepler_E
    '''
    #...Data declaration:
    e = 0.37255
    M = 3.6029
    #...

    #...Pass the input data to the function kepler_E, which returns E:
    E = kepler_E.kepler_E(e, M)

    #...Echo the input data and output to the console:
    print('---------------------------------------------------')
    print(f'\n Eccentricity                = {e}')
    print(f'\n Mean anomaly (radians)      = {M}')
    print(f'\n Eccentric anomaly (radians) = {E}')
    print('---------------------------------------------------')

if __name__ == '__main__':
    kepler_E_test()