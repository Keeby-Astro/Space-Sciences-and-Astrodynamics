import numpy as np
import planet_elements_and_sv
import month_planet_names

def compute_planet_elements_and_sv():
    '''
    This program uses Algorithm 8.1 to compute the orbital elements
    and state vector of the earth at the date and time specified.
    To obtain the same results for Mars, set planet_id = 4.

    mu        - gravitational parameter of the sun (km^3/s^2)
    deg       - conversion factor between degrees and radians
    pi        - 3.1415926...

    coe       - vector of heliocentric orbital elements
                [h  e  RA  incl  w  TA  a  w_hat  L  M  E],
                where
                h     = angular momentum                    (km^2/s)
                e     = eccentricity
                RA    = right ascension                     (deg)
                incl  = inclination                         (deg)
                w     = argument of perihelion              (deg)
                TA    = true anomaly                        (deg)
                a     = semimajor axis                      (km)
                w_hat = longitude of perihelion ( = RA + w) (deg)
                L     = mean longitude ( = w_hat + M)       (deg)
                M     = mean anomaly                        (deg)
                E     = eccentric anomaly                   (deg)

    r         - heliocentric position vector (km)
    v         - heliocentric velocity vector (km/s) 

    planet_id - planet identifier:
                1 = Mercury
                2 = Venus
                3 = Earth
                4 = Mars
                5 = Jupiter
                6 = Saturn
                7 = Uranus
                8 = Neptune
                9 = Pluto

    year      - range: 1901 - 2099
    month     - range: 1 - 12
    day       - range: 1 - 31
    hour      - range: 0 - 23
    minute    - range: 0 - 60
    second    - range: 0 - 60 

    User py-functions required: planet_elements_and_sv, month_planet_names
    '''
    global mu
    mu = 1.327124e11
    deg = np.pi / 180

    # Input data
    planet_id = 3
    year      = 2003
    month     = 8
    day       = 27
    hour      = 12
    minute    = 0
    second    = 0

    # Algorithm 8.1:
    coe, r, v, jd = planet_elements_and_sv.planet_elements_and_sv(planet_id, year, month, day, 
                                                                  hour, minute, second)

    # Convert the planet_id and month numbers into names for output:
    month_name, planet_name = month_planet_names.month_planet_names(month, planet_id)

    # Echo the input data and output the solution
    print('-----------------------------------------------------')
    print('\n Input data:')
    print(f'\n   Planet: {planet_name}')
    print(f'   Year  : {year}')
    print(f'   Month : {month_name}')
    print(f'   Day   : {day}')
    print(f'   Hour  : {hour}')
    print(f'   Minute: {minute}')
    print(f'   Second: {second}')
    print(f'\n\n   Julian day: {jd:.3f}')

    print('\n\n Orbital elements:')

    print(f'\n  Angular momentum (km^2/s)                   = {coe[0]}')
    print(f'  Eccentricity                                = {coe[1]}')
    print(f'  Right ascension of the ascending node (deg) = {coe[2]}')
    print(f'  Inclination to the ecliptic (deg)           = {coe[3]}')
    print(f'  Argument of perihelion (deg)                = {coe[4]}')
    print(f'  True anomaly (deg)                          = {coe[5]}')
    print(f'  Semimajor axis (km)                         = {coe[6]}')

    print(f'\n  Longitude of perihelion (deg)               = {coe[7]}')
    print(f'  Mean longitude (deg)                        = {coe[8]}')
    print(f'  Mean anomaly (deg)                          = {coe[9]}')
    print(f'  Eccentric anomaly (deg)                     = {coe[10]}')

    print('\n\n State vector:')

    print(f'\n  Position vector (km) = [{r[0]}  {r[1]}  {r[2]}]')
    print(f'  Magnitude            = {np.linalg.norm(r):.5f}')
    print(f'\n  Velocity (km/s)      = [{v[0]}  {v[1]}  {v[2]}]')
    print(f'  Magnitude            = {np.linalg.norm(v):.5f}')

    print('-----------------------------------------------------')

if __name__ == '__main__':
    compute_planet_elements_and_sv()