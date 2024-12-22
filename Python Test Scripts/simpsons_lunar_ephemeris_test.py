import simpsons_lunar_ephemeris

def simpsons_lunar_ephemeris_test():
    '''
    Outputs the state vector of the moon at a given time
    relative to the earth's geocentric equatorial frame using a curve fit
    to JPL's DE200 (1982) ephemeris model.

    jd      - julian date (days)
    pos     - position vector (km)
    vel     - velocity vector (km/s)
    a       - matrix of amplitudes (km)
    b       - matrix of frequencies (rad/century)
    c       - matrix of phase angles (rad)
    t       - time in centuries since J2000
    tfac    - no. of seconds in a Julian century (36525 days)

    User py-functions required: simpsons_lunar_ephemeris
    '''
    # Define a sample Julian date (e.g., January 1, 2024, at 12:00 TT)
    jd_sample = 2460271.0  # Julian date

    # Call the simpsons_lunar_ephemeris function
    pos, vel = simpsons_lunar_ephemeris.simpsons_lunar_ephemeris(jd_sample)

    # Display the Juilan date
    print('\nJulian Date:')
    print(f'{jd_sample:.1f}\n')

    # Display the position vector
    print('Position vector (km):')
    print(f'X: {pos[0]:.6f} km')
    print(f'Y: {pos[1]:.6f} km')
    print(f'Z: {pos[2]:.6f} km')

    # Display the velocity vector
    print('\nVelocity vector (km/s):')
    print(f'VX: {vel[0]:.6f} km/s')
    print(f'VY: {vel[1]:.6f} km/s')
    print(f'VZ: {vel[2]:.6f} km/s')

if __name__ == "__main__":
    simpsons_lunar_ephemeris_test()