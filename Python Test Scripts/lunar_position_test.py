import lunar_position

def lunar_position_test():
    '''
    Outputs the geocentric equatorial position vector of the moon
    given the Julian day.

    User M-functions required: None
    '''
    jd_1 = 2451545.0  # Julian Date for J2000 epoch
    result_1 = lunar_position.lunar_position(jd_1)
    jd_2 = 2459200.5  # Julian Date for September 22, 2020
    result_2 = lunar_position.lunar_position(jd_2)
    jd_3 = 2440587.5  # Julian Date for January 1, 1970
    result_3 = lunar_position.lunar_position(jd_3)

    # Display results
    print(f'Julian Date = {jd_1:.1f}')
    print('Geocentric position vector (km):')
    print(f'  X = {result_1[0]:.3f} km')
    print(f'  Y = {result_1[1]:.3f} km')
    print(f'  Z = {result_1[2]:.3f} km\n')

    print(f'Julian Date = {jd_2:.1f}')
    print('Geocentric position vector (km):')
    print(f'  X = {result_2[0]:.3f} km')
    print(f'  Y = {result_2[1]:.3f} km')
    print(f'  Z = {result_2[2]:.3f} km\n')

    print(f'Julian Date = {jd_3:.1f}')
    print('Geocentric position vector (km):')
    print(f'  X = {result_3[0]:.3f} km')
    print(f'  Y = {result_3[1]:.3f} km')
    print(f'  Z = {result_3[2]:.3f} km\n')

if __name__ == "__main__":
    lunar_position_test()