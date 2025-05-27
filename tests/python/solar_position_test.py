import solar_position

def solar_position_test():
    '''
    Outputs the geocentric equatorial position vector
    of the sun, given the Julian date.

    User py-functions required: None
    '''
    # Julian date for 2023-12-21 00:00:00 UTC (Winter Solstice)
    jd_test = 2460271.5  # Julian date

    # Call the solar_position function
    lamda, eps, r_S = solar_position.solar_position(jd_test)

    # Display the results
    print(f'Test Case: Julian Date = {jd_test:.2f}')
    print(f'Apparent Ecliptic Longitude (lamda): {lamda:.6f} degrees')
    print(f'Obliquity of the Ecliptic (eps): {eps:.6f} degrees')
    print('Geocentric Position Vector (r_S):')
    print(f'  x = {r_S[0]:.6f} km')
    print(f'  y = {r_S[1]:.6f} km')
    print(f'  z = {r_S[2]:.6f} km')

if __name__ == '__main__':
    solar_position_test()