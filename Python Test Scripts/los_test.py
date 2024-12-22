import los

def los_test():
    '''
    Uses the ECI position vectors of the satellite (r_sat) and the sun (r_sun)
    and evaluates whether the earth is in the line of sight between the two.

    User py-functions required: None
    '''
    # Scenario 1: Satellite in sunlight
    r_sat1 = [7000, 0, 0]    # Satellite position vector (ECI, km)
    r_sun1 = [1.496e8, 0, 0] # Sun position vector (ECI, km)
    light_switch1 = los.los(r_sat1, r_sun1)

    if light_switch1 == 1:
        print('Satellite 1 is in sunlight')
    else:
        print('Satellite 1 is not in sunlight')

    # Scenario 2: Satellite in Earth\'s shadow
    r_sat2 = [7000, 0, 0]     # Satellite position vector (ECI, km)
    r_sun2 = [-1.496e8, 0, 0] # Sun position vector (ECI, km)
    light_switch2 = los.los(r_sat2, r_sun2)

    if light_switch2 == 0:
        print('Satellite 2 is in Earth\'s shadow')
    else:
        print('Satellite 2 is not in Earth\'s shadow')

if __name__ == "__main__":
    los_test()