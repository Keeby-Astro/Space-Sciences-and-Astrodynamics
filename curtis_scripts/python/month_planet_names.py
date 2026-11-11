def month_planet_names(month_id, planet_id):
    '''
    This function returns the name of the month and the planet
    corresponding, respectively, to the numbers "month_id" and
    "planet_id".

    months    - a vector containing the names of the 12 months
    planets   - a vector containing the names of the 9 planets
    month_id  - the month number (1 - 12)
    planet_id - the planet number (1 - 9)

    User py-functions required: none
    '''
    months = ['January  ',
              'February ',
              'March    ',
              'April    ',
              'May      ',
              'June     ',
              'July     ',
              'August   ',
              'September',
              'October  ',
              'November ',
              'December ']

    planets = ['Mercury',
               'Venus  ',
               'Earth  ',
               'Mars   ',
               'Jupiter',
               'Saturn ',
               'Uranus ',
               'Neptune',
               'Pluto  ']

    if not isinstance(month_id, int) or not 1 <= month_id <= 12 or not isinstance(planet_id, int) or not 1 <= planet_id <= 9:
        raise ValueError("Month and planet identifiers must be integers in 1..12 and 1..9.")
    month = months[month_id - 1]
    planet = planets[planet_id - 1]

    return month, planet