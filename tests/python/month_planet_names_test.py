import sys
from pathlib import Path

# Locate companion scripts relative to this test file
_SCRIPTS_DIR = str(Path(__file__).resolve().parents[2] / "curtis_scripts" / "python")
if _SCRIPTS_DIR not in sys.path:
    sys.path.append(_SCRIPTS_DIR)

import month_planet_names

def month_planet_names_test():
    '''
    Outputs the name of the month and the planet corresponding, respectively,
    to the numbers "month_id" and "planet_id".

    month_id  - the month number (1 - 12)
    planet_id - the planet number (1 - 9)

    User functions required: month_planet_names
    '''
    print('Month and Planet Names Test')
    print('---------------------------')

    # Test month and planet IDs
    for month_id in range(1, 13):  # Months 1 to 12
        for planet_id in range(1, 10):  # Planets 1 to 9
            month, planet = month_planet_names.month_planet_names(month_id, planet_id)
            
            # Output
            print(f'Month ID: {month_id:2d}, Month: {month.strip()} | Planet ID: {planet_id:2d}, Planet: {planet.strip()}')

if __name__ == "__main__":
    month_planet_names_test()