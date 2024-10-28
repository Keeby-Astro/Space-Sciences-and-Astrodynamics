# Algorithm 10.3: Determine whether or not a satellite is in earth’s shadow.
import numpy as np

def los(r_sat, r_sun):
    """
    This function uses the ECI position vectors of the satellite (r_sat)
    and the sun (r_sun) to determine whether the earth is in the line of
    sight between the two.
    
    User py-functions required: None
    """
    RE = 6378  # Earth's radius (km)
    rsat = np.linalg.norm(r_sat)
    rsun = np.linalg.norm(r_sun)

    # Angle between sun and satellite position vectors:
    theta = np.degrees(np.arccos(np.dot(r_sat, r_sun) / (rsat * rsun)))

    # Angle between the satellite position vector and the radial to the point
    # of tangency with the earth of a line from the satellite:
    theta_sat = np.degrees(np.arccos(RE / rsat))

    # Angle between the sun position vector and the radial to the point
    # of tangency with the earth of a line from the sun:
    theta_sun = np.degrees(np.arccos(RE / rsun))

    # Determine whether a line from the sun to the satellite
    # intersects the earth:
    if theta_sat + theta_sun <= theta:
        return 0  # yes
    else:
        return 1  # no
