# ALGORITHM 4.5: CALCULATION OF THE STATE VECTOR FROM THE ORBITAL ELEMENTS

import numpy as np

def sv_from_coe(coe, mu):
    '''
    This function computes the state vector (r,v) from the
    classical orbital elements (coe).

    mu   - gravitational parameter (km^3/s^2)
    coe  - orbital elements [h e RA incl w TA]
           where
               h = angular momentum (km^2/s)
               e = eccentricity
               RA = right ascension of the ascending node (rad)
               incl = inclination of the orbit (rad)
               w = argument of perigee (rad)
               TA = true anomaly (rad)
    R3_w - Rotation matrix about the z-axis through the angle w
    R1_i - Rotation matrix about the x-axis through the angle i
    R3_W - Rotation matrix about the z-axis through the angle RA
    Q_pX - Matrix of the transformation from perifocal to geocentric
           equatorial frame
    rp   - position vector in the perifocal frame (km)
    vp   - velocity vector in the perifocal frame (km/s)
    r    - position vector in the geocentric equatorial frame (km)
    v    - velocity vector in the geocentric equatorial frame (km/s)
    
    User py-functions required: none
    ''' 
    h    = coe[0]
    e    = coe[1]
    RA   = coe[2]
    incl = coe[3]
    w    = coe[4]
    TA   = coe[5]

    #...Equations 4.45 and 4.46 (rp and vp are column vectors):
    rp = (h**2/mu) * (1/(1 + e*np.cos(TA))) * (np.cos(TA)*np.array([1, 0, 0]) + np.sin(TA)*np.array([0, 1, 0]))
    vp = (mu/h) * (-np.sin(TA)*np.array([1, 0, 0]) + (e + np.cos(TA))*np.array([0, 1, 0]))

    #...Equation 4.34:
    R3_W = np.array([[ np.cos(RA),  np.sin(RA), 0],
                     [-np.sin(RA),  np.cos(RA), 0],
                     [           0,           0, 1]])

    #...Equation 4.32:
    R1_i = np.array([[1,            0,            0],
                     [0,  np.cos(incl),  np.sin(incl)],
                     [0, -np.sin(incl),  np.cos(incl)]])

    #...Equation 4.34:
    R3_w = np.array([[ np.cos(w),  np.sin(w), 0],
                     [-np.sin(w),  np.cos(w), 0],
                     [          0,           0, 1]])

    #...Equation 4.49:
    Q_pX = (R3_w @ R1_i @ R3_W).T

    #...Equations 4.51 (r and v are column vectors):
    r = Q_pX @ rp
    v = Q_pX @ vp

    #...Convert r and v into row vectors:
    r = r.T
    v = v.T

    return r, v