% sv_from_coe_test
%{
    This program uses Algorithm 4.5 to obtain the state vector from
    the orbital elements.

    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    coe - orbital elements [h e RA incl w TA a]
          where h = angular momentum (km^2/s)
          e       = eccentricity
          RA      = right ascension of the ascending node (rad)
          incl    = orbit inclination (rad)
          w       = argument of perigee (rad)
          TA      = true anomaly (rad)
          a       = semimajor axis (km)
    r   - position vector (km) in geocentric equatorial frame
    v   - velocity vector (km) in geocentric equatorial frame

    User M-function required: sv_from_coe
%}
% ----------------------------------------------
clear all; clc
deg = pi/180;
mu  = 398600;

%...Data declaration (angles in degrees):
h    = 80000;
e    = 1.4;
RA   = 40;
incl = 30;
w    = 60;
TA   = 30;
%...

coe = [h, e, RA*deg, incl*deg, w*deg, TA*deg];

%...Algorithm 4.5 (requires angular elements be in radians):
[r, v] = sv_from_coe(coe, mu);

%...Echo the input data and output the results to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Gravitational parameter (km^3/s^2)  = %g\n', mu)
fprintf('\n Angular momentum (km^2/s)           = %g', h)
fprintf('\n Eccentricity                        = %g', e)
fprintf('\n Right ascension (deg)               = %g', RA)
fprintf('\n Argument of perigee (deg)           = %g', w)
fprintf('\n True anomaly (deg)                  = %g', TA)
fprintf('\n\n State vector:')
fprintf('\n   r (km)   = [%g %g %g]', r(1), r(2), r(3))
fprintf('\n   v (km/s) = [%g %g %g]', v(1), v(2), v(3))
fprintf('\n-----------------------------------------------------\n')
