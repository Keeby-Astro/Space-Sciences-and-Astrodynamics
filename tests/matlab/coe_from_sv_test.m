% coe_from_sv_test
%{
    This program uses Algorithm 4.2 to obtain the orbital
    elements from the state vector.
    pi  - 3.1415926...
    deg - factor for converting between degrees and radians
    mu  - gravitational parameter (km^3/s^2)
    r   - position vector (km) in the geocentric equatorial frame
    v   - velocity vector (km/s) in the geocentric equatorial frame
    coe - orbital elements [h e RA incl w TA a]
          where h    = angular momentum (km^2/s)
                e    = eccentricity
                RA   = right ascension of the ascending node (rad)
                incl = orbit inclination (rad)
                w    = argument of perigee (rad)
                TA   = true anomaly (rad)
                a    = semimajor axis (km)
    T   - Period of an elliptic orbit (s)

    User M-function required: coe_from_sv
%}
% ----------------------------------------------
clear all; clc
deg = pi/180;
mu = 398600;

%...Data declaration:
r = [ -6045 -3490 2500];
v = [-3.457 6.618 2.533];
%...

%...Algorithm 4.2:
coe = coe_from_sv(r,v,mu);

%...Echo the input data and output results to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu)
fprintf('\n State vector:\n')
fprintf('\n r (km)                             = [%g %g %g]', ...
                                                 r(1), r(2), r(3))
fprintf('\n v (km/s)                           = [%g %g %g]', ...
                                                 v(1), v(2), v(3))
disp(' ')
fprintf('\n Angular momentum (km^2/s)          = %g', coe(1))
fprintf('\n Eccentricity                       = %g', coe(2))
fprintf('\n Right ascension (deg)              = %g', coe(3)/deg)
fprintf('\n Inclination (deg)                  = %g', coe(4)/deg)
fprintf('\n Argument of perigee (deg)          = %g', coe(5)/deg)
fprintf('\n True anomaly (deg)                 = %g', coe(6)/deg)
fprintf('\n Semimajor axis (km):               = %g', coe(7))

%...if the orbit is an ellipse, output its period (Equation 2.73):
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n Period:')
    fprintf('\n Seconds = %g', T)
    fprintf('\n Minutes = %g', T/60)
    fprintf('\n Hours   = %g', T/3600)
    fprintf('\n Days    = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------------\n')