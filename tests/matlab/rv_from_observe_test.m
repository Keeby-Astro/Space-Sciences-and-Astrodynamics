% rv_from_observe_test
%
% This program uses Algorithms 5.4 and 4.2 to obtain the orbital
% elements from the observational data.
%
% deg    - conversion factor between degrees and radians
% pi     - 3.1415926...
% mu     - gravitational parameter (km^3/s^2)

% Re     - equatorial radius of the earth (km)
% f      - earth's flattening factor
% wE     - angular velocity of the earth (rad/s)
% omega  - earth's angular velocity vector (rad/s) in the
%          geocentric equatorial frame

% rho    - slant range of object (km)
% rhodot - range rate (km/s)
% A      - azimuth (deg) of object relative to observation site
% Adot   - time rate of change of azimuth (deg/s)
% a      - elevation angle (deg) of object relative to observation site
% adot   - time rate of change of elevation angle (degrees/s)

% theta  - local sidereal time (deg) of tracking site
% phi    - geodetic latitude (deg) of site
% H      - elevation of site (km)

% r      - geocentric equatorial position vector of object (km)
% v      - geocentric equatorial velocity vector of object (km)

% coe    - orbital elements [h e RA incl w TA a]
%          where
%              h    = angular momentum (km^2/s)
%              e    = eccentricity
%              RA   = right ascension of the ascending node (rad)
%              incl = inclination of the orbit (rad)
%              w    = argument of perigee (rad)
%              TA   = true anomaly (rad)
%              a    = semimajor axis (km)
% rp     - perigee radius (km)
% T      - period of elliptical orbit (s)
%
% User M-functions required: rv_from_observe, coe_from_sv
% --------------------------------------------------------------------
clear all; clc
global f Re wE

deg = pi/180;
f   = 1/298.256421867;
Re  = 6378.13655;
wE  = 7.292115e-5;
mu  = 398600.4418;

%...Data declaration:
rho    = 2551;
rhodot = 0;
A      = 90;
Adot   = 0.1130;
a      = 30;
adot   = 0.05651;
theta  = 300;
phi    = 60;
H      = 0;
%...

%...Algorithm 5.4:
[r,v] = rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H);

%...Algorithm 4.2:
coe = coe_from_sv(r,v,mu);

h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);
a    = coe(7);

%...Equation 2.40:
rp = h^2/mu/(1 + e);

%...Echo the input data and output the solution to
% the command window:
fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n');
fprintf('\n Slant range (km)              = %g', rho);
fprintf('\n Slant range rate (km/s)       = %g', rhodot);
fprintf('\n Azimuth (deg)                 = %g', A);
fprintf('\n Azimuth rate (deg/s)          = %g', Adot);
fprintf('\n Elevation (deg)               = %g', a);
fprintf('\n Elevation rate (deg/s)        = %g', adot);
fprintf('\n Local sidereal time (deg)     = %g', theta);
fprintf('\n Latitude (deg)                = %g', phi);
fprintf('\n Altitude above sea level (km) = %g', H);
fprintf('\n\n');

fprintf(' Solution:')

fprintf('\n\n State vector:\n');
fprintf('\n r (km)                        = [%g, %g, %g]', ...
                                        r(1), r(2), r(3));
fprintf('\n v (km/s)                      = [%g, %g, %g]', ...
                                        v(1), v(2), v(3));

fprintf('\n\n Orbital elements:\n')
fprintf('\n   Angular momentum (km^2/s)   = %g', h)
fprintf('\n   Eccentricity                = %g', e)
fprintf('\n   Inclination (deg)           = %g', incl/deg)
fprintf('\n   RA of ascending node (deg)  = %g', RA/deg)
fprintf('\n   Argument of perigee (deg)   = %g', w/deg)
fprintf('\n   True anomaly (deg)          = %g\n', TA/deg)
fprintf('\n   Semimajor axis (km)         = %g', a)
fprintf('\n   Perigee radius (km)         = %g', rp)
%...If the orbit is an ellipse, output its period:
if e < 1
    T = 2*pi/sqrt(mu)*a^1.5;
    fprintf('\n   Period:')
    fprintf('\n     Seconds = %g', T)
    fprintf('\n     Minutes = %g', T/60)
    fprintf('\n     Hours   = %g', T/3600)
    fprintf('\n     Days    = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------------\n')
