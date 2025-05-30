% gauss_test
%{
    This program uses Algorithms 5.5 and 5.6 (Gauss's method) to compute
    the state vector from the data provided.

    deg          - factor for converting between degrees and radians
    pi           - 3.1415926...
    mu           - gravitational parameter (km^3/s^2)
    Re           - earth's radius (km)
    f            - earth's flattening factor
    H            - elevation of observation site (km)
    phi          - latitude of site (deg)
    t            - vector of observation times t1, t2, t3 (s)
    ra           - vector of topocentric equatorial right ascensions
                    at t1, t2, t3 (deg)
    dec          - vector of topocentric equatorial right declinations
                    at t1, t2, t3 (deg)
    theta        - vector of local sidereal times for t1, t2, t3 (deg)
    R            - matrix of site position vectors at t1, t2, t3 (km)
    rho          - matrix of direction cosine vectors at t1, t2, t3
    fac1, fac2   - common factors
    r_old, v_old - the state vector without iterative improvement (km, km/s)
    r, v         - the state vector with iterative improvement (km, km/s)
    coe          - vector of orbital elements for r, v: 
                    [h, e, RA, incl, w, TA, a]
                    where h    = angular momentum (km^2/s)
                        e    = eccentricity
                        incl = inclination (rad)
                        w    = argument of perigee (rad)
                        TA   = true anomaly (rad)
                        a    = semimajor axis (km)
    coe_old     - vector of orbital elements for r_old, v_old       

    User M-functions required: gauss, coe_from_sv
%}
% ---------------------------------------------
 
clear all; clc
 
global mu
 
deg = pi/180;
mu  = 398600;
Re  = 6378;
f   = 1/298.26;
 
%...Data declaration for Example 5.11:
H     = 1;
phi   = 40*deg;
t     = [       0   118.104   237.577];
ra    = [ 43.5365   54.4196   64.3178]*deg;
dec   = [-8.78334  -12.0739  -15.1054]*deg;
theta = [ 44.5065    45.000   45.4992]*deg;
%...
 
%...Equations 5.64, 5.76 and 5.79:
fac1 = Re/sqrt(1-(2*f - f*f)*sin(phi)^2);
fac2 = (Re*(1-f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*sin(phi);
for i = 1:3
    R(i,1) = (fac1 + H)*cos(phi)*cos(theta(i));
    R(i,2) = (fac1 + H)*cos(phi)*sin(theta(i));
    R(i,3) = fac2;
    rho(i,1) = cos(dec(i))*cos(ra(i));
    rho(i,2) = cos(dec(i))*sin(ra(i));
    rho(i,3) = sin(dec(i));
end
 
%...Algorithms 5.5 and 5.6:
[r, v, r_old, v_old] = gauss(rho(1,:), rho(2,:), rho(3,:), ...
                               R(1,:),    R(2,:),  R(3,:), ...
                               t(1),      t(2),    t(3));
 
%...Algorithm 4.2 for the initial estimate of the state vector
%   and for the iteratively improved one:
coe_old = coe_from_sv(r_old,v_old,mu);
coe     = coe_from_sv(r,v,mu);
 
%...Echo the input data and output the solution to
%   the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.11: Orbit determination by the Gauss method\n')
fprintf('\n Radius of earth (km)               = %g', Re)
fprintf('\n Flattening factor                  = %g', f)
fprintf('\n Gravitational parameter (km^3/s^2) = %g', mu)
fprintf('\n\n Input data:\n');
fprintf('\n Latitude (deg)                = %g', phi/deg);
fprintf('\n Altitude above sea level (km) = %g', H);
fprintf('\n\n Observations:')
fprintf('\n               Right')
fprintf('                                     Local')
fprintf('\n   Time (s)   Ascension (deg)   Declination (deg)') 
fprintf('   Sidereal time (deg)')
for i = 1:3
    fprintf('\n %9.4g %11.4f %19.4f %20.4f', ...
                 t(i), ra(i)/deg, dec(i)/deg, theta(i)/deg)
end
 
fprintf('\n\n Solution:\n')
 
fprintf('\n Without iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km)                          = [%g, %g, %g]', ...
                                   r_old(1), r_old(2), r_old(3))
fprintf('\n v (km/s)                        = [%g, %g, %g]', ...
                                   v_old(1), v_old(2), v_old(3))
fprintf('\n');
 
fprintf('\n   Angular momentum (km^2/s)     = %g', coe_old(1))
fprintf('\n   Eccentricity                  = %g', coe_old(2))
fprintf('\n   RA of ascending node (deg)    = %g', coe_old(3)/deg)
fprintf('\n   Inclination (deg)             = %g', coe_old(4)/deg)
fprintf('\n   Argument of perigee (deg)     = %g', coe_old(5)/deg)
fprintf('\n   True anomaly (deg)            = %g', coe_old(6)/deg)
fprintf('\n   Semimajor axis (km)           = %g', coe_old(7))
fprintf('\n   Periapse radius (km)          = %g', coe_old(1)^2 ...
                                             /mu/(1 + coe_old(2)))
%...If the orbit is an ellipse, output the period:
if coe_old(2)<1
    T = 2*pi/sqrt(mu)*coe_old(7)^1.5; 
    fprintf('\n   Period:')
    fprintf('\n     Seconds                     = %g', T) 
    fprintf('\n     Minutes                     = %g', T/60)
    fprintf('\n     Hours                       = %g', T/3600)
    fprintf('\n     Days                        = %g', T/24/3600)
end
 
fprintf('\n\n With iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km)                          = [%g, %g, %g]', ...
                                               r(1), r(2), r(3))
fprintf('\n v (km/s)                        = [%g, %g, %g]', ...
                                               v(1), v(2), v(3))
fprintf('\n');
fprintf('\n   Angular momentum (km^2/s)     = %g', coe(1))
fprintf('\n   Eccentricity                  = %g', coe(2))
fprintf('\n   RA of ascending node (deg)    = %g', coe(3)/deg)
fprintf('\n   Inclination (deg)             = %g', coe(4)/deg)
fprintf('\n   Argument of perigee (deg)     = %g', coe(5)/deg)
fprintf('\n   True anomaly (deg)            = %g', coe(6)/deg)
fprintf('\n   Semimajor axis (km)           = %g', coe(7))
fprintf('\n   Periapse radius (km)          = %g', coe(1)^2 ...
                                             /mu/(1 + coe(2)))
%...If the orbit is an ellipse, output the period:
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5; 
    fprintf('\n   Period:')
    fprintf('\n     Seconds                     = %g', T) 
    fprintf('\n     Minutes                     = %g', T/60)
    fprintf('\n     Hours                       = %g', T/3600)
    fprintf('\n     Days                        = %g', T/24/3600)
end
fprintf('\n-----------------------------------------------------\n')
