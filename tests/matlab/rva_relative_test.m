%rva_relative_test
%{
  This program uses the data to calculate the position,
  velocity and acceleration of an orbiting chaser B relative to an
  orbiting target A.
 
  mu               - gravitational parameter (km^3/s^2)
  deg              - conversion factor from degrees to radians
 
                   Spacecraft A & B:
  h_A, h_B         -   angular momentum (km^2/s)
  e_A, E_B         -   eccentricity
  i_A, i_B         -   inclination (radians)
  RAAN_A, RAAN_B   -   right ascension of the ascending node (radians)
  omega_A, omega_B -   argument of perigee (radians)
  theta_A, theta_A -   true anomaly (radians)
 
  rA, vA           - inertial position (km) and velocity (km/s) of A
  rB, vB           - inertial position (km) and velocity (km/s) of B
  r                - position (km) of B relative to A in A's
                     co-moving frame
  v                - velocity (km/s) of B relative to A in A's
                     co-moving frame
  a                - acceleration (km/s^2) of B relative to A in A's
                     co-moving frame
 
  User M-function required:   sv_from_coe, rva_relative
%}
% -----------------------------------------------------------
 
clear all; clc
global mu
mu  = 398600;
deg = pi/180;
 
%...Input data:
 
%   Spacecraft A:
h_A     = 52059;
e_A     = 0.025724;
i_A     = 60*deg;
RAAN_A  = 40*deg;
omega_A = 30*deg;
theta_A = 40*deg;
 
%   Spacecraft B:
h_B     = 52362;
e_B     = 0.0072696;
i_B     = 50*deg;
RAAN_B  = 40*deg;
omega_B = 120*deg;
theta_B = 40*deg;
 
%...End input data
 
%...Compute the initial state vectors of A and B using Algorithm 4.5:
[rA,vA] = sv_from_coe([h_A e_A RAAN_A i_A omega_A theta_A],mu);
[rB,vB] = sv_from_coe([h_B e_B RAAN_B i_B omega_B theta_B],mu);
 
%...Compute relative position, velocity and acceleration using
%   Algorithm 7.1:
[r,v,a] = rva_relative(rA,vA,rB,vB);
 
%...Output
fprintf('\n\n--------------------------------------------------------\n\n')
fprintf('\nOrbital parameters of spacecraft A:')
fprintf('\n   angular momentum    = %g (km^2/s)', h_A)
fprintf('\n   eccentricity        = %g'         , e_A)
fprintf('\n   inclination         = %g (deg)'   , i_A/deg)
fprintf('\n   RAAN                = %g (deg)'   , RAAN_A/deg)
fprintf('\n   argument of perigee = %g (deg)'   , omega_A/deg)
fprintf('\n   true anomaly        = %g (deg)\n' , theta_A/deg)
 
fprintf('\nState vector of spacecraft A:')
fprintf('\n   r = [%g, %g, %g]', rA(1), rA(2), rA(3))
fprintf('\n       (magnitude = %g)', norm(rA))
fprintf('\n   v = [%g, %g, %g]', vA(1), vA(2), vA(3))
fprintf('\n       (magnitude = %g)\n', norm(vA))
 
fprintf('\nOrbital parameters of spacecraft B:')
fprintf('\n   angular momentum    = %g (km^2/s)', h_B)
fprintf('\n   eccentricity        = %g'         , e_B)
fprintf('\n   inclination         = %g (deg)'   , i_B/deg)
fprintf('\n   RAAN                = %g (deg)'   , RAAN_B/deg)
fprintf('\n   argument of perigee = %g (deg)'   , omega_B/deg)
fprintf('\n   true anomaly        = %g (deg)\n' , theta_B/deg)
 
fprintf('\nState vector of spacecraft B:')
fprintf('\n   r = [%g, %g, %g]', rB(1), rB(2), rB(3))
fprintf('\n       (magnitude = %g)', norm(rB))
fprintf('\n   v = [%g, %g, %g]', vB(1), vB(2), vB(3))
fprintf('\n       (magnitude = %g)\n', norm(vB))
 
fprintf('\nIn the co-moving frame attached to A:')
fprintf('\n   Position of B relative to A = [%g, %g, %g]', ...
                           r(1), r(2), r(3))
fprintf('\n      (magnitude = %g)\n', norm(r))
fprintf('\n   Velocity of B relative to A = [%g, %g, %g]', ...
                           v(1), v(2), v(3))
fprintf('\n      (magnitude = %g)\n', norm(v))
fprintf('\n   Acceleration of B relative to A = [%g, %g, %g]', ...
                           a(1), a(2), a(3))
fprintf('\n      (magnitude = %g)\n', norm(a))
fprintf('\n\n--------------------------------------------------------\n\n')