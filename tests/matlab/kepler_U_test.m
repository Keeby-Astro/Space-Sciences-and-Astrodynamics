% kepler_U_test
%{
    This program uses Algorithm 3.3 and the data
    to solve the universal Kepler's equation.

    mu  - gravitational parameter (km^3/s^2)
    x   - the universal anomaly (km^0.5)
    dt  - time since x = 0 (s)
    ro  - radial position when x = 0 (km)
    vro - radial velocity when x = 0 (km/s)
    a   - semimajor axis (km)

    User M-function required: kepler_U
%}
% ---------------------------------------------

clear all; clc
global mu
mu = 398600;

%...Data declaration:
ro = 10000;
vro = 3.0752;
dt = 3600;
a = -19655;
%...

%...Pass the input data to the function keplfr_U, which returns x
%...(Universal Kepler's requires the reciprocal of semimajor axis):
x = kepler_U(dt, ro, vro, 1/a);

...Echo the input data and output the results to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Initial radial coordinate (km) = %g',ro)
fprintf('\n Initial radial velocity (km/s) = %g',vro)
fprintf('\n Elapsed time (seconds)         = %g',dt)
fprintf('\n Semimajor axis (km)            = %g\n',a)
fprintf('\n Universal anomaly (km^0.5)     = %g',x)
fprintf('\n-----------------------------------------------------\n')
