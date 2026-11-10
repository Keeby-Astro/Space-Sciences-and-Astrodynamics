% kepler_H_test
%{
    This program uses Algorithm 3.2 and the data to solve 
    Kepler's equation for the hyperbola.

    e - eccentricity
    M - hyperbolic mean anomaly (dimensionless)
    F - hyperbolic eccentric anomaly (dimensionless)

    User M-function required: kepler_H
%}
% -------------------------------------------

% Locate companion functions relative to this test file
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'curtis_scripts', 'matlab'));

clear

%...Data declaration:
e = 2.7696;
M = 40.69;
%...

%...Pass the input data to the function kepler_H, which returns F:
F = kepler_H(e, M);

%...Echo the input data and output to the command window:
fprintf('----------------------------------------------')
fprintf('\n Eccentricity                 = %g',e)
fprintf('\n Hyperbolic mean anomaly      = %g\n',M)
fprintf('\n Hyperbolic eccentric anomaly = %g',F)
fprintf('\n----------------------------------------------\n')