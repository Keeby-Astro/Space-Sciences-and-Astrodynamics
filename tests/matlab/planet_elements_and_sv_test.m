%planet_elements_and_sv_test
%
% This program uses Algorithm 8.1 to compute the orbital elements
% and state vector of the earth at the date and time specified.
% To obtain the same results for Mars, set planet_id = 4.
%
% mu        - gravitational parameter of the sun (km^3/s^2)
% deg       - conversion factor between degrees and radians
% pi        - 3.1415926...
%
% coe       - vector of heliocentric orbital elements
%             [h  e  RA  incl  w  TA  a  w_hat  L  M  E],
%             where
%              h     = angular momentum                    (km^2/s)
%              e     = eccentricity
%              RA    = right ascension                     (deg)
%              incl  = inclination                         (deg)
%              w     = argument of perihelion              (deg)
%              TA    = true anomaly                        (deg)
%              a     = semimajor axis                      (km)
%              w_hat = longitude of perihelion ( = RA + w) (deg)
%              L     = mean longitude ( = w_hat + M)       (deg)
%              M     = mean anomaly                        (deg)
%              E     = eccentric anomaly                   (deg)
%
% r         - heliocentric position vector (km)
% v         - heliocentric velocity vector (km/s) 
%
% planet_id - planet identifier:
%              1 = Mercury
%              2 = Venus
%              3 = Earth
%              4 = Mars
%              5 = Jupiter
%              6 = Saturn
%              7 = Uranus
%              8 = Neptune
%              9 = Pluto
%
% year      - range: 1901 - 2099
% month     - range: 1 - 12
% day       - range: 1 - 31
% hour      - range: 0 - 23
% minute    - range: 0 - 60
% second    - range: 0 - 60 
%
% User M-functions required: planet_elements_and_sv, month_planet_names
% --------------------------------------------------------------------

global mu
mu  = 1.327124e11;
deg = pi/180;

%...Input data
planet_id = 3;
year      = 2003;
month     = 8;
day       = 27;
hour      = 12;
minute    = 0;
second    = 0;
%...

%...Algorithm 8.1:
[coe, r, v, jd] = planet_elements_and_sv ...
              (planet_id, year, month, day, hour, minute, second);

%...Convert the planet_id and month numbers into names for output:
[month_name, planet_name] = month_planet_names(month, planet_id);

%...Echo the input data and output the solution to
%   the command window:
fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n');
fprintf('\n   Planet: %s', planet_name)
fprintf('\n   Year  : %g', year)
fprintf('\n   Month : %s', month_name)
fprintf('\n   Day   : %g', day)
fprintf('\n   Hour  : %g', hour)
fprintf('\n   Minute: %g', minute)
fprintf('\n   Second: %g', second)
fprintf('\n\n   Julian day: %11.3f', jd)

fprintf('\n\n');
fprintf(' Orbital elements:')
fprintf('\n');

fprintf('\n  Angular momentum (km^2/s)                   = %g', coe(1));
fprintf('\n  Eccentricity                                = %g', coe(2));
fprintf('\n  Right ascension of the ascending node (deg) = %g', coe(3));
fprintf('\n  Inclination to the ecliptic (deg)           = %g', coe(4));
fprintf('\n  Argument of perihelion (deg)                = %g', coe(5));
fprintf('\n  True anomaly (deg)                          = %g', coe(6));
fprintf('\n  Semimajor axis (km)                         = %g', coe(7));

fprintf('\n');

fprintf('\n  Longitude of perihelion (deg)               = %g', coe(8));
fprintf('\n  Mean longitude (deg)                        = %g', coe(9));
fprintf('\n  Mean anomaly (deg)                          = %g', coe(10));
fprintf('\n  Eccentric anomaly (deg)                     = %g', coe(11));

fprintf('\n\n');
fprintf(' State vector:')
fprintf('\n');

fprintf('\n  Position vector (km) = [%g  %g  %g]', r(1), r(2), r(3))
fprintf('\n  Magnitude            = %g\n', norm(r))
fprintf('\n  Velocity (km/s)      = [%g  %g  %g]', v(1), v(2), v(3))
fprintf('\n  Magnitude            = %g', norm(v))

fprintf('\n-----------------------------------------------------\n')