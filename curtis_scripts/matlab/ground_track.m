% ALGORITHM 4.6: CALCULATE THE GROUND TRACK OF A SATELLITE FROM ITS ORBITAL ELEMENTS

function ground_track
%{
    This program plots the ground track of an earth satellite
    for which the orbital elements are specified

    mu          - gravitational parameter (km^3/s^2)
    deg         - factor that converts degrees to radians
    J2          - second zonal harmonic
    Re          - earth s radius (km)
    we          - earth s angular velocity (rad/s)
    rP          - perigee of orbit (km)
    rA          - apogee of orbit (km)
    TA, TAo     - true anomaly, initial true anomaly of satellite (rad)
    RA, RAo     - right ascension, initial right ascension of the node (rad)
    incl        - orbit inclination (rad)
    wp, wpo     - argument of perigee, initial argument of perigee (rad)
    n_periods   - number of periods for which ground track is to be plotted
    a           - semimajor axis of orbit (km)
    T           - period of orbit (s)
    e           - eccentricity of orbit
    h           - angular momentum of orbit (km^2/s)
    E, Eo       - eccentric anomaly, initial eccentric anomaly (rad)
    M, Mo       - mean anomaly, initial mean anomaly (rad)
    to, tf      - initial and final times for the ground track (s)
    fac         - common factor in Equations 4.53 and 4.53
    RAdot       - rate of regression of the node (rad/s)
    wpdot       - rate of advance of perigee (rad/s)
    times       - times at which ground track is plotted (s)
    ra          - vector of right ascensions of the spacecraft (deg)
    dec         - vector of declinations of the spacecraft (deg)
    TA          - true anomaly (rad)
    r           - perifocal position vector of satellite (km)
    R           - geocentric equatorial position vector (km)
    R1          - DCM for rotation about z through RA
    R2          - DCM for rotation about x through incl
    R3          - DCM for rotation about z through wp
    QxX         - DCM for rotation from perifocal to geocentric equatorial
    Q           - DCM for rotation from geocentric equatorial
                  into earth-fixed frame
    r_rel       - position vector in earth-fixed frame (km)
    alpha       - satellite right ascension (deg)
    delta       - satellite declination (deg)
    n_curves    - number of curves comprising the ground track plot
    RA          - cell array containing the right ascensions for each of
                  the curves comprising the ground track plot
    Dec         - cell array containing the declinations for each of
                  the curves comprising the ground track plot

    User M-functions required: sv_from_coe, kepler_E, ra_and_dec_from_r
%}
% 
clear all; close all; clc
global ra dec n_curves RA Dec

%...Constants
deg       = pi/180;
mu        = 398600;
J2        = 0.00108263;
Re        = 6378;
we        = (2*pi + 2*pi/365.26)/(24*3600);

%...Data declaration for Example 4.12:
rP        = 6700;
rA        = 10000;
TAo       = 230*deg;
Wo        = 270*deg;
incl      = 60*deg;
wpo       = 45*deg;
n_periods = 3.25;
%...End data declaration

%...Compute the initial time (since perigee) and
% the rates of node regression and perigee advance
a         = (rA + rP)/2;
T         = 2*pi/sqrt(mu)*a^(3/2);
e         = (rA - rP)/(rA + rP);
h         = sqrt(mu*a*(1 - e^2));
Eo        = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
Mo        = Eo - e*sin(Eo);
to        = Mo*(T/2/pi);
tf        = to + n_periods*T;
fac       = -3/2*sqrt(mu)*J2*Re^2/(1-e^2)^2/a^(7/2);
Wdot      = fac*cos(incl);
wpdot     = fac*(5/2*sin(incl)^2 - 2);

find_ra_and_dec
form_separate_curves
plot_ground_track
print_orbital_data

return

% ----------------------
function find_ra_and_dec
% ----------------------
% Propagates the orbit over the specified time interval, transforming
% the position vector into the earth-fixed frame and, from that,
% computing the right ascension and declination histories
% ----------------------
%
times = linspace(to,tf,1000);
ra    = [];
dec   = [];
theta = 0;
for i               = 1:length(times)
    t               = times(i);
    M               = 2*pi/T*t;
    E               = kepler_E(e, M);
    TA              = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    r               = h^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]';

    W               = Wo + Wdot*t;
    wp              = wpo + wpdot*t;

    R1              = [ cos(W) sin(W) 0
                       -sin(W) cos(W) 0
                           0      0   1];

    R2              = [1     0         0
                       0  cos(incl) sin(incl)
                       0 -sin(incl) cos(incl)];

    R3              = [ cos(wp) sin(wp) 0
                       -sin(wp) cos(wp) 0
                           0       0    1];

    QxX = (R3*R2*R1)';
    R = QxX*r;

    theta           = we*(t - to);
    Q               = [ cos(theta) sin(theta) 0
                       -sin(theta) cos(theta) 0
                            0         0       1];
    r_rel           = Q*R;

    [alpha delta]   = ra_and_dec_from_r(r_rel);

    ra              = [ra; alpha];
    dec             = [dec; delta];
end

end %find_ra_and_dec

% ---------------------------
function form_separate_curves
% ---------------------------
% Breaks the ground track up into separate curves which start
% and terminate at right ascensions in the range [0,360 deg].
% ---------------------------
tol = 100;
curve_no = 1;
n_curves = 1;
k        = 0;
ra_prev  = ra(1);
for i = 1:length(ra)
    if abs(ra(i) - ra_prev) > tol
        curve_no = curve_no + 1;
        n_curves = n_curves + 1;
        k = 0;
    end
    k                = k + 1;
    RA{curve_no}(k)  = ra(i);
    Dec{curve_no}(k) = dec(i);
    ra_prev          = ra(i);
end

end %form_separate_curves

% -------------------
function plot_ground_track
% -------------------
hold on
xlabel('East longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
grid on
for i = 1:n_curves
    plot(RA{i}, Dec{i})
end

axis ([0 360 -90 90])
text( ra(1), dec(1), 'o Start')
text(ra(end), dec(end), 'o Finish')
line([min(ra) max(ra)],[0 0], 'Color','k') % the equator

end %plot_ground_track

    % -------------------------
    function print_orbital_data
    % -------------------------
    coe = [h e Wo incl wpo TAo];
    [ro, vo] = sv_from_coe(coe, mu);
    fprintf('\n ----------------------------------------------\n')
    fprintf('\n Angular momentum     = %g km^2/s' , h)
    fprintf('\n Eccentricity         = %g' , e)
    fprintf('\n Semimajor axis       = %g km' , a)
    fprintf('\n Perigee radius       = %g km' , rP)
    fprintf('\n Apogee radius        = %g km' , rA)
    fprintf('\n Period               = %g hours' , T/3600)
    fprintf('\n Inclination          = %g deg' , incl/deg)
    fprintf('\n Initial true anomaly = %g deg' , TAo/deg)
    fprintf('\n Time since perigee   = %g hours' , to/3600)
    fprintf('\n Initial RA           = %g deg' , Wo/deg)
    fprintf('\n RA_dot               = %g deg/period' , Wdot/deg*T)
    fprintf('\n Initial wp           = %g deg' , wpo/deg)
    fprintf('\n wp_dot               = %g deg/period' , wpdot/deg*T)
    fprintf('\n')
    fprintf('\n r0 = [%12g, %12g, %12g] (km)', ro(1), ro(2), ro(3))
    fprintf('\n magnitude = %g km\n', norm(ro))
    fprintf('\n v0 = [%12g, %12g, %12g] (km)', vo(1), vo(2), vo(3))
    fprintf('\n magnitude = %g km\n', norm(vo))
    fprintf('\n ----------------------------------------------\n')
    end %print_orbital_data

end %ground_track