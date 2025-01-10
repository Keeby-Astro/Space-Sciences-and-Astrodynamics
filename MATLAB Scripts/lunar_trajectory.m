%lunar_trajectory
%{
    This program presents the graphical solution of the motion of a
    spacecraft in the gravity fields of both the earth and the moon for
    the initial data provided in the input declaration below.
    MATLAB's ode45 Runge-Kutta solver is used.

    deg                         - conversion factor, degrees to radians
    days                        - conversion factor, days to seconds
    Re, Rm                      - radii of earth and moon, respectively (km)
    m_e, m_m                    - masses of earth and moon, respectively (kg)
    mu_e, mu_m                  - gravitational parameters of earth and moon,
                                  respectively (km^3/s^2)
    D                           - semimajor axis of moon's orbit (km)
    I_, J_, K_                  - unit vectors of ECI frame
    RS                          - radius of moon's sphere of influence (km)
                                  year, month, hour, minute
    second                      - Date and time of spacecraft's lunar arrival
    t0                          - initial time on the trajectory (s)
    z0                          - initial altitude of the trajectory (km)
    alpha0, dec0                - initial right ascension and declination of
                                  spacecraft (deg)
    gamma0                      - initial flight path angle (deg)
    fac                         - ratio of spaccraft's initial speed to the
                                  escape speed.
    ttt                         - predicted time to perilune (s)
    tf                          - time at end of trajectory (s)
    jd0                         - julian date of lunar arrival
    rm0, vm0                    - state vector of the moon at jd0 (km, km/s)
    RA, Dec                     - right ascension and declination of the moon
                                  at jd0 (deg)
    hmoon_, hmoon               - moon's angular momentum vector and magnitude
                                  at jd0 (km^2/s)
    inclmoon                    - inclination of moon's orbit earth's
                                  equatorial plane (deg)
    r0                          - initial radius from earth's center to
                                  probe (km)
    r0_                         - initial ECI position vector of probe (km)
    vesc                        - escape speed at r0 (km/s)
    v0                          - initial ECI speed of probe (km/s)
    w0_                         - unit vector normal to plane of translunar
                                  orbit at time t0
    ur_                         - radial unit vector to probe at time t0
    uperp_                      - transverse unit vector at time t0
    vr                          - initial radial speed of probe (km/s)
    vperp                       - initial transverse speed of probe (km/s)
    v0_                         - initial velocity vector of probe (km/s)
    uv0_                        - initial tangential unit vector
    y0                          - initial state vector of the probe (km, km/s)
    t                           - vector containing the times from t0 to tf at
                                  which the state vector is evaluated (s)
    y                           - a matrix whose 6 columns contain the inertial
                                  position and velocity components evaluated
                                  at the times t(:) (km, km/s)
    X, Y, Z                     - the probe's inertial position vector history
    vX, vY, VZ                  - the probe's inertial velocity history
    x, y, z                     - the probe's position vector history in the
                                  Moon-fixed frame
    Xm, Ym, Zm                  - the Moon's inertial position vector history
    vXm, vYm, vZm               - the Moon's inertial velocity vector history
    ti                          - the ith time of the set [t0,tf] (s)
    r_                          - probe's inertial position vector at time ti
                                  (km)
    r                           - magnitude of r_ (km)
    jd                          - julian date of corresponding to ti (days)
    rm_, vm_                    - the moon's state vector at time ti (km,km/s)
    x_, y_, z_                  - vectors along the axes of the rotating
                                  moon-fixed at time ti (km)
    i_, j_, k_                  - unit vectors of the moon-fixed rotating frame
                                  at time ti
    Q                           - DCM of transformation from ECI to moon-fixed
                                  frame at time ti
    rx_                         - probe's inertial position vector in moon-
                                  fixed coordinates at time ti (km)
    rmx_                        - Moon's inertial position vector in moon-
                                  fixed coordinates at time ti (km)
    dist_                       - position vector of probe relative to the moon
                                  at time ti (km)
    dist                        - magnitude of dist_ (km)
    dist_min                    - perilune of trajectory (km)
    rmTLI_                      - Moon's position vector at TLI
    RATLI, DecTLI               - Moon's right ascension and declination at
                                  TKI (deg)
    v_atdmin_                   - Probe's velocity vector at perilune (km/s)
    rm_perilume, vm_perilune    - Moon's state vector when the probe is at
                                  perilune (km, km/s)
    rel_speed                   - Speed of probe relative to the Moon at
                                  perilune (km/s)
    RA_at_perilune              - Moon's RA at perilune arrival (deg)
    Dec_at_perilune             - Moon's Dec at perilune arrival (deg)
    target_error                - Distance between Moon's actual position at
                                  perilune arrival and its position after the
                                  predicted flight time, ttt (km).
    rms_                        - position vector of moon relative to
                                  spacecraft (km)
    rms                         - magnitude of rms_ (km)
    aearth_                     - acceleration of spacecraft due to
                                  earth (km/s^2)
    amoon_                      - acceleration of spacecraft due to
                                  moon (km/s^2)
    atot_                       - aearth_ + amoon_ (km/s^2)
    binormal_                   - unit vector normal to the osculating plane
    incl                        - angle between inertial Z axis and the
                                  binormal (deg)
    rend_                       - Position vector of end point of trajectory
                                  (km)
    alt_end                     - Altitude of end point of trajectory (km)
    ra_end, dec_end             - Right ascension and declination of end point
                                  of trajectory (km)

    User M-functions required: none
    User subfunctions required: rates, plotit_XYZ, plotit_xyz
%}
% --------------------------------------------------------------------

clear all; close all; clc
fprintf('\nlunar_restricted_threebody.m %s\n\n', datestr(now))

global jd0 days ttt mu_m mu_e Re Rm rm_ rm0_

%...general data
deg  = pi/180;
days = 24*3600;
Re   = 6378;
Rm   = 1737;
m_e  = 5974.e21;
m_m  = 73.48e21;
mu_e = 398600.4;
mu_m = 4902.8;
D    = 384400;
RS   = D*(m_m/m_e)^(2/5);
%...

%...Data declaration
Title = 'Lunar Trajectory';
%   Date and time of lunar arrival:
  year   = 2020;
  month  = 5;
  day    = 4;
  hour   = 12;
  minute = 0;
  second = 0;
t0       = 0;
z0       = 320;
alpha0   = 90;
dec0     = 15;
gamma0   = 40;
fac      = .9924; %Fraction of Vesc
ttt      = 3*days;
tf       = ttt + 2.667*days;
%...End data declaration

%...State vector of moon at target date:
jd0         = juliandate(year, month, day, hour, minute, second);
[rm0_,vm0_] = simpsons_lunar_ephemeris(jd0);
%[rm0_,vm0] = planetEphemeris(jd0, 'Earth', 'Moon', '430');
[RA, Dec]   = ra_and_dec_from_r(rm0_);
distance    = norm(rm0_);
hmoon_      = cross(rm0_,vm0_);
hmoon       = norm(hmoon_);
inclmoon    = acosd(hmoon_(3)/hmoon);

%...Initial position vector of probe:
I_      = [1;0;0];
J_      = [0;1;0];
K_      = cross(I_,J_);
r0      = Re + z0;
r0_     = r0*(cosd(alpha0)*cosd(dec0)*I_ + ...
              sind(alpha0)*cosd(dec0)*J_ + ...
              sind(dec0)*K_);
vesc    = sqrt(2*mu_e/r0);
v0      = fac*vesc;
w0_     = cross(r0_,rm0_)/norm(cross(r0_,rm0_));

%...Initial velocity vector of probe:
ur_     = r0_/norm(r0_);
uperp_  = cross(w0_,ur_)/norm(cross(w0_,ur_));
vr      = v0*sind(gamma0);
vperp   = v0*cosd(gamma0);
v0_     = vr*ur_ + vperp*uperp_;
uv0_    = v0_/v0;

%...Initial state vector of the probe:
y0      = [r0_(1) r0_(2) r0_(3) v0_(1) v0_(2) v0_(3)]';

%...Pass the initial conditions and time interval to ode45, which
%   calculates the position and velocity of the spacecraft at discrete
%   times t, returning the solution in the column vector y. ode45 uses
%   the subfunction 'rates' below to evaluate the spacecraft acceleration
%   at each integration time step.

options = odeset('RelTol', 1.e-10, 'AbsTol', 1.e-10,'Stats', 'off');
[t,y]   = ode45(@rates, [t0 tf], y0, options);

%...Spacecraft trajectory
% in ECI frame:
X   = y(:,1); Y   = y(:,2); Z   = y(:,3);
vX  = y(:,4); vY  = y(:,5); vZ  = y(:,6);
% in Moon-fixed frame:
x   = [];     y   = [];     z   = [];

%...Moon trajectory
% in ECI frame:
Xm  = [];     Ym  = [];     Zm  = [];
vXm = [];     vYm = [];     vZm = [];
% in Moon-fixed frame:
xm  = [];     ym  = [];     zm  = [];

%...Compute the Moon's trajectory from an ephemeris, find perilune of the
%   probe's trajectory, and project the probe's trajectory onto the axes
%   of the Moon-fixed rotating frame:

dist_min = 1.e30; %Starting value in the search for perilune
for i = 1:length(t)
    ti = t(i);

    %...Probe's inertial position vector at time ti:
    r_ = [X(i) Y(i) Z(i)]';

    %...Moon's inertial position and velocity vectors at time ti:
    jd        = jd0 - (ttt - ti)/days;
    [rm_,vm_] = simpsons_lunar_ephemeris(jd);

    %...Moon's inertial state vector at time ti:
    Xm  = [ Xm;rm_(1)]; Ym  = [ Ym;rm_(2)]; Zm  = [ Zm;rm_(3)];
    vXm = [vXm;vm_(1)]; vYm = [vYm;vm_(2)]; vZm = [vZm;vm_(3)];

    %...Moon-fixed rotating xyz frame:
    x_ = rm_;
    z_ = cross(x_,vm_);
    y_ = cross(z_,x_);
    i_ = x_/norm(x_);
    j_ = y_/norm(y_);
    k_ = z_/norm(z_);

    %...DCM of transformation from ECI to moon-fixed frame:
    Q = [i_'; j_'; k_'];

    %...Components of probe's inertial position vector in moon-fixed frame:
    rx_ = Q*r_;
    x   = [x;rx_(1)]; y = [y;rx_(2)]; z = [z;rx_(3)];

    %...Components of moon's inertial position vector in moon-fixed frame:
    rmx_ = Q*rm_;
    xm   = [xm;rmx_(1)]; ym = [ym;rmx_(2)]; zm = [zm;rmx_(3)];

    %...Find perilune of the probe:
    dist_ = r_ - rm_;
    dist  = norm(dist_);
    if dist < dist_min
        imin     = i;
        ist_min_ = dist_;
        dist_min = dist;
    end
end

%...Location of the Moon at TLI:
rmTLI_          = [Xm(1); Ym(1); Zm(1)];
[RATLI, DecTLI] = ra_and_dec_from_r(rmTLI_);

%...Spacecraft velocity at perilune:
v_atdmin_ = [vX(imin);vY(imin);vZ(imin)];

%...State vector and celestial position of moon when probe is at perilune:
rm_perilune_                      = [Xm(imin) Ym(imin) Zm(imin)]';
vm_perilune_                      = [vXm(imin) vYm(imin) vZm(imin)]';
[RA_at_perilune, Dec_at_perilune] = ra_and_dec_from_r(rm_perilune_);
target_error                      = norm(rm_perilune_ - rm0_);

%...Speed of probe relative to Moon at perilune:
rel_speed = norm(v_atdmin_ - vm_perilune_);

%...End point of trajectory:
rend_             = [X(end); Y(end); Z(end)];
alt_end           = norm(rend_) - Re;
[ra_end, dec_end] = ra_and_dec_from_r(rend_);

%...Find the history of the trajectory's binormal:
for i = 1:imin
    time(i) = t(i);
    r_          = [X(i) Y(i) Z(i)]';
    r           = norm(r_);
    v_          = [vX(i) vY(i) vZ(i)]';
    rm_         = [Xm(i) Ym(i) Zm(i)]';
    rm          = norm(rm_);
    rms_        = rm_ - r_;
    rms(i)      = norm(rms_);
    aearth_     = -mu_e*r_/r^3;
    amoon_      = mu_m*(rms_/rms(i)^3 - rm_/rm^3);
    atot_       = aearth_ + amoon_;
    binormal_   = cross(v_,atot_)/norm(cross(v_,atot_));
    binormalz   = binormal_(3);
    incl(i)     = acosd(binormalz);
end

%...Output:
fprintf('\n\n%s\n\n', Title)
fprintf('Date and time of arrival at moon: ')
fprintf('%s/%s/%s %s:%s:%s', ...
         num2str(month), num2str(day), num2str(year), ...
         num2str(hour), num2str(minute), num2str(second))
fprintf('\nMoon''s position: ')
fprintf('\n Distance                        = %11g km'      , distance)
fprintf('\n Right Ascension                 = %11g deg'     , RA)
fprintf('\n Declination                     = %11g deg '    , Dec)
fprintf('\nMoon''s orbital inclination      = %11g deg\n'   , inclmoon)

fprintf('\nThe probe at earth departure (t = %g sec):'      , t0)
fprintf('\n Altitude                        = %11g km'      , z0)
fprintf('\n Right ascension                 = %11g deg'     , alpha0)
fprintf('\n Declination                     = %11g deg'     , dec0)
fprintf('\n Flight path angle               = %11g deg'     , gamma0)
fprintf('\n Speed                           = %11g km/s'    , v0)
fprintf('\n Escape speed                    = %11g km/s'    , vesc)
fprintf('\n v/vesc                          = %11g'         , v0/vesc)
fprintf('\n Inclination of translunar orbit = %11g deg\n'   , ...
                                                           acosd(w0_(3)))

fprintf('\nThe moon when the probe is at TLI:')
fprintf('\n Distance                        = %11g km'  , norm(rmTLI_))
fprintf('\n Right ascension                 = %11g deg' , RATLI)
fprintf('\n Declination                     = %11g deg' , DecTLI)

fprintf('\nThe moon when the probe is at perilune: ')
fprintf('\n Distance                        = %11g km'  , ...
norm(rm_perilune_))
fprintf('\n Speed                           = %11g km/s', ...
norm(vm_perilune_))
fprintf('\n Right ascension                 = %11g deg' ,RA_at_perilune)
fprintf('\n Declination                     = %11g deg' ,Dec_at_perilune)
fprintf('\n Target error                    = %11g km'  , target_error)

fprintf('\n\nThe probe at perilune:')
fprintf('\n Altitude                        = %11g km'  , dist_min - Rm)
fprintf('\n Speed                           = %11g km/s',norm(v_atdmin_))
fprintf('\n Relative speed                  = %11g km/s',rel_speed)
fprintf('\n Inclination of osculating plane = %11g deg' ,incl(imin))
fprintf('\n Time from TLI to perilune       = %11g hours (%g days)', ...
                                                  abs(t(imin))/3600, ...
                                                  abs(t(imin))/3600/24)

fprintf('\n\nTotal time of flight           = %11g days'  , t(end)/days)
fprintf('\nTime to target point             = %11g days'  , ttt/days)
fprintf('\nFinal earth altitude             = %11g km'    , alt_end)
fprintf('\nFinal right ascension            = %11g deg'   , ra_end)
fprintf('\nFinal declination                = %11g deg\n' , dec_end)
%...End output

%...Graphical output"
%   Plot the trajectory relative to the inertial frame:
plotit_XYZ(X,Y,Z,Xm,Ym,Zm,imin)

%   Plot inclination of the osculating plane vs distance from the Moon
figure
hold on
plot(rms/RS,incl)
line([0 6],[90 90],'Linestyle','-','color','red')
title('Osculating Plane Inclination vs Distance from Moon')
xlabel('r_{ms}/R_s')
ylabel('Inclination deg)')
grid on
grid minor

%   Plot the trajectory relative to the rotating Moon-fixed frame:
plotit_xyz(x,y,z,xm,ym,zm,imin)
%...End graphical output

return

function dydt = rates(t,y)
%{
    This function evaluates the 3D acceleration of the spacecraft in a
    restricted 3-body system at time t from its position and velocity
    and the position of the moon at that time.

    t           - time (s)
    ttt         - flight time, TLI to target point (s)
    jd0         - Julian Date on arrival at target (days)
    jd          - Julian Date at time t (days)
    X, Y, Z     - Components of spacecraft's geocentric position vector (km)
    vX, vY, vZ  - Components of spacecraft's geocentric velocity vector (km/s)
    aX, aY, aZ  - Components of spacecraft's geocentric acceleration
                  vector (km/s^2)
    y           - column vector containing the geocentric position and
                  velocity components of the spacecraft at time t
    r_          - geocentric position vector [X Y Z] of the spacecraft
    rm_         - geocentric position vector of the moon
    rms_        - rm_ - r_, the position of the moon relative to the
                  spacecraft
    aearth_     - spacecraft acceleration vector due to earth's gravity
    amoon_      - spacecraft acceleration vector due to lunar gravity
    a_          - total spacececraft acceleration vector
    dydt        - column vector containing the geocentric velocity and
                  acceleration components of the spacecraft at time t
%}
% --------------------------------------------------------------------

global jd0 days mu_m mu_e ttt

jd      = jd0 - (ttt - t)/days;

X       = y(1);
Y       = y(2);
Z       = y(3);

vX      = y(4);
vY      = y(5);
vZ      = y(6);

r_      = [X Y Z]';
r       = norm(r_);
[rm_,]  = simpsons_lunar_ephemeris(jd);
%[rm_,] = planetEphemeris(jd0, 'Earth', 'Moon', '430');
rm      = norm(rm_);
rms_    = rm_ - r_;
rms     = norm(rms_);
aearth_ = -mu_e*r_/r^3;
amoon_  = mu_m*(rms_/rms^3 - rm_/rm^3);
a_      = aearth_ + amoon_;
aX      = a_(1);
aY      = a_(2);
aZ      = a_(3);

dydt    = [vX vY vZ aX aY aZ]';

end %rates

% --------------------------------------
function plotit_XYZ(X,Y,Z,Xm,Ym,Zm,imin)
% --------------------------------------
global Re Rm

figure ('Name','Trajectories of Spacecraft (red) and Moon (green)', ...
        'Color', [1 1 1])
[xx, yy, zz] = sphere(128);
hold on

%...Geocentric inertial coordinate axes:
L = 20*Re;
line([0 L], [0 0], [0 0],'color','k')
text(L,0,0, 'X', 'FontSize',12, 'FontAngle','italic','FontName','Palatino')
line([0 0], [0 L], [0 0],'color','k')
text(0,L,0, 'Y', 'FontSize',12, 'FontAngle','italic','FontName','Palatino')
line([0 0], [0 0], [0 L],'color','k')
text(0,0,L, 'Z', 'FontSize',12, 'FontAngle','italic','FontName','Palatino')

%...Earth:
Earth = surfl(Re*xx, Re*yy, Re*zz);
set(Earth, 'FaceAlpha', 0.5);
shading interp

%...Spacecraft at TLI
plot3(X(1), Y(1), Z(1), 'o', ...
'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',3)

%...Spacecraft at closest approach
plot3(X(imin), Y(imin), Z(imin), 'o', ...
'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',2)

%...Spacecraft at tf
plot3(X(end), Y(end), Z(end), 'o', ...
'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',3)

%...Moon at TLI:
text(Xm(1), Ym(1), Zm(1), 'Moon at TLI')
Moon = surfl(Rm*xx + Xm(1), Rm*yy + Ym(1), Rm*zz + Zm(1));
set(Moon, 'FaceAlpha', 0.99)
shading interp

%...Moon at closest approach:
Moon = surfl(Rm*xx + Xm(imin), Rm*yy + Ym(imin), Rm*zz + Zm(imin));
set(Moon, 'FaceAlpha', 0.99)
shading interp

%...Moon at end of simulation:
Moon = surfl(Rm*xx + Xm(end), Rm*yy + Ym(end), Rm*zz + Zm(end));
set(Moon, 'FaceAlpha', 0.99)
shading interp

%...Spacecraft trajectory
plot3( X, Y, Z, 'r', 'LineWidth', 1.5)

%...Moon trajectory
plot3(Xm, Ym, Zm, 'g', 'LineWidth', 0.5)

axis image
axis off
axis vis3d
view([1,1,1])

end %plotit_XYZ

% --------------------------------------
function plotit_xyz(x,y,z,xm,ym,zm,imin)
% --------------------------------------
global Re Rm Rm0_ Q0

figure ('Name','Spacecraft trajectory in Moon-fixed rotating frame', ...
        'Color', [1 1 1])
[xx, yy, zz] = sphere(128);
hold on

%...Spacecraft trajectory:
plot3( x, y, z, 'r', 'LineWidth', 2.0)

%...Moon trajectory:
plot3(xm, ym, zm, 'g', 'LineWidth', 0.5)

%...Earth:
Earth = surfl(Re*xx, Re*yy, Re*zz);
set(Earth, 'FaceAlpha', 0.5);
shading interp

%...Geocentric moon-fixed coordinate axes:
L1 = 63*Re; L2 = 20*Re; L3 = 29*Re;
line([0 L1], [0 0], [0 0],'color','k')
text(L1, 0, 0, 'x','FontSize', 12, 'FontAngle', 'italic','FontName','Palatino')
line([0 0], [0 L2], [0 0],'color','k')
text(0, L2, 0, 'y','FontSize', 12, 'FontAngle', 'italic','FontName','Palatino')
line([0 0], [0 0], [0 L3],'color','k')
text(0, 0, L3, 'z','FontSize', 12, 'FontAngle', 'italic','FontName','Palatino')

%...Spacecraft at TLI
plot3(x(1), y(1), z(1), 'o', ...
'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',3)

%...Spacecraft at closest approach
plot3(x(imin), y(imin), z(imin), 'o', ...
'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',2)

%...Spacecraft at tf
plot3(x(end), y(end), z(end), 'o', ...
'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',3)

%...Moon at TLI:
text(xm(1), ym(1), zm(1),'Moon at TLI')
Moon = surfl(Rm*xx + xm(1), Rm*yy + ym(1), Rm*zz + zm(1));
set(Moon, 'FaceAlpha', 0.99)
shading interp

%...Moon at spacecraft closest approach:
Moon = surfl(Rm*xx + xm(imin), Rm*yy + ym(imin), Rm*zz + zm(imin));
set(Moon, 'FaceAlpha', 0.99)
shading interp

%...Moon at end of simulation:
Moon = surfl(Rm*xx + xm(end), Rm*yy + ym(end), Rm*zz + zm(end));
set(Moon, 'FaceAlpha', 0.99)
shading interp

axis image
axis vis3d
axis off
view([1,1,1])
end %plotit_xyz

%end lunar_trajectory
