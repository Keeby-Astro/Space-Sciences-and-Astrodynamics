function spinning_top
%{
  This program numerically integrates Euler's equations of motion
  for the spinning top. The  quaternion is used to obtain the time history
  of the top's orientation.
 
  User M-functions required: rkf45, q_from_dcm, dcm_from_q, dcm_to_euler
  User subfunction required: rates
%}
%----------------------------------------------
clear all; close all; clc
 
%...Data declaration:
g     = 9.807;         % Acceleration of gravity (m/s^2)
m     = 0.5;           % Mass in kg
d     = 0.05;          % Distance of center of mass from pivot point (m)  
A     = 12.e-4;        % Moment of inertia about body x (kg-m^2)
B     = 12.e-4;        % Moment of inertia about body y (kg-m^2)
C     = 4.5e-4;        % Moment of inertia about body z (kg-m^2)
ws0   = 1000*2*pi/60;  % Spin rate (rad/s)
 
wp0   = 51.93*2*pi/60; % Precession rate (rad/s) Use to obtain Fig. 11.33
wp0   = 0;             %                         Use to obtain Fig, 11.34
 
wn0   = 0;             % Nutation rate (deg/s)
theta = 60;            % Initial nutation angle (deg)
z     = [0 -sind(theta) cosd(theta)];  % Initial z-axis direction:
p     = [1 0 0];                       % Initial x-axis direction
                                       %   (or a line defining x-z plane)
%...
y     = cross(z,p);     % y-axis direction (normal to x-z plane)
x     = cross(y,z);     % x-axis direction (normal to y-z plane)
i     = x/norm(x);      % Unit vector along x axis
j     = y/norm(y);      % Unit vector along y axis
k     = z/norm(z);      % Unit vector along z axis
QXx   = [i; j; k];      % Initial direction cosine matrix
 
%...Initial precession, nutation, and spin angles (deg):
[phi0 theta0 psi0]  = dcm_to_euler(QXx); 
 
%...Initial quaternion (column vector):
q0    = q_from_dcm(QXx);
 
%...Initial body-frame angular velocity, column vector (rad/s):
w0    = [wp0*sind(theta0)*sin(psi0) + wn0*cosd(psi0), ...
         wp0*sind(theta0)*cos(psi0) - wn0*sind(psi0), ...
         ws0 + wp0*cosd(theta0)]';
   
t0    = 0;              % Initial time (s)
tf    = 1.153;          % Final time (s)  (for 360 degrees of precession)
f0    = [q0; w0];       % Initial conditions vector (quaternion & angular velocities)
                     
%...RKF4(5) numerical ODE solver. Time derivatives computed in
%   function 'rates' below.
[t,f] = rkf45(@rates, [t0,tf], f0); 
                                    
%...Solutions for quaternion and angular velocities at 'nsteps' times
%   from t0 to tf 
q     = f(:,1:4); 
wx    = f(:,5); 
wy    = f(:,6);
wz    = f(:,7);
 
%...Obtain the direction cosine matrix, the Euler angles and the Euler 
%   angle rates at each solution time:
for m = 1:length(t)
    %...DCM from the quaternion:
    QXx       = dcm_from_q(q(m,:));
    %...Euler angles (deg) from DCM:
    [prec(m) ...
     nut(m)  ...
     spin(m)] = dcm_to_euler(QXx);
    %...Euler rates from Eqs. 11.116:
    wp(m)     = (wx(m)*sind(spin(m)) + wy(m)*cosd(spin(m)))/sind(nut(m));
    wn(m)     =  wx(m)*cosd(spin(m)) - wy(m)*sind(spin(m));
    ws(m)     = -wp(m)*cosd(nut(m)) + wz(m);
end
 
plotit

function dfdt = rates(t,f)
%
q      = f(1:4);          % components of quaternion
wx     = f(5);            % angular velocity along x
wy     = f(6);            % angular velocity along y
wz     = f(7);            % angular velocity along z
 
q      = q/norm(q);       % normalize the quaternion
 
Q      = dcm_from_q(q);   % DCM from quaternion
 
%...Body frame components of the moment of the weight vector
%   about the pivot point:
M      = Q*[-m*g*d*Q(3,2)
             m*g*d*Q(3,1)
                       0];
%...Skew-symmetric matrix of angular velocities:
Omega  = [  0   wz  -wy   wx
          -wz    0   wx   wy
           wy  -wx    0   wz
          -wx  -wy  -wz    0];                              
q_dot  = Omega*q/2;                % time derivative of quaternion
 
%...Euler's equations:  
wx_dot = M(1)/A - (C - B)*wy*wz/A; % time derivative of wx
wy_dot = M(2)/B - (A - C)*wz*wx/B; % time derivative of wy
wz_dot = M(3)/C - (B - A)*wx*wy/C; % time derivative of wz
 
%...Return the rates in a column vector:
dfdt   = [q_dot; wx_dot; wy_dot; wz_dot];
 
end %rates
 
function plotit
 
figure('Name', 'Euler angles and their rates', 'color', [1 1 1])
 
subplot(321)
plot(t, prec )
xlabel('time (s)')
ylabel('Precession angle (deg)')
axis([-inf, inf, -inf, inf])
axis([-inf, inf, -inf, inf])
grid
 
subplot(322)
plot(t, wp*60/2/pi)
xlabel('time (s)')
ylabel('Precession rate (rpm)')
axis([0, 1.153, 51, 53])
axis([-inf, inf, -inf, inf])
grid
 
subplot(323)
plot(t, nut)
xlabel('time (s)')
ylabel('Nutation angle (deg)')
axis([0, 1.153, 59, 61])
axis([-inf, inf, -inf, inf])
grid
 
subplot(324)
plot(t, wn*180/pi)
xlabel('time (s)')
ylabel('Nutation rate (deg/s)')
axis([-inf, inf, -inf, inf])
grid
 
subplot(325)
plot(t, spin)
xlabel('time (s)')
ylabel('Spin angle (deg)')
 
axis([-inf, inf, -inf, inf])
grid
 
subplot(326)
plot(t, ws*60/2/pi)
xlabel('time (s)')
ylabel('Spin rate (rpm)')
axis([-inf, inf, -inf, inf])
grid
 
end %plotit
 
end %spinning_top
