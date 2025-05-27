function gravity_turn
%{
  This program numerically integrates Equations 13.6 through
  13.8 for a gravity turn trajectory.
 
  M-functions required:      atmosisa
  User M-functions required: rkf45
  User subfunction requred:  rates
%}
% ----------------------------------------------
clear all;close all;clc
 
deg    = pi/180;        % Convert degrees to radians
g0     = 9.81;          % Sea-level acceleration of gravity (m/s)
Re     = 6378e3;        % Radius of the earth (m)
hscale = 7.5e3;         % Density scale height (m)
rho0   = 1.225;         % Sea level density of atmosphere (kg/m^3)
 
diam   = 196.85/12 ...
         *0.3048;       % Vehicle diameter (m)
A      = pi/4*(diam)^2; % Frontal area (m^2)
CD     = 0.5;           % Drag coefficient (assumed constant)
m0     = 149912*.4536;  % Lift-off mass (kg)
n      = 7;             % Mass ratio
T2W    = 1.4;           % Thrust to weight ratio
Isp    = 390;           % Specific impulse (s)
 
mfinal = m0/n;          % Burnout mass (kg)
Thrust = T2W*m0*g0;     % Rocket thrust (N)
m_dot  = Thrust/Isp/g0; % Propellant mass flow rate (kg/s)
mprop  = m0 - mfinal;   % Propellant mass (kg)
tburn  = mprop/m_dot;   % Burn time (s)
hturn  = 130;           % Height at which pitchover begins (m)
 
t0     = 0;             % Initial time for the numerical integration
tf     = tburn;         % Final time for the numerical integration
tspan  = [t0,tf];       % Range of integration
 
% ...Initial conditions:
v0     = 0;             % Initial velocity (m/s)
gamma0 = 89.85*deg;     % Initial flight path angle (rad)
x0     = 0;             % Initial downrange distance (km)
h0     = 0;             % Initial altitude (km)
vD0    = 0;             % Initial value of velocity loss due
                        % to drag (m/s)
vG0    = 0;             % Initial value of velocity loss due
                        % to gravity (m/s)
 
%...Initial conditions vector:
f0 = [v0; gamma0; x0; h0; vD0; vG0];
 
%...Call to Runge-Kutta numerical integrator 'rkf45'
%   rkf45 solves the system of equations df/dt = f(t):
 
[t,f] = rkf45(@rates, tspan, f0);
   
% ...Solution f(t) returned on the time interval [t0 tf]: 
v      =  f(:,1)*1.e-3;  % Velocity (km/s)
gamma  =  f(:,2)/deg;    % Flight path angle (degrees)
x      =  f(:,3)*1.e-3;  % Downrange distance (km)
h      =  f(:,4)*1.e-3;  % Altitude (km)
vD     = -f(:,5)*1.e-3;  % Velocity loss due to drag (km/s)
vG     = -f(:,6)*1.e-3;  % Velocity loss due to gravity (km/s)
 
%...Dynamic pressure vs time:
for i = 1:length(t)
    Rho  = rho0 * exp(-h(i)*1000/hscale);     % Air density (kg/m^3)
    q(i) = 1/2*Rho*(v(i)*1.e3)^2;             % Dynamic pressure (Pa)
    [dum a(i) dum dum] = atmosisa(h(i)*1000); % Speed of sound (m/s)
    M(i) = 1000*v(i)/a(i);                    % Mach number
end
 
%...Maximum dynamic pressure and corresponding time, speed, altitude and
%   Mach number:
[maxQ,imax]      = max(q);                 % qMax
tQ               = t(imax);                % Time
vQ               = v(imax);                % Speed
hQ               = h(imax);                % Altitude
[dum aQ dum dum] = atmosisa(h(imax)*1000); % Speed of sound at altitude
MQ               = 1000*vQ/aQ;
 
output
 
return
  
function dydt = rates(t,y)
%
% Calculates the time rates dy/dt of the variables y(t) 
% in the equations of motion of a gravity turn trajectory.
%-------------------------
 
%...Initialize dydt as a column vector:
dydt  = zeros(6,1);
 
v     = y(1);                      % Velocity
gamma = y(2);                      % Flight path angle
x     = y(3);                      % Downrange distance
h     = y(4);                      % Altitude
vD    = y(5);                      % Velocity loss due to drag
vG    = y(6);                      % Velocity loss due to gravity
 
%...When time t exceeds the burn time, set the thrust
%   and the mass flow rate equal to zero:
if t < tburn
    m = m0 - m_dot*t;              % Current vehicle mass
    T = Thrust;                    % Current thrust
else
    m = m0 - m_dot*tburn;          % Current vehicle mass
    T = 0;                         % Current thrust
end
    
g     = g0/(1 + h/Re)^2;           % Gravitational variation
                                   % with altitude h
rho   = rho0*exp(-h/hscale);       % Exponential density variation
                                   % with altitude
D     = 0.5*rho*v^2*A*CD;          % Drag [Equation 13.1]
 
%...Define the first derivatives of v, gamma, x, h, vD and vG
%   ("dot" means time derivative):
%v_dot = T/m - D/m - g*sin(gamma); % Equation 13.6
 
%...Start the gravity turn when h = hturn:
if h <= hturn
    gamma_dot = 0;
    v_dot     = T/m - D/m - g;
    x_dot     = 0;
    h_dot     = v;
    vG_dot    = -g;
else
    v_dot = T/m - D/m - g*sin(gamma);
    gamma_dot = -1/v*(g - v^2/(Re + h))*cos(gamma);% Equation 13.7
    x_dot    = Re/(Re + h)*v*cos(gamma);           % Equation 13.8(1)
    h_dot    = v*sin(gamma);                       % Equation 13.8(2)
    vG_dot   = -g*sin(gamma);                      % Gravity loss rate
end                                                % Equation 13.27(1)
 
vD_dot   = -D/m;                                   % Drag loss rate
                                                   % Equation 13.27(2)
%...Load the first derivatives of y(t) into the vector dydt:
dydt(1)  = v_dot;
dydt(2)  = gamma_dot;
dydt(3)  = x_dot;
dydt(4)  = h_dot;
dydt(5)  = vD_dot;
dydt(6)  = vG_dot;
end %rates
 
function output
%
fprintf('\n\n -----------------------------------\n')
fprintf('\n Initial flight path angle = %10.3f deg ',gamma0/deg)
fprintf('\n Pitchover altitude        = %10.3f m   ',hturn)
fprintf('\n Burn time                 = %10.3f s   ',tburn)
fprintf('\n Maximum dynamic pressure  = %10.3f atm ',maxQ*9.869e-6)
fprintf('\n    Time                   = %10.3f min ',tQ/60)
fprintf('\n    Speed                  = %10.3f km/s',vQ)
fprintf('\n    Altitude               = %10.3f km  ',hQ)
fprintf('\n    Mach Number            = %10.3f     ',MQ)
fprintf('\n At burnout:')
fprintf('\n    Speed                  = %10.3f km/s',v(end))
fprintf('\n    Flight path angle      = %10.3f deg ',gamma(end))
fprintf('\n    Altitude               = %10.3f km  ',h(end))
fprintf('\n    Downrange distance     = %10.3f km  ',x(end))
fprintf('\n    Drag loss              = %10.3f km/s',vD(end))
fprintf('\n    Gravity loss           = %10.3f km/s',vG(end))
fprintf('\n\n -----------------------------------\n')
 
 
figure('Name','Trajectory and Dynamic Pressure')
subplot(2,1,1)
plot(x,h)
title('(a) Altitude vs Downrange Distance')
axis equal
xlabel('Downrange Distance (km)')
ylabel('Altitude (km)')
axis([-inf, inf, -inf, inf])
grid
 
subplot(2,1,2)
plot(h, q*9.869e-6)
title('(b) Dynamic Pressure vs Altitude')
xlabel('Altitude (km)')
ylabel('Dynamic pressure (atm)')
axis([-inf, inf, -inf, inf])
xticks([0:10:120])
grid
 
end %output
 
end %gravity_turn