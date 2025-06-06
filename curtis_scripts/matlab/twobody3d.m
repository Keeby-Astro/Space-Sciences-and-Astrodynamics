% ALGORITHM 2.1: NUMERICAL SOLUTION OF THE TWO-BODY PROBLEM
% RELATIVE TO AN INERTIAL FRAME

function twobody3d
%{
    This function solves the inertial two-body problem in three dimensions
    numerically using the RKF4(5) method.

    G               - universal gravitational constant (km^3/kg/s^2)
    m1,m2           - the masses of the two bodies (kg)
    m               - the total mass (kg)
    t0              - initial time (s)
    tf              - final time (s)
    R1_0,V1_0       - 3 by 1 column vectors containing the components of tbe
                      initial position (km) and velocity (km/s) of m1
    R2_0,V2_0       - 3 by 1 column vectors containing the components of the
                      initial position (km) and velocity (km/s) of m2
    y0              - 12 by 1 column vector containing the initial values
                      of the state vectors of the two bodies:
                      [R1_0; R2_0; V1_0; V2_0]
    t               - column vector of the times at which the solution is found
    X1,Y1,Z1        - column vectors containing the X,Y and Z coordinates (km)
                      of m1 at the times in t
    X2,Y2,Z2        - column vectors containing the X,Y and Z coordinates (km)
                      of m2 at the times in t
    VX1, VY1, VZ1   - column vectors containing the X,Y and Z components
                      of the velocity (km/s) of m1 at the times in t
    VX2, VY2, VZ2   - column vectors containing the X,Y and Z components
                      of the velocity (km/s) of m2 at the times in t
    y               - a matrix whose 12 columns are, respectively,
                      X1,Y1,Z1; X2,Y2,Z2; VX1,VY1,VZ1; VX2,VY2,VZ2
    XG,YG,ZG        - column vectors containing the X,Y and Z coordinates (km)
                      the center of mass at the times in t

    User M-function required: rkf45
    User subfunctions required: rates, output
%}
% -------------------------------------------------------------------------
clc; clear all; close all
G = 6.67259e-20;

%...Input data:
m1  = 1.e26;
m2  = 1.e26;
t0  = 0;
tf  = 480;

R1_0 = [   0;   0;   0];
R2_0 = [3000;   0;   0];

V1_0 = [  10;  20;  30];
V2_0 = [   0;  40;   0];
%...End input data

y0 = [R1_0; R2_0; V1_0; V2_0];

%...Integrate the equations of motion:
[t,y] = rkf45(@rates, [t0 tf], y0);

%...Output the results:
output

return

% ------------------------
function dydt = rates(t,y)
% ------------------------
%{
    This function calculates the accelerations in Equations 2.19.

    t       - time
    y       - column vector containing the position and velocity vectors
              of the system at time t
    R1, R2  - position vectors of m1 & m2
    V1, V2  - velocity vectors of m1 & m2
    r       - magnitude of the relative position vector
    A1, A2  - acceleration vectors of m1 & m2
    dydt    - column vector containing the velocity and acceleration
              vectors of the system at time t
%}
% ------------------------
R1 = [y(1); y(2); y(3)];
R2 = [y(4); y(5); y(6)];

V1 = [y(7); y(8); y(9)];
V2 = [y(10); y(11); y(12)];

r = norm(R2 - R1);

A1 = G*m2*(R2 - R1)/r^3;
A2 = G*m1*(R1 - R2)/r^3;

dydt = [V1; V2; A1; A2];
end %rates

    % -------------
    function output
    % -------------
    %{
        This function calculates the trajectory of the center of mass and
        plots
        (a) the motion of m1, m2 and G relative to the inertial frame
        (b) the motion of m2 and G relative to m1
        (c) the motion of m1 and m2 relative to G
    
        User subfunction required: common_axis_settings
    %}
    % -------------
    
    %...Extract the particle trajectories:
    X1 = y(:,1); Y1 = y(:,2); Z1 = y(:,3);
    X2 = y(:,4); Y2 = y(:,5); Z2 = y(:,6);
    
    %...Locate the center of mass at each time step:
    XG = []; YG = []; ZG = [];
    for i = 1:length(t)
        XG = [XG; (m1*X1(i) + m2*X2(i))/(m1 + m2)];
        YG = [YG; (m1*Y1(i) + m2*Y2(i))/(m1 + m2)];
        ZG = [ZG; (m1*Z1(i) + m2*Z2(i))/(m1 + m2)];
    end
    
    %...Plot the trajectories:
    figure (1)
    title('Figure 2.3: Motion relative to the inertial frame')
    hold on
    plot3(X1, Y1, Z1, '-r')
    plot3(X2, Y2, Z2, '-g')
    plot3(XG, YG, ZG, '-b')
    common_axis_settings
    
    figure (2)
    title('Figure 2.4a: Motion of m2 and G relative to m1')
    hold on
    plot3(X2 - X1, Y2 - Y1, Z2 - Z1, '-g')
    plot3(XG - X1, YG - Y1, ZG - Z1, '-b')
    common_axis_settings
    
    figure (3)
    title('Figure 2.4b: Motion of m1 and m2 relative to G')
    hold on
    plot3(X1 - XG, Y1 - YG, Z1 - ZG, '-r')
    plot3(X2 - XG, Y2 - YG, Z2 - ZG, '-g')
    common_axis_settings
    
        % ---------------------------
        function common_axis_settings
        % ---------------------------
        %{
            This function establishes axis properties common to the several plots.
        %}
        % ---------------------------
        text(0, 0, 0, 'o')
        axis('equal')
        view([2,4,1.2])
        grid on
        axis equal
        xlabel('X (km)')
        ylabel('Y (km)')
        zlabel('Z (km)')
        end %common_axis_settings
    
    end %output

end %twobody3d