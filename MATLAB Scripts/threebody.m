function threebody
    %{
        This program presents the graphical solution of the motion of three
        bodies in the plane for data provided in the input definitions below.
    
        MATLAB's ode45 Runge-Kutta solver is used.
    
        G                           - gravitational constant (km^3/kg/s^2)
        t0, tf                      - initial and final times (s)
        m1, m2, m3                  - masses of the three bodies (kg)
        m                           - total mass (kg)
        X1,Y1; X2,Y2; X3,Y3         - coordinates of the three masses (km)
        VX1,VY1; VX2,VY2; VX3,VY3   - velocity components of the three
                                      masses (km/s)
        XG, YG                      - coordinates of the center of mass (km)
        y0                          - column vector of the initial conditions
        t                           - column vector of times at which the solution
                                      was computed
        y                           - matrix, the columns of which contain the
                                      position and velocity components evaluated at
                                      the times t(:):
                                        y(:,1) , y(:, 2) = X1(:), Y1(:)
                                        y(:,3) , y(:, 4) = X2(:), Y2(:)
                                        y(:,5) , y(:, 6) = X3(:), Y3(:)
                                        y(:,7) , y(:, 8) = VX1(:), VY1(:)
                                        y(:,9) , y(:,10) = VX2(:), VY2(:)
                                        y(:,11), y(:,12) = VX3(:), VY3(:)
        
        User M-functions required: none
        User subfunctions required: rates, plotit
    %}
    % -------------------------------------------------------------------------
    clear all
    close all
    clc
    
    G = 6.67259e-20;
    
    %...Input data:
    m1 = 1.e29; m2 = 1.e29; m3 = 1.e29;
    
    t0 = 0; tf = 67000;
    
    X1 = 0; Y1 = 0;
    X2 = 300000; Y2 = 0;
    X3 = 2*X2; Y3 = 0;
    
    VX1 = 0; VY1 = 0;
    VX2 = 250; VY2 = 250;
    VX3 = 0; VY3 = 0;
    %...End input data
    
    m = m1 + m2 + m3;
    y0 = [X1 Y1 X2 Y2 X3 Y3 VX1 VY1 VX2 VY2 VX3 VY3]';
    
    %...Pass the initial conditions and time interval to ode45, which
    % calculates the position and velocity of each particle at discrete
    % times t, returning the solution in the column vector y. ode45 uses
    % the subfunction 'rates' below to evaluate the accelerations at each
    % integration time step.
    [t,y] = ode45(@rates, [t0 tf], y0);
    
    X1 = y(:,1); Y1 = y(:,2);
    X2 = y(:,3); Y2 = y(:,4);
    X3 = y(:,5); Y3 = y(:,6);
    
    %...Locate the center of mass at each time step:
    XG = []; YG = [];
    for i = 1:length(t)
        XG = [XG; (m1*X1(i) + m2*X2(i) + m3*X3(i))/m];
        YG = [YG; (m1*Y1(i) + m2*Y2(i) + m3*Y3(i))/m];
    end
    
    %...Coordinates of each particle relative to the center of mass:
    X1G = X1 - XG; Y1G = Y1 - YG;
    X2G = X2 - XG; Y2G = Y2 - YG;
    X3G = X3 - XG; Y3G = Y3 - YG;
    
    plotit
    
    return
    % -------------------------------------------------------------------------
    function dydt = rates(t,y)
    %{
        This function evaluates the acceleration of each member of a planar
        3-body system at time t from their positions and velocities
        at that time.

        t                           - time (s)
        y                           - column vector containing the position and
                                      velocity components of the three masses
                                      at time t
        R12                         - cube of the distance between m1 and m2 (km^3)
        R13                         - cube of the distance between m1 and m3 (km^3)
        R23                         - cube of the distance between m2 and m3 (km^3)
        AX1,AY1; AX2,AY2; AX3,AY3   - acceleration components of each mass (km/s^2)
        dydt                        - column vector containing the velocity and
                                      acceleration components of the three
                                      masses at time t
    %}
    % -------------------------------------------------------------------------
    X1 = y( 1);
    Y1 = y( 2);
    X2 = y( 3);
    Y2 = y( 4);
    X3 = y( 5);
    Y3 = y( 6);
    VX1 = y( 7);
    VY1 = y( 8);
    VX2 = y( 9);
    VY2 = y(10);
    VX3 = y(11);
    VY3 = y(12);

    %...Equations C.8:
    R12 = norm([X2 - X1, Y2 - Y1])^3;
    R13 = norm([X3 - X1, Y3 - Y1])^3;
    R23 = norm([X3 - X2, Y3 - Y2])^3;

    %...Equations C.7:
    AX1 = G*m2*(X2 - X1)/R12 + G*m3*(X3 - X1)/R13;
    AY1 = G*m2*(Y2 - Y1)/R12 + G*m3*(Y3 - Y1)/R13;
    AX2 = G*m1*(X1 - X2)/R12 + G*m3*(X3 - X2)/R23;
    AY2 = G*m1*(Y1 - Y2)/R12 + G*m3*(Y3 - Y2)/R23;
    AX3 = G*m1*(X1 - X3)/R13 + G*m2*(X2 - X3)/R23;
    AY3 = G*m1*(Y1 - Y3)/R13 + G*m2*(Y2 - Y3)/R23;

    dydt = [VX1 VY1 VX2 VY2 VX3 VY3 AX1 AY1 AX2 AY2 AX3 AY3]';
    end %rates
    % -------------------------------------------------------------------------
    function plotit
    %...Plot the motions relative to the inertial frame (Figure 2.4):
    figure(1)
    title('Figure 2.4: Motion relative to the inertial frame', ...
           'Fontweight', 'bold', 'FontSize', 12)
    hold on
    plot(XG, YG, '--k', 'LineWidth', 0.25)
    plot(X1, Y1, 'r', 'LineWidth', 0.5)
    plot(X2, Y2, 'g', 'LineWidth', 0.75)
    plot(X3, Y3, 'b', 'LineWidth', 1.00)
    xlabel('X(km)'); ylabel('Y(km)')
    grid on
    axis('equal')

    %...Plot the motions relative to the center of mass (Figure 2.5):
    figure(2)
    title('Figure 2.5: Motion relative to the center of mass', ...
          'Fontweight', 'bold', 'FontSize', 12)
    hold on
    plot(X1G, Y1G, 'r', 'LineWidth', 0.5)
    plot(X2G, Y2G, '--g', 'LineWidth', 0.75)
    plot(X3G, Y3G, 'b', 'LineWidth', 1.00)
    xlabel('X(km)'); ylabel('Y(km)')
    grid on
    axis('equal')
    end %plotit
    
end %threebody