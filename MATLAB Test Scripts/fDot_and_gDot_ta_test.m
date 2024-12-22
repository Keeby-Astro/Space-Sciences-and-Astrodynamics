function fDot_and_gDot_ta_test()
%{
    This script tests the function fDot_and_gDot_ta by passing example
    inputs and verifying the outputs. It computes the time derivatives of
    the Lagrange f and g coefficients for a given position and velocity
    vector, and a change in true anomaly.

    r0      - position vector at time t0 (km)
    v0      - velocity vector at time t0 (km/s)
    dt      - change in true anomaly (degrees)
    mu      - gravitational parameter (km^3/s^2)
    fdot    - time derivative of the Lagrange f coefficient (1/s)
    gdot    - time derivative of the Lagrange g coefficient (dimensionless)

    User M-functions required: fDot_and_gDot_ta
%}
% --------------------------------------------------------

% Inputs
r0 = [7000, 0, 0];        % Initial position vector (km)
v0 = [0, 7.546, 0];       % Initial velocity vector (km/s)
dt = 30;                  % Change in true anomaly (degrees)
mu = 398600;              % Gravitational parameter (km^3/s^2)

% Call the function
[fdot, gdot] = fDot_and_gDot_ta(r0, v0, dt, mu);

% Display the results
disp('Test Results for fDot_and_gDot_ta:');
disp(['fdot: ', num2str(fdot), ' (1/s)']);
disp(['gdot: ', num2str(gdot), ' (dimensionless)']);

end