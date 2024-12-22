function f_and_g_ta_test()
%{
    The script defines input parameters for r0, v0, dt, and mu, and then
    calculates the f and g coefficients using the f_and_g_ta function.
    It also prints the results for inspection.

    Inputs:
    r0  - Initial position vector (km)
    v0  - Initial velocity vector (km/s)
    dt  - Change in true anomaly (degrees)
    mu  - Gravitational parameter (km^3/s^2)

    Outputs:
    f   - Lagrange f coefficient (dimensionless)
    g   - Lagrange g coefficient (s)

    User M-functions required: f_and_g_ta
%}
% -------------------------------------------

% Define input parameters:
mu = 398600;               % Gravitational parameter for Earth (km^3/s^2)
r0 = [7000, 0, 0];         % Initial position vector (km)
v0 = [0, 7.5, 0];          % Initial velocity vector (km/s)
dt = 45;                   % Change in true anomaly (degrees)

% Call the f_and_g_ta function:
[f, g] = f_and_g_ta(r0, v0, dt, mu);

% Display the results:
fprintf('Results:\n');
fprintf('Lagrange f coefficient: %f\n', f);
fprintf('Lagrange g coefficient: %f\n', g);

end