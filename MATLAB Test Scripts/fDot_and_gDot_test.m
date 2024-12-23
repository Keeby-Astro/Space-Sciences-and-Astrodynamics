% fDot_and_gDot_test
%{
    Outputs the time derivatives of the Lagrange f and g coefficients.

    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)

    User M-functions required: fDot_and_gDot
%}
% --------------------------------------------------
global mu;
mu = 398600;

% Inputs
x = 1.5;    % Universal anomaly (km^0.5)
r = 7000;   % Radial position after time t (km)
ro = 6878;  % Radial position at time to (km)
a = 1/8000; % Reciprocal of semimajor axis (1/km)

[fdot, gdot] = fDot_and_gDot(x, r, ro, a);

% Display results
fprintf('Results:\n');
fprintf('f_dot = %.6e 1/s\n', fdot);
fprintf('g_dot = %.6f (dimensionless)\n', gdot);