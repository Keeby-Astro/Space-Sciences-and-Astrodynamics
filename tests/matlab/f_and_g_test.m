function f_and_g_test()
%{
    This function tests the f_and_g function by providing sample inputs
    and verifying the outputs. 

    User M-functions required: f_and_g, stumpC, stumpS
%}

% Define the global variable mu (gravitational parameter)
global mu;
mu = 398600; % km^3/s^2, standard for Earth

%...Inputs
ro = 7000;   % km, radial position at time to
t = 500;     % s, time elapsed
x = 1.5;     % universal anomaly
a = 1/10000; % 1/km, reciprocal of the semimajor axis

% Call the f_and_g function
[f, g] = f_and_g(x, t, ro, a);

%...Results
fprintf('Results:\n');
fprintf('Lagrange f: %.6f\n', f);
fprintf('Lagrange g: %.6f\n', g);

end
