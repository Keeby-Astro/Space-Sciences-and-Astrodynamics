function atmosphere_test
%{
    Tests the atmosphere function by plotting the density
    as a function of altitude from sea level through 1000 km.
%}

% Locate companion functions relative to this test file
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'curtis_scripts', 'matlab'));

%...Geometric altitudes (km):
z = linspace(0, 1000, 500); % Generate 500 points from 0 to 1000 km

%...Calculate densities:
density = arrayfun(@atmosphere, z);

%...Plot the results:
figure;
semilogy(z, density, 'LineWidth', 1.5);
grid on;
title('Atmospheric Density vs. Altitude');
xlabel('Altitude (km)');
ylabel('Density (kg/m^3)');
legend('Density', 'Location', 'northeast');

end