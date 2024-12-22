%lunar_position_test
%{
    Outputs the geocentric equatorial position vector of the moon
    given the Julian day.
    
    User M-functions required: None
%}
% -------------------------------------------------------------------------

jd_1 = 2451545.0; % Julian Date for J2000 epoch
result_1 = lunar_position(jd_1);
jd_2 = 2459200.5; % Julian Date for September 22, 2020
result_2 = lunar_position(jd_2);
jd_3 = 2440587.5; % Julian Date for January 1, 1970
result_3 = lunar_position(jd_3);

% Display results
fprintf('Julian Date = %.1f\n', jd_1);
fprintf('Geocentric position vector (km):\n');
fprintf('  X = %.3f km\n', result_1(1));
fprintf('  Y = %.3f km\n', result_1(2));
fprintf('  Z = %.3f km\n\n', result_1(3));

fprintf('Julian Date = %.1f\n', jd_2);
fprintf('Geocentric position vector (km):\n');
fprintf('  X = %.3f km\n', result_2(1));
fprintf('  Y = %.3f km\n', result_2(2));
fprintf('  Z = %.3f km\n\n', result_2(3));
fprintf('Julian Date = %.1f\n', jd_3);

fprintf('Geocentric position vector (km):\n');
fprintf('  X = %.3f km\n', result_3(1));
fprintf('  Y = %.3f km\n', result_3(2));
fprintf('  Z = %.3f km\n\n', result_3(3));