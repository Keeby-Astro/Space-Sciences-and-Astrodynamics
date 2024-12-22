%simpsons_lunar_ephemeris_test
%{
    Outputs the state vector of the moon at a given time
    relative to the earth's geocentric equatorial frame using a curve fit
    to JPL's DE200 (1982) ephemeris model.
    
    jd      - julian date (days)
    pos     - position vector (km)
    vel     - velocity vector (km/s)
    a       - matrix of amplitudes (km)
    b       - matrix of frequencies (rad/century)
    c       - matrix of phase angles (rad)
    t       - time in centuries since J2000
    tfac    - no. of seconds in a Julian century (36525 days)

    User M-functions required: simpsons_lunar_ephemeris
%}

% Define a sample Julian date (e.g., January 1, 2024, at 12:00 TT)
jd_sample = 2460271.0; % Julian date

% Call the simpsons_lunar_ephemeris function
[pos, vel] = simpsons_lunar_ephemeris(jd_sample);

% Display the Julian date
fprintf('Julian date:\n');
fprintf('JD: %.6f\n', jd_sample);

% Display the position vector
fprintf('Position vector (km):\n');
fprintf('X: %.6f km\n', pos(1));
fprintf('Y: %.6f km\n', pos(2));
fprintf('Z: %.6f km\n', pos(3));

% Display the velocity vector
fprintf('\nVelocity vector (km/s):\n');
fprintf('VX: %.6f km/s\n', vel(1));
fprintf('VY: %.6f km/s\n', vel(2));
fprintf('VZ: %.6f km/s\n', vel(3));