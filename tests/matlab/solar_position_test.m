function solar_position_test()
%{
    Ouputs the geocentric equatorial position vector
    of the sun, given the julian date.
    
    User M-functions required: None
%}
% -------------------------------------------------------------------------

    % Julian date for 2023-12-21 00:00:00 UTC (Winter Solstice)
    jd_test = 2460271.5; % Julian date

    % Call the solar_position function
    [lamda, eps, r_S] = solar_position(jd_test);

    % Display the results
    fprintf('Julian Date = %.2f\n', jd_test);
    fprintf('Apparent Ecliptic Longitude (lamda): %.6f degrees\n', lamda);
    fprintf('Obliquity of the Ecliptic (eps): %.6f degrees\n', eps);
    fprintf('Geocentric Position Vector (r_S):\n');
    fprintf('  x = %.6f km\n', r_S(1));
    fprintf('  y = %.6f km\n', r_S(2));
    fprintf('  z = %.6f km\n', r_S(3));

end
