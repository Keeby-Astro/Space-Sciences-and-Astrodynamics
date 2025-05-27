function los_test()
%{
    Uses the ECI position vectors of the satellite (r_sat) and the sun (r_sun)
    and evaluates whether the earth is in the line of sight between the two.

    User M-functions required: None
%}
% -----------------------------------------------------------------------------

% Scenario 1: Satellite in sunlight
r_sat1 = [7000, 0, 0]; % Satellite position vector (ECI, km)
r_sun1 = [1.496e8, 0, 0]; % Sun position vector (ECI, km)
light_switch1 = los(r_sat1, r_sun1);

if light_switch1 == 1
    fprintf('Satellite is in sunlight\n');
else
    fprintf('Satellite is not in sunlight\n');
end

% Scenario 2: Satellite in Earth\'s shadow
r_sat2 = [7000, 0, 0];     % Satellite position vector (ECI, km)
r_sun2 = [-1.496e8, 0, 0]; % Sun position vector (ECI, km)
light_switch2 = los(r_sat2, r_sun2);

if light_switch2 == 0
    fprintf('Satellite is in Earth''s shadow\n');
else
    fprintf('Satellite is not in Earth''s shadow\n');
end

% Scenario 3: Satellite on Earth\'s surface
RE = 6378; % Earth\'s radius (km)
r_sat3 = [RE, 0, 0];      % Satellite on Earth\'s surface
r_sun3 = [1.496e8, 0, 0]; % Sun position vector (ECI, km)
light_switch3 = los(r_sat3, r_sun3);

if light_switch3 == 1
    fprintf('Satellite on Earth''s surface is in sunlight\n');
else
    fprintf('Satellite on Earth''s surface is not in sunlight\n');
end

end
