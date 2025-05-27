%month_planet_names_test
%{
    Outputs the name of the month and the planet corresponding, respectively,
     to the numbers "month_id" and "planet_id".
    
    month_id  - the month number (1 - 12)
    planet_id - the planet number (1 - 9)

    User M-functions required: month_planet_names
%}
% -------------------------------------------------------------------

disp('Month and Planet Names Test');
disp('---------------------------');

% Test month and planet IDs
for month_id = 1:12
    for planet_id = 1:9
        [month, planet] = month_planet_names(month_id, planet_id);
        
        % Output
        fprintf('Month ID: %2d, Month: %s | Planet ID: %2d, Planet: %s\n', ...
                month_id, strtrim(month), planet_id, strtrim(planet));
    end
end