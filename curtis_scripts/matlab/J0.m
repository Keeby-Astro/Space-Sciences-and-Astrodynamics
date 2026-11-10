function j0 = J0(year, month, day)
%{
    This function computes the Julian day number at 0 UT for any year
    between 1900 and 2100 using Equation 5.48.

    j0      - Julian day at 0 hr UT (Universal Time)
    year    - range: 1901 - 2099
    month   - range: 1 - 12
    day     - range: 1 - 31

    User M-functions required: none
%}
% ----------------------------------

if any(~isfinite([year month day])) || any(fix([year month day]) ~= [year month day]) ...
        || year < 1901 || year > 2099 || month < 1 || month > 12
    error('Use an integer calendar date in 1901..2099.')
end
month_days = [31 28+(mod(year,4)==0) 31 30 31 30 31 31 30 31 30 31];
if day < 1 || day > month_days(month)
    error('Invalid calendar day.')
end
j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
    + fix(275*month/9) + day + 1721013.5;

end