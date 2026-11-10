% ALGORITHM 5.3: CALCULATION OF LOCAL SIDEREAL TIME

function lst = LST(y, m, d, ut, EL)
%{
    This function calculates the local sidereal time.
    lst - local sidereal time (degrees)

    y   - year
    m   - month
    d   - day
    ut  - Universal Time (hours)
    EL  - east longitude (degrees)
    j0  - Julian day number at 0 hr UT
    j   - number of centuries since J2000
    g0  - Greenwich sidereal time (degrees) at 0 hr UT
    gst - Greenwich sidereal time (degrees) at the specified UT

    User M-function required: J0
    User subfunction required: zeroTo360
%}
% ----------------------------------------------

%...Equation 5.48;
j0 = J0(y, m, d);

%...Equation 5.49:
j = (j0 - 2451545)/36525;

%...Equation 5.50:
g0 = 100.4606184 + 36000.77004*j + 0.000387933*j^2 - 2.583e-8*j^3;

%...Reduce g0 so it lies in the range 0 - 360 degrees
g0 = zeroTo360(g0);

%...Equation 5.51:
gst = g0 + 360.98564724*ut/24;

%...Equation 5.52:
lst = gst + EL;

%...Reduce lst to the range 0 - 360 degrees:
lst = mod(lst,360);

return

    % -----------------------
    function y = zeroTo360(x)
    % -----------------------
    %{
        This subfunction reduces an angle to the range 0 - 360 degrees.
    
        x - The angle (degrees) to be reduced
        y - The reduced value
    %}
    % -------------------------
    y = mod(x,360);
    end %zeroTo360

end %LST