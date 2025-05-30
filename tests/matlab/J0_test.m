% J0_test
%{
    This program computes J0 and the Julian day number using the data.
    
    year   - range: 1901 - 2099
    month  - range: 1 - 12
    day    - range: 1 - 31
    hour   - range: 0 - 23 (Universal Time)
    minute - rage: 0 - 60
    second - range: 0 - 60
    ut     - universal time (hr)
    j0     - Julian day number at 0 hr UT
    jd     - Julian day number at specified UT

    User M-function required: J0
%}
% ----------------------------------------------

clear all; clc

%...Data declaration:
year  = 2004;
month = 5;
day   = 12;

hour   = 14;
minute = 45;
second = 30;
%...

ut = hour + minute/60 + second/3600;

%...Equation 5.46:
j0 = J0(year, month, day);

%...Equation 5.47:
jd = j0 + ut/24;

%...Echo the input data and output the results to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.4: Julian day calculation\n')
fprintf('\n Input data:\n');
fprintf('\n   Year            = %g', year)
fprintf('\n   Month           = %g', month)
fprintf('\n   Day             = %g', day)
fprintf('\n   Hour            = %g', hour)
fprintf('\n   Minute          = %g', minute)
fprintf('\n   Second          = %g\n', second)

fprintf('\n Julian day number = %11.3f', jd);
fprintf('\n-----------------------------------------------------\n')
