% gibbs_test
%{
    This program uses Algorithm 5.1 (Gibbs method) and Algorithm 4.2
    to obtain the orbital elements from the data.

    deg        - factor for converting between degrees and radians
    pi         - 3.1415926...
    mu         - gravitational parameter (km^3/s^2)
    r1, r2, r3 - three coplanar geocentric position vectors (km)
    ierr       - 0 if r1, r2, r3 are found to be coplanar
                 1 otherwise
    v2         - the velocity corresponding to r2 (km/s)
    coe        - orbital elements [h e RA incl w TA a]
                 where h = angular momentum (km^2/s)
                 e       = eccentricity
                 RA      = right ascension of the ascending node (rad)
                 incl    = orbit inclination (rad)
                 w       = argument of perigee (rad)
                 TA      = true anomaly (rad)
                 a       = semimajor axis (km)
    T          - period of elliptic orbit (s)

    User M-functions required: gibbs, coe_from_sv
%}
% ----------------------------------------------

clear all; clc
deg = pi/180;
global mu

%...Data declaration for Example 5.1:
mu = 398600;
r1 = [-294.32 4265.1 5986.7];
r2 = [-1365.5 3637.6 6346.8];
r3 = [-2940.3 2473.7 6555.8];
%...

%...Echo the input data to the command window:
fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n')
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu)
fprintf('\n r1 (km) = [%g %g %g]', r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', r2(1), r2(2), r2(3))
fprintf('\n r3 (km) = [%g %g %g]', r3(1), r3(2), r3(3))
fprintf('\n\n');

%...Algorithm 5.1:
[v2, ierr] = gibbs(r1, r2, r3);

%...If the vectors r1, r2, r3, are not coplanar, abort:
if ierr == 1
    fprintf('\n These vectors are not coplanar.\n\n')
    return
end

%...Algorithm 4.2:
coe = coe_from_sv(r2,v2,mu);

h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);
a    = coe(7);

%...Output the results to the command window:
fprintf(' Solution:')
fprintf('\n');
fprintf('\n v2 (km/s) = [%g %g %g]', v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:');
fprintf('\n Angular momentum (km^2/s)  = %g', h)
fprintf('\n Eccentricity               = %g', e)
fprintf('\n Inclination (deg)          = %g', incl/deg)
fprintf('\n RA of ascending node (deg) = %g', RA/deg)
fprintf('\n Argument of perigee (deg)  = %g', w/deg)
fprintf('\n True anomaly (deg)         = %g', TA/deg)
fprintf('\n Semimajor axis (km)        = %g', a)
%...If the orbit is an ellipse, output the period:
if e < 1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n Period (s)             = %g', T)
end
fprintf('\n-----------------------------------------------------\n')