function [f, g] = f_and_g(x, t, ro, a)
%{
    This function calculates the Lagrange f and g coefficients.

    mu - the gravitational parameter (km^3/s^2)
    a - reciprocal of the semimajor axis (1/km)
    ro - the radial position at time to (km)
    t - the time elapsed since ro (s)
    x - the universal anomaly after time t (km^0.5)
    f - the Lagrange f coefficient (dimensionless)
    g - the Lagrange g coefficient (s)

    User M-functions required: stumpC, stumpS
%}
% --------------------------------------------------

global mu

z = a*x^2;

%...Equation 3.69a:
f = 1 - x^2/ro*stumpC(z);

%...Equation 3.69b:
g = t - 1/sqrt(mu)*x^3*stumpS(z);

end
