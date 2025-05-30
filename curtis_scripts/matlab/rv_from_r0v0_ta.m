% ALGORITHM 2.3: CALCULATE THE STATE VECTOR FROM THE INITIAL
% STATE VECTOR AND THE CHANGE IN TRUE ANOMALY

function [r,v] = rv_from_r0v0_ta(r0, v0, dt, mu)
%{
    This function computes the state vector (r,v) from the
    initial state vector (r0,v0) and the change in true anomaly.

    mu - gravitational parameter (km^3/s^2)
    r0 - initial position vector (km)
    v0 - initial velocity vector (km/s)
    dt - change in true anomaly (degrees)
    r  - final position vector (km)
    v  - final velocity vector (km/s)

    User M-functions required: f_and_g_ta, fDot_and_gDot_ta
%}
% -------------------------------------------------------------------------

%global mu

%...Compute the f and g functions and their derivatives:
[f, g]       =       f_and_g_ta(r0, v0, dt, mu);
[fdot, gdot] = fDot_and_gDot_ta(r0, v0, dt, mu);

%...Compute the final position and velocity vectors:
r =    f*r0 +g*v0;
v = fdot*r0 + gdot*v0;

end