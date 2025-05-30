% ALGORITHM 7.1: FIND THE POSITION, VELOCITY, AND ACCELERATION
% OF B RELATIVE TO A'S LVLH FRAME

function [r_rel_x, v_rel_x, a_rel_x] = rva_relative(rA,vA,rB,vB)
%{
    This function uses the state vectors of spacecraft A and B
    to find the position, velocity and acceleration of B relative
    to A in the LVLH frame attached to A (see Figure 7.1).

    rA,vA       - state vector of A (km, km/s)
    rB,vB       - state vector of B (km, km/s)
    mu          - gravitational parameter (km^3/s^2)
    hA          - angular momentum vector of A (km^2/s)
    i, j, k     - unit vectors along the x, y and z axes of A's
                  LVLH frame
    QXx         - DCM of the LVLH frame relative to the geocentric
                  equatorial frame (GEF)
    Omega       - angular velocity of the LVLH frame (rad/s)
    Omega_dot   - angular acceleration of the LVLH frame (rad/s^2)
    aA, aB      - absolute accelerations of A and B (km/s^2)
    r_rel       - position of B relative to A in GEF (km)
    v_rel       - velocity of B relative to A in GEF (km/s)
    a_rel       - acceleration of B relative to A in GEF (km/s^2)
    r_rel_x     - position of B relative to A in the LVLH frame
    v_rel_x     - velocity of B relative to A in the LVLH frame
    a_rel_x     - acceleration of B relative to A in the LVLH frame

    User M-functions required: None
%}
% ------------------------------------------------------------

global mu

%...Calculate the vector hA:
hA = cross(rA, vA);

%...Calculate the unit vectors i, j and k:
i = rA/norm(rA);
k = hA/norm(hA);
j = cross(k,i);

%...Calculate the transformation matrix Qxx:
QXx = [i; j; k];

%...Calculate Omega and Omega_dot:
Omega     = hA/norm(rA)^2; % Equation 7.5
Omega_dot = -2*dot(rA,vA)/norm(rA)^2*Omega;% Equation 7.6

%...Calculate the accelerations aA and aB:
aA = -mu*rA/norm(rA)^3;
aB = -mu*rB/norm(rB)^3;

%...Calculate r_rel:
r_rel = rB - rA;

%...Calculate v_rel:
v_rel = vB - vA - cross(Omega,r_rel);

%...Calculate a_rel:
a_rel = aB - aA - cross(Omega_dot,r_rel)...
       - cross(Omega,cross(Omega,r_rel))...
       - 2*cross(Omega,v_rel);

%...Calculate r_rel_x, v_rel_x and a_rel_x:
r_rel_x = QXx*r_rel';
v_rel_x = QXx*v_rel';
a_rel_x = QXx*a_rel';

end