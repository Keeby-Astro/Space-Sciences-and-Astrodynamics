function coe = coe_from_sv(R,V,mu)
%{
    This function computes the classical orbital elements (coe)
    from the state vector (R,V) using Algorithm 4.1.

    mu      - gravitational parameter (km^3/s^2)
    R       - position vector in the geocentric equatorial frame (km)
    V       - velocity vector in the geocentric equatorial frame (km)
    r, v    - the magnitudes of R and V
    vr      - radial velocity component (km/s)
    H       - the angular momentum vector (km^2/s)
    h       - the magnitude of H (km^2/s)
    incl    - inclination of the orbit (rad)
    N       - the node line vector (km^2/s)
    n       - the magnitude of N
    cp      - cross product of N and R
    RA      - right ascension of the ascending node (rad)
    E       - eccentricity vector
    e       - eccentricity (magnitude of E)
    tol     - a small number below which the eccentricity is considered
              to be zero
    w       - argument of perigee (rad)
    TA      - true anomaly (rad)
    a       - semimajor axis (km)
    pi      - 3.1415926...
    coe     - vector of orbital elements [h e RA incl w TA a]

    User M-functions required: none
%}
% ---------------------------------------------

tol = 1.e-10;
R = R(:)';
V = V(:)';
if numel(R) ~= 3 || numel(V) ~= 3 || any(~isfinite([R V mu])) || ~isscalar(mu) || mu <= 0
    error('Finite three-vectors and a positive mu are required.')
end
r = norm(R);
v = norm(V);
H = cross(R,V);
h = norm(H);
if r == 0 || h == 0
    error('Classical elements require nonzero position and angular momentum.')
end
vr = dot(R,V)/r;

%...Equation 4.7:
incl = atan2(hypot(H(1),H(2)),H(3));

%...Equation 4.8:
N = cross([0 0 1],H);
n = norm(N);

%...Equation 4.10:
E = ((v^2 - mu/r)*R - r*vr*V)/mu;
e = norm(E);

%...Equation 4.9:
RA = 0;
if n > tol*h
    RA = mod(atan2(N(2),N(1)),2*pi);
end

%...Undefined angles are set to zero; the remaining angle retains the orientation.

%...Equation 4.12 (incorporating the case e = 0):
if e > tol
    if n > tol*h
        w = mod(atan2(dot(cross(N,E),H)/h,dot(N,E)),2*pi);
    else
        w = mod(atan2(sign(H(3))*E(2),E(1)),2*pi);
    end

    %...Equation 4.13a:
    TA = mod(atan2(dot(cross(E,R),H)/h,dot(E,R)),2*pi);
else
    w = 0;
    if n > tol*h
        TA = mod(atan2(dot(cross(N,R),H)/h,dot(N,R)),2*pi);
    else
        TA = mod(atan2(sign(H(3))*R(2),R(1)),2*pi);
    end
end

%...Equation 4.62 (a < 0 for a hyperbola):
den = 1 - e^2;
if abs(den) <= 32*eps
    a = Inf;
else
    a = h^2/mu/den;
end
coe = [h e RA incl w TA a];
end