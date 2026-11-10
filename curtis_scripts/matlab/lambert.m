function [V1, V2] = lambert(R1, R2, t, string)
%{
    This function solves the zero-revolution Lambert problem.

    mu      - gravitational parameter (km^3/s^2)
    R1, R2  - initial and final position vectors (km)
    r1, r2  - magnitudes of R1 and R2
    t       - the time of flight from R1 to R2 (a constant) (s)
    V1, V2  - initial and final velocity vectors (km/s)
    c12     - cross product of R1 into R2
    theta   - angle between R1 and R2
    string  - 'pro' if the orbit is prograde
              'retro' if the orbit is retrograde
    A       - a constant given by Equation 5.35
    z       - alpha*x^2, where alpha is the reciprocal of the
              semimajor axis and x is the universal anomaly
    y(z)    - a function of z given by Equation 5.38
    F(z,t)  - a function of the variable z and constant t,
            - given by Equation 5.40
    f, g    - Lagrange coefficients
    gdot    - time derivative of g
              C(z), S(z) - Stumpff functions
    dum     - a dummy variable

    User M-functions required: stumpC, stumpS
%}
% ----------------------------------------------

global mu
if nargin < 4
    string = 'pro';
end
if isempty(mu) || ~isscalar(mu) || ~isfinite(mu) || mu <= 0 || ~isfinite(t) || t <= 0
    error('Set global mu to a positive value and supply a positive flight time.')
end
R1 = R1(:)';
R2 = R2(:)';
if numel(R1) ~= 3 || numel(R2) ~= 3 || any(~isfinite([R1 R2]))
    error('Lambert requires finite position three-vectors.')
end
if ~strcmp(string,'pro') && ~strcmp(string,'retro')
    error('Orbit type must be pro or retro.')
end

%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1,R2);
if r1 == 0 || r2 == 0 || norm(c12) <= 1.e-12*r1*r2
    error('Collinear or zero vectors require a separate Lambert formulation.')
end

%...Determine whether the orbit is prograde or retrograde:
theta = acos(max(-1,min(1,dot(R1,R2)/(r1*r2))));
if (strcmp(string,'pro') && c12(3) <= 0) || (strcmp(string,'retro') && c12(3) >= 0)
    theta = 2*pi-theta;
end

%...Equation 5.35:
A = sin(theta)*sqrt(r1*r2/(1-cos(theta)));

%...Bracket the zero-revolution branch, z < 4*pi^2.
lo = -1;
hi = 4*pi^2*(1-1.e-7);
for n = 1:100
    value = F(lo);
    if ~isfinite(value)
        error('Lambert lower bracket overflowed.')
    end
    if value <= 0
        break
    end
    lo = 2*lo;
end
if value > 0
    error('Could not bracket Lambert solution.')
end
if ~isfinite(F(hi)) || F(hi) <= 0
    error('Flight time exceeds this solver zero-revolution bracket.')
end

%...Bisect the bracket until the time-of-flight residual satisfies the tolerance:
for n = 1:200
    z = (lo+hi)/2;
    value = F(z);
    if ~isfinite(value)
        error('Lambert iteration is nonfinite.')
    end
    if abs(value) <= 1.e-11*sqrt(mu)*t
        break
    end
    if value > 0
        hi = z;
    else
        lo = z;
    end
end
if n == 200 && abs(value) > 1.e-11*sqrt(mu)*t
    error('Lambert solver failed to converge.')
end
yz = y(z);
if yz <= 0
    error('Lambert solution is outside the physical domain.')
end

%...Equation 5.46a:
f = 1-yz/r1;

%...Equation 5.46b:
g = A*sqrt(yz/mu);

%...Equation 5.46d:
gdot = 1-yz/r2;

%...Equation 5.28:
V1 = (R2-f*R1)/g;

%...Equation 5.29:
V2 = (gdot*R2-R1)/g;

    %...Equation 5.38:
    function yz = y(z)
        yz = r1+r2+A*(z*stumpS(z)-1)/sqrt(stumpC(z));
    end

    %...Equation 5.40:
    function value = F(z)
        yz = y(z);
        if yz <= 0
            value = -sqrt(mu)*t;
        else
            value = (yz/stumpC(z))^1.5*stumpS(z)+A*sqrt(yz)-sqrt(mu)*t;
        end
    end
end %lambert