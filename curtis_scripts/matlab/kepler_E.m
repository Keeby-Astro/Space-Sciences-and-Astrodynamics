% ALGORITHM 3.1: SOLUTION OF KEPLER'S EQUATION BY NEWTON'S METHOD

function E = kepler_E(e, M)
%{
    This function uses Newton's method to solve Kepler's
    equation E - e*sin(E) = M for the eccentric anomaly,
    given the eccentricity and the mean anomaly.

    E  - eccentric anomaly (radians)
    e  - eccentricity, passed from the calling program
    M  - mean anomaly (radians), passed from the calling program
    pi - 3.1415926...

    User M-functions required: none
%}
% ----------------------------------------------

%...Check the supplied orbital parameters:
if ~all(isfinite([e M])) || e < 0 || e >= 1
    error('Elliptic Kepler equation requires 0 <= e < 1 and finite M.')
end

%...Select a starting value for E and bracket the root:
m = mod(M+pi,2*pi)-pi;
lo = m-e;
hi = m+e;
E = m + sign(m)*0.85*e;

%...Iterate on Equation 3.17, using bisection when a Newton step leaves
%...the bracket. Limit the iteration count and test the error tolerance:
for n = 1:100
    F = E-e*sin(E)-m;
    if F == 0
        E = E+(M-m);
        return
    end
    if F > 0
        hi = E;
    else
        lo = E;
    end
    trial = E-F/(1-e*cos(E));
    if trial <= lo || trial >= hi
        trial = (lo+hi)/2;
    end
    if abs(trial-E) <= 1.e-13*max(1,abs(E))
        E = trial+(M-m);
        return
    end
    E = trial;
end
error('Elliptic Kepler equation failed to converge.')
end