% ALGORITHM 3.2: SOLUTION OF KEPLER'S EQUATION FOR THE
% HYPERBOLA USING NEWTON'S METHOD

function F = kepler_H(e, M)
%{
    This function uses Newton's method to solve Kepler's equation
    for the hyperbola e*sinh(F) - F = M for the hyperbolic
    eccentric anomaly, given the eccentricity and the hyperbolic
    mean anomaly.

    F - hyperbolic eccentric anomaly (radians)
    e - eccentricity, passed from the calling program
    M - hyperbolic mean anomaly (radians), passed from the
    calling program

    User M-functions required: none
%}
% --------------------------------------------

%...Check the supplied orbital parameters:
if ~all(isfinite([e M])) || e <= 1
    error('Hyperbolic Kepler equation requires e > 1 and finite M.')
end
if M == 0
    F = 0;
    return
end

%...Starting value for F and an upper bound on the root:
m = abs(M);
lo = 0;
hi = asinh(m/e)+1;
while e*sinh(hi)-hi < m
    hi = 2*hi;
end
F = asinh(m/e);

%...Iterate on Equation 3.45, using bisection when a Newton step leaves
%...the bracket. Limit the iteration count and test the error tolerance:
for n = 1:100
    residual = e*sinh(F)-F-m;
    if residual == 0
        F = sign(M)*F;
        return
    end
    if residual > 0
        hi = F;
    else
        lo = F;
    end
    trial = F-residual/(e*cosh(F)-1);
    if trial <= lo || trial >= hi
        trial = (lo+hi)/2;
    end
    if abs(trial-F) <= 1.e-13*max(1,abs(F))
        F = sign(M)*trial;
        return
    end
    F = trial;
end
error('Hyperbolic Kepler equation failed to converge.')
end