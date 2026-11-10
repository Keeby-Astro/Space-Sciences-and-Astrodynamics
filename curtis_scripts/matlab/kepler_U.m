% ALGORITHM 3.3: SOLUTION OF THE UNIVERSAL KEPLER'S EQUATION
% USING NEWTON'S METHOD

function x = kepler_U(dt, ro, vro, a)
%{
    This function uses Newton's method to solve the universal
    Kepler equation for the universal anomaly.

    mu   - gravitational parameter (km^3/s^2)
    x    - the universal anomaly (km^0.5)
    dt   - time since x = 0 (s)
    ro   - radial position (km) when x = 0
    vro  - radial velocity (km/s) when x = 0
    a    - reciprocal of the semimajor axis (1/km)
    z    - auxiliary variable (z = a*x^2)
    C    - value of Stumpff function C(z)
    S    - value of Stumpff function S(z)
    n    - number of iterations for convergence
    lo, hi - lower and upper bounds on the signed anomaly magnitude

    User M-functions required: stumpC, stumpS
%}
% ---------------------------------------------
global mu

%...Check the supplied orbital parameters:
if ~all(isfinite([dt ro vro a mu])) || ro <= 0 || isempty(mu) || mu <= 0
    error('Universal propagation requires finite inputs and positive ro and mu.')
end
if dt == 0
    x = 0;
    return
end

%...Bracket the universal anomaly in the direction of propagation:
direction = sign(dt);
target = sqrt(mu)*abs(dt);
lo = 0;
if a < 0
    scale = sqrt(-1/a);
    hi = max(1,scale*log1p(target/(ro*scale)));
else
    hi = max(1,target*max(a,1/ro));
end
for n = 1:100
    [value,~] = equation(hi);
    if ~isfinite(value)
        error('Universal anomaly bracket overflowed.')
    end
    if value >= 0
        break
    end
    hi = 2*hi;
end
if value < 0
    error('Could not bracket the universal anomaly.')
end

%...Iterate on Equation 3.65 with a safeguarded Newton step
%...until the residual or bracket width satisfies the error tolerance:
u = (lo+hi)/2;
for n = 1:200
    [value,derivative] = equation(u);
    if ~all(isfinite([value derivative]))
        error('Universal anomaly iteration is nonfinite.')
    end
    if abs(value) <= 1.e-12*max(target,1)
        x = direction*u;
        return
    end
    if value > 0
        hi = u;
    else
        lo = u;
    end
    trial = u-value/derivative;
    if derivative <= 0 || ~isfinite(trial) || trial <= lo || trial >= hi
        trial = (lo+hi)/2;
    end
    if hi-lo <= 1.e-12*max(1,abs(u))
        x = direction*trial;
        return
    end
    u = trial;
end
error('Universal Kepler equation failed to converge.')

    %...Universal Kepler equation and its derivative:
    function [F,dFdx] = equation(u)
        xx = direction*u;
        C     = stumpC(a*xx^2);
        S     = stumpS(a*xx^2);
        F     = direction*(ro*vro/sqrt(mu)*xx^2*C + (1-a*ro)*xx^3*S + ro*xx - sqrt(mu)*dt);
        dFdx  = ro*vro/sqrt(mu)*xx*(1-a*xx^2*S) + (1-a*ro)*xx^2*C + ro;
    end
end