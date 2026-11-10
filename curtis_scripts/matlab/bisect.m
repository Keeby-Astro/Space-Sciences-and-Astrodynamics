% ALGORITHM 2.4: FIND THE ROOT OF A FUNCTION USING THE
% BISECTION METHOD

function root = bisect(fun, xl, xu)
%{
    This function evaluates a root of a function using
    the bisection method.

    tol  - error to within which the root is computed
    n    - number of iterations
    xl   - low end of the interval containing the root
    xu   - upper end of the interval containing the root
    i    - loop index
    xm   - mid-point of the interval from xl to xu
    fun  - name of the function whose root is being found
    fxl  - value of fun at xl
    fxm  - value of fun at xm
    root - the computed root

    User M-functions required: none
%}
% ----------------------------------------------

tol = 1.e-6;
if ~all(isfinite([xl xu])) || xl > xu
    error('Bounds must be finite and ordered.')
end
fxl = fun(xl);
fxu = fun(xu);
if ~all(isfinite([fxl fxu]))
    error('Function values must be finite.')
end
if fxl == 0
    root = xl;
    return
end
if fxu == 0
    root = xu;
    return
end
if sign(fxl) == sign(fxu)
    error('Bounds must bracket a root.')
end
while xu - xl > 2*tol
    xm = xl + (xu - xl)/2;
    fxm = fun(xm);
    if ~isfinite(fxm)
        error('Function values must be finite.')
    end
    if fxm == 0
        root = xm;
        return
    end
    if sign(fxl) == sign(fxm)
        xl = xm;
        fxl = fxm;
    else
        xu = xm;
    end
end
root = xl + (xu - xl)/2;
end