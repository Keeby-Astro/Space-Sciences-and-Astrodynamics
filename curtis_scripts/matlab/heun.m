% ALGORITHM 1.2: NUMERICAL INTEGRATION BY HEUN'S PREDICTOR-CORRECTOR METHOD

function [tout, yout] = heun(ode_function, tspan, y0, h)
%{
    This function uses the predictor-corrector method to integrate a system
    of first-order differential equations dy/dt = f(t,y).

    y               - column vector of solutions
    f               - column vector of the derivatives dy/dt
    ode_function    - handle for the user M-function in which the derivatives
                      f are computed
    t               - time
    t0              - initial time
    tf              - final time
    tspan           - the vector [t0 tf] giving the time interval for the
                      solution
    h               - time step
    y0              - column vector of initial values of the vector y
    tout            - column vector of the times at which y was evaluated
    yout            - a matrix, each row of which contains the components of y
                      evaluated at the correponding time in tout
    feval           - a built-in MATLAB function which executes 'ode_function'
                      at the arguments t and y
    tol             - Maximum allowable relative error for determining
                      convergence of the corrector
    itermax         - maximum allowable number of iterations for corrector
                      convergence
    iter            - iteration number in the corrector convergence loop
    t1              - time at the beginning of a time step
    y1              - value of y at the beginning of a time step
    f1              - derivative of y at the beginning of a time step
    f2              - derivative of y at the end of a time step
    favg            - average of f1 and f2
    y2p             - predicted value of y at the end of a time step
    y2              - corrected value of y at the end of a time step
    err             - maximum relative error (for all components) between y2p
                      and y2 for given iteration
    eps             - unit roundoff error (the smallest number for which
                      1 + eps > 1). Used to avoid a zero denominator.

    User M-function required: ode_function
%}
% ----------------------------------------------------------------------

tol = 1.e-6;
itermax = 100;

t0 = tspan(1);
tf = tspan(2);
if ~all(isfinite(tspan)) || tf < t0 || isempty(y0) || ~isvector(y0) || any(~isfinite(y0(:)))
    error('Finite forward times and a finite initial vector are required.')
end
if ~isscalar(h) || ~isfinite(h) || h <= 0
    error('Step must be finite and positive.')
end

t = t0;
y = y0(:);
tout = t;
yout = y';

while t < tf
    h = min(h, tf-t);
    if t+h == t
        error('Integration step cannot advance time.')
    end
    t1 = t;
    y1 = y;
    f1 = feval(ode_function, t1, y1);
    f1 = f1(:);
    if numel(f1) ~= numel(y) || any(~isfinite(f1))
        error('Derivative must match the state and be finite.')
    end
    y2 = y1 + f1*h;
    t2 = t1 + h;
    err = tol + 1;
    iter = 0;
    while err > tol && iter <= itermax
        y2p = y2;
        f2 = feval(ode_function, t2, y2p);
        f2 = f2(:);
        if numel(f2) ~= numel(y) || any(~isfinite(f2))
            error('Derivative must match the state and be finite.')
        end
        favg = (f1 + f2)/2;
        y2 = y1 + favg*h;
        err = max(abs((y2 - y2p)./(y2 + eps)));
        iter = iter + 1;
    end

    if iter > itermax
        error('Heun corrector failed to converge.')
    end

    t = t + h;
    y = y2;
    tout = [tout;t]; % adds t to the bottom of the column vector tout
    yout = [yout;y']; % adds y' to the bottom of the matrix yout
end