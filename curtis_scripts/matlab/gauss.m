% ALGORITHMS 5.5 AND 5.6: GAUSS' METHOD OF PRELIMINARY ORBIT
% DETERMINATION WITH ITERATIVE IMPROVEMENT

function [r, v, r_old, v_old] = ...
gauss(Rho1, Rho2, Rho3, R1, R2, R3, t1, t2, t3)
%{
    This function uses the Gauss method with iterative improvement
    (Algorithms 5.5 and 5.6) to calculate the state vector of an
    orbiting body from angles-only observations at three
    closely spaced times.

    mu                  - the gravitational parameter (km^3/s^2)
    t1, t2, t3          - the times of the observations (s)
    tau, tau1, tau3     - time intervals between observations (s)
    R1, R2, R3          - the observation site position vectors
                          at t1, t2, t3 (km)
    Rho1, Rho2, Rho3    - the direction cosine vectors of the
                          satellite at t1, t2, t3
    p1, p2, p3          - cross products among the three direction
                          cosine vectors
    Do                  - scalar triple product of Rho1, Rho2 and Rho3
    D                   - Matrix of the nine scalar triple products
                          of R1, R2 and R3 with p1, p2 and p3
    E                   - dot product of R2 and Rho2
    A, B                - constants in the expression relating slant range
                          to geocentric radius
    a,b,c               - coefficients of the 8th order polynomial
                          in the estimated geocentric radius x
    x                   - positive root of the 8th order polynomial
    rho1, rho2, rho3    - the slant ranges at t1, t2, t3
    r1, r2, r3          - the position vectors at t1, t2, t3 (km)
    r_old, v_old        - the estimated state vector at the end of
                          Algorithm 5.5 (km, km/s)
    rho1_old,
    rho2_old, and
    rho3_old            - the values of the slant ranges at t1, t2, t3
                          at the beginning of iterative improvement
                          (Algorithm 5.6) (km)
    diff1, diff2,
    and diff3           - the magnitudes of the differences between the
                          old and new slant ranges at the end of
                          each iteration
    tol                 - the error tolerance determining
                          convergence
    n                   - number of passes through the
                          iterative improvement loop
    nmax                - limit on the number of iterations
    ro, vo              - magnitude of the position and
                          velocity vectors (km, km/s)
    vro                 - radial velocity component (km)
    a                   - reciprocal of the semimajor axis (1/km)
    v2                  - computed velocity at time t2 (km/s)
    r, v                - the state vector at the end of Algorithm 5.6
                          (km, km/s)

    User m-functions required: kepler_U, f_and_g
    User subfunctions required: posroot
%}
% --------------------------------------------------

global mu

%...Equations 5.98:
tau1 = t1 - t2;
tau3 = t3 - t2;

%...Equation 5.101:
tau = tau3 - tau1;

%...Independent cross products among the direction cosine vectors:
p1 = cross(Rho2,Rho3);
p2 = cross(Rho1,Rho3);
p3 = cross(Rho1,Rho2);

%...Equation 5.108:
Do = dot(Rho1,p1);

%...Equations 5.109b, 5.110b and 5.111b:
D = [[dot(R1,p1) dot(R1,p2) dot(R1,p3)]
     [dot(R2,p1) dot(R2,p2) dot(R2,p3)]
     [dot(R3,p1) dot(R3,p2) dot(R3,p3)]];

%...Equation 5.115b:
E = dot(R2,Rho2);

%...Equations 5.112b and 5.112c:
A = 1/Do*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);
B = 1/6/Do*(D(1,2)*(tau3^2 - tau^2)*tau3/tau ...
            + D(3,2)*(tau^2 - tau1^2)*tau1/tau);

%...Equations 5.117:
a = -(A^2 + 2*A*E + norm(R2)^2);
b = -2*mu*B*(A + E);
c = -(mu*B)^2;

%...Calculate the roots of Equation 5.116 using MATLAB's
% polynomial 'roots' solver:
Roots = roots([1 0 a 0 0 b 0 0 c]);

%...Find the positive real root:
x = posroot(Roots);

%...Equations 5.99a and 5.99b:
f1 = 1 - 1/2*mu*tau1^2/x^3;
f3 = 1 - 1/2*mu*tau3^2/x^3;

%...Equations 5.100a and 5.100b:
g1 = tau1 - 1/6*mu*(tau1/x)^3;
g3 = tau3 - 1/6*mu*(tau3/x)^3;

%...Equation 5.112a:
rho2 = A + mu*B/x^3;

%...Equation 5.113:
rho1 = 1/Do*((6*(D(3,1)*tau1/tau3 + D(2,1)*tau/tau3)*x^3 ...
        + mu*D(3,1)*(tau^2 - tau1^2)*tau1/tau3) ...
        /(6*x^3 + mu*(tau^2 - tau3^2)) - D(1,1));

%...Equation 5.114:
rho3 = 1/Do*((6*(D(1,3)*tau3/tau1 - D(2,3)*tau/tau1)*x^3 ...
        + mu*D(1,3)*(tau^2 - tau3^2)*tau3/tau1) ...
        /(6*x^3 + mu*(tau^2 - tau1^2)) - D(3,3));

%...Equations 5.86:
r1 = R1 + rho1*Rho1;
r2 = R2 + rho2*Rho2;
r3 = R3 + rho3*Rho3;

%...Equation 5.118:
v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

%...Save the initial estimates of r2 and v2:
r_old = r2;
v_old = v2;

%...End of Algorithm 5.5

%...Use Algorithm 5.6 to improve the accuracy of the initial estimates.
%...Initialize the iterative improvement loop and set error tolerance:
rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;
diff1    = 1;    diff2    = 1;    diff3    = 1;
n    = 0;
nmax = 1000;
tol  = 1.e-8;

%...Iterative improvement loop:
while ((diff1 > tol) & (diff2 > tol) & (diff3 > tol)) & (n < nmax)
    n = n+1;

    %...Compute quantities required by universal kepler's equation:
    ro = norm(r2);
    vo = norm(v2);
    vro = dot(v2,r2)/ro;
    a   = 2/ro - vo^2/mu;

    %...Solve universal Kepler's equation at times tau1 and tau3 for
    % universal anomalies x1 and x3:
    x1 = kepler_U(tau1, ro, vro, a);
    x3 = kepler_U(tau3, ro, vro, a);

    %...Calculate the Lagrange f and g coefficients at times tau1
    % and tau3:
    [ff1, gg1] = f_and_g(x1, tau1, ro, a);
    [ff3, gg3] = f_and_g(x3, tau3, ro, a);

    %...Update the f and g functions at times tau1 and tau3 by
    % averaging old and new:
    f1 = (f1 + ff1)/2;
    f3 = (f3 + ff3)/2;
    g1 = (g1 + gg1)/2;
    g3 = (g3 + gg3)/2;

    %...Equations 5.96 and 5.97:
    c1 = g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);

    %...Equations 5.109a, 5.110a and 5.111a:
    rho1 = 1/Do*( -D(1,1) + 1/c1*D(2,1) - c3/c1*D(3,1));
    rho2 = 1/Do*( -c1*D(1,2) +D(2,2) - c3*D(3,2));
    rho3 = 1/Do*(-c1/c3*D(1,3) + 1/c3*D(2,3) - D(3,3));

    %...Equations 5.86:
    r1 = R1 + rho1*Rho1;
    r2 = R2 + rho2*Rho2;
    r3 = R3 + rho3*Rho3;

    %...Equation 5.118:
    v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

    %...Calculate differences upon which to base convergence:
    diff1 = abs(rho1 - rho1_old);
    diff2 = abs(rho2 - rho2_old);
    diff3 = abs(rho3 - rho3_old);

    %...Update the slant ranges:
    rho1_old = rho1; rho2_old = rho2; rho3_old = rho3;
end
%...End iterative improvement loop

fprintf('\n( **Number of Gauss improvement iterations = %g)\n\n',n)

if n >= nmax
    fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax);
end

%...Return the state vector for the central observation:
r = r2;
v = v2;

return

    % -------------------------
    function x = posroot(Roots)
    % -------------------------
    %{
        This subfunction extracts the positive real roots from
        those obtained in the call to MATLAB's 'roots' function.
        If there is more than one positive root, the user is
        prompted to select the one to use.
    
        x        - the determined or selected positive root
        Roots    - the vector of roots of a polynomial
        posroots - vector of positive roots
    
        User M-functions required: none
    %}
    
    %...Construct the vector of positive real roots:
    posroots = Roots(find(Roots>0 & ~imag(Roots)));
    npositive = length(posroots);
    
    %...Exit if no positive roots exist:
    if npositive == 0
        fprintf('\n\n ** There are no positive roots. \n\n')
        return
    end
    
    %...If there is more than one positive root, output the
    % roots to the command window and prompt the user to
    % select which one to use:
    if npositive == 1
        x = posroots;
    else
        fprintf('\n\n ** There are two or more positive roots.\n')
        for i = 1:npositive
            fprintf('\n root #%g = %g',i,posroots(i))
        end
        fprintf('\n\n Make a choice:\n')
        nchoice = 0;
        while nchoice < 1 | nchoice > npositive
            nchoice = input(' Use root #? ');
        end
        x = posroots(nchoice);
        fprintf('\n We will use %g .\n', x)
    end

    end %posroot

end %gauss