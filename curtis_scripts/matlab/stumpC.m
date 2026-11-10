function c = stumpC(z)
%{
    This function evaluates the Stumpff function C(z) according
    to Equation 3.53.

    z - input argument
    c - value of C(z)

    User M-functions required: none
%}
% ----------------------------------------------
if abs(z) < 1.e-3
    c = 1/2 - z/24 + z^2/720 - z^3/40320 + z^4/3628800;
elseif z > 0
    c = (1 - cos(sqrt(z)))/z;
elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
else
    c = 1/2;
end