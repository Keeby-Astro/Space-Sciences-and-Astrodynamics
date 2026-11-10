function s = stumpS(z)
%{
    This function evaluates the Stumpff function S(z) according
    to Equation 3.52.

    z - input argument
    s - value of S(z)

    User M-functions required: none
%}
% ----------------------------------------------
if abs(z) < 1.e-3
    s = 1/6 - z/120 + z^2/5040 - z^3/362880 + z^4/39916800;
elseif z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
    s = 1/6;
end