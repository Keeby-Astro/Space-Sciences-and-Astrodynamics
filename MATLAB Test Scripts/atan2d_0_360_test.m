% atan2d_0_360_test
%{
    This script tests the atan2d_0_360 function for various inputs.

    t - angle in degrees

    User M-function required: atan2d_0_360
%}
% ----------------------------------------------

%...Inputs
inputs = [ 0, 0; 1, 0;   0, 1; -1, 0;
          0, -1; 1, 1; -1, -1; 1, -1;
          -1, 1; sqrt(3), 1; -sqrt(3), -1;
          sqrt(3), -1; -sqrt(3), 1];

for i = 1:size(inputs, 1)
    y = inputs(i, 1);
    x = inputs(i, 2);
    result = atan2d_0_360(y, x);

    fprintf('Test case %d: atan2d_0_360(%f, %f)\n', i, y, x);
end
