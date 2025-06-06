% dcm_to_euler_test 
%{
    This script tests the dcm_to_euler function using various input
    direction cosine matrices and verifies the output Euler angles.

    Q     - direction cosine matrix

    User M-function required: dcm_to_euler
%}
% ---------------------------------------------------

% Input Matrices
input_matrices = {
    eye(3), ...                             % Identity matrix (no rotation)
    [ 0 -1  0; 1  0  0; 0  0  1], ...      % 90-degree rotation about Z-axis
    [ 0  0 -1; 0  1  0; 1  0  0], ...      % 90-degree rotation about X-axis
    [ 0  1  0; -1  0  0; 0  0  1], ...     % -90-degree rotation about Z-axis
    [ cosd(45) -sind(45) 0; ...            % 45-degree rotation about Z-axis
      sind(45)  cosd(45) 0; ...
      0         0        1]
};

for i = 1:length(input_matrices)
    Q = input_matrices{i};

    % Compute Euler Angles
    [alpha, beta, gamma] = dcm_to_euler(Q);

    % Results
    fprintf('Input DCM:\n');
    disp(Q);
    fprintf('Output Euler Angles: alpha = %.2f, beta = %.2f, gamma = %.2f (degrees)\n', alpha, beta, gamma);
    fprintf('\n');
end