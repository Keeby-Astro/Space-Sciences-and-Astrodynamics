% quat_rotate_test
%{
    Evaluates known quaternions and vectors to compute
    that the rotation properly.
    
    q - Quaternion to rotate by 1 x 4 matrix
    v - Vector to rotate by 1 x 3 matrix
    r - Rotated vector by 1 x 3 matrix

    User M-functions required: quat_rotate
%}
q1 = [sqrt(2)/2, 0, 0, sqrt(2)/2]; % Quaternion for 90-degree Z-rotation
q2 = [0, 0, 1, 0];                 % Quaternion for 180-degree Y-rotation
q3 = [1, 0, 0, 0];                 % Identity quaternion (no rotation)
q4 = [sqrt(2)/2, sqrt(2)/2, 0, 0]; % Quaternion for 90-degree X-rotation

% Vector to rotate
v1 = [1, 0, 0]; 
v2 = [1, 0, 0];
v3 = [1, 2, 3]; 
v4 = [0, 1, 0]; 

% Rotated vector
r1 = quat_rotate(q1, v1); 
r2 = quat_rotate(q2, v2); 
r3 = quat_rotate(q3, v3);
r4 = quat_rotate(q4, v4);

% Results
fprintf('90-degree rotation around Z-axis\n');
fprintf('Input Vector: [%f, %f, %f]\n', v1);
fprintf('Rotated Vector: [%f, %f, %f]\n\n', r1);

fprintf('180-degree rotation around Y-axis\n');
fprintf('Input Vector: [%f, %f, %f]\n', v2);
fprintf('Rotated Vector: [%f, %f, %f]\n\n', r2);

fprintf('No rotation\n');
fprintf('Input Vector: [%f, %f, %f]\n', v3);
fprintf('Rotated Vector: [%f, %f, %f]\n\n', r3);

fprintf('90-degree rotation around X-axis\n');
fprintf('Input Vector: [%f, %f, %f]\n', v4);
fprintf('Rotated Vector: [%f, %f, %f]\n\n', r4);