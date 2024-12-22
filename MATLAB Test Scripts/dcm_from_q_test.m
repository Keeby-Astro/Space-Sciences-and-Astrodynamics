function dcm_from_q_test()
%{
    Verifies the direction cosine matrix from known quaternions.

    q - quaternion (where q(4) is the scalar part)
    Q - direction cosine matrix

    User M-functions required: dcm_from_q
%}
% ----------------------------------------------
q1 = [1, 0, 0, 0];  % Identity quaternion
q2 = [0, 1, 0, 0];  % 180-degree rotation about x-axis
q3 = [0, 0, 1, 0];  % 180-degree rotation about y-axis
q4 = [0, 0, 0, 1];  % 180-degree rotation about z-axis

% Calculate the direction cosine matrices using the function
Q1 = dcm_from_q(q1);
Q2 = dcm_from_q(q2);
Q3 = dcm_from_q(q3);
Q4 = dcm_from_q(q4);

% Display the results
fprintf('Results for quaternion q1: [%d, %d, %d, %d]\n', q1);
fprintf('Direction Cosine Matrix:\n');
fprintf('%f %f %f\n', Q1.');

fprintf('Results for quaternion q2: [%d, %d, %d, %d]\n', q2);
fprintf('Direction Cosine Matrix:\n');
fprintf('%f %f %f\n', Q2.');

fprintf('Results for quaternion q3: [%d, %d, %d, %d]\n', q3);
fprintf('Direction Cosine Matrix:\n');
fprintf('%f %f %f\n', Q3.');

fprintf('Results for quaternion q4: [%d, %d, %d, %d]\n', q4);
fprintf('Direction Cosine Matrix:\n');
fprintf('%f %f %f\n', Q4.');

end