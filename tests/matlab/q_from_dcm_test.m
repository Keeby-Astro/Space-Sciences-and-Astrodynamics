function q_from_dcm_test()
%{
    Verifies the calculation of the quaternion from the direction
    cosine matrix using known inputs.

    Q - direction cosine matrix
    q - quaternion (where q(4) is the scalar part)

    User M-functions required: q_from_dcm
%}
% ----------------------------------------------

% Define test direction cosine matrices
Q1 = eye(3); % Identity matrix
Q2 = [1,  0,  0;
      0, -1,  0;
      0,  0, -1]; % 180-degree rotation about x-axis
Q3 = [-1,  0,  0;
       0,  1,  0;
       0,  0, -1]; % 180-degree rotation about y-axis
Q4 = [-1,  0,  0;
       0, -1,  0;
       0,  0,  1]; % 180-degree rotation about z-axis

% Calculate the quaternions using the function
q1 = q_from_dcm(Q1);
q2 = q_from_dcm(Q2);
q3 = q_from_dcm(Q3);
q4 = q_from_dcm(Q4);

% Display the results
fprintf('Results for direction cosine matrix Q1:\n');
fprintf('Quaternion: [%f, %f, %f, %f]\n', q1);

fprintf('Results for direction cosine matrix Q2:\n');
fprintf('Quaternion: [%f, %f, %f, %f]\n', q2);

fprintf('Results for direction cosine matrix Q3:\n');
fprintf('Quaternion: [%f, %f, %f, %f]\n', q3);

fprintf('Results for direction cosine matrix Q4:\n');
fprintf('Quaternion: [%f, %f, %f, %f]\n', q4);

end