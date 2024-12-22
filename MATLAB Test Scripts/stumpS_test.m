function stumpS_test()
%{
    This function tests the stumpS function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User M-functions required: stumpS
%}
% ----------------------------------------------

    % Define a range of test values for z
    z_values = [-10, -1, -0.1, 0, 0.1, 1, 10];

    % Initialize an array to store results
    s_values = zeros(size(z_values));

    % Loop through test values and compute S(z)
    for i = 1:length(z_values)
        z = z_values(i);
        s_values(i) = stumpS(z);
    end

    % Display the results
    fprintf('z\t\tS(z)\n');
    fprintf('-------------------\n');
    for i = 1:length(z_values)
        fprintf('%.2f\t\t%.6f\n', z_values(i), s_values(i));
    end

    % Verify specific cases with known outcomes
    assert(abs(stumpS(0) - 1/6) < 1e-6, 'Test failed for z = 0');
    assert(abs(stumpS(1) - (sqrt(1) - sin(sqrt(1)))/(sqrt(1))^3) < 1e-6, 'Test failed for z > 0');
    assert(abs(stumpS(-1) - (sinh(sqrt(1)) - sqrt(1))/(sqrt(1))^3) < 1e-6, 'Test failed for z < 0');
end