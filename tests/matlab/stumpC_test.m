function stumpC_test()
%{
    This function tests the stumpC function for various cases, 
    including positive, negative, and zero values of z. 

    The function also verifies specific cases with known outcomes.

    User M-functions required: stumpC
%}
% ----------------------------------------------

    % Define a range of test values for z
    z_values = [-10, -1, -0.1, 0, 0.1, 1, 10];

    % Initialize an array to store results
    c_values = zeros(size(z_values));

    % Loop through test values and compute C(z)
    for i = 1:length(z_values)
        z = z_values(i);
        c_values(i) = stumpC(z);
    end

    % Display the results
    fprintf('z\t\tC(z)\n');
    fprintf('-------------------\n');
    for i = 1:length(z_values)
        fprintf('%.2f\t\t%.6f\n', z_values(i), c_values(i));
    end

    % Verify specific cases with known outcomes
    assert(abs(stumpC(0) - 0.5) < 1e-6, 'Test failed for z = 0');
    assert(abs(stumpC(1) - (1 - cos(sqrt(1)))/1) < 1e-6, 'Test failed for z > 0');
    assert(abs(stumpC(-1) - (cosh(sqrt(1)) - 1)/1) < 1e-6, 'Test failed for z < 0');
end