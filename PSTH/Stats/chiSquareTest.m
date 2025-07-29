function [chi_square_stat, p_value, df] = chiSquareTest(observed, alpha)
%CHISQUARETEST Perform a Chi-square test.
%
%   [CHI_SQUARE_STAT, P_VALUE, DF] = chiSquareTest(OBSERVED, ALPHA)) performs a 
%   Chi-square test on the contingency table OBSERVED. It computes the 
%   Chi-square statistic, p-value, and degrees of freedom. The null hypothesis 
%   assumes that the variables are independent (for 2D tables). It allows 
%   specification of the significance level ALPHA.
%
%   Inputs:
%       OBSERVED - Matrix of observed frequencies (contingency table).
%       ALPHA    - Scalar, significance level (optional).
%
%   Outputs:
%       CHI_SQUARE_STAT - Chi-square test statistic (float)
%       P_VALUE         - P-value of the test (scalar, in range [0,1])
%       DF              - Degrees of freedom (integer)
%
%   Example:
%       % Test independence in a 2x3 contingency table
%       observed = [20 10 15; 15 25 10];
%       [stat, pVal, df] = chiSquareTest(observed, 0.01);
%
%   Notes:
%       - Expected frequencies less than 5 in many cells may invalidate 
%         the Chi-square approximation.
%
%  Malek Aouadi, Laboratory of Systems Neuroscience
%  Institute of Experimental Medicine, Budapest, Hungary
%  2025

    % Calculate the row and column sums
    row_sums = sum(observed, 2);
    col_sums = sum(observed, 1);
    total_sum = sum(row_sums);

    % Calculate the expected frequencies
    expected = (row_sums * col_sums) / total_sum;

    % Calculate the Chi-square statistic
    chi_square_stat = sum(sum((observed - expected).^2 ./ expected));

    % Degrees of freedom
    [num_rows, num_cols] = size(observed);
    df = (num_rows - 1) * (num_cols - 1);

    % Calculate the p-value
    p_value = 1 - chi2cdf(chi_square_stat, df);

    % Display the results 
    fprintf('Chi-square statistic: %.3f\n', chi_square_stat);
    fprintf('Degrees of freedom: %d\n', df);
    fprintf('P-value: %.2e\n', p_value);

    % Interpret the result 
    if p_value <= alpha
        fprintf('Result: Reject the null hypothesis. There is a significant association.\n');
    else
        fprintf('Result: Fail to reject the null hypothesis. No significant association.\n');
    end
end
