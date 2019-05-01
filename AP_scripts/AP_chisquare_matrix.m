% Chi-square test for whether two proportions are statistically significantly different.
%
% Modified to work on matricies:
% input [x,y]
% x and y matricies of observations(rows) by dimensions (columns)
% must have the same number of dimensions (i.e. equal cell number)

function p = AP_chisquare_matrix(x,y)
if size(x,2) ~= size(y,2)
    error('Number of dimensions not the same')
end

n1 = nansum(x,1);
n2 = nansum(y,1);

N1 = sum(~isnan(x),1);
N2 = sum(~isnan(y),1);

% Pooled estimate of proportion
p0 = (n1+n2) ./ (N1+N2);
% Expected counts under H0 (null hypothesis)
n10 = N1 .* p0;
n20 = N2 .* p0;
% Chi-square test, by hand
observed = [n1;N1-n1;n2;N2-n2];
expected = [n10;N1-n10;n20;N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected,1);
p = 1 - chi2cdf(chi2stat,1);
 