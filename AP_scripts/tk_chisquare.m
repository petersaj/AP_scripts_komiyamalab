% Chi-square test for whether two proportions are statistically significantly different.
% input [n1, N1, n2, N2]:
% n1/N1 vs n2/N2 fractions
function p = tk_chisquare(n1, N1, n2, N2)
if isrow(n1)
    n1 = n1';
end
if isrow(N1)
    N1 = N1';
end
if isrow(n2)
    n2 = n2';
end
if isrow(N2)
    N2 = N2';
end
% Pooled estimate of proportion
p0 = (n1+n2) ./ (N1+N2);
% Expected counts under H0 (null hypothesis)
n10 = N1 .* p0;
n20 = N2 .* p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected,2);
p = 1 - chi2cdf(chi2stat,1);
 