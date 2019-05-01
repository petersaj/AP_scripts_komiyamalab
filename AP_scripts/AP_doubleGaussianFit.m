function R_error = AP_doubleGaussianFit(params, R_measure)
% Hiroshi: double gaussian fit

a0 = params(1);
a1 = params(2);
a2 = params(3);
theta0 = params(4);
sigma = params(5);

theta = 1:1:360;

R_fit = a0+a1*exp(-(theta-theta0).^2/(2*sigma^2))+a2*exp(-(theta-(theta0+180)).^2/(2*sigma^2));
R_error = sum((R_measure-R_fit).^2)/sum((R_measure-mean(R_measure)).^2);