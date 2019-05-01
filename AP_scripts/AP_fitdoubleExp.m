function sse = AP_fitdoubleExp(params,t,actual_trace,t0_est,tau_on_est)
% function fitted_curve = AP_fitdoubleExp(params,input)
% embed previously gotten values, since don't know how to pass through
%t0 = params(5);
%tau_on = params(6);
t0 = t0_est;
tau_on = tau_on_est;
A1 = params(1);
tau1 = params(2);
A2 = params(3);
tau2 = params(4);

fitted_curve = zeros(length(t));
fitted_curve =(1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1) + A2*exp(-(t-t0)./tau2));
fitted_curve(t < t0) = 0;
error_vector = fitted_curve(t > t0) - actual_trace(t > t0);
sse = sum(error_vector.^2);