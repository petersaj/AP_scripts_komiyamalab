function sse = AP_fitExp(params,t,actual_trace)
%function fitted_curve = AP_fitExp(params,t)
A = params(1);
t0 = params(2);
tau_on = params(3);
tau1 = params(4);

fitted_curve = zeros(length(t));
fitted_curve = A*(1-exp(-(t-t0)/tau_on)).*exp(-(t-t0)./tau1);
fitted_curve(t < t0) = 0;
error_vector = fitted_curve - actual_trace;
sse = sum(error_vector.^2);