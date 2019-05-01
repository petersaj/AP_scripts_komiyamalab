function sse = AP_deconvSpike(params,oopsi_trace,trace)

t0 = 1;
amplitude = params(1);
tau_up = params(2);
tau_down = params(3);
t = [1:100];

event_template = amplitude*exp(t/tau_up);
%event_template(t0:end) = event_template(t0)*exp(-(t(t0:end)-t0)/tau_down);
event_template(t0:end) = amplitude*exp(-(t(t0:end)-t0)/tau_down);
event_template = fliplr(event_template);

fitted_trace = conv(oopsi_trace,event_template,'same');
sse = sum((trace - fitted_trace).^2);