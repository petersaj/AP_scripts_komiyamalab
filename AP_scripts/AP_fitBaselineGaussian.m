function error = AP_fitBaselineGaussian(estimated_baseline,estimated_std,trace)
% find baseline closest to estimated noise level

% threshold the trace to get below what's potential baseline, normalize
threshold_trace = trace(trace < estimated_baseline)-estimated_baseline;
% flip the thresholded trace to simulate even distribution
mirrored_trace = [threshold_trace -threshold_trace];

error = abs(estimated_std - std(mirrored_trace));

% % new approach: minimize std AND mean about estimated baseline
% 
% % threshold the trace to get below what's potential baseline, normalize
% below_threshold_trace = trace(trace < estimated_baseline)-estimated_baseline;
% std_error = abs(estimated_std - std(below_threshold_trace));
% 
% baseline_ub = estimated_baseline + 2*estimated_std;
% above_threshold_trace = trace(trace > estimated_baseline & ...
%     trace < baseline_ub);
% mean_error = abs(mean(below_threshold_trace) - mean(above_threshold_trace));
% 
% error = std_error + mean_error;
% 
% % new approach: just minimize squared error with a straight line
% % threshold the trace to get below what's potential baseline, normalize
% threshold_trace = trace(trace < estimated_baseline);
% % flip the thresholded trace to simulate even distribution
% mirrored_trace = [threshold_trace -threshold_trace];
% 
% baseline_ub = estimated_baseline + 2*estimated_std;
% baseline_lb = estimated_baseline - 2*estimated_std;
% baseline_trace = trace(trace > baseline_lb &...
%     trace < baseline_ub);
% 
% straight_line = estimated_baseline*ones(1,length(threshold_trace));
% 
% error = sum((threshold_trace-straight_line).^2);