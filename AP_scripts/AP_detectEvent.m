function [spike_frames spike_amplitudes spike_estimates spike_t0s] = AP_detectEvent(roi_trace_long)

% Subtract exponential decay from transients
% [spike_frames spike_amplitudes spike_estimates spike_t0s] = AP_detectEvent(roi_trace_long)


% try
%     matlabpool
% catch me
% end
    
num_rois = size(roi_trace_long,1);
residuals_all = cell(num_rois,1);
spike_frames = cell(num_rois,1);
spike_amplitudes = cell(num_rois,1);
curr_peak_all = cell(num_rois,1);
spike_t0s = cell(num_rois,1);
spike_estimates = cell(num_rois,1);

parfor roi_num = 1:num_rois;
%for roi_num = 2;
% optional smoothing, can't go back at this point
%roi_trace_long(roi_num,:,:) = smooth(roi_trace_long(roi_num,:,:),2);

% set parameters
curr_roi = roi_num; 
sse_criterion = 0.1; % high threshold for sse acceptance
area_criterion = 0.4; % fraction of area fit cannot deviate from
frame_search = 30; % frames to look past peak
frame_back = 20; % frames to look behind peak
baseline_search = 200; % surrounding frames to estimate local baseline
max_std = 100; % keep events over this std amount, even if bad fit
peak_back = 10; % frames within Ca jump has to occur
peakstd = [2 2 2 2 2]; % definitions for max trigger, in stds

% get total std
if all(roi_trace_long(roi_num) == 0);
    continue
end
try
    norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',2); % there are two hidden curves
catch me
    norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',1);
end
baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve(1)));
% initialize event and baseline vectors
%%%peak_frames = roi_peak_max{curr_roi}(:,1); % don't use this anymore
event_detect = zeros(1,size(roi_trace_long,2));
baseline_noEvent = roi_trace_long(curr_roi,:);
event_estimate = zeros(1,size(roi_trace_long,2));
event_estimate_bad = zeros(1,size(roi_trace_long,2));
event_localBaseline = zeros(1,size(roi_trace_long,2));
% residuals_all{roi_num} = [];
% spike_frames{roi_num} = [];
% spike_amplitudes{roi_num} = [];
% curr_peak_all{roi_num} = [];

curr_point = peak_back + 1;
% go through all points
while curr_point < (size(roi_trace_long,2) - frame_search+1) && ...
    curr_point > peak_back;
    
    % find peaks, only use first
    peak_frame = [];
    peak_amplitude = [];
    curr_peak = [];
    minpeak = trace_std*peakstd;
    [peak_frame, peak_amplitude] = AP_schmittTrigger(baseline_noEvent(curr_point-peak_back:end),minpeak,peak_back,frame_search,trace_std);
    curr_peak = (curr_point-1) + peak_frame - peak_back;
    curr_point = curr_peak + 1;

    if isempty(curr_peak)
        curr_point = size(roi_trace_long,2);
        continue
    end
    
    if  curr_peak+frame_search-1 < size(roi_trace_long,2) && ...
            curr_peak - frame_back > 0;
        curr_peak_all{roi_num} = [curr_peak_all{roi_num} curr_peak];
        % amplitude of the current peak
        curr_peak_amp = roi_trace_long(curr_roi,curr_peak);
        % frames to compare event with Ca transient
        event_frames = [curr_peak-frame_back:curr_peak+frame_search-1];
        % correct for local baseline
        % get 50 frames, or as close as possible, for baseline estimation
        local_frames = [curr_peak-baseline_search:curr_peak+baseline_search];
        local_frames = local_frames(local_frames > 0 ...
            & local_frames < size(roi_trace_long,2));
        % estimate local (100 frames) baseline to correct for baseline
        % drift
        try
            norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,local_frames)',2);
        catch me
            continue
        end
        baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
        baseline_local = norm_fit_obj.mu(baseline_curve);
        % create template for local-baseline event, fit with double exp
        event_trace = baseline_noEvent(event_frames)-baseline_local(1);  %roi_trace_long(curr_roi,event_frames) - baseline_local;
        % initial fit to get t0
        A_guess = max(event_trace);
        starting = [A_guess 200 30 600];
        t = (0:frame_search+frame_back-1)*(1.24*128);
        
        %%% find fit
        
        options = optimset('Algorithm','interior-point');
        lb = [0 -frame_back*1.24*128 0 0];
        ub = [2*A_guess frame_back*1.24*128 1000 5000];
        % estimate transform of orig -> curr
        display 'Finding best fit...'
        estimates = fmincon(@(x) AP_fitExp(x,t,event_trace),starting,[],[],[],[],lb,ub,[],options);
        A = estimates(1);
        t0 = estimates(2);
        t0_frame = event_frames(find((t-t0) == min(t-t0),1)); 
        tau_on = estimates(3);
        tau1 = estimates(4);
        % entire transient fit
        starting = [A_guess,300,A_guess,800];%,t0,tau_on];
        t = (0:frame_search+frame_back-1)*(1.24*128);
        % embed values just gotten to pass to function
        lb = [0,0,0,0];%,-frame_back*1.24*128,0];
        ub = [5,5000,5,5000];%,1000,500];
        estimates = fmincon(@(x) AP_fitdoubleExp(x,t,event_trace,t0,tau_on),starting,[],[],[],[],lb,ub,[],options);
        %t0 = estimates(5);
        %tau_on = estimates(6);
        A1 = estimates(1);
        tau1 = estimates(2);
        A2 = estimates(3);
        tau2 = estimates(4);
        estimate_fit = (1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1) + A2*exp(-(t-t0)./tau2));
        estimate_fit(t < t0) = 0;
        estimate_fit(estimate_fit < 0) = 0;
        %compare fit, exclude if above sse criterion
        error_vector = estimate_fit(t > t0) - event_trace(t > t0);
        sse = sum(error_vector.^2);
        
        % want SSE criterion to be scaled by square root of peak
        est_peak_frame = event_frames(estimate_fit == max(estimate_fit));
        sse_criterion = sse_criterion/baseline_noEvent(est_peak_frame);
        
        if peak_amplitude > max_std*trace_std
            est_peak_frame = event_frames(estimate_fit == max(estimate_fit));
            spike_frames{roi_num} = [spike_frames{roi_num}; est_peak_frame];
            spike_amplitudes{roi_num} = [spike_amplitudes{roi_num}; baseline_noEvent(est_peak_frame)'];
            spike_t0s{roi_num} = [spike_t0s{roi_num}; t0_frame];
            residuals_all{roi_num} = [residuals_all{roi_num}; sse];
        elseif sse > sse_criterion || all(estimate_fit == 0) || sum(estimate_fit) < (1-area_criterion)*sum(event_trace) ...
                || sum(estimate_fit) > (1+area_criterion)*sum(event_trace);
            event_estimate_bad(event_frames) = estimate_fit;
            curr_point = curr_peak+1;
            continue
        else
            est_peak_frame = event_frames(estimate_fit == max(estimate_fit));
            spike_frames{roi_num} = [spike_frames{roi_num}; est_peak_frame];
            spike_amplitudes{roi_num} = [spike_amplitudes{roi_num}; baseline_noEvent(est_peak_frame)'];
            spike_t0s{roi_num} = [spike_t0s{roi_num}; t0_frame];
            residuals_all{roi_num} = [residuals_all{roi_num}; sse];
        end
        % extrapolate long event to remove (50 frames or as near)
        est_frames = local_frames(local_frames >= curr_peak-frame_back);
        t_est = (est_frames-est_frames(1)).*(1.24*128);
        estimate_fit_long = (1-exp(-(t_est-t0)./tau_on)).*(A1*exp(-(t_est-t0)./tau1) + A2*exp(-(t_est-t0)./tau2));
        estimate_fit_long(t_est < t0) = 0;
        estimate_fit_long(estimate_fit_long < 0) = 0;
        event_estimate(est_frames) = event_estimate(est_frames)+ estimate_fit_long;
        % subtract along template except for peak, place in event_detect,
        % correct for local baseline
        event_detect(est_frames) = (roi_trace_long(curr_roi,est_frames)-baseline_local(1)) ...
            - [0 estimate_fit_long(2:end)];
        % subtract events from baseline, store in baseline_noEvent, correct
        % for local baseline
        baseline_noEvent(est_frames) = ((baseline_noEvent(est_frames)) ...
            - estimate_fit_long); % keeps orig baseline
        event_localBaseline(est_frames) = roi_trace_long(curr_roi,est_frames)-baseline_local(1);
    
    else
        curr_point = curr_peak+1;
    end

end
    spike_estimates{roi_num} = event_estimate;
    disp(['Finished ROI ' num2str(curr_roi)]);
end










