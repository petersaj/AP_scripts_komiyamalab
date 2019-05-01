% Subtract exponential decay from transients

% set parameters
curr_roi = 81; 
sse_criterion = 0.1; % high threshold for sse acceptance
area_criterion = 0.4; % percent of area fit cannot deviate from
frame_search = 6; % frames to look past peak
frame_back = 3; % frames to look behind peak

% get total std
norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',2); % there are two hidden curves
baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve));
% initialize event and baseline vectors
peak_frames = roi_peak_max{curr_roi}(:,1);
event_detect = zeros(1,size(roi_trace_long,2));
baseline_noEvent = roi_trace_long(curr_roi,:);
event_estimate = zeros(1,size(roi_trace_long,2));
event_estimate_bad = zeros(1,size(roi_trace_long,2));
event_localBaseline = zeros(1,size(roi_trace_long,2));
residuals_all = zeros(1,length(peak_frames));
spike_frames = [];

w = waitbar(0,'Fitting events');
curr_point = 1;
for curr_peak = 1:length(peak_frames)
    if  peak_frames(curr_peak)+frame_search-1 < size(roi_trace_long,2) && ...
            peak_frames(curr_peak)- frame_back > 0;
        if curr_peak == 118
            pause(0.1)
        end
        % amplitude of the current peak
        curr_peak_amp = roi_trace_long(curr_roi,peak_frames(curr_peak));
        % frames to compare event with Ca transient
        event_frames = [peak_frames(curr_peak)-frame_back:peak_frames(curr_peak)+frame_search-1];
        % correct for local baseline
        % get 50 frames, or as close as possible, for baseline estimation
        local_frames = [peak_frames(curr_peak)-50:peak_frames(curr_peak)+50];
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
        event_trace = baseline_noEvent(event_frames)-baseline_local;  %roi_trace_long(curr_roi,event_frames) - baseline_local;
        % initial fit to get t0
        A_guess = max(event_trace);
        starting = [A_guess 200 30 600];
        t = (0:frame_search+frame_back-1)*(1.24*128);
        
        %%% find fit
        
        options = optimset('Algorithm','interior-point');
        lb = [0 -frame_back*1.24*128 0 0];
        ub = [2*A_guess frame_back*1.24*128 300 2000];
        % estimate transform of orig -> curr
        display 'Finding best fit...'
        estimates = fmincon(@(x) AP_fitExp(x,t,event_trace),starting,[],[],[],[],lb,ub,[],options);
        A = estimates(1);
        t0 = estimates(2);
        tau_on = estimates(3);
        tau1 = estimates(4);
        % entire transient fit
        starting = [A_guess,300,A_guess,800];%,t0,tau_on];
        t = (0:frame_search+frame_back-1)*(1.24*128);
        % embed values just gotten to pass to function
        lb = [0,0,0,0];%,-frame_back*1.24*128,0];
        ub = [5,2000,5,2000];%,1000,500];
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
        residuals_all(curr_peak) = sse;
        if sse > sse_criterion || all(estimate_fit == 0) || sum(estimate_fit) < (1-area_criterion)*sum(event_trace) ...
                || sum(estimate_fit) > (1+area_criterion)*sum(event_trace);
            event_estimate_bad(event_frames) = estimate_fit;
            continue
        else
            est_peak_frame = event_frames(estimate_fit == max(estimate_fit));
            spike_frames = [spike_frames; est_peak_frame];
        end
        % extrapolate long event to remove (50 frames or as near)
        est_frames = local_frames(local_frames >= peak_frames(curr_peak)-frame_back);
        t_est = (est_frames-est_frames(1)).*(1.24*128);
        estimate_fit_long = (1-exp(-(t_est-t0)./tau_on)).*(A1*exp(-(t_est-t0)./tau1) + A2*exp(-(t_est-t0)./tau2));
        estimate_fit_long(t_est < t0) = 0;
        estimate_fit_long(estimate_fit_long < 0) = 0;
        event_estimate(est_frames) = event_estimate(est_frames)+ estimate_fit_long;
        % subtract along template except for peak, place in event_detect,
        % correct for local baseline
        event_detect(est_frames) = (roi_trace_long(curr_roi,est_frames)-baseline_local) ...
            - [0 estimate_fit_long(2:end)];%decay_template.*(roi_trace_long(curr_roi,event_frames-1)-baseline_local);
        % subtract events from baseline, store in baseline_noEvent, correct
        % for local baseline
        baseline_noEvent(est_frames) = ((baseline_noEvent(est_frames)) ...
            - estimate_fit_long); % keeps orig baseline
        event_localBaseline(est_frames) = roi_trace_long(curr_roi,est_frames)-baseline_local;
    end
    waitbar(curr_peak/length(peak_frames),w,'Fitting events');
end
close(w);











