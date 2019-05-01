% Subtract exponential decay from transients

% create template to subtract
% manual tau for now, single exponential
tau = 400; % tau in ms
tau_frame = tau./(1.24*128); % tau in frame
curr_roi = 3;
frame_search = 7;
frame_back = 1;
% template for entire event
event_template = exp(-(0:frame_search-1)/tau_frame);
% template leaving first point intact
decay_template = [0 event_template(2:end)];
% get total std
norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',2); % there are two hidden curves
baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve));
% initialize event and baseline vectors
peak_frames = roi_peak_max{curr_roi}(:,1);
event_detect = zeros(1,size(roi_trace_long,2));
baseline_noEvent = roi_trace_long(curr_roi,:);
event_estimate = zeros(1,size(roi_trace_long,2));
event_localBaseline = zeros(1,size(roi_trace_long,2));

w = waitbar(0,'Fitting events');
curr_point = 1;
while curr_point < length(baseline_noEvent)
    % detect next peak
    % 1) define criterion for detection
    % be > +1std +2std -1std -1std
    detect_min = [1*trace_std 2*trace_std -1*trace_std];
    % 2) find the next point that satisfies this criterion from peeled
    curr_detect = find([zeros(1,curr_point-1) diff(baseline_noEvent(curr_point:end))] > detect_min(1),1);
    if isempty(curr_detect) || ...
            curr_detect + length(detect_min) > length(baseline_noEvent);
        curr_point = curr_point + 1;
        continue
    else
        comp_detect = baseline_noEvent(curr_detect+1:curr_detect+length(detect_min))-baseline_noEvent(curr_detect);
        if all(comp_detect > detect_min)
            curr_peak = curr_detect+1;
                if  curr_peak+frame_search-1 < size(roi_trace_long,2) && ...
                    curr_peak- frame_back > 0;
                % amplitude of the current peak
                curr_peak_amp = roi_trace_long(curr_roi,curr_peak);
                % frames to compare event with Ca transient
                event_frames = [curr_peak-frame_back:curr_peak+frame_search-1];
                % correct for local baseline
                % get 50 frames, or as close as possible, for baseline estimation
                local_frames = [curr_peak-50:curr_peak+50];
                local_frames = local_frames(local_frames > 0 ...
                    & local_frames < size(roi_trace_long,2));
                % estimate local (100 frames) baseline to correct for baseline
                % drift by taking average of everything under abs(1std)
                try
                    norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,local_frames)',2);
                catch me
                    continue
                end
                baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
                baseline_local = norm_fit_obj.mu(baseline_curve);
                % create template for local-baseline event, fit with double exp
                event_trace = baseline_noEvent(event_frames)-baseline_local;%roi_trace_long(curr_roi,event_frames) - baseline_local;%
                % initial fit to get t0
                starting = [1 0 800];
                t = (0:frame_search+frame_back-1)*(1.24*128);
                options = [];
                estimates = fminsearch(@AP_fitExp,starting,options,t,event_trace);
                A = estimates(1);
                t0 = estimates(2);
                tau_on = estimates(3);
                % entire transient fit
                starting = [800,1,800,1,800];
                t = (0:frame_search+frame_back-1)*(1.24*128);
                options = [];
                estimates = fminsearch(@AP_fitdoubleExp,starting,options,t,event_trace,t0);
                tau_on = estimates(1);
                A1 = estimates(2);
                tau1 = estimates(3);
                A2 = estimates(4);
                tau2 = estimates(5);
                estimate_fit = (1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1 + A2*exp(-(t-t0)./tau2)));
                estimate_fit(t < t0) = 0;
                estimate_fit(estimate_fit < 0) = 0;
                % extrapolate long event to remove (50 frames or as near)
                est_frames = local_frames(local_frames >= curr_peak-frame_back);
                t_est = (est_frames-est_frames(1)).*(1.24*128);
                estimate_fit_long = (1-exp(-(t_est-t0)./tau_on)).*(A1*exp(-(t_est-t0)./tau1 + A2*exp(-(t_est-t0)./tau2)));
                estimate_fit_long(t_est < t0) = 0;
                estimate_fit_long(estimate_fit_long < 0) = 0;
                event_estimate(est_frames) = estimate_fit_long;
                % subtract along template except for peak, place in event_detect,
                % correct for local baseline
                event_detect(est_frames) = (roi_trace_long(curr_roi,est_frames)-baseline_local) ...
                    - [0 estimate_fit_long(2:end)];%decay_template.*(roi_trace_long(curr_roi,event_frames-1)-baseline_local);
                % subtract events from baseline, store in baseline_noEvent, correct
                % for local baseline
                baseline_noEvent(est_frames) = ((baseline_noEvent(est_frames)) ...
                    - estimate_fit_long + baseline_local); % to keep orig
                event_localBaseline(est_frames) = roi_trace_long(curr_roi,est_frames)-baseline_local;
                end
                curr_point = curr_peak;
        else
            curr_point = curr_detect + 1;
            continue
        end
        waitbar(curr_point/length(baseline_noEvent),w,'Fitting events');
    end
end
    close(w);
    
    
    
    
    
    
    
    
    
    
    
