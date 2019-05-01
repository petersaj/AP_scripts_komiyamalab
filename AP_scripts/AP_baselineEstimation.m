function [roi_trace_df roi_trace_baseline] = AP_baselineEstimation(roi_trace_long,framerate)
% [roi_trace_df roi_trace_baseline] = AP_baselineEstimation(roi_trace_long,framerate)
% roi_trace_df: the df/f normalized trace
% roi_trace_baseline (optional): the estimated baseline of the trace

disp('Baseline correcting traces (%):')
fprintf('%2d',0);

baseline_est_runs = 2;

[num_rois num_frames] = size(roi_trace_long);
roi_trace_df = nan(size(roi_trace_long));
roi_trace_baseline = nan(size(roi_trace_long));
for curr_cell = 1:num_rois;
        
    curr_trace = roi_trace_long(curr_cell,:);
    
    % if there's a NaN in the current trace, which happens for instance
    % when the ROI is off of the screen during this day, ignore this ROI
    if any(isnan(curr_trace));
        roi_trace_df(curr_cell,:) = nan(1,size(roi_trace_df,2));
        continue
    end
    
    % Estimate baseline, get rid of active zones, run baseine again.
    for baseline_est = 1:baseline_est_runs;
        
        % 1) Take long moving average to get general moving baseline, and
        % produce the first df estimate.
        
        % amount to smooth for baseline estimation, in minutes
        baseline_smooth_min = 1;
        baseline_smooth_frames = round(baseline_smooth_min*60*framerate);
        if size(curr_trace,2) <= baseline_smooth_frames%TK: if too short, consider 1/3 as window
            baseline_smooth_frames = round(size(curr_trace,2) / 3);
        end
        % smooth by moving average to flatten activity
        baseline_smooth = ...
            smooth(curr_trace,baseline_smooth_frames)';
        % copy the ends for reliable smooth points, not near edges
        smooth_reliable = ceil(baseline_smooth_frames/2);
        baseline_smooth(1:smooth_reliable) = ...
            baseline_smooth(smooth_reliable);
        baseline_smooth(end-smooth_reliable:end) = ...
            baseline_smooth(end-smooth_reliable);
        
        % create initial df estimate to find the mode
        baseline_smooth_df = curr_trace - baseline_smooth;
        
        % 2) Find mode or gaussian fit of loess smoothed trace and get the
        % average noise by raw-smoothed
        
        % amount to loess smooth, usually one second looks pretty good
        % -------> actually, it was that 30 frames looks pretty good,
        % regardless of framerate (which makes sense)
        noise_smooth_sec = 1;
        %noise_smooth_frames = round(noise_smooth_sec*framerate);
        noise_smooth_frames = 30;
        % loess smooth on the original trace
        noise_smooth = ...
            smooth(curr_trace,noise_smooth_frames,'loess')';
        % loess smooth on the df estimate
        noise_smooth_df= ...
            smooth(baseline_smooth_df,noise_smooth_frames,'loess')';
        % estimate baseline by mode (since drifts should already be fixed)
        noise_smooth_mode = mode(round(noise_smooth_df));
        
        % try to find baseline by gaussian fitting trace (even though
        % usually not gaussian, and definitely not on the first run)
        norm_fit_obj = gmdistribution.fit(noise_smooth_df',1);
        dx = (max(noise_smooth_df')-min(noise_smooth_df'))/length(noise_smooth_df');
        x = min(noise_smooth_df'):dx:max(noise_smooth_df');
        x = x(1:end-1); % it adds an extra
        norm_fit = pdf(norm_fit_obj,x');
        noise_smooth_mode2 = x(norm_fit == max(norm_fit));
        if abs(noise_smooth_mode - noise_smooth_mode2)<100
            noise_smooth_mode = noise_smooth_mode2;
        end
        
        noise_std = std(baseline_smooth_df-noise_smooth_df);
        
        % 3) Get rid of active times, and re-estimate the baseline.
        % Take where the trace goes above and below some lever as activity in some
        % direction, and restimate mode ignoring things around these points.
            
        % take 2x std above and below as "activity" bounds. if this is one
        % of the preliminary estimates, then only worry about things going
        % over threshold.
        lb = noise_smooth_mode + baseline_smooth - 2*noise_std;
        ub = noise_smooth_mode + baseline_smooth + 2*noise_std;
        if baseline_est < baseline_est_runs
            bound_indx = noise_smooth > ub;
        else
            bound_indx = noise_smooth > ub | noise_smooth < lb;
        end
        % index everything within n seconds of active zones
        activity_surround_sec = 5;
        activity_surround_frames = round(activity_surround_sec*framerate);
        % do this by convolving with filter the size of surround, taking >0's
        activity_surround_kernel = ones(1,round(activity_surround_frames/2));
        active_surround_frames = conv2(+bound_indx,activity_surround_kernel,'same');
        active_indx_surround_curr = bound_indx;
        active_indx_surround_curr(active_surround_frames > 0) = true;
        
        % define which frames are being looked at - first iter is all
        % frames, also initialize the full active frame index
        if baseline_est == 1
            curr_trace_indx = [1:length(curr_trace)];
            active_indx_surround = false(1,size(roi_trace_long,2));
        end
        
        % update the final index for active frames
        active_indx_surround(curr_trace_indx) = active_indx_surround_curr;
                   
        % update which frames will be looked at in the next iteration
        curr_trace_indx = curr_trace_indx(~active_indx_surround_curr);
       
        % truncate out the active parts, use that for next iteration    
        inactive_trace_trunk = curr_trace(~active_indx_surround_curr);
        curr_trace = inactive_trace_trunk;
                        
%         set the current trace to be only inactive parts with interpolated
%         fluorescence between them, use this for future estimations
%         inactive_trace = NaN(size(curr_trace));
%         inactive_trace(~active_indx_surround) = ...
%             curr_trace(~active_indx_surround);
%         make first and last frames equal first and last values
%         inactive_trace(1) = ...
%             curr_trace(find(~active_indx_surround,1));
%         inactive_trace(end) = ...
%             curr_trace(find(~active_indx_surround,1,'last'));
%         active_indx_surround(1) = 0;
%         active_indx_surround(end) = 0;
%         interpolate between active regions
%         inactive_trace =  ...
%             interp1(find(~active_indx_surround), ...
%             inactive_trace(~active_indx_surround), ...
%             [1:length(curr_trace)]);
%         curr_trace = inactive_trace;
        
    end
    
    % if (incredibly rare - and the cells are probably filled anyway)
    % active 95% of the time, skip the cell
    if length(curr_trace) <= 0.05*length(roi_trace_long(curr_cell,:));
        roi_trace_df(curr_cell,:) = nan(1,size(roi_trace_df,2));
        continue
    end
    
    % the current trace is once again the raw trace
    curr_trace = roi_trace_long(curr_cell,:);
    
    whole_noise_smooth = ...
        smooth(curr_trace,noise_smooth_frames,'loess')';
    
    % with the active parts indexed, fill a vector with only inactive
    % frames smoothed by previous large-scale amount    
    inactive_baseline_smooth_part = [];
    inactive_baseline_smooth_part = ...
        smooth(whole_noise_smooth(~active_indx_surround), ...
        baseline_smooth_frames);
    
    % copy the ends for reliable smooth points, not near edges
    smooth_reliable = ceil(baseline_smooth_frames/2);
%      if size(baseline_smooth,2) <= smooth_reliable%TK: if too short, consider 1/3 as edge
%          smooth_reliable = round(size(baseline_smooth,2) / 4);
%      end

     if smooth_reliable >= length(inactive_baseline_smooth_part)%TK: if too short, consider 1/3 as edge
         smooth_reliable = round(size(inactive_baseline_smooth_part,1) / 3);
     end
    inactive_baseline_smooth_part(1:smooth_reliable) = ...
        inactive_baseline_smooth_part(smooth_reliable);
    inactive_baseline_smooth_part(end-smooth_reliable:end) = ...
        inactive_baseline_smooth_part(end-smooth_reliable);
    inactive_baseline_smooth = NaN(size(curr_trace));
    inactive_baseline_smooth(~active_indx_surround) = ...
        inactive_baseline_smooth_part;
     
    % make the first and last frames equal the first and last values
    inactive_baseline_smooth(1) = inactive_baseline_smooth_part(1);
    inactive_baseline_smooth(end) = inactive_baseline_smooth_part(end);
    active_indx_surround(1) = 0;
    active_indx_surround(end) = 0;
    
    % interpolate between active sites
    inactive_baseline_full = [];
    inactive_baseline_full = ...
        interp1(find(~active_indx_surround), ...
        inactive_baseline_smooth(~active_indx_surround), ...
        [1:length(inactive_baseline_smooth)]);
    
    % make another df estimate based on inactive baseline
    inactive_df = curr_trace - inactive_baseline_full;
    
    % 4) At this point, all drift should be compensated for and all active
    % zones should have been ignored for baselinining. Do one final global
    % baseline estimation by again finding active zones and finding the
    % mode as the baseline.
    
    inactive_df_smooth = smooth(inactive_df,noise_smooth_frames,'loess')';
    noise_std2 = std(abs(inactive_df - inactive_df_smooth));
    lb = -2*noise_std2;
    ub = 2*noise_std2;
    baseline_indx = inactive_df_smooth > ub | inactive_df_smooth < lb;
    
    % index everything within n seconds of active zones
    activity_surround_sec = 5;
    activity_surround_frames = round(activity_surround_sec*framerate);
    % do this by convolving with filter the size of surround, taking >0's
    activity_surround_kernel = ones(1,round(activity_surround_frames/2));
    active_surround_frames = conv2(+baseline_indx,activity_surround_kernel,'same');
    active_indx_surround = baseline_indx;
    active_indx_surround(active_surround_frames > 0) = true;
    
    baseline_offset = ...
        mode(round(inactive_df_smooth(~active_indx_surround)));
    
    % if (incredibly rare - and the cells are probably filled anyway)
    % active 95% of the time, skip the cell
    if sum(active_indx_surround) >= 0.95*length(active_indx_surround);
        roi_trace_df(curr_cell,:) = nan(1,size(roi_trace_df,2));
        continue
    end
    
    % try gaussian fitting (it really should be a gaussian by now)
    norm_fit_obj = gmdistribution.fit(inactive_df_smooth(~active_indx_surround)',1);
    dx = (max(inactive_df_smooth(~active_indx_surround)')-min(inactive_df_smooth(~active_indx_surround)'))/length(inactive_df_smooth(~active_indx_surround)');
    x = min(inactive_df_smooth(~active_indx_surround)'):dx:max(inactive_df_smooth(~active_indx_surround)');
    x = x(1:end-1); % it adds an extra
    norm_fit = pdf(norm_fit_obj,x');
    baseline_offset2 = x(norm_fit == max(norm_fit));
    if abs(baseline_offset2 - baseline_offset)<100
        baseline_offset = baseline_offset2;
    end
    
        
    % 5) Create the final baseline trace, then do the final normalization
    
    final_baseline = inactive_baseline_full+baseline_offset;
    roi_trace_df(curr_cell,:) = ...
        (roi_trace_long(curr_cell,:) - final_baseline) ./ final_baseline;
    roi_trace_baseline(curr_cell,:) = final_baseline;
    
    fprintf('%c%c%2d',8,8,floor(100*curr_cell/num_rois));

end

fprintf('%c%c%c',8,8,8);
disp('Finished.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Everything below this point was an old attempt that wasn't as good

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % Calculate moving baseline by finding gaussian fit within window
% tic
% 
% baseline_window = 6000;
% baseline_window_step = 2000;
% min_baseline_length = 300;
% min_below_thresh = 100;
% max_consec = 20;
% size_file = 4000;
% 
% size_smooth = 90;
% smooth_type = 'loess';
% 
% 
% % This was to do it loop by loop, ignoring for now
% % Estimate noise through smoothing
% roi_trace_long_smooth = zeros(size(roi_trace_long));
% for i = 1:size(roi_trace_long,1);
%     for j = 1:size(roi_trace_long,2)/size_file;
%         roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
%     end
% end
% smooth_std = zeros(size(roi_trace_long,1),1);
% smooth_std = sqrt(sum((abs(roi_trace_long - roi_trace_long_smooth).^2)/size(roi_trace_long_smooth,2),2));
% 
% % Estimate baseline: threshold trace, mirror below threshold to above
% % threshold, match standard deviation with estimated noise
% roi_trace_df = zeros(size(roi_trace_long));
% roi_trace_smooth_df = zeros(size(roi_trace_long));
% 
% % make a baseline matrix: all points not estimated are NaN
% estimate_baseline_all = nan(size(roi_trace_long));
% baseline_mode_all = nan(size(roi_trace_long));
% baseline_trace_indx_all = cell(size(roi_trace_long,1),1);
% 
% for roi = 1:size(roi_trace_long,1);
%     
%     curr_frame = baseline_window/2+1;
%     
%     while curr_frame < size(roi_trace_long,2)-baseline_window/2;
%         
%         curr_frame_range = [curr_frame-baseline_window/2:curr_frame+baseline_window/2];
%         curr_trace = roi_trace_long(roi,curr_frame_range);
%         
%         options = optimset('Algorithm','interior-point');
%         options = optimset('MaxFunEvals',2000);
%         lb = min(roi_trace_long(roi,:));
%         ub = max(roi_trace_long(roi,:));
%         starting = lb+(ub-lb)/2;
%         
%         %     % this function minimizes the absolute difference between the estimated
%         %     % standard deviation, and the standard deviation of the mirrored trace
%         %     % below threshold
%         %     estimate_baseline = fminsearch(@(x) ...
%         %         abs(smooth_std(i)- ...
%         %         std([(roi_trace_long(i,(roi_trace_long(i,:) < x))-x) ...
%         %         -(roi_trace_long(i,(roi_trace_long(i,:) < x))-x)])),starting);
%         %
%         %%% THIS IS FOR THE MODE RESTIMATION
%         % this function minimizes the absolute difference between the estimated
%         % standard deviation, and the standard deviation of the mirrored trace
%         % below threshold
%         estimate_baseline = fminsearch(@(x) ...
%             AP_fitBaselineGaussian(x,smooth_std(roi),curr_trace),starting);
%         
%         % correct the baseline by looking for the mode of the cutoff assumed to
%         % be just baseline (get better mode estimate by rounding baseline_trace
%         % to the nearest interger)
%         baseline_ub = estimate_baseline + 2*smooth_std(roi);
%         baseline_lb = estimate_baseline - 2*smooth_std(roi);
%         baseline_trace = curr_trace(curr_trace > baseline_lb &...
%             curr_trace < baseline_ub);
%         baseline_trace_indx = find(curr_trace > baseline_lb &...
%             curr_trace < baseline_ub) + curr_frame-baseline_window/2;
%         baseline_mode = mode(round(baseline_trace));
%         
%         baseline_trace_below_indx = find(curr_trace < baseline_mode);
%         
%         % find longest consecutive run
%         consec_below = [0 diff(baseline_trace_below_indx) == 1 0];
%         p = find(~consec_below);%find 0
%         longest_consec = max(diff(p)-1);
%         
%         % 1) ignore if number of baseline points below minimum
%         % 2) ignore if based off of too many consecutive samples below thresh
%         % 3) ignore if number of below threshold points is below min
%         if length(baseline_trace) > min_baseline_length &&...
%                 longest_consec < max_consec && ...
%                 length(curr_trace(curr_trace < estimate_baseline)) > min_below_thresh
%             baseline_mode_all(roi,curr_frame) = baseline_mode;
%             estimate_baseline_all(roi,curr_frame) = estimate_baseline;
%             baseline_trace_indx_all{roi} = [baseline_trace_indx_all{roi} baseline_trace_indx];
%         end
%                 
%         curr_frame = curr_frame + baseline_window_step;
%         
%     end
% end
% 
% % if no good estimate, then just try to estimate by the whole trace
% for roi = 1:size(roi_trace_long,1);
%     if ~all(isnan(estimate_baseline_all(roi,:)))
%         continue
%     end
%     disp(['WARNING: No good baseline for ' num2str(roi) ', estimating from whole trace']);
%      options = optimset('Algorithm','interior-point');
%         options = optimset('MaxFunEvals',2000);
%         lb = min(roi_trace_long(roi,:));
%         ub = max(roi_trace_long(roi,:));
%         starting = lb+(ub-lb)/2;
%         
%         curr_trace = roi_trace_long(roi,:);
%         curr_frame = round(size(roi_trace_long,2)/2);
%         %%% THIS IS FOR THE MODE RESTIMATION
%         % this function minimizes the absolute difference between the estimated
%         % standard deviation, and the standard deviation of the mirrored trace
%         % below threshold
%         estimate_baseline = fminsearch(@(x) ...
%             AP_fitBaselineGaussian(x,smooth_std(roi),curr_trace),starting);
%         
%         % correct the baseline by looking for the mode of the cutoff assumed to
%         % be just baseline (get better mode estimate by rounding baseline_trace
%         % to the nearest interger)
%         baseline_ub = estimate_baseline + 2*smooth_std(roi);
%         baseline_lb = estimate_baseline - 2*smooth_std(roi);
%         baseline_trace = curr_trace(curr_trace > baseline_lb &...
%             curr_trace < baseline_ub);
%         baseline_trace_indx = find(curr_trace > baseline_lb &...
%             curr_trace < baseline_ub);
%         baseline_trace_indx_all{roi} = [baseline_trace_indx_all{roi} baseline_trace_indx];
%         baseline_mode = mode(round(baseline_trace));
%         baseline_mode_all(roi,curr_frame) = baseline_mode;
%         estimate_baseline_all(roi,curr_frame) = estimate_baseline;
% end
% 
% % interpolate baseline between estimated baseline points
% estimate_baseline_discrete = estimate_baseline_all;
% for roi = 1:size(roi_trace_long,1);
%     if all(isnan(estimate_baseline_all(roi,:)))
%         disp(['WARNING: No good estimate for ROI ' num2str(roi)]);
%         continue
%     end
%     frame_estimates_x = [];
%     frame_estimates_x = find(isnan(estimate_baseline_all(roi,:)) == 0);
%     frame_estimates_y = [];
%     frame_estimates_y = estimate_baseline_all(roi,frame_estimates_x);
%     % assume baseline start and end are the same as the first and last
%     % estimates of baseline
%     frame_estimates_x = [1 frame_estimates_x size(roi_trace_long,2)];
%     frame_estimates_y = [frame_estimates_y(1) frame_estimates_y frame_estimates_y(end)];
%     estimate_baseline_all(roi,:) = ...
%         interp1(frame_estimates_x,frame_estimates_y,1:size(roi_trace_long,2));
% end
% 
% % interpolate baseline between estimated baseline mode points
% baseline_mode_discrete = baseline_mode_all;
% for roi = 1:size(roi_trace_long,1);
%     if all(isnan(baseline_mode_all(roi,:)))
%         continue
%     end
%     frame_estimates_x = [];
%     frame_estimates_x = find(isnan(baseline_mode_all(roi,:)) == 0);
%     frame_estimates_y = [];
%     frame_estimates_y = baseline_mode_all(roi,frame_estimates_x);
%     % assume baseline start and end are the same as the first and last
%     % estimates of baseline
%     frame_estimates_x = [1 frame_estimates_x size(roi_trace_long,2)];
%     frame_estimates_y = [frame_estimates_y(1) frame_estimates_y frame_estimates_y(end)];
%     baseline_mode_all(roi,:) = ...
%         interp1(frame_estimates_x,frame_estimates_y,1:size(roi_trace_long,2));
% end
% 
% % this is to actually correct the raw trace
% roi_trace_df = (roi_trace_long-baseline_mode_all)./baseline_mode_all;
% roi_trace_smooth_df = (roi_trace_long_smooth-baseline_mode_all)./baseline_mode_all;
% 
% toc

