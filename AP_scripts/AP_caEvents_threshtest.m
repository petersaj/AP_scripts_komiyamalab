function [caEventIdxs caEvents caStds] = AP_caEvents_threshtest(df_matrix,thresh)
% [caEventIdxs caEvents caStds] = AP_caEvents_threshtest(df_matrix,thresh)
% Create a matrix of calcium events from raw df/f traces
% 'thresh' - noise multiplier threshold
%
% THIS IS A FUNCTION TO TEST VARIOUS THRESHOLDS: THIS IS NOT INTENDED FOR
% ACTUAL DATA ANALYSIS

num_cells = size(df_matrix,1);

% Discard cells that contain all NaNs and zero occasional NaNs (is this ok?)
discard_cells = false(num_cells,1);
discard_cells(all(isnan(df_matrix),2)) = true;
df_matrix(isnan(df_matrix)) = 0;

use_cells = find(~discard_cells);

% Initialize calcium events matrix
caEventIdxs = cell(num_cells,1);
caEvents = nan(size(df_matrix));
caStds = cell(num_cells,1);
for curr_cell = use_cells'
    
    curr_trace = df_matrix(curr_cell,:);
    
    temp_smooth_30 = smooth(curr_trace,30,'loess')';
    
    %noise_est = mean(abs(df_matrix(curr_cell,:) - temp_smooth));
    % NEW NOISE EST: the std of (mirrored) below-zero df/f values
    below_zero_trace = curr_trace(curr_trace < 0);
    noise_est = std([below_zero_trace -below_zero_trace]);
    
    thresh_lo = temp_smooth_30 > noise_est*1;
    thresh_hi = temp_smooth_30 > noise_est*thresh;
    
    % find edges of long-smooth above-thresh periods
    thresh_lo_start = diff([0 thresh_lo 0]) == 1;
    thresh_lo_stop = diff([0 thresh_lo 0]) == -1;
    thresh_hi_start = find(diff([0 thresh_hi 0]) == 1);
    
    thresh_lo_hi_smooth_idx = arrayfun(@(x) x-find(thresh_lo_start(x:-1:1),1)+1: ...
        x+find(thresh_lo_stop(x:end),1)-1, ...
        thresh_hi_start,'uni',false);
    
    % exclude activity that overlaps with start/end of imaging
    exclude_activity_bounds = cellfun(@(x) any(x < 1 | x > ...
        length(curr_trace)),thresh_lo_hi_smooth_idx);
    
    % refine start times of activity to when df/f goes about hi thresh
    thresh_lo_hi_raw_idx = cellfun(@(x) x(1) +...
        find(curr_trace(x:end) > noise_est*thresh,1) - ...
        find(curr_trace(x(1)+find(curr_trace(x:end) > noise_est*thresh,1):-1:1) < noise_est*1,1): ...
        x(end),thresh_lo_hi_smooth_idx(~exclude_activity_bounds),'uni',false);
    
    % exclude cases where there active frames are active (happens when
    % activity is close to beginning/end when doesn't have a bordering
    % baseline time)
    exclude_activity = cellfun(@(x) any(x(x <= 1 | ...
        x >= length(curr_trace))) | isempty(x),thresh_lo_hi_raw_idx);
    
    thresh_lo_hi_raw_idx_use = thresh_lo_hi_raw_idx(~exclude_activity);
    
    % get the maximum point (in noise stds) for each active portion
    caStds{curr_cell} = cellfun(@(x) max(temp_smooth_30(x))./noise_est, ...
        thresh_lo_hi_raw_idx_use);

    over_thresh_frames = horzcat(thresh_lo_hi_raw_idx_use{:});
    
    thresh_lo_hi = false(1,size(caEvents,2));
    thresh_lo_hi(over_thresh_frames) = true;
        
    caEventIdxs{curr_cell} = thresh_lo_hi_raw_idx_use;
    
    caEvents(curr_cell,~thresh_lo_hi) = 0;
    caEvents(curr_cell,thresh_lo_hi) = curr_trace(:,thresh_lo_hi);
end










