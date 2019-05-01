function caEvents = AP_caEvents_thresh(df_matrix,thresh,method)
% caEvent_matrix = AP_caEvents_thresh(df_matrix,thresh,type)
% Create a matrix of calcium events from raw df/f traces
% 'thresh' - noise multiplier threshold
%
% method: 
% 0 (default) = zero inactive portions of the trace
% 1 = zero the trace when df/f is falling
% 2 = retain max-min df of increasing active df/f portions

if nargin < 3
    method = 0;
end

num_cells = size(df_matrix,1);

% Discard cells that contain all NaNs and zero occasional NaNs (is this ok?)
discard_cells = false(num_cells,1);
discard_cells(all(isnan(df_matrix),2)) = true;
df_matrix(isnan(df_matrix)) = 0;

use_cells = find(~discard_cells);

% Initialize calcium events matrix
caEvents = nan(size(df_matrix));

for curr_cell = use_cells'
    
    curr_trace = df_matrix(curr_cell,:);
    
    temp_smooth_30 = smooth(curr_trace,30,'loess')';
    
    %noise_est = mean(abs(df_matrix(curr_cell,:) - temp_smooth));
    % NEW NOISE EST: the std of (mirrored) below-zero df/f values
    below_zero_trace = curr_trace(curr_trace < 0);
    noise_est = std([below_zero_trace -below_zero_trace]);
    
    thresh_lo = temp_smooth_30 > noise_est*1;
    thresh_hi = temp_smooth_30 > noise_est*thresh;
    
    % fill in hi-threshold portions where lo threshold is not crossed (dips
    % down but not to baseline, so continuously active portion)
    
    % find edges of long-smooth above-thresh periods
    thresh_lo_start = diff([0 thresh_lo 0]) == 1;
    thresh_lo_stop = diff([0 thresh_lo 0]) == -1;
    thresh_hi_start = diff([0 thresh_hi 0]) == 1;
    thresh_hi_stop = diff([0 thresh_hi 0]) == -1;
    
    thresh_hi_start_idx = find(thresh_hi_start);
    thresh_hi_stop_idx = find(thresh_hi_stop);
    
    thresh_lo_hi_smooth_idx = arrayfun(@(x,y) ...
        x-find(thresh_lo_start(x:-1:1),1)+1:y, ...
        thresh_hi_start_idx,thresh_hi_stop_idx,'uni',false);
    
    % don't bother trying to find activity that starts before imaging or
    % continues after last frame
    exclude_activity = cellfun(@(x) any(x(x <= 1 | ...
        x >= length(curr_trace))),thresh_lo_hi_smooth_idx);
    
    % refine start times of activity to when df/f goes above hi thresh
    thresh_lo_hi_raw_idx = cellfun(@(x) x(1) +...
        find(curr_trace(x:end) > noise_est*thresh,1) - ...
        find(curr_trace(x(1)+find(curr_trace(x:end) > noise_est*thresh,1):-1:1) < noise_est*1,1): ...
        x(end),thresh_lo_hi_smooth_idx(~exclude_activity),'uni',false);
    
    % again, filter out events that overlap with beginning or end
    exclude_activity_2 = cellfun(@(x) any(x(x <= 1 | ...
        x >= length(curr_trace))),thresh_lo_hi_raw_idx);
    thresh_lo_hi_raw_idx(exclude_activity_2) = [];
    
    % find continuous active portions after this process      
    active_trace = zeros(size(curr_trace));
    active_trace(horzcat(thresh_lo_hi_raw_idx{:})) = true;
    
    active_idx = arrayfun(@(x,y) x:y,find(diff([0 active_trace 0]) == 1), ...
        find(diff([0 active_trace 0]) == -1),'uni',false);
    
    % Create trace from chosen method
    switch method
        % 0 (default) = retain over-threshold, zero others
        case 0
            over_thresh_frames = horzcat(active_idx{:});
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = curr_trace(over_thresh_frames);
            
        % 1 = retain over-threshold & increasing, zero others
        case 1
            thresh_lo_hi_rising = cellfun(@(x) x(diff([0 temp_smooth_30(x)]) > 0), ...
                active_idx,'uni',false);
            
            over_thresh_frames = horzcat(thresh_lo_hi_rising{:});
            
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = curr_trace(:,over_thresh_frames);
            
        % 2 = retain ddf of over-threshold & increasing portions
        case 2
            % find rising portions
            thresh_lo_hi_rising = cellfun(@(x) x(diff([0 temp_smooth_30(x)]) > 0), ...
                active_idx,'uni',false);
            % split rising portions
            active_rising_trace = zeros(size(curr_trace));
            active_rising_trace(horzcat(thresh_lo_hi_rising{:})) = true;            
            active_rising_idx = arrayfun(@(x,y) x:y,find(diff([0 active_rising_trace 0]) == 1), ...
                find(diff([0 active_rising_trace 0]) == -1),'uni',false);
    
            thresh_lo_hi_rising_ddf = cellfun(@(x) repmat((max(curr_trace(x)) - ...
                min(curr_trace(x))),1,length(x)),active_rising_idx,'uni',false);
            
            over_thresh_frames = horzcat(active_rising_idx{:});
            over_thresh_values = horzcat(thresh_lo_hi_rising_ddf{:});
            
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = over_thresh_values;
    end    
end










