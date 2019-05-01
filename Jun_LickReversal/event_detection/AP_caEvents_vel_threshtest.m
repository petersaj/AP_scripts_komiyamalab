function [caEventIdxs caEvents caStds_vel caStds_df] = ...
    AP_caEvents_vel_threshtest(df_matrix,thresh,use_cells)
% [caEventIdxs caEvents caStds] = AP_caEvents_threshtest(df_matrix,thresh)
% Create a matrix of calcium events from raw df/f traces
% 'thresh' - noise multiplier threshold * NOW OF VELOCITY
%
% THIS IS A FUNCTION TO TEST VARIOUS THRESHOLDS: THIS IS NOT INTENDED FOR
% ACTUAL DATA ANALYSIS
%
% New version (vel): 
% test VELOCITY thresholds 

num_cells = size(df_matrix,1);

% Discard cells that contain all NaNs and zero occasional NaNs (is this ok?)
discard_cells = false(num_cells,1);
discard_cells(all(isnan(df_matrix),2)) = true;
% Discard cells instructed not to use
discard_cells(~use_cells) = true;

df_matrix(isnan(df_matrix)) = 0;

use_cells = find(~discard_cells);

num_cells = size(df_matrix,1);

% Smooth traces
size_smooth = 30;
smooth_type = 'loess';
df_matrix_smooth = nan(size(df_matrix));
for curr_cell = 1:num_cells
    df_matrix_smooth(curr_cell,:) = smooth(df_matrix(curr_cell,:),size_smooth,smooth_type);
end

% Get first and second derivatives for detecting peaks
ddf_matrix_smooth = [zeros(num_cells,1) diff(df_matrix_smooth,[],2)];
d2df_matrix_smooth = [zeros(num_cells,1) diff(ddf_matrix_smooth,[],2)];

% Get noise estimates for each cell based on smoothing
noise_std = nan(num_cells,1);
noise_thresh = nan(num_cells,1);
for curr_cell = 1:num_cells
    noise_std(curr_cell) = std([df_matrix(curr_cell,df_matrix(curr_cell,:) < 0) ...
        -df_matrix(curr_cell,df_matrix(curr_cell,:) < 0)]);
    noise_thresh(curr_cell) = 2*noise_std(curr_cell);
end

% Initialize calcium events matrix
caEventIdxs = cell(num_cells,1);
caEvents = nan(size(df_matrix));
caStds_vel = cell(num_cells,1);
caStds_df = cell(num_cells,1);
for curr_cell = use_cells'
    
    % Get velocity noise for current cell
    velocity_noise = std(ddf_matrix_smooth(curr_cell, ...
        df_matrix_smooth(curr_cell,:) < noise_thresh(curr_cell)));
    velocity_noise_thresh = thresh*velocity_noise;
    
    % Find rising edges of where that threshold is crossed
    velocity_thresh_rise = ...
        [0 diff(sign(ddf_matrix_smooth(curr_cell,:) - velocity_noise_thresh))] == 2 & ...
        d2df_matrix_smooth(curr_cell,:) > 0;
    
    % Find the peak values associated with threshold crossings
    zero_vel_down = [0 diff(sign(ddf_matrix_smooth(curr_cell,:)))] == -2 & ...
        d2df_matrix_smooth(curr_cell,:) < 0;
    zero_vel_up = [0 diff(sign(ddf_matrix_smooth(curr_cell,:)))] == 2 & ...
        d2df_matrix_smooth(curr_cell,:) > 0;
    
    % Index the fast rises
    velocity_thresh_rise_idx = find(velocity_thresh_rise);
            
    % Ignore velocity threshold crossings at end not followed by peak or if
    % there isn't a previous zero-velocity crossing
    velocity_thresh_rise_idx(velocity_thresh_rise_idx > ...
        find(zero_vel_down,1,'last') | velocity_thresh_rise_idx < ...
        find(zero_vel_up,1,'first')) = [];
    
    peak_frames = arrayfun(@(x) ...
        velocity_thresh_rise_idx(x) - 1 + ...
        find(zero_vel_down(velocity_thresh_rise_idx(x):end),1), ...
        1:length(velocity_thresh_rise_idx));
        
    % Find maximum raw df/f value within +/- 5 frames of smoothed peak
    peak_frame_search = arrayfun(@(x) ...
        peak_frames(x)-5:peak_frames(x)+5,1:length(peak_frames),'uni',false);
    peak_frame_search_cutoff = cellfun(@(x) x(x > 0 & x < ...
        size(df_matrix,2)),peak_frame_search,'uni',false);
    [peak_df peak_df_frame_offset] = ...
        cellfun(@(x) max(df_matrix(curr_cell,x)),peak_frame_search_cutoff);
    peak_df_frame = peak_frames+(peak_df_frame_offset-6);
    
    % Find event start time and fluorescence
    vel_start_frames = arrayfun(@(x) ...
        find(zero_vel_up(1:velocity_thresh_rise_idx(x)),1,'last'), ...
        1:length(velocity_thresh_rise_idx));
    low_df = df_matrix(curr_cell,vel_start_frames);
    
    % Need to define start from low_df (in case riding on another spike)
    % NOTE: 141008 AP - changed peak_frames to peak_df_frame here to make
    % sure that the start is before the fixed peak
    low_df_thresh = low_df + noise_std(curr_cell);
    event_start_frames = arrayfun(@(x) find(df_matrix(curr_cell,...
        (1:peak_df_frame(x))) < low_df_thresh(x),1,'last'), ...
        1:length(peak_df_frame),'uni',false);
    
    % Correct peak df values by fluorescence at start of event
    peak_df_corrected = peak_df - low_df;
    
    % Pare down false positive events
    %%%%
    eliminate_events = false(1,length(event_start_frames));
    
    % Eliminate any events which are below the noise threshold
    eliminate_events = eliminate_events | peak_df < noise_thresh(curr_cell);
    
    % Eliminiate events which aren't above noise for 2 frames
    eliminate_events = eliminate_events | df_matrix(curr_cell,peak_df_frame+1) - ...
        low_df < 2*noise_std(curr_cell) | ...
        df_matrix(curr_cell,peak_df_frame) < 3*noise_std(curr_cell);
    
    % Eliminate events which have a "baseline" of < 2 noise std
    eliminate_events = eliminate_events | low_df < -2*noise_std(curr_cell);
    
    % Eliminate events that have no start frame, i.e. at beginning of movie
    eliminate_events = eliminate_events | cellfun(@isempty,event_start_frames);
    
    % Cut bad events out
    event_start_frames(eliminate_events) = [];
    peak_df_corrected(eliminate_events) = [];
    peak_df_frame(eliminate_events) = [];
    
    event_start_frames = cell2mat(event_start_frames);
  
    %%%% Collect data to save
    
    % 1) event trace, 2) max velocity (stds) of event, 3) start/end frames
    
    caStds_vel{curr_cell} = nan(length(event_start_frames),1);
    caEventIdxs{curr_cell} = cell(length(event_start_frames),1);
    for curr_event = 1:length(event_start_frames);
        
        caEvents(curr_cell,event_start_frames(curr_event):...
            peak_df_frame(curr_event)) = peak_df_corrected(curr_event);
        
        % Get max velocity surrounding event (-5 start to +5 end)
        curr_search_frames = event_start_frames(curr_event)-5: ...
            peak_df_frame(curr_event)+5;
        curr_search_frames(curr_search_frames < 1 | ...
            curr_search_frames > size(df_matrix,2)) = [];
        
        caStds_vel{curr_cell}(curr_event) = ...
            max(ddf_matrix_smooth(curr_cell,curr_search_frames))/velocity_noise;
        
        caStds_df{curr_cell}(curr_event) = ...
            max(df_matrix_smooth(curr_cell,curr_search_frames))/noise_std(curr_cell);
       
        caEventIdxs{curr_cell}{curr_event} = ...
            event_start_frames(curr_event):peak_df_frame(curr_event);
        
    end    
    
    % Set current cell event trace NaN's to zeros
    caEvents(curr_cell,isnan(caEvents(curr_cell,:))) = 0;
  
end










