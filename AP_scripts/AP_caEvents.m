function caEvent_matrix = AP_caEvents(df_matrix,pyr_cells,gad_cells)
% caEvent_matrix = AP_caEvents(df_matrix,pyr_cells,gad_cells)
% Create a matrix of calcium events from raw df/f traces
% Pyr cells = peak detection, gad cells = start-end thresholding

% If pyr/gad not specified, assume all are pyramidal
if nargin == 1
    pyr_cells = true(size(df_matrix,1),1);
    gad_cells = ~pyr_cells;
end

% If logical, change to index
if islogical(pyr_cells)
    pyr_cells = find(pyr_cells); 
end
if islogical(gad_cells)
    gad_cells = find(gad_cells);
end

num_cells = size(df_matrix,1);

% Discard cells that contain all NaNs and fix occasional NaNs (is this ok?)
discard_cells = false(num_cells,1);
discard_cells(all(isnan(df_matrix),2)) = true;
df_matrix(isnan(df_matrix)) = 0;

pyr_use = setdiff(pyr_cells,find(discard_cells));
gad_use = setdiff(gad_cells,find(discard_cells));

% Initialize calcium events matrix
caEvent_matrix = nan(size(df_matrix));

% Find pyramidal events
if ~isempty(pyr_use)
    caEvent_matrix(pyr_use,:) = caEvents_pyr(df_matrix(pyr_use,:));
end

% Find gad events
if ~isempty(gad_use)
    caEvent_matrix(gad_use,:) = caEvents_gad(df_matrix(gad_use,:));
end

disp('Done.')

end








function pyr_caEvent_matrix = caEvents_pyr(pyr_matrix)
disp('Finding pyramidal calcium events...')
fprintf('%2d%%',0);

% Initialize pyramidal event matrix
pyr_caEvent_matrix = zeros(size(pyr_matrix));

num_cells = size(pyr_matrix,1);

% Smooth trace for easier detection
size_smooth = 30;
smooth_type = 'loess';
df_matrix_smooth = nan(size(pyr_matrix));
for curr_cell = 1:num_cells
    df_matrix_smooth(curr_cell,:) = smooth(pyr_matrix(curr_cell,:),size_smooth,smooth_type);
end

% Get first and second derivatives for detecting peaks
ddf_matrix_smooth = [zeros(num_cells,1) diff(df_matrix_smooth,[],2)];
d2df_matrix_smooth = [zeros(num_cells,1) diff(ddf_matrix_smooth,[],2)];

% Get noise estimates for each cell based on smoothing
noise_std = std(pyr_matrix - df_matrix_smooth,[],2);
noise_thresh = 3*noise_std;

% Loop through cells and find peaks
for curr_cell = 1:size(pyr_matrix,1);
    
    % Get velocity noise for current cell
    velocity_noise_thresh = 5*std(ddf_matrix_smooth(curr_cell, ...
        df_matrix_smooth(curr_cell,:) < noise_thresh(curr_cell)));
    
    % Find rising edges of where that threshold is crossed
    velocity_thresh_rise = ...
        [0 diff(sign(ddf_matrix_smooth(curr_cell,:) - velocity_noise_thresh))] == 2 & ...
        d2df_matrix_smooth(curr_cell,:) > 0;
    
    % Find the peak values associated with threshold crossings
    zero_vel_down = [0 diff(sign(ddf_matrix_smooth(curr_cell,:)))] == -2 & ...
        d2df_matrix_smooth(curr_cell,:) < 0;
    zero_vel_up = [0 diff(sign(ddf_matrix_smooth(curr_cell,:)))] == 2 & ...
        d2df_matrix_smooth(curr_cell,:) > 0;
    
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
        size(pyr_matrix,2)),peak_frame_search,'uni',false);
    [peak_df peak_df_frame_offset] = ...
        cellfun(@(x) max(pyr_matrix(curr_cell,x)),peak_frame_search_cutoff);
    peak_df_frame = peak_frames+(peak_df_frame_offset-6);
    
    % Find event start time and fluorescence
    vel_start_frames = arrayfun(@(x) ...
        find(zero_vel_up(1:velocity_thresh_rise_idx(x)),1,'last'), ...
        1:length(velocity_thresh_rise_idx));
    low_df = pyr_matrix(curr_cell,vel_start_frames);
    % Need to define start from low_df (in case riding on another spike)
    % NOTE: 141008 AP - changed peak_frames to peak_df_frame here to make
    % sure that the start is before the fixed peak
    low_df_thresh = low_df + noise_std(curr_cell);
    event_start_frames = arrayfun(@(x) find(pyr_matrix(curr_cell,...
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
    eliminate_events = eliminate_events | pyr_matrix(curr_cell,peak_df_frame+1) - ...
        low_df < 2*noise_std(curr_cell) | ...
        pyr_matrix(curr_cell,peak_df_frame) < 3*noise_std(curr_cell);
    
    % Eliminate events which have a "baseline" of < 2 noise std
    eliminate_events = eliminate_events | low_df < -2*noise_std(curr_cell);
    
    % Eliminate events that have no start frame, i.e. at beginning of movie
    eliminate_events = eliminate_events | cellfun(@isempty,event_start_frames);
    
    % Cut bad events out
    event_start_frames(eliminate_events) = [];
    peak_df_corrected(eliminate_events) = [];
    peak_df_frame(eliminate_events) = [];
    %%%%
    
    event_start_frames = cell2mat(event_start_frames);
    for record_event = 1:length(event_start_frames);
        pyr_caEvent_matrix(curr_cell,event_start_frames(record_event):...
            peak_df_frame(record_event)) = peak_df_corrected(record_event);
    end
    
    fprintf('%c%c%c%2d%%',8,8,8,round(100*curr_cell/num_cells));
end
fprintf('%c%c%c%c',8,8,8,8)

end




function gad_caEvent_matrix = caEvents_gad(gad_matrix)
disp('Finding gad calcium events...')
fprintf('%2d%%',0);

num_cells = size(gad_matrix,1);

% Initialize gad event matrix
gad_caEvent_matrix = zeros(size(gad_matrix));
for curr_cell = 1:size(gad_matrix,1)
    temp_smooth = smooth(gad_matrix(curr_cell,:),30,'loess')';
    noise_est = mean(abs(gad_matrix(curr_cell,:) - temp_smooth));
    
    thresh_1 = temp_smooth > noise_est*1;
    thresh_3 = temp_smooth > noise_est*3;
    
    % find edges of above-thresh periods
    thresh_1_start = diff([0 thresh_1 0]) == 1;
    thresh_3_start = find(diff([0 thresh_3 0]) == 1);
    thresh_3_stop = find(diff([0 thresh_3 0]) == -1);
    thresh_1_3_idx = arrayfun(@(x,y) x-find(thresh_1_start(x:-1:1),1)+1:y-1, ...
        thresh_3_start,thresh_3_stop,'uni',false);
    
    thresh_1_3 = false(1,size(gad_matrix,2));
    thresh_1_3(horzcat(thresh_1_3_idx{:})) = true;
    
    gad_caEvent_matrix(curr_cell,thresh_1_3) = gad_matrix(curr_cell,thresh_1_3);
    fprintf('%c%c%c%2d%%',8,8,8,round(100*curr_cell/num_cells));
end
fprintf('%c%c%c%c',8,8,8,8)
end







