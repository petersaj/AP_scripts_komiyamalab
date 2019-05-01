% Corticospinal Analysis toolbox
%
% 
% This script contains general use analysis which is usually developed in
% dated analysis scripts
%
%

%% Group ROIs into likely same cell through correlation

d_grp_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    dendrite_corr = [];
    
    analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
    analysis_files = dir([analysis_path filesep '*.mat']);
    num_sessions = length(analysis_files);
    
    for curr_session = 1:num_sessions;
        
        curr_file = [analysis_path filesep analysis_files(curr_session).name];
        curr_data = load(curr_file);
        
        % Smooth data
        roi_trace_df_smooth = nan(size(curr_data.im.roi_trace_df));
        for curr_roi = 1:size(roi_trace_df_smooth,1);
            roi_trace_df_smooth(curr_roi,:) = smooth(curr_data.im.roi_trace_df(curr_roi,:),30,'loess');
        end
        
        curr_corr = corrcoef(roi_trace_df_smooth');
        dendrite_corr(:,:,curr_session) = curr_corr;
        
    end
    
    median_dendrite_corr = tril(median(dendrite_corr,3),-1);
    
    [d1,d2] = find(median_dendrite_corr >= 0.8);
    d_pairs = [d1 d2];
    
    % Group by multi-way combinations
    d_grp = {};
    all_d = unique([d1;d2]);
    for curr_d = all_d'
        
        old_grp = cellfun(@(x) ismember(curr_d,x),d_grp);
        
        if ~any(old_grp)
            curr_grp = length(d_grp) + 1;
            d_grp{curr_grp,1} = [];
        else
            % NOTE: this shouldn't yield more than 1 grp, but sometimes does
            % so do something about this later
            curr_grp = find(old_grp,1);
        end
        
        curr_pair = d_pairs == curr_d;
        paired_d = setdiff(reshape(d_pairs(any(curr_pair,2),:),[],1),curr_d);
        
        d_grp{curr_grp,1} = unique([d_grp{curr_grp};curr_d;paired_d]);
        
    end
    ungrouped = setdiff([1:size(dendrite_corr,1)]',unique(vertcat(d_grp{:})));
    d_grp = [d_grp;num2cell(ungrouped)];
    
    d_grp_all{curr_animal} = d_grp;
    
    disp(['Grouped animal ' animal]);
    
end

%% (Paired with above) Apply dendrite grouping to aligned analysis
% Requires analysis structure from AP_corticospinal_prepare_processed
% Requires above grouping script be run
for curr_animal = 1:length(analysis)
   for curr_session = 1:length(analysis(curr_animal).im)
       
       if isempty(analysis(curr_animal).im(curr_session).move_onset_aligned)
           continue
       end
       
        curr_grp_act = cell2mat(cellfun(@(x) ...
            nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,x),3), ...
            permute(d_grp_all{curr_animal},[2 3 1]),'uni',false));
    
        analysis(curr_animal).im(curr_session).move_onset_aligned = ...
            curr_grp_act;
        
   end
end


%% Get activity around any movements 

extracted_movements = struct('movements',cell(length(animals),1), ...
    'activity',cell(length(animals),1),'move_time',cell(length(animals),1));
    
for curr_animal = 1:length(animals);
    
    animal = animals{curr_animal};
    
       
    % Use this if data not loaded in
    data_dir = '/usr/local/lab/People/Andy/Data';
    curr_analysis_dir = [data_dir filesep animal filesep ...
        animal '_batch_thresh_roi'];
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    n_sessions = length(curr_analysis_files);
    for curr_session = 1:n_sessions
        curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
        curr_data = load(curr_file);
        
        % Threshold data
        event_thresh = 3;
        if exist('d_grp_all','var')
            % If ROIs grouped, then average together and threshold
            curr_df_grp = cell2mat(cellfun(@(x) nanmean(curr_data.im.roi_trace_df(x,:),1), ...
                d_grp_all{curr_animal},'uni',false));
            curr_df_thresh = AP_caEvents_thresh(curr_df_grp,event_thresh);
        else
            % If not, just threshold all
            curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
        end
        
        [lever_active,lever_force_resample] = ...
            AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 500;
        max_move_time = 3000;
        min_pre_move_iti = 1000;
        min_post_move_iti = 1000;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        criteria_movements = movement_starts(use_movements);
        
        [~,use_movement_start_frames] = ...
            arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
            criteria_movements(x))), 1:length(criteria_movements));
        
        
        frames_back = 0;
        frames_forward = round(max_move_time./curr_data.bhv.framerate);
        valid_frames = (use_movement_start_frames - frames_back > 0) & ...
            (use_movement_start_frames + frames_forward < size(curr_data.im.roi_trace_df,2));
        
        
        extracted_movements(curr_animal).activity{curr_session} = cell2mat(cellfun(@(x) curr_df_thresh(:, ...
            x-frames_back:x+frames_forward)',permute(num2cell( ...
            use_movement_start_frames(valid_frames)),[1 3 2]),'uni',false));
        
        extracted_movements(curr_animal).movements{curr_session} = cell2mat(cellfun(@(x) lever_force_resample(x:x+max_move_time), ...
            num2cell(criteria_movements(valid_frames))','uni',false));
        
        extracted_movements(curr_animal).move_time{curr_session} = movement_durations(use_movements);
        
        disp(curr_session)
    end
    
end
























