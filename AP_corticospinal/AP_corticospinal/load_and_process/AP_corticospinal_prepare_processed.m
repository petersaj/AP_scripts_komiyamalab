function analysis = AP_corticospinal_prepare_processed(animals)
% analysis = AP_corticospinal_prepare_processed(animals,sessions)
%
% If 'sessions' isn't specified, all sessions loaded
%
% Prepare the processed data to be analyzed
% data from AP_corticospinal_prepare_processed

analysis = struct;

% % Pick analysis files
% all_analysis_path = cell(size(animals));
% all_analysis_files = cell(size(animals));
% for curr_animal = 1:length(animals)
%     [all_analysis_files{curr_animal} all_analysis_path{curr_animal}] = ...
%         uigetfile(['/usr/local/lab/People/Andy/Data' filesep animals{curr_animal}], ...
%         'multiselect','on',['Analysis files: ' animals{curr_animal}]);
% end

% Grab anaylsis files (assume dendrites in specific folder)
all_analysis_path = cell(size(animals));
all_analysis_files = cell(size(animals));
data_dir = '/usr/local/lab/People/Andy/Data';
for curr_animal = 1:length(animals)
    
    curr_analysis_dir = [data_dir filesep animals{curr_animal} filesep ...
        animals{curr_animal} '_batch_thresh_roi'];
    
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    
    all_analysis_path{curr_animal} = curr_analysis_dir;
    all_analysis_files{curr_animal} = sort({curr_analysis_files.name});
    
end

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    %% Load processed files
    
    %%%%% FIX THIS FOR LATER: MAKE SURE THAT SESSIONS CORRESPOND?
    
    %analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
    analysis_path = all_analysis_path{curr_animal};
    %analysis_files = dir([analysis_path filesep '*.mat']);
    analysis_files = sort(all_analysis_files{curr_animal});
    num_sessions = length(analysis_files);
           
%     if nargin == 1;
%         get_sessions = 1:num_sessions;
%     elseif nargin == 2;
%         get_sessions = sessions;
%     end
    
    for curr_session = 1:num_sessions
        
        %curr_file = [analysis_path filesep analysis_files(curr_session).name];
        curr_file = [analysis_path filesep analysis_files{curr_session}];
        curr_data = load(curr_file);
        
        % Skip if no behavior file (forgot to save)
        if ~isfield(curr_data.bhv,'bhv_times')
            continue
        end
        
        num_trials = length(curr_data.bhv.bhv_times);
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),curr_data.bhv.bhv_times);
        
        %%% Lever properties
        
        % Initialize
        analysis(curr_animal).lever(curr_session).reaction_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).move_to_reward_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement = cell(num_trials,1);
        
        % Parse movement - NOTE: this doesn't seem to be doing as good of a job
        % as it used to?
        [lever_active,lever_force_resample] = ...
            AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        curr_move_onset_frames = nan(num_trials,1);
        
        % Only go through trials which are rewarded and imaged
        for curr_trial = intersect(find(rewarded_trials)',curr_data.bhv.trial_list(:,2)');
            
            % Get the offset between xsg and bhv
            curr_list_idx = find(curr_data.bhv.trial_list(:,2) == curr_trial);
            curr_xsg_bhv_offset = curr_data.bhv.trial_list(curr_list_idx,1) - ...
                curr_data.bhv.bhv_times{curr_trial}.states.bitcode(1);
            
            % Get the reward in xsg time (in ms)
            if isfield(curr_data.bhv.bhv_times{curr_trial}.states,'reward_delay')
                curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward_delay(1) + ...
                    curr_xsg_bhv_offset)*1000);
            else
                curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward(1) + ...
                    curr_xsg_bhv_offset)*1000);
            end
            
            % Get the boundaries of the movement that overlaps with reward
            curr_rewarded_movement_onset = curr_reward_time - find( ...
                ~lever_active(curr_reward_time:-1:1),1) + 2;
            curr_rewarded_movement_offset = curr_reward_time + find( ...
                ~lever_active(curr_reward_time:end),1) - 2;
            
            [~,curr_move_onset_frames(curr_trial)] = ...
                min(abs(curr_data.bhv.frame_times*1000 - ...
                curr_rewarded_movement_onset));
            
            % Get reaction time (cue to movement onset, can be negative if
            % movement is animal moved early)
            curr_reaction_time = curr_rewarded_movement_onset - ...
                (curr_data.bhv.bhv_times{curr_trial}.states.cue(1) + ...
                curr_xsg_bhv_offset)*1000;
            
            % Get movement onset to reward time
            curr_move_to_reward_time = curr_reward_time - curr_rewarded_movement_onset;
            
            % Get total time of rewarded movement
            curr_rewarded_movement_time = curr_rewarded_movement_offset - ...
                curr_rewarded_movement_onset;
            
            % Get the rewarded movement from onset to 3 seconds after movement
            rewarded_movement_back = 1000;
            rewarded_movement_forward = 5000;
            curr_rewarded_movement = lever_force_resample(curr_rewarded_movement_onset-rewarded_movement_back: ...
                curr_rewarded_movement_onset+rewarded_movement_forward);
            
            
            % Package lever data into structure
            analysis(curr_animal).lever(curr_session).reaction_time(curr_trial) = curr_reaction_time;
            analysis(curr_animal).lever(curr_session).move_to_reward_time(curr_trial) = curr_move_to_reward_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement_time(curr_trial) = curr_rewarded_movement_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement{curr_trial} = curr_rewarded_movement;
            
        end
        
        % Define cued movement trials where reaction time is positive
        cued_movement_trials = analysis(curr_animal).lever(curr_session).reaction_time > 0;
        analysis(curr_animal).lever(curr_session).cued_movement_trials = cued_movement_trials;
        
        %%% Image alignment
        
        % Threshold data
        event_thresh = 3;
        curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
        
        % Get frames for cued rewarded movements
        use_move_onset_frames = curr_move_onset_frames(cued_movement_trials);
        use_reward_frames = cellfun(@(x) x.states.reward(1),curr_data.bhv.bhv_frames(cued_movement_trials));
        
        % Make matricies for pulling out frames
        surrounding_frames = -30:90;
        analysis(curr_animal).surrounding_frames = surrounding_frames;
        
        move_onset_frames_grab = repmat(use_move_onset_frames,1,length(surrounding_frames)) + ...
            repmat(surrounding_frames,length(use_move_onset_frames),1);
        move_onset_frames_grab(any(move_onset_frames_grab < 1,2),:) = [];
        
        reward_frames_grab = repmat(use_reward_frames,1,length(surrounding_frames)) + ...
            repmat(surrounding_frames,length(use_reward_frames),1);
        reward_frames_grab(any(reward_frames_grab < 1,2),:) = [];
        
        % Pull out frames
        num_rois = size(curr_data.im.roi_trace_df,1);
        
        curr_move_onset_aligned = permute(reshape(curr_df_thresh(:,move_onset_frames_grab'), ...
            num_rois,length(surrounding_frames),[]),[3,2,1]);
        
        curr_reward_aligned = permute(reshape(curr_df_thresh(:,reward_frames_grab'), ...
            num_rois,length(surrounding_frames),[]),[3,2,1]);
        
        % Package aligned images into structure
        analysis(curr_animal).im(curr_session).move_onset_aligned = curr_move_onset_aligned;
        analysis(curr_animal).im(curr_session).reward_aligned = curr_reward_aligned;
        
        % Display which file was processed
        disp(['Processed ' curr_file]);
        
    end
    
    disp('Finished processing');
    
    % Get ROI labels if available
    roilabels_file = dir([analysis_path filesep '*.roilabel']);
    if ~isempty(roilabels_file)
        curr_labels = load([analysis_path filesep roilabels_file.name],'-MAT');
        gad_cells = cellfun(@(x) any(strcmp('gad',x)),curr_labels.roi_labels)';
        pyr_cells = ~gad_cells;
        filled_cells = cellfun(@(x) any(strcmp('filled',x)),curr_labels.roi_labels)';
        
        analysis(curr_animal).labels.gad_cells = gad_cells;
        analysis(curr_animal).labels.filled_cells = filled_cells;
    end
    
end




