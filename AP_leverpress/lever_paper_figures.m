%% PUT IN MEMORY: all lever data (including Simon's animals)

%animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' ...
%    'AP78' 'AP79' 'SC027' 'SC028' 'SC029' ...
%    'AP91' 'AP92' 'AP93' 'AP94'};
animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' ...
    'AP78' 'AP79' 'SC027' 'SC028' 'SC029'};
lever_param_all = cell(size(animals));

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    filename_delimiter = cellfun(@(x) strfind(x,filesep),data_days,'uni',false);
    days = cellfun(@(x,y) x(y(end)+1:end),data_days,filename_delimiter, ...
        'UniformOutput',false);
    
    lever_param.pre_cue_movement= cell(size(days));
    lever_param.response_movement = cell(size(days));
    lever_param.reward_reaction = cell(size(days));
    lever_param.reaction_time = cell(size(days));
    lever_param.movement_reward_time = cell(size(days));
    lever_param.cued_movement = cell(size(days));
    lever_param.catch_trial = cell(size(days));
    lever_param.iti_movement = cell(size(days));
    lever_param.rewarded_lever_force = cell(size(days));
    lever_param.rewarded_lever_force_fixtime = cell(size(days));
    lever_param.move_reward_time = cell(size(days));
    lever_param.reward_move_ratio = cell(size(days));
    lever_param.reward_move_postcue_time = cell(size(days));
    lever_param.reward_move_precue_time = cell(size(days));
    lever_param.reward_move_postcue = cell(size(days));
    lever_param.iti_time = cell(size(days));
    lever_param.early_movement_ratio = cell(size(days));
    lever_param.rewarded_lever_force = cell(size(days));
    lever_param.rewarded_movement_distance = cell(size(days));
    lever_param.cue_quiescence = cell(size(days));
    lever_param.iti_quiescence = cell(size(days));
    lever_param.iti_quiescence_rewardedmove = cell(size(days));
    lever_param.iti_quiescence_middlemove = cell(size(days));
    lever_param.iti_quiescence_endmove = cell(size(days));
    lever_param.iti_quiescence_nonrewardedmove = cell(size(days));
    lever_param.reward_aligned_lever = cell(size(days));
    lever_param.reward_lever_peak = cell(size(days));
    
    
    for curr_day = 1:length(days)
        % Initialize
        clearvars -except animals curr_animal animal days ...
            curr_day lever_param data_path lever_param_all
        
        day = num2str(days{curr_day});
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep animal(1:2) repmat('0',1,4-length(animal(3:end))) animal(3:end)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % if no xsg (because accidentally not recorded), skip
        if isempty(xsg_fullfilenames)
            continue
        end
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        lever_param.pre_cue_movement{curr_day} = nan(num_trials,1);
        lever_param.response_movement{curr_day} = nan(num_trials,1);
        lever_param.reward_reaction{curr_day} = nan(num_trials,1);
        lever_param.reaction_time{curr_day} = nan(num_trials,1);
        lever_param.movement_reward_time{curr_day} = nan(num_trials,1);
        lever_param.cued_movement{curr_day} = false(num_trials,1);
        lever_param.catch_trial{curr_day} = false(num_trials,1);
        lever_param.iti_movement{curr_day} = nan(num_trials,1);
        lever_param.rewarded_lever_force{curr_day} = cell(num_trials,1);
        lever_param.rewarded_lever_force_fixtime{curr_day} = cell(num_trials,1);
        lever_param.move_reward_time{curr_day} = nan(num_trials,1);
        lever_param.reward_move_ratio{curr_day} = nan(num_trials,1);
        lever_param.reward_move_postcue_time{curr_day} = nan(num_trials,1);
        lever_param.reward_move_precue_time{curr_day} = nan(num_trials,1);
        lever_param.reward_move_postcue{curr_day} = nan(num_trials,1);
        lever_param.iti_time{curr_day} = nan(num_trials,1);
        lever_param.early_movement_ratio{curr_day} = nan(num_trials,1);
        lever_param.rewarded_movement_distance{curr_day} = nan(num_trials,1);
        lever_param.cue_quiescence{curr_day} = nan(num_trials,1);
        lever_param.iti_quiescence{curr_day} = nan(num_trials,1);
        lever_param.iti_quiescence_rewardedmove{curr_day} = nan(num_trials,1);
        lever_param.iti_quiescence_middlemove{curr_day} = nan(num_trials,1);
        lever_param.iti_quiescence_endmove{curr_day} = nan(num_trials,1);
        lever_param.iti_quiescence_nonrewardedmove{curr_day} = nan(num_trials,1);
        lever_param.reward_lever_peak{curr_day} = nan(num_trials,1);
        % Process
        
        % Only look at rewarded trials
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));

        % record catch trials, if they exist (remember: can include up to
        % two unrecorded trials, 1 for not finished, 2 for not started)
        if isfield(saved_history,'CatchSection_catch_trial')
            lever_param.catch_trial{curr_day}(1:length( ...
                saved_history.CatchSection_catch_trial)) = ...
                cell2mat(saved_history.CatchSection_catch_trial);
        end
        
        % Initializing saving raw traces around reward
        time_surround = 1;
        % NOTE! CHANGED TEMP FROM -0.5:4 to -1:5 (131202)
        %sample_surround = [-1*xsg_samplerate:5*xsg_samplerate];
        sample_surround = [-0.5*xsg_samplerate:4*xsg_samplerate];
        lever_param.reward_aligned_lever{curr_day} = ...
            nan(num_trials,length(sample_surround));
        use_samples = 1:length(sample_surround);

        
        % go through all xsg files, grab movement properties
        plots_on = 0;
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            
            % for Simon: sometimes XSG has no trials, if so, skip
            if isempty(curr_trial_list)
                continue
            end
            
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            % create movement velocity, active times
            % quantify movement before
            
            [lever_active lever_force_resample lever_force_smooth ...
                lever_velocity_envelope_smooth] = AP_parseLeverMovement(xsg_data);
            
            lever_stopstart = [0;diff(lever_active)];
            
            
            % get cue and reward times in resampled xsg
            cue_sample = nan(length(curr_trial_list(:,2)),1);
            reward_start_sample = nan(length(curr_trial_list(:,2)),1);
            reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
            iti_end_sample = nan(length(curr_trial_list(:,2)),1);
            for curr_trial_idx = 1:length(curr_trial_list(:,2));
                curr_trial = curr_trial_list(curr_trial_idx,2);
                             
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials) 
                    continue
                end

                % get trial start in dispatcher
                curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
                % get cue and reward sample in xsg (downsampled to ms)
                curr_bhv_cue_sample_rel =  (bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                cue_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_cue_sample_rel);
                
                curr_bhv_reward_start_sample_rel =  (bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_start_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_start_sample_rel);     
                
                curr_bhv_reward_stop_sample_rel =  (bhv{curr_trial}.states.reward(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_stop_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_stop_sample_rel); 
                
                curr_bhv_iti_end_sample_rel =  (bhv{curr_trial}.states.iti(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                iti_end_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_iti_end_sample_rel);
                                              
                
            end
        
            
            for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
                
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials) || ...
                    length(lever_force_resample) < cue_sample(curr_trial_idx+1)
                    continue
                end
               
                % ignore if out of bounds
                if cue_sample(curr_trial_idx+1) > length( ...
                        lever_force_resample)
                    continue
                end
                
                % ignore if reward doesn't come during classified movement
                % (usually because it's at the end, incomplete movement
                if ~lever_active(reward_start_sample(curr_trial_idx))
                    continue
                end
                
                % classify trial cue/movement overlap
                % rewarded movement has to start after cue, 100ms to cue
                % time has to be quiescent
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                cued_movement = rewarded_move_start > cue_sample(curr_trial_idx) && ...
                    ~any(lever_active(cue_sample(curr_trial_idx)-100: ...
                    cue_sample(curr_trial_idx)));
                lever_param.cued_movement{curr_day}(curr_trial) = cued_movement;
                
                % find movement parameters
                
                % get lever peak following movement
                % NOTE: using here approx threshold = 1.6
                prereward_thresh_cross = find(lever_force_resample(1: ...
                    reward_start_sample(curr_trial_idx)) < 1.6 & ...
                    diff([0;lever_force_smooth( ...
                    1:reward_start_sample(curr_trial_idx))]) > 0,1,'last');
                rewarded_peak = find(diff([0;lever_force_smooth( ...
                    prereward_thresh_cross:end)]) <= 0,1) + ...
                    prereward_thresh_cross;
                lever_param.reward_lever_peak{curr_day}(curr_trial) = ...
                    lever_force_resample(rewarded_peak);
                
                % get movement-onset aligned movemenet
                %curr_reward_sample = reward_start_sample(curr_trial_idx)*10;
                curr_reward_sample = rewarded_move_start*10;
                curr_samples = curr_reward_sample + sample_surround;
                valid_samples = curr_samples > 0 & curr_samples < ...
                    length(xsg_data.data.acquirer.trace_2);
                lever_param.reward_aligned_lever{curr_day} ...
                    (curr_trial,use_samples(valid_samples)) = ...
                    xsg_data.data.acquirer.trace_2(curr_samples(valid_samples));
                
                % get movement before cues
                pre_cue_sec = 1;
                pre_cue_samples = round(pre_cue_sec*(xsg_samplerate/10));
                if cue_sample(curr_trial_idx)-pre_cue_samples > 0
                    lever_param.pre_cue_movement{curr_day}(curr_trial) = ...
                        sum(lever_velocity_envelope_smooth( ...
                        cue_sample(curr_trial_idx)-pre_cue_samples: ...
                        cue_sample(curr_trial_idx)));
                end
                
                % get movement between cue and reward
                lever_param.response_movement{curr_day}(curr_trial) = ...
                    sum(lever_velocity_envelope_smooth( ...
                    cue_sample(curr_trial_idx):reward_stop_sample(curr_trial_idx))); %...
                    %/(reward_sample(curr_trial_idx)-cue_sample(curr_trial_idx));
                
                % get time between cue and reward
                lever_param.reward_reaction{curr_day}(curr_trial) = ...
                    reward_start_sample(curr_trial_idx) - cue_sample(curr_trial_idx);
                
                % get time from rewarded movement onset to reward
                lever_param.movement_reward_time{curr_day}(curr_trial) = ...
                    find(~lever_active(reward_start_sample(curr_trial_idx):-1:1),1);
                
                % get time to first movement after cue
                lever_param.reaction_time{curr_day}(curr_trial) = ...
                    find(lever_active(cue_sample(curr_trial_idx):...
                    reward_start_sample(curr_trial_idx)),1);
                if lever_param.reaction_time{curr_day}(curr_trial) == 1 && cued_movement
                    keyboard
                end
                lever_param.cue_quiescence{curr_day}(curr_trial) = ...
                    sum(lever_active(cue_sample(curr_trial_idx):...
                    reward_start_sample(curr_trial_idx)))/...
                    (reward_start_sample(curr_trial_idx)- ...
                    cue_sample(curr_trial_idx)+1);
                
                % break this down more: movement continued from reward,
                % movements in the middle, movements that don't end when
                % the trial ends
                lever_param.iti_quiescence{curr_day}(curr_trial) = ...
                    sum(lever_active(reward_stop_sample(curr_trial_idx):...
                    iti_end_sample(curr_trial_idx)))/...
                    (iti_end_sample(curr_trial_idx)- ...
                    reward_stop_sample(curr_trial_idx)+1);
                
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                rewarded_move_end = reward_start_sample(curr_trial_idx) + ...
                    find(lever_active(...
                    reward_start_sample(curr_trial_idx):end) == 0,1);
                % only used for breaking down iti
                endtrial_move_start = find(lever_active(1: ...
                    iti_end_sample(curr_trial_idx)) == 0,1,'last');
                
                if rewarded_move_end > iti_end_sample(curr_trial_idx)
                   rewarded_move_end = iti_end_sample(curr_trial_idx); 
                end
                
                if rewarded_move_end < reward_stop_sample(curr_trial_idx);
                    rewarded_move_end = reward_stop_sample(curr_trial_idx);
                end

                
                lever_param.iti_quiescence_rewardedmove{curr_day}(curr_trial) = ...
                    sum(lever_active(reward_stop_sample(curr_trial_idx):...
                    rewarded_move_end))/...
                    (iti_end_sample(curr_trial_idx)- ...
                    reward_stop_sample(curr_trial_idx)+1);

                lever_param.iti_quiescence_middlemove{curr_day}(curr_trial) = ...
                    sum(lever_active(rewarded_move_end:...
                    endtrial_move_start))/...
                    (iti_end_sample(curr_trial_idx)- ...
                    reward_stop_sample(curr_trial_idx)+1);
                
                lever_param.iti_quiescence_endmove{curr_day}(curr_trial) = ...
                    sum(lever_active(endtrial_move_start:...
                    iti_end_sample(curr_trial_idx)))/...
                    (iti_end_sample(curr_trial_idx)- ...
                    reward_stop_sample(curr_trial_idx)+1);
                
                lever_param.iti_quiescence_nonrewardedmove{curr_day}(curr_trial) = ...
                    sum(lever_active(rewarded_move_end:...
                    iti_end_sample(curr_trial_idx)))/...
                    (iti_end_sample(curr_trial_idx)- ...
                    rewarded_move_end+1);
                                             
                % get movement during the iti
                if curr_trial_idx ~= length(curr_trial_list(:,2)) && ...
                        cue_sample(curr_trial_idx+1) < length( ...
                        lever_velocity_envelope_smooth);
                    lever_param.iti_movement{curr_day}(curr_trial) = ...
                        sum(lever_velocity_envelope_smooth( ...
                        reward_stop_sample(curr_trial_idx): ...
                        cue_sample(curr_trial_idx+1))); % ...
                    % /(cue_sample(curr_trial_idx+1) - ...
                    % reward_sample(curr_trial_idx));
                    
                    % get iti time
                    lever_param.iti_time{curr_day}(curr_trial) = ...
                        cue_sample(curr_trial_idx+1) - ...
                        reward_stop_sample(curr_trial_idx);
                end
                    
                % extract the rewarded movement
                % get start of rewarded movement bout
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                rewarded_move_end = reward_start_sample(curr_trial_idx) + ...
                    find(lever_active(...
                    reward_start_sample(curr_trial_idx):end) == 0,1);
                % if the movement goes until the end of the trace, just
                % take what's available
                if isempty(rewarded_move_end)
                    rewarded_move_end = length(lever_force_resample);
                end
                lever_param.rewarded_lever_force{curr_day}{curr_trial} =  ...
                    lever_force_resample(rewarded_move_start: ...
                    rewarded_move_end);
                
                lever_param.rewarded_movement_distance{curr_day}(curr_trial) =  ...
                    sum(lever_velocity_envelope_smooth(rewarded_move_start: ...
                    rewarded_move_end));
                
                % get lever press for fixed time around reward
                if reward_start_sample(curr_trial_idx)-1000>0 && ...
                        reward_start_sample(curr_trial_idx)+300<length(lever_force_resample);
                    lever_param.rewarded_lever_force_fixtime{curr_day}{curr_trial} = ...
                        lever_force_resample(reward_start_sample(curr_trial_idx)-1000: ...
                        reward_start_sample(curr_trial_idx)+300);
                end
                
                lever_param.move_reward_time{curr_day}(curr_trial) = ...
                    reward_start_sample(curr_trial_idx) - rewarded_move_start;
                
                % get the ratio of post-cue rewarded movements/time to all
                % other movement/time
                if rewarded_move_start > cue_sample(curr_trial_idx)
                    temp_start = rewarded_move_start;
                else
                    temp_start = cue_sample(curr_trial_idx);                   
                end
                if curr_trial_idx < length(cue_sample)
                    if ~isnan(cue_sample(curr_trial_idx+1))
                        temp_reward_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (temp_start:...
                            rewarded_move_end));%/(temp_start-cue_sample(curr_trial_idx));
                        temp_all_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx):...
                            cue_sample(curr_trial_idx+1))); %/(cue_sample(curr_trial_idx+1)- ...
                        %cue_sample(curr_trial_idx));
                        lever_param.reward_move_ratio{curr_day}(curr_trial) = ...
                            temp_reward_move/temp_all_move;
                    end
                end
                
                
                % get movement ratio from fixed early time to all time
                % (median response time of all animals, 490ms, +400 ms water
                % time)
                
                pre_cue_time = 0;
                post_cue_time = 2000; % was previously 890
                if curr_trial_idx > 1 && curr_trial_idx < length(cue_sample)
                    if ~isnan(cue_sample(curr_trial_idx+1))
                        temp_early_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx)-pre_cue_time: ...
                            cue_sample(curr_trial_idx)+post_cue_time));
                        temp_all_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx): ...
                            cue_sample(curr_trial_idx+1)));
                        
                        lever_param.early_movement_ratio{curr_day}(curr_trial) = ...
                            temp_early_move/temp_all_move;
                    end
                end
                
                % If this is needed in the future: decide reward start/end
%                 % get amount of movement before cue/after reward
%                 lever_param.reward_move_postcue{curr_day}(curr_trial) = ...
%                     sum(lever_velocity_envelope_smooth( ...
%                     rewarded_move_end - reward_sample(curr_trial_idx)));
%                 
%                 lever_param.reward_move_postcue_time{curr_day}(curr_trial) = ...
%                     rewarded_move_end - reward_sample(curr_trial_idx);
%                 lever_param.reward_move_precue_time{curr_day}(curr_trial) = ...
%                     cue_sample(curr_trial_idx) - rewarded_move_start;
                          
                
                if plots_on
                    subplot(h1)
                    line([reward_sample_start(curr_trial_idx) reward_sample_start(curr_trial_idx)],ylim);
                    line([reward_sample_stop(curr_trial_idx) reward_sample_stop(curr_trial_idx)],ylim);
                    if cued_movement
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
                    else
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
                    end
                    subplot(h2)
                    line([reward_sample_start(curr_trial_idx) reward_sample_start(curr_trial_idx)],ylim);
                    line([reward_sample_stop(curr_trial_idx) reward_sample_stop(curr_trial_idx)],ylim);
                    if cued_movement
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
                    else
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
                    end
                end
            end
        end
        disp(['Finished day: ' day]);
    end
    disp('Finished all')
    
lever_param_all{curr_animal} = lever_param;

end
disp('Finished all animals')

%% PUT IN MEMORY: activity/classification, movements/dispatcher

%animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79', ...
%    'AP95' 'AP96' 'AP97' 'AP98' 'AP99' 'AP100'};
animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

gad_activity_thresh_all = cell(length(animals),15);
gad_class_all = cell(length(animals),15);

pyr_activity_all = cell(length(animals),15);
pyr_class_all = cell(length(animals),15);

pyr_class_rewardmove = cell(length(animals),15);
gad_class_rewardmove = cell(length(animals),15);

movement_trace_all = cell(length(animals),15);
move_epoch_frames = cell(length(animals),15);
rewarded_movements = cell(length(animals),15);
cued_rewarded_movements = cell(length(animals),15);
catch_rewarded_movements = cell(length(animals),15);

lick_times_all = cell(length(animals),15);
lick_trace_all = cell(length(animals),15);

cue_frames_all = cell(length(animals),15);
reward_frames_all = cell(length(animals),15);
cued_reward_frames_all = cell(length(animals),15);
catch_reward_frames_all = cell(length(animals),15);

trial_frames_all = cell(length(animals),15);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    % load ROIs and get ROI centroids and distances
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    roi_dir = dir(roi_path);
    roi_filenames = {roi_dir(:).name};
    roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
    roi_files_sort = sort(roi_filenames(roi_files));
    load([roi_path filesep roi_files_sort{end}],'-MAT');
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for curr_roi  = 1:num_rois
        [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
            polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
    end
    roi_dist = nan(length(roi_center_x));
    for k = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
                (roi_center_y(k) - roi_center_y(j))^2);
            roi_dist(k,j) = curr_dist;
        end
    end
    
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    % load movement classifications for cells
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        'final_class/' animal '_classified_cells_caEvents_10k'])
        
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);

    % Classification based on rewarded moves only
%            
%     clear classified_cells
%     load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
%         'final_class/' animal '_classified_cells_caEvents_rewardedmove_10k'])
%         
%     pyr_class_rewardmove(curr_animal,sessions) = ...
%         cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
%     gad_class_rewardmove(curr_animal,sessions) = ...
%         cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
%         
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);


    if strcmp(animal,'AP98') || strcmp(animal,'AP99') ||strcmp(animal,'AP100')
        days = days(8:end);
        sessions = 1:7;
    end
    
    for curr_day = sessions
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        %%%% Find and store activity
        pyr_activity_all{curr_animal,curr_day} = ...
            AP_caEvents(im.roi_trace_df(pyr_unfilled,:));
        
        gad_activity_thresh_all{curr_animal,curr_day} = ...
        AP_caEvents(im.roi_trace_df(gad_unfilled,:),[],1:sum(gad_unfilled));
       
        %%%% Find and store movements
        
        % get lick frames (in licks per frame);
        lick_trace_all{curr_animal,curr_day} = zeros(1,size(im.roi_trace_df,2));
        lick_trace_all{curr_animal,curr_day}(round( ...
            bhv.lick_frames(~isnan(bhv.lick_frames)))) = 1;
        lick_times_all{curr_animal,curr_day} = bhv.lick_frames(~isnan(bhv.lick_frames));
        
        num_frames = size(im.roi_trace_df,2);
        
        % load split image file to get loop sizes        
        %loop_frames = unique(bhv.add_frames);
        %loop_size = [diff(loop_frames);num_frames-loop_frames(end)];
        framerate = bhv.framerate;
        clear bhv;
        
        %if loading in image anyway, this is probably safer
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        % get behavior files/data
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));        
        
        catch_trials = [];
        if isfield(saved_history,'CatchSection_catch_trial')
            catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
            % there's one or more uncompleted trials at the end
            catch_trials = find(catch_trials(1:length(bhv)));
        end
                
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP' repmat('0',1,4-length(animal(3:end))) animal(3:end)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
             
        lever_active_full = cell(length(xsg_filenames),1);
        lever_velocity_full = cell(length(xsg_filenames),1);
        cue_frames_full =  cell(length(xsg_filenames),1);
        reward_frames_full =  cell(length(xsg_filenames),1);
        cued_reward_frames_full = cell(length(xsg_filenames),1);
        catch_reward_frames_full = cell(length(xsg_filenames),1);
        
        trial_frames = cellfun(@(x) nan(1,x),num2cell(loop_size),'uni',false);
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            [lever_active lever_force_smooth lever_force_velocity] = ...
                AP_parseLeverMovement(xsg_data);
            
            % find which movements are (cued) rewarded
            cue_sample = nan(length(curr_trial_list(:,2)),1);
            reward_start_sample = nan(length(curr_trial_list(:,2)),1);
            reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
            cued_movement = nan(length(curr_trial_list(:,2)),1);
            catch_movement = nan(length(curr_trial_list(:,2)),1);
            for curr_trial_idx = 1:length(curr_trial_list(:,2));
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials)
                    continue
                end
                
                % get trial start in dispatcher
                curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
                % get cue and reward sample in xsg (downsampled to ms)
                curr_bhv_cue_sample_rel =  (bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                cue_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_cue_sample_rel);
                
                curr_bhv_reward_start_sample_rel =  (bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_start_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_start_sample_rel);
                
                curr_bhv_reward_stop_sample_rel =  (bhv{curr_trial}.states.reward(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_stop_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_stop_sample_rel);
                
                if reward_start_sample(curr_trial_idx) > length(lever_active)
                    continue
                end
                
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                cued_movement(curr_trial_idx) = ...
                    rewarded_move_start > cue_sample(curr_trial_idx) && ...
                    ~any(lever_active(cue_sample(curr_trial_idx)-100: ...
                    cue_sample(curr_trial_idx)));
                                
                if ismember(curr_trial,catch_trials)
                    catch_movement(curr_trial_idx) = true;
                end
                
            end      
            cue_frames = round((cue_sample/1000)*framerate);
            reward_frames = round((reward_start_sample/1000)*framerate);
            cued_reward_frames = reward_frames(cued_movement == 1);
            catch_reward_frames = reward_frames(catch_movement == 1);
            
            cue_frames_full{curr_xsg} = cue_frames( ...
                cue_frames > 0 & ...
                cue_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);
            reward_frames_full{curr_xsg} = reward_frames( ...
                reward_frames > 0 & ...
                reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);            
            cued_reward_frames_full{curr_xsg} = cued_reward_frames( ...
                cued_reward_frames > 0 & ...
                cued_reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);            
            catch_reward_frames_full{curr_xsg} = catch_reward_frames( ...
                catch_reward_frames > 0 & ...
                catch_reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);
           
            %save the active lever portions in frames
            curr_lever_active_frames = round((find(lever_active)/1000) ...
                *framerate);
            lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
                curr_lever_active_frames <= loop_size(curr_xsg) & ...
                curr_lever_active_frames > 0)) + loop_frames(curr_xsg);
            
            [n d] = rat(loop_size(curr_xsg)/length(lever_force_velocity));
            lever_force_velocity_resample = resample(lever_force_velocity,n,d);
            lever_velocity_full{curr_xsg} = ...
                lever_force_velocity_resample(1:loop_size(curr_xsg));
            
            % save the trial numbers that each frame falls into
            % convert trial start times to frames
            curr_trial_list_frames = curr_trial_list;
            curr_trial_list_frames(:,1) = round( ...
                curr_trial_list_frames(:,1)*framerate);
            % save trials for each frame
            for curr_trial = 1:size(curr_trial_list_frames,1)-1;
                trial_frames{curr_xsg}(curr_trial_list_frames(curr_trial,1) : ...
                    curr_trial_list_frames(curr_trial+1,1)-1) = ...
                    curr_trial_list_frames(curr_trial,2);                
            end
            % fill in the end
            trial_frames{curr_xsg}(curr_trial_list_frames(end,1):end) = ...
                curr_trial_list_frames(end,2);            
            
        end
        trial_frames_all{curr_animal,curr_day} = horzcat(trial_frames{:});
        
        lever_velocity = vertcat(lever_velocity_full{:});
        
        lever_movement = vertcat(lever_active_full{:});
        movement_trace = zeros(1,num_frames);
        movement_trace(lever_movement) = 1;
        
        movement_trace_all{curr_animal,curr_day} = movement_trace;
        move_starts = find(diff([0 movement_trace 0]) == 1);
        move_stops =find(diff([0 movement_trace 0]) == -1);
        curr_move_epoch_frames = arrayfun(@(x) move_starts(x): ...
            move_stops(x),1:length(move_starts),'uni',false);
        
        move_epoch_frames{curr_animal,curr_day} = curr_move_epoch_frames;
        rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(reward_frames_full{:}))),curr_move_epoch_frames);
        cued_rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(cued_reward_frames_full{:}))),curr_move_epoch_frames);
        catch_rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(catch_reward_frames_full{:}))),curr_move_epoch_frames);
        
        cue_frames_all{curr_animal,curr_day} = vertcat(cue_frames_full{:});
        reward_frames_all{curr_animal,curr_day} = vertcat(reward_frames_full{:});
        cued_reward_frames_all{curr_animal,curr_day} = vertcat(cued_reward_frames_full{:});
        catch_reward_frames_all{curr_animal,curr_day} = vertcat(cued_reward_frames_full{:});
        
    end
  
    clc
    curr_animal
end


%% Figure 1: lever timing learning

% Cue to movement onset (cued movements only)
use_animals = 1:10;

curr_param = cell(7,15);
for i = use_animals
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reaction_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_param(i,1:length(curr_data)) = curr_data;  
end
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(curr_param{:,day_combine{i}});
    day_combine_data{i} = curr_data(:);
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
num_trials = cellfun(@length,curr_param);
curr_param_median(num_trials < 10) = NaN;
plot(curr_param_median(use_animals,1:14)','k');
%plot(nanmean(curr_param_median(use_animals,1:14)),'r','linewidth',2);
plot(cellfun(@nanmedian,day_combine_data),'r','linewidth',2)
%errorbar(nanmean(curr_param_median), ...
%    nanstd(curr_param_median)./sqrt(length(curr_param_median)),'r','linewidth',2);
% figure;errorbar(cellfun(@nanmean,day_combine_data), ...
%     cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);
xlabel('Session')
ylabel('Time (ms)')
title('Cue to movement onset')
% plot the median alone for this so that it can be blown up for inset
figure;
plot(cellfun(@nanmedian,day_combine_data),'r','linewidth',2)
xlabel('Session')
ylabel('Time (ms)')
title('Cue to movement onset (mean of medians only)')

% stats
p = anova1(curr_param_median,[],'off');


% Movement onset to reward
use_animals = 1:10;

curr_param = cell(7,15);
for i = use_animals
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reaction_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_data2 = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reward_reaction, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_data3 = cellfun(@(x,y) y-x,curr_data,curr_data2,'uni',false);
    curr_param(i,1:length(curr_data3)) = curr_data3;  
end
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(curr_param{:,day_combine{i}});
    day_combine_data{i} = curr_data(:);
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
num_trials = cellfun(@length,curr_param);
curr_param_median(num_trials < 10) = NaN;
plot(curr_param_median(use_animals,1:14)','k');
%plot(nanmean(curr_param_median(use_animals,1:14)),'r','linewidth',2);
plot(cellfun(@nanmedian,day_combine_data),'r','linewidth',2)
%errorbar(nanmean(curr_param_median), ...
%    nanstd(curr_param_median)./sqrt(length(curr_param_median)),'r','linewidth',2);
% figure;errorbar(cellfun(@nanmean,day_combine_data), ...
%     cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);
xlabel('Session')
ylabel('Time (ms)')
title('Movement onset to reward') 

% stats
p = anova1(curr_param_median,[],'off');

% Rewarded movement onset to reward (cued movements only)
curr_param = cell(7,15);
for i = use_animals
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.movement_reward_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_param(i,1:length(curr_data)) = curr_data;  
end
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(curr_param{:,day_combine{i}});
    day_combine_data{i} = curr_data(:);
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
num_trials = cellfun(@length,curr_param);
curr_param_median(num_trials < 10) = NaN;
plot(curr_param_median(use_animals,1:14)','k');
plot(nanmean(curr_param_median(use_animals,1:14)),'r','linewidth',2);
%errorbar(nanmean(curr_param_median), ...
%    nanstd(curr_param_median)./sqrt(length(curr_param_median)),'r','linewidth',2);
% figure;errorbar(cellfun(@nanmean,day_combine_data), ...
%     cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);
xlabel('Session')
ylabel('Time (ms)')
title('Rewarded movement onset to reward')

% stats
p = anova1(curr_param_median,[],'off');

% Fraction movement during ITI
curr_param = cell(7,15);
for i = use_animals
    curr_data = lever_param_all{i}.iti_quiescence_nonrewardedmove;

    curr_param(i,1:length(curr_data)) = curr_data;  
end
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(curr_param{:,day_combine{i}});
    day_combine_data{i} = curr_data(:);
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
num_trials = cellfun(@length,curr_param);
curr_param_median(num_trials < 10) = NaN;
plot(curr_param_median(use_animals,1:14)','k');
%plot(nanmean(curr_param_median(use_animals,1:14)),'r','linewidth',2);
plot(cellfun(@nanmedian,day_combine_data),'r','linewidth',2)
%errorbar(nanmean(curr_param_median), ...
%    nanstd(curr_param_median)./sqrt(length(curr_param_median)),'r','linewidth',2);
% figure;errorbar(cellfun(@nanmean,day_combine_data), ...
%     cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);
xlabel('Session')
ylabel('Fraction')
title('Fraction of time moving during ITI')

% stats
p = anova1(curr_param_median,[],'off');


% Fraction movement during ITI



%% Figure 1: lever average-average correlation and trial-trial correlation

lever_downsamp = 100;
use_samples = 49:350; % in 10 ms groups
lever_used = cell(11,15);
for curr_animal = 1:10
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    % extracting all lever
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        
        % use cued and uncued trials
        %curr_cued = true(size(curr_cued));
        
        temp_lever = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp_lever),2) & ...
            ~all(temp_lever == 0,2);
        temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
        total_curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);

    end
    lever_used(curr_animal,sessions) = total_curr_cat_lever(sessions);
    
    disp(curr_animal)
end

% plot trial-trial correlation
lever_trial_corr = cellfun(@(x) AP_itril(corrcoef(x'),-1),lever_used,'uni',false);
lever_trial_corr_median = cellfun(@median,lever_trial_corr);
lever_trial_corr_cat = cell(14,1);
use_animals = 1:10;
for i = 1:14
   %lever_trial_corr_cat{i} = vertcat(lever_trial_corr{use_animals,i}); 
   lever_trial_corr_cat{i} = vertcat(lever_trial_corr_median(use_animals,i)); 
end
lever_trial_corr_mean = cellfun(@nanmean,lever_trial_corr_cat); 
lever_trial_corr_sem = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), ...
    lever_trial_corr_cat);
figure;
errorbar(lever_trial_corr_mean,lever_trial_corr_sem,'k','linewidth',2)
xlabel('Session')
ylabel('Correlation')
title('Mean trial-trial correlation')

% stats
sessions = repmat(1:13,10,1);
lever_trial_corr_median_use = lever_trial_corr_median(1:10,1:13);
[r p] = corrcoef(lever_trial_corr_median_use(~isnan(lever_trial_corr_median_use)), ...
    sessions(~isnan(lever_trial_corr_median_use)))

% % plot avg-avg corrlation
% lever_avg = cellfun(@nanmean,lever_used,'uni',false);
% lever_avg_corr = nan(15,15,10);
% for curr_animal = 1:10
%    sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
%    curr_lever = vertcat(lever_avg{curr_animal,sessions});
%    lever_avg_corr(sessions,sessions,curr_animal) = corrcoef(curr_lever');  
% end
% use_animals = 1:10;
% figure;
% imagesc(nanmean(lever_avg_corr(1:14,1:14,use_animals),3))
% colormap(hot);caxis([0 1])
% xlabel('Session')
% ylabel('Session')
% title('Mean lever press correlation')
% 
% % plot avg-avg next-day correlation
% lever_avg_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_avg_corr(:,:,x),-1), ...
%     use_animals,'uni',false))';
% figure;errorbar(nanmean(lever_avg_corr_next_session),nanstd(...
%     lever_avg_corr_next_session)./sqrt(sum(~isnan(lever_avg_corr_next_session))), ...
%     'k','linewidth',2);
% 
% % plot avg trial-trial correlation across sessions
% lever_trial_corr = cell(15,15,10);
% for curr_animal = 1:10
%     sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
%     for curr_session_1 = sessions
%         for curr_session_2 = sessions           
%            if curr_session_1 == curr_session_2
%                curr_lever_cc = corrcoef(lever_used{curr_animal,curr_session_1}');
%                lever_trial_corr{curr_session_1,curr_session_2,curr_animal} = ...
%                    AP_itril(curr_lever_cc,-1);
%                continue
%            end
%            curr_lever_cc = corrcoef([lever_used{curr_animal,curr_session_1}; ...
%                lever_used{curr_animal,curr_session_2}]');
%            lever_trial_corr{curr_session_1,curr_session_2,curr_animal} = ...
%                curr_lever_cc(1:size(lever_used{curr_animal,curr_session_1},1), ...
%                size(lever_used{curr_animal,curr_session_1},1)+1:end);
%         end
%     end
% end
% use_animals = 1:7;
% lever_trial_corr_median = cellfun(@(x) nanmedian(x(:)),lever_trial_corr);
% figure;imagesc(nanmean(lever_trial_corr_median(1:14,1:14,use_animals),3));colormap(hot);
% title('Trial-trial correlation')
% 
% lever_trial_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_trial_corr_median(1:14,1:14,x),-1), ...
%     use_animals,'uni',false))';
% figure;errorbar(nanmean(lever_trial_corr_next_session),nanstd(...
%     lever_trial_corr_next_session)./sqrt(sum(~isnan(lever_trial_corr_next_session))), ...
%     'k','linewidth',2);
% title('Trial-trial next-session correlation')

% plot example lever traces
lever_used_mean = cellfun(@nanmean,lever_used,'uni',false);
% figure;
% for i = 1:11
%     for j = find(cellfun(@(x) ~isempty(x),lever_used(i,:)))
%         subplot(11,15,(i-1)*15+j)
%         plot(lever_used_mean{i,j},'k');
%         ylim([1 2])
%     end
% end

figure; hold on
plot(lever_used{6,1}(randperm(size(lever_used{6,1},1),10),:)','color',[0.5 0.5 0.5]);
plot(lever_used_mean{6,1},'k','linewidth',2)
title('AP76 Session 1')

figure; hold on
plot(lever_used{6,5}(randperm(size(lever_used{6,5},1),10),:)','color',[0.5 0.5 0.5]);
plot(lever_used_mean{6,5},'k','linewidth',2)
title('AP76 Session 5')

figure; hold on
plot(lever_used{6,14}(randperm(size(lever_used{6,14},1),10),:)','color',[0.5 0.5 0.5]);
plot(lever_used_mean{6,14},'k','linewidth',2);
title('AP76 Session 14')


%% Figure 1: example lever/pyr/gad traces


load('/usr/local/lab/People/Andy/Data/AP75/AP75_bgROI/120906_AP75_bgROI_analysis.mat' ...
    ,'im','bhv');
load('/usr/local/lab/People/Andy/Data/AP75/AP75_roilabels.roilabel','-MAT');

gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
pyr_cells = ~gad_cells;
filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';

pyr_unfilled = find(pyr_cells & ~filled_cells);
gad_unfilled = find(gad_cells & ~filled_cells);

use_frames = [20000:24000];

% plot bhv
curr_cue = intersect(round(bhv.cue_frames),use_frames)-use_frames(1)+1;
curr_reward = intersect(round(bhv.reward_frames),use_frames)-use_frames(1)+1;

cue_img = zeros(1,length(use_frames));
reward_img = zeros(1,length(use_frames));
move_img = zeros(1,length(use_frames));

for i = 1:length(curr_cue);
    cue_img(curr_cue(i):curr_reward(i)) = 1;
end

for i = 1:length(curr_reward);
    reward_img(curr_reward(i):curr_reward(i)+round(28*0.5)) = 1;
end

move_img(logical(movement_trace_all{5,12}(use_frames))) = 1;


figure; colormap(gray);
cue_plot = subplot(4,1,1);
imagesc(cue_img); title('cue')
reward_plot = subplot(4,1,2);
imagesc(reward_img); title('reward')
move_plot = subplot(4,1,3);
imagesc(move_img);title('move binary')
move_force_plot = subplot(4,1,4);
plot(-bhv.lever_force_resample(use_frames),'k');title('lever force')

linkaxes([cue_plot reward_plot move_plot move_force_plot],'x');

figure;
% plot only move-related, sort by timing
curr_timing = permute(mean(moveonset_aligned_pyr{5,12},1),[3 2 1]);
curr_move_idx = find(pyr_class_all{5,12});
[c ci] = max(curr_timing(curr_move_idx,:),[],2);
[asdf sort_idx] = sort(ci);
imagesc(im.roi_trace_df(pyr_unfilled(curr_move_idx(sort_idx)),use_frames));
caxis([0 0.5])
pyr_cmap = zeros(255,3);
pyr_cmap(:,2) = linspace(0,1,255);
pyr_cmap(128:end,1) = linspace(0,1,128);
pyr_cmap(128:end,3) = linspace(0,1,128);
colormap(pyr_cmap);
figure;plot(nanmean(im.roi_trace_df(pyr_unfilled,use_frames)),'color',[0 0.8,0]);
figure;plot(im.roi_trace_df(pyr_unfilled(69),use_frames),'color',[0 0.8,0]);


figure;
% plot only move-related, sort by timing
curr_timing = permute(mean(moveonset_aligned_gad{5,12},1),[3 2 1]);
curr_move_idx = find(gad_class_all{5,12});
[c ci] = max(curr_timing(curr_move_idx,:),[],2);
[asdf sort_idx] = sort(ci);
imagesc(im.roi_trace_df(gad_unfilled(curr_move_idx(sort_idx)),use_frames));
caxis([0 0.5])
gad_cmap = zeros(255,3);
gad_cmap(:,1) = linspace(0,1,255);
gad_cmap(128:end,2) = linspace(0,1,128);
gad_cmap(128:end,3) = linspace(0,1,128);
colormap(gad_cmap);
figure;plot(nanmean(im.roi_trace_df(gad_unfilled,use_frames)),'color',[0.8 0,0]);
figure;plot(im.roi_trace_df(gad_unfilled(14),use_frames),'color',[0.8 0,0]);

%% Figure 1: plot ROIs for example imaged field

load('/usr/local/lab/People/Andy/Data/AP71/AP71_bgROI/120619_AP71_bgROI.roi','-MAT');
load('/usr/local/lab/People/Andy/Data/AP71/AP71_roilabels.roilabel','-MAT')

gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels);
pyr_cells = ~gad_cells;
filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
pyr_unfilled = find(pyr_cells & ~filled_cells');
gad_unfilled = find(gad_cells & ~filled_cells');

pyr_idx = find(pyr_cells);
gad_idx = find(gad_cells);

figure;

for i = 1:length(pyr_idx)
    curr_cell = pyr_idx(i);
    curr_roi = polygon.ROI{curr_cell};
    
    smooth_filt = ones(10,1)/10;
    if size(curr_roi,1) > 20
        curr_roi = conv2(repmat(curr_roi,3,1),smooth_filt,'same');
        curr_roi = curr_roi(size(curr_roi,1)/3+1:end-size(curr_roi,1)/3,:);
    end
    
    patch(curr_roi(:,1),-curr_roi(:,2),[0 0.8 0]);
end

for i = 1:length(gad_idx)
    curr_cell = gad_idx(i);
    curr_roi = polygon.ROI{curr_cell};
    
    smooth_filt = ones(10,1)/10;
    if size(curr_roi,1) > 20
        curr_roi = conv2(repmat(curr_roi,3,1),smooth_filt,'same');
        curr_roi = curr_roi(size(curr_roi,1)/3+1:end-size(curr_roi,1)/3,:);
    end
    
    patch(curr_roi(:,1),-curr_roi(:,2),[0.8 0 0]);
end

ylim([-512*(um_px_y/um_px_x) 0]);
xlim([0 512]);

% make scalebar
% set pixels to microns ratio
um_px_x = 472/512;
um_px_y = 507/512;
figure;
line([0 0],[0 100/um_px_x],'linewidth',3,'color','k');
ylim([-512*(um_px_y/um_px_x) 0]);
xlim([0 512]);  

%% Figure 1: spine dynamics

% load spine data
load('/usr/local/lab/Projects/M1/LeverPress/Morphology/SC_allspinedata.mat')

figure; hold on;
SCmean_baseline = nanmean(SCbaseline_norm)/100;
SCsem_baseline = (nanstd(SCbaseline_norm)./sqrt(sum(~isnan(SCbaseline_norm))))/100;
errorbar(-6:0,SCmean_baseline,SCsem_baseline,'k','linewidth',2)
errorbar(1:14,SCmean, SCsem,'k','linewidth',2);
title('SC summary')
xlabel('Session');ylabel('Spines (fraction)');
xlim([-7 15]);ylim([0.95,1.15]);

figure; hold on
spine_add = cellfun(@(x) nanmean(x),Addition,'uni',false);
spine_elim = cellfun(@(x) nanmean(x),Elimination,'uni',false);
add_mean = nanmean(vertcat(spine_add{:}));
add_sem = nanstd(vertcat(spine_add{:}))/sqrt(3);
elim_mean = nanmean(-vertcat(spine_elim{:}));
elim_sem = nanstd(vertcat(spine_elim{:}))/sqrt(3);
add_bar = bar(add_mean);
elim_bar = bar(elim_mean);
set(add_bar,'FaceColor','k','linewidth',2);
set(elim_bar,'FaceColor',[0.5 0.5 0.5],'linewidth',2,'EdgeColor',[0.5 0.5 0.5]);
errorbar(add_mean,add_sem,'.k','linewidth',2);
errorbar(elim_mean,elim_sem,'.','color',[0.5 0.5 0.5],'linewidth',2);
title('Spine turnover')
xlabel('Session');ylabel('Spines')

% for lumping baseline sessions into one bar
add_mean_baseline = nanmean(baseline_Addition(:))/100;
elim_mean_baseline = nanmean(-baseline_Subtraction(:))/100;
add_sem_baseline = nanstd(baseline_Addition(:))/sqrt(sum(~isnan(baseline_Addition(:))))/100;
elim_sem_baseline = nanstd(baseline_Subtraction(:))/sqrt(sum(~isnan(baseline_Addition(:))))/100;
add_bar_baseline = bar(-1,add_mean_baseline);
elim_bar_baseline = bar(-1,elim_mean_baseline);
set(add_bar_baseline,'FaceColor','k','linewidth',2);
set(elim_bar_baseline,'FaceColor',[0.5 0.5 0.5],'linewidth',2,'EdgeColor',[0.5 0.5 0.5]);
errorbar(-1,add_mean_baseline,add_sem_baseline,'.k','linewidth',2);
errorbar(-1,elim_mean_baseline,elim_sem_baseline,'.','color',[0.5 0.5 0.5],'linewidth',2);

% for plotting all baseline bars
add_mean_baseline = nanmean(baseline_Addition)/100;
elim_mean_baseline = nanmean(-baseline_Subtraction)/100;
add_sem_baseline = nanstd(baseline_Addition)./sqrt(sum(~isnan(baseline_Addition)))/100;
elim_sem_baseline = nanstd(baseline_Subtraction)./sqrt(sum(~isnan(baseline_Addition)))/100;
add_bar_baseline = bar(-7:-1,add_mean_baseline);
elim_bar_baseline = bar(-7:-1,elim_mean_baseline);
set(add_bar_baseline,'FaceColor','k','linewidth',2);
set(elim_bar_baseline,'FaceColor',[0.5 0.5 0.5],'linewidth',2,'EdgeColor',[0.5 0.5 0.5]);
errorbar(-7:-1,add_mean_baseline,add_sem_baseline,'.k','linewidth',2);
errorbar(-7:-1,elim_mean_baseline,elim_sem_baseline,'.','color',[0.5 0.5 0.5],'linewidth',2);


% statistics
SCnorm = cell2mat(cellfun(@(x) nansum(x)/nansum(x(:,1)),SCcorr,'uni',false)');

early_spines = SCnorm(:,1:5);
early_sessions = ([1:5]'*ones(1,3))';
[cc_early p_early] = corrcoef(early_spines(:),early_sessions(:))
late_spines = SCnorm(:,5:14);
late_sessions = ([5:14]'*ones(1,3))';
[cc_late p_late] = corrcoef(late_spines(:),late_sessions(:))

%% Figure 2: movement pyr-gad balance

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity with lever 
pyr_move_active = cell(7,15);
gad_move_active = cell(7,15);
movement_lengths = cell(7,15);
act_trials_used = cell(7,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        %curr_moves = move_epoch_frames{curr_animal,curr_session}( ...
        %    rewarded_movements{curr_animal,curr_session});
        curr_moves = move_epoch_frames{curr_animal,curr_session};
      
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x) [leeway ...
            leeway],curr_moves,'uni',false);
        curr_moves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_moves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_moves_leeway);
        curr_moves_leeway(eliminate_movements) = [];
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_moves_leeway);
        
    end

    curr_animal
    
end

% plot pyr:gad ratio with all concatenated

use_animals = 1:7;
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
gc = horzcat(g{use_animals,:});
bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
figure;
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement based activity')




%% Figure 2: pyr-gad distance correlation

% set pixels to microns ratio
um_px_x = 472/512;
um_px_y = 507/512;

roi_distances = cell(7,1);
pyr_unfilled_idx = cell(7,1);
gad_unfilled_idx = cell(7,1);

animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for curr_animal = 1:7
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    pyr_unfilled_idx{curr_animal} = find(pyr_cells & ~filled_cells);
    gad_unfilled_idx{curr_animal} = find(gad_cells & ~filled_cells);
    
    % load ROIs and get ROI centroids and distances
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    roi_dir = dir(roi_path);
    roi_filenames = {roi_dir(:).name};
    roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
    roi_files_sort = sort(roi_filenames(roi_files));
    load([roi_path filesep roi_files_sort{end}],'-MAT');
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for curr_roi  = 1:num_rois
        [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
            polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
    end
    roi_center_x = roi_center_x*um_px_x;
    roi_center_y = roi_center_y*um_px_y;
    
    roi_dist = nan(length(roi_center_x));
    for k = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
                (roi_center_y(k) - roi_center_y(j))^2);
            roi_dist(k,j) = curr_dist;
        end
    end
    
    roi_distances{curr_animal} = roi_dist;
    
    disp(['Distances: animal ' animal]) 
    
end

% get pyr-gad correlations
pyr_gad_correlation_grid = cellfun(@(x,y) corrcoef([x' y']), ...
    pyr_activity_all,gad_activity_thresh_all,'uni',false);
pyr_gad_correlation = cellfun(@(x,y) x(1:size(y,1),size(y,1)+1:end), ...
    pyr_gad_correlation_grid,pyr_activity_all,'uni',false);

roi_distances_unfilled = cellfun(@(x,p,g) x(p,g), ...
    roi_distances,pyr_unfilled_idx,gad_unfilled_idx,'uni',false);

dist_bins = [0:50:300 Inf];
[n roi_dist_bins] = cellfun(@(x) histc(x,dist_bins), ...
    roi_distances_unfilled,'uni',false);

roi_dist_bin_rep = repmat(roi_dist_bins,1,15);

pyr_gad_correlation_vect = cellfun(@(x) x(:),pyr_gad_correlation,'uni',false);
roi_dist_bin_rep_vect = cellfun(@(x) x(:),roi_dist_bin_rep,'uni',false);

use_animals = 1:7;
animalday_used = cellfun(@(x) ~isempty(x),pyr_gad_correlation);
animalday_used(setdiff(1:7,use_animals),:) = 0;

pyr_gad_correlation_cat = vertcat(pyr_gad_correlation_vect{animalday_used});
roi_dist_bin_cat = vertcat(roi_dist_bin_rep_vect{animalday_used});

nonzero_bin = roi_dist_bin_cat ~= 0 & roi_dist_bin_cat ~= length(dist_bins);
pyr_gad_correlation_cat = pyr_gad_correlation_cat(nonzero_bin);
roi_dist_bin_cat = roi_dist_bin_cat(nonzero_bin);

corr_mean = grpstats(pyr_gad_correlation_cat,roi_dist_bin_cat,'mean');
corr_sem = grpstats(pyr_gad_correlation_cat,roi_dist_bin_cat,'sem');
figure;errorbar([dist_bins(2:end-1)-25 375],corr_mean,corr_sem,'k','linewidth',2);
xlabel('Distance (\mum)')
ylabel('Correlation');
title('Pyramidal-gad correlation by distance')

% stats
under_150_bins = roi_dist_bin_cat <= 3;
over_150_bins = roi_dist_bin_cat > 3;

under_150_corr = pyr_gad_correlation_cat(under_150_bins);
over_150_corr = pyr_gad_correlation_cat(over_150_bins);
p = ranksum(under_150_corr(~isnan(under_150_corr)), ...
    over_150_corr(~isnan(over_150_corr)));

%% Figure 2: gad/pyr dynamic populations

% plot class grids
use_animals = 1:7;

figure;
persistence = cell(8,2);
for curr_animal = use_animals   
    
    curr_gad_class = horzcat(gad_class_all{curr_animal,:});
    day_center = nanmean(curr_gad_class.*meshgrid(1:size(curr_gad_class,2), ...
        1:size(curr_gad_class,1)),2);
    gad_class_cells = find(any(curr_gad_class,2));
    [asdf sort_idx] = sort(day_center(gad_class_cells));
    
    subplot(8,2,curr_animal*2-1);
    imagesc(curr_gad_class(gad_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Gad')
    end
    
    curr_pyr_class = horzcat(pyr_class_all{curr_animal,:});
    day_center = nanmean(curr_pyr_class.*meshgrid(1:size(curr_pyr_class,2), ...
        1:size(curr_pyr_class,1)),2);
    pyr_class_cells = find(any(curr_pyr_class,2));
    [asdf sort_idx] = sort(day_center(pyr_class_cells));
    
    subplot(8,2,curr_animal*2);
    imagesc(curr_pyr_class(pyr_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Pyr')
    end
    
    persistence{curr_animal,1} = nanmean(curr_gad_class(gad_class_cells,:),2);
    persistence{curr_animal,2} = nanmean(curr_pyr_class(pyr_class_cells,:),2);
    
end

avg_persistence = cellfun(@nanmean,persistence);
figure; hold on;
bar(nanmean(avg_persistence),'FaceColor','w','linewidth',2,'BarWidth',0.5);
errorbar(nanmean(avg_persistence), ...
    nanstd(avg_persistence)./sqrt(sum(~isnan(avg_persistence))),'.k','linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Gad' 'Pyr'})
ylabel('Persistence')
xlim([0 3])

persistence_cat = cell(1,2);
persistence_cat{1} = vertcat(persistence{use_animals,1});
persistence_cat{2} = vertcat(persistence{use_animals,2});
figure;hold on
plot(sort(persistence_cat{2}),linspace(0,1, ...
    length(persistence_cat{2})),'color',[0 0.8 0],'linewidth',2)
plot(sort(persistence_cat{1}),linspace(0,1, ...
    length(persistence_cat{1})),'color',[0.8 0 0],'linewidth',2)
ylabel('Cumulative fraction of cells')
xlabel('Fraction of days active')
title('Gad cells are more persistent')
legend({'Pyr' 'Gad'},'location','se')
xlim([-0.1 1.1]);

% statistics
[h p] = kstest2(persistence_cat{1},persistence_cat{2});


%% Figure 3: fraction of movement-related cells

use_animals = 1:7;
pyr_class_mean_use = cellfun(@nanmean,pyr_class_all(use_animals,1:14));
pyr_class_mean_use_norm = bsxfun(@times,pyr_class_mean_use,1./max(pyr_class_mean_use,[],2));
figure; hold on;
errorbar(nanmean(pyr_class_mean_use),nanstd(pyr_class_mean_use)./ ...
    sqrt(sum(~isnan(pyr_class_mean_use))),'k','linewidth',2);
xlabel('Session')
ylabel('Fraction of cells')
title('Movement-related cells')
xlim([0 15]);

% statistics
[p_anova anova_stats] = anovan(pyr_class_mean_use_norm(~isnan(pyr_class_mean_use_norm)), ...
{g1(~isnan(pyr_class_mean_use_norm))},'display','off');

early_frac = pyr_class_mean_use_norm(:,1:3);
early_sessions = ([1:3]'*ones(1,7))';
early_nan_frac = isnan(early_frac);
[cc_early p_early] = corrcoef(early_frac(~early_nan_frac),early_sessions(~early_nan_frac))
late_frac = pyr_class_mean_use_norm(:,4:14);
late_sessions = ([3:14]'*ones(1,7))';
late_nan_frac = isnan(late_frac);
[cc_late p_late] = corrcoef(late_frac(~late_nan_frac),late_sessions(~late_nan_frac))

% use_animals = 1:7;
% class_cat = cell(14,1);
% for i = 1:14
%    class_cat{i} = vertcat(pyr_class_all{use_animals,i});
% end
% figure; hold on;
% plot(cellfun(@nanmean,class_cat),'k','linewidth',2);
% xlabel('Session')
% ylabel('Fraction of cells')
% title('Movement-related cells')
% xlim([0 15]);
% ylim([0.08 0.2]);

% use_animals = 1:7;
% class_cat = cell(14,1);
% for i = 1:14
%    class_cat{i} = vertcat(gad_class_all{use_animals,i});
% end
% figure; hold on;
% plot(cellfun(@nanmean,class_cat),'k','linewidth',2);
% xlabel('Session')
% ylabel('Fraction of cells')
% title('Movement-related cells')

%% Figure 3: examples of movement-related ROIs

% set pixels to microns ratio
um_px_x = 472/512;
um_px_y = 507/512;

animal = 'AP79';

load(['/usr/local/lab/People/Andy/Data/' ...
    animal '/' animal '_roilabels.roilabel'],'-MAT');

% load ROIs and get ROI centroids and distances
roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
roi_dir = dir(roi_path);
roi_filenames = {roi_dir(:).name};
roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
roi_files_sort = sort(roi_filenames(roi_files));
load([roi_path filesep roi_files_sort{end}],'-MAT');

gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels);
pyr_cells = ~gad_cells;

filled = cellfun(@(x) any(strcmp('filled',x)),roi_labels);

pyr_unfilled_idx = find(pyr_cells & ~filled);

% plot movement-related cell rois
figure;
plot_days = [1 3 10];

for curr_day = 1:length(plot_days)
    subplot(1,length(plot_days),curr_day);
    
    cell_color = ones(length(pyr_unfilled_idx),3)*0.8;
    cell_color(pyr_class_all{8,plot_days(curr_day)},1) = 0;
    cell_color(pyr_class_all{8,plot_days(curr_day)},2) = 0.5;
    cell_color(pyr_class_all{8,plot_days(curr_day)},3) = 0;
    
    roi_handle = [];
    for i = 1:length(pyr_unfilled_idx)
        curr_cell = pyr_unfilled_idx(i);
        roi_handle(i) = patch(polygon.ROI{curr_cell}(:,1),-polygon.ROI{curr_cell}(:,2), ...
            cell_color(i,:),'EdgeColor','none');
    end
    ylim([-512*(um_px_y/um_px_x) 0]);
    xlim([0 512]);   
    
    title(plot_days(curr_day));
end

% plot scalebar

figure;
subplot(1,length(plot_days),1);
line([0 0],[0 100/um_px_y],'linewidth',3,'color','k');
ylim([-512*(um_px_y/um_px_x) 0]);
xlim([0 512]);   
title('100 um')


%% Figure 3: pyr move cell population correlation

pyr_pop_corr = nan(15,15,7);
for i = 1:7
    sessions = find(cellfun(@(x) ~isempty(x),pyr_class_all(i,:)));
    curr_pop = horzcat(pyr_class_all{i,:});
    pyr_pop_corr(sessions,sessions,i) = corrcoef(curr_pop);
end
use_animals = 1:7;
figure;imagesc(nanmean(pyr_pop_corr(1:14,1:14,use_animals),3));
caxis([0 1]); colormap(hot);
xlabel('Session')
ylabel('Session')
title('Pyramidal move cell population correlation')
colorbar

% plot the second diagonal
pyr_pop_corr_diag = cell2mat(arrayfun(@(x) diag(pyr_pop_corr(:,:,x),-1), ...
    1:size(pyr_pop_corr,3),'uni',false))';
figure;errorbar(nanmean(pyr_pop_corr_diag(use_animals,1:13)), ...
    nanstd(pyr_pop_corr_diag(use_animals,1:13))./sqrt(sum(~isnan( ...
    pyr_pop_corr_diag(use_animals,1:13)))),'k','linewidth',2);
xlabel('Session')
ylabel('Correlation')
title('Pyramidal move cell population correlation, i, i+1')



%%%% try doing this but with the mean of rewarded movement activity (FINAL)

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity with lever 
pyr_move_active = cell(7,15);
gad_move_active = cell(7,15);
movement_lengths = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        curr_moves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        %curr_moves = move_epoch_frames{curr_animal,curr_session};
      
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x) [leeway ...
            leeway],curr_moves,'uni',false);
        curr_moves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_moves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_moves_leeway);
        curr_moves_leeway(eliminate_movements) = [];
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_moves_leeway);
        
    end

    curr_animal
    
end

pyr_move_active_mean = cellfun(@(x) nanmean(x>0,2),pyr_move_active,'uni',false);
pyr_pop_corr = nan(15,15,7);
for i = 1:7
    sessions = find(cellfun(@(x) ~isempty(x),pyr_move_active_mean(i,:)));
    curr_pop = horzcat(pyr_move_active_mean{i,:});
    pyr_pop_corr(sessions,sessions,i) = corrcoef(curr_pop);
end
use_animals = 1:7;
figure;imagesc(nanmean(pyr_pop_corr(1:14,1:14,use_animals),3));
caxis([0 1]); colormap(hot);
xlabel('Session')
ylabel('Session')
title('Pyramidal move cell population correlation (reliability)')
colorbar

% plot the second diagonal
pyr_pop_corr_diag = cell2mat(arrayfun(@(x) diag(pyr_pop_corr(:,:,x),-1), ...
    1:size(pyr_pop_corr,3),'uni',false))';
figure;errorbar(nanmean(pyr_pop_corr_diag(use_animals,1:13)), ...
    nanstd(pyr_pop_corr_diag(use_animals,1:13))./sqrt(sum(~isnan( ...
    pyr_pop_corr_diag(use_animals,1:13)))),'k','linewidth',2);
xlabel('Session')
ylabel('Correlation')
title('Pyramidal move cell population correlation, i, i+1')

% stats
early_corr = arrayfun(@(x) AP_itril(pyr_pop_corr(1:4,1:4,x),-1),1:7,'uni',false);
early_corr_cat = vertcat(early_corr{:});
late_corr = arrayfun(@(x) AP_itril(pyr_pop_corr(10:14,10:14,x),-1),1:7,'uni',false);
late_corr_cat = vertcat(late_corr{:});
p = ranksum(early_corr_cat(~isnan(early_corr_cat)),late_corr_cat(~isnan(late_corr_cat)));



%% Figure 3: bar plot of population composition

% plot class grids
pyr_oldfrac = nan(7,15);
pyr_firstfrac = nan(7,15);
pyr_newsessions = nan(7,15);
pyr_oldsessions = nan(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    curr_pyr_class = horzcat(pyr_class_all{curr_animal,:});
    curr_pyr_old = ~diff(curr_pyr_class,[],2) & curr_pyr_class(:,2:end);
    curr_pyr_new = curr_pyr_class(:,2:end) & ~curr_pyr_old;
    
    pyr_oldfrac(curr_animal,sessions(2:end)) = ...
        sum(curr_pyr_old)./sum(curr_pyr_class(:,2:end));
    
    curr_pyr_first = zeros(size(curr_pyr_class));
    [c ci] = max(curr_pyr_class,[],2);
    for i = 1:length(c)
        curr_pyr_first(i,ci(i)) = c(i);
    end
    
    pyr_firstfrac(curr_animal,sessions(2:end)) = ...
        sum(curr_pyr_first(:,2:end))./sum(curr_pyr_class(:,2:end));
    
    pyr_numsessions = nanmean(curr_pyr_class,2);
    pyr_newsessions(curr_animal,sessions(2:end)) = ...
        nansum(repmat(pyr_numsessions,1,size(curr_pyr_new,2)).* ...
        curr_pyr_new)./nansum(curr_pyr_new);
    pyr_oldsessions(curr_animal,sessions(2:end)) = ...
        nansum(repmat(pyr_numsessions,1,size(curr_pyr_old,2)).* ...
        curr_pyr_old)./nansum(curr_pyr_old);
    
end
use_animals = 1:7;

mean_oldfrac = nanmean(pyr_oldfrac(use_animals,2:14));
std_oldfrac = nanstd(pyr_oldfrac(use_animals,2:14))./ ...
    sqrt(sum(~isnan(pyr_oldfrac(use_animals,2:14))));
mean_firstfrac = nanmean(pyr_firstfrac(use_animals,2:14));
std_firstfrac = nanstd(pyr_firstfrac(use_animals,2:14))./ ...
    sqrt(sum(~isnan(pyr_firstfrac(use_animals,2:14))));

figure; hold on;
bar(2:14,[mean_oldfrac' mean_firstfrac' 1-(mean_oldfrac'+mean_firstfrac')], ...
    'stack','linewidth',2)
errorbar(2:14,nanmean(pyr_oldfrac(use_animals,2:14)), ...
    std_oldfrac,'.k','linewidth',2)

ylabel('Fraction of movement cells')
xlabel('Session')
legend({'Old' 'New'})
xlim([1 15]);
ylim([-0.1 1.1])



figure; hold on;
errorbar(nanmean(pyr_newsessions),nanstd(pyr_newsessions)./sqrt(sum( ...
    ~isnan(pyr_newsessions))),'k','linewidth',2);
errorbar(nanmean(pyr_oldsessions),nanstd(pyr_oldsessions)./sqrt(sum( ...
    ~isnan(pyr_oldsessions))),'r','linewidth',2);


% stats
early_oldfrac = reshape(pyr_oldfrac(:,1:3),[],1);
late_oldfrac = reshape(pyr_oldfrac(:,10:14),[],1);
p = ranksum(early_oldfrac(~isnan(early_oldfrac)), ...
    late_oldfrac(~isnan(late_oldfrac)));



% % get the random first appearance
% num_rep = 20;
% temp = nan(num_rep,13);
% for curr_rep = 1:num_rep
%     pyr_oldfrac = nan(7,15);
%     pyr_firstfrac = nan(7,15);
%     for curr_animal = 1:8
%         sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
%         
%         curr_pyr_class = shake(horzcat(pyr_class_all{curr_animal,:}),2);
%         curr_pyr_old = ~diff(curr_pyr_class,[],2) & curr_pyr_class(:,2:end);
%         
%         pyr_oldfrac(curr_animal,sessions(2:end)) = ...
%             sum(curr_pyr_old)./sum(curr_pyr_class(:,2:end));
%         
%         curr_pyr_first = zeros(size(curr_pyr_class));
%         [c ci] = max(curr_pyr_class,[],2);
%         for i = 1:length(c)
%             curr_pyr_first(i,ci(i)) = c(i);
%         end
%         
%         pyr_firstfrac(curr_animal,sessions(2:end)) = ...
%             sum(curr_pyr_first(:,2:end))./sum(curr_pyr_class(:,2:end));
%         
%     end
%     use_animals = 1:7;
%     
%     temp(curr_rep,:) = nanmean(pyr_firstfrac(use_animals,2:14));
%     curr_rep
% end

%% Figure 3: activity timing standard deviation / correlation

% get rewarded-movement onset activity timing
animaldays = find(cellfun(@(x) ~isempty(x),pyr_class_all))';
movement_timing_std = cell(7,15);
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    curr_pyrclass = pyr_class_all{curr_animalday};

    %curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_pyr_activity = pyr_activity_onset(curr_pyrclass,:);
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % find rewarded move onset aligned timings
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr_cat = cat(3,curr_moveonset_aligned_pyr{:});

    [c ci] = max(moveonset_aligned_pyr_cat>0,[],2);    
    
    use_cells = sum(c,3) > 5;
    
    ci_use = ci(use_cells,:,:);
    ci_use(c(use_cells,:,:) == 0) = NaN;
    movement_timing_std{curr_animalday} = nanstd(ci_use,[],3);
        
    curr_animalday/max(animaldays)

end

movement_timing_std_cat = cell(14,1);
use_animals = 1:7;
for i = 1:14
   movement_timing_std_cat{i} = vertcat(movement_timing_std{use_animals,i})/28;
end

m = cellfun(@nanmean,movement_timing_std_cat);
s = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),movement_timing_std_cat);
figure;errorbar(m,s,'k','linewidth',2);
ylabel('Standard deviation (s)')
xlabel('Session')
title('Standard deviation of pyr move cell activity timing')

% statistics
session_x = cellfun(@(x,y) x*ones(size(y)),num2cell(1:14)',movement_timing_std_cat,'uni',false);
session_x_cat = vertcat(session_x{:});
movement_timing_std_cat_cat = vertcat(movement_timing_std_cat{:});
[r p] = corrcoef(session_x_cat,movement_timing_std_cat_cat);


%%%%%%%%%%% the rest currently no in paper
% get change in std from cells which are present over sessions

animaldays = find(cellfun(@(x) ~isempty(x),pyr_class_all))';
movement_timing_std = cell(7,15);
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
    
     % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    curr_pyrclass = pyr_class_all{curr_animalday};

    %curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_pyr_activity = pyr_activity_onset(:,:);
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % find rewarded move onset aligned timings
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr_cat = cat(3,curr_moveonset_aligned_pyr{:});

    [c ci] = max(moveonset_aligned_pyr_cat>0,[],2);    
    
    use_cells = 1:size(c,1);
    rarely_active_cells = sum(c,3) > 10;
    
    ci_use = ci(use_cells,:,:);
    ci_use(c(use_cells,:,:) == 0) = NaN;
    movement_timing_std{curr_animalday} = nanstd(ci_use,[],3);
    movement_timing_std{curr_animalday}(rarely_active_cells) = NaN;
    
    curr_animalday/max(animaldays)

end

std_change = cell(7,1);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_class = horzcat(pyr_class_all{curr_animal,:});
    use_cells = find(sum(curr_class,2) > 3);
    std_change{curr_animal} = nan(length(use_cells),2);
    for curr_cell = 1:length(use_cells)
        first_session = sessions(find(curr_class(use_cells(curr_cell),:),1));
        last_session = sessions(find(curr_class(use_cells(curr_cell),:),1,'last'));
        
        std_change{curr_animal}(curr_cell,1) = ...
            movement_timing_std{curr_animal,first_session}(use_cells(curr_cell));
        std_change{curr_animal}(curr_cell,2) = ...
            movement_timing_std{curr_animal,last_session}(use_cells(curr_cell));
    end
end
use_animals = 1:7;
std_median = cell2mat(cellfun(@nanmedian,std_change(use_animals),'uni',false));
figure;plot(std_median','color',[0.5 0.5 0.5])
hold on
plot(nanmean(std_median),'k','linewidth',2)
xlim([0 3]);
ylabel('Standard deviation of activity timing')
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'First session classified' 'Last session classified'})

% CORRELATION 


% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

%%%%%%%% Trying out correlation types


%%% cell-trial mean-subtracted correlation
moveonset_aligned_pyr_meansub_trial = cellfun(@(x) reshape(permute( ...
    bsxfun(@minus,x,nanmean(x,2)),[2 3 1]),[],size(x,1)), ...
    moveonset_aligned_pyr,'uni',false);

moveonset_aligned_pyr_meansub_trial_nonan = ...
    moveonset_aligned_pyr_meansub_trial;
for i = find(cellfun(@(x) ~isempty(x),moveonset_aligned_pyr_meansub_trial))'
    moveonset_aligned_pyr_meansub_trial_nonan{i}(isnan( ...
        moveonset_aligned_pyr_meansub_trial_nonan{i})) = 0;
end

meansub_trial_corr = cellfun(@(x) AP_itril(corrcoef(x),-1), ...
    moveonset_aligned_pyr_meansub_trial_nonan,'uni',false);

use_animals = 1:7;

meansub_trial_corr_mean = cellfun(@nanmean,meansub_trial_corr(use_animals,1:14));
meansub_trial_corr_m = nanmean(meansub_trial_corr_mean);
meansub_trial_corr_s = nanstd(meansub_trial_corr_mean)./ ...
    sqrt(sum(~isnan(meansub_trial_corr_mean)));
figure;errorbar(meansub_trial_corr_m,meansub_trial_corr_s,'k','linewidth',2)

meansub_trial_corr_cat = cell(14,1);
for i = 1:14
   meansub_trial_corr_cat{i} = vertcat(meansub_trial_corr{use_animals,i});
end
meansub_trial_corr_m = cellfun(@nanmean,meansub_trial_corr_cat);
meansub_trial_corr_s = cellfun(@(x) nanstd(x)./ ...
    sqrt(sum(~isnan(x))),meansub_trial_corr_cat);
figure;errorbar(meansub_trial_corr_m,meansub_trial_corr_s,'k','linewidth',2)

%%%%%%%%

%%% get trial-average correlation

mean_aligned_activity = cellfun(@(x) reshape(permute(nanmean(x),[2 3 1]),[],1), ...
    moveonset_aligned_pyr,'uni',false);
temporal_act_avg_corr = cell(7,15);
temporal_act_trial_corr = cell(7,15);
for i = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
    for j = sessions
        curr_data = reshape(permute(moveonset_aligned_pyr{i,j},[2 3 1]), ...
            [],size(moveonset_aligned_pyr{i,j},1));           
        % ignore nan values for correlation
        curr_data_mean = repmat(nanmean(curr_data),size(curr_data,1),1);
        curr_data(isnan(curr_data)) = curr_data_mean(isnan(curr_data));
        curr_template = mean_aligned_activity{i,j};
        curr_template(isnan(curr_template)) = nanmean(curr_template);
        
        curr_corr = corrcoef([curr_template curr_data]);
        temporal_act_avg_corr{i,j} = curr_corr(1,2:end);
        temporal_act_trial_corr{i,j} = AP_itril(curr_corr(2:end,2:end),-1);
    end
end

% shuffle time to subtract off non-temporal component
num_rep = 100;
use_animals = 1:7;
temporal_act_avg_corr_shuff_all = cell(8,15,num_rep);
temporal_act_trial_corr_shuff_all = cell(8,15,num_rep);
for curr_rep = 1:num_rep
    temporal_act_avg_corr_shuff = cell(7,15);
    temporal_act_trial_corr_shuff = cell(7,15);
    for i = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
        for j = sessions
            
            % randomly circle shift each cell/trial
            curr_data = permute(moveonset_aligned_pyr{i,j},[2 3 1]);
            curr_data_shifted = nan(size(curr_data));
            
            [m n p] = size(curr_data);
            
            for curr_p = 1:p
                D = randi([0, m - 1], [1, n]);       
                for ii = (1 : n)
                    curr_data_shifted(:,ii,curr_p) = ...
                        [curr_data((m - D(ii) + 1 : m),ii,curr_p); ...
                        curr_data((1 : m - D(ii) ),ii,curr_p)];
                end
            end
            
            curr_data = reshape(curr_data_shifted, ...
                [],size(moveonset_aligned_pyr{i,j},1));
            % ignore nan values for correlation
            curr_data_mean = repmat(nanmean(curr_data),size(curr_data,1),1);
            curr_data(isnan(curr_data)) = curr_data_mean(isnan(curr_data));
            curr_template = mean_aligned_activity{i,j};
            curr_template(isnan(curr_template)) = nanmean(curr_template);
            
            curr_corr = corrcoef([curr_template curr_data]);
            temporal_act_avg_corr_shuff{i,j} = curr_corr(1,2:end);
            temporal_act_trial_corr_shuff{i,j} = AP_itril(curr_corr(2:end,2:end),-1);
        end
        i
    end
    temporal_act_avg_corr_shuff_all(:,:,curr_rep) = ...
        temporal_act_avg_corr_shuff;
    temporal_act_trial_corr_shuff_all(:,:,curr_rep) = ...
        temporal_act_trial_corr_shuff;
    curr_rep
end

% trial-trial correlation
temporal_act_trial_corr_shuff_all_cat = cell(7,15);
for i = 1:8
    for j = 1:15
        temporal_act_trial_corr_shuff_all_cat{i,j} = ...
            horzcat(temporal_act_trial_corr_shuff_all{i,j,:});
    end
end

act_trial_shuff_sub = cellfun(@(x,y) nanmean(x,2)-nanmean(y,2), ...
    temporal_act_trial_corr, temporal_act_trial_corr_shuff_all_cat,'uni',false);
act_trial_shuff_sub_cat = cell(14,1);
use_animals = 1:7;
for i = 1:14
   act_trial_shuff_sub_cat{i} = vertcat(act_trial_shuff_sub{use_animals,i}); 
end
% get total shuff subtracted mean
act_trial_shuff_sub_m = cellfun(@nanmean,act_trial_shuff_sub);
figure;errorbar(nanmean(act_trial_shuff_sub_m), ...
    nanstd(act_trial_shuff_sub_m)./sqrt(sum(~isnan(act_trial_shuff_sub_m))),'k','linewidth',2)
title('trial-trial correlation - shifted subtract')
figure;errorbar(cellfun(@nanmean,act_trial_shuff_sub_cat), ...
    cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),act_trial_shuff_sub_cat), ...
    'k','linewidth',2);
title('trial-trial correlation - shifted subtract')


% trial-average correlation
temporal_act_avg_corr_shuff_all_cat = cell(7,15);
for i = 1:8
    for j = 1:15
        temporal_act_avg_corr_shuff_all_cat{i,j} = ...
            vertcat(temporal_act_avg_corr_shuff_all{i,j,:});
    end
end

act_avg_shuff_sub = cellfun(@(x,y) nanmean(x,1)-nanmean(y,1), ...
    temporal_act_avg_corr, temporal_act_avg_corr_shuff_all_cat,'uni',false);
act_avg_shuff_sub_cat = cell(14,1);
use_animals = 1:7;
for i = 1:14
   act_avg_shuff_sub_cat{i} = horzcat(act_avg_shuff_sub{use_animals,i}); 
end
% get total shuff subtracted mean
act_avg_shuff_sub_m = cellfun(@nanmean,act_avg_shuff_sub);
figure;errorbar(nanmean(act_avg_shuff_sub_m), ...
    nanstd(act_avg_shuff_sub_m)./sqrt(sum(~isnan(act_avg_shuff_sub_m))),'k','linewidth',2)
title('trial-average correlation')
figure;errorbar(cellfun(@nanmean,act_avg_shuff_sub_cat), ...
    cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),act_avg_shuff_sub_cat), ...
    'k','linewidth',2);
title('trial-average correlation')





%% Figure 3: temporal alignment example

% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 28;fs_f = 28*3;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

% this is to look for good examples

% plot class grid
curr_animal = 5;
figure;imagesc(horzcat(pyr_class_all{curr_animal,:})');colormap(gray)

% plot aligned activity of all high-reliable cells
reliability_cutoff = 0.2;
reliability = cellfun(@(x) permute(nanmean(any(x>0,2),1),[3 1 2]), ...
    moveonset_aligned_pyr,'uni',false);
high_reliability_cells = cell(7,1);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_reliability_cat = horzcat(reliability{curr_animal,sessions});
    curr_class_cat = horzcat(pyr_class_all{curr_animal,sessions});
    
    curr_reliability_cat(~curr_class_cat) = NaN;
    high_reliability_cells{curr_animal} = find(nanmean(curr_reliability_cat,2) > ...
        reliability_cutoff);
end
curr_animal = 4;
for curr_cell = 1:length(high_reliability_cells{curr_animal})
    figure;
    curr_pyrcell = high_reliability_cells{curr_animal}(curr_cell);
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for session = sessions
        subplot(4,4,session);colormap(gray);
        imagesc(moveonset_aligned_pyr{curr_animal,session}(:,:,curr_pyrcell)>0);
        axis off;
        if session == sessions(1)
            title(['Cell ' num2str(curr_pyrcell)]);
        end
    end
end
set(gcf,'Name',['Animal ' num2str(curr_animal)]);

% plot example mean avg day plot
curr_animal = 5;
curr_session = 1;
a = permute(nanmean(moveonset_aligned_pyr{curr_animal,curr_session}),[3 2 1]);
figure;imagesc(a)

% plot the trials from animal/cell
curr_animal = 6;
curr_pyrcell = 37;
figure; 
sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
for session = sessions
    subplot(4,4,session);colormap(gray);
    imagesc(moveonset_aligned_pyr{curr_animal,session}(:,:,curr_pyrcell)>0);
    axis off;
end
set(gcf,'Name',['Animal ' num2str(curr_animal) ', cell ' num2str(curr_pyrcell)]);

% FINAL SINGLE TEMPORAL EXAMPLE plot single cell example for figure
curr_animal = 6;
curr_pyrcell = 37;
use_sessions = [4 8 14];
use_num_trials = min(cellfun(@(x) sum(any(x(:,:,curr_pyrcell),2)), ...
    moveonset_aligned_pyr(curr_animal,use_sessions)));
figure;
for session = 1:length(use_sessions)
    subplot(2,3,session);colormap(gray);
    active_trials = find(any(moveonset_aligned_pyr{curr_animal, ...
        use_sessions(session)}(:,:,curr_pyrcell),2));
    use_trials = sort(active_trials(randperm(length(active_trials),use_num_trials)));
    plot_activity = moveonset_aligned_pyr{curr_animal, ...
        use_sessions(session)}(use_trials,:,curr_pyrcell)>0;
    imagesc(plot_activity);
    
    subplot(2,3,session+3);
    plot(nanmean(moveonset_aligned_pyr{curr_animal, ...
        use_sessions(session)}(active_trials,:,curr_pyrcell)>0),'k','linewidth',2)
    xlim([0 size(plot_activity,2)+1]);ylim([0 1])
end
set(gcf,'Name',['Animal ' num2str(curr_animal) ', cell ' num2str(curr_pyrcell)]);

    



% make matricies of first onsets
moveonset_aligned_pyr_first_onset = cell(size(moveonset_aligned_pyr));
for i = 1:8
    for j = 1:15
        [c first_peak] = max(moveonset_aligned_pyr{i,j},[],2);
        x = zeros(size(moveonset_aligned_pyr{i,j}));
        for curr_cell = 1:size(x,3)
            for trial = 1:size(first_peak,1)
                x(trial,first_peak(trial,1,curr_cell),curr_cell) = c(trial,1,curr_cell);
            end
        end
        moveonset_aligned_pyr_first_onset{i,j} = x;
    end
end

% EXAMPLE 1

% plot examples
curr_animal = 1;
curr_cells = [79 32 241];
curr_days = [5 7 11];
figure
for i = 1:length(curr_days)
    for j = 1:length(curr_cells)
        subplot(length(curr_days),length(curr_cells), ...
            (i-1)*length(curr_cells)+j) 
        
        imagesc(moveonset_aligned_pyr_first_onset{curr_animal,curr_days(i)} ...
            (:,:,curr_cells(j)) > 0);
        
        colormap(gray);
        
        line([61 61],ylim,'color','r','linewidth',2);
        
        if j == 1
            ylabel(['Day ' num2str(curr_days(i))])
        end
        
        if i == 1
            title(['Animal ' num2str(curr_animal) ', cell ' num2str(curr_cells(j))])
        end
        
    end    
end

% plot example scalebar
figure;
subplot(length(curr_days),length(curr_cells),1);
line([0 3*28],[0 0],'linewidth',3,'color','k');
title('3 seconds');
xlim([0 fs_b+fs_f+1]);

% plot both cells together in one block for first and last day used
a = moveonset_aligned_pyr_onset{curr_animal,curr_days(end)}(:,:,curr_cells(1))>0;
a(:,:,2) = moveonset_aligned_pyr_onset{curr_animal,curr_days(end)}(:,:,curr_cells(2))>0;
a(:,:,3) = moveonset_aligned_pyr_onset{curr_animal,curr_days(end)}(:,:,curr_cells(3))>0;
figure;imagesc(a)
line([61 61],ylim,'color','y','linewidth',2);

[c a1_idx] = max(a(:,:,1),[],2);
a1_idx_long = cellfun(@(x) x-5:x+50,num2cell(a1_idx),'uni',false);
use_trials = find(cellfun(@(x) all(x < size(a,2)) & ...
    all(x > 0),a1_idx_long) & sum(any(a,2),3) > 1);
a_1start = cell2mat(arrayfun(@(x) a(use_trials(x),a1_idx_long{use_trials(x)},:), ...
    1:length(use_trials),'uni',false)');
[c ci] = max(a_1start(:,:,2),[],2);
[asdf sort_idx] = sort(ci);
figure;imagesc(a_1start(sort_idx,:,:));

a = moveonset_aligned_pyr_onset{curr_animal,curr_days(1)}(:,:,curr_cells(1))>0;
a(:,:,2) = moveonset_aligned_pyr_onset{curr_animal,curr_days(1)}(:,:,curr_cells(2))>0;
a(:,:,3) = moveonset_aligned_pyr_onset{curr_animal,curr_days(1)}(:,:,curr_cells(3))>0;
figure;imagesc(a)
line([61 61],ylim,'color','y','linewidth',2);

[c a1_idx] = max(a(:,:,1),[],2);
a1_idx_long = cellfun(@(x) x-5:x+50,num2cell(a1_idx),'uni',false);
use_trials = find(cellfun(@(x) all(x < size(a,2)) & ...
    all(x > 0),a1_idx_long) & sum(any(a,2),3) > 1);
a_1start = cell2mat(arrayfun(@(x) a(use_trials(x),a1_idx_long{use_trials(x)},:), ...
    1:length(use_trials),'uni',false)');
[c ci] = max(a_1start(:,:,2),[],2);
[asdf sort_idx] = sort(ci);
figure;imagesc(a_1start(sort_idx,:,:));



use_trials = find(sum(any(a,2),3) > 1);
[c ci] = max(a(:,:,1),[],2);
[asdf sort_idx] = sort(ci(use_trials));
figure;imagesc(a(use_trials(sort_idx),:,:));


% EXAMPLE 2

% plot examples
curr_animal = 4;
curr_cells = [75 14];
curr_days = [3 7 12];
figure
for i = 1:length(curr_days)
    for j = 1:length(curr_cells)
        subplot(length(curr_days),length(curr_cells), ...
            (i-1)*length(curr_cells)+j) 
        
        imagesc(moveonset_aligned_pyr{curr_animal,curr_days(i)} ...
            (:,:,curr_cells(j)) > 0);
        
        colormap(gray);
        
        line([61 61],ylim,'color','r','linewidth',2);
        
        if j == 1
            ylabel(['Day ' num2str(curr_days(i))])
        end
        
        if i == 1
            title(['Animal ' num2str(curr_animal) ', cell ' num2str(curr_cells(j))])
        end
        
    end    
end

% plot example scalebar
figure;
subplot(length(curr_days),length(curr_cells),1);
line([0 3*28],[0 0],'linewidth',3,'color','k');
title('3 seconds');
xlim([0 fs_b+fs_f+1]);

% plot both cells together in one block for last day used
a = moveonset_aligned_pyr{curr_animal,curr_days(end)}(:,:,curr_cells(1))>0;
a(:,:,2) = moveonset_aligned_pyr{curr_animal,curr_days(end)}(:,:,curr_cells(2))>0;
a(:,:,3) = zeros(size(a(:,:,1)));
figure;imagesc(a)
line([61 61],ylim,'color','y','linewidth',2);


%%%%%%%%% FINAL RASTERPLOT FIGURE
% find move onset to +3s
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
reward_aligned_pyr = cell(size(movement_trace_all));
reward_aligned_pyr_onset = cell(size(movement_trace_all));
movement_lengths = cell(7,15);
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 28;fs_f = 28*3;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
    curr_rewardmoves(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
   curr_reward_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_reward_surround,'uni',false);
   curr_reward_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_reward_surround,'uni',false);
   
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    
    reward_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_reward_aligned_pyr{:}), [3 2 1]);
    reward_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_reward_aligned_pyr_onsets{:}), [3 2 1]);

    
    movement_lengths{curr_animalday} = cellfun(@length,curr_rewardmoves);      
    
    curr_animalday/max(animaldays)

end
% make matricies of first onsets
moveonset_aligned_pyr_first_onset = cell(size(moveonset_aligned_pyr_onset));
for i = 1:7
    for j = 1:15
        [c first_peak] = max(moveonset_aligned_pyr_onset{i,j}>0,[],2);
        x = zeros(size(moveonset_aligned_pyr_onset{i,j}));
        for curr_cell = 1:size(x,3)
            for trial = 1:size(first_peak,1)
                x(trial,first_peak(trial,1,curr_cell),curr_cell) = c(trial,1,curr_cell);
            end
        end
        moveonset_aligned_pyr_first_onset{i,j} = x;
    end
end
reward_aligned_pyr_first_onset = cell(size(reward_aligned_pyr_onset));
for i = 1:7
    for j = 1:15
        [c first_peak] = max(reward_aligned_pyr_onset{i,j}>0,[],2);
        x = zeros(size(reward_aligned_pyr_onset{i,j}));
        for curr_cell = 1:size(x,3)
            for trial = 1:size(first_peak,1)
                x(trial,first_peak(trial,1,curr_cell),curr_cell) = c(trial,1,curr_cell);
            end
        end
        reward_aligned_pyr_first_onset{i,j} = x;
    end
end


% plot all instances of a cell's firing on top of each other

%reliability = cellfun(@(x) permute(nanmean(any(x>0,2),1),[3 1 2]), ...
%    moveonset_aligned_pyr,'uni',false);

reliability_cutoff = 0.1;

curr_animal = 2;
sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));

all_sessions = {[1 2 3] [5 6 7] [12 13 14]};
for curr_session_group = 1:length(all_sessions)
    
    curr_sessions = all_sessions{curr_session_group};
    
    curr_data = vertcat(moveonset_aligned_pyr{curr_animal,curr_sessions});
    curr_data_onset = vertcat(moveonset_aligned_pyr_first_onset{curr_animal,curr_sessions});
    curr_reward_data_onset = vertcat(reward_aligned_pyr_first_onset ...
        {curr_animal,curr_sessions});
    
    curr_movelengths = horzcat(movement_lengths{curr_animal,curr_sessions});
    
    curr_mean_act = permute(nanmean(curr_data_onset > 0),[3 2 1]);
    [c ci] = max(curr_mean_act(:,fs_b:end),[],2);
    curr_reliability = permute(nanmean(any(curr_data_onset(:,fs_b-5:end,:),2)),[3 1 2]);
    reliable_cells = find(curr_reliability > reliability_cutoff);
    [asdf sort_idx] = sortrows([ci(reliable_cells) c(reliable_cells)],[1 -2]);
    sort_cells = reliable_cells(sort_idx);
    
    trial_active = arrayfun(@(x) find(any(curr_data_onset ...
        (:,:,x),2)),sort_cells,'uni',false);
    
    use_movelengths = cellfun(@(x) curr_movelengths(x), ...
        trial_active,'uni',false);
    
    % only use some number of trials if > 20 trials 
    %plot_active = arrayfun(@(x) curr_data_onset ...
    %    (trial_active{x}(randperm(length(trial_active{x}),20)), ...
    %    :,sort_cells(x))>0,1:length(sort_cells),'uni',false);
   
     plot_active = arrayfun(@(x) curr_data_onset ...
         (trial_active{x},:,sort_cells(x))>0,1:length(sort_cells),'uni',false);
    
    %col = lines(length(sort_cells));
    %plot_active_colored = arrayfun(@(x) repmat(plot_active{x},[1 1 3]).* ...
    %    repmat(permute(col(x,:),[1 3 2]),size(plot_active{x},1), ...
    %    size(plot_active{x},2)),1:length(sort_cells),'uni',false);
    %plot_active_cat = vertcat(plot_active_colored{:});    
    %subplot(1,max(sessions),curr_session);
    %figure
    %imagesc(plot_active_cat);
    
    % plot as points
    [c activity_times] = cellfun(@(x) max(x,[],2),plot_active,'uni',false);
    trial_lengths = [0 cumsum(cellfun(@length,activity_times(1:end-1)))];
    col = lines(length(sort_cells));
    col = col(circshift(1:size(col,1),[0 randi(size(col,1))]),:);
    figure; hold on;
    for i = 1:length(activity_times)
        for j = 1:length(activity_times{i})            
            plot(activity_times{i}(j),j+trial_lengths(i),'.','color', ...
                col(i,:),'linewidth',2)
            % to rescale time
%             plot((activity_times{i}(j)-5)/use_movelengths{i}(j), ...
%                 j+trial_lengths(i),'.','color', ...
%                 col(i,:),'linewidth',2)
        end
    end
    ylim([0 sum(cellfun(@length,activity_times))+1])
    title(num2str(curr_sessions));
    line([fs_b+1 fs_b+1],ylim,'color','k','linewidth',2)
end


% Plot mean activity of many cells
use_animals = 1:7;
curr_sessions = [12:14];

reliability_cutoff = 0.1;
curr_onset_cat = arrayfun(@(x) vertcat( ...
    moveonset_aligned_pyr_first_onset{x,curr_sessions}),use_animals,'uni',false);
reliability_cells = cellfun(@(x) permute(nanmean(any(x,2)),[3 1 2]) > ...
    reliability_cutoff,curr_onset_cat,'uni',false);

curr_onset_mean = arrayfun(@(x) permute(nanmean( ...
    curr_onset_cat{x}(:,:,reliability_cells{x})),[3 2 1]), ...
    1:length(use_animals),'uni',false);
curr_onset_mean_cat = vertcat(curr_onset_mean{:});

[c ci] = max(curr_onset_mean_cat,[],2);
[asdf sort_idx] = sort(ci);
curr_onset_mean_cat_norm = bsxfun(@times,curr_onset_mean_cat, ...
    1./max(curr_onset_mean_cat,[],2));
curr_onset_mean_cat_conv = conv2(curr_onset_mean_cat,h,'same');
curr_onset_mean_cat_conv_norm = bsxfun(@times,curr_onset_mean_cat_conv, ...
    1./max(curr_onset_mean_cat_conv,[],2));
figure;imagesc(curr_onset_mean_cat_conv_norm(sort_idx,:));colormap(hot)



%% Figure 3: heatmap of normalized movement related mean activity

% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 20;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

mean_moveonset_pyr = cellfun(@(x,y) permute(nanmean(x(:,:,y),1),[3 2 1]), ...
    moveonset_aligned_pyr,pyr_class_all,'uni',false);
use_animals = 1:7;
early_moveonset_cat = vertcat(mean_moveonset_pyr{use_animals,2});
early_moveonset_cat_norm = bsxfun(@times,early_moveonset_cat,1./ ...
    max(early_moveonset_cat,[],2));
[c ci] = max(early_moveonset_cat_norm,[],2);
[asdf early_sort] = sort(ci);
exclude_cells = ~any(early_moveonset_cat_norm(early_sort,:),2) | ci(early_sort) > 104;
early_sort(exclude_cells) = [];

late_moveonset_cat = vertcat(mean_moveonset_pyr{use_animals,14});
late_moveonset_cat_norm = bsxfun(@times,late_moveonset_cat,1./ ...
    max(late_moveonset_cat,[],2));
[c ci] = max(late_moveonset_cat_norm,[],2);
[asdf late_sort] = sort(ci);
exclude_cells = ~any(late_moveonset_cat_norm(late_sort,:),2) | ci(late_sort) > 104;
late_sort(exclude_cells) = [];

figure; colormap(hot)
subplot(1,2,1);imagesc(early_moveonset_cat_norm(early_sort,:));
line([fs_b+1 fs_b+1],ylim,'color','b','linewidth',2);
subplot(1,2,2);imagesc(late_moveonset_cat_norm(late_sort,:));
line([fs_b+1 fs_b+1],ylim,'color','b','linewidth',2);

figure; colormap(hot)
subplot(1,2,1);imagesc(early_moveonset_cat(early_sort,:));
subplot(1,2,2);imagesc(late_moveonset_cat(late_sort,:));
%% Figure 4: template-based lever/activity

% get activity with lever 
movement_activity = cell(7,15);
act_trials_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        

        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);
  

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% get relevant lever traces
use_samples = 10001:100:40001;
lever_trials_used = cell(7,15);
lever_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);

    % extracting all lever
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        %curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
        %    reward_aligned_lever{session},1));
        % try using cued and uncued trials
        %curr_cued = true(size(curr_cued));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        total_curr_cat_lever{session} = temp(use_trials,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    lever_used(curr_animal,sessions) = total_curr_cat_lever;
    % creating template
    for session = sessions(end-2:end);
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        % try using cued and uncued trials
        %curr_cued = true(size(curr_cued));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        curr_cat_lever{session} = temp(use_trials,use_samples);
    end
    disp(curr_animal)
end


% use only trials with both activity and lever traces
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);

movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

% since random trials are picked for template: do this many times and take
% the mean results
num_rep = 1;
corr_bins = [-1:0.2:1-0.2 Inf];
m = nan(num_rep,length(corr_bins)-1,2);
s = nan(num_rep,length(corr_bins)-1,2);
r = nan(num_rep,2);
p = nan(num_rep,2);
anova_p = nan(num_rep,2);
bin_ranksum = cell(num_rep,1);

corr_fig = figure; hold on;
act_fig = figure;
lever_fig = figure;

for curr_rep = 1:num_rep
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    
    % For this purpose - template is just half of the trials in last 3 days
    lever_avg_template_trials = cell(7,15);
    lever_avg_template = cell(7,1);
    lever_avg_template_corr = cell(7,15);
    lever_avg_avg_template_corr = nan(7,15);
    activity_template_all = cell(7,1);
    activity_template_corr = cell(7,15);
    for curr_animal = 1:7
        
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        
        % pick random subset of trials from last 3 days to make lever template
        template_sessions = sessions(end-2:end);
        curr_lever_template_trials = cell(1,max(sessions));
        curr_lever_template_trials(template_sessions) = ...
            cellfun(@(x) randsample(size(x,1),round(size(x,1)/2)), ...
            lever_used_crop(curr_animal,template_sessions),'uni',false);
        lever_avg_template_trials(curr_animal,1:max(sessions)) = ...
            curr_lever_template_trials;
        
        % make lever template
        lever_avg_use = cellfun(@(x,y) x(y,:), ...
            lever_used_crop(curr_animal,1:max(sessions)), ...
            curr_lever_template_trials,'uni',false);
        curr_template = nanmean(vertcat(lever_avg_use{:}));
        lever_avg_template{curr_animal} = curr_template;
        
        % get the lever correlation for trials not used to make template
        lever_avg_template_corr_grid = cellfun(@(x) ...
            corrcoef([x' curr_template']), ...
            lever_used_crop(curr_animal,sessions),'uni',false);
        lever_avg_template_corr(curr_animal,sessions) = cellfun(@(x,y) ...
            x(end,setdiff(1:size(x,2)-1,y)),lever_avg_template_corr_grid, ...
            curr_lever_template_trials(sessions),'uni',false);
        
        % get lever correlation avg of day with template
        curr_lever_avg = cellfun(@(x,y) nanmean(x(setdiff(1:size(x,1),y),:)), ...
            lever_used_crop(curr_animal,sessions), ...
            curr_lever_template_trials(sessions),'uni',false);
        curr_lever_avg_cat = vertcat(curr_lever_avg{:});
        lever_avg_avg_template_corr_grid = corrcoef( ...
            [curr_lever_avg_cat' curr_template']);
        lever_avg_avg_template_corr(curr_animal,sessions) = ...
            lever_avg_avg_template_corr_grid(end,1:end-1);
        
        % make activity template from same trials used for lever template
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            curr_lever_template_trials(sessions),'uni',false);
        
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        activity_template_all{curr_animal} = activity_template;
        
        % get the activity correlation for trials not used to make template
        act_avg_template_corr_grid = cellfun(@(x) ...
            corrcoef([x activity_template]), ...
            movement_activity_crop(curr_animal,sessions),'uni',false);
        activity_template_corr(curr_animal,sessions) = cellfun(@(x,y) ...
            x(end,setdiff(1:size(x,2)-1,y)),act_avg_template_corr_grid, ...
            curr_lever_template_trials(sessions),'uni',false);
        
        curr_animal
    end
    
    % combine all animals to get total activity and correlations
    use_animals = 1:7;
    curr_cc_data = cell(2,1);
    curr_bin_data = cell(2,1);
    for earlylate = 1:2
        
        % they need to be sorted for down the road
        use_animals = sort(use_animals);
        
        if earlylate == 1
            curr_sessions = 1:3;
            plot_col = 'k';
        elseif earlylate == 2
            curr_sessions = 10:14;
            plot_col = 'r';
        end
        
        total_lever_cc = cell(7,1);
        total_act_cc = cell(7,1);
        total_act = cell(7,1);
        total_lever = cell(7,1);
        for curr_animal = use_animals;
            sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
            use_sessions = intersect(curr_sessions,sessions);
            
            cat_lever_cc = horzcat(lever_avg_template_corr{curr_animal,use_sessions});
            cat_act_cc = horzcat(activity_template_corr{curr_animal,use_sessions});
            
            curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
                movement_activity_crop(curr_animal,use_sessions), ...
                lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
            curr_act_cat = horzcat(curr_act{:});
            
            curr_lever = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
                lever_used_crop(curr_animal,use_sessions), ...
                lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
            curr_lever_cat = vertcat(curr_lever{:});
            
            total_lever_cc{curr_animal} = cat_lever_cc;
            total_act_cc{curr_animal} = cat_act_cc;
            total_act{curr_animal} = curr_act_cat;
            total_lever{curr_animal} = curr_lever_cat;
        end
        
        total_lever_cc_cat = horzcat(total_lever_cc{:});
        total_act_cc_cat = horzcat(total_act_cc{:});
        
        
        [n bins] = histc(total_lever_cc_cat,corr_bins);
        curr_cc_data{earlylate} = total_act_cc_cat;
        curr_bin_data{earlylate} = bins;
        
        m(curr_rep,unique(bins),earlylate) = grpstats(total_act_cc_cat,bins,'nanmean');
        s(curr_rep,unique(bins),earlylate) = grpstats(total_act_cc_cat,bins,'sem');
        [curr_r curr_p] = robustfit(bins,total_act_cc_cat);
        r(curr_rep,earlylate) = curr_r(2);
        p(curr_rep,earlylate) = curr_p.p(2);

        % plot the average activity of cells in bins (only on the last rep)
        
        if curr_rep == num_rep
            %act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
            act_corr_bins = [-1 -0.1667 0.1667 0.5 Inf];
            
            [n bin_split] = cellfun(@(x) histc(x,act_corr_bins),total_lever_cc,'uni',false);
            act_bins = horzcat(bin_split{:});
            act_grp_all = cell(7,1);
            lever_grp_all = cell(7,1);
            for curr_animal = use_animals
                act_grp = nan(size(total_act{curr_animal},1),length(act_corr_bins));
                act_grp(:,unique(bin_split{curr_animal})) = ...
                    grpstats([total_act{curr_animal}]',bin_split{curr_animal})';
                act_grp_reshape = reshape(act_grp,move_time_framespan+2*leeway+1, ...
                    [],length(act_corr_bins));
                act_grp_all{curr_animal} = act_grp_reshape;
                
                lever_grp = nan(size(lever_avg_template{curr_animal},2),length(act_corr_bins));
                lever_grp(:,unique(bin_split{curr_animal})) = ...
                    grpstats([total_lever{curr_animal}],bin_split{curr_animal})';
                lever_grp_all{curr_animal} = lever_grp;
            end
            
            act_grp_cat = cat(2,act_grp_all{:});
            
            activity_template_cat = vertcat(activity_template_all{use_animals});
            activity_template_reshape = reshape(activity_template_cat, ...
                move_time_framespan+2*leeway+1,[])';
            [c ci] = max(activity_template_reshape,[],2);
            use_cells = find(c > 0);
            
            % restrict use cells only to those movement-related on the
            % later sessions
            cat_class = cell(7,1);
            for i = use_animals
                cat_class{i} = any(horzcat(pyr_class_all{i,10:14}),2);
            end
            cat_class_idx = find(vertcat(cat_class{:}));
            
            use_cells = intersect(use_cells,cat_class_idx);

            [asdf sort_idx] = sort(ci(use_cells));
            
            activity_template_reshape_norm = bsxfun(@times,activity_template_reshape,1./c);
            
            figure(act_fig);
            for i = 1:length(act_corr_bins)-1
                subplot(2,length(act_corr_bins),i+((earlylate-1)*length(act_corr_bins)));
                
                imagesc(act_grp_cat(:,use_cells(sort_idx),i)');
                
                colormap(hot);
                caxis([0 0.2]);
                axis off
                title([num2str(act_corr_bins(i)) ...
                    ':' num2str(act_corr_bins(i+1))]);
            end
            subplot(2,length(act_corr_bins),length(act_corr_bins)+ ...
                ((earlylate-1)*length(act_corr_bins)));
            imagesc(activity_template_reshape(use_cells(sort_idx),:));
            caxis([0 0.2]);
            colormap(hot);
            axis off
            title('Template')
            colormap(hot);
            
            set(gcf,'Name',['Animals ' num2str(use_animals)]);
            
            figure(lever_fig);
            for plot_animal = use_animals
                for plot_bin = 1:length(act_corr_bins)
                    subplot(8,length(act_corr_bins),(length(act_corr_bins))*(plot_animal-1) + plot_bin)
                    hold on;
                    plot(lever_grp_all{plot_animal}(:,plot_bin),'color',plot_col,'linewidth',2);
                end
                subplot(8,length(act_corr_bins),(length(act_corr_bins)* ...
                    (plot_animal-1)) + length(act_corr_bins))
                plot(lever_avg_template{plot_animal},'color',plot_col,'linewidth',2);
            end
        end
        
    end

    curr_cc_data_cat = horzcat(curr_cc_data{:});
    % anova_p saved as 1) bin effect, 2) group effect
    g1 = horzcat(curr_bin_data{:});
    g2 = [ones(size(curr_cc_data{1})) 2*ones(size(curr_cc_data{2}))];
    anova_p(curr_rep,:) = anovan(curr_cc_data_cat(~isnan(curr_cc_data_cat)), ...
        {g1(~isnan(curr_cc_data_cat)) g2(~isnan(curr_cc_data_cat))},'display','off');

    bin_ranksum{curr_rep} = nan(length(corr_bins)-1,1);
    for curr_bin = intersect(unique(curr_bin_data{1}),unique(curr_bin_data{2}))
        curr_early_data = curr_cc_data{1}(curr_bin_data{1} == curr_bin);
        curr_late_data = curr_cc_data{2}(curr_bin_data{2} == curr_bin);
        bin_ranksum{curr_rep}(curr_bin) = ranksum(curr_early_data(~isnan(curr_early_data)), ...
            curr_late_data(~isnan(curr_late_data)));
    end
    

    curr_rep
end

figure(corr_fig);
m_plot = nanmean(m,1);
s_plot = nanmean(s,1);
errorbar(corr_bins(1:end-1),m_plot(1,:,1),s_plot(1,:,1),'k','linewidth',2);
errorbar(corr_bins(1:end-1),m_plot(1,:,2),s_plot(1,:,2),'r','linewidth',2);
xlabel('Lever template correlation')
ylabel('Activity template correlation')
title(['Animals ' num2str(use_animals)])
legend({'Early' 'Late'})

% statistics 
mean_r = nanmean(r);
mean_p = nanmean(p);
mean_anova_p = nanmean(anova_p);

% 
% % plot the average lever trace for each used animal in each bin
% use_animals = 1:8;
% for earlylate = 1:2
%     
%     if earlylate == 1
%         curr_sessions = 1:3;
%         figure;
%         col = [0 0 1];
%     elseif earlylate == 2
%         curr_sessions = 8:14;
%         col = [1 0 0];
%     end
%         
%     total_lever_cc = cell(7,1);
%     for curr_animal = use_animals;
%         sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
%         use_sessions = intersect(curr_sessions,sessions);
% 
%         cat_lever_cc = horzcat(lever_avg_template_corr{curr_animal,use_sessions});
% 
%         total_lever_cc{curr_animal} = cat_lever_cc;      
%     end
%       
%     total_lever_cc_cat = horzcat(total_lever_cc{:});
%        
%     act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
%     
%     [n act_bins] = histc(total_lever_cc_cat,act_corr_bins);
%     
%     bin_split = mat2cell(act_bins,1,cellfun(@length,total_lever_cc));
%     
%     lever_used_bin = nan(length(act_corr_bins)-1,302,8);
%     lever_used_bin_sem = nan(length(act_corr_bins)-1,302,8);
%     for curr_animal = use_animals
%         
%         sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
%         use_sessions = intersect(curr_sessions,sessions);
%         
%         lever_used_crop_nontemplate = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
%             lever_used_crop(curr_animal,use_sessions), ...
%             lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
%         lever_cat = vertcat(lever_used_crop_nontemplate{:});
%         
%         lever_used_bin(unique(bin_split{curr_animal}),:,curr_animal) = ...
%             grpstats(lever_cat,bin_split{curr_animal});
%         lever_used_bin_sem(unique(bin_split{curr_animal}),:,curr_animal) = ...
%             grpstats(lever_cat,bin_split{curr_animal},'sem');
%         
%     end
%     
%     for curr_animal = use_animals
%         for i = 1:length(act_corr_bins)-1
%             subplot(length(act_corr_bins),8,8*(i-1)+curr_animal);hold on;
%             plot(lever_used_bin(i,:,curr_animal),'color',col,'linewidth',2);
%             
%             if earlylate == 1
%                 plot(lever_avg_template{curr_animal}, ...
%                     '--k','linewidth',2)
%             end           
%             
%             if i == 1
%                 title(curr_animal)
%             end
%             
%             if curr_animal == 1
%                 ylabel(num2str(act_corr_bins(i)));
%             end
% 
%             axis off
%         end
%     end
% end

%% Supplementary ?: fraction of rewarded trials

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' ...
    'AP78' 'AP79' 'SC027' 'SC028' 'SC029'};
rewarded_trials_all = cell(length(animals),15);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    filename_delimiter = cellfun(@(x) strfind(x,filesep),data_days,'uni',false);
    days = cellfun(@(x,y) x(y(end)+1:end),data_days,filename_delimiter, ...
        'UniformOutput',false);
    
    for curr_day = 1:length(days)
     
        day = num2str(days{curr_day});
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep animal(1:2) '00' animal(end-1:end)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % if no xsg (because accidentally not recorded), skip
        if isempty(xsg_fullfilenames)
            continue
        end
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        
        rewarded_trials_all{curr_animal,curr_day} = ...
            cellfun(@(x) ~isempty(x.states.reward), bhv);
        disp(curr_day)
    end
    disp(['animal:' animal])
end

use_animals = 1:10;
rewarded_trials_frac_use = cellfun(@nanmean,rewarded_trials_all(use_animals,1:14));
figure; hold on;
plot(rewarded_trials_frac_use','color',[0.5 0.5 0.5]);
errorbar(nanmean(rewarded_trials_frac_use), ...
    nanstd(rewarded_trials_frac_use)./sqrt(sum( ...
    ~isnan(rewarded_trials_frac_use))),'k','linewidth',2)

p = anova(rewarded_trials_frac_use,[],'off');






%% Supplementary 1: Spatial clustering

% get all roi distances and gad/pyr unfilled indicies

% set pixels to microns ratio
um_px_x = 590/512;
um_px_y = 500/512;

roi_distances = cell(7,1);
pyr_unfilled_idx = cell(7,1);
gad_unfilled_idx = cell(7,1);

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    pyr_unfilled_idx{curr_animal} = find(pyr_cells & ~filled_cells);
    gad_unfilled_idx{curr_animal} = find(gad_cells & ~filled_cells);
    
    % load ROIs and get ROI centroids and distances
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    roi_dir = dir(roi_path);
    roi_filenames = {roi_dir(:).name};
    roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
    roi_files_sort = sort(roi_filenames(roi_files));
    load([roi_path filesep roi_files_sort{end}],'-MAT');
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for curr_roi  = 1:num_rois
        [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
            polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
    end
    roi_center_x = roi_center_x*um_px_x;
    roi_center_y = roi_center_y*um_px_y;
    
    roi_dist = nan(length(roi_center_x));
    for k = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
                (roi_center_y(k) - roi_center_y(j))^2);
            roi_dist(k,j) = curr_dist;
        end
    end
    
    roi_distances{curr_animal} = roi_dist;
    
    disp(['Distances: animal ' animal]) 
    
end

% check clustering of movement cells

pyr_move_dist = cell(7,15);
pyr_move_time_diff = cell(7,15);
pyr_move_dist_p = nan(7,15);
aligned_max_grid = cell(7,15);
cell_dist_p = cell(7,15);
shuff_dist = cell(7,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    
    for session = sessions
        
        curr_move_cells = pyr_unfilled_idx{curr_animal}( ...
            pyr_class_all{curr_animal,session});
        curr_dist = AP_itril(roi_distances{curr_animal}( ...
            curr_move_cells,curr_move_cells),-1);
        pyr_move_dist{curr_animal,session} = curr_dist;
        
        curr_move_cell_dist_mean = nanmean(curr_dist);
        
        pyr_move_time_diff{curr_animal,session} = ... 
            AP_itril(pyr_time_diff_grid{curr_animal,session}( ...
            pyr_class_all{curr_animal,session}, ...
            pyr_class_all{curr_animal,session}),-1);
    
        % get probability distribution for random cells
        num_shuff = 10000;
        rand_pyr_unfilled = pyr_unfilled_idx{curr_animal}( ...
            randi(length(pyr_unfilled_idx{curr_animal}), ...
            length(curr_move_cells),num_shuff));
        
        rand_pyr_reshape = permute(rand_pyr_unfilled,[1 3 2]);
        rand_pyr_grid_1 = repmat(rand_pyr_reshape,1, ...
            size(rand_pyr_reshape,1));
        rand_pyr_grid_2 = permute(rand_pyr_grid_1,[2 1 3]);
        
        rand_pyr_idx = sub2ind(size(roi_distances{curr_animal}), ...
            rand_pyr_grid_1,rand_pyr_grid_2);
        rand_pyr_idx_tril = arrayfun(@(x) AP_itril( ...
            rand_pyr_idx(:,:,x),-1),1:num_shuff,'uni',false);
        
        rand_pyr_dist = cellfun(@(x) nanmean(roi_distances{curr_animal}(x)) ,...
            rand_pyr_idx_tril);
        
        shuff_dist{curr_animal,session} = rand_pyr_dist;
                
        curr_dist_rank = tiedrank([rand_pyr_dist curr_move_cell_dist_mean]);
        
        pyr_move_dist_p(curr_animal,session) = curr_dist_rank(end)/(num_shuff+1);
        
        cell_dist_p{curr_animal,session} = nan(size(curr_dist));
        for j = 1:length(curr_dist)
            curr_cell_rank = tiedrank([rand_pyr_dist curr_dist(j)]);
            cell_dist_p{curr_animal,session}(j) = ...
                curr_cell_rank(end)/(num_shuff+1);
        end           
        
    end
    
    curr_animal
    
end

% plot mean/std and real values of distances
figure;
for curr_animal = [1 2 4 5 6 7 8];
    curr_sessions = find(cellfun(@(x) ~isempty(x), shuff_dist(curr_animal,1:14)));
    curr_shuff_mean = cellfun(@nanmean,shuff_dist(curr_animal,curr_sessions));
    curr_shuff_ci_low = cellfun(@(x) prctile(x,2.5), ...
        shuff_dist(curr_animal,curr_sessions));
    curr_shuff_ci_high = cellfun(@(x) prctile(x,97.5), ...
        shuff_dist(curr_animal,curr_sessions));
    subplot(3,3,curr_animal); hold on;
    plot(curr_sessions,curr_shuff_mean,'k','linewidth',2)
    plot(curr_sessions,curr_shuff_ci_low,'--k','linewidth',2)
    plot(curr_sessions,curr_shuff_ci_high,'--k','linewidth',2)
    plot(curr_sessions,cellfun(@nanmean,pyr_move_dist(curr_animal,curr_sessions)), ...
        'or','MarkerSize',4,'MarkerFaceColor','r');
    xlim([0 15])
    ylim([100 400])
    set(gca,'XTick',2:2:14);
    set(gca,'YDir','reverse');
end

%% Supplementary 2: pyr/gad balance over days

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity with lever 
pyr_move_active = cell(7,15);
gad_move_active = cell(7,15);
movement_lengths = cell(7,15);
act_trials_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        %curr_moves = move_epoch_frames{curr_animal,curr_session}( ...
        %    rewarded_movements{curr_animal,curr_session});
        curr_moves = move_epoch_frames{curr_animal,curr_session};
      
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x) [leeway ...
            leeway],curr_moves,'uni',false);
        curr_moves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_moves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_moves_leeway);
        curr_moves_leeway(eliminate_movements) = [];
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_moves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_moves_leeway);
        
    end

    curr_animal
    
end

% get pyr:gad ratio for movements over days

use_animals = 1:7;
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);

figure; hold on
day_combine = {1:4 5:8 9:14};
col = jet(length(day_combine));
data_all = cell(size(day_combine));
bins_all = cell(size(day_combine));
for i = 1:length(day_combine)

    pc = horzcat(p{use_animals,day_combine{i}});
    gc = horzcat(g{use_animals,day_combine{i}});
    
    bin_edges = [0:0.02:0.25];
    %bin_edges = [0:0.1:1-0.1 Inf];
    [n bins] = histc(pc,bin_edges);
    gc_mean = grpstats(gc,bins);
    gc_sem = grpstats(gc,bins,'sem');
    use_bins = unique(bins) ~= 0;
    use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
    errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins), ...
        'color',col(i,:),'linewidth',2);
        
    data_all{i} = gc;
    bins_all{i} = bins;
end

pc = horzcat(p{use_animals,:});
gc = horzcat(g{use_animals,:});
bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'--k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement based activity')

legend([cellfun(@num2str,day_combine,'uni',false) 'All'])

% stats
g1 = horzcat(bins_all{:});
g2 = cellfun(@(x,y) x*ones(size(y)),num2cell(1:length(bins_all)),bins_all,'uni',false); 
g2_cat = horzcat(g2{:});
data_cat = horzcat(data_all{:});
[anova_p anova_stats] = anovan(data_cat(~isnan(data_cat)), ...
    {g1(~isnan(data_cat)) g2_cat(~isnan(data_cat))},'display','off');

p_g_ratio = cellfun(@(x,y) nanmedian(x./y),p,g);
sessions = repmat(1:14,7,1);
p_g_ratio_use = p_g_ratio(:,1:14);
[corr_r corr_p] = corrcoef(sessions(:),p_g_ratio_use(:),'rows','complete');
[anova_p anova_stats] = anovan(p_g_ratio_use(~isnan(p_g_ratio_use)), ...
    {sessions(~isnan(p_g_ratio_use))},'display','off');


[r p] = corrcoef(data_all{1},bins_all{1})
[r p] = corrcoef(data_all{2},bins_all{2})
[r p] = corrcoef(data_all{3},bins_all{3})


%% Supplementary 3: classification grid for all animals / cells

use_animals = 1:7;

figure;
for curr_animal = use_animals   
    
    curr_gad_class = horzcat(gad_class_all{curr_animal,:});
    day_center = nanmean(curr_gad_class.*meshgrid(1:size(curr_gad_class,2), ...
        1:size(curr_gad_class,1)),2);
    gad_class_cells = find(any(curr_gad_class,2));
    [asdf sort_idx] = sort(day_center(gad_class_cells));
    
    subplot(8,2,curr_animal*2-1);
    imagesc(curr_gad_class(gad_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Gad')
    end
    
    curr_pyr_class = horzcat(pyr_class_all{curr_animal,:});
    day_center = nanmean(curr_pyr_class.*meshgrid(1:size(curr_pyr_class,2), ...
        1:size(curr_pyr_class,1)),2);
    pyr_class_cells = find(any(curr_pyr_class,2));
    [asdf sort_idx] = sort(day_center(pyr_class_cells));
    
    subplot(8,2,curr_animal*2);
    imagesc(curr_pyr_class(pyr_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Pyr')
    end

end

%% Supplementary 4: fraction of cells active during movement

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity with lever 
pyr_move_active = cell(7,15);
gad_move_active = cell(7,15);
movement_lengths = cell(7,15);
act_trials_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        curr_reward_moves = move_epoch_frames{curr_animal,curr_session}( ...
            cued_rewarded_movements{curr_animal,curr_session});
      
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x) [leeway ...
            leeway],curr_reward_moves,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_reward_moves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        
        
%         %%% instead of using actual movements, to use 3s after
%         
%         % look at 3 seconds following movement onset (to go with lever)
%         move_time_framespan = 28*3;
%         move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
%             move_epoch_frames{curr_animal,curr_session},'uni',false);
%         
%         curr_rewardmoves_timed = move_time_frames( ...
%             cued_rewarded_movements{curr_animal,curr_session});
% 
%         % add leeway to movements (fixed frames)
%         leeway = 5;
%         leeway_frames = cellfun(@(x) [leeway ...
%             leeway],curr_rewardmoves_timed,'uni',false);
%         curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
%             curr_rewardmoves_timed,leeway_frames,'uni',false);
%         eliminate_movements = cellfun(@(x) ...
%             any(x < 1 | x > length(movement_trace_all ...
%             {curr_animal,curr_session})),curr_rewardmoves_leeway);
%         curr_rewardmoves_leeway(eliminate_movements) = [];

        
        %%%
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_rewardmoves_leeway);
        
    end

    curr_animal
    
end
use_animals = [1:7];
frac_cells_trial_active = cellfun(@nanmean,pyr_move_active(use_animals,1:14),'uni',false);
median_frac_cells_trial_active = cellfun(@nanmedian,frac_cells_trial_active);
figure; 
errorbar(nanmean(median_frac_cells_trial_active), nanstd( ...
    median_frac_cells_trial_active)./sqrt(sum(~isnan( ...
    median_frac_cells_trial_active))),'k','linewidth',2)
ylim([0 0.08]);
set(gca,'XTick',2:2:14)
xlabel('Session')
ylabel('Fraction of active excitatory neurons per movement')


% statistics
g1 = repmat(1:14,[14 1]);
[p anova_stats] = anovan(median_frac_cells_trial_active(~isnan(median_frac_cells_trial_active)), ...
    {g1(~isnan(median_frac_cells_trial_active))},'display','off');
[r p] = corrcoef(g1(~isnan(median_frac_cells_trial_active)), ...
    median_frac_cells_trial_active(~isnan(median_frac_cells_trial_active)));

%% Supplementary 5: changing for time preference
time_pref = cell(7,15);

for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        movetime_trace = nan(size(movement_trace_all{curr_animal,curr_session}));
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            pyr_activity_all{curr_animal,curr_session}(:,x), ...
            curr_rewardmoves_leeway,'uni',false);
        
        curr_activity_cat = cat(3,curr_activity{:});
        [c ci] = max(curr_activity_cat,[],2);
        curr_onset_timing = permute((c > 0).*ci,[1,3,2]);
        curr_onset_timing(~curr_onset_timing) = NaN;
        %time_pref{curr_animal,curr_session} = nanmean(curr_onset_timing,2);
        
        %time_pref{curr_animal,curr_session} = ...
        %    curr_onset_timing(pyr_class_all{curr_animal,curr_session},:);
        %temp_timing = curr_onset_timing(pyr_class_all{curr_animal,curr_session},:);
        %time_pref{curr_animal,curr_session} = temp_timing(~isnan(temp_timing));  
        
        % to get grid of all onsets       
        curr_onset_act = nan(size(curr_activity_cat));
        [x_grid y_grid z_grid] = meshgrid(1:size(c,2),1:size(c,1),1:size(c,3));
        onset_idx = sub2ind(size(curr_onset_act),y_grid(:),ci(:),z_grid(:));
        curr_onset_act(onset_idx(c(:) ~= 0)) = 1;
        time_pref{curr_animal,curr_session} = ...
            curr_onset_act(pyr_class_all{curr_animal,curr_session},:,:);       
                
    end
    disp(curr_animal);
end

% using activity onsets, find distribution
use_animals = 1:7;
%day_combine = {[1 2 3] [4 5 6] [7 8 9] [10 11 12] [13 14]};
day_combine = {1:3 4:10 11:14};
cat_act_onset = cell(length(day_combine),1);
total_onset_all = cell(length(day_combine));
col = [linspace(0,1,length(day_combine))' zeros(length(day_combine),1) ...
    zeros(length(day_combine),1)];
figure; hold on;
for i = 1:length(day_combine)
    curr_data = cellfun(@(x) permute(x,[2 3 1]), ...
        time_pref(use_animals,day_combine{i}),'uni',false);
    curr_data_reshape = cellfun(@(x) reshape(x,size(x,1),[]),curr_data,'uni',false);    
    cat_act_onset{i} = horzcat(curr_data_reshape{:}); 
    
    total_onset = nansum(cat_act_onset{i}');
    plot(linspace(0,length(total_onset)/28,length(total_onset))-5/28, ...
        cumsum(total_onset)/sum(total_onset),'color',col(i,:),'linewidth',2)
    
    total_onset_all{i} = total_onset;
end
xlabel('Time (s)')
ylabel('Fraction of activity onsets')
legend(cellfun(@(x) num2str(x),sessions_combine,'uni',false))

[h p] = kstest2(total_onset_all{1},total_onset_all{2})
[h p] = kstest2(total_onset_all{1},total_onset_all{3})
[h p] = kstest2(total_onset_all{2},total_onset_all{3})

% statistics
% get all onsets
[onset_max onset_times] = cellfun(@(x) max(x,[],1),cat_act_onset,'uni',false);
early_use_onsets = ~isnan(onset_max{1});
late_use_onsets = ~isnan(onset_max{5});
[h p] = kstest2(onset_times{1}(early_use_onsets),onset_times{5}(late_use_onsets));


%% TESTING: ranking from timing

% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end
% make matricies of first onsets
moveonset_aligned_pyr_first_onset = cell(size(moveonset_aligned_pyr_onset));
for i = 1:8
    for j = 1:15
        [c first_peak] = max(moveonset_aligned_pyr_onset{i,j}>0,[],2);
        x = zeros(size(moveonset_aligned_pyr_onset{i,j}));
        for curr_cell = 1:size(x,3)
            for trial = 1:size(first_peak,1)
                x(trial,first_peak(trial,1,curr_cell),curr_cell) = c(trial,1,curr_cell);
            end
        end
        moveonset_aligned_pyr_first_onset{i,j} = x;
    end
end


all_rank_increase_order = nan(size(moveonset_aligned_pyr_first_onset));
for curr_animalday = animaldays
    curr_activity = permute(moveonset_aligned_pyr_first_onset ...
        {curr_animalday},[3 2 1]);
    curr_mean_activity = nanmean(curr_activity,3);
    timing_matrix = repmat(1:size(curr_mean_activity,2), ...
        size(curr_mean_activity,1),1);
    curr_mean_timing = nansum(curr_mean_activity.*timing_matrix,2)./ ...
        nansum(curr_mean_activity,2);
    [asdf curr_ranks] = sort(curr_mean_timing);
    
    trial_ranks = repmat(curr_ranks,[1, ...
        size(curr_activity,2),size(curr_activity,3)]).*curr_activity; 
    
    % not sure what to do about same timing at the moment, for now just use
    % the one with the max rank
    trial_ranks_max = permute(max(trial_ranks),[3 2 1]);
    trial_rank_vectors = arrayfun(@(x) trial_ranks_max(x, ...
        trial_ranks_max(x,:) > 0),1:size(trial_ranks_max,1),'uni',false);
    
    rank_increase_order = nan(size(trial_rank_vectors));
    for curr_trial = 1:length(trial_rank_vectors)
        num_cells = length(trial_rank_vectors{curr_trial});
        % for now, timing's sake: ignore trials over 10 cells active
        if num_cells > 10
            continue
        end
        for i = length(trial_rank_vectors{curr_trial}):-1:1
            rank_perms = any(all(diff(nchoosek(trial_rank_vectors{curr_trial},i),[],2) > 0,2));
            if rank_perms
                rank_increase_order(curr_trial) = i;
                break
            end
            %disp([num2str(i) '/' num2str(length(trial_rank_vectors{curr_trial}))]);
        end
        %disp(curr_trial);
    end
    
    all_rank_increase_order(curr_animalday) = nanmean(rank_increase_order);
    disp(['curr_animalday ' num2str(curr_animalday)]);
end


num_shuff = 100;
all_rank_increase_order_shuff = cell(size(moveonset_aligned_pyr_first_onset));
for curr_shuff = 1:num_shuff
    for curr_animalday = animaldays
        
        curr_activity = permute(moveonset_aligned_pyr_first_onset ...
            {curr_animalday},[3 2 1]);
        
        % shuffle activity of active trials
        curr_activity_shuff = zeros(size(curr_activity));
        for curr_cell = 1:size(curr_activity,1)
            active_trials = find(any(curr_activity(curr_cell,:,:),2));
            curr_activity_shuff(curr_cell,:,active_trials) = ...
                curr_activity(curr_cell,:, ...
                active_trials(randperm(length(active_trials))));
        end
        
        curr_mean_activity = nanmean(curr_activity_shuff,3);
        timing_matrix = repmat(1:size(curr_mean_activity,2), ...
            size(curr_mean_activity,1),1);
        curr_mean_timing = nansum(curr_mean_activity.*timing_matrix,2)./ ...
            nansum(curr_mean_activity,2);
        [asdf curr_ranks] = sort(curr_mean_timing);
        
        trial_ranks = repmat(curr_ranks,[1, ...
            size(curr_activity_shuff,2),size(curr_activity_shuff,3)]).*curr_activity_shuff;
        
        % not sure what to do about same timing at the moment, for now just use
        % the one with the max rank
        trial_ranks_max = permute(max(trial_ranks),[3 2 1]);
        trial_rank_vectors = arrayfun(@(x) trial_ranks_max(x, ...
            trial_ranks_max(x,:) > 0),1:size(trial_ranks_max,1),'uni',false);
        
        rank_increase_order = nan(size(trial_rank_vectors));
        for curr_trial = 1:length(trial_rank_vectors)
            num_cells = length(trial_rank_vectors{curr_trial});
            if num_cells > 10
                continue
            end
            for i = length(trial_rank_vectors{curr_trial}):-1:1
                rank_perms = any(all(diff(nchoosek( ...
                    trial_rank_vectors{curr_trial},i),[],2) > 0,2));
                if rank_perms
                    rank_increase_order(curr_trial) = i;
                    break
                end
                %disp([num2str(i) '/' num2str(length(trial_rank_vectors{curr_trial}))]);
            end
            %disp(curr_trial);
        end
        
        all_rank_increase_order_shuff{curr_animalday}(curr_shuff) = ...
            nanmean(rank_increase_order);
    end
    disp(['curr_shuff ' num2str(curr_shuff)]);
end



%% TESTING: shuffle-subtracted correlation

% CORRELATION 


animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

%%%%%%%% Trying out correlation types


%%% cell-trial mean-subtracted correlation
moveonset_aligned_pyr_meansub_trial = cellfun(@(x) reshape(permute( ...
    bsxfun(@minus,x,nanmean(x,2)),[2 3 1]),[],size(x,1)), ...
    moveonset_aligned_pyr,'uni',false);

moveonset_aligned_pyr_meansub_trial_nonan = ...
    moveonset_aligned_pyr_meansub_trial;
for i = find(cellfun(@(x) ~isempty(x),moveonset_aligned_pyr_meansub_trial))'
    moveonset_aligned_pyr_meansub_trial_nonan{i}(isnan( ...
        moveonset_aligned_pyr_meansub_trial_nonan{i})) = 0;
end

meansub_trial_corr = cellfun(@(x) AP_itril(corrcoef(x),-1), ...
    moveonset_aligned_pyr_meansub_trial_nonan,'uni',false);

use_animals = 1:7;

meansub_trial_corr_mean = cellfun(@nanmean,meansub_trial_corr(use_animals,1:14));
meansub_trial_corr_m = nanmean(meansub_trial_corr_mean);
meansub_trial_corr_s = nanstd(meansub_trial_corr_mean)./ ...
    sqrt(sum(~isnan(meansub_trial_corr_mean)));
figure;errorbar(meansub_trial_corr_m,meansub_trial_corr_s,'k','linewidth',2)

meansub_trial_corr_cat = cell(14,1);
for i = 1:14
   meansub_trial_corr_cat{i} = vertcat(meansub_trial_corr{use_animals,i});
end
meansub_trial_corr_m = cellfun(@nanmean,meansub_trial_corr_cat);
meansub_trial_corr_s = cellfun(@(x) nanstd(x)./ ...
    sqrt(sum(~isnan(x))),meansub_trial_corr_cat);
figure;errorbar(meansub_trial_corr_m,meansub_trial_corr_s,'k','linewidth',2)

%%%%%%%%

%%% get trial-average correlation

mean_aligned_activity = cellfun(@(x) reshape(permute(nanmean(x),[2 3 1]),[],1), ...
    moveonset_aligned_pyr,'uni',false);
temporal_act_avg_corr = cell(7,15);
temporal_act_trial_corr = cell(7,15);
for i = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
    for j = sessions
        curr_data = reshape(permute(moveonset_aligned_pyr{i,j},[2 3 1]), ...
            [],size(moveonset_aligned_pyr{i,j},1));           
        % ignore nan values for correlation
        curr_data_mean = repmat(nanmean(curr_data),size(curr_data,1),1);
        curr_data(isnan(curr_data)) = curr_data_mean(isnan(curr_data));
        curr_template = mean_aligned_activity{i,j};
        curr_template(isnan(curr_template)) = nanmean(curr_template);
        
        curr_corr = corrcoef([curr_template curr_data]);
        temporal_act_avg_corr{i,j} = curr_corr(1,2:end);
        temporal_act_trial_corr{i,j} = AP_itril(curr_corr(2:end,2:end),-1);
    end
end

% shuffle time to subtract off non-temporal component
num_rep = 10000;
use_animals = 1:7;
temporal_act_avg_corr_shuff_all = cell(8,15,num_rep);
temporal_act_trial_corr_shuff_all = cell(8,15,num_rep);
for curr_rep = 1:num_rep
    temporal_act_avg_corr_shuff = cell(7,15);
    temporal_act_trial_corr_shuff = cell(7,15);
    for i = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
        for j = sessions
            
            % randomly circle shift each cell/trial
            curr_data = permute(moveonset_aligned_pyr{i,j},[2 3 1]);
            curr_data_shifted = nan(size(curr_data));
            
            [m n p] = size(curr_data);
            
            for curr_p = 1:p
                D = randi([0, m - 1], [1, n]);       
                for ii = (1 : n)
                    curr_data_shifted(:,ii,curr_p) = ...
                        [curr_data((m - D(ii) + 1 : m),ii,curr_p); ...
                        curr_data((1 : m - D(ii) ),ii,curr_p)];
                end
            end
            
            curr_data = reshape(curr_data_shifted, ...
                [],size(moveonset_aligned_pyr{i,j},1));
            % ignore nan values for correlation
            curr_data_mean = repmat(nanmean(curr_data),size(curr_data,1),1);
            curr_data(isnan(curr_data)) = curr_data_mean(isnan(curr_data));
            curr_template = mean_aligned_activity{i,j};
            curr_template(isnan(curr_template)) = nanmean(curr_template);
            
            curr_corr = corrcoef([curr_template curr_data]);
            temporal_act_avg_corr_shuff{i,j} = curr_corr(1,2:end);
            temporal_act_trial_corr_shuff{i,j} = AP_itril(curr_corr(2:end,2:end),-1);
        end
        i
    end
    temporal_act_avg_corr_shuff_all(:,:,curr_rep) = ...
        temporal_act_avg_corr_shuff;
    temporal_act_trial_corr_shuff_all(:,:,curr_rep) = ...
        temporal_act_trial_corr_shuff;
    curr_rep
end

% trial-trial correlation
temporal_act_trial_corr_shuff_all_cat = cell(7,15);
for i = 1:8
    for j = 1:15
        temporal_act_trial_corr_shuff_all_cat{i,j} = ...
            horzcat(temporal_act_trial_corr_shuff_all{i,j,:});
    end
end

act_trial_shuff_sub = cellfun(@(x,y) nanmean(x,2)-nanmean(y,2), ...
    temporal_act_trial_corr, temporal_act_trial_corr_shuff_all_cat,'uni',false);
act_trial_shuff_sub_cat = cell(14,1);
use_animals = 1:7;
for i = 1:14
   act_trial_shuff_sub_cat{i} = vertcat(act_trial_shuff_sub{use_animals,i}); 
end
% get total shuff subtracted mean
act_trial_shuff_sub_m = cellfun(@nanmean,act_trial_shuff_sub);
figure;errorbar(nanmean(act_trial_shuff_sub_m), ...
    nanstd(act_trial_shuff_sub_m)./sqrt(sum(~isnan(act_trial_shuff_sub_m))),'k','linewidth',2)
title('trial-trial correlation - shifted subtract')
figure;errorbar(cellfun(@nanmean,act_trial_shuff_sub_cat), ...
    cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),act_trial_shuff_sub_cat), ...
    'k','linewidth',2);
title('trial-trial correlation - shifted subtract')


% trial-average correlation
temporal_act_avg_corr_shuff_all_cat = cell(7,15);
for i = 1:8
    for j = 1:15
        temporal_act_avg_corr_shuff_all_cat{i,j} = ...
            vertcat(temporal_act_avg_corr_shuff_all{i,j,:});
    end
end

act_avg_shuff_sub = cellfun(@(x,y) nanmean(x,1)-nanmean(y,1), ...
    temporal_act_avg_corr, temporal_act_avg_corr_shuff_all_cat,'uni',false);
act_avg_shuff_sub_cat = cell(14,1);
use_animals = 1:7;
for i = 1:14
   act_avg_shuff_sub_cat{i} = horzcat(act_avg_shuff_sub{use_animals,i}); 
end
% get total shuff subtracted mean
act_avg_shuff_sub_m = cellfun(@nanmean,act_avg_shuff_sub);
figure;errorbar(nanmean(act_avg_shuff_sub_m), ...
    nanstd(act_avg_shuff_sub_m)./sqrt(sum(~isnan(act_avg_shuff_sub_m))),'k','linewidth',2)
title('trial-average correlation')
figure;errorbar(cellfun(@nanmean,act_avg_shuff_sub_cat), ...
    cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),act_avg_shuff_sub_cat), ...
    'k','linewidth',2);
title('trial-average correlation')



%% TESTING: trial-by-trial rank correlation (FIGURE 3j) 

% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_pyr_onset = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        cued_rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 5;fs_f = 3*28;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_pyr_onsets = cellfun(@(x) pyr_activity_onset(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_pyr_onset{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr_onsets{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

mean_moveonset_pyr = cellfun(@(x,y) permute(nanmean(x(:,:,y),1),[3 2 1]), ...
    moveonset_aligned_pyr,pyr_class_all,'uni',false);


use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);

pyr_trial_rank_corr = cell(size(use_animaldays));

num_cells = cell(size(use_animaldays));
all_reliability = cell(size(use_animaldays));
pyr_move_active_relative = cell(size(use_animaldays));
for i = find(use_animaldays)'
    
    curr_animal = mod(i,8);
    if curr_animal == 0;
       curr_animal = 8; 
    end
    %use_cells = pyr_class_all{i} & ...
    %    (move_onset_cells_all{i}' &...
    %    ~reward_cells_all{i}');
    %use_cells = pyr_class_all{i} & ...
    %    ~any(horzcat(move_onset_cells_all{curr_animal,:}),2)' &...
    %    any(horzcat(reward_cells_all{curr_animal,:}),2)';
    %use_cells = pyr_class_all{i} & ...
    %    any(horzcat(reward_cells_all{curr_animal,:}),2)';
    use_cells = pyr_class_all{i};
    use_cells_idx = find(use_cells);
    
    curr_rewardmoves = move_epoch_frames{i}( ...
        cued_rewarded_movements{i});   
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{i}) - x(1) + 1,curr_rewardmoves,'uni',false); 
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);    
     % add leeway to movements
    leeway = 0.25;
    %leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
    %    round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
    leeway_frames = cellfun(@(x,y) [10 ...
        10],curr_rewardmoves,curr_rewardframes,'uni',false);
    curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
        curr_rewardmoves,leeway_frames,'uni',false);
    eliminate_movements = cellfun(@(x) ...
        any(x < 1 | x > length(movement_trace_all{i})),curr_rewardmoves_leeway);
    curr_rewardmoves_leeway(eliminate_movements) = [];
    
    [curr_activity curr_time] = cellfun(@(x) ...
        max(pyr_activity_all{i}(use_cells,x)>0,[],2), ...
        curr_rewardmoves_leeway,'uni',false);
    
    curr_movement_lengths = cellfun(@length,...
        move_epoch_frames{i}(cued_rewarded_movements{i}));
    
    pyr_move_active{i} = cell2mat(curr_activity);
    pyr_move_time{i} = cell2mat(curr_time);
    
    % get timings/rank of mean and trial-trial, correlate
    [asdf mean_max_time] = max(mean_moveonset_pyr{i},[],2);
    all_rank_template = tiedrank(mean_max_time);
    
    [trial_max trial_max_time] = arrayfun(@(x) ...
        max(moveonset_aligned_pyr{i}(x,:,pyr_class_all{i}),[],2),1:size( ...
        moveonset_aligned_pyr{i},1),'uni',false);
    
    trial_rank = cellfun(@(x,y) tiedrank(permute(x(:,:,y > 0),[3 1 2])), ...
        trial_max_time,trial_max,'uni',false);
    
    trial_mean_rank = cellfun(@(x) all_rank_template(x > 0), trial_max,'uni',false);
      
    trial_rank_corr = cellfun(@(x,y) corrcoef(x,y), ...
        trial_rank,trial_mean_rank,'uni',false);
    
    % only look at trials with at least 3 cells active
    pyr_trial_rank_corr{i} = cellfun(@(x) x(2), trial_rank_corr( ...
        cellfun(@length,trial_rank) > 3));

end

pyr_corr_mean = cellfun(@nanmedian,pyr_trial_rank_corr);
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
use_animals = 1:7;
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
data_cat = cell(size(day_combine));
for i = 1:length(day_combine)
    %curr_data = horzcat(pyr_trial_rank_corr{use_animals,day_combine{i}});
    curr_data = vertcat(pyr_corr_mean(use_animals,day_combine{i}));
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
    data_cat{i} = curr_data;
end
figure;
errorbar(data_mean,data_sem,'k','linewidth',2);
ylabel('Rank correlation')
xlabel('Session')

% statistics
session_x = cellfun(@(x,y) x*ones(size(y)),num2cell(1:length(day_combine)),data_cat,'uni',false);
session_x_cat = vertcat(session_x{:});
data_cat_cat = vertcat(data_cat{:});
[r p] = corrcoef(session_x_cat(~isnan(data_cat_cat)),data_cat_cat(~isnan(data_cat_cat)));

%% TESTING: Movement vs. activity similarity within session

% get activity with lever 
movement_activity = cell(7,15);
act_trials_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            cued_rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% get relevant lever traces
use_samples = 10001:100:40001;
lever_trials_used = cell(7,15);
lever_used = cell(7,15);
for curr_animal = 1:7
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);

    % extracting all lever
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        % try using cued and uncued trials
        %curr_cued = true(size(curr_cued));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        total_curr_cat_lever{session} = temp(use_trials,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    lever_used(curr_animal,sessions) = total_curr_cat_lever;
    % creating template
    for session = sessions(end-2:end);
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        % try using cued and uncued trials
        %curr_cued = true(size(curr_cued));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        curr_cat_lever{session} = temp(use_trials,use_samples);
    end
    disp(curr_animal)
end


% use only trials with both activity and lever traces
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);

movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);


movement_activity_crop_cc = cellfun(@(x) AP_itril(corrcoef(x),-1),movement_activity_crop,'uni',false);
lever_used_crop_cc = cellfun(@(x) AP_itril(corrcoef(x'),-1),lever_used_crop,'uni',false);

corr_bins = [-1:0.1:1-0.1 Inf];
m = nan(length(corr_bins),14,7);
for curr_animal = 1:7;    
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        [n bins] = histc(lever_used_crop_cc{curr_animal,curr_session},corr_bins);
        
        m(unique(bins),curr_session,curr_animal) = ...
            grpstats(movement_activity_crop_cc{curr_animal,curr_session},bins);
        s(unique(bins),curr_session,curr_animal) = ...
            grpstats(movement_activity_crop_cc{curr_animal,curr_session},bins,'sem');
    end
    disp(curr_animal)
end

use_animals = 1:7;
figure; hold on;
session_combine = {1:3 4:8 9:14};
col = summer(length(session_combine));
slope_all = nan(length(session_combine),1);
curr_data_all = cell(size(session_combine));
for i = 1:length(session_combine);
   curr_data = reshape(m(:,session_combine{i},use_animals), ...
       length(corr_bins),[]);
   curr_data_all{i} = curr_data;
   %plot(nanmean(curr_data,2),'color',col(i,:),'linewidth',2);
   errorbar(nanmean(curr_data,2),nanstd(curr_data,[],2)./ ...
       sqrt(sum(~isnan(curr_data),2)),'color',col(i,:),'linewidth',2);
end
set(gca,'XTick',1:5:21);
set(gca,'XTickLabel',[corr_bins(1:5:end-5) 1]);
xlabel('Pairwise lever correlation');
ylabel('Pairwise activity correlation');
title('Activity within similar movements');
legend(cellfun(@num2str,session_combine,'uni',false),'location','nw')

% statistics
curr_data_cat = horzcat(curr_data_all{:});
g1 = repmat([1:21]',1,size(curr_data_cat,2));
g2_split = cellfun(@(x,y) x*ones(size(y)),num2cell(1:length(session_combine)),curr_data_all,'uni',false);
g2 = horzcat(g2_split{:});
[p anova_stats] = anovan(curr_data_cat(~isnan(curr_data_cat)), ...
    {g1(~isnan(curr_data_cat)) g2(~isnan(curr_data_cat))},'display','off');

curr_compare = [1 3];
curr_data_cat = horzcat(curr_data_all{curr_compare});
g3 = horzcat(g2_split{curr_compare});
[p anova_stats] = anovan(curr_data_cat(~isnan(curr_data_cat)), ...
    {g3(~isnan(curr_data_cat))},'display','off');

curr_corr = 3;
g4 = cellfun(@(x) repmat([1:size(x,1)]',1,size(x,2)),curr_data_all,'uni',false);
[r p] = corrcoef(curr_data_all{curr_corr}(~isnan(curr_data_all{curr_corr})), ...
    g4{curr_corr}(~isnan(curr_data_all{curr_corr})));








