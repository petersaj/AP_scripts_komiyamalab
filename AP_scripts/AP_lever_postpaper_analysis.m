%% PUT IN MEMORY: all lever data (including Simon's animals)

animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
lever_param_all = cell(size(animals));

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days (new way: go by bgROI files)
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_dir = dir([analysis_path filesep '*bgROI_analysis.mat']);
    analysis_names = {analysis_dir.name};
    days = sort(cellfun(@(x) x(1:6),analysis_names,'UniformOutput',false));
    
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
        
        % get behavior files/data
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        bhv_file = dir([data_path filesep animal '_dispatcher' ...
            filesep 'data_@lever2p_experimenter_' animal '_' day '*.mat']);
        bhv_filename = bhv_file.name;
        
        % get xsg files
        xsg_folder = [data_path filesep animal '_xsg' filesep day '_' animal '_xsg'];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = sort({xsg_dir.name});
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
            load([data_path filesep animal '_dispatcher' filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        % double check the dates are correct
        if datenum(saved.SavingSection_SaveTime(1:11)) ~= datenum([day(3:4) '/' day(5:6) '/' day(1:2)])
            warning([animal ' ' day ' bhv dates not correct']);
            keyboard;
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
    
    % find days (new way: go by bgROI files)
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_dir = dir([analysis_path filesep '*bgROI_analysis.mat']);
    analysis_names = {analysis_dir.name};
    days = sort(cellfun(@(x) x(1:6),analysis_names,'UniformOutput',false));

    if strcmp(animal,'AP98') || strcmp(animal,'AP99') ||strcmp(animal,'AP100')
        days = days(8:end);
        sessions = 1:7;
    end
    
    for curr_day = sessions
        day = days{curr_day};
        
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
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        bhv_file = dir([data_path filesep animal '_dispatcher' ...
            filesep 'data_@lever2p_experimenter_' animal '_' day '*.mat']);
        bhv_filename = bhv_file.name;
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep animal '_dispatcher' filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        % double check the dates are correct
        if datenum(saved.SavingSection_SaveTime(1:11)) ~= datenum([day(3:4) '/' day(5:6) '/' day(1:2)])
            warning([animal ' ' day ' bhv dates not correct']);
            keyboard;
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
        xsg_folder = [data_path filesep animal '_xsg' filesep day '_' animal '_xsg'];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = sort({xsg_dir.name});
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


%% Get 3s data around movement

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
            pyr_activity_all{curr_animal,curr_session}(:,x), ...
            curr_rewardmoves_leeway,'uni',false);
  

        cat_activity = permute(cat(3,curr_activity{:}),[3 2 1]);
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

%% Get corresponding lever traces

% get relevant lever traces
use_samples = 5001:100:35001;
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

%% Crop the activity and lever traces so they're corresponding

% use only trials with both activity and lever traces
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);

movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);


%% Get PCA of lever traces

check_pc = 5;

lever_curr_zscore = cellfun(@(x) zscore(x'),lever_used_crop,'uni',false);

[coeff score latent] = cellfun(@princomp,lever_curr_zscore,'uni',false);

curr_pc_varexp = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);

x_pc_varexp = nan(size(curr_pc_varexp));
use = cellfun(@(x) length(x) > check_pc,curr_pc_varexp);
x_pc_varexp(use) = cellfun(@(x) x(check_pc),curr_pc_varexp(use));
% normalize each animal by zscore
%x_pc_varexp_minusmean = bsxfun(@minus,x_pc_varexp,nanmean(x_pc_varexp,2));
%x_pc_varexp_zscore = bsxfun(@times,x_pc_varexp_minusmean,1./nanstd(x_pc_varexp_minusmean,[],2));
%figure;errorbar(nanmean(x_pc_varexp_zscore(:,1:14)),nanstd(x_pc_varexp_zscore(:,1:14))./ ...
%    sqrt(sum(~isnan(x_pc_varexp_zscore(:,1:14)))),'k','linewidth',2);
figure;errorbar(nanmean(x_pc_varexp(:,1:14)),nanstd(x_pc_varexp(:,1:14))./ ...
    sqrt(sum(~isnan(x_pc_varexp(:,1:14)))),'k','linewidth',2);
title('Lever')
ylabel(['Normalized variance up to PC ' num2str(check_pc)]);
xlabel('Session')


%% Get PCA of mean/normalized data

check_pc = 5;

movement_activity_curr = cellfun(@(x) permute(nanmean(x,1),[2 3 1]), ...
    movement_activity,'uni',false);

movement_activity_curr_zscore = cellfun(@zscore,movement_activity_curr,'uni',false);

[coeff score latent] = cellfun(@princomp,movement_activity_curr_zscore,'uni',false);

curr_pc_varexp = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);


x_pc_varexp = nan(size(curr_pc_varexp));
use = cellfun(@(x) length(x) ~= 1,curr_pc_varexp);
x_pc_varexp(use) = cellfun(@(x) x(check_pc),curr_pc_varexp(use));
% normalize each animal by zscore
x_pc_varexp_minusmean = bsxfun(@minus,x_pc_varexp,nanmean(x_pc_varexp,2));
x_pc_varexp_zscore = bsxfun(@times,x_pc_varexp_minusmean,1./nanstd(x_pc_varexp_minusmean,[],2));
figure;errorbar(nanmean(x_pc_varexp_zscore(:,1:14)),nanstd(x_pc_varexp_zscore(:,1:14))./ ...
    sqrt(sum(~isnan(x_pc_varexp_zscore(:,1:14)))),'k','linewidth',2);

title('Mean activity')
ylabel(['Normalized variance up to PC ' num2str(check_pc)]);
xlabel('Session')


%% Get PCA of trial/normalized data

check_pc = 5;

movement_activity_curr = cellfun(@(x) reshape(permute(x,[2 1 3]),[],size(x,3)), ...
    movement_activity,'uni',false);

movement_activity_curr_zscore = cellfun(@zscore,movement_activity_curr,'uni',false);

[coeff score latent] = cellfun(@princomp,movement_activity_curr_zscore,'uni',false);

curr_pc_varexp = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);


x_pc_varexp = nan(size(curr_pc_varexp));
use = cellfun(@(x) length(x) ~= 1,curr_pc_varexp);
x_pc_varexp(use) = cellfun(@(x) x(check_pc),curr_pc_varexp(use));
% normalize each animal by zscore
x_pc_varexp_minusmean = bsxfun(@minus,x_pc_varexp,nanmean(x_pc_varexp,2));
x_pc_varexp_zscore = bsxfun(@times,x_pc_varexp_minusmean,1./nanstd(x_pc_varexp_minusmean,[],2));
figure;errorbar(nanmean(x_pc_varexp_zscore(:,1:14)),nanstd(x_pc_varexp_zscore(:,1:14))./ ...
    sqrt(sum(~isnan(x_pc_varexp_zscore(:,1:14)))),'k','linewidth',2);

title('All trial activity')
ylabel(['Normalized variance up to PC ' num2str(check_pc)]);
xlabel('Session')


%% Get PCA of trial collapse mean/normalized data

check_pc = 5;

movement_activity_curr = cellfun(@(x) reshape(permute(nanmean(x,2),[1 3 2]),[],size(x,3)), ...
    movement_activity,'uni',false);

movement_activity_curr_zscore = cellfun(@zscore,movement_activity_curr,'uni',false);

[coeff score latent] = cellfun(@princomp,movement_activity_curr_zscore,'uni',false);

curr_pc_varexp = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);


x_pc_varexp = nan(size(curr_pc_varexp));
use = cellfun(@(x) length(x) ~= 1,curr_pc_varexp);
x_pc_varexp(use) = cellfun(@(x) x(check_pc),curr_pc_varexp(use));
% normalize each animal by zscore
x_pc_varexp_minusmean = bsxfun(@minus,x_pc_varexp,nanmean(x_pc_varexp,2));
x_pc_varexp_zscore = bsxfun(@times,x_pc_varexp_minusmean,1./nanstd(x_pc_varexp_minusmean,[],2));
figure;errorbar(nanmean(x_pc_varexp_zscore(:,1:14)),nanstd(x_pc_varexp_zscore(:,1:14))./ ...
    sqrt(sum(~isnan(x_pc_varexp_zscore(:,1:14)))),'k','linewidth',2);

title('Trial time-averaged activity')
ylabel(['Normalized variance up to PC ' num2str(check_pc)]);
xlabel('Session')



%% Get similarity of PCs across days (within the first n PCs)
num_pcs = 5;
pc_similarity = nan(15,15,length(animals));
for curr_animal = 1:length(animals);
    
    curr_sessions = cellfun(@(x) ~isempty(x),score(curr_animal,:));
    curr_pcs = cellfun(@(x) x(:,1:num_pcs),score(curr_animal,curr_sessions),'uni',false);
    cat_pcs = horzcat(curr_pcs{:});
    pc_cc = corrcoef(cat_pcs);
    pc_cc_split = mat2cell(pc_cc,num_pcs*ones(length(curr_pcs),1),num_pcs*ones(1,length(curr_pcs)));
    curr_pc_similarity = cellfun(@(x) ...
        nanmean(max(max(abs(x),[],2),[],2)),pc_cc_split);
    pc_similarity(curr_sessions,curr_sessions,curr_animal) = curr_pc_similarity;
    
end
curr_diag = 1;
compare_pc_cc = arrayfun(@(x) diag(pc_similarity(:,:,x),curr_diag),1:length(animals),'uni',false);
cat_compare_pc_cc = horzcat(compare_pc_cc{:})';
figure;errorbar(nanmean(cat_compare_pc_cc(:,1:14-abs(curr_diag))), ...
    nanstd(cat_compare_pc_cc(:,1:14-abs(curr_diag)))./ ...
    sqrt(sum(~isnan(cat_compare_pc_cc(:,1:14-abs(curr_diag))))),'k','linewidth',2);

ylabel('Average PC match correlations')
xlabel('Sessions')
title(['Session compare: ' num2str(curr_diag)])


compare_pc_grid = nanmean(pc_similarity(1:14,1:14,:),3);
compare_pc_grid(logical(eye(size(compare_pc_grid)))) = NaN;
figure;imagesc(compare_pc_grid);colormap(gray);


% trying out plotting examples
curr_animal = 6;

curr_sessions = cellfun(@(x) ~isempty(x),score(curr_animal,:));
a = cellfun(@(x) x(:,1:num_pcs),score(curr_animal,curr_sessions),'uni',false);

% figure; hold off
% col = jet(length(a));
% for i = 1:length(a)
%     plot(a{i},'color',col(i,:));
%     waitforbuttonpress
% end


curr_sessions = cellfun(@(x) ~isempty(x),score(curr_animal,:));
a = cellfun(@(x) x(:,1:num_pcs),score(curr_animal,curr_sessions),'uni',false);
b = horzcat(a{:});
c = corrcoef(b);
c_split = mat2cell(c,num_pcs*ones(length(a),1),num_pcs*ones(1,length(a)));

[mm idx] = cellfun(@(x) max(abs(x),[],2),c_split(end,:),'uni',false);

figure;
for curr_pc = 1:num_pcs;
r = cellfun(@(x) x(curr_pc),idx,'uni',false);
r2 = cellfun(@(x,y) x(:,y),score(curr_animal,curr_sessions),r,'uni',false);
a = horzcat(r2{:});
curr_cc = corrcoef(a);
sign_flip = sign(curr_cc(end,:));
curr_compare_pcs = bsxfun(@times,a,sign_flip);
subplot(1,num_pcs,curr_pc);imagesc(curr_compare_pcs);colormap(gray);
end



%% Masamizu/Matsuzaki 2014: binary movement linear regression

% Max normalize pyr activity
% WARNING! This normalization turns all 0 rows into all nan rows
pyr_activity_all_norm = cellfun(@(df) bsxfun(@times,df,1./ ...
    (max(df,[],2))),pyr_activity_all,'uni',false);

% Build up linear regression in random 1/10 of timeseries
regression_parts = 10;
rand_timepoints = cellfun(@(x) randperm(size(x,2)),movement_trace_all,'uni',false);
rand_timepoints_split = cellfun(@(x) mat2cell(x,1, ...
    repmat(length(x)/regression_parts,regression_parts,1)),rand_timepoints,'uni',false);

regressed_movement = cellfun(@(x) nan(size(x)),movement_trace_all,'uni',false);
for curr_animal = 1:size(pyr_activity_all_norm,1);
    for curr_session = 1:size(pyr_activity_all_norm,2);
        for curr_regression = 1:regression_parts
            if ~isempty(pyr_activity_all_norm{curr_animal,curr_session})
                warning off
                weights = pyr_activity_all_norm{curr_animal,curr_session}(~all( ...
                    isnan(pyr_activity_all_norm{curr_animal,curr_session}),2), ...
                    horzcat(rand_timepoints_split{curr_animal,curr_session} ...
                    {setdiff(1:regression_parts,curr_regression)}))'\ ...
                    movement_trace_all{curr_animal,curr_session}( ...
                    horzcat(rand_timepoints_split{curr_animal,curr_session} ...
                    {setdiff(1:regression_parts,curr_regression)}))';
                warning on
                
                regressed_movement{curr_animal,curr_session} ...
                    (rand_timepoints_split{curr_animal,curr_session}{curr_regression}) = ...
                    pyr_activity_all_norm{curr_animal,curr_session}(... 
                    ~all(isnan(pyr_activity_all_norm{curr_animal,curr_session}),2), ...
                    horzcat(rand_timepoints_split{curr_animal,curr_session} ...
                    {curr_regression}))'*weights;
            end
        end
        disp(curr_session);
    end
    disp(['Animal: ' num2str(curr_animal)]);
end


regression_correlation_sq = cellfun(@(x,y) corrcoef(x,y), ...
    regressed_movement,movement_trace_all,'uni',false);
regression_correlation = cellfun(@(x) x(2), regression_correlation_sq);

early_corr = nanmean(regression_correlation(:,1:4),2);
late_corr = nanmean(regression_correlation(:,11:14),2);
figure; hold on
plot([early_corr late_corr]','color',[0.5 0.5 0.5]);
plot([nanmean(early_corr) nanmean(late_corr)],'k','linewidth',2);
p = signrank(early_corr,late_corr);
ylim([0 0.7]);
xlim([0 3]);
ylabel('Correlation')
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Early' 'Late'});
title(['Movement/nonmovement linear regression (p = ' num2str(p) ')'])
















