%% Copy all XSG/BHV from imaging folder to data folder (non-imaging expts)

animal = 'SC064';
days = {'131125' '131126' '131128' '131129' '131130' '131201' ...
    '131202' '131203' '131204' '131205' '131206' '131207' '131208'};

img_path = ['/usr/local/lab/Data/ImagingRig3/'];
data_path = ['/usr/local/lab/People/Andy/Data/' animal];

if ~exist(data_path,'dir')
    mkdir(data_path)
end


for curr_day = 1:length(days)
    
    day = days{curr_day};
    curr_source = [img_path filesep day filesep animal];
    curr_destination = [data_path filesep day];
    if ~exist(curr_destination,'dir')
        mkdir(curr_destination)
    end
    
    % bhv files
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = (dir_filenames(bhv_file == 1));
    
    for i = 1:length(bhv_filename)
        copyfile([curr_source filesep bhv_filename{i}], ...
            [curr_destination filesep bhv_filename{i}]);
    end
    
    % xsg files
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    dir_isfolders = cellfun(@(x) isdir([curr_source filesep x]), ...
        dir_filenames);
    xsg_regexp = regexp(dir_filenames(dir_isfolders),'\w\w\d\d\d\d');
    xsg_folders = cellfun(@(x) ~isempty(x),xsg_regexp);
    dir_folders = dir_filenames(dir_isfolders);
    curr_xsg_folders = dir_folders(xsg_folders);
    
    for i = 1:length(curr_xsg_folders)
        copyfile([curr_source filesep curr_xsg_folders{i}], ...
            [curr_destination filesep curr_xsg_folders{i}]);
    end
    
    
    disp(['Finished day ' num2str(curr_day) '/' num2str(length(days))])

end


%% PUT IN MEMORY: all lever data
animals = {'SC063'};
%animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' ...
%    'AP78' 'AP79' 'SC027' 'SC028' 'SC029' ...
%    'AP91' 'AP92' 'AP93' 'AP94'};
%animals = {'AP91' 'AP92' 'AP93' 'AP94'};
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
    
    lever_param.opto_trial= cell(size(days));
    lever_param.pre_cue_movement= cell(size(days));
    lever_param.response_movement = cell(size(days));
    lever_param.reward_reaction = cell(size(days));
    lever_param.reaction_time = cell(size(days));
    lever_param.movement_reward_time = cell(size(days));
    lever_param.cued_movement = cell(size(days));
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
    lever_param.global_reward_times = cell(size(days));
    
    
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
        
        % get which trials had light stim (can be 1-2 uncompleted trials at
        % end)
        if isfield(saved_history,'OptoSection_opto_trial');
            lever_param.opto_trial{curr_day} = cell2mat(saved_history.OptoSection_opto_trial(1:num_trials));
        end
        lever_param.pre_cue_movement{curr_day} = nan(num_trials,1);
        lever_param.response_movement{curr_day} = nan(num_trials,1);
        lever_param.reward_reaction{curr_day} = nan(num_trials,1);
        lever_param.reaction_time{curr_day} = nan(num_trials,1);
        lever_param.movement_reward_time{curr_day} = nan(num_trials,1);
        lever_param.cued_movement{curr_day} = false(num_trials,1);
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
        lever_param.global_reward_times{curr_day} = nan(num_trials,1);
        
        % Process
        
        % Only look at rewarded trials
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));

        % get times for rewards in those trials
        lever_param.global_reward_times{curr_day}(rewarded_trials) = cellfun(@(x)  ...
        x.states.reward(1), bhv(rewarded_trials));
        
        % Initializing saving raw traces around reward
        time_surround = 1;
        sample_surround = [-0.5*xsg_samplerate:4*xsg_samplerate];
        lever_param.reward_aligned_lever{curr_day} = ...
            nan(num_trials,length(sample_surround));
        use_samples = 1:length(sample_surround);

        
        % go through all xsg files, grab movement properties
        plots_on = false;
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

            end
            if plots_on
                figure;hold on;
                %h1 = subplot(2,1,1);
                %h2 = subplot(2,1,2);
                title('Green = cued, Red = unued, Cyan = opto')
                for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
                    curr_trial = curr_trial_list(curr_trial_idx,2);
                    %subplot(h1)
                    lever_force_resample_active = lever_force_resample;
                    lever_force_resample_active(~lever_active) = NaN;
                    lever_force_resample_inactive = lever_force_resample;
                    lever_force_resample_inactive(lever_active) = NaN;
                    
                    plot(lever_force_resample_active,'r');
                    plot(lever_force_resample_inactive,'k');
                    line([reward_start_sample(curr_trial_idx) reward_start_sample(curr_trial_idx)],ylim);
                    if lever_param.cued_movement{curr_day}(curr_trial)
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
                    else
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
                    end
                    if lever_param.opto_trial{curr_day}(curr_trial)
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim, ...
                            'color','w','linestyle','--');
                    end
%                     subplot(h2)
%                     line([reward_start_sample(curr_trial_idx) reward_start_sample(curr_trial_idx)],ylim);
%                     line([reward_stop_sample(curr_trial_idx) reward_stop_sample(curr_trial_idx)],ylim);
%                     if cued_movement
%                         line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
%                     else
%                         line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
%                     end
                end
                keyboard
            end
        end
        disp(['Finished day: ' day]);
    end
    disp('Finished all')
    
lever_param_all{curr_animal} = lever_param;

end
disp('Finished all animals')

%% Get correlation of levers
max_days = max(cellfun(@(x) length(x.cued_movement),lever_param_all));

lever_downsamp = 100;
use_samples = 49:350; % in 10 ms groups
lever_used = cell(length(lever_param_all),max_days);
for curr_animal = 1:length(lever_param_all)
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
lever_trial_corr_cat = cell(max_days,1);
for i = 1:length(lever_param_all)
   %lever_trial_corr_cat{i} = vertcat(lever_trial_corr{use_animals,i}); 
   lever_trial_corr_cat{i} = vertcat(lever_trial_corr_median(i,1)); 
end
lever_trial_corr_mean = cellfun(@nanmean,lever_trial_corr_cat); 
lever_trial_corr_sem = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), ...
    lever_trial_corr_cat);
figure;
errorbar(lever_trial_corr_mean,lever_trial_corr_sem,'k','linewidth',2)
xlabel('Session')
ylabel('Correlation')
title('Mean trial-trial correlation')

% plot avg-avg corrlation
lever_avg = cellfun(@nanmean,lever_used,'uni',false);
lever_avg_corr = nan(max_days,max_days,length(lever_param_all));
for curr_animal = 1:length(lever_param_all)
   sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
   curr_lever = vertcat(lever_avg{curr_animal,sessions});
   lever_avg_corr(sessions,sessions,curr_animal) = corrcoef(curr_lever');  
end
figure;
imagesc(nanmean(lever_avg_corr(1:max_days,1:max_days,:),3))
colormap(hot);caxis([0 1])
xlabel('Session')
ylabel('Session')
title('Mean lever press correlation')

% plot avg-avg next-day correlation
lever_avg_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_avg_corr(:,:,x),-1), ...
    1:length(lever_param_all),'uni',false))';
figure;errorbar(nanmean(lever_avg_corr_next_session),nanstd(...
    lever_avg_corr_next_session)./sqrt(sum(~isnan(lever_avg_corr_next_session))), ...
    'k','linewidth',2);

% plot avg trial-trial correlation across sessions
lever_trial_corr = cell(max_days,max_days,length(lever_param_all));
for curr_animal = 1:11
    sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
    for curr_session_1 = sessions
        for curr_session_2 = sessions           
           if curr_session_1 == curr_session_2
               curr_lever_cc = corrcoef(lever_used{curr_animal,curr_session_1}');
               lever_trial_corr{curr_session_1,curr_session_2,curr_animal} = ...
                   AP_itril(curr_lever_cc,-1);
               continue
           end
           curr_lever_cc = corrcoef([lever_used{curr_animal,curr_session_1}; ...
               lever_used{curr_animal,curr_session_2}]');
           lever_trial_corr{curr_session_1,curr_session_2,curr_animal} = ...
               curr_lever_cc(1:size(lever_used{curr_animal,curr_session_1},1), ...
               size(lever_used{curr_animal,curr_session_1},1)+1:end);
        end
    end
end
lever_trial_corr_median = cellfun(@(x) nanmedian(x(:)),lever_trial_corr);
figure;imagesc(nanmean(lever_trial_corr_median(1:max_days,1:max_days,:),3));colormap(hot);
title('Trial-trial correlation')

lever_trial_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_trial_corr_median(1:max_days,1:max_days,x),-1), ...
    1:length(lever_param_all),'uni',false))';
figure;errorbar(nanmean(lever_trial_corr_next_session),nanstd(...
    lever_trial_corr_next_session)./sqrt(sum(~isnan(lever_trial_corr_next_session))), ...
    'k','linewidth',2);
title('Trial-trial next-session correlation')

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


%% M1 lesion

normal_animals = 1:10;
lesion_animals = 11:14;

% Lever trial-trial correlation
max_days = max(cellfun(@(x) length(x.cued_movement),lever_param_all));

lever_downsamp = 100;
use_samples = 49:350; % in 10 ms groups
lever_used = cell(length(lever_param_all),max_days);
for curr_animal = 1:length(lever_param_all)
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
figure; hold on;
errorbar(nanmean(lever_trial_corr_median(normal_animals,1:14)), nanstd( ...
    lever_trial_corr_median(normal_animals,1:14))./sqrt(sum(~isnan( ...
    lever_trial_corr_median(normal_animals,1:14)))),'k','linewidth',2);
errorbar(nanmean(lever_trial_corr_median(lesion_animals,1:14)), nanstd( ...
    lever_trial_corr_median(lesion_animals,1:14))./sqrt(sum(~isnan( ...
    lever_trial_corr_median(lesion_animals,1:14)))),'r','linewidth',2);

% combine all normal/lesion correlations
lever_trial_corr_combine = cell(2,14);
for i = 1:14
    lever_trial_corr_combine{1,i} = vertcat(lever_trial_corr{normal_animals,i});
    lever_trial_corr_combine{2,i} = vertcat(lever_trial_corr{lesion_animals,i});
end
p = arrayfun(@(x) ranksum(lever_trial_corr_combine{1,x}, ...
    lever_trial_corr_combine{2,x}),1:14);
% Behavior times

% Reaction time (cued movements only)
curr_param = cell(length(lever_param_all),max_days);
for i = 1:length(lever_param_all);
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reaction_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_param(i,1:length(curr_data)) = curr_data;  
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
errorbar(nanmean(curr_param_median(normal_animals,1:14)), nanstd( ...
    curr_param_median(normal_animals,1:14))./sqrt(sum(~isnan( ...
    curr_param_median(normal_animals,1:14)))),'k','linewidth',2);
errorbar(nanmean(curr_param_median(lesion_animals,1:14)), nanstd( ...
    curr_param_median(lesion_animals,1:14))./sqrt(sum(~isnan( ...
    curr_param_median(lesion_animals,1:14)))),'r','linewidth',2);
title('Reaction time')

% Movement onset to reward (cued movements only)
curr_param = cell(length(lever_param_all),max_days);
for i = 1:length(lever_param_all)
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reaction_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_data2 = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reward_reaction, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_data3 = cellfun(@(x,y) y-x,curr_data,curr_data2,'uni',false);
    curr_param(i,1:length(curr_data3)) = curr_data3;  
end
curr_param_median = cellfun(@nanmedian,curr_param);
figure; hold on;
errorbar(nanmean(curr_param_median(normal_animals,1:14)), nanstd( ...
    curr_param_median(normal_animals,1:14))./sqrt(sum(~isnan( ...
    curr_param_median(normal_animals,1:14)))),'k','linewidth',2);
errorbar(nanmean(curr_param_median(lesion_animals,1:14)), nanstd( ...
    curr_param_median(lesion_animals,1:14))./sqrt(sum(~isnan( ...
    curr_param_median(lesion_animals,1:14)))),'r','linewidth',2);
title('Movement onset to reward time')

% Rate of rewards (only day 3+ b/c consistent times from then on)
rewarded_trials = cellfun(@(y) cellfun(@(x) ~isnan(x),y.global_reward_times, ...
    'uni',false), lever_param_all,'uni',false);
reward_fit = cellfun(@(y,z) cellfun(@(x,w) polyfit(x(w),cumsum(w(w)),1), ...
    y.global_reward_times(3:end),z(3:end),'uni',false),lever_param_all,rewarded_trials,'uni',false);
reward_rate = cellfun(@(y) cellfun(@(x) x(1),y),reward_fit,'uni',false);
reward_rate_matrix = nan(length(lever_param_all),max_days-2);
for curr_animal = 1:length(lever_param_all)
    sessions = length(reward_rate{curr_animal});
    reward_rate_matrix(curr_animal,1:sessions) = ...
        reward_rate{curr_animal};
end


%% PV-ChR2: comparison of light v non-light trials

% just get rewarded fraction
animals = {'SC063','SC064'};

rewarded_trials = cell(length(animals),2);
opto_trials = cell(length(animals),2);
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
        
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
        
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        
        rewarded_trials{curr_animal,curr_day} = ...
            cellfun(@(x) ~isempty(x.states.reward), bhv);
        opto_trials{curr_animal,curr_day} = ...
            cell2mat(saved_history.OptoSection_opto_trial(1:num_trials));
 
    end
    curr_animal;
end
rewarded_opto = cellfun(@(x,y) x(logical(y)),rewarded_trials,opto_trials,'uni',false);
rewarded_nonopto = cellfun(@(x,y) x(~logical(y)),rewarded_trials,opto_trials,'uni',false);

rewarded_opto_mean = cellfun(@nanmean,rewarded_opto);
rewarded_nonopto_mean = cellfun(@nanmean,rewarded_nonopto);
figure;
subplot(2,1,1);hold on;
plot(rewarded_nonopto_mean','k','linewidth',2);
plot(rewarded_opto_mean','color',[0.5 0.8 1],'linewidth',2)

subplot(2,1,2);
light_on = vertcat(rewarded_opto{:,8:12});
light_off = vertcat(rewarded_nonopto{:,8:12});
ctrl_light_on = vertcat(rewarded_opto{:,13:14});
ctrl_light_off = vertcat(rewarded_nonopto{:,13:14});
plot_mean = [nanmean(light_off(:)) nanmean(light_on(:)); ...
    nanmean(ctrl_light_off(:)) nanmean(ctrl_light_on(:))];
plot_sem = [nanmean(light_off(:))/sqrt(sum(~isnan(light_off(:)))) ...
    nanmean(light_on(:))/sqrt(sum(~isnan(light_on(:)))); ...
    nanmean(ctrl_light_off(:))/sqrt(sum(~isnan(ctrl_light_off(:))))...
    nanmean(ctrl_light_on(:))/sqrt(sum(~isnan(ctrl_light_on(:))))];
h = bar(plot_mean);
set(gca,'XTickLabel',{'Window clear' 'Window blocked'})
set(h(1),'FaceColor',[0 0 0])
set(h(2),'FaceColor',[0.5 0.8 1])
legend({'Light off' 'Light on'})
xlim([0 3]);
%hold on
%x1 = mean(reshape(unique(get(get(h(1),'children'),'xdata')),2,2));
%x2 = mean(reshape(unique(get(get(h(2),'children'),'xdata')),2,2));
%errorbar([x1' x2'],plot_mean,plot_sem,'.k','linewidth',2);

% % Look at lever params
% use_animals = 1:length(animals);
% max_days = max(cellfun(@(x) length(x.cued_movement),lever_param_all));
% 
% % current expt: combine days past 7
% opto_day_combine = structfun(@(x) [x(1:7) {vertcat(x{8:end})}],lever_param_all{1},'uni',false);
% lever_param_all{1} = opto_day_combine;
% 
% % Reaction time (cued movements only)
% curr_param_opto = cell(length(use_animals),max_days);
% for i = use_animals
%     curr_data = cellfun(@(x,y,z) x(y & z), ...
%         lever_param_all{i}.reaction_time, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_param_opto(i,1:length(curr_data)) = curr_data;  
% end
% 
% curr_param_nonopto = cell(length(use_animals),max_days);
% for i = use_animals
%     curr_data = cellfun(@(x,y,z) x(y & ~z), ...
%         lever_param_all{i}.reaction_time, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_param_nonopto(i,1:length(curr_data)) = curr_data;  
% end
% 
% %p = cellfun(@(x,y) ranksum(x,y),curr_param_opto,curr_param_nonopto)
% curr_param_nonopto_median = cellfun(@nanmedian,curr_param_nonopto);
% curr_param_opto_median = cellfun(@nanmedian,curr_param_opto);
% figure; hold on;
% errorbar(nanmean(curr_param_nonopto_median(use_animals,:),1), nanstd( ...
%     curr_param_nonopto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_nonopto_median(use_animals,:)))),'k','linewidth',2);
% errorbar(nanmean(curr_param_opto_median(use_animals,:),1), nanstd( ...
%     curr_param_opto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_opto_median(use_animals,:)))),'color',[0.5 0.8 1],'linewidth',2);
% title('Reaction time')
% ylabel('Time (ms)')
% xlabel('Sessions')
% 
% % Movement onset to reward (cued movements only)
% curr_param_opto = cell(length(use_animals),max_days);
% for i = use_animals
%     curr_data = cellfun(@(x,y,z) x(y & z), ...
%         lever_param_all{i}.reaction_time, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_data2 = cellfun(@(x,y,z) x(y & z), ...
%         lever_param_all{i}.reward_reaction, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_data3 = cellfun(@(x,y) y-x,curr_data,curr_data2,'uni',false);
%     curr_param_opto(i,1:length(curr_data)) = curr_data3;  
% end
% 
% curr_param_nonopto = cell(length(use_animals),max_days);
% for i = use_animals
%     curr_data = cellfun(@(x,y,z) x(y & ~z), ...
%         lever_param_all{i}.reaction_time, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_data2 = cellfun(@(x,y,z) x(y & ~z), ...
%         lever_param_all{i}.reward_reaction, ...
%         lever_param_all{i}.cued_movement, ...
%         lever_param_all{i}.opto_trial,'uni',false);
%     curr_data3 = cellfun(@(x,y) y-x,curr_data,curr_data2,'uni',false);
%     curr_param_nonopto(i,1:length(curr_data)) = curr_data3;
% end
% 
% %p = ranksum(curr_param_opto{1},curr_param_nonopto{1})
% curr_param_nonopto_median = cellfun(@nanmedian,curr_param_nonopto);
% curr_param_opto_median = cellfun(@nanmedian,curr_param_opto);
% figure; hold on;
% errorbar(nanmean(curr_param_nonopto_median(use_animals,:),1), nanstd( ...
%     curr_param_nonopto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_nonopto_median(use_animals,:)))),'k','linewidth',2);
% errorbar(nanmean(curr_param_opto_median(use_animals,:),1), nanstd( ...
%     curr_param_opto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_opto_median(use_animals,:)))),'color',[0.5 0.8 1],'linewidth',2);
% title('Movement onset to reward time')
% ylabel('Time (ms)')
% xlabel('Sessions')
% 
% % trial-trial correlation
% max_days = max(cellfun(@(x) length(x.cued_movement),lever_param_all));
% lever_downsamp = 100;
% use_samples = 49:350; % in 10 ms groups
% 
% lever_used_opto = cell(length(lever_param_all),max_days);
% for curr_animal = 1:length(lever_param_all)
%     sessions = find(cellfun(@(x) ~isempty(x), ...
%         lever_param_all{curr_animal}.reward_aligned_lever));
%     total_curr_cat_lever = cell(max(sessions),1);
%     curr_cat_lever = cell(length(sessions),1);
%     % extracting all lever
%     for session = sessions;
%         curr_cued = lever_param_all{curr_animal}.cued_movement{session};
%         curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session},1));
%         % check all trials
%         %curr_cued(:) = true;
%         
%         curr_opto = lever_param_all{curr_animal}.opto_trial{session};
%         curr_opto = curr_opto(1:size(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session},1));
%         
%         temp_lever = vertcat(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session}(curr_cued & curr_opto,:));   
%         use_trials = ~any(isnan(temp_lever),2) & ...
%             ~all(temp_lever == 0,2);
%         if ~isempty(temp_lever);
%             temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
%             total_curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);
%         end
%         
%     end
%     lever_used_opto(curr_animal,sessions) = total_curr_cat_lever(sessions);
%     
%     disp(curr_animal)
% end
% lever_corr_opto = cellfun(@(x) AP_itril(corrcoef(x'),-1),lever_used_opto,'uni',false);
% 
% lever_used_nonopto = cell(length(lever_param_all),max_days);
% for curr_animal = 1:length(lever_param_all)
%     sessions = find(cellfun(@(x) ~isempty(x), ...
%         lever_param_all{curr_animal}.reward_aligned_lever));
%     total_curr_cat_lever = cell(max(sessions),1);
%     curr_cat_lever = cell(length(sessions),1);
%     % extracting all lever
%     for session = sessions;
%         curr_cued = lever_param_all{curr_animal}.cued_movement{session};
%         curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session},1));
%         % check all trials
%         curr_cued(:) = true;
%         
%         curr_opto = lever_param_all{curr_animal}.opto_trial{session};
%         curr_opto = curr_opto(1:size(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session},1));
%         
%         temp_lever = vertcat(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session}(curr_cued & ~curr_opto,:));   
%         use_trials = ~any(isnan(temp_lever),2) & ...
%             ~all(temp_lever == 0,2);
%         if ~isempty(temp_lever);
%             temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
%             total_curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);
%         end       
% 
%     end
%     lever_used_nonopto(curr_animal,sessions) = total_curr_cat_lever(sessions);
%     
%     disp(curr_animal)
% end
% lever_corr_nonopto = cellfun(@(x) AP_itril(corrcoef(x'),-1),lever_used_nonopto,'uni',false);
% 
% %p = ranksum(lever_corr_opto{9},lever_corr_nonopto{9})
% curr_param_nonopto_median = cellfun(@nanmedian,lever_corr_nonopto);
% curr_param_opto_median = cellfun(@nanmedian,lever_corr_opto);
% figure; hold on;
% errorbar(nanmean(curr_param_nonopto_median(use_animals,:),1), nanstd( ...
%     curr_param_nonopto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_nonopto_median(use_animals,:)))),'k','linewidth',2);
% errorbar(nanmean(curr_param_opto_median(use_animals,:),1), nanstd( ...
%     curr_param_opto_median(use_animals,:),[],1)./sqrt(sum(~isnan( ...
%     curr_param_opto_median(use_animals,:)))),'color',[0.5 0.8 1],'linewidth',2);
% title('Within-day lever correlation')
% ylabel('Correlation')
% xlabel('Sessions')
% 
% % % At some point, could pull out lever releases by resonant freq (~35 hz)
% % [S F T P] = spectrogram(lever_force_resample,100,5,1000,1000,'yaxis');
% % a = 10*log10(P);
% % j = [find(F > 31,1) find(F > 35,1)];
% % b = nanmean(a(j(1):j(2),:));
% 
% % fraction of rewarded trials



%% Get correlation of levers
max_days = max(cellfun(@(x) length(x.cued_movement),lever_param_all));

lever_downsamp = 100;
use_samples = 49:350; % in 10 ms groups
lever_used = cell(length(lever_param_all),max_days);
for curr_animal = 1:length(lever_param_all)
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




%% (for iacuc test)
bhv = cell(4,1);

bhv{1} = load('/usr/local/lab/People/Andy/Data/AP87/130906/data_@lever2p_experimenter_AP87_130906b.mat');
bhv{2} = load('/usr/local/lab/People/Andy/Data/AP87/130907/data_@lever2p_experimenter_AP87_130907b.mat');
bhv{3} = load('/usr/local/lab/People/Andy/Data/AP87/130908/data_@lever2p_experimenter_AP87_130908b.mat');
bhv{4} = load('/usr/local/lab/People/Andy/Data/AP87/130925/data_@lever2p_experimenter_AP87_130925b.mat');

reward_trials = cell(4,1);
reward_times = cell(4,1);
for i = 1:4
    reward_trials{i} = find(cellfun(@(x) ~isempty( ...
        x.states.reward), bhv{i}.saved_history.ProtocolsSection_parsed_events));
    
    reward_times{i} = cellfun(@(x)  ...
        x.states.reward(1), bhv{i}.saved_history.ProtocolsSection_parsed_events(reward_trials{i}));
end





















