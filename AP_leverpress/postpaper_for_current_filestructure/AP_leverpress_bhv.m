%%  FOR POSTPAPER

% The data was restructure after the paper, this is to access the new
% structure (was restructured to send to Taro Toyoizumi, instead of having
% separate days in folders, bhv/xsg/data in separate folders, then
% separated by day



%% PUT IN MEMORY: all lever data (including Simon's animals)

%animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' ...
%    'AP78' 'AP79' 'SC027' 'SC028' 'SC029' ...
%    'AP91' 'AP92' 'AP93' 'AP94'};
animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' ...
    'AP78' 'AP79'};
lever_param_all = cell(size(animals));

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    
    bhv_path = [data_path filesep animal '_dispatcher'];
    xsg_path = [data_path filesep animal '_xsg'];
    
    bhv_all = dir(bhv_path);
    xsg_all = dir(xsg_path);
    
    days = cellfun(@(x) x(1:6),{xsg_all(3:end).name}, ...
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
            curr_day lever_param data_path lever_param_all ...
            bhv_path xsg_path bhv_all xsg_all
        
        
        day = num2str(days{curr_day});
        
        % get behavior file
        bhv_file = cellfun(@(x) ~isempty(strfind([bhv_path filesep x],day)),{bhv_all.name});
        bhv_filename = [bhv_path filesep bhv_all(bhv_file).name];
        
        % get xsg files
        xsg_folder = [xsg_path filesep xsg_all(curr_day+2).name];
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
            load(bhv_filename,'-MAT');
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



%% Figure 1: lever average-average correlation and trial-trial correlation

lever_downsamp = 100;
use_samples = 49:350; % in 10 ms groups
lever_used = cell(11,15);
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
lever_trial_corr_cat = cell(14,1);
use_animals = 1:length(lever_param_all);
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
    sessions(~isnan(lever_trial_corr_median_use)));


% plot avg-avg corrlation
lever_avg = cellfun(@nanmean,lever_used,'uni',false);
lever_avg_corr = nan(15,15,10);
for curr_animal = 1:10
   sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
   curr_lever = vertcat(lever_avg{curr_animal,sessions});
   lever_avg_corr(sessions,sessions,curr_animal) = corrcoef(curr_lever');  
end
use_animals = 1:length(lever_param_all);
figure;
imagesc(nanmean(lever_avg_corr(1:14,1:14,use_animals),3))
colormap(hot);caxis([0 1])
xlabel('Session')
ylabel('Session')
title('Mean lever press correlation')

% plot avg-avg next-day correlation
lever_avg_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_avg_corr(:,:,x),-1), ...
    use_animals,'uni',false))';
figure;errorbar(nanmean(lever_avg_corr_next_session),nanstd(...
    lever_avg_corr_next_session)./sqrt(sum(~isnan(lever_avg_corr_next_session))), ...
    'k','linewidth',2);

% plot avg trial-trial correlation across sessions
lever_trial_corr = cell(15,15,10);
for curr_animal = 1:10
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
use_animals = 1:length(lever_param_all);
lever_trial_corr_median = cellfun(@(x) nanmedian(x(:)),lever_trial_corr);
figure;imagesc(nanmean(lever_trial_corr_median(1:14,1:14,use_animals),3));colormap(hot);
title('Trial-trial correlation')

lever_trial_corr_next_session = cell2mat(arrayfun(@(x) diag(lever_trial_corr_median(1:14,1:14,x),-1), ...
    use_animals,'uni',false))';
figure;errorbar(nanmean(lever_trial_corr_next_session),nanstd(...
    lever_trial_corr_next_session)./sqrt(sum(~isnan(lever_trial_corr_next_session))), ...
    'k','linewidth',2);
title('Trial-trial next-session correlation')











