%% Quantify lever activity over days

% Lever paper

clear all

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
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
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);
    
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
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
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

        
        % Initializing saving raw traces around reward
        time_surround = 1;
        sample_surround = [-0.5*xsg_samplerate:4*xsg_samplerate];
        lever_param.reward_aligned_lever{curr_day} = ...
            nan(num_trials,length(sample_surround));
        use_samples = 1:length(sample_surround);

        
        % go through all xsg files, grab movement properties
        plots_on = 0;
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            % create movement velocity, active times
            % quantify movement before
            
            [lever_active lever_force_resample lever_force_smooth ...
                lever_velocity_envelope_smooth] = AP_parseLeverMovement(xsg_data);
            
            lever_stopstart = [0;diff(lever_active)];
            
            %             move_start_frames = find(lever_stopstart == 1);
            %             move_stop_frames = find(lever_stopstart == -1);
            %             % get rid of unpaired movement bouts at start and end
            %             if move_stop_frames(1) < move_start_frames(1)
            %                 move_stop_frames(1) = [];
            %             end
            %
            %             if move_start_frames(end) > move_stop_frames(end)
            %                 move_start_frames(end) = [];
            %             end
%             if plots_on
%                 figure;
%                 h1 = subplot(2,1,1);
%                 plot(lever_force_resample,'k')
%                 temp = lever_force_resample;
%                 temp(~lever_active) = NaN;
%                 hold on; plot(temp,'r');
%                 title('Position')
%                 h2 = subplot(2,1,2); hold on;
%                 plot(lever_velocity_resample_smooth,'m')
%                 plot(lever_velocity_envelope_smooth,'k')
%                 temp = lever_velocity_envelope_smooth;
%                 temp(~lever_active) = NaN;
%                 hold on; plot(temp,'r');
%                 hold on; plot(temp,'r');
%                 line(xlim,[movethresh movethresh],'color','m');
%                 title('Velocity')
%                 linkaxes([h1 h2],'x');
%             end
            
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
                % rewarded-movement/cue overlap can't be > 150ms
                % NEW: also, lever has to be inactive for 300ms before
                % cue/cue buffer time
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                cued_movement = rewarded_move_start > cue_sample(curr_trial_idx)-150 && ...
                    ~any(lever_active(cue_sample(curr_trial_idx)-450: ...
                    cue_sample(curr_trial_idx)-150));
%                 cued_movement = ~any(lever_active( ...
%                     cue_sample(curr_trial_idx)-5:cue_sample(curr_trial_idx)));
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
%                 move_thresh = lever_velocity_envelope_smooth ...
%                     (cue_sample(curr_trial_idx):reward_start_sample(curr_trial_idx)) ...
%                     > movethresh;
%                 if any(move_thresh)
%                     lever_param.reaction_time{curr_day}(curr_trial) = find(move_thresh,1);
%                 end
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
 

%% plot behavior 

clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' ...
    'AP76' 'AP78' 'AP79'};
curr_path = '/usr/local/lab/People/Andy/Data/lever_data';
animal_data = cell(length(animals),15);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    animal_data{curr_animal} = {};
    load([curr_path filesep animal 'lever_params.mat']);
        
    curr_param = cellfun(@(x,y) x(y), ...
        lever_param.reward_reaction, lever_param.cued_movement,'uni',false);
    
    %curr_param = lever_param.cued_movement;
    
    % normalize lever trace
    %curr_param = cellfun(@(x) cellfun(@(y) y./sqrt(sum(y.^2)),x,'uni',false), ...
    %    curr_param,'uni',false);
    %curr_param = cellfun(@(x) mean(reshape(horzcat(x{:})'*horzcat(x{:}), ...
    %    size(horzcat(x{:}),2)^2,1)),curr_param,'uni',false);
    %curr_param = cellfun(@(x) mean(horzcat(x{:}),2),curr_param,'uni',false);
    %curr_param = cellfun(@(x) x./sqrt(sum(x.^2)),curr_param,'uni',false);
    %curr_param = horzcat(curr_param{:})'*horzcat(curr_param{:});
    days = length(curr_param);
    
    % go through days, split into fractions of day
    day_frac = 1;
    num_params_day = cellfun(@length,curr_param,'uni',false);
    curr_param_cell = cellfun(@(x,y) ...
        mat2cell(x(1:y-mod(y,day_frac)), ...
        ones(1,day_frac)*floor(y/day_frac), ...
        1),curr_param,num_params_day,'uni',false);
  
   animal_data(curr_animal,1:days) = curr_param_cell;  
   %animal_data(curr_animal,1:days) = curr_param;
   %animal_data{curr_animal} = curr_param;
end

% %%%%
% % temporary: find number of consecutive very high reaction times
% temp_data = cellfun(@cell2mat,animal_data,'uni',false);
% temp_diff = cellfun(@(x) ...
% diff(find(x > (nanmedian(x) + 1.5*nanstd(x)))), ...
% temp_data,'uni',false);
% 
% a = cellfun(@(x) ~isempty(x),temp_data);
% 
% n_rep = 1000;
% temp_diff_shuffle = cellfun(@(x,y) ...
% diff(sort(randi(length(x),[n_rep length(y)]),2),[],2), ...
% temp_data(a),temp_diff(a),'uni',false);
% 
% temp_diff_conf = cellfun(@(x) ...
%     [mean(median(x,2)) prctile(median(x,2),[5 95])], ...
%     temp_diff_shuffle,'uni',false);
% % reshape
% temp_diff_shuffle_rs = cell(size(temp_data));
% temp_diff_shuffle_rs(a) = temp_diff_conf;
% 
% % plot
% figure;
% for i = 1:8
%    subplot(3,3,i); hold on;
%    curr_data = temp_diff_shuffle_rs(i,:);
%    curr_days = sum(cellfun(@(x) ~isempty(x), curr_data));
%    y1 = cellfun(@(x) x(2),curr_data(1:curr_days));
%    y2 = cellfun(@(x) x(3),curr_data(1:curr_days));
%    plot(y1,'r','linewidth',2);
%    plot(y2,'r','linewidth',2);
%    
%    plot(cellfun(@(x) nanmedian(x),temp_diff(1,1:curr_days)),'.b','MarkerSize',10)
%   
% end
% 
% %%%%

% get mean, sem of all data
curr_param_means = cell(1,14);
curr_param_stds = cell(1,14);
curr_param_sems = cell(1,14);

curr_param_bootcis = cell(1,14);
a = [];
for i = 1:14
    curr_data = animal_data(:,i);
    curr_data_full = horzcat(curr_data{:});
    curr_data_full_means = cellfun(@nanmean,curr_data_full);
    clear cat_data
    cat_data = vertcat(curr_data_full{:});
    
    curr_param_means{i} = nanmean(cat_data);
    
    curr_param_sems{i} = nanstd(vertcat(cat_data))./ ...
        sqrt(length(cat_data));
    
    %curr_param_bootcis{i} = bootci(1000,@nanmean,cat_data');
    
    %curr_param_means{i} = nanmean(curr_data_full_means,2);
    %curr_param_sems{i} = nanstd(curr_data_full_means,[],2)/ ...
    %    sqrt(size(curr_data_full_means,2));
    %curr_param_bootcis{i} = bootci(1000,@nanmean,curr_data_full_means');
    %curr_param_stds{i} = nanstd(curr_data_full_means,[],2);
    a = [a;vertcat(curr_data_full{:})];
end

% prepare data for each animal individually
animal_data_separate = cell(size(animal_data,1),1);
for i = 1:size(animal_data,1)
    animal_data_separate{i} = cellfun(@nanmedian, ...
        vertcat(animal_data{i,:}));
end

% plot mean/sem quarters of all days, nans in between
nan_breaks = num2cell(nan(1,14));
mean_plot = [curr_param_means(1:14);nan_breaks];
sem_plot = [curr_param_sems(1:14);nan_breaks];
std_plot = [curr_param_stds(1:14);nan_breaks];
% bootci_plot_l = [cellfun(@(x) x(1,:)', ...
%     curr_param_bootcis(1:14),'uni',false);nan_breaks];
% bootci_plot_u = [cellfun(@(x) x(2,:)', ...
%     curr_param_bootcis(1:14),'uni',false);nan_breaks];

mean_plot_cat = vertcat(mean_plot{:});
sem_plot_cat = vertcat(sem_plot{:});
std_plot_cat = vertcat(std_plot{:});
% bootci_plot_l_cat = vertcat(bootci_plot_l{:});
% bootci_plot_u_cat = vertcat(bootci_plot_u{:});

figure;
%errorbar([1:14+14*day_frac], ...
%    mean_plot_cat,bootci_plot_l_cat-mean_plot_cat, ...
%    bootci_plot_u_cat - mean_plot_cat,'k','linewidth',2)

errorbar(mean_plot_cat,sem_plot_cat,'k','linewidth',2);


%% Analyze/break down of movements
clear all

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

lever_movements_all = cell(size(animals));

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
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);
    
    lever_movements = struct;
    for curr_day = 1:length(days)
        % Initialize
        clearvars -except animals curr_animal animal days curr_day ...
            lever_movements data_path lever_movements_all
        
        day = num2str(days{curr_day});
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
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
        
        % Process
        
        % Only look at rewarded trials
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));
        
        time_surround = 1; % in seconds
        sample_surround = time_surround*xsg_samplerate;
        rewarded_lever_trace = nan(num_trials,sample_surround*2+1);
        rewarded_lever_trace_vel = nan(num_trials,sample_surround*2+1);
        
        % initialize vars
        lever_move_epoch_xsg = [];
        lever_move_epoch_samples_xsg = [];
        
        lever_move_epoch_downsamp_xsg = [];
        lever_move_epoch_downsamp_vel_xsg = [];
        lever_move_epoch_samples_downsamp_xsg = [];
        
        lever_still_epoch_downsamp_xsg = [];
        lever_still_epoch_samples_downsamp_xsg = [];
        
        rewarded_movements_xsg = logical([]);
        
        lever_reward_surround = [];
        lever_vel_reward_surround = [];
        
        lever_vel = [];
        
        % get lengths of all xsgs
        xsg_lengths = 0;
        for curr_xsg = 1:length(xsg_fullfilenames)-1
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            xsg_lengths(curr_xsg+1) = length(xsg_data.data.acquirer.trace_2);
        end
        xsg_lengths = cumsum(xsg_lengths);
        
        % go through all xsg files, grab movement properties
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            % create movement velocity, active times
            % quantify movement before
            
            [lever_active lever_force_resample lever_velocity_envelope_smooth]...
                = AP_parseLeverMovement(xsg_data);
            
            lever_stopstart = [0;diff(lever_active)];
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
            
            % save the lever trace surrounding the reward
            time_surround = 500; % ms
            reward_surround = cellfun(@(x) x*10-time_surround*10:x*10+time_surround*10, ...
                num2cell(reward_start_sample( ...
                ~isnan(reward_start_sample) & ...
                reward_start_sample*10-time_surround*10 > 0 & ...
                reward_start_sample*10+time_surround*10 < length(...
                xsg_data.data.acquirer.trace_2))),'uni',false);
            lever_reward_surround = [lever_reward_surround;cell2mat(cellfun(@(x) ...
                xsg_data.data.acquirer.trace_2(x)',reward_surround,'uni',false))];
            
            lever_vel_curr = [];
            lever_vel_curr = diff([0;lever_force_resample])';
            reward_surround_vel = cellfun(@(x) x-time_surround:x+time_surround, ...
                num2cell(reward_start_sample( ...
                ~isnan(reward_start_sample) & ...
                reward_start_sample-time_surround > 0 & ...
                reward_start_sample+time_surround < length(...
                lever_vel_curr))),'uni',false);        
            lever_vel_reward_surround = [lever_vel_reward_surround;cell2mat(cellfun(@(x) ...
                lever_vel_curr(x),reward_surround_vel,'uni',false))];
            
            lever_vel = [lever_vel lever_vel_curr];
            
            % break down lever into movement epochs
            lever_stopstart = diff(lever_active);
            lever_epochs = [0 find(lever_stopstart ~= 0)' length(lever_active)];
            lever_epoch_size = diff(lever_epochs)*10; %*10 for resampled            
            lever_epoch = mat2cell(xsg_data.data.acquirer.trace_2', ...
                1,lever_epoch_size);
            lever_epoch_samples = mat2cell( ...
                [1+xsg_lengths(curr_xsg): ...
                length(xsg_data.data.acquirer.trace_2)+xsg_lengths(curr_xsg)], ...
                1,lever_epoch_size);
   
            % the movement ones start from 2 and are every other epoch
            lever_move_epoch = lever_epoch(2:2:end);
            lever_move_epoch_samples = lever_epoch_samples(2:2:end);
            
            % split up again, but the downsampled one
            lever_epoch_downsamp = mat2cell(lever_force_resample', ...
                1,lever_epoch_size/10);
            lever_epoch_samples_downsamp = mat2cell( ...
                [1+xsg_lengths(curr_xsg)/10: ...
                length(lever_force_resample)+xsg_lengths(curr_xsg)/10], ...
                1,lever_epoch_size/10);
            
            lever_epoch_downsamp_vel = mat2cell(diff([0;lever_force_resample])', ...
                1,lever_epoch_size/10);
            
            % the movement ones start from 2 and are every other epoch
            lever_move_epoch_downsamp = lever_epoch_downsamp(2:2:end);
            lever_move_epoch_samples_downsamp = lever_epoch_samples_downsamp(2:2:end);
            
            lever_move_epoch_downsamp_vel = lever_epoch_downsamp_vel(2:2:end);
            
            lever_still_epoch_downsamp = lever_epoch_downsamp(1:2:end);
            lever_still_epoch_samples_downsamp = lever_epoch_samples_downsamp(1:2:end);
            
         
            % mark which movements were rewarded (/10 from resample)
            rewarded_movements = cellfun(@(x) ...
                any(ismember(round((x-xsg_lengths(curr_xsg))/10),reward_start_sample)), ...
                lever_move_epoch_samples);
            
            lever_move_epoch_xsg = [lever_move_epoch_xsg lever_move_epoch];
            lever_move_epoch_samples_xsg = [lever_move_epoch_samples_xsg ...
                lever_move_epoch_samples];
            
            lever_move_epoch_downsamp_xsg = [lever_move_epoch_downsamp_xsg ...
                lever_move_epoch_downsamp];
            lever_move_epoch_samples_downsamp_xsg = [lever_move_epoch_samples_downsamp_xsg ...
                lever_move_epoch_samples_downsamp];
            
            lever_still_epoch_downsamp_xsg = [lever_still_epoch_downsamp_xsg ...
                lever_still_epoch_downsamp];
            lever_still_epoch_samples_downsamp_xsg = [lever_still_epoch_samples_downsamp_xsg ...
                lever_still_epoch_samples_downsamp];
            
            lever_move_epoch_downsamp_vel_xsg = [lever_move_epoch_downsamp_vel_xsg ...
                lever_move_epoch_downsamp_vel];
            
            rewarded_movements_xsg = [rewarded_movements_xsg rewarded_movements];
            
        end
        
        lever_movements.lever_move_epoch_xsg{curr_day} = lever_move_epoch_xsg;
        lever_movements.lever_move_epoch_samples_xsg{curr_day} = lever_move_epoch_samples_xsg;
        lever_movements.rewarded_movements{curr_day} = rewarded_movements_xsg;
        
        lever_movements.lever_move_epoch_downsamp_xsg{curr_day} = lever_move_epoch_downsamp_xsg;
        lever_movements.lever_move_epoch_samples_downsamp_xsg{curr_day} = ...
            lever_move_epoch_samples_downsamp_xsg;
        lever_movements.lever_still_epoch_downsamp_xsg{curr_day} = lever_still_epoch_downsamp_xsg;
        lever_movements.lever_still_epoch_samples_downsamp_xsg{curr_day} = ...
            lever_still_epoch_samples_downsamp_xsg;
        
        lever_movements.lever_move_epoch_downsamp_vel_xsg{curr_day} = ...
            lever_move_epoch_downsamp_vel_xsg;
        
        lever_movements.lever_reward_surround{curr_day} = lever_reward_surround;
        lever_movements.lever_vel_reward_surround{curr_day} = lever_vel_reward_surround;
        lever_movements.lever_vel{curr_day} = lever_vel;

        
        
        disp(['Finished day: ' day]);
    end
    disp('Finished all')
    %save([animal 'lever_movements.mat'],'lever_movements');
    lever_movements_all{curr_animal} = lever_movements;
end
disp('Finished all animals')



% plot the max amplitude of each movement
figure;
for i = 1:8
    subplot(3,3,i);
    hold on;
    a = cellfun(@length,horzcat(lever_movements_all{i}.lever_move_epoch_downsamp_vel_xsg{:}));
    b = horzcat(lever_movements_all{i}.rewarded_movements{:});
    c = cumsum(cellfun(@length,lever_movements_all{i}.lever_move_epoch_xsg));
    
    plot(find(~b),a(~b),'.k');
    plot(find(b),a(b),'.r');
    for j = 1:length(c)
       line([c(j) c(j)],ylim) 
    end
    
    xlim([0 c(end)]);
    
end




%% Plot lever params (new)

% probably most important: 
% iti_quiescence_nonrewardedmove (fraction of iti spent moving, not during
% the rewarded movement)
% movement_reward_time: time from movement onset to reward
% reaction time: time from cue onset to movement

curr_param = cell(8,15);
for i = [1 2 3 4 5 6 7 8] 
    %curr_data = lever_param_all{i}.iti_quiescence_nonrewardedmove;
    curr_data = cellfun(@(x,y) x(y), ...
        lever_param_all{i}.reaction_time, ...
        lever_param_all{i}.cued_movement,'uni',false);
    curr_param(i,1:length(curr_data)) = curr_data;  
end
curr_param_mean = cellfun(@nanmean,curr_param);
%day_combine = {[1 2],[3 4],[5 6],[7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    %curr_data = vertcat(curr_param{:,day_combine{i}});
    curr_data = vertcat(curr_param_mean(:,day_combine{i}));
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);


%% Plot lever correlation at movement onset for rewarded movements

lever_template = cell(8,1);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    for session = sessions
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session});
        temp = temp(curr_cued(1:size(temp,1)),4801:100:10001);
        lever_template{curr_animal,session} = nanmean(temp);
    end
end
trial_template_corr = cell(15,15,8);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    for session = sessions
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        temp = lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(:,4801:100:10001);
        m = corrcoef([vertcat(lever_template{curr_animal,:})' ...
            temp(curr_cued(1:size(temp,1)),:)']);
        session_corr = mat2cell(m(1:length(sessions),length(sessions)+1:end), ...
            ones(length(sessions),1),sum(curr_cued(1:size(temp,1))));
        trial_template_corr(sessions,session,curr_animal) = session_corr;
    end
end
use_animals = [1 2 4 5 6 7 8];
cat_data = cell(15,15);
for i = 1:15
    for j = 1:15
        cat_data{i,j} = horzcat(trial_template_corr{i,j,use_animals});
    end
end


% Plot correlation of mean movement between days
use_animals = [1 2 4 5 6 7 8];
movecorr = nan(15,15,8);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    m = cell2mat(cellfun(@(x,y) nanmean(x(y(1:size(x,1)),:)), ...
        lever_param_all{curr_animal}.reward_aligned_lever, ...
        lever_param_all{curr_animal}.cued_movement,'uni',false)');
    m = m(:,5001:100:15001)';
    movecorr(sessions,sessions,curr_animal) = corrcoef(m);
end
use_data = movecorr(1:14,1:14,use_animals);
temp_data1 = cell2mat(arrayfun(@(x) diag(use_data(:,:,x),-1),1:size(use_data,3),'uni',false));
temp_data2 = cell2mat(arrayfun(@(x) diag(use_data(:,:,x),-2),1:size(use_data,3),'uni',false));
temp_data3 = cell2mat(arrayfun(@(x) diag(use_data(:,:,x),-3),1:size(use_data,3),'uni',false));

curr_data = [temp_data1(1:end-2,:) temp_data2(1:end-1,:) temp_data3(1:end,:)]';
figure;errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2);


% Plot example traces
example_animal = 8;
m = cellfun(@(x) nanmean(bsxfun( ...
    @minus,x(:,1:100:end),mean(x(:,1:100:3000),2))), ...
    lever_param_all{example_animal}.reward_aligned_lever,'uni',false);
figure; set(gcf,'Name',num2str(example_animal));
subplot(3,1,1);
plot(m{1},'k','linewidth',2);
xlim([0 252]);ylim([-0.1 0.3])
line([51 51],ylim,'linewidth',2,'color','r','linestyle','--');
title('Day 1')
subplot(3,1,2);
plot(m{3},'k','linewidth',2);
xlim([0 252]);ylim([-0.1 0.3])
line([51 51],ylim,'linewidth',2,'color','r','linestyle','--');
title('Day 3')
subplot(3,1,3);
plot(m{10},'k','linewidth',2);
xlim([0 252]);ylim([-0.1 0.3])
line([51 51],ylim,'linewidth',2,'color','r','linestyle','--');
title('Day 10')




%% Find best timing for correlation with lever


start_time = 4801;
end_times = [6001:1000:45001];
end_corrs = nan(length(end_times),10);

for curr_end = 1:length(end_times)
    
    lever_template = cell(8,1);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), ...
            lever_param_all{curr_animal}.reward_aligned_lever));
        for session = sessions
            curr_cued = lever_param_all{curr_animal}.cued_movement{session};
            temp = vertcat(lever_param_all{curr_animal}. ...
                reward_aligned_lever{session});
            temp = temp(curr_cued(1:size(temp,1)),start_time:100:end_times(curr_end));
            lever_template{curr_animal,session} = nanmean(temp);
        end
    end
    trial_template_corr = cell(15,15,8);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), ...
            lever_param_all{curr_animal}.reward_aligned_lever));
        for session = sessions
            curr_cued = lever_param_all{curr_animal}.cued_movement{session};
            temp = lever_param_all{curr_animal}. ...
                reward_aligned_lever{session}(:,start_time:100:end_times(curr_end));
            m = corrcoef([vertcat(lever_template{curr_animal,:})' ...
                temp(curr_cued(1:size(temp,1)),:)']);
            session_corr = mat2cell(m(1:length(sessions),length(sessions)+1:end), ...
                ones(length(sessions),1),sum(curr_cued(1:size(temp,1))));
            trial_template_corr(sessions,session,curr_animal) = session_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    cat_data = cell(15,15);
    for i = 1:15
        for j = 1:15
            cat_data{i,j} = horzcat(trial_template_corr{i,j,use_animals});
        end
    end
       
    cat_data_mean = cellfun(@nanmean,cat_data);
    %curr_end_corr = diag(cellfun(@nanmean,cat_data));
    curr_end_corr = cat_data_mean(10,1:10);
    end_corrs(curr_end,:) = curr_end_corr;
    disp(curr_end/length(end_times));
    
    
end

figure;hold on;
col = jet(size(end_corrs,1));
for i = 1:size(end_corrs,1)
    plot(end_corrs(i,:),'color',col(i,:))
end

%% Combining with activity: find correlations of lever with end

use_samples = 4901:100:35001;
%use_samples = 3000:100:7000;

% TK: Make, plot properties of n k-means groups of total lever
k_clusters = 3;
lever_k_templates = cell(8,1);
inter_cluster_cc = cell(8,1);
template_cluster_cc = cell(8,15);
move_reward_time = cell(8,1);
lever_trials_used = cell(8,15);
lever_used = cell(8,15);
template_trials = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    use_move_reward_time = cell(length(session),1);
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
        
        use_move_reward_time{session} = lever_param_all{curr_animal}. ...
            move_reward_time{session}(use_trials);
    end
    %all_cat_lever = vertcat(curr_cat_lever{:});  
    all_cat_lever_downsamp = vertcat(curr_cat_lever{:});  
    
    all_cat_lever_use = all_cat_lever_downsamp;
    kidx = kmeans(all_cat_lever_use,k_clusters);
    lever_k_templates{curr_animal} = grpstats(all_cat_lever_use,kidx);
    
    kidx_split = cell(max(sessions),1);
    kidx_split(sessions) = mat2cell(kidx,cellfun(@(x) size(x,1),curr_cat_lever),1);
    template_trials(curr_animal,sessions) = kidx_split;
    
    for session = sessions
        temp_cc = corrcoef([lever_k_templates{curr_animal}' ...
            total_curr_cat_lever{session}']);
        template_cluster_cc{curr_animal,session} = temp_cc( ...
            1:k_clusters,k_clusters+1:end);
    end
    
    for curr_k_cluster = 1:k_clusters
        inter_cluster_cc{curr_animal}(curr_k_cluster) = ...
            nanmean(AP_itril(corrcoef(all_cat_lever_use( ...
            kidx == curr_k_cluster,:)'),-1));
    end
    
    cat_move_reward_time = vertcat(use_move_reward_time{:});
    move_reward_time{curr_animal} = grpstats(cat_move_reward_time,kidx);
    
end

use_animals = [1 2 4 5 6 7 8];
for curr_animal = use_animals;
    figure;
    for curr_k_cluster = 1:k_clusters
        subplot(4,k_clusters,curr_k_cluster);
        plot(lever_k_templates{curr_animal}(curr_k_cluster,:),'k','linewidth',2)
        axis off
        
        subplot(4,k_clusters,2*k_clusters+curr_k_cluster);
        sessions = find(cellfun(@(x) ~isempty(x),template_cluster_cc(curr_animal,:)));
        curr_template_cluster_cc = cellfun(@(x) nanmean(x(curr_k_cluster,:),2), ...
            template_cluster_cc(curr_animal,sessions));
        plot(sessions,curr_template_cluster_cc,'k','linewidth',2);
        if curr_k_cluster == 1
            ylabel('Template cc');
        end
    end
    subplot(4,k_clusters,curr_k_cluster+1:curr_k_cluster*2);
    plot(inter_cluster_cc{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Inter-cluster cc');
    subplot(4,k_clusters,3*curr_k_cluster+1:curr_k_cluster*4);
    plot(move_reward_time{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Move-reward time');
    set(gcf,'Name',num2str(curr_animal));
end

% Haven't been using below

% % Template based on average of last three days
% lever_template = cell(8,1);
% for curr_animal = 1:8
%     sessions = find(cellfun(@(x) ~isempty(x), ...
%         lever_param_all{curr_animal}.reward_aligned_lever));
%     curr_cat_lever = cell(3,1);
%     for end_session = 1:3
%         session = sessions(end - (end_session - 1));
%         curr_cued = lever_param_all{curr_animal}.cued_movement{session};
%         temp = vertcat(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session});
%         temp = temp(curr_cued(1:size(temp,1)),:);
%         curr_cat_lever{end_session} = temp;
%     end
%     lever_template{curr_animal} = nanmean(vertcat(curr_cat_lever{:}));
% end
% 
% % Template based on k-means of last three days
% lever_template_1 = cell(8,1);
% lever_template_2 = cell(8,1);
% learned_corr = nan(8,1);
% nonlearned_corr = nan(8,1);
% template_trials = cell(8,2);
% for curr_animal = 1:8
%     sessions = find(cellfun(@(x) ~isempty(x), ...
%         lever_param_all{curr_animal}.reward_aligned_lever));
%     curr_cat_lever = cell(3,1);
%     for end_session = 1:3
%         session = sessions(end - (end_session - 1));
%         curr_cued = lever_param_all{curr_animal}.cued_movement{session};
%         curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session},1));
%         temp = vertcat(lever_param_all{curr_animal}. ...
%             reward_aligned_lever{session}(curr_cued,:));      
%         curr_cat_lever{end_session} = temp;
%     end
%     all_cat_lever = vertcat(curr_cat_lever{:});
%     
%     % TRYING OUT: align instead by reward
%     m = cellfun(@(x,y) x(y),...
%         lever_param_all{curr_animal}.move_reward_time(sessions(end-2):sessions(end)), ...
%         lever_param_all{curr_animal}.cued_movement(sessions(end-2):sessions(end)),'uni',false);
%     m1 = vertcat(m{:});
%     m_fix = m1*10+5000;
%     all_cat_lever_realign = nan(size(all_cat_lever,1),10001);
%     for i = 1:size(all_cat_lever,1);
%         if m_fix(i)+5000 > size(all_cat_lever,2) || isnan(m_fix(i))
%             continue
%         end
%         all_cat_lever_realign(i,:) = all_cat_lever(i,m_fix(i)-5000:m_fix(i)+5000);       
%     end
%     
%     use_trials = ~any(isnan(all_cat_lever),2) & ~all(all_cat_lever == 0,2);
%     all_cat_lever_use = all_cat_lever(use_trials,use_samples);
%     kidx = kmeans(all_cat_lever_use,2);
%     
%     % pick out more inter-correlated ("learned") press
%     kidx_1_corr = nanmean(AP_itril(corrcoef(all_cat_lever_use(kidx == 1,:)')));
%     kidx_2_corr = nanmean(AP_itril(corrcoef(all_cat_lever_use(kidx == 2,:)')));
%     
%     if kidx_1_corr > kidx_2_corr
%         learned_kidx = 1;
%         learned_corr(curr_animal) = kidx_1_corr;
%         nonlearned_corr(curr_animal) = kidx_2_corr;
%     elseif kidx_2_corr > kidx_1_corr
%         learned_kidx = 2;
%         learned_corr(curr_animal) = kidx_2_corr;
%         nonlearned_corr(curr_animal) = kidx_1_corr;
%     end
%     nonlearned_kidx = setdiff([1 2],learned_kidx);       
%     
%     lever_template_1{curr_animal} = nanmean(all_cat_lever_use(kidx == learned_kidx,:));
%     lever_template_2{curr_animal} = nanmean(all_cat_lever_use(kidx == nonlearned_kidx,:));
%     
%     
%     % pick out which trials are used to form template (for use in activity)
%     kidx_full = zeros(length(use_trials),1);
%     use_trials_idx = find(use_trials);
%     kidx_full(use_trials_idx(kidx == 1)) = 1;
%     kidx_full(use_trials_idx(kidx == 2)) = 2;
%     curr_trial_num = cellfun(@(x) ...
%         find(~any(isnan(x),2) & ~all(x == 0,2)),curr_cat_lever,'uni',false);
%     template_trials{curr_animal,1} = cellfun(@(x,y) x(y == learned_kidx),curr_trial_num, ...
%         mat2cell(kidx_full,cellfun(@length,curr_trial_num),1),'uni',false);
%     template_trials{curr_animal,2} = cellfun(@(x,y) x(y == nonlearned_kidx),curr_trial_num, ...
%         mat2cell(kidx_full,cellfun(@length,curr_trial_num),1),'uni',false);
% end
% 
% lever_corr_1 = cell(8,15);
% lever_corr_2 = cell(8,15);
% lever_corr_trials = cell(8,15);
% lever_use = cell(8,15);
% use_dot = cell(8,15);
% for curr_animal = 1:8
%     sessions = find(cellfun(@(x) ~isempty(x),lever_param_all{curr_animal}.cued_movement));
%     for curr_session = sessions    
%         
%         curr_cued = lever_param_all{curr_animal}.cued_movement{session};
%         
%         curr_lever = lever_param_all{curr_animal}. ...
%             reward_aligned_lever{curr_session};
%         
%         use_trials = 1:size(curr_lever,1);%trial_template_trials{curr_animal,curr_session};
%         end_trials = use_trials > size(curr_lever,1);
%         use_trials(end_trials) = [];
%         
%         % use only cued movement trials
%         use_trials_cued = intersect(find(curr_cued),use_trials);
%                 
%         curr_lever_trials = curr_lever(use_trials,:);
%         
% %         lever_corr_grid = corrcoef([lever_template{curr_animal}(use_samples)' ...
% %             curr_lever_trials(:,use_samples)']);
% %         lever_corr{curr_animal,curr_session} = lever_corr_grid(1,2:end);  
%         lever_corr_grid_1 = corrcoef([lever_template_1{curr_animal}' ...
%             curr_lever_trials(:,use_samples)']);
%         lever_corr_1{curr_animal,curr_session} = lever_corr_grid_1(1,2:end);  
%         
%         lever_corr_grid_2 = corrcoef([lever_template_2{curr_animal}' ...
%             curr_lever_trials(:,use_samples)']);
%         lever_corr_2{curr_animal,curr_session} = lever_corr_grid_2(1,2:end);  
%         
%         lever_corr_trials{curr_animal,curr_session} = use_trials;
%         lever_use{curr_animal,curr_session} = curr_lever_trials(:,use_samples);
% %         use_trials_dot = ~end_trials; %& ismember( ...
% %             %trial_template_trials{curr_animal,curr_session},find(curr_cued));
% %         use_dot{curr_animal,curr_session} = ...
% %             trial_template_corr{curr_animal,curr_session}(use_trials_dot);
%     end
%     curr_animal
% end

%% Cluster movements to find template, align by xcorr

lever_downsamp = 100;
use_samples = 49:350; % new: in 10 ms groups
%use_samples = 4901:100:35001;
%use_samples = 3000:100:7000;

% TK: Make, plot properties of n k-means groups of total lever
k_clusters = 3;
lever_k_templates = cell(8,1);
inter_cluster_cc = cell(8,1);
template_cluster_cc = cell(8,15);
move_reward_time = cell(8,1);
lever_trials_used = cell(8,15);
lever_used = cell(8,15);
template_trials = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    use_move_reward_time = cell(length(session),1);
    % extracting all lever
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        
        % try using cued and uncued trials
        curr_cued = true(size(curr_cued));
        
        temp_lever = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp_lever),2) & ...
            ~all(temp_lever == 0,2);
        temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
        total_curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    lever_used(curr_animal,sessions) = total_curr_cat_lever;
    % creating template
    for session = sessions(1:5);
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        
        % try using cued and uncued trials
        curr_cued = true(size(curr_cued));
        
        temp_lever = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp_lever),2) & ...
            ~all(temp_lever == 0,2);
        temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
        
        curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);
        
        use_move_reward_time{session} = lever_param_all{curr_animal}. ...
            move_reward_time{session}(use_trials);
    end
    %all_cat_lever = vertcat(curr_cat_lever{:});  
    all_cat_lever_downsamp = vertcat(curr_cat_lever{:});  
    
    all_cat_lever_use = all_cat_lever_downsamp;
    kidx = kmeans(all_cat_lever_use,k_clusters);
    lever_k_templates{curr_animal} = grpstats(all_cat_lever_use,kidx);
    
    kidx_split = cell(max(sessions),1);
    kidx_split(sessions) = mat2cell(kidx,cellfun(@(x) size(x,1),curr_cat_lever),1);
    template_trials(curr_animal,sessions) = kidx_split;
    
    for session = sessions
        temp_cc = corrcoef([lever_k_templates{curr_animal}' ...
            total_curr_cat_lever{session}']);
        template_cluster_cc{curr_animal,session} = temp_cc( ...
            1:k_clusters,k_clusters+1:end);
    end
    
    for curr_k_cluster = 1:k_clusters
        inter_cluster_cc{curr_animal}(curr_k_cluster) = ...
            nanmean(AP_itril(corrcoef(all_cat_lever_use( ...
            kidx == curr_k_cluster,:)'),-1));
    end
    
    cat_move_reward_time = vertcat(use_move_reward_time{:});
    move_reward_time{curr_animal} = grpstats(cat_move_reward_time,kidx);
    
end

use_animals = [1 2 4 5 6 7 8];
for curr_animal = use_animals;
    figure;
    for curr_k_cluster = 1:k_clusters
        subplot(4,k_clusters,curr_k_cluster);
        plot(lever_k_templates{curr_animal}(curr_k_cluster,:),'k','linewidth',2)
        axis off
        
        subplot(4,k_clusters,2*k_clusters+curr_k_cluster);
        sessions = find(cellfun(@(x) ~isempty(x),template_cluster_cc(curr_animal,:)));
        curr_template_cluster_cc = cellfun(@(x) nanmean(x(curr_k_cluster,:),2), ...
            template_cluster_cc(curr_animal,sessions));
        plot(sessions,curr_template_cluster_cc,'k','linewidth',2);
        if curr_k_cluster == 1
            ylabel('Template cc');
        end
    end
    subplot(4,k_clusters,curr_k_cluster+1:curr_k_cluster*2);
    plot(inter_cluster_cc{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Inter-cluster cc');
    subplot(4,k_clusters,3*curr_k_cluster+1:curr_k_cluster*4);
    plot(move_reward_time{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Move-reward time');
    set(gcf,'Name',num2str(curr_animal));
end

keyboard

%%%%%%%%
% CHOOSE TEMPLATES BASED ON ABOVE
%%%%%%%%

% find best alignment for each lever trace to template

up_templates = [3 3 nan 1 3 3 2 2];
use_animals = [1 2 4 5 6 7 8];

lever_aligned = cell(8,15);
lever_aligned_cc = cell(8,15);
align_shift = cell(8,15);
for curr_animal = use_animals
    sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
    curr_template = lever_k_templates{curr_animal}(up_templates(curr_animal),:);
    for curr_session = sessions
        lever_xcorr = arrayfun(@(x) xcov(lever_used{curr_animal,curr_session}(x,:), ...
            curr_template),1:size(lever_used{curr_animal,curr_session},1),'uni',false);
        [c ci] = max(vertcat(lever_xcorr{:}),[],2);
        ci_fix = -(ci - size(curr_template,2));
        align_shift{curr_animal,curr_session} = ci_fix;
        curr_lever_aligned = ...
            arrayfun(@(x) circshift(lever_used{curr_animal,curr_session}(x,:), ...
            [0 ci_fix(x)]),1:length(ci_fix),'uni',false);
        lever_aligned{curr_animal,curr_session} = vertcat(curr_lever_aligned{:});
        
        curr_lever_aligned_cc_grid = corrcoef( ...
            [lever_k_templates{curr_animal}(up_templates(curr_animal),:)' ...
            vertcat(curr_lever_aligned{:})']);
        lever_aligned_cc{curr_animal,curr_session} = ...
            curr_lever_aligned_cc_grid(1,2:end);
    end
end

%% Same as above, but split early/late days for kmeans
use_samples = 4901:100:35001;
%use_samples = 3000:100:7000;

% TK: Make, plot properties of n k-means groups of total lever

% templates from early days
k_clusters_early = 3;
lever_k_templates_early = cell(8,1);
inter_cluster_cc_early = cell(8,1);
template_cluster_cc_early = cell(8,15);
move_reward_time_early = cell(8,1);
lever_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    use_move_reward_time = cell(length(session),1);
    % extracting all activity
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        total_curr_cat_lever{session} = temp(:,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    % creating template
    for session = sessions(1:3);
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        curr_cat_lever{session} = temp(:,use_samples);
        
        use_move_reward_time{session} = lever_param_all{curr_animal}. ...
            move_reward_time{session}(use_trials);
    end
    %all_cat_lever = vertcat(curr_cat_lever{:});  
    all_cat_lever_downsamp = vertcat(curr_cat_lever{:});  
    
    all_cat_lever_use = all_cat_lever_downsamp;
    kidx = kmeans(all_cat_lever_use,k_clusters_early);
    lever_k_templates_early{curr_animal} = grpstats(all_cat_lever_use,kidx);
    
    kidx_split = cell(max(sessions),1);
    kidx_split(sessions) = mat2cell(kidx,cellfun(@(x) size(x,1),curr_cat_lever),1);
    
    for session = sessions
        temp_cc = corrcoef([lever_k_templates_early{curr_animal}' ...
            total_curr_cat_lever{session}']);
        template_cluster_cc_early{curr_animal,session} = temp_cc( ...
            1:k_clusters_early,k_clusters_early+1:end);
    end
    
    for curr_k_cluster = 1:k_clusters_early
        inter_cluster_cc_early{curr_animal}(curr_k_cluster) = ...
            nanmean(AP_itril(corrcoef(all_cat_lever_use( ...
            kidx == curr_k_cluster,:)'),-1));
    end
    
    cat_move_reward_time = vertcat(use_move_reward_time{:});
    move_reward_time_early{curr_animal} = grpstats(cat_move_reward_time,kidx);
    
end

% templates from late days
k_clusters_late = 3;
lever_k_templates_late = cell(8,1);
inter_cluster_cc_late = cell(8,1);
template_cluster_cc_late = cell(8,15);
move_reward_time_late = cell(8,1);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    use_move_reward_time = cell(length(session),1);
    % extracting all activity
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        total_curr_cat_lever{session} = temp(:,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    % creating template
    for session = sessions(end-2:end);
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        temp = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp),2) & ...
            ~all(temp == 0,2);
        curr_cat_lever{session} = temp(:,use_samples);
        
        use_move_reward_time{session} = lever_param_all{curr_animal}. ...
            move_reward_time{session}(use_trials);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    %all_cat_lever = vertcat(curr_cat_lever{:});  
    all_cat_lever_downsamp = vertcat(curr_cat_lever{:});  
    
    all_cat_lever_use = all_cat_lever_downsamp;
    kidx = kmeans(all_cat_lever_use,k_clusters_late);
    lever_k_templates_late{curr_animal} = grpstats(all_cat_lever_use,kidx);
    
    kidx_split = cell(max(sessions),1);
    kidx_split(sessions) = mat2cell(kidx,cellfun(@(x) size(x,1),curr_cat_lever),1);
    
    for session = sessions
        temp_cc = corrcoef([lever_k_templates_late{curr_animal}' ...
            total_curr_cat_lever{session}']);
        template_cluster_cc_late{curr_animal,session} = temp_cc( ...
            1:k_clusters_late,k_clusters_late+1:end);
    end
    
    for curr_k_cluster = 1:k_clusters_late
        inter_cluster_cc_late{curr_animal}(curr_k_cluster) = ...
            nanmean(AP_itril(corrcoef(all_cat_lever_use( ...
            kidx == curr_k_cluster,:)'),-1));
    end
    
    cat_move_reward_time = vertcat(use_move_reward_time{:});
    move_reward_time_late{curr_animal} = grpstats(cat_move_reward_time,kidx);
    
end

k_clusters = k_clusters_early + k_clusters_late;
lever_k_templates = cellfun(@(x,y) vertcat(x,y), ...
    lever_k_templates_early,lever_k_templates_late,'uni',false);
inter_cluster_cc = cellfun(@(x,y) horzcat(x,y), ...
    inter_cluster_cc_early,inter_cluster_cc_late,'uni',false);
template_cluster_cc = cellfun(@(x,y) vertcat(x,y), ...
    template_cluster_cc_early,template_cluster_cc_late,'uni',false);
move_reward_time = cellfun(@(x,y) vertcat(x,y), ...
    move_reward_time_early,move_reward_time_late,'uni',false);

% plot params

use_animals = [1 2 4 5 6 7 8];
for curr_animal = use_animals;
    figure;
    for curr_k_cluster = 1:k_clusters
        subplot(4,k_clusters,curr_k_cluster);
        plot(lever_k_templates{curr_animal}(curr_k_cluster,:),'k','linewidth',2)
        axis off
        
        subplot(4,k_clusters,2*k_clusters+curr_k_cluster);
        sessions = find(cellfun(@(x) ~isempty(x),template_cluster_cc(curr_animal,:)));
        curr_template_cluster_cc = cellfun(@(x) nanmean(x(curr_k_cluster,:),2), ...
            template_cluster_cc(curr_animal,sessions));
        plot(sessions,curr_template_cluster_cc,'k','linewidth',2);
        if curr_k_cluster == 1
            ylabel('Template cc');
        end
    end
    subplot(4,k_clusters,curr_k_cluster+1:curr_k_cluster*2);
    plot(inter_cluster_cc{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Inter-cluster cc');
    subplot(4,k_clusters,3*curr_k_cluster+1:curr_k_cluster*4);
    plot(move_reward_time{curr_animal},'ok','MarkerSize',10, ...
        'MarkerFaceColor',[1 0 0]);
    ylabel('Move-reward time');
    set(gcf,'Name',num2str(curr_animal));
end



%% Find correlation of 3s trace


lever_downsamp = 100;
use_samples = 49:350; % new: in 10 ms groups
lever_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), ...
        lever_param_all{curr_animal}.reward_aligned_lever));
    total_curr_cat_lever = cell(length(sessions),1);
    curr_cat_lever = cell(length(sessions),1);
    use_move_reward_time = cell(length(session),1);
    % extracting all lever
    for session = sessions;
        curr_cued = lever_param_all{curr_animal}.cued_movement{session};
        curr_cued = curr_cued(1:size(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session},1));
        
        % try using cued and uncued trials
        curr_cued = true(size(curr_cued));
        
        temp_lever = vertcat(lever_param_all{curr_animal}. ...
            reward_aligned_lever{session}(curr_cued,:));   
        use_trials = ~any(isnan(temp_lever),2) & ...
            ~all(temp_lever == 0,2);
        temp_lever_downsamp = resample(temp_lever',1,lever_downsamp)';
        total_curr_cat_lever{session} = temp_lever_downsamp(use_trials,use_samples);
        
        curr_cued_idx = find(curr_cued);
        lever_trials_used{curr_animal,session} = curr_cued_idx(use_trials);
    end
    lever_used(curr_animal,sessions) = total_curr_cat_lever;
    
    disp(curr_animal)
end


% plot trial-trial correlation
lever_trial_corr = cellfun(@(x) AP_itril(corrcoef(x'),-1),lever_used,'uni',false);
lever_trial_corr_cat = cell(14,1);
use_animals = [1 2 4 5 6 7 8];
for i = 1:14
   lever_trial_corr_cat{i} = vertcat(lever_trial_corr{use_animals,i}); 
end
lever_trial_corr_mean = cellfun(@nanmedian,lever_trial_corr_cat); 
lever_trial_corr_sem = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), ...
    lever_trial_corr_cat);
figure;
errorbar(lever_trial_corr_mean,lever_trial_corr_sem,'k','linewidth',2)

% plot avg-avg corrlation
lever_avg = cellfun(@nanmean,lever_used,'uni',false);
lever_avg_corr = nan(15,15,8);
for curr_animal = 1:8
   sessions = find(cellfun(@(x) ~isempty(x),lever_used(curr_animal,:)));
   curr_lever = vertcat(lever_avg{curr_animal,sessions});
   lever_avg_corr(sessions,sessions,curr_animal) = corrcoef(curr_lever');  
end
use_animals = [1 2 4 5 6 7 8];
figure;
imagesc(nanmean(lever_avg_corr(1:14,1:14,use_animals),3))
colormap(hot);caxis([0 1])
xlabel('Session')
ylabel('Session')
title('Mean lever press correlation')














