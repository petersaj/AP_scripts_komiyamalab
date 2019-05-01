function [data,analysis] = AP_corticospinal_load_layer23
% [data,analysis,classified_rois] = AP_corticospinal_load_layer23 
% Load and process layer 2/3 data using updated protocols and methods for
% comparison to corticospinal data
% (the output structures can be passed to
% AP_classify_movement_cells_continuous)

%% L2/3 animals

animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

data = struct('animal',cell(length(animals),1), 'bhv',cell(length(animals),1), ...
    'im',cell(length(animals),1));

analysis = struct('lever',cell(length(animals),1), ...
    'surrounding_frames',cell(length(animals),1));


%% Load and process data (only for unfilled pyr cells)

movement_trace_all = cell(length(animals),15);
move_epoch_frames = cell(length(animals),15);
rewarded_movements = cell(length(animals),15);
cued_rewarded_movements = cell(length(animals),15);
catch_rewarded_movements = cell(length(animals),15);

cue_frames_all = cell(length(animals),15);
reward_frames_all = cell(length(animals),15);
cued_reward_frames_all = cell(length(animals),15);
catch_reward_frames_all = cell(length(animals),15);

trial_frames_all = cell(length(animals),15);

disp('Loading and processing L2/3 data')
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
    
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    
    bhv_path = [data_path filesep animal '_dispatcher'];
    xsg_path = [data_path filesep animal '_xsg'];
    
    bhv_all = dir(bhv_path);
    xsg_all = dir(xsg_path);
    
    days = cellfun(@(x) x(1:6),{xsg_all(3:end).name}, ...
        'UniformOutput',false);
    
    for curr_day = 1:length(days);
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        % Re-do baseline estimation: allow for more nans than before
        % Background subtraction
        [~,bg_baseline] = AP_baselineEstimation(im.roi_trace_long_bg,30);
        bg_sub_trace = im.roi_trace_long - (im.roi_trace_long_bg - bg_baseline);
        
        im.roi_trace_df = AP_baselineEstimation(bg_sub_trace,30);
        
        %%%% Find and store activity
        pyr_activity = ...
            AP_caEvents_thresh(im.roi_trace_df(pyr_unfilled,:),3,2);
       
        %%%% Find and store movements
        
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
                
        % get behavior data
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
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));        
        
        catch_trials = [];
        if isfield(saved_history,'CatchSection_catch_trial')
            catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
            % there's one or more uncompleted trials at the end
            catch_trials = find(catch_trials(1:length(bhv)));
        end
                       
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
             
        lever_active_full = cell(length(xsg_filenames),1);
        lever_position_full = cell(length(xsg_filenames),1);
        lever_velocity_full = cell(length(xsg_filenames),1);
        cue_frames_full =  cell(length(xsg_filenames),1);
        reward_frames_full =  cell(length(xsg_filenames),1);
        cued_reward_frames_full = cell(length(xsg_filenames),1);
        catch_reward_frames_full = cell(length(xsg_filenames),1);
        
        trial_frames = cellfun(@(x) nan(1,x),num2cell(loop_size),'uni',false);
        
        lever_use_trials = false(num_trials,1);
        analysis(curr_animal).lever(curr_day). ...
            cued_movement_trials = false(num_trials,1);
        analysis(curr_animal).lever(curr_day). ...
            rewarded_movement = cell(num_trials,1);
        analysis(curr_animal).lever(curr_day). ...
            rewarded_movement_fixtime = cell(num_trials,1);
            
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            [lever_active,lever_force_resample,~,lever_force_velocity] = ...
                AP_parseLeverMovement_updated(xsg_data);
            
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
                
                rewarded_move_stop = find(lever_active( ...
                    reward_start_sample(curr_trial_idx):end) == 0,1) + ...
                    reward_start_sample(curr_trial_idx);
                
                cued_movement(curr_trial_idx) = ...
                    rewarded_move_start > cue_sample(curr_trial_idx) && ...
                    ~any(lever_active(cue_sample(curr_trial_idx)-100: ...
                    cue_sample(curr_trial_idx)));
                                
                if ismember(curr_trial,catch_trials)
                    catch_movement(curr_trial_idx) = true;
                end
                
                % Store cued movement status
                 analysis(curr_animal).lever(curr_day). ...
                    cued_movement_trials(curr_trial,1) = cued_movement(curr_trial_idx);              
                
                % Store rewarded movements
                rewarded_movement_back = 1000;
                rewarded_movement_forward = 4000;
                
                % don't store movement if not enough surrounding samples
                if rewarded_move_start - rewarded_movement_back < 1 || ...
                        rewarded_move_start + rewarded_movement_forward > ...
                        length(lever_force_resample)
                    continue
                end
                
                curr_rewarded_movement = lever_force_resample( ...
                    rewarded_move_start:rewarded_move_stop);
                
                curr_rewarded_movement_fixtime = lever_force_resample( ...
                    rewarded_move_start - rewarded_movement_back: ...
                    rewarded_move_start + rewarded_movement_forward);
                
                analysis(curr_animal).lever(curr_day). ...
                    rewarded_movement{curr_trial,1} = curr_rewarded_movement;
                
                analysis(curr_animal).lever(curr_day). ...
                    rewarded_movement_fixtime{curr_trial,1} = curr_rewarded_movement_fixtime; 
                
                % store trial that has full lever data
                lever_use_trials(curr_trial) = true;
                
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
            
            % Resample lever position and speed to frames
            [n d] = rat(loop_size(curr_xsg)/length(lever_force_resample));
            lever_force_frame_resample = resample(lever_force_resample,n,d);
            lever_position_full{curr_xsg} = ...
                lever_force_frame_resample(1:loop_size(curr_xsg));
            
            lever_force_velocity_frame_resample = resample(lever_force_velocity,n,d);
            lever_velocity_full{curr_xsg} = ...
                lever_force_velocity_frame_resample(1:loop_size(curr_xsg));
                      
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
        
        lever_position = vertcat(lever_position_full{:});
        lever_speed = vertcat(lever_velocity_full{:});
        
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
        
        % Align and store activity to movement on/offset (CR movements)       
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_day}( ...
            cued_rewarded_movements{curr_animal,curr_day});
        
        %%% Make sure lever and activity trials match
        
        % (only use trials with lever data)
        curr_rewardtrials = cellfun(@(x) trial_frames_all{curr_animal,curr_day} ...
            (x(1)),curr_rewardmoves);
        curr_rewardmoves(~ismember(curr_rewardtrials,find(lever_use_trials))) = [];
        
        % (the other direction: only keep lever with activity)
        cued_lever_trials = find(lever_use_trials & analysis(curr_animal).lever(curr_day).cued_movement_trials);
        lever_noact = setdiff(cued_lever_trials,curr_rewardtrials);
        analysis(curr_animal).lever(curr_day).rewarded_movement(lever_noact) = ...
            cell(length(lever_noact),1);
        analysis(curr_animal).lever(curr_day).rewarded_movement_fixtime(lever_noact) = ...
            cell(length(lever_noact),1);
               
        curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
        curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
                      
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_day}),curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
                   
        surrounding_frames = -90:90;
        curr_moveonset_surround = cellfun(@(x) x+surrounding_frames,curr_moveonsets,'uni',false)';
        curr_reward_surround = cellfun(@(x) x+surrounding_frames,curr_rewardframes,'uni',false)';
        curr_moveoffset_surround = cellfun(@(x) x+surrounding_frames,curr_moveoffsets,'uni',false)';
        
        elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
            length(movement_trace_all{curr_animal,curr_day})),curr_moveonset_surround, ...
            curr_reward_surround,curr_moveoffset_surround);
        curr_moveonset_surround(elim_moves) = [];
        curr_reward_surround(elim_moves) = [];
        curr_moveoffset_surround(elim_moves) = [];
        analysis(curr_animal).lever(curr_day).rewarded_movement(curr_rewardtrials(elim_moves)) = ...
            cell(sum(elim_moves),1);
        analysis(curr_animal).lever(curr_day).rewarded_movement_fixtime(curr_rewardtrials(elim_moves)) = ...
            cell(sum(elim_moves),1);
        
        curr_moveonset_aligned_pyr = cellfun(@(x) pyr_activity(:,x), ...
            curr_moveonset_surround,'uni',false);
        curr_moveoffset_aligned_pyr = cellfun(@(x) pyr_activity(:,x), ...
            curr_moveoffset_surround,'uni',false);
        
        analysis(curr_animal).surrounding_frames = surrounding_frames;
        analysis(curr_animal).im(curr_day).move_onset_aligned = ...
            permute(cat(3,curr_moveonset_aligned_pyr{:}),[3 2 1]);
        analysis(curr_animal).im(curr_day).move_offset_aligned = ...
            permute(cat(3,curr_moveoffset_aligned_pyr{:}),[3 2 1]);
             
        %%%%
                
        % Store in data/analysis structures
        data(curr_animal).im(curr_day).roi_trace_df = im.roi_trace_df(pyr_unfilled,:);
        data(curr_animal).im(curr_day).roi_trace_thresh = pyr_activity;
        analysis(curr_animal).lever(curr_day).lever_move_frames = movement_trace';    
               
        % Store lever position, velocity, movement in frames
        median_lever = median(smooth(lever_position,30,'loess'));

        analysis(curr_animal).lever(curr_day).lever_position_frames = ...
            lever_position-median_lever;
        
        analysis(curr_animal).lever(curr_day).lever_speed_frames = ...
            lever_speed;       
 
        
        disp(['Loaded ' animal ' day ' num2str(curr_day)]);
        
    end
    
    data(curr_animal).animal = animal;
  
end

disp('done.')









