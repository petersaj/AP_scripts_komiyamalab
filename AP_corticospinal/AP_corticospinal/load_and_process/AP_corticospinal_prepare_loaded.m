function analysis = AP_corticospinal_prepare_loaded(data)
% analysis = AP_corticospinal_prepare_loaded(data)
%
% Extract lever/task information
% Extract df traces locked to cued-rewarded movements

analysis = struct; 


%% Lever and task behavior 
clearvars -except data analysis

for curr_animal = 1:length(data)
    
    num_sessions = length(data(curr_animal).im);
    
    for curr_session = 1:num_sessions
        
        
        % Skip if no behavior file (forgot to save)
        if isempty(data(curr_animal).bhv(curr_session).bhv_times)
            continue
        end
        
        num_trials = length(data(curr_animal).bhv(curr_session).bhv_times);
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),data(curr_animal).bhv(curr_session).bhv_times);
        
        num_frames = min(size(data(curr_animal).im(curr_session).roi_trace_df,2), ...
            length(data(curr_animal).bhv(curr_session).frame_times));
        
        
        % Parse movement
        [lever_active,lever_force_resample, ...
            lever_force_smooth,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated(data(curr_animal).bhv(curr_session).lever_force);
        
        %% 1) Get lever params for each frame
        
        lever_velocity = [0;smooth(diff(lever_force_smooth),5)];
        
        frame_times = (data(curr_animal).bhv(curr_session).frame_times( ...
            1:num_frames))*1000;
        
        imaged_lever_idx = ((1:length(lever_active)) >= ...
            frame_times(1)) & ...
            ((1:length(lever_active)) <= ...
            ceil(frame_times(end)+median(diff(frame_times))));
        
        lever_active_imaged = lever_active(imaged_lever_idx);
        lever_position_resample_imaged = lever_force_resample(imaged_lever_idx);
        lever_velocity_resample_imaged = lever_velocity(imaged_lever_idx);
        lever_speed_resample_imaged = lever_velocity_envelope_smooth(imaged_lever_idx);
        
        relative_frame_times = ceil(frame_times - frame_times(1));
        lever_time_frames = nan(1,length(lever_position_resample_imaged));
        for curr_frame = 1:length(relative_frame_times);
            if curr_frame ~= length(relative_frame_times)
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    relative_frame_times(curr_frame+1)) = curr_frame;
            else
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    end) = curr_frame;
            end
        end
              
        % Split lever by frames
        frame_boundaries = [0,find(lever_time_frames(2:end) ~= ...
            lever_time_frames(1:end-1)),length(lever_time_frames)];
        lever_position_frame_split = mat2cell(lever_position_resample_imaged, ...
            diff(frame_boundaries)',1);
        lever_velocity_frame_split = mat2cell(lever_velocity_resample_imaged, ...
            diff(frame_boundaries)',1);
        lever_speed_frame_split = mat2cell(lever_speed_resample_imaged, ...
            diff(frame_boundaries)',1);
        
        lever_active_split = mat2cell(lever_active_imaged, ...
            diff(frame_boundaries)',1);        

        % Get average lever parameter in each frame
        lever_position_frame = cellfun(@mean,lever_position_frame_split);
        lever_velocity_frame = cellfun(@mean,lever_velocity_frame_split);
        lever_speed_frame = cellfun(@mean,lever_speed_frame_split);
        lever_active_frames = cellfun(@any,lever_active_split);

        % Get median lever of entire day
        median_lever = median(lever_force_smooth);
        
        % Store lever position, velocity, movement in frames
        analysis(curr_animal).lever(curr_session).lever_position_frames = ...
            lever_position_frame-median_lever;
        
        analysis(curr_animal).lever(curr_session).lever_velocity_frames = ...
            lever_velocity_frame;
        
        analysis(curr_animal).lever(curr_session).lever_speed_frames = ...
            lever_speed_frame;
        
        analysis(curr_animal).lever(curr_session).lever_move_frames = ...
            lever_active_frames;
        
        %% 2) Get values related to rewarded movements
                        
        % Initialize
        analysis(curr_animal).lever(curr_session).reaction_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).move_to_reward_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement = cell(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement_fixtime = cell(num_trials,1);
        
        curr_move_onset_frames = nan(num_trials,1);
        curr_move_offset_frames = nan(num_trials,1);
        
        % Only go through trials which are rewarded and imaged
        for curr_trial = intersect(find(rewarded_trials)',data(curr_animal).bhv(curr_session).trial_list(:,2)');
            
            % Get the offset between xsg and bhv
            curr_list_idx = data(curr_animal).bhv(curr_session).trial_list(:,2) == curr_trial;
            curr_xsg_bhv_offset = data(curr_animal).bhv(curr_session).trial_list(curr_list_idx,1) - ...
                data(curr_animal).bhv(curr_session).bhv_times{curr_trial}.states.bitcode(1);
            
            % Get the reward in xsg time (in ms)
            if isfield(data(curr_animal).bhv(curr_session).bhv_times{curr_trial}.states,'reward_delay')
                curr_reward_time = round((data(curr_animal).bhv(curr_session).bhv_times{curr_trial}.states.reward_delay(1) + ...
                    curr_xsg_bhv_offset)*1000);
            else
                curr_reward_time = round((data(curr_animal).bhv(curr_session).bhv_times{curr_trial}.states.reward(1) + ...
                    curr_xsg_bhv_offset)*1000);
            end
            
            % Get the boundaries of the movement that overlaps with reward
            curr_rewarded_movement_onset = curr_reward_time - find( ...
                ~lever_active(curr_reward_time:-1:1),1) + 2;
            curr_rewarded_movement_offset = curr_reward_time + find( ...
                ~lever_active(curr_reward_time:end),1) - 2;
            
            [~,curr_move_onset_frames(curr_trial)] = ...
                min(abs(data(curr_animal).bhv(curr_session).frame_times*1000 - ...
                curr_rewarded_movement_onset));
            
            [~,curr_move_offset_frames(curr_trial)] = ...
                min(abs(data(curr_animal).bhv(curr_session).frame_times*1000 - ...
                curr_rewarded_movement_offset));
            
            % Get reaction time (cue to movement onset, can be negative if
            % movement is animal moved early)
            curr_reaction_time = curr_rewarded_movement_onset - ...
                (data(curr_animal).bhv(curr_session).bhv_times{curr_trial}.states.cue(1) + ...
                curr_xsg_bhv_offset)*1000;
            
            % Get movement onset to reward time
            curr_move_to_reward_time = curr_reward_time - curr_rewarded_movement_onset;
            
            % Get total time of rewarded movement
            curr_rewarded_movement_time = curr_rewarded_movement_offset - ...
                curr_rewarded_movement_onset;
            
            % Get rewarded movement
            curr_rewarded_movement = lever_force_resample(curr_rewarded_movement_onset: ...
                curr_rewarded_movement_offset);
            
            rewarded_movement_back = 1000;
            rewarded_movement_forward = 4000;
            curr_rewarded_movement_fixtime = lever_force_resample( ...
                curr_rewarded_movement_onset-rewarded_movement_back: ...
                curr_rewarded_movement_onset + rewarded_movement_forward);
                      
            % Package lever data into structure
            analysis(curr_animal).lever(curr_session).reaction_time(curr_trial) = curr_reaction_time;
            analysis(curr_animal).lever(curr_session).move_to_reward_time(curr_trial) = curr_move_to_reward_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement_time(curr_trial) = curr_rewarded_movement_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement{curr_trial} = curr_rewarded_movement;
            analysis(curr_animal).lever(curr_session).rewarded_movement_fixtime{curr_trial} = curr_rewarded_movement_fixtime;
        end
        
        % Define cued movement trials where reaction time is positive
        cued_movement_trials = analysis(curr_animal).lever(curr_session).reaction_time > 0;
        
        % Store cued-rewarded trial info
        analysis(curr_animal).lever(curr_session).cued_movement_trials = cued_movement_trials;
        analysis(curr_animal).lever(curr_session).move_onset_frames = curr_move_onset_frames;
        analysis(curr_animal).lever(curr_session).move_offset_frames = curr_move_offset_frames;
    end
    
    disp(['Prepared lever/bhv, animal ' data(curr_animal).animal]);
end
        

%% Align activity to events
clearvars -except data analysis

for curr_animal = 1:length(data)
    
    num_sessions = length(data(curr_animal).im);
    
    for curr_session = 1:num_sessions
        
        curr_move_onset_frames = analysis(curr_animal).lever(curr_session).move_onset_frames;
        curr_move_offset_frames = analysis(curr_animal).lever(curr_session).move_offset_frames;
        
        curr_cued_movement_trials = analysis(curr_animal).lever(curr_session).cued_movement_trials;
        n_cued_movement_trials = sum(curr_cued_movement_trials);
        
        %%% Align to cued-movement onset
        
        % Get frames for cued rewarded movements
        use_move_onset_frames = curr_move_onset_frames( ...
            curr_cued_movement_trials);
        
        % Make matricies for pulling out frames
        surrounding_frames = -90:90;
        analysis(curr_animal).surrounding_frames = surrounding_frames;
        
        move_onset_frames_grab = repmat(use_move_onset_frames,1,length(surrounding_frames)) + ...
            repmat(surrounding_frames,length(use_move_onset_frames),1);
        move_onset_oob = any(move_onset_frames_grab < 1,2);
        move_onset_frames_grab(move_onset_oob,:) = [];
        
        % Pull out frames
        num_rois = size(data(curr_animal).im(curr_session).roi_trace_df,1);
                        
        curr_move_onset_aligned = nan(n_cued_movement_trials, ...
            length(surrounding_frames),num_rois);
        
        curr_move_onset_aligned(~move_onset_oob,:,:) = permute(reshape( ...
            data(curr_animal).im(curr_session).roi_trace_thresh(:,move_onset_frames_grab'), ...
            num_rois,length(surrounding_frames),[]),[3,2,1]);

        
        %%% Align to cued-movement offset
        
        % Get frames for cued rewarded movements
        use_move_offset_frames = curr_move_offset_frames( ...
            curr_cued_movement_trials);
        
        % Make matricies for pulling out frames
        surrounding_frames = -90:90;
        analysis(curr_animal).surrounding_frames = surrounding_frames;
        
        move_offset_frames_grab = repmat(use_move_offset_frames,1,length(surrounding_frames)) + ...
            repmat(surrounding_frames,length(use_move_offset_frames),1);
        move_offset_oob = any(move_offset_frames_grab < 1,2);
        move_offset_frames_grab(move_offset_oob,:) = [];
        
        % Pull out frames
        num_rois = size(data(curr_animal).im(curr_session).roi_trace_df,1);
        
        curr_move_offset_aligned = nan(n_cued_movement_trials, ...
            length(surrounding_frames),num_rois);
        
        curr_move_offset_aligned(~move_offset_oob,:,:) = permute(reshape( ...
            data(curr_animal).im(curr_session).roi_trace_thresh(:,move_offset_frames_grab'), ...
            num_rois,length(surrounding_frames),[]),[3,2,1]);
             

        %%% Store aligned images into structure
        analysis(curr_animal).im(curr_session).move_onset_aligned = curr_move_onset_aligned;
        analysis(curr_animal).im(curr_session).move_offset_aligned = curr_move_offset_aligned;
        
    end
    
    disp(['Prepared aligned activity, animal ' data(curr_animal).animal]);
    
end

disp('Finished all');












        