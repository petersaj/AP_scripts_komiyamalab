function analysis = AP_corticospinal_process_bhv_only(animals)

% Initialize saving strucure
analysis = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Get sessions from day folders
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_path_dir = dir(data_path);
    session_folders = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d\d\d','once')), ...
        {data_path_dir.name});
    
    sessions = {data_path_dir(session_folders).name};
    
    for curr_session = 1:length(sessions)
        
        % Get day and corresponding TIFF path
        day = sessions{curr_session};
        tiff_path = [data_path filesep day];
        
        % Get behavior data
        disp('Getting behavior data...')
        curr_data = struct;
        curr_data.bhv = AP_getLeverBehavior_continuous(animal,tiff_path);
        
        % Analyze lever/trial properties
        num_trials = length(curr_data.bhv.bhv_times);
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),curr_data.bhv.bhv_times);
        
        % If optogenetics, get trials with stim on
        opto_trials = cellfun(@(x) isfield(x.states,'opto_on'),curr_data.bhv.bhv_times);
        opto_session = any(opto_trials);
        if opto_session && ~isfield(curr_data.bhv,'opto')
            error('Opto trials but no TTL trace')
        end
        
        analysis(curr_animal).bhv(curr_session).rewarded_trials = rewarded_trials;
        if opto_session
            analysis(curr_animal).bhv(curr_session).opto_trials = opto_trials;
        end
        
        %%% Lever properties
        
        % Initialize
        analysis(curr_animal).lever(curr_session).cue_to_reward_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).reaction_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).move_to_reward_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).move_to_water_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement_time = nan(num_trials,1);
        analysis(curr_animal).lever(curr_session).rewarded_movement = cell(num_trials,1);
        if opto_session
            analysis(curr_animal).lever(curr_session).opto_trial = nan(num_trials,1);
        end
        
        % Parse movement
        disp('Getting lever movement...')
        [lever_active,lever_force_resample] = AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        curr_move_onset_frames = nan(num_trials,1);
        
        % Only go through trials which are rewarded and imaged
        disp('Getting trial properties...')
        for curr_trial = intersect(find(rewarded_trials)',curr_data.bhv.trial_list(:,2)');
            
            % Get the offset between xsg and bhv
            curr_list_idx = find(curr_data.bhv.trial_list(:,2) == curr_trial);
            curr_xsg_bhv_offset = curr_data.bhv.trial_list(curr_list_idx,1) - ...
                curr_data.bhv.bhv_times{curr_trial}.states.bitcode(1);
            
            % If optogenetics, check for TTL during trial
            if opto_session
                curr_cue_on_time_fulltime = round((curr_data.bhv.bhv_times{curr_trial}.states.cue(1) + ...
                    curr_xsg_bhv_offset)*10000);
                curr_cue_off_time_fulltime = round((curr_data.bhv.bhv_times{curr_trial}.states.cue(2) + ...
                    curr_xsg_bhv_offset)*10000);
                
                curr_opto_trial = any(curr_data.bhv.opto(curr_cue_on_time_fulltime: ...
                    curr_cue_off_time_fulltime) > 2.5);
                
                % If mismatch between TTL and trial, error
                if curr_opto_trial ~= opto_trials(curr_trial)
                    error('Opto trial/TTL mismatch')
                end
            end
            
            curr_cue_on_time = round((curr_data.bhv.bhv_times{curr_trial}.states.cue(1) + ...
                curr_xsg_bhv_offset)*1000);
            curr_cue_off_time = round((curr_data.bhv.bhv_times{curr_trial}.states.cue(2) + ...
                curr_xsg_bhv_offset)*1000);
            curr_cue_to_reward_time = curr_cue_off_time - curr_cue_on_time;
            
            % Get the reward in xsg time (in ms)
            if isfield(curr_data.bhv.bhv_times{curr_trial}.states,'reward_delay')
                % If it's the delayed reward paradigm
                curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward_delay(1) + ...
                    curr_xsg_bhv_offset)*1000);
            else
                % If it's the instant reward paradigm
                curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward(1) + ...
                    curr_xsg_bhv_offset)*1000);
            end
            
            curr_water_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward(1) + ...
                    curr_xsg_bhv_offset)*1000);
            
            % Get the boundaries of the movement that overlaps with reward
            curr_rewarded_movement_onset = curr_reward_time - find( ...
                ~lever_active(curr_reward_time:-1:1),1) + 2;
            curr_rewarded_movement_offset = curr_reward_time + find( ...
                ~lever_active(curr_reward_time:end),1) - 2;
            
            % Get reaction time (cue to movement onset, can be negative if
            % movement is animal moved early)
            curr_reaction_time = curr_rewarded_movement_onset - ...
                (curr_data.bhv.bhv_times{curr_trial}.states.cue(1) + ...
                curr_xsg_bhv_offset)*1000;
            
            % Get movement onset to reward time
            curr_move_to_reward_time = curr_reward_time - curr_rewarded_movement_onset;
            
            % Get movement onset to water time (in case theres delay)
            curr_move_to_water_time = curr_water_time - curr_rewarded_movement_onset;
            
            % Get total time of rewarded movement
            curr_rewarded_movement_time = curr_rewarded_movement_offset - ...
                curr_rewarded_movement_onset;
            
            % Get the rewarded movement from onset to 3 seconds after movement
            rewarded_movement_back = 500;
            rewarded_movement_forward = 3000;
            curr_rewarded_movement = lever_force_resample(curr_rewarded_movement_onset-rewarded_movement_back: ...
                curr_rewarded_movement_onset+rewarded_movement_forward);
            
            
            % Package lever data into structure
            analysis(curr_animal).lever(curr_session).cue_to_reward_time(curr_trial) = curr_cue_to_reward_time;
            analysis(curr_animal).lever(curr_session).reaction_time(curr_trial) = curr_reaction_time;
            analysis(curr_animal).lever(curr_session).move_to_reward_time(curr_trial) = curr_move_to_reward_time;
            analysis(curr_animal).lever(curr_session).move_to_water_time(curr_trial) = curr_move_to_water_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement_time(curr_trial) = curr_rewarded_movement_time;
            analysis(curr_animal).lever(curr_session).rewarded_movement{curr_trial} = curr_rewarded_movement;
            if opto_session
                analysis(curr_animal).lever(curr_session).opto_trial(curr_trial) = curr_opto_trial;
            end
            
        end
        
        % Define cued movement trials where reaction time is positive
        cued_movement_trials = analysis(curr_animal).lever(curr_session).reaction_time > 0;
        analysis(curr_animal).lever(curr_session).cued_movement_trials = cued_movement_trials;
            
        disp(['Finished analyzing bhv: ' animal ' day ' day])
        
    end
    
end