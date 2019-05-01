n_channels = 1; % true at the moment for all motion corrected files
n_frames = 1000; %get this later

bhv = load(bhv_filename,'-MAT');

saved_history = [];
saved_history = bhv.saved_history; % this in for now, somewhere below fix

w = waitbar(0, 'Getting and concatenating events');

left_down_all = [];
left_up_all = [];
reward_all = [];
unrewarded_all = [];
cue_onset_all = [];
frames_left_down_all = [];
frames_left_up_all = [];
frames_reward_all = [];
frames_unrewarded_all = [];
frames_cue_onset_all = [];
timeDiff = [];
reward_trial = [];
unrewarded_trial = [];

for i = 1:length(xsg_filenames);
    
    % get current acq filename
    acq_filename = [xsg_filenames{i}];
    
    % this is hard coded until the header problem is fixed
    % % %     % pull out time for each for each frame from header
    % % %     img_info = imfinfo([img_pathname img_filename]);
    % % %     msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
    % % %     msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
    % % %     msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
    % % %     lines = img_info(1).Height;
    % % %     fpms = msPerLine*lines;
    
    msPerLine = 1.24;
    lines = 128;
    fpms = msPerLine*lines;
    
    % initialize vars
    left_down = [];
    left_up = [];
    reward = [];
    unrewarded = [];
    cue_onset = [];
    
    % Get trials in sec since started
    TrialNumList = TK_ReadBitCode([acq_filename]);
    close(gcf)
    
    if isempty(TrialNumList);
        continue
    end

    % put trial time in ms
    TrialNumList(2,:) = TrialNumList(2,:)*1000;
       
    for trial = TrialNumList(1,:)
        if trial > length(bhv.saved_history.RewardsSection_LastTrialEvents);
            continue
        end
        bhv_current = bhv.saved_history.RewardsSection_LastTrialEvents{trial,1};
        % fixes possible bug of putting in old timepoints at the end
        LastGoodTrial = find(diff(bhv_current(:,3)) < 0,1);
        if ~isempty(LastGoodTrial)
            bhv_current = bhv_current(1:LastGoodTrial,:);
        end
        % get time lag of behavior to acq (t_bhv + timeDiff = t_acq) in ms
        trialStart_bhv = bhv_current(find(bhv_current(:,1) == 101 & bhv_current(:,2) == 0,1,'last'),3)*1000;
        lag_time = TrialNumList(2,find(TrialNumList(1,:) == trial)) - trialStart_bhv;
        timeDiff(trial,1) = lag_time;
        % make a second column for which acquisition you're in to frame shift
        timeDiff(trial,2) = i*ones(length(lag_time),1);
        % time for all lever up and down (left down = 5, left up = 6) in ms
        %left_down = [left_down; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff(trial)];
        %left_up = [left_up; bhv_current(find(bhv_current(:,2) == 6),3)*1000 + timeDiff(trial)];
        reward = [reward; bhv_current(find(bhv_current(:,1) == 44,1),3)*1000 + timeDiff(trial,1)];
        unrewarded = [unrewarded; bhv_current(find(bhv_current(:,1) == 45,1),3)*1000 + timeDiff(trial,1)];
        if ~isempty(find(bhv_current(:,1) == 44,1))
            reward_trial = [reward_trial trial];
        end
        if ~isempty(find(bhv_current(:,1) == 45,1))
            unrewarded_trial = [unrewarded_trial trial];
        end        
        cue_onset = [cue_onset; bhv_current(find(bhv_current(:,1) == 41,1),3)*1000 + timeDiff(trial,1)];
    end
    
    % concatenate all events, multiply by img number to get total time
    
    %left_down_all = [left_down_all; left_down];
    %left_up_all = [left_up_all; left_up];
    unrewarded_all = [unrewarded_all; unrewarded];
    reward_all = [reward_all; reward];
    cue_onset_all = [cue_onset_all; cue_onset];
    
    %frames_left_down = ceil(left_down(~isnan(left_down))/fpms)';
    %frames_left_up = ceil(left_up(~isnan(left_up))/fpms)';
    frames_reward = ceil(reward(~isnan(reward))/fpms)';
    frames_unrewarded = ceil(unrewarded(~isnan(unrewarded))/fpms)';
    frames_cue_onset = ceil(cue_onset(~isnan(cue_onset))/fpms)';
    
    % round up for number of channels
    %frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
    %frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
    frames_cue_onset = (1-mod(frames_cue_onset,n_channels))+frames_cue_onset;
    frames_reward = (1-mod(frames_reward,n_channels))+frames_reward;
    frame_unrewarded = (1-mod(frames_unrewarded,n_channels))+frames_unrewarded;
    
    % Get frame for each event
    %frames_left_down_all = [frames_left_down_all frames_left_down+n_frames*(i-1)];
    %frames_left_up_all = [frames_left_up_all frames_left_up+n_frames*(i-1)];
    frames_reward_all = [frames_reward_all frames_reward+n_frames*(i-1)];
    frames_unrewarded_all = [frames_unrewarded_all frames_unrewarded+n_frames*(i-1)];    
    frames_cue_onset_all = [frames_cue_onset_all frames_cue_onset+n_frames*(i-1)];
    
    waitbar(i/length(xsg_filenames),w,'Getting and concatenating events');
    
    
end

% get all lever timing at once to get rid of artifacts correctly
AP_getLeftLeverTimesImg
% get rid of lever presses that weren't imaged, maybe only at end?
unimaged_presses_down = trial_down > size(timeDiff,1);
trial_down(unimaged_presses_down) = [];
left_down(unimaged_presses_down) = [];
states_down(unimaged_presses_down) = [];

unimaged_presses_up = trial_up > size(timeDiff,1);
trial_up(unimaged_presses_up) = [];
left_up(unimaged_presses_up) = [];
states_up(unimaged_presses_up) = [];

% correct lever times for time delay, keep old ones
left_down_acq = left_down.*1000 + timeDiff(trial_down,1);
trial_down_acq = trial_down;
states_down_acq = states_down;
left_up_acq = left_up.*1000 + timeDiff(trial_up,1);
trial_up_acq = trial_up;
states_up_acq = states_up;
% get rid of lever times when trial not recorded
trial_notRecorded = find(timeDiff(:,1) == 0);
elim_left_down = ismember(trial_down,trial_notRecorded);
elim_left_up = ismember(trial_up,trial_notRecorded);

left_down_acq(elim_left_down) = [];
trial_down_acq(elim_left_down) = [];
states_down_acq(elim_left_down) = [];
left_up_acq(elim_left_up) = [];
trial_up_acq(elim_left_up) = [];
states_up_acq(elim_left_up) = [];
% pick out the left_down that were in the right place (state 41)
left_down_air = left_down_acq(states_down_acq == 41);

% convert lever times to frames
frames_left_down = ceil(left_down_acq(~isnan(left_down_acq))/fpms)';
frames_left_up = ceil(left_up_acq(~isnan(left_up_acq))/fpms)';
% round frames based on channels
frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
% account for concatenation
% get amount of frames to shift based on trial number
frame_shift_left_down = [];
frame_shift_left_up = [];
frame_shift_left_down = (timeDiff(trial_down_acq,2)-1)*n_frames;
frame_shift_left_up = (timeDiff(trial_up_acq,2)-1)*n_frames;
frames_left_down_all = [frames_left_down_all frames_left_down'+frame_shift_left_down];
frames_left_up_all = [frames_left_up_all frames_left_up'+frame_shift_left_up];

% get the frames where left lever down during correct time
frames_left_down_air_all = frames_left_down_all(states_down_acq == 41);

close(w)
cd ..