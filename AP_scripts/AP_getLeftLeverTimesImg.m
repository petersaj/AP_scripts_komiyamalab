%% Get times, states, trials for left lever down

times = bhv.saved_history.RewardsSection_LastTrialEvents;
trial_down = [];
trial_up = [];
states_down = [];
states_up = [];
left_down = [];
left_up = [];
left_times = [];

% start on 2nd trial - 1st can get weird
for curr_trial = 2:length(times);

    % fixes possible bug of putting in old timepoints at the end
    LastGoodTrial = find(diff(times{curr_trial}(:,3)) < 0,1);
    if ~isempty(LastGoodTrial)
        times{curr_trial} = times{curr_trial}(1:LastGoodTrial,:);
    end
    
    % ignore jitter state 
    left_down_indx = find(times{curr_trial}(:,2) == 5 &...
        times{curr_trial}(:,1) ~= 46);
    left_up_indx = find(times{curr_trial}(:,2) == 6 &...
        times{curr_trial}(:,1) ~= 46);
    
    trial_down_curr = curr_trial*ones(length(left_down_indx),1);
    trial_up_curr = curr_trial*ones(length(left_up_indx),1);
    states_down_curr = times{curr_trial}(left_down_indx,1);
    states_up_curr = times{curr_trial}(left_up_indx,1);
    
    states_down = [states_down;states_down_curr];
    states_up = [states_up;states_up_curr];
    trial_down = [trial_down;trial_down_curr];
    trial_up = [trial_up;trial_up_curr];
    left_down = [left_down;times{curr_trial}(left_down_indx,3)];
    left_up = [left_up;times{curr_trial}(left_up_indx,3)];
end

% find bad lever press matches, should make equal lengths
% first for double down/ups
bad_up = [];
bad_down = [];

left_all = [left_down 1*ones(length(left_down),1);...
left_up 2*ones(length(left_up),1)]; % cat all rows, mark down 1 up 2
left_sorted = sortrows(left_all); % sort by time, keep marks
left_alt = diff(left_sorted(:,2));
false_indx = find(left_alt == 0); % look for 2 same in a row
false_dir = left_sorted(false_indx,2);
% 1 is down, 2 is up
for curr_false = 1:length(false_dir);
    if false_dir(curr_false) == 1;
        bad_down = [bad_down; find(left_down == left_sorted(false_indx(curr_false),1),1,'first')];
    elseif false_dir(curr_false) == 2;
        bad_up = [bad_up; find(left_up == left_sorted(false_indx(curr_false),1),1,'last')];
    end
end

% correct ends for start/ending lever position
if left_up(1) - left_down(1) < 0; % lever was down at start
    bad_up = [bad_up;1];
end
if left_up(end) - left_down(end) < 0; % lever down at end
    bad_down = [bad_down;length(left_down)];
end

% eliminate the bad times from doubles or start/end
left_up(bad_up)= [];
left_down(bad_down) = [];
states_up(bad_up) = [];
states_down(bad_down) = [];
trial_up(bad_up) = [];
trial_down(bad_down) = [];

% calculate lever hold times
left_times = left_up - left_down;

% lever times < 2 ms
nonartifact_times = left_times >= 0.002;

left_down = left_down(nonartifact_times);
left_up = left_up(nonartifact_times);
left_times = left_times(nonartifact_times);
states_up = states_up(nonartifact_times);
states_down = states_down(nonartifact_times);
trial_up = trial_up(nonartifact_times);
trial_down = trial_down(nonartifact_times);