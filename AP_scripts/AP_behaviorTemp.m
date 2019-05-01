%% Get times, states, trials for left lever down

times = saved_history.RewardsSection_LastTrialEvents;
trial_num = [];
states_down = [];
states_up = [];
left_down = [];
left_up = [];
left_times = [];
for i = 1:length(times)
    
    % fixes possible bug of putting in old timepoints at the end
    LastGoodTrial = find(diff(times{i}(:,3)) < 0,1);
    if ~isempty(LastGoodTrial)
        times{i} = times{i}(1:LastGoodTrial,:);
    end
    
    left_down_indx = find(times{i}(:,2) == 5);
    left_up_indx = find(times{i}(:,2) == 6);
    
    trial_curr = i*ones(length(left_down_indx),1);
    states_down_curr = times{i}(left_down_indx,1);
    states_up_curr = times{i}(left_up_indx,1);
    
    states_down = [states_down;states_down_curr];
    states_up = [states_up;states_up_curr];
    trial_num = [trial_num;trial_curr];
    left_down = [left_down;times{i}(left_down_indx,3)];
    left_up = [left_up;times{i}(left_up_indx,3)];
end

% correct for unpaired
if length(left_up) ~= length(left_down);
    if left_up(1) - left_down(1) > 0; % extra down at end
        left_times = left_up - left_down(1:end-1);
    elseif left_up(1) - left_down(1) < 0; % lever down at start
        left_times = left_up - left_down(2:end);
    end
else
    left_times = left_up - left_down;
end
    
% elimitate lever times > 1 s or < 2 ms
nonartifact_times = left_times < 1 & left_times > 0.002;
left_down = left_down(nonartifact_times);
left_up = left_up(nonartifact_times);
left_times = left_times(nonartifact_times);
states_up = states_up(nonartifact_times);
states_down = states_down(nonartifact_times);
trial_num = trial_num(nonartifact_times);

%%
% 35 = trial reset, 40 = ITI, 41 = Odor, waiting, 42 = Lever held + Odor
% 43 = no air, wait, 44 = reward, 45 = extra ITI

% Get only rewarded trials
rewarded = find(states_up == 43);

%% plot histogram of lever hold times
half = floor(length(left_times))./2;
left_times_early = left_times(1:half);
left_times_late = left_times(half+1:end);

hist(left_times_early,100)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(left_times_late,100)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
ylabel('Frequency');
xlabel('Lever hold time (s)');
title('Change in frequency of hold time across session');
legend('First half of session','Second half of session');

%% plot states over lever press number
down_43 = find(states_down == 43);
down_45 = find(states_down == 45);
down_35 = find(states_down == 35);
down_40 = find(states_down == 40);
down_41 = find(states_down == 41);
down_42 = find(states_down == 42);
figure
hold on;
plot(trial_num(down_45),[1:length(down_45)]./length(down_45),'.')
plot(trial_num(down_35),[1:length(down_35)]./length(down_35),'r.')
plot(trial_num(down_40),[1:length(down_40)]./length(down_40),'g.')
plot(trial_num(down_41),[1:length(down_41)]./length(down_41),'m<')
plot(trial_num(down_42),[1:length(down_42)]./length(down_42),'k.')
plot(trial_num(down_43),[1:length(down_43)]./length(down_43),'c.')
legend('45','35','40','41','42','43','location','SE')
title('States during lever depression')
xlabel('Trial #');
ylabel('Lever press %');

%% plot reaction time from odor delivery to lever press
left_correct_indx = find(states_down == 41);
left_correct_trials = trial_num(left_correct_indx);
[left_correct_trials, left_correct_trials_indx] = unique(left_correct_trials,'first');

times = saved_history.RewardsSection_LastTrialEvents;
trials_begin_41 = [];
for curr_trial = left_correct_trials';
    begin_41 = find(times{curr_trial}(:,1) == 41,1);
    trials_begin_41 = [trials_begin_41; times{curr_trial}(begin_41,3)];
end

odor_rxn_times = left_down(left_correct_indx) - trials_begin_41;
plot(odor_rxn_times,'k.');
title('Reaction time to odor')
xlabel('Lever depression #')
ylabel('Time between odor and lever depression (s)')

%% plot percent press v trial

figure
subplot(2,1,1)
plot(left_down,[1:length(left_down)]./length(left_down),'k.');
title('Lever depression timing')
xlabel('Time (s)')
ylabel('Normalized lever press')
subplot(2,1,2)
left_down_diff = diff(left_down);
plot([1:length(left_down_diff)]./length(left_down_diff),left_down_diff,'k.')
xlabel('Normalized lever press')
ylabel('Time between lever depressions (s)')


