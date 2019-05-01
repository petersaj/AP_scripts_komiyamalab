%% Perform analysis for OdorLeverHold

AP_getLeftLeverTimes;


%% plot raw lever hold timesstd
times_fig = figure;
subplot(4,4,[1 2])
plot(left_times,[1:length(left_times)],'k.','MarkerSize',1);
xlabel('Lever hold time (s)');
ylabel('Lever press #');
title('Lever hold times');
ylim([1 length(left_times)])

%% plot histogram of lever hold times
half = floor(length(left_times)/2);
left_times_early = left_times(1:half);
left_times_late = left_times(half+1:end);
left_times_early_med = median(left_times_early);
left_times_late_med = median(left_times_late);

halfmax_time = length(times)/2;
half_session = find(trial_down <= halfmax_time,1,'last');
left_times_first = left_times(1:half_session);
left_times_second = left_times(half_session+1:end);
left_times_first_med = median(left_times_first);
left_times_second_med = median(left_times_second);

subplot(4,4,5)
hist(left_times_early,100)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(left_times_late,100)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
ylabel('Frequency');
xlabel('Lever hold time (s)');
title('Across lever presses');
legend('1/2 presses','2/2 presses');
% cutoff for displaying, and lines for medians
xlim([0 1]);
line([left_times_early_med left_times_early_med],ylim,'color','r');
line([left_times_late_med left_times_late_med],ylim,'color','b');

subplot(4,4,6)
hist(left_times_first,100)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(left_times_second,100)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
title('Across time');
ylabel('Frequency');
xlabel('Lever hold time (s)');
legend('1/2 time','2/2 time');
% cutoff for displaying, and lines for medians
xlim([0 1]);
line([left_times_first_med left_times_first_med],ylim,'color','r');
line([left_times_second_med left_times_second_med],ylim,'color','b');

%% plot reaction time from odor delivery to lever press

times = saved_history.RewardsSection_LastTrialEvents;
odor_rxn_time_curr = [];
odor_rxn_times = [];
odor_rxn_trials = [];
for curr_trial = 1:length(times);
    begin_41 = [];
    down_41 = [];
    
    begin_41 = find(times{curr_trial}(:,1) == 41,1);
    down_41 = find(times{curr_trial}(:,1) == 41 &...
        times{curr_trial}(:,2) == 5,1);
    if ~isempty(down_41);
        odor_rxn_time_curr = times{curr_trial}(down_41,3)...
            - times{curr_trial}(begin_41,3);
        odor_rxn_times = [odor_rxn_times;odor_rxn_time_curr];
        odor_rxn_trials = [odor_rxn_trials; curr_trial];
    end
end

subplot(4,4,[3 4])
plot(odor_rxn_times,[1:length(odor_rxn_times)],'k.','MarkerSize',1);
title('Reaction time to air')
xlabel('Reaction time (s)')
ylabel('Lever press #')
ylim([1 length(odor_rxn_times)])

%% plot histogram of reaction times

half = floor(length(odor_rxn_times)/2);
odor_rxn_times_early = odor_rxn_times(1:half);
odor_rxn_times_late = odor_rxn_times(half+1:end);
odor_rxn_times_early_med = median(odor_rxn_times_early);
odor_rxn_times_late_med = median(odor_rxn_times_late);

halfmax_time = length(times)/2;
half_session = find(odor_rxn_trials <= halfmax_time,1,'last');
odor_rxn_times_first = odor_rxn_times(1:half_session);
odor_rxn_times_second = odor_rxn_times(half_session+1:end);
odor_rxn_times_first_med = median(odor_rxn_times_first);
odor_rxn_times_second_med = median(odor_rxn_times_second);

subplot(4,4,7)
hist(odor_rxn_times_early,50)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(odor_rxn_times_late,50)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
ylabel('Frequency');
xlabel('Air reaction time (s)');
title('Across lever presses');
legend('1/2 presses','2/2 presses');
% cutoff for displaying, and lines for medians
xlim([0 2]);
line([odor_rxn_times_early_med odor_rxn_times_early_med],ylim,'color','r');
line([odor_rxn_times_late_med odor_rxn_times_late_med],ylim,'color','b');

subplot(4,4,8)
hist(odor_rxn_times_first,50)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(odor_rxn_times_second,50)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
title('Across time');
ylabel('Frequency');
xlabel('Air reaction time (s)');
legend('1/2 time','2/2 time');
% cutoff for displaying, and lines for medians
xlim([0 2]);
line([odor_rxn_times_first_med odor_rxn_times_first_med],ylim,'color','r');
line([odor_rxn_times_second_med odor_rxn_times_second_med],ylim,'color','b');


%% plot first hold time after air

times = saved_history.RewardsSection_LastTrialEvents;
air_triggered_hold = [];
for curr_trial = 1:length(times);
    first_down_41 = [];
    first_up_post41 = [];
    curr_air_triggered_hold = [];
    
    first_down_41 = find(times{curr_trial}(:,1) == 41 ...
        & times{curr_trial}(:,2) == 5,1);
    
    first_up_post41 = find((times{curr_trial}(:,1) == 42 ...
        | times{curr_trial}(:,1) == 43) ...
        & times{curr_trial}(:,2) == 6,1);
    
    curr_air_triggered_hold = times{curr_trial}(first_up_post41,3)...
        - times{curr_trial}(first_down_41,3);
    
    if ~isempty(curr_air_triggered_hold);
        air_triggered_hold = [air_triggered_hold;...
            curr_trial curr_air_triggered_hold];
    end
end

subplot(4,4,[9 10])
plot(air_triggered_hold(:,2),air_triggered_hold(:,1),'.k','MarkerSize',1);
xlim([0 1])
title('First hold time after air on')
ylabel('Trial #')
xlabel('Hold time (s)')
ylim([1 length(times)]);

%% plot histogram of first hold time

half = floor(length(air_triggered_hold)/2);
air_triggered_hold_early = air_triggered_hold(1:half,2);
air_triggered_hold_late = air_triggered_hold(half+1:end,2);
air_triggered_hold_early_med = median(air_triggered_hold_early);
air_triggered_hold_late_med = median(air_triggered_hold_late);

halfmax_time = length(times)/2;
half_session = find(air_triggered_hold(:,1) <= halfmax_time,1,'last');
air_triggered_hold_first = air_triggered_hold(1:half_session,2);
air_triggered_hold_second = air_triggered_hold(half_session+1:end,2);
air_triggered_hold_first_med = median(air_triggered_hold_first);
air_triggered_hold_second_med = median(air_triggered_hold_second);

subplot(4,4,13)
bin = max(air_triggered_hold_early)/.02;
hist(air_triggered_hold_early,bin)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
bin = max(air_triggered_hold_late)/.02;
hist(air_triggered_hold_late,bin)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
ylabel('Frequency');
xlabel('Lever hold time (s)');
title('Across lever presses');
legend('1/2 presses','2/2 presses');
% cutoff for displaying, and lines for medians
xlim([0 1]);
line([air_triggered_hold_early_med air_triggered_hold_early_med],ylim,'color','r');
line([air_triggered_hold_late_med air_triggered_hold_late_med],ylim,'color','b');

subplot(4,4,14)
bin = max(air_triggered_hold_first)/.02;
hist(air_triggered_hold_first,bin)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
bin = max(air_triggered_hold_second)/.02;
hist(air_triggered_hold_second,bin)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
title('Across time');
ylabel('Frequency');
xlabel('Lever hold time (s)');
legend('1/2 time','2/2 time');
% cutoff for displaying, and lines for medians
xlim([0 1]);
line([air_triggered_hold_first_med air_triggered_hold_first_med],ylim,'color','r');
line([air_triggered_hold_second_med air_triggered_hold_second_med],ylim,'color','b');

%% plot raster of lever press times for each trial

rewarded = find(states_up == 43);
subplot(4,4,[11 12])
hold on;

% get start times for each trial
begin_40 = [];
for curr_trial = 1:length(times);
    curr_trial_40 = [];
    curr_trial_40 = times{curr_trial}(find(times{curr_trial}(:,1) == 40,1),3);  
    begin_40 = [begin_40 curr_trial_40];
end

% plot raster of lever presses
left_down_trial_time = [];
left_down_trials = [];
for curr_trial = 1:length(times)
    curr_trial_left = find(trial_down == curr_trial);
    for curr_press = 1:length(curr_trial_left);
        left_down_relative_time = left_down(curr_trial_left(curr_press))...
            - begin_40(curr_trial);
        line([left_down_relative_time left_down_relative_time],...
            [curr_trial+0.1, curr_trial+0.9],...
            'color','black');
        left_down_trial_time = [left_down_trial_time; left_down_relative_time];
        left_down_trials = [left_down_trials; ...
            curr_trial*ones(length(left_down_relative_time))];
    end
end

line([5 5], [ylim],'LineStyle','--','color','black')
title('Lever presses of each trial')
xlabel('Time (s)')
ylabel('Trial')
xlim([-1 10])
ylim([1 length(times)])

%% plot histogram of lever press times within trial

half = floor(length(left_down_trial_time)/2);
left_down_trial_time_early = left_down_trial_time(1:half);
left_down_trial_time_late = left_down_trial_time(half+1:end);
left_down_trial_time_early_med = median(left_down_trial_time_early);
left_down_trial_time_late_med = median(left_down_trial_time_late);

halfmax_time = length(times)/2;
half_session = find(left_down_trials <= halfmax_time,1,'last');
left_down_trial_time_first = left_down_trial_time(1:half_session);
left_down_trial_time_second = left_down_trial_time(half_session+1:end);
left_down_trial_time_first_med = median(left_down_trial_time_first);
left_down_trial_time_second_med = median(left_down_trial_time_second);

subplot(4,4,15)
bin = max(left_down_trial_time_early)/.2;
hist(left_down_trial_time_early,bin)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
bin = max(left_down_trial_time_late)/.2;
hist(left_down_trial_time_late,bin)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
ylabel('Frequency');
xlabel('Lever hold time (s)');
title('Across lever presses');
legend('1/2 presses','2/2 presses');
% cutoff for displaying, and lines for medians
xlim([-1 10]);
line([left_down_trial_time_early_med left_down_trial_time_early_med],ylim,'color','r');
line([left_down_trial_time_late_med left_down_trial_time_late_med],ylim,'color','b');

subplot(4,4,16)
bin = max(left_down_trial_time_first)/.2;
hist(left_down_trial_time_first,bin)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
bin = max(left_down_trial_time_second)/.2;
hist(left_down_trial_time_second,bin)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
title('Across time');
ylabel('Frequency');
xlabel('Lever hold time (s)');
legend('1/2 time','2/2 time');
% cutoff for displaying, and lines for medians
xlim([-1 10]);
line([left_down_trial_time_first_med left_down_trial_time_first_med],ylim,'color','r');
line([left_down_trial_time_second_med left_down_trial_time_second_med],ylim,'color','b');


%% plot histogram of rewarded v unrewarded

times = saved_history.RewardsSection_LastTrialEvents;
num_rewarded = 0;
rewarded_trials = [];
for curr_trial = 1:length(times);
    r_indx = [];
    r_indx = find(times{curr_trial}(:,1) == 44,1);
    if ~isempty(r_indx)
        num_rewarded = num_rewarded + 1;
        rewarded_trials = [rewarded_trials curr_trial];
    end
end

unrewarded_trials = setdiff([1:length(times)],rewarded_trials);

inactive_41 = [];
no_41 = [];
for curr_trial = 1:length(times);
    start_41 = [];
    start_41 = find(times{curr_trial}(:,1) == 41 & ...
        times{curr_trial}(:,2) == 7);
    if ~isempty(start_41);
        inactive_41 = [inactive_41 curr_trial];
    end
    
    is_41 = [];
    is_41 = find(times{curr_trial}(:,1) == 41,1);
    if isempty(is_41);
        no_41 = [no_41 curr_trial];
    end
end

trials_fig = figure;
subplot(3,1,1)
plot(rewarded_trials,0,'.g')
hold on
plot(unrewarded_trials,1,'.r');
% got this from earlier process
plot(inactive_41,1,'.k');
ylim([-10 10])
title('Trial Outcomes')
xlabel('Trial #');
xlim([1 length(times)]);

subplot(3,1,2)
hist(unrewarded_trials,50)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','red','EdgeColor','k','facealpha',0.5)
hold on;
hist(rewarded_trials,50)
h1 = findobj(gca,'Type','patch');
set(h1(1),'FaceColor','green','facealpha',0.5);
ylabel('Frequency');
xlabel('Trial #');
xlim([1 length(times)]);
title('Frequency of Correct Trials');
legend('Inorrect','Correct');

%% Get times of trials

trials_start = [];
for curr_trial = 1:length(times);
    curr_start = [];
    curr_start_indx = find(times{curr_trial}(:,1) == 40,1);
    curr_start_time = times{curr_trial}(curr_start_indx,3);
    trials_start = [trials_start; curr_start_time];
end

subplot(3,2,5)
plot(trials_start,[1:length(trials_start)],'k')
hold on
line([0 max(trials_start)],[0 length(trials_start)],'color','red','linestyle','--')
title('Trial times')
xlabel('Time (s)');
ylabel('Trial #');

%% plot percent press v trial, interpress time

% plot lever press timing in general
subplot(3,2,6)
plot(left_down,[1:length(left_down)],'k.','MarkerSize',1);
title('Lever press timing')
xlabel('Time (s)')
ylabel('Lever press #')


%% plot histogram of lever press states

% turn this on or off, haven't been using it
hist_flag = 0;

if hist_flag == 1;
up_43 = trial_up(find(states_up == 43));
down_45 = trial_down(find(states_down == 45));
down_35 = trial_down(find(states_down == 35));
down_40 = trial_down(find(states_down == 40));
down_41 = trial_down(find(states_down == 41));
up_42 = trial_up(find(states_up == 42));
%--------------------------------------------------
% get times where inactive 41 and no 41
times = saved_history.RewardsSection_LastTrialEvents;
inactive_41 = [];
no_41 = [];
for curr_trial = 1:length(times);
    start_41 = [];
    start_41 = find(times{curr_trial}(:,1) == 41 & ...
        times{curr_trial}(:,2) == 7);
    if ~isempty(start_41);
        inactive_41 = [inactive_41 curr_trial];
    end
    
    is_41 = [];
    is_41 = find(times{curr_trial}(:,1) == 41,1);
    if isempty(is_41);
        no_41 = [no_41 curr_trial];
    end
end
%---------------------------------------------------
states_fig = figure;
hold on;
subplot(2,4,1)
num_trials = length(times);
bin = num_trials/10; % group 10 trials per bar
[d45 d45out] = hist(down_45,bin); d45 = bar(d45out,d45,'r')
title('down 45')
subplot(2,4,2)
[d35 d35out] = hist(down_35,bin); d35 = bar(d35out,d35,'r')
title('down 35')
subplot(2,4,3)
[d40 d40out] = hist(down_40,bin); d40 = bar(d40out,d40,'g')
title('down 40')
subplot(2,4,4)
[d41 d41out] = hist(down_41,bin); d41 = bar(d41out,d41,'m')
title('down 41')
subplot(2,4,5)
[u42 u42out] = hist(up_42,bin); u42 = bar(u42out,u42,'c')
title('up 42')
subplot(2,4,6)
[u43 u43out] = hist(up_43,bin); u43 = bar(u43out,u43,'y')
title('up 43')
subplot(2,4,7)
[i41 i41out] = hist(inactive_41,bin); i41 = bar(i41out,i41,'k')
title('inactive 41')
subplot(2,4,8)
[n41 n41out] = hist(no_41,bin); n41 = bar(n41out,n41,'b')
title('no 41')

xlabel('Trial #');
xlim([1 length(times)]);
ylabel('Frequency');
end

%% plot raw states over lever press number
up_43 = find(states_up == 43);
down_45 = find(states_down == 45);
down_35 = find(states_down == 35);
down_40 = find(states_down == 40);
down_41 = find(states_down == 41);
up_42 = find(states_up == 42);
%--------------------------------------------------
% get times where inactive 41 and no 41
times = saved_history.RewardsSection_LastTrialEvents;
inactive_41 = [];
no_41 = [];
for curr_trial = 1:length(times);
    start_41 = [];
    start_41 = find(times{curr_trial}(:,1) == 41 & ...
        times{curr_trial}(:,2) == 7);
    if ~isempty(start_41);
        inactive_41 = [inactive_41 curr_trial];
    end
    
    is_41 = [];
    is_41 = find(times{curr_trial}(:,1) == 41,1);
    if isempty(is_41);
        no_41 = [no_41 curr_trial];
    end
end
%---------------------------------------------------
figure
hold on;
% if ~isempty(down_45) && ~isempty(down_35) && ~isempty(down_40)...
%         && ~isempty(down_41) && ~isempty(up_42)...
%         && ~isempty(up_43) && ~isempty(inactive_41)...
%         && ~isempty(no_41)
    d45 = plot(trial_down(down_45),[1:length(down_45)]./length(down_45),'.');
    d35 = plot(trial_down(down_35),[1:length(down_35)]./length(down_35),'r.');
    d40 = plot(trial_down(down_40),[1:length(down_40)]./length(down_40),'g.');
    d41 = plot(trial_down(down_41),[1:length(down_41)]./length(down_41),'m<');
    u42 = plot(trial_up(up_42),[1:length(up_42)]./length(up_42),'k*');
    u43 = plot(trial_up(up_43),[1:length(up_43)]./length(up_43),'c*');
    i41 = plot(inactive_41,[1:length(inactive_41)]./length(inactive_41),'py');
    %n41 = plot(no_41,[1:length(no_41)]./length(no_41),'pr');
    legend([d45 d35 d40 d41 u42 u43 i41],...
        '45','35','40','41','42 (up)','43 (up)','41 (inact)','location','NW');
    line([0 length(times)],[0 1],'color','k','linestyle','--');
    title('States during lever press')
    xlabel('Trial #');
    ylabel('Normalized lever press');
    xlim([1 length(times)])
% end
