%% Perform analysis for OdorLeverHold

AP_getLeftLeverTimes;


%% plot histogram of lever hold times
half = floor(length(left_times)/2);
left_times_early = left_times(1:half);
left_times_late = left_times(half+1:end);
left_times_early_med = median(left_times_early);
left_times_late_med = median(left_times_late);

halfmax_time = floor(max(left_down)/2);
half_session = find(left_down <= halfmax_time,1,'last');
left_times_first = left_times(1:half_session);
left_times_second = left_times(half_session+1:end);
left_times_first_med = median(left_times_first);
left_times_second_med = median(left_times_second);

figure
subplot(2,4,1)
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

subplot(2,4,2)
hist(left_times_first,100)
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','k','facealpha',0.5)
hold on;
hist(left_times_second,100)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5);
title('Across time');
legend('1/2 time','2/2 time');
% cutoff for displaying, and lines for medians
xlim([0 1]);
line([left_times_first_med left_times_first_med],ylim,'color','r');
line([left_times_second_med left_times_second_med],ylim,'color','b');


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
subplot(2,2,2)
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


%% plot reaction time from odor delivery to lever press

times = saved_history.RewardsSection_LastTrialEvents;
odor_rxn_time_curr = [];
odor_rxn_times = [];
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
    end
end

subplot(4,4,[9])
plot(odor_rxn_times,'k.','MarkerSize',1);
title('Reaction time to odor')
xlabel('Lever press #')
ylabel('Reaction time (s)')
xlim([1 length(odor_rxn_times)])

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

subplot(4,4,[10])
plot(air_triggered_hold(:,1),air_triggered_hold(:,2),'.k','MarkerSize',1);
ylim([0 1])
title('First hold time after air on')
xlabel('Trial #')
ylabel('Hold time (s)')
xlim([1 length(times)]);

%% plot time between rewarded lever presses

rewarded = find(states_up == 43);
subplot(4,4,13)
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
for curr_trial = 1:length(times)
    curr_trial_left = find(trial_down == curr_trial);
    for curr_press = 1:length(curr_trial_left);
        left_down_relative_time = left_down(curr_trial_left(curr_press))...
            - begin_40(curr_trial);
        line([left_down_relative_time left_down_relative_time],...
            [curr_trial+0.1, curr_trial+0.9],...
            'color','black');
        left_down_trial_time = [left_down_trial_time left_down_relative_time];
    end
end

line([5 5], [ylim],'LineStyle','--','color','black')
title('Lever presses of each trial')
xlabel('Time (s)')
ylabel('Trial')
xlim([-1 10])
ylim([1 length(times)])

%% plot raw lever hold timesstd
subplot(4,4,14)
plot(left_times,'k.','MarkerSize',1);
ylabel('Lever hold time (s)');
xlabel('Lever press #');
title('Lever hold times');

%% plot percent press v trial, interpress time

% plot lever press timing in general
subplot(4,2,6)
plot(left_down,[1:length(left_down)],'k.','MarkerSize',1);
title('Lever press timing')
xlabel('Time (s)')
ylabel('Lever press #')
% plot all lever press interpress times
subplot(4,4,15)
left_down_diff = diff(left_down);
plot([1:length(left_down_diff)],left_down_diff,'k.','MarkerSize',1)
title('Time between lever press')
xlabel('Lever press #')
ylabel('Interpress time (s)')
ylim([0 10])
% plot mean and variability of every 1/4 session
subplot(4,4,16)
fourth_trial = floor(length(times)/4);

left_down_1_4 = left_down(find(trial_down <= fourth_trial));
left_down_2_4 = left_down(find(trial_down > fourth_trial & trial_down <= fourth_trial*2));
left_down_3_4 = left_down(find(trial_down > fourth_trial*2 & trial_down <= fourth_trial*3));
left_down_4_4 = left_down(find(trial_down > fourth_trial*3));

% get interpress times for all and only those over 1 second

diff_left_down_1_4 = diff(left_down_1_4);
diff_left_down2_1_4 = diff_left_down_1_4(diff_left_down_1_4 > 1);
diff_left_down_2_4 = diff(left_down_2_4);
diff_left_down2_2_4 = diff_left_down_2_4(diff_left_down_2_4 > 1);
diff_left_down_3_4 = diff(left_down_3_4);
diff_left_down2_3_4 = diff_left_down_3_4(diff_left_down_3_4 > 1);
diff_left_down_4_4 = diff(left_down_4_4);
diff_left_down2_4_4 = diff_left_down_4_4(diff_left_down_4_4 > 1);

median_diff_left_down_1_4 = median(diff_left_down_1_4);
median_diff_left_down2_1_4 = median(diff_left_down2_1_4);
median_diff_left_down_2_4 = median(diff_left_down_2_4);
median_diff_left_down2_2_4 = median(diff_left_down2_2_4);
median_diff_left_down_3_4 = median(diff_left_down_3_4);
median_diff_left_down2_3_4 = median(diff_left_down2_3_4);
median_diff_left_down_4_4 = median(diff_left_down_4_4);
median_diff_left_down2_4_4 = median(diff_left_down2_4_4);

std_diff_left_down_1_4 = std(diff_left_down_1_4);
std_diff_left_down2_1_4 = std(diff_left_down2_1_4);
std_diff_left_down_2_4 = std(diff_left_down_2_4);
std_diff_left_down2_2_4 = std(diff_left_down2_2_4);
std_diff_left_down_3_4 = std(diff_left_down_3_4);
std_diff_left_down2_3_4 = std(diff_left_down2_3_4);
std_diff_left_down_4_4 = std(diff_left_down_4_4);
std_diff_left_down2_4_4 = std(diff_left_down2_4_4);

median_diff = [median_diff_left_down_1_4,median_diff_left_down_2_4,...
    median_diff_left_down_3_4,median_diff_left_down_4_4,];

median2_diff = [median_diff_left_down2_1_4,median_diff_left_down2_2_4,...
    median_diff_left_down2_3_4,median_diff_left_down2_4_4,];

std_diff = [std_diff_left_down_1_4,std_diff_left_down_2_4,...
    std_diff_left_down_3_4,std_diff_left_down_4_4,];

std2_diff = [std_diff_left_down2_1_4,std_diff_left_down2_2_4,...
    std_diff_left_down2_3_4,std_diff_left_down2_4_4,];

plot([1 2 3 4],median_diff)
hold on;
plot([1 2 3 4],median2_diff,'red')


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

figure
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

subplot(3,1,3)
plot(trials_start,[1:length(trials_start)],'k')
hold on
line([0 max(trials_start)],[0 length(trials_start)],'color','red','linestyle','--')
title('Trial times')
xlabel('Time (s)');
ylabel('Trial #');


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
figure
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