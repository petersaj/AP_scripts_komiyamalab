%% Get relationship between event/no event for rewards and presses

mod_cells = 9;
means = {};
sems = {};
conf95 = {};

for curr_mod_cell = 1:length(mod_cells);
    
    %mod_cells = mod_up_cells(curr_mod_cell);
    
    % reaction time from air onset to lever press
    
    times = saved_history.RewardsSection_LastTrialEvents;
    odor_rxn_time_curr = [];
    odor_rxn_times = [];
    odor_rxn_trials = [];
    for curr_trial = reward_trial;
        if ~ismember(curr_trial, trial_notRecorded)
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
    end
    
    % Hold time for the press that counted
    
    times = saved_history.RewardsSection_LastTrialEvents;
    air_triggered_hold = [];
    for curr_trial = reward_trial;
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
    
    % Get number of presses during the trial before the air came on
    AP_getLeftLeverTimes;
    
    times = saved_history.RewardsSection_LastTrialEvents;
    pre_air_presses = [];
    for curr_trial = reward_trial
        begin_40 = [];
        down_41 = [];
        
        begin_40 = find(times{curr_trial}(:,1) == 40,1);
        down_41 = find(times{curr_trial}(:,1) == 41 &...
            times{curr_trial}(:,2) == 5,1);
        
        begin_40 = times{curr_trial}(begin_40,3);
        down_41 = times{curr_trial}(down_41,3);
        
        pre_air_presses_curr = length(find(left_down > begin_40 & ...
            left_down < down_41));
        
        pre_air_presses = [pre_air_presses;pre_air_presses_curr];
    end
    
    
    % Compare behavior from when cell was active around reward time
    
    reward_event = find(surround_reward(mod_cells,:) > 0);
    reward_noEvent = find(surround_reward(mod_cells,:) == 0);
    
    reward_event_rxn = odor_rxn_times(reward_event);
    reward_noEvent_rxn = odor_rxn_times(reward_noEvent);
    
    reward_event_hold = air_triggered_hold(reward_event,2);
    reward_noEvent_hold = air_triggered_hold(reward_noEvent,2);
    
    reward_event_prepresses = pre_air_presses(reward_event);
    reward_noEvent_prepresses = pre_air_presses(reward_noEvent);
    
    % plot means and stds
    means{curr_mod_cell} = [median(reward_event_rxn) median(reward_noEvent_rxn); ...
        median(reward_event_prepresses) median(reward_noEvent_prepresses); ...
        median(reward_event_hold) median(reward_noEvent_hold)];
    
    num_reward_event = length(reward_event);
    num_reward_noEvent = length(reward_noEvent);
    
    sems{curr_mod_cell} = [std(reward_event_rxn)/sqrt(num_reward_event) std(reward_noEvent_rxn)/sqrt(num_reward_noEvent); ...
        std(reward_event_prepresses)/sqrt(num_reward_event) std(reward_noEvent_prepresses)/sqrt(num_reward_noEvent);
        std(reward_event_hold)/sqrt(num_reward_event) std(reward_event_hold)/sqrt(num_reward_noEvent)];
    
    conf95{curr_mod_cell} = sems{curr_mod_cell}.*1.96;
    
end

figure;
subplot(1,3,1); hold on
bar(means{curr_mod_cell}(1,:))
errorb(means{curr_mod_cell}(1,:),sems{curr_mod_cell}(1,:),'top')
colormap(gray)
title('Reaction Time')
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{'Event','No Event'});
ylabel('Seconds')

subplot(1,3,2); hold on
bar(means{curr_mod_cell}(2,:))
errorb(means{curr_mod_cell}(2,:),sems{curr_mod_cell}(2,:),'top')
colormap(gray)
title('Inter-trial Presses')
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{'Event','No Event'});
ylabel('Presses')

subplot(1,3,3); hold on
bar(means{curr_mod_cell}(3,:))
errorb(means{curr_mod_cell}(3,:),sems{curr_mod_cell}(3,:),'top')
colormap(gray)
title('Hold Time')
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{'Event','No Event'});
ylabel('Seconds')

%% Do that, but across days

mod_cells = 1;

animal = 35;

listing_cell = {...
    '111101' ...
    '111102' ...
    '111103' ...
    '111104' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    '111109' ...
    '111110' ...
    };

reward_trace_days = {};
grab_var = [];
for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file mod_cells
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    curr_data2 = load([curr_loop_folder filesep 'bhv.mat']);
    unstruct(curr_data2);
    
    % can't do this if spikes weren't found yet
    if ~exist('spike_amplitudes','var')
        continue
    end
    
    medians = {};
    sems = {};
    conf95 = {};
    
    %%%% get events around rewards
    frame_back = 10;
    frame_forward = 10;
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        for curr_roi = 1:length(spike_t0s)
            event_detect = [];
            event_detect = ismember(spike_t0s{curr_roi},curr_frame_search);
            if any(event_detect)
                % get the largest event which occured around the searched time
                event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
                surround_reward(curr_roi,curr_reward) = event_detect_amp;
            end
        end
    end
    %%%%
    
    for curr_mod_cell = 1:length(mod_cells);
        
        %mod_cells = mod_up_cells(curr_mod_cell);
        
        % reaction time from air onset to lever press
        
        times = saved_history.RewardsSection_LastTrialEvents;
        odor_rxn_time_curr = [];
        odor_rxn_times = [];
        odor_rxn_trials = [];
        for curr_trial = reward_trial;
            if ~ismember(curr_trial, trial_notRecorded)
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
        end
        
        % Hold time for the press that counted
        
        times = saved_history.RewardsSection_LastTrialEvents;
        air_triggered_hold = [];
        for curr_trial = reward_trial;
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
        
        % Get number of presses during the trial before the air came on
        AP_getLeftLeverTimes;
        
        times = saved_history.RewardsSection_LastTrialEvents;
        pre_air_presses = [];
        for curr_trial = reward_trial
            begin_40 = [];
            down_41 = [];
            
            begin_40 = find(times{curr_trial}(:,1) == 40,1);
            down_41 = find(times{curr_trial}(:,1) == 41 &...
                times{curr_trial}(:,2) == 5,1);
            
            begin_40 = times{curr_trial}(begin_40,3);
            down_41 = times{curr_trial}(down_41,3);
            
            pre_air_presses_curr = length(find(left_down > begin_40 & ...
                left_down < down_41));
            
            pre_air_presses = [pre_air_presses;pre_air_presses_curr];
        end
        
        
        % Compare behavior from when cell was active around reward time
        
        reward_event = find(surround_reward(mod_cells,:) > 0);
        reward_noEvent = find(surround_reward(mod_cells,:) == 0);
        
        reward_event_rxn = odor_rxn_times(reward_event);
        reward_noEvent_rxn = odor_rxn_times(reward_noEvent);
        
        reward_event_hold = air_triggered_hold(reward_event,2);
        reward_noEvent_hold = air_triggered_hold(reward_noEvent,2);
        
        reward_event_prepresses = pre_air_presses(reward_event);
        reward_noEvent_prepresses = pre_air_presses(reward_noEvent);
        
        grab_var.event{trace_loop_file} = [reward_event_rxn reward_event_hold reward_event_prepresses];
        grab_var.noEvent{trace_loop_file} = [reward_noEvent_rxn reward_noEvent_hold reward_noEvent_prepresses];
    end
    
    disp(['Finished ' listing_cell{trace_loop_file} ', ' num2str(trace_loop_file) '/' num2str(length(listing_cell))]);
end

%%%% Plot everything

% plot reaction times to air
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,1)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Reaction time (s)')
xlabel('Session')
title('Reaction time to air - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,1)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Reaction time (s)')
xlabel('Session')
title('Reaction time to air - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,1)) - median(grab_var.noEvent{i}(:,1));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,1));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,1));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Reaction time difference (s)')
xlabel('Session')
title('Reaction time difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);

% plot hold times that count
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,2)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Hold time (s)')
xlabel('Session')
title('Hold time to air - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,2)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Hold time (s)')
xlabel('Session')
title('Hold time to air - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,2)) - median(grab_var.noEvent{i}(:,2));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,2));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,2));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Hold time difference (s)')
xlabel('Session')
title('Hold time difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);

% plot number of pre-presses
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,3)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,3)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,3)) - median(grab_var.noEvent{i}(:,3));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,3));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,3));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);


%% Cross-day behavior with cell, but do it while finding spikes the dumb way

mod_cells = 1;

animal = 35;

listing_cell = {...
    '111101' ...
    '111102' ...
    '111103' ...
    '111104' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    '111109' ...
    '111110' ...
    };

reward_trace_days = {};
grab_var = [];
for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file mod_cells
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    curr_data2 = load([curr_loop_folder filesep 'bhv.mat']);
    unstruct(curr_data2);
    
    % get spike information for current roi
    peakstd = [3*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
    peak_back = 5;
    local_max = 20;
    drop_percentage = .5;
    spike_frames = {};
    movie_frames = 1000;
    
    num_rois = size(roi_trace_long,1);
    for i = mod_cells
        % get estimate of baseline variance
        try
            norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',2); % there are two hidden curves
        catch me
            norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',1); % there are two hidden curves
        end
        baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
        % I think this is std? not sure
        curr_trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve));
        minpeak = curr_trace_std*peakstd;
        spike_frames{i} = ap_schmittDetect(roi_trace_long(i,:),minpeak,peak_back,local_max,drop_percentage);
        disp([num2str(i) '/' num2str(num_rois)])
    end
    
    medians = {};
    sems = {};
    conf95 = {};
    
    %%%% get events around rewards
    frame_back = 10;
    frame_forward = 10;
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        for curr_roi = 1:length(spike_frames)
            event_detect = [];
            event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
            if any(event_detect)
                surround_reward(curr_roi,curr_reward) = 1;
            end
        end
    end
    %%%%
    
    for curr_mod_cell = 1:length(mod_cells);
        
        %mod_cells = mod_up_cells(curr_mod_cell);
        
        % reaction time from air onset to lever press
        
        times = saved_history.RewardsSection_LastTrialEvents;
        odor_rxn_time_curr = [];
        odor_rxn_times = [];
        odor_rxn_trials = [];
        for curr_trial = reward_trial;
            if ~ismember(curr_trial, trial_notRecorded)
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
        end
        
        % Hold time for the press that counted
        
        times = saved_history.RewardsSection_LastTrialEvents;
        air_triggered_hold = [];
        for curr_trial = reward_trial;
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
        
        % Get number of presses during the trial before the air came on
        AP_getLeftLeverTimes;
        
        times = saved_history.RewardsSection_LastTrialEvents;
        pre_air_presses = [];
        for curr_trial = reward_trial
            begin_40 = [];
            down_41 = [];
            
            begin_40 = find(times{curr_trial}(:,1) == 40,1);
            down_41 = find(times{curr_trial}(:,1) == 41 &...
                times{curr_trial}(:,2) == 5,1);
            
            begin_40 = times{curr_trial}(begin_40,3);
            down_41 = times{curr_trial}(down_41,3);
            
            pre_air_presses_curr = length(find(left_down > begin_40 & ...
                left_down < down_41));
            
            pre_air_presses = [pre_air_presses;pre_air_presses_curr];
        end
        
        
        % Compare behavior from when cell was active around reward time
        
        reward_event = find(surround_reward(mod_cells,:) > 0);
        reward_noEvent = find(surround_reward(mod_cells,:) == 0);
        
        reward_event_rxn = odor_rxn_times(reward_event);
        reward_noEvent_rxn = odor_rxn_times(reward_noEvent);
        
        reward_event_hold = air_triggered_hold(reward_event,2);
        reward_noEvent_hold = air_triggered_hold(reward_noEvent,2);
        
        reward_event_prepresses = pre_air_presses(reward_event);
        reward_noEvent_prepresses = pre_air_presses(reward_noEvent);
        
        grab_var.event{trace_loop_file} = [reward_event_rxn reward_event_hold reward_event_prepresses];
        grab_var.noEvent{trace_loop_file} = [reward_noEvent_rxn reward_noEvent_hold reward_noEvent_prepresses];
    end
    
    disp(['Finished ' listing_cell{trace_loop_file} ', ' num2str(trace_loop_file) '/' num2str(length(listing_cell))]);
end

%%%% Plot everything

% plot reaction times to air
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,1)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Reaction time (s)')
xlabel('Session')
title('Reaction time to air - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,1)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Reaction time (s)')
xlabel('Session')
title('Reaction time to air - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,1)) - median(grab_var.noEvent{i}(:,1));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,1));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,1));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Reaction time difference (s)')
xlabel('Session')
title('Reaction time difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);

% plot hold times that count
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,2)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Hold time (s)')
xlabel('Session')
title('Hold time to air - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,2)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Hold time (s)')
xlabel('Session')
title('Hold time to air - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,2)) - median(grab_var.noEvent{i}(:,2));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,2));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,2));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Hold time difference (s)')
xlabel('Session')
title('Hold time difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);

% plot number of pre-presses
figure('Name',['ROI ' num2str(mod_cells)]);
hold on;
subplot(3,1,1);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i})
        continue
    end
    box_prep = [box_prep; grab_var.event{i}(:,3)];
    box_group = [box_group; ones(size(grab_var.event{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses - spike')

subplot(3,1,2);
box_prep = [];
box_group = [];
for i = 1:length(grab_var.noEvent)
    if isempty(grab_var.noEvent{i})
        continue
    end
    box_prep = [box_prep; grab_var.noEvent{i}(:,3)];
    box_group = [box_group; ones(size(grab_var.noEvent{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','compact','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses - no spike')

linkaxes;

subplot(3,1,3)
event_median_diff = [];
p = [];
num_rep = 1000;
for i = 1:length(grab_var.event)
    if isempty(grab_var.event{i}) || isempty(grab_var.noEvent{i})
        continue
    end
    event_median_diff(i) = median(grab_var.event{i}(:,3)) - median(grab_var.noEvent{i}(:,3));
    bootstrp_event = [];
    bootstrp_noEvent = [];
    bootstr_diff = [];
    rank = [];
    bootstrp_event = bootstrp(num_rep,@median,grab_var.event{i}(:,3));
    bootstrp_noEvent = bootstrp(num_rep,@median,grab_var.noEvent{i}(:,3));
    bootstrp_diff = bootstrp_event - bootstrp_noEvent;
    rank = [bootstrp_diff;0];
    rank = tiedrank(rank);
    rank = rank(end)/size(rank,1);
    p(i) = rank;
end
plot(event_median_diff,'k','LineWidth',3);
ylabel('Pre-air presses')
xlabel('Session')
title('Pre-air presses difference')
sig_increase_x = find(p < 0.025 | p > 0.975);
curr_ylim = ylim;
sig_increase_y = ones(size(sig_increase_x))*(curr_ylim(2)-0.2);
text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14);
