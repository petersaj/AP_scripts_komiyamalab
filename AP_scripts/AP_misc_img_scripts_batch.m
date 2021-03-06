%% Get long traces from many sessions
clear all
animal = 35;
listing_cell = {...
    '111102' ...
    '111103' ...
    '111104' ...
    '111106' ...
    '111107' ...
    '111108' ...
    '111109' ...
    '111110' ...
    }'

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    if ~exist('roi_trace_long','var')
        roi_trace_long = [];
        roi_trace_long = ap_getConcatTrace;
        save('Manual_ROI.mat');
    end
end

% not finished yet

%% Get concatenated traces with behavior for batch files
clear all
animal = 38;
listing_cell = {...
    '111102' ...
    '111103' ...
    '111104' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    }'

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI_Adapted110928.mat']);
    unstruct(curr_data);
end

% not finished yet

%% find calcium events for batch files
clear all
animal = 38;
listing_cell = {...
    '111102' ...
    '111103' ...
    '111104' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    }'

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    if ~exist('spike_amplitudes','var')
        [spike_frames spike_amplitudes spike_estimates spike_t0s] = AP_detectEvent(roi_trace_long);
    end
    save([curr_loop_folder filesep 'Manual_ROI.mat']);
    disp(listing_cell{trace_loop_file})
end


%% The point of the following is to find the cross correlations between good cells to get relative timing

% get and concatenate traces from cells normalized by largest spike
good_cells_trace = {};
good_cells = [77 101 104 123 179 180 177 16 36 43 63 64 74 46 54 55 62 56];

listing = dir;

for loop_file = 3:length(listing)
    curr_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP27\' listing(loop_file).name];
    cd(curr_folder)
    clearvars -except listing curr_folder loop_file good_cells good_cells_trace
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load('Manual_ROI_Adapted110929.mat');
    unstruct(curr_data);
    
    roi_trace_norm = roi_trace_long;
    for i = 1:size(roi_trace_norm,1)
        if ~isempty(spike_amplitudes{i})
            max_spike = max(spike_amplitudes{i});
            roi_trace_norm(i,:) = ((roi_trace_norm(i,:))/(max_spike));
        end
    end
    
    good_cells_trace{loop_file-2} = roi_trace_norm(:,:);
end


%% get correlation coefficients of all traces
sqr_num = ceil(sqrt(length(good_cells_trace)));
figure;
for i = 1:length(good_cells_trace)
    subplot(sqr_num,sqr_num,i);
    imagesc(corrcoef(good_cells_trace{i}'));
    caxis([0 1]);
    title(num2str(i));
end

sample = tril(corrcoef(good_cells_trace{i}'));
num_curves = sum(sum(sample ~= 0));
corrcoef_grid = [];
for i = 1:length(good_cells_trace)
    curr_corrcoef = [];
    curr_corrcoef = tril(corrcoef(good_cells_trace{i}'),-1);
    curr_corrcoef(curr_corrcoef == 0) = [];
    corrcoef_grid(:,i) = curr_corrcoef(:);
end
figure;imagesc(corrcoef_grid);

%% get cross correlations of all traces
spike_bin = {};
for i = 1:size(roi_trace_long,1)
    spike_bin{i} = zeros(size(roi_trace_long,2),1);
    spike_bin{i}(spike_frames{i}(:)) = 1;
end

xcorr_all = zeros(size(roi_trace_long,1));
xcorr_time = zeros(size(roi_trace_long,1));
for a = 1:size(roi_trace_long,1)
    for b = 1:a
        %xcorr_curr = xcorr(spike_bin{a},spike_bin{b},50);
        xcorr_curr = xcorr(roi_trace_long(a,:),roi_trace_long(b,:),50);
        %xcorr_curr = smooth(xcorr_curr,5)*5;
        xcorr_curr_max = max(xcorr_curr);%xcorr_curr(size(roi_trace_long,2));%max(xcorr_curr);
        xcorr_all(a,b) = xcorr_curr_max;
        xcorr_time(a,b) = find(xcorr_curr == xcorr_curr_max,1) - 50;
    end
end

%% create full xcorr grid, just add tril to triu, sort into clusters
xcorr_grid = xcorr_grid + tril(xcorr_grid,-1)';
xcorr_grid(logical(eye((size(xcorr_grid))))) = 0;% get rid of diagonal
for i = 1:size(xcorr_grid,1)
    xcorr_grid(i,:) = xcorr_grid(i,:)./max(abs(xcorr_grid(i,:)));
end
[temp1 centers] = kmeans(xcorr_grid,5);
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

% sort roi_trace_long by clusters
% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    temp_trace(i,:) = roi_trace_long(sort_temp1(i,2),:);
end

% also try doing this while keeping amplitude, checking the difference




%% THIS IS TO CREATE FIGURE OF AVG REWARD TRACE OVER DAYS
clear all

trace_days.animal = 35;
trace_days.plot_cells_all = [1];
trace_days.save_flag = 0;

%trace_days.listing_cell = {110917 110918 110919 110920 110922 110923 110924 110926 110928 110929 111007 111010 111011 111012};
%trace_days.listing_cell = cellfun(@(x) num2str(x),trace_days.listing_cell,'UniformOutput',0);
trace_days.listing_cell = {...
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
    }'

trace_days.roi_filename = 'Manual_ROI.mat';

trace_days.reward_trace_days = {};
trace_days.reward_trace_days_full = [];
trace_days.plot_cells = [];

for all_cell_trace_loop = 1:length(trace_days.plot_cells_all)
    
    trace_days.plot_cells = trace_days.plot_cells_all(all_cell_trace_loop);
    
    for trace_loop_file = 1:length(trace_days.listing_cell)
        
        clearvars -except trace_days trace_loop_file
        curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' ...
            num2str(trace_days.animal) '\' trace_days.listing_cell{trace_loop_file}];
        % load into structure then unpack to avoid overwriting metaloop vars
        curr_data = load([curr_loop_folder filesep trace_days.roi_filename], ...
            'roi_trace_long');
        unstruct(curr_data);
        
        % load behavior data
        curr_data = load([curr_loop_folder filesep 'bhv.mat'], ...
            'frames_reward_all');
        unstruct(curr_data);
        
        frame_sample = 50;
        reward_trace_concat = zeros(length(frames_reward_all),frame_sample*2+1,length(trace_days.plot_cells));
        
        %for curr_mod = 1:size(roi_trace_long,1)
        for curr_mod = 1:length(trace_days.plot_cells);
            
            curr_roi = trace_days.plot_cells(curr_mod);
            %curr_roi = curr_mod;
            
            % get traces time-locked to events
            reward_trace = zeros(length(frames_reward_all),frame_sample*2+1);
            for reward_num = 1:length(frames_reward_all);
                if frames_reward_all(reward_num)-frame_sample > 1 && ...
                        frames_reward_all(reward_num)+frame_sample < size(roi_trace_long,2);
                    reward_trace(reward_num,:) = roi_trace_long...
                        (curr_roi,frames_reward_all(reward_num)-frame_sample:frames_reward_all(reward_num)+frame_sample);
                else
                    reward_trace(reward_num,:) = NaN;
                end
            end
            reward_trace_concat(:,:,curr_mod) = reward_trace;
        end
        
        %     %%%% grab only the trials where there was a spike
        %     if exist('spike_amplitudes','var')
        %     frame_back = 10;
        %     frame_forward = 10;
        %     % events around reward/left up
        %     surround_reward = zeros(1,length(frames_reward_all));
        %     for curr_reward = 1:length(frames_reward_all)
        %         curr_frame = frames_reward_all(curr_reward);
        %         curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        %         event_detect = [];
        %         event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        %         if any(event_detect)
        %             % get the largest event which occured around the searched time
        %             event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
        %             surround_reward(curr_reward) = event_detect_amp;
        %         end
        %     end
        %     reward_event = find(surround_reward > 0);
        %     reward_noEvent = find(surround_reward == 0);
        %     %%%%
        %
        %     else
        %         reward_event = [1:size(reward_trace_concat,1)];
        %     end
        
        trace_days.reward_trace_days{trace_loop_file} = reward_trace_concat;
        trace_days.reward_trace_days_full{trace_loop_file} = roi_trace_long(curr_roi,:);
        %trace_days.spike_frames{trace_loop_file} = spike_frames{curr_roi};
        %trace_days.spike_amplitudes{trace_loop_file} = spike_amplitudes{curr_roi};
        
        % Get number of times there was a spike
        if exist('spike_amplitudes','var')
            frame_back = 10;
            frame_forward = 10;
            % events around reward/left up
            surround_reward = zeros(1,length(frames_reward_all));
            for curr_reward = 1:length(frames_reward_all)
                curr_frame = frames_reward_all(curr_reward);
                curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
                event_detect = [];
                event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
                if any(event_detect)
                    % get the largest event which occured around the searched time
                    event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
                    surround_reward(curr_reward) = event_detect_amp;
                end
            end
            trace_days.num_reward_event{trace_loop_file} = length(find(surround_reward > 0));
            trace_days.num_reward_noEvent{trace_loop_file} = length(find(surround_reward == 0));
        end
        
    end
    clearvars -except save_flag trace_days
    
    % plot the traces over days
    
    % get ylimit
    ymax = max(cellfun(@(x) max(nanmean(x)+1.96*(nanstd(x)/sqrt(size(x,1)))),trace_days.reward_trace_days))+0.1;
    ymin = min(cellfun(@(x) min(nanmean(x)-1.96*(nanstd(x)/sqrt(size(x,1)))),trace_days.reward_trace_days));
    
    % plot all long traces
    figure('Name',num2str(trace_days.plot_cells));
    hold on
    for i = 1:length(trace_days.reward_trace_days_full);
        subplot(length(trace_days.reward_trace_days_full),1,i)
        hold on
        plot(trace_days.reward_trace_days_full{i})
        %plot(trace_days.spike_frames{i},trace_days.spike_amplitudes{i},'.r');
        axis off
        ylim([ymin ymax]);
    end
    
    % plot reward_traces for each day
    
    sqr_num = ceil(sqrt(length(trace_days.listing_cell)+1));
    figure('Name',num2str(trace_days.plot_cells)); hold on
    for i = 1:length(trace_days.reward_trace_days)
        subplot(sqr_num,sqr_num,i);
        imagesc(trace_days.reward_trace_days{i});
        colormap(gray);
        hold on
        line([size(trace_days.reward_trace_days{i},2)/2+1 size(trace_days.reward_trace_days{i},2)/2+1],...
            ylim,'linestyle','--','color','k');
        axis off
        title(trace_days.listing_cell(i));
    end
    
    % plot mean traces for each day
    sqr_num = ceil(sqrt(length(trace_days.listing_cell)+1));
    figure('Name',num2str(trace_days.plot_cells)); hold on
    for i = 1:length(trace_days.reward_trace_days)
        subplot(sqr_num,sqr_num,i);
        plot(nanmean(trace_days.reward_trace_days{i}),'linewidth',2,'color','k')
        hold on
        nanstd_reward_trace_u = nanmean(trace_days.reward_trace_days{i})+...
            nanstd(trace_days.reward_trace_days{i});
        nanstd_reward_trace_l = nanmean(trace_days.reward_trace_days{i})-...
            nanstd(trace_days.reward_trace_days{i});
        nansem_reward_trace_u = nanmean(trace_days.reward_trace_days{i}) +...
            1.96*((nanstd(trace_days.reward_trace_days{i})/sqrt...
            (size(trace_days.reward_trace_days{i},1))));
        nansem_reward_trace_l = nanmean(trace_days.reward_trace_days{i}) -...
            1.96*((nanstd(trace_days.reward_trace_days{i})/sqrt...
            (size(trace_days.reward_trace_days{i},1))));
        x_reward = [1:1:size(trace_days.reward_trace_days{i},2)];
        jbfill(x_reward,nansem_reward_trace_u,...
            nansem_reward_trace_l,'blue','none',0,0.5);
        set(gca,'xlim',[0 size(trace_days.reward_trace_days{i},2)]);
        ylim([ymin ymax])
        line([size(trace_days.reward_trace_days{i},2)/2+1 size(trace_days.reward_trace_days{i},2)/2+1],...
            ylim,'linestyle','--','color','k');
        axis off
        title(trace_days.listing_cell(i));
        
        if trace_days.save_flag == 1;
            if ~exist(['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(trace_days.animal) filesep 'Traces_Cell' num2str(curr_roi)],'dir')
                mkdir(['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(trace_days.animal) filesep 'Traces_Cell' num2str(curr_roi)])
                cd ..
            end
            save_filename = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(trace_days.animal) filesep 'Traces_Cell' num2str(curr_roi) filesep 'Session' num2str(i)];
            set(gcf,'PaperPositionMode','auto')
            print('-dpng', save_filename)
        end
        
    end
    
    % plot scale bar
    subplot(sqr_num,sqr_num,sqr_num^2);
    line([0 12.6008],[0 0],'linewidth',3,'color','k') % 2 seconds in frames
    line([0 0],[0 0.5],'linewidth',3,'color','k') % 50% DF/F0
    ylim([ymin ymax])
    xlim([-size(trace_days.reward_trace_days{1},2)/2+1 size(trace_days.reward_trace_days{1},2)/2+1])
    axis off
    
end
% frame_sample = 30;
% figure;
% plot(([start_event{:}]-(frame_sample+1))*(1.24*128/1000),'k','linewidth',3);title('Event start time');
% ylabel('Start time after reward')
% xlabel('Session')
% figure; hold on;
% plot([trace_days.num_reward_event{:}],'linewidth',3);
% plot([trace_days.num_reward_noEvent{:}],'--r','linewidth',3);
% title('Number of rewards with spikes (blue) or without spikes (red)')

%% get and create behavioral data from each day
clear all

animal = 38;

listing_cell = {...
    '111102' ...
    '111103' ...
    '111104' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    }'

for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    
    % get the behavior filename - look for 'data_@' start
    dir_currfolder = dir(curr_loop_folder);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [curr_loop_folder filesep cell2mat(dir_filenames(bhv_file == 1))];
    
    if exist([curr_loop_folder filesep 'bhv.mat'],'file')
        continue
    end
    
    % get the XSG filenames - look for (Initials)(Month)(Day)
    xsg_dirname = ['AP' num2str(listing_cell{trace_loop_file}(3:6))];
    dir_xsg = dir([curr_loop_folder filesep xsg_dirname]);
    xsgdir_filenames = {dir_xsg.name};
    % this is incredibly dumb, but that's because matlab didn't make this
    % function flexible enough: you have to look for the beginning of
    % strings, so flip the strings and loop for backwards extension
    xsgdir_filenames_flipped = cellfun(@(x) x(end:-1:1),xsgdir_filenames,'UniformOutput',0);
    xsg_files = strncmp('gsx.',xsgdir_filenames_flipped,4);
    xsg_part_filenames = xsgdir_filenames([xsg_files == 1]);
    xsg_filenames = cellfun(@(x) [curr_loop_folder filesep xsg_dirname filesep x],xsg_part_filenames,'UniformOutput',0);
    
    AP_getBehavior
    
    save([curr_loop_folder filesep 'bhv.mat']);
end

%% plot behavioral data from each day
clear all

animal = 27;

% just another format, quicker to edit
% listing_cell = {110916 110917 110918 110920 110922 110923 110924 110926 110928 110929 111006 111007 111010 111011 111012 111013 111014 111015 111016 111017};
% listing_cell = cellfun(@(x) num2str(x),listing_cell,'UniformOutput',0);

listing_cell = {...
    '110916' ...
    '110917' ...
    '110918' ...
    '110919' ...
    '110920' ...
    '110922' ...
    '110923' ...
    '110924' ...
    '110926' ...
    '110928' ...
    '110929' ...
    '111006' ...
    '111007' ...
    '111010' ...
    '111011' ...
    '111012' ...
    '111013' ...
    }'

reward_trace_days = {};
grab_var = [];
for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'bhv.mat']);
    unstruct(curr_data);
    
    %%%% get ITI times
    times = saved_history.RewardsSection_LastTrialEvents;
    ITI_times = [];
    ITI_trials = [];
    for curr_trial = 1:length(times)
        if ~ismember(curr_trial, trial_notRecorded)
            % ITI is time from state 40 to 41
            begin_40 = [];
            begin_41 = [];
            
            begin_40= find(times{curr_trial}(:,1) == 40,1);
            begin_41= find(times{curr_trial}(:,1) == 41,1);
            ITI_times = [ITI_times times{curr_trial}(begin_41,3)-times{curr_trial}(begin_40,3)];
            ITI_trials = [ITI_trials curr_trial];
        end
    end
    grab_var.ITI_times{trace_loop_file} = [ITI_times];
    
    %%%% get air reaction times
    times = saved_history.RewardsSection_LastTrialEvents;
    odor_rxn_time_curr = [];
    odor_rxn_times = [];
    odor_rxn_trials = [];
    for curr_trial = 1:length(times);
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
    
    grab_var.rxn_times{trace_loop_file} = [odor_rxn_times];
    
    %%%% Get hold time for the press that counted
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
    
    grab_var.hold_cued{trace_loop_file} = [air_triggered_hold];
    
    %%%% Get all hold times
    grab_var.left_times{trace_loop_file} = [left_times];
    
    %%%% Get lever timing within each trial
    % get start times for each trial
    begin_40 = [];
    for curr_trial = 1:length(times);
        curr_trial_40 = [];
        curr_trial_40 = times{curr_trial}(find(times{curr_trial}(:,1) == 40,1),3);
        begin_40 = [begin_40 curr_trial_40];
    end
    
    % get relative timing for lever presses
    left_down_trial_time = [];
    left_down_trials = [];
    for curr_trial = 1:length(times)
        curr_trial_left = find(trial_down == curr_trial);
        for curr_press = 1:length(curr_trial_left);
            left_down_relative_time = left_down(curr_trial_left(curr_press))...
                - begin_40(curr_trial);
            left_down_trial_time = [left_down_trial_time; curr_trial left_down_relative_time];
        end
    end
    grab_var.left_down_trial_time{trace_loop_file} = [left_down_trial_time];
    
    % get number of successfull and missed trials
    reward_trial = [];
    unrewarded_trial = [];
    for curr_trial = 1:length(times)
        check_trial = [];
        check_trial = find((times{curr_trial}(:,1) == 44));
        if isempty(check_trial)
            unrewarded_trial = [unrewarded_trial curr_trial];
        elseif ~isempty(check_trial)
            reward_trial = [reward_trial curr_trial];
        end
    end
    grab_var.reward_trial{trace_loop_file} = reward_trial;
    grab_var.unrewarded_trial{trace_loop_file} = unrewarded_trial;
    
    disp(['Finished ' listing_cell{trace_loop_file} ', ' num2str(trace_loop_file) '/' num2str(length(listing_cell))]);
    
end

clearvars -except grab_var

%%%% plot everything out
% plot number of correct trials
figure;
bar([cellfun(@length, grab_var.reward_trial)' cellfun(@length, grab_var.unrewarded_trial)'],'stacked');
colormap(gray)
xlabel('Session')
ylabel('# Trials')

temp = [cellfun(@length, grab_var.reward_trial)' cellfun(@length, grab_var.unrewarded_trial)'];


% plot reaction times to air
figure; hold on;
box_prep = [];
box_group = [];
for i = 1:length(grab_var.rxn_times)
    box_prep = [box_prep; grab_var.rxn_times{i}];
    box_group = [box_group; ones(size(grab_var.rxn_times{i},1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','traditional','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Reaction time (s)')
xlabel('Session')

% plot cued hold time
figure; hold on;
box_prep = [];
box_group = [];
for i = 1:length(grab_var.hold_cued)
    box_prep = [box_prep; grab_var.hold_cued{i}(:,2)];
    box_group = [box_group; ones(size(grab_var.hold_cued{i}(:,2),1),1)*i];
end
box_group = num2str(box_group);
boxplot(box_prep,box_group,'notch','on','plotstyle','traditional','whisker',0,'colors',[0 0 0],'symbol','.k')
ylabel('Cued hold time (s)')
xlabel('Session')

% plot relative lever times
figure;hold on
for i = 1:length(grab_var.left_down_trial_time)
    subplot(length(grab_var.left_down_trial_time),1,i)
    plot(grab_var.left_down_trial_time{i}(:,2),grab_var.left_down_trial_time{i}(:,1),'k.','MarkerSize',3)
    xlim([0 15]);
end

% make median +- SEM plots
hold_median = [];
hold_sem = [];
for a = 1:length(grab_var.hold_cued)
    hold_median(a) = median(grab_var.hold_cued{a}(:,2));
    hold_sem(a) = std(grab_var.hold_cued{a}(:,2))/sqrt(size(grab_var.hold_cued{a},1));
end
figure; errorbar(hold_median,hold_sem,'linewidth',2,'color','k');
ylabel('Hold time (s)')
xlabel('Session')

rxn_median = [];
rxn_sem = [];
for a = 1:length(grab_var.rxn_times)
    rxn_median(a) = median(grab_var.rxn_times{a});
    rxn_sem(a) = std(grab_var.rxn_times{a})/sqrt(length(grab_var.rxn_times{a}));
end
figure; errorbar(rxn_median,rxn_sem,'linewidth',2,'color','k');
ylabel('Reaction time (s)')
xlabel('Session')

%% grab a variable from each day

animal = 29;

listing_cell = {...
    '110928' ...
    '110929' ...
    '111006' ...
    '111007' ...
    '111010' ...
    '111011' ...
    '111012' ...
    };

reward_trace_days = {};
grab_var = [];
for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI_Adapted110928.mat']);
    unstruct(curr_data);
    
    grab_var = [grab_var;mod_up_cells;mod_down_cells];
    
end

clearvars -except grab_var


%% show up/down/non-modulated cells across days
clear all

animal = 35;

listing_cell = {111101 111102 111103 111104 111105 111106 111107 111108 111109 111110};
listing_cell = cellfun(@(x) num2str(x),listing_cell,'UniformOutput',0);

reward_trace_days = {};
grab_var = [];
figure
for trace_loop_file = 1:length(listing_cell)
    clearvars -except grab_var animal listing_cell trace_loop_file
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid matlabpool overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    curr_data2 = load([curr_loop_folder filesep 'Manual_ROI.roi'],'-mat');
    unstruct(curr_data2);
    if ~exist('mod_up_cells','var')
        continue
    end
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %  Find modulated cells if not %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     if ~exist('mod_up_cells','var')
    %         disp(['Finding modulated cells, ' listing_cell{trace_loop_file} '...'])
    %         frame_back = 10;
    %         frame_forward = 10;
    %         % events around reward/left up
    %         surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    %         for curr_reward = 1:length(frames_reward_all)
    %             curr_frame = frames_reward_all(curr_reward);
    %             curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    %             for curr_roi = 1:length(spike_frames)
    %                 event_detect = [];
    %                 event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
    %                 if any(event_detect)
    %                     % get the largest event which occured around the searched time
    %                     event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
    %                     surround_reward(curr_roi,curr_reward) = event_detect_amp;
    %                 end
    %             end
    %         end
    %
    %         tic
    %         n_sample = 10000;
    %
    %         frame_back = 10;
    %         frame_forward = 10;
    %
    %         surround_random_all = zeros(size(roi_trace_long,1),n_sample);
    %         frames_random = zeros(n_sample, length(frames_reward_all));
    %         % create matrix for random reward times
    %         for num_rand = 1:n_sample;
    %             frames_random(num_rand,:) = ...
    %                 randsample(size(roi_trace_long,2),length(frames_reward_all));
    %         end
    %
    %         try
    %             matlabpool
    %         catch me
    %         end
    %
    %         parfor curr_roi = 1:length(spike_frames)
    %             frames_random_spike = [];
    %             frames_random_spike = repmat(frames_random,[1 1 length(spike_frames{curr_roi})]) ;
    %             % get difference between random frames and spikes
    %             for curr_spike = 1:length(spike_frames{curr_roi});
    %                 frames_random_spike(:,:,curr_spike) = frames_random_spike(:,:,curr_spike) - spike_frames{curr_roi}(curr_spike);
    %             end
    %             % find any of that matrix that are between frame_back and frame_forward
    %             frames_random_spike = frames_random_spike > -frame_back & frames_random_spike < frame_forward;
    %             surround_random_all(curr_roi,:) = sum(any(frames_random_spike,3),2);
    %             disp(['ROI' num2str(curr_roi) '/' num2str(length(spike_frames))]);
    %         end
    %
    %         % make rewarded events binary
    %         surround_reward_binary = surround_reward;
    %         surround_reward_binary(surround_reward_binary > 0) = 1;
    %         surround_reward_binary = sum(surround_reward_binary,2);
    %
    %         p = zeros(size(roi_trace_long,1),1);
    %         % loop through each ROI, get percentile of events around reward from random
    %         for curr_roi = 1:size(roi_trace_long,1)
    %             % tack on test value to end of random distribution
    %             distr_test = [surround_random_all(curr_roi,:) surround_reward_binary(curr_roi)];
    %             % rank all values, nonparametric way to get percentile
    %             distr_rank = tiedrank(distr_test);
    %             % percentile is the ratio of the test value to all total ranks
    %             p(curr_roi) = distr_rank(end)/length(distr_rank);
    %         end
    %
    %         mod_up_cells = [];
    %         mod_down_cells = [];
    %
    %         mod_up_cells = find(p > 0.975);
    %         mod_down_cells = find(p < 0.025);
    %         toc
    %         disp(['Finished finding modulated cells, saving ' listing_cell{trace_loop_file} '...'])
    %         save([curr_loop_folder filesep 'Manual_ROI.mat']);
    %         disp(['Finished modulated cells, ' listing_cell{trace_loop_file}])
    %     end
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sqr_num = ceil(sqrt(length(listing_cell)));
    
    subplot(sqr_num,sqr_num,trace_loop_file);
    
    hold on
    for i = 1:length(polygon.ROI)
        patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
    end
    
    for i = 1:length(mod_up_cells)
        patch(polygon.ROI{mod_up_cells(i)}(:,1),-polygon.ROI{mod_up_cells(i)}(:,2),[0 1 0]);
    end
    
    for i = 1:length(mod_down_cells)
        patch(polygon.ROI{mod_down_cells(i)}(:,1),-polygon.ROI{mod_down_cells(i)}(:,2),[1 0 0]);
    end
    ylim([-128 0]);
    xlim([0 512]);
    title(listing_cell{trace_loop_file})
    disp(['Finished ' listing_cell{trace_loop_file}])
end

clear all

%% Make pie chart of how many cells are modulated

clear all

animal = 35;

listing_cell = {111101 111102 111103 111104 111105 111106 111107 111108 111109 111110};
listing_cell = cellfun(@(x) num2str(x),listing_cell,'UniformOutput',0);

reward_trace_days = {};
grab_var = [];
figure
for trace_loop_file = 1:length(listing_cell)
    clearvars -except grab_var animal listing_cell trace_loop_file
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid matlabpool overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    curr_data2 = load([curr_loop_folder filesep 'Manual_ROI.roi'],'-mat');
    unstruct(curr_data2);
    if ~exist('mod_up_cells','var')
        continue
    end
    sqr_num = ceil(sqrt(length(listing_cell)));
    
    subplot(sqr_num,sqr_num,trace_loop_file);
    
    hold on
    pie([length(mod_down_cells) (length(polygon.ROI)-(length(mod_down_cells)+length(mod_up_cells))) length(mod_up_cells)])
    colormap gray
    title(listing_cell{trace_loop_file})
    disp(['Finished ' listing_cell{trace_loop_file}])
end

clear all


%% Get and plot the times for trials where a cell spiked

clear all

movie_flag = 0;

mod_cells = 61;

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
grab_var.surround_reward_trace = [];
grab_var.surround_reward_trace_all = [];
grab_var.event_detect_t0 = [];
grab_var.event_detect_maxtime = [];
grab_var.session_markers = [];
grab_var.session_markers_all = [];

for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file mod_cells movie_flag
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    curr_data2 = load([curr_loop_folder filesep 'bhv.mat']);
    unstruct(curr_data2);
    
    medians = {};
    sems = {};
    conf95 = {};
    
    %%%% get events around rewards
    frame_back = 20;
    frame_forward = 20;
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        
        % get all of the traces
        if any(curr_frame_search > size(roi_trace_long,2));
            curr_frame_search(curr_frame_search > size(roi_trace_long,2)) = size(roi_trace_long,2);
        end
        grab_var.surround_reward_trace_all = [grab_var.surround_reward_trace_all; roi_trace_long(mod_cells,curr_frame_search)];
        
        if exist('spike_amplitudes','var')
            % get only the traces with events
            event_detect = [];
            event_detect = ismember(spike_t0s{mod_cells},curr_frame_search);
            if any(event_detect)
                grab_var.surround_reward_trace = [grab_var.surround_reward_trace;roi_trace_long(mod_cells,curr_frame_search)];
                % get the time 0 of the trace
                grab_var.event_detect_t0 = [grab_var.event_detect_t0; (min(spike_t0s{mod_cells}(event_detect)) - curr_frame)];
                % get the time of the max
                event_detect_maxtime = find(grab_var.surround_reward_trace(end,:) == max(grab_var.surround_reward_trace(end,:)));
                grab_var.event_detect_maxtime = [grab_var.event_detect_maxtime; event_detect_maxtime];
            end
        end
        
    end
    
    grab_var.session_markers = [grab_var.session_markers; size(grab_var.surround_reward_trace,1)];
    grab_var.session_markers_all = [grab_var.session_markers_all; size(grab_var.surround_reward_trace_all,1)];
    
end

% figure;plot(grab_var.event_detect_t0,'.')
% figure;imagesc(grab_var.surround_reward_trace)
figure;plot(grab_var.event_detect_maxtime,'.')
% insert lines for session markers
for i = 1:length(grab_var.session_markers)
    hold on
    line([grab_var.session_markers(i) grab_var.session_markers(i)],ylim,'color','r','linestyle','--');
end

figure;plot(smooth(grab_var.event_detect_maxtime,50),'.')

% show the traces with events
figure;
imagesc(grab_var.surround_reward_trace)
colormap(gray);
% insert lines for session markers
for i = 1:length(grab_var.session_markers)
    hold on
    line(xlim,[grab_var.session_markers(i) grab_var.session_markers(i)],'color','r','linestyle','--');
end
title('Traces around rewards with events')

% show all traces with rewards
figure;
imagesc(grab_var.surround_reward_trace_all)
colormap(gray);
% insert lines for session markers
for i = 1:length(grab_var.session_markers_all)
    hold on
    line(xlim,[grab_var.session_markers_all(i) grab_var.session_markers_all(i)],'color','r','linestyle','--');
end
title('All traces around rewards')

% if you want to mark where the peaks are
% hold on;
% for i = 1:size(grab_var.surround_reward_trace,1)
%     plot(grab_var.event_detect_maxtime(i),i,'.r')
% end

if movie_flag == 1;
    figure;
    clear movie_frames
    for i = 1:size(grab_var.surround_reward_trace,1)
        plot(grab_var.surround_reward_trace(i,:))
        ylim([0 max(max(grab_var.surround_reward_trace))]);
        ylim([0 max(grab_var.surround_reward_trace(i,:))]);
        xlim([0 size(grab_var.surround_reward_trace,2)]);
        hold on
        line([grab_var.event_detect_t0(i)+frame_back+1 grab_var.event_detect_t0(i)+frame_back+1],ylim,'color','red');
        movie_frames(i) = getframe(gcf);
        hold off
    end
end


%% Get and plot the times for trials where a cell spiked, but find events the dumb way

clear all

movie_flag = 0;

mod_cells = 51;

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
grab_var.surround_reward_trace = [];
grab_var.surround_reward_trace_all = [];
grab_var.event_detect_t0 = [];
grab_var.event_detect_maxtime = [];
grab_var.session_markers = [];
grab_var.session_markers_all = [];

for trace_loop_file = 1:length(listing_cell)
    
    clearvars -except grab_var animal listing_cell trace_loop_file mod_cells movie_flag
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    curr_data2 = load([curr_loop_folder filesep 'bhv.mat']);
    unstruct(curr_data2);
    
    medians = {};
    sems = {};
    conf95 = {};
    
    %%%% find events
    
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
        
        %     % eliminate peaks at file junctions because they're mistimed
        %     if ~isempty(peak_frames)
        %         file_stitch = [];
        %         file_stitch = find(mod(peak_frames-1,movie_frames) == 0);
        %         peak_frames(file_stitch,:) = [];
        %     end
        disp([num2str(i) '/' num2str(num_rois)])
    end
    %%%%
    
    %%%% get events around rewards
    frame_back = 10;
    frame_forward = 30;
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        
        % get all of the traces
        if any(curr_frame_search > size(roi_trace_long,2));
            curr_frame_search(curr_frame_search > size(roi_trace_long,2)) = size(roi_trace_long,2);
        end
        grab_var.surround_reward_trace_all = [grab_var.surround_reward_trace_all; roi_trace_long(mod_cells,curr_frame_search)];
        % get only the traces with events
        event_detect = [];
        event_detect = ismember(spike_frames{mod_cells},curr_frame_search);
        if any(event_detect)
            grab_var.surround_reward_trace = [grab_var.surround_reward_trace;roi_trace_long(mod_cells,curr_frame_search)];
            % get the time 0 of the trace
            grab_var.event_detect_t0 = [grab_var.event_detect_t0; (min(spike_frames{mod_cells}(event_detect)) - curr_frame)];
        end
    end
    
    grab_var.session_markers = [grab_var.session_markers; size(grab_var.surround_reward_trace,1)];
    grab_var.session_markers_all = [grab_var.session_markers_all; size(grab_var.surround_reward_trace_all,1)];
    
end

figure;plot(grab_var.event_detect_t0,'.')
% figure;imagesc(grab_var.surround_reward_trace)
% insert lines for session markers
for i = 1:length(grab_var.session_markers)
    hold on
    line([grab_var.session_markers(i) grab_var.session_markers(i)],ylim,'color','r','linestyle','--');
end

figure;plot(smooth(grab_var.event_detect_t0,50),'.')

% show the traces with events
figure;
imagesc(grab_var.surround_reward_trace)
colormap(gray);
% insert lines for session markers
for i = 1:length(grab_var.session_markers)
    hold on
    line(xlim,[grab_var.session_markers(i) grab_var.session_markers(i)],'color','r','linestyle','--');
end
title('Traces around rewards with events')

%if you want to mark where the peaks are
hold on;
for i = 1:size(grab_var.surround_reward_trace,1)
    plot(grab_var.event_detect_t0(i),i,'.r')
end

% show all traces with rewards
figure;
imagesc(grab_var.surround_reward_trace_all)
colormap(gray);
% insert lines for session markers
for i = 1:length(grab_var.session_markers_all)
    hold on
    line(xlim,[grab_var.session_markers_all(i) grab_var.session_markers_all(i)],'color','r','linestyle','--');
end
title('All traces around rewards')

if movie_flag == 1;
    figure;
    clear movie_frames
    for i = 1:size(grab_var.surround_reward_trace,1)
        plot(grab_var.surround_reward_trace(i,:))
        ylim([0 max(max(grab_var.surround_reward_trace))]);
        ylim([0 max(grab_var.surround_reward_trace(i,:))]);
        xlim([0 size(grab_var.surround_reward_trace,2)]);
        hold on
        line([grab_var.event_detect_t0(i)+frame_back+1 grab_var.event_detect_t0(i)+frame_back+1],ylim,'color','red');
        movie_frames(i) = getframe(gcf);
        hold off
    end
end