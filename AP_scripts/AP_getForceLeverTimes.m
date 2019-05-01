%% Get events for force lever
% This is with up-to-date current setup

% Time for events:
% bhv.saved_history.ProtocolsSection_parsed_events{curr_trial}.states. PUT EVENT HERE (1 = starttime, 2 = endtime);

lick_raster_flag = 1;
check_acq_delay_flag = 1;

[bhv_filename bhv_path] = uigetfile('*','Behavior file');
bhv = load([bhv_path bhv_filename],'-MAT');
parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;

num_trials = length(parsed_events);

all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);

% get trial starts and ends
all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);

reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];

% NOTE ON SETUP: C = lickport, L = lever
% Index if/which lever press ended the cue state
lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);

% Get all lick times, rasterplot them
if lick_raster_flag == 1
    lick_times = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    % rasterplot
    figure; hold on
    for i = 1:length(lick_times)
        if ~isempty(lick_times{i}) && reward_trials(i) == 1;
            plot(lick_times{i}-all_reward{i}(1),i,'.k');
            plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
            plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
            plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
        end
    end
    title('Lick Raster')
end

% Get lever force

[acq_filename_all acq_path] = uigetfile('*.xsg','Choose XSG files','Multiselect','on');

lever_force_split = cell(length(acq_filename_all),1);

lever_force = [];
trial_list = [];
lever_force = [];
bitcode_trace = [];
trial_stitch = [];
% At the moment, concatenate all and pretend no interruption
for i = 1:length(acq_filename_all);
    xsg = load([acq_path acq_filename_all{i}],'-MAT');
    
    % Get all trial offsets in samples
    xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
    
    % Get trials in raw samples since started
    curr_trial_list = [];
    curr_trial_list = AP_ReadBitCode([acq_path acq_filename_all{i}]);
    if ~isempty(curr_trial_list)
        curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
        trial_list = [trial_list; curr_trial_list];
    end
    bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
    lever_force_split{i} = xsg.data.acquirer.trace_2;
    lever_force = [lever_force; xsg.data.acquirer.trace_2];
    trial_stitch = [trial_stitch trial_list(end,2)];
end

% put curr_trial_list time in raw samples
trial_list(:,1) = trial_list(:,1)*xsg_sample_rate;

xsg_bhv_offset = zeros(num_trials,1);
for curr_trial = 1:num_trials
    xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
    if ~isempty(xsg_trial_indx)
        xsg_bhv_offset(curr_trial) = ...
            trial_list(xsg_trial_indx,1) - ...
            all_trial_start{curr_trial}(1)*(xsg_sample_rate);
    else
        xsg_bhv_offset(curr_trial) = NaN;
    end
end

% Get lever presses when rewarded
reward_trials_num = find(reward_trials == 1);
lever_surround_sec = 2;
lever_surround_sample = lever_surround_sec*xsg_sample_rate;
lever_force_reward = nan(num_trials,lever_surround_sample+1);
force_point = nan(num_trials,1);
for curr_reward_trial = reward_trials_num';
    % Skip this if trial start wasn't recorded
    if ~isnan(xsg_bhv_offset(curr_reward_trial)) && ...
            ~ismember(curr_reward_trial,trial_stitch)
        % Current xsg trial time
        xsg_trial_indx = find(trial_list(:,2) == curr_reward_trial);
        lever_sample_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(xsg_sample_rate);
        % xsg lever time = relative lever sample trial start + xsg trial start
        lever_sample_xsg = lever_sample_bhv + trial_list(xsg_trial_indx,1);
        
        if round((lever_sample_xsg-(lever_surround_sample/2))) > 0 && ...
                round((lever_sample_xsg+(lever_surround_sample/2)))
        lever_force_reward(curr_reward_trial,:) = lever_force(round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
        force_point(curr_reward_trial) = lever_sample_xsg;
        end
    end
end

figure; hold on;
plot(lever_force_reward');
line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
title('Rewarded lever force')
% figure; imagesc(lever_force_reward);

% % make stacked plot of lever traces
% lever_force_reward_nonan = lever_force_reward;
% lever_force_reward_nonan(any(isnan(lever_force_reward),2),:) = [];
% figure; hold on
% for i = 1:size(lever_force_reward_nonan,1)
%     plot(lever_force_reward_nonan(i,:)-1.37+i,'k')
% end


figure; hold on;
lever_mean = nanmean(lever_force_reward);
lever_std = nanstd(lever_force_reward);
plot(nanmean(lever_force_reward));
x = 1:size(lever_force_reward,2);
%jbfill(x,lever_mean+lever_std,lever_mean-lever_std)
plot(lever_mean+lever_std,'--r');
plot(lever_mean-lever_std,'--r');
title('Mean rewarded lever force')

% Dot product/similarity stuff

% % plot all the dot products - does this make sense?
% lever_force_reward_nonan = lever_force_reward;
% lever_force_reward_nonan(any(isnan(lever_force_reward),2),:) = [];
% lever_force_reward_nonan_dot = lever_force_reward_nonan*lever_force_reward_nonan';
% figure;imagesc(lever_force_reward*lever_force_reward');
% figure; imagesc(lever_force_reward_nonan_dot);
%
% % plot the dot products with the average
% lever_force_reward_avg = nanmean(lever_force_reward);
% lever_force_reward_nonan_avg = lever_force_reward_nonan*lever_force_reward_avg';
% figure;plot(lever_force_reward_nonan_avg);
%
% % this is the normalized dot I did for oopsi, I think it makes sense...
% n_best = lever_force_reward_nonan;
% % dot product of lever forces
% spike_dot = [];
% spike_dot = n_best*n_best';
% % dot products of lever force on a trial with itself
% self_dot = [];
% self_dot = diag(spike_dot);
% % make rows of the self dots
% self_dot_rep = [];
% self_dot_rep = repmat(self_dot,1,length(self_dot));
% % get the averages of all self dot combinations
% dot_avg = [];
% dot_avg = (self_dot_rep+self_dot_rep')/2;
% % get area under curve for lever forces
% spike_integral = [];
% spike_integral = sum(n_best,2)/size(n_best,2);
% % get expected
% spike_expected = [];
% spike_expected = size(n_best,2)*(spike_integral*spike_integral');
%
% norm_spike_dot = [];
% norm_spike_dot = (spike_dot-spike_expected)./dot_avg;
% figure;imagesc(norm_spike_dot);
%
% % plot dot products normalized by average pairs of self dot products
% figure;imagesc(spike_dot./dot_avg);

% % CHECK ACQ DELAY
% if check_acq_delay_flag == 1
%     
%     baseline = 1.37;
%     thresh_l = baseline + 0.3;
%     thresh_u = baseline + 0.5;
%     rise_time = 100;
%     cross_thresh = nan(num_trials,1);
%     
%     for curr_reward_trial = reward_trials_num';
%         % Skip this if trial start wasn't recorded
%         if ~isnan(xsg_bhv_offset(curr_reward_trial))
%             % Current xsg trial time
%             xsg_trial_indx = find(trial_list(:,2) == curr_reward_trial);
%             for j = trial_list(xsg_trial_indx,1):length(lever_force)-200
%                 curr_sample = lever_force(j:j+200);
%                 min_indx = find(curr_sample == min(curr_sample),1);
%                 thresh_indx = find(curr_sample(min_indx:end) > thresh_u,1);
%                 if ~isempty(thresh_indx)
%                     thresh_indx = thresh_indx + min_indx;
%                     if min(curr_sample) < thresh_l;
%                         cross_thresh(curr_reward_trial) = thresh_indx+j;
%                         break
%                     end
%                 else
%                     continue
%                 end
%             end
%         end
%         disp(curr_reward_trial)
%     end
%     
%     % Plot trials aligned by xsg threshold crossing point
%     reward_trials_num = find(reward_trials == 1);
%     lever_surround_sec = 0.5;
%     lever_surround_sample = lever_surround_sec*xsg_sample_rate;
%     lever_force_reward = nan(num_trials,lever_surround_sample+1);
%     for curr_reward_trial = reward_trials_num';
%         % Skip this if trial start wasn't recorded
%         if ~isnan(cross_thresh(curr_reward_trial))
%             lever_sample_xsg = cross_thresh(curr_reward_trial);
%             lever_force_reward(curr_reward_trial,:) = lever_force(round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
%         end
%     end
%     figure; hold on;
%     plot(lever_force_reward');
%     line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
%     title('Rewarded lever force - XSG')
%     
% end


%% Get events for force lever (batch)
% at the moment: the path is set to run on my PC
clear all
mouse = 'AP56';
days = {120316 120317 120318 120319 120320 120321 120322 120324 120325 120326 120327 120328};
days = cellfun(@(x) num2str(x), days,'UniformOutput',false);
save_path = ['Z:\People\Andy\Data\bhv_data\AP56\'];

for curr_day = 1:length(days)
    clearvars -except mouse days save_path curr_day
    lick_raster_flag = 1;
    
    curr_path = ['Z:\Data\ImagingRig3\' days{curr_day} '\' mouse];
    
    % get the behavior filename - look for 'data_@' start
    dir_currfolder = dir(curr_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    if ~any(bhv_file)
        disp(['Warning: No bhv file for ' days{curr_day}]);
        continue
    end
    bhv_filename = [curr_path filesep cell2mat(dir_filenames(bhv_file == 1))];
    warning off;
    bhv = load(bhv_filename,'-MAT');
    warning on;
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get all lick times, rasterplot them
    if lick_raster_flag == 1
        lick_times = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
        % rasterplot
        figure; hold on
        for i = 1:length(lick_times)
            if ~isempty(lick_times{i}) && reward_trials(i) == 1;
                plot(lick_times{i}-all_reward{i}(1),i,'.k');
                plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
                plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
            end
        end
        title([days{curr_day} ': Lick Raster'])
    end
    
    % Get lever force
    % get the XSG filenames - look for (Initials)(Month)(Day)
    xsg_dirname = [curr_path filesep mouse(1:2) '00' mouse(3:4)];
    dir_xsg = dir(xsg_dirname);
    xsgdir_filenames = {dir_xsg.name};
    % this is incredibly dumb, but that's because matlab didn't make this
    % function flexible enough: you have to look for the beginning of
    % strings, so flip the strings and loop for backwards extension
    xsgdir_filenames_flipped = cellfun(@(x) x(end:-1:1),xsgdir_filenames,'UniformOutput',0);
    xsg_files = strncmp('gsx.',xsgdir_filenames_flipped,4);
    xsg_part_filenames = xsgdir_filenames([xsg_files == 1]);
    xsg_filenames = cellfun(@(x) [xsg_dirname filesep x],xsg_part_filenames,'UniformOutput',0);
    % make sure it's in increasing order
    xsg_filenames = sort(xsg_filenames);
    
    lever_force_split = cell(length(xsg_filenames),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    for i = 1:length(xsg_filenames);
        xsg = load([xsg_filenames{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([xsg_filenames{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; xsg.data.acquirer.trace_2];
        trial_stitch = [trial_stitch trial_list(end,2)];
    end
    
    % put curr_trial_list time in raw samples
    trial_list(:,1) = trial_list(:,1)*xsg_sample_rate;
    
    xsg_bhv_offset = zeros(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1)*(xsg_sample_rate);
        else
            xsg_bhv_offset(curr_trial) = NaN;
        end
    end
    
    % Get lever presses when rewarded
    reward_trials_num = find(reward_trials == 1);
    lever_surround_sec = 1;
    lever_surround_sample = lever_surround_sec*xsg_sample_rate;
    lever_force_reward = nan(num_trials,lever_surround_sample+1);
    force_point = nan(num_trials,1);
    for curr_reward_trial = reward_trials_num';
        % Skip this if trial start wasn't recorded
        if ~isnan(xsg_bhv_offset(curr_reward_trial)) & ...
                ~ismember(curr_reward_trial,trial_stitch)
            % Current xsg trial time
            xsg_trial_indx = find(trial_list(:,2) == curr_reward_trial);
            lever_sample_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(xsg_sample_rate);
            % xsg lever time = relative lever sample trial start + xsg trial start
            lever_sample_xsg = lever_sample_bhv + trial_list(xsg_trial_indx,1);
            lever_force_reward(curr_reward_trial,:) = lever_force(round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
            force_point(curr_reward_trial) = lever_sample_xsg;
        end
    end
    
    figure; hold on;
    plot(lever_force_reward');
    line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
    title([days{curr_day} ': Rewarded lever force'])
    
    
    figure; hold on;
    lever_mean = nanmean(lever_force_reward);
    lever_std = nanstd(lever_force_reward);
    plot(nanmean(lever_force_reward));
    x = 1:size(lever_force_reward,2);
    jbfill(x,lever_mean+lever_std,lever_mean-lever_std);
    title([days{curr_day} ': Mean rewarded lever press +/- std'])
    
    % go through open figures, save (based on assumed order/plot)
    open_figs = findall(0,'type','figure');
    open_figs = sort(open_figs);
    for i = 1:length(open_figs)
       if i == 1;
            save_filename = [days{curr_day} '_' mouse '_' 'lick'];
       elseif i == 2;
           save_filename = [days{curr_day} '_' mouse '_' 'all_lever'];
       elseif i == 3;
           save_filename = [days{curr_day} '_' mouse '_' 'avg_lever'];
       end
       saveas(open_figs(i),[save_path save_filename],'fig');
    end
    close all
end
disp('Finished.')