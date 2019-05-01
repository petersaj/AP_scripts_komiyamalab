% This is to get data and figures for the PSEM experiments
clear all

animal = 'AP80';

% find all days
data_path = ['/usr/local/lab/People/Andy/Data/' animal];
data_dir = dir(data_path);
data_names = {data_dir.name};
data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
    'UniformOutput',false);
data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
data_folders = data_fullnames(data_isfolders);
data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);

days = {'121202' '121203' '121204' '121205'};

% get behavior from each day
for curr_day = 1:length(days)
    clearvars -except animal days curr_day psem
    
    day = num2str(days{curr_day});
    
    data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
    
    % Get behavior filename
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    bhv_fullfilename = [data_path filesep bhv_filename];
    warning off
    load(bhv_fullfilename)
    warning on
    bhv = saved_history.ProtocolsSection_parsed_events;
    
    %%%% get start time of all trials
    all_bitcode = cellfun(@(x) x.states.bitcode(1), bhv);
    
    %%%% get number of rewarded trials
    all_reward = cellfun(@(x) x.states.reward, bhv,'uni',0);
    reward_trial = cellfun(@(x) ~isempty(x),all_reward);
    psem.freq_reward(curr_day) = sum(reward_trial)/length(bhv);

    %%%% get number of licks
    all_licks = cellfun(@(x) x.pokes.C(:,1), bhv,'uni',0);
    lick_trial = cellfun(@length,all_licks);
    psem.freq_licks(curr_day) = sum(lick_trial)/length(bhv);

    %%%% get total movement of lever
    
    % Get XSG filenames
    xsg_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_xsg = dir(xsg_path);
    dir_xsg_filenames = {dir_currfolder_xsg.name};
    xsg_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_xsg_filenames);
    xsg_filename_all = dir_xsg_filenames(xsg_filename_indx);
    % make sure they're in order
    xsg_filename_all = sort(xsg_filename_all);xsg_filename_all = dir_xsg_filenames(xsg_filename_indx);
    % make sure they're in order
    xsg_filename_all = sort(xsg_filename_all);
    xsg_fullfilenames = cellfun(@(x) [xsg_path filesep x],xsg_filename_all,'UniformOutput',false);
    
    lever_vel_all = cell(size(xsg_fullfilenames));
    lever_time_all = cell(size(xsg_fullfilenames));
    
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];    
        [lever_active lever_force_smooth lever_force_velocity] = ...
            AP_parseLeverMovement(xsg_data);  
        
        lever_vel_all{curr_xsg} = lever_force_velocity;
        xsg_start = all_bitcode(curr_trial_list(1,2)) - ...
            curr_trial_list(1,1);
        % resampled lever force velocity is at 1kHz
        lever_time_all{curr_xsg} = xsg_start:0.001: ...
            (xsg_start-0.001)+length(lever_force_velocity)*0.001;
    end
    
    lever_vel_concat = vertcat(lever_vel_all{:});
    lever_time_concat = horzcat(lever_time_all{:})';
    
    psem.avg_lever(curr_day) = mean(lever_vel_concat);
    
    % break data into 5 minute segments for averaging
    bin_edges = lever_time_all{1}(1):300:lever_time_all{end}(end);
    [n,trial_bin] = histc(all_bitcode,bin_edges);
    [n,lever_bin] = histc(lever_time_concat,bin_edges);

    reward_bin_mean = grpstats(reward_trial, trial_bin,'sum');
    lick_bin_mean = grpstats(lick_trial, trial_bin,'sum');
    lever_bin_mean = grpstats(lever_vel_concat, lever_bin,'sum');
    
    % ignore things that aren't binned
    psem.reward_bin_mean{curr_day} = reward_bin_mean(2:end);
    psem.lick_bin_mean{curr_day} = lick_bin_mean(2:end);
    psem.lever_bin_mean{curr_day} = lever_bin_mean(2:end);
    
    disp(['Finished ' num2str(curr_day/length(days))]);
    
end


%% Plot PSEM data

reward_all = cell(2,1);
lick_all = cell(2,1);
lever_all = cell(2,1);
for mouse = 1:2;
    switch mouse
        case 1
            load('AP80_PSEM')
            time_idx = {[1 2],[1 2],[1 2],[1 2],[1 2],[3 4],[1 2]};
        case 2
            load('AP81_PSEM')
            time_idx = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2]};
    end
    
    reward_all{mouse} = cell2mat(cellfun(@(x,y) x(y),psem.reward_bin_mean,time_idx,'uni',false));
    reward_all{mouse} = reward_all{mouse}/mean(reward_all{mouse}(1:6));
    
    lick_all{mouse} = cell2mat(cellfun(@(x,y) x(y),psem.lick_bin_mean,time_idx,'uni',false));
    lick_all{mouse} = lick_all{mouse}/mean(lick_all{mouse}(1:6));
    
    lever_all{mouse} = cell2mat(cellfun(@(x,y) x(y),psem.lever_bin_mean,time_idx,'uni',false));
    lever_all{mouse} = lever_all{mouse}/mean(lever_all{mouse}(1:6));
    
end

figure; hold on;
reward_concat = cell2mat(reward_all);
errorbar(mean(reward_concat),std(reward_concat)/sqrt(2),'k','linewidth',2);

lick_concat = cell2mat(lick_all);
errorbar(mean(lick_concat),std(lick_concat)/sqrt(2),'r','linewidth',2);

lever_concat = cell2mat(lever_all);
errorbar(mean(lever_concat),std(lever_concat)/sqrt(2),'b','linewidth',2);

legend({'Rewards' 'Licks' 'Lever'},'location','SW')
set(gca,'XTick',[1:7]);
set(gca,'XTickLabel',{'Day 12' 'Day 13' 'Day 14' 'PSEM' 'Saline' 'PSEM' 'Saline'})
























