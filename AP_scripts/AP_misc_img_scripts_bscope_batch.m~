%%%%% NOTE!!!!!!
%
%
%   all of these scripts are a bit outdated, probably don't use




%% Get long trace in batch
% requires NO OTHER TIFF FILES in folder
clear all
animal = 60;
listing_cell = [120409 ...
    120410 120411 120412 120413 120414 120415 120416 120417 120418];

listing_cell = num2str((listing_cell'));

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['/usr/local/lab/People/Andy/Data/AP' num2str(animal) filesep listing_cell(trace_loop_file,:)];
    roi_filename = 'AP60_template_roi';
    
    roi_trace_long = [];
    roi_trace_long_split = [];
    [roi_trace_long roi_trace_long_split] = AP_getConcatTrace_batch(curr_loop_folder,roi_filename);
    save([curr_loop_folder filesep roi_filename '.mat']);
end


%% Get behavior times (commented out version of next cell)

% commented out all movie-related stuff

animal  = 'AP56';

listing_cell = {120316 120317 120318 120319 120320 120321 120322 ...
    120324 120325 120326 120327 120328};
listing_cell = cellfun(@(x) num2str(x),listing_cell,'UniformOutput',0);

bhv_all = struct;

for curr_day = 1:length(listing_cell);
    
    % clean up shop
    clearvars -except animal listing_cell bhv_all curr_day
    
    day = listing_cell{curr_day};
    
    %     % load roi data
    %     roi_struct = [];
    %     roi_path = ['/usr/local/lab/People/Andy/Data/' animal];
    %     roi_struct = load([roi_path filesep day filesep animal '_template_roi.mat']);
    %     unstruct(roi_struct);
    
    data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
    
    % Get behavior
    % [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    
    warning off
    bhv = load([data_path filesep bhv_filename],'-MAT');
    warning on
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_licks = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    all_press = cellfun(@(x) x.pokes.L(:,1), parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get lever force
    
    % find and load the XSG files
    acq_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_acq = dir(acq_path);
    dir_acq_filenames = {dir_currfolder_acq.name};
    acq_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_acq_filenames);
    
    
    acq_filename_all = dir_acq_filenames(acq_filename_indx);
    % make sure they're in order
    acq_filename_all = sort(acq_filename_all);
    
    lever_force_split = cell(length(acq_filename_all),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    disp('ASSUMING: 8 FULL FILES FOR EACH LOOP')
    
    % Get framerate
    
    % Woah woah WOAH woah woah - this isn't always correct! This depends on
    % hitting the 'measure framerate' button at the beginning, and even then
    % it may not be right. Here's what we'll do instead: above, get the
    % framerate by dividing the number of frames by the length of the xsg in
    % seconds. Note that this is not robust to dropped frames.
    % Edit: at the moment, I know how many frames there should be, so I'll just
    % hardcode it.
    % 2nd Edit: This doesn't make sense - the XSG is slightly longer than the
    % acq, the ONLY information is 'framerate' in the header, which is really
    % bullshit that it's not. For now, assume it's right, but fix for AP60
    % which I know didn't have measured framerate for the first couple days
    % disp('ASSUMING 4000 FRAMES per XSG!');
    % framerate = 4000/(mean(cellfun(@length,lever_force_split))/xsg_sample_rate);
    
    % [img_filename img_path] = uigetfile('*.tif','Pick sample file','Multiselect','off');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    
    img_filename = dir_filenames(tiff_filename_indx);
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    %     framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    %     framerate = str2num(img_value{framerate_indx});
    %
    %     if strcmp(animal,'AP60')
    %         framerate = 28.7224;
    %         disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
    %     end
    
    for i = 1:length(acq_filename_all);
        xsg = load([acq_path filesep acq_filename_all{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([acq_path filesep acq_filename_all{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        
        % Cut out pieces of the lever trace that weren't imaged (dropped
        % frames) - Assume they were clipped at the end
        %         curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
        %         curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
        curr_lever_trace = xsg.data.acquirer.trace_2(1:end); % changed to include all from (1:curr_imaged_lever_samples)
        
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; curr_lever_trace];
        trial_stitch = [trial_stitch trial_list(end,2)];
    end
    
    %     % resample lever force to be 1:1 with numframes
    %     numframes = size(roi_trace_long,2);
    %     % to make resampling work: numframes must be even number
    %     % if it's not, cut out one frame at the end
    %     if mod(numframes,2) ~= 0
    %         roi_trace_long(:,end) = [];
    %         numframes = size(roi_trace_long,2);
    %     end
    %     % resample lever force to so # samples = # frames
    %     [n d] = rat(numframes/length(lever_force));
    %     lever_force_resample = resample(lever_force,n,d);
    %     % add or delete last sample if necessary
    %     if length(lever_force_resample) < numframes
    %         lever_force_resample(end+1:numframes) = 0;
    %     elseif length(lever_force_resample) > numframes
    %         lever_force_resample = lever_force_resample(1:numframes);
    %     end
    
    % Get trial offsets (seconds)
    xsg_bhv_offset = nan(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            continue;
        end
    end
    
    %     % Get dispatcher times in frames
    %     cue_frames = [];
    %     reward_frames = [];
    %     lick_frames = [];
    %     for curr_trial = 1:length(xsg_bhv_offset);
    %         if isnan(xsg_bhv_offset(curr_trial))
    %             continue
    %         end
    %
    %         % relative cue time
    %         cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
    %         % relative reward time
    %         if ~isempty(all_reward{curr_trial})
    %             reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
    %         end
    %         % relative lick time
    %         lick_frames = [lick_frames;(all_licks{curr_trial} + xsg_bhv_offset(curr_trial))*framerate];
    %     end
    
    
    %%%% Plot behavior
    % assuming last cell run
    % Reaction times, licking, press stereotypy
    
    % reaction times: reward - cue
    rxn_time = zeros(length(all_cue),1);
    for i = 1:length(all_cue)
        if ~isempty(all_reward{i})
            rxn_time(i) = all_reward{i}(1) - all_cue{i}(1);
        end
    end
    
    % relative lick times to reward
    lick_times = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    for i = 1:length(lick_times)
        if ~isempty(lick_times{i}) && reward_trials(i) == 1;
            lick_times{i} = lick_times{i} - all_reward{i}(1);
        end
    end
    
    % relative press times to reward
    press_times = cellfun(@(x) x.pokes.L(:,1), parsed_events,'UniformOutput',0);
    for i = 1:length(press_times)
        if ~isempty(press_times{i}) && reward_trials(i) == 1;
            press_times{i} = press_times{i} - all_reward{i}(1);
        end
    end
    
    
    % figure; hold on
    % for i = 1:length(lick_times)
    %     if ~isempty(lick_times{i}) && reward_trials(i) == 1;
    %         plot(lick_times{i},i,'.k');
    %         plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
    %         plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
    %         plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
    %     end
    % end
    % title('Lick Raster')
    
    % press stereotypy
    
    % put curr_trial_list time in raw samples
    trial_list_samples = [];
    trial_list_samples(:,1) = trial_list(:,1)*xsg_sample_rate;
    trial_list_samples(:,2) = trial_list(:,2);
    
    xsg_bhv_offset = zeros(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list_samples(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list_samples(xsg_trial_indx,1) - ...
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
    lever_force_diff_reward = nan(num_trials,lever_surround_sample+1);
    force_point = nan(num_trials,1);
    % smooth by 50 ms, take the velocity
    lever_force_smooth = smooth(lever_force,500);
    lever_force_smooth_diff = diff([lever_force_smooth;0]);
    for curr_reward_trial = reward_trials_num';
        % Skip this if trial start wasn't recorded
        if ~isnan(xsg_bhv_offset(curr_reward_trial)) && ...
                ~ismember(curr_reward_trial,trial_stitch)
            % Current xsg trial time
            xsg_trial_indx = find(trial_list_samples(:,2) == curr_reward_trial);
            lever_sample_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(xsg_sample_rate);
            % xsg lever time = relative lever sample trial start + xsg trial start
            lever_sample_xsg = lever_sample_bhv + trial_list_samples(xsg_trial_indx,1);
            
            if round((lever_sample_xsg-(lever_surround_sample/2))) > 0 && ...
                    round((lever_sample_xsg+(lever_surround_sample/2)))
                lever_force_reward(curr_reward_trial,:) = lever_force...
                    (round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
                lever_force_diff_reward(curr_reward_trial,:) = lever_force_smooth_diff...
                    (round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
                
                force_point(curr_reward_trial) = lever_sample_xsg;
            end
        end
    end
    
    % figure; hold on;
    % plot(lever_force_reward');
    % line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
    % title('Rewarded lever force')
    
    avg_lever_force_reward = nanmean(lever_force_reward)';
    lever_force_reward_norm = bsxfun(@times,lever_force_reward,1./sqrt(sum(lever_force_reward.^2,2)));
    lever_force_reward_dot = lever_force_reward_norm*avg_lever_force_reward;
    
    
    
    % save the bhv data from this day
    % figure; hold on
    % for i = 1:length(lick_times)
    %     if ~isempty(lick_times{i}) && reward_trials(i) == 1;
    %         plot(lick_times{i},i,'.k');
    %         plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
    %         plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
    %         plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
    %     end
    % end
    % title('Lick Raster')
    
    bhv_all.lick{curr_day} = lick_times;
    bhv_all.trial_start{curr_day} = all_trial_start;
    bhv_all.trial_end{curr_day} = all_trial_end;
    bhv_all.cue {curr_day}= all_cue;
    bhv_all.reward {curr_day} = all_reward;
    bhv_all.rxn_time{curr_day} = rxn_time;
    bhv_all.lever_force_reward{curr_day} = lever_force_reward;
    bhv_all.lever_force_diff_reward{curr_day} = lever_force_diff_reward;
    bhv_all.press{curr_day} = press_times;
    
    disp(['Finished bhv' num2str(curr_day/length(listing_cell))]);
end


a = cellfun(@transpose, bhv_all.lever_force_reward,'UniformOutput',false);
b = [a{:}];
b = b';

c = bsxfun(@times,b,1./sqrt(sum(b.^2,2)));
d = nanmean(c)';

e = c*d;


mean_trace = cellfun(@nanmean, bhv_all.lever_force_reward,'UniformOutput',false);
mean_trace = cellfun(@transpose, mean_trace,'UniformOutput',false);
mean_trace = [mean_trace{:}]';
mean_trace_norm = bsxfun(@times,mean_trace,1./sqrt(sum(mean_trace.^2,2)));
figure; hold on;

for i = 1:size(mean_trace,1)
    plot(mean_trace(i,:),'color',[0 i/size(mean_trace,1) 0]);
end
plot(nanmean(b),'linewidth',2);
title(animal);

% compare dots of lever press with the day's average
lever_dot = [];
lever_dot_norm = [];
for i = 1:size(mean_trace_norm,1)
    % try for the last day's average
    lever_dot = [lever_dot;bhv_all.lever_force_reward{i}*mean_trace(end,:)'];
    
    temp_norm = bsxfun(@times,bhv_all.lever_force_reward{i},1./sqrt(sum(bhv_all.lever_force_reward{i}.^2,2)));
    lever_dot_norm = [lever_dot_norm;temp_norm*mean_trace_norm(i,:)'];
end

% press raster, for simon
for i = 1:length(bhv_all.press)
    figure; hold on
    press_times = bhv_all.press{i};
    for i = 1:length(press_times)
        if ~isempty(press_times{i}) && reward_trials(i) == 1;
            plot(press_times{i},i,'.k');
            plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
            plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
            plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
        end
    end
    title('Press Raster')
end

%% Get trace and behavior times

% at the moment this does stuff with frames, don't remember exact what

animal  = 'AP61';

listing_cell = {120414  120419  ...
     120426 };
listing_cell = cellfun(@(x) num2str(x),listing_cell,'UniformOutput',0);

bhv_all = struct;

for curr_day = 1:length(listing_cell);
    
    % clean up shop
    clearvars -except animal listing_cell bhv_all curr_day
    
    day = listing_cell{curr_day};
    
%     % load roi data
%     roi_struct = [];
%     roi_path = ['/usr/local/lab/People/Andy/Data/' animal];
%     roi_struct = load([roi_path filesep day filesep animal '_template_roi.mat']);
%     unstruct(roi_struct);
    
    data_path = ['/usr/local/lab/Data/ImagingRig3/120108-120525/' day filesep animal];
    
    % Get behavior
    % [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    
    warning off
    bhv = load([data_path filesep bhv_filename],'-MAT');
    warning on
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_licks = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get lever force
    
    % find and load the XSG files
    acq_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_acq = dir(acq_path);
    dir_acq_filenames = {dir_currfolder_acq.name};
    acq_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_acq_filenames);
    
    
    acq_filename_all = dir_acq_filenames(acq_filename_indx);
    % make sure they're in order
    acq_filename_all = sort(acq_filename_all);
    
    lever_force_split = cell(length(acq_filename_all),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    disp('ASSUMING: 8 FULL FILES FOR EACH LOOP')
    
    % Get framerate
    
    % Woah woah WOAH woah woah - this isn't always correct! This depends on
    % hitting the 'measure framerate' button at the beginning, and even then
    % it may not be right. Here's what we'll do instead: above, get the
    % framerate by dividing the number of frames by the length of the xsg in
    % seconds. Note that this is not robust to dropped frames.
    % Edit: at the moment, I know how many frames there should be, so I'll just
    % hardcode it.
    % 2nd Edit: This doesn't make sense - the XSG is slightly longer than the
    % acq, the ONLY information is 'framerate' in the header, which is really
    % bullshit that it's not. For now, assume it's right, but fix for AP60
    % which I know didn't have measured framerate for the first couple days
    % disp('ASSUMING 4000 FRAMES per XSG!');
    % framerate = 4000/(mean(cellfun(@length,lever_force_split))/xsg_sample_rate);
    
    % [img_filename img_path] = uigetfile('*.tif','Pick sample file','Multiselect','off');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    
    img_filename = dir_filenames(tiff_filename_indx);
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    
    if strcmp(animal,'AP60')
        framerate = 28.7224;
        disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
    end
    
    for i = 1:length(acq_filename_all);
        xsg = load([acq_path filesep acq_filename_all{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([acq_path filesep acq_filename_all{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        
        % Cut out pieces of the lever trace that weren't imaged (dropped
        % frames) - Assume they were clipped at the end
        curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
        curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
        curr_lever_trace = xsg.data.acquirer.trace_2(1:curr_imaged_lever_samples);
        
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; curr_lever_trace];
        trial_stitch = [trial_stitch trial_list(end,2)];
    end
    
    % resample lever force to be 1:1 with numframes
    numframes = size(roi_trace_long,2);
    % to make resampling work: numframes must be even number
    % if it's not, cut out one frame at the end
    if mod(numframes,2) ~= 0
        roi_trace_long(:,end) = [];
        numframes = size(roi_trace_long,2);
    end
    % resample lever force to so # samples = # frames
    [n d] = rat(numframes/length(lever_force));
    lever_force_resample = resample(lever_force,n,d);
    % add or delete last sample if necessary
    if length(lever_force_resample) < numframes
        lever_force_resample(end+1:numframes) = 0;
    elseif length(lever_force_resample) > numframes
        lever_force_resample = lever_force_resample(1:numframes);
    end
    
    % Get trial offsets (seconds)
    xsg_bhv_offset = nan(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            continue;
        end
    end
    
    % Get dispatcher times in frames
    cue_frames = [];
    reward_frames = [];
    lick_frames = [];
    for curr_trial = 1:length(xsg_bhv_offset);
        if isnan(xsg_bhv_offset(curr_trial))
            continue
        end
        
        % relative cue time
        cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        % relative reward time
        if ~isempty(all_reward{curr_trial})
            reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        end
        % relative lick time
        lick_frames = [lick_frames;(all_licks{curr_trial} + xsg_bhv_offset(curr_trial))*framerate];
    end
    
    
    %%%% Plot behavior
    % assuming last cell run
    % Reaction times, licking, press stereotypy
    
    % reaction times: reward - cue
    rxn_time = zeros(length(all_cue),1);
    for i = 1:length(all_cue)
        if ~isempty(all_reward{i})
            rxn_time(i) = all_reward{i}(1) - all_cue{i}(1);
        end
    end
    
    % relative lick times to reward
    lick_times = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    for i = 1:length(lick_times)
        if ~isempty(lick_times{i}) && reward_trials(i) == 1;
            lick_times{i} = lick_times{i} - all_reward{i}(1);
        end
    end
    
    % figure; hold on
    % for i = 1:length(lick_times)
    %     if ~isempty(lick_times{i}) && reward_trials(i) == 1;
    %         plot(lick_times{i},i,'.k');
    %         plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
    %         plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
    %         plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
    %     end
    % end
    % title('Lick Raster')
    
    % press stereotypy
    
    % put curr_trial_list time in raw samples
    trial_list_samples = [];
    trial_list_samples(:,1) = trial_list(:,1)*xsg_sample_rate;
    trial_list_samples(:,2) = trial_list(:,2);
    
    xsg_bhv_offset = zeros(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list_samples(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list_samples(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1)*(xsg_sample_rate);
        else
            xsg_bhv_offset(curr_trial) = NaN;
        end
    end
    
    % Get lever presses when rewarded
    reward_trials_num = find(reward_trials == 1);
    lever_surround_sec = 0.5;
    lever_surround_sample = lever_surround_sec*xsg_sample_rate;
    lever_force_reward = nan(num_trials,lever_surround_sample+1);
    force_point = nan(num_trials,1);
    % smooth by 50 ms, take the velocity
    lever_force_smooth = smooth(lever_force,500);
    lever_force_smooth_diff = diff([lever_force_smooth;0]);
    for curr_reward_trial = reward_trials_num';
        % Skip this if trial start wasn't recorded
        if ~isnan(xsg_bhv_offset(curr_reward_trial)) && ...
                ~ismember(curr_reward_trial,trial_stitch)
            % Current xsg trial time
            xsg_trial_indx = find(trial_list_samples(:,2) == curr_reward_trial);
            lever_sample_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(xsg_sample_rate);
            % xsg lever time = relative lever sample trial start + xsg trial start
            lever_sample_xsg = lever_sample_bhv + trial_list_samples(xsg_trial_indx,1);
            
            if round((lever_sample_xsg-(lever_surround_sample/2))) > 0 && ...
                    round((lever_sample_xsg+(lever_surround_sample/2)))
                lever_force_reward(curr_reward_trial,:) = lever_force_smooth_diff...
                    (round((lever_sample_xsg-(lever_surround_sample/2))):round((lever_sample_xsg+(lever_surround_sample/2))));
                force_point(curr_reward_trial) = lever_sample_xsg;
            end
        end
    end
    
    % figure; hold on;
    % plot(lever_force_reward');
    % line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
    % title('Rewarded lever force')
    
    avg_lever_force_reward = nanmean(lever_force_reward)';
    lever_force_reward_norm = bsxfun(@times,lever_force_reward,1./sqrt(sum(lever_force_reward.^2,2)));
    lever_force_reward_dot = lever_force_reward_norm*avg_lever_force_reward;
    
    
    
    % save the bhv data from this day
    bhv_all.rxn_time{curr_day} = rxn_time;
    bhv_all.lever_force_reward{curr_day} = lever_force_reward;
    
    disp(['Finished bhv' num2str(curr_day/length(listing_cell))]);
end

% 
% a = cellfun(@transpose, bhv_all.lever_force_reward,'UniformOutput',false);
% b = [a{:}];
% b = b';
% 
% c = bsxfun(@times,b,1./sqrt(sum(b.^2,2)));
% d = nanmean(c)';
% 
% e = c*d;
% 
% 
% mean_trace = cellfun(@nanmean, bhv_all.lever_force_reward,'UniformOutput',false);
% mean_trace = cellfun(@transpose, mean_trace,'UniformOutput',false);
% mean_trace = [mean_trace{:}]';
% %mean_trace = bsxfun(@times,mean_trace,1./sqrt(sum(mean_trace.^2,2)));
% figure; hold on;
% 
% for i = 1:size(mean_trace,1)
%     plot(mean_trace(i,:),'color',[0 i/size(mean_trace,1) 0]);
% end
% plot(nanmean(b),'linewidth',2);




%% Get weighted force for a bunch of days

%%%%%% WARNING!!!
%%%%%% As it is, OOPSI picks up way too many movement artifacts

clear all
animal = 'AP60';
listing_cell = [120330 120331 120401 120402 120403 120404 120409 ...
    120410 120411 120412 120413 120414 120415 120416 120417 120418];

listing_cell = num2str((listing_cell'));

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['/usr/local/lab/People/Andy/Data/' animal filesep listing_cell(trace_loop_file,:)];
    data_filename = 'AP60_template_roi.mat';
    
    data_struct = load([curr_loop_folder filesep data_filename]);
    unstruct(data_struct);
    
    %%%%
    %%%% First: get bhv
    %%%%
    day = listing_cell(trace_loop_file,:);
    
    data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
    
    % Get behavior
    % [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    
    warning off
    bhv = load([data_path filesep bhv_filename],'-MAT');
    warning on
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_licks = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get lever force
    
    % find and load the XSG files
    acq_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_acq = dir(acq_path);
    dir_acq_filenames = {dir_currfolder_acq.name};
    acq_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_acq_filenames);
    
    
    acq_filename_all = dir_acq_filenames(acq_filename_indx);
    % make sure they're in order
    acq_filename_all = sort(acq_filename_all);
    
    lever_force_split = cell(length(acq_filename_all),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    disp('ASSUMING: 8 FULL FILES FOR EACH LOOP')
    
    % Get framerate
    
    % Woah woah WOAH woah woah - this isn't always correct! This depends on
    % hitting the 'measure framerate' button at the beginning, and even then
    % it may not be right. Here's what we'll do instead: above, get the
    % framerate by dividing the number of frames by the length of the xsg in
    % seconds. Note that this is not robust to dropped frames.
    % Edit: at the moment, I know how many frames there should be, so I'll just
    % hardcode it.
    % 2nd Edit: This doesn't make sense - the XSG is slightly longer than the
    % acq, the ONLY information is 'framerate' in the header, which is really
    % bullshit that it's not. For now, assume it's right, but fix for AP60
    % which I know didn't have measured framerate for the first couple days
    % disp('ASSUMING 4000 FRAMES per XSG!');
    % framerate = 4000/(mean(cellfun(@length,lever_force_split))/xsg_sample_rate);
    
    % [img_filename img_path] = uigetfile('*.tif','Pick sample file','Multiselect','off');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    
    img_filename = dir_filenames(tiff_filename_indx);
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    
    if strcmp(img_filename{1}(8:11),'AP60')
        framerate = 28.7224;
        disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
    end
    
    for i = 1:length(acq_filename_all);
        xsg = load([acq_path filesep acq_filename_all{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([acq_path filesep acq_filename_all{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        
        % Cut out pieces of the lever trace that weren't imaged (dropped
        % frames) - Assume they were clipped at the end
        curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
        curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
        curr_lever_trace = xsg.data.acquirer.trace_2(1:curr_imaged_lever_samples);
        
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; curr_lever_trace];
        trial_stitch = [trial_stitch trial_list(end,2)];
    end
    
    % resample lever force to be 1:1 with numframes
    numframes = size(roi_trace_long,2);
    % to make resampling work: numframes must be even number
    % if it's not, cut out one frame at the end
    if mod(numframes,2) ~= 0
        roi_trace_long(:,end) = [];
        numframes = size(roi_trace_long,2);
    end
    % resample lever force to so # samples = # frames
    [n d] = rat(numframes/length(lever_force));
    lever_force_resample = resample(lever_force,n,d);
    % add or delete last sample if necessary
    if length(lever_force_resample) < numframes
        lever_force_resample(end+1:numframes) = 0;
    elseif length(lever_force_resample) > numframes
        lever_force_resample = lever_force_resample(1:numframes);
    end
    
    % Get trial offsets (seconds)
    xsg_bhv_offset = nan(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            continue;
        end
    end
    
    % Get dispatcher times in frames
    cue_frames = [];
    reward_frames = [];
    lick_frames = [];
    for curr_trial = 1:length(xsg_bhv_offset);
        if isnan(xsg_bhv_offset(curr_trial))
            continue
        end
        
        % relative cue time
        cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        % relative reward time
        if ~isempty(all_reward{curr_trial})
            reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        end
        % relative lick time
        lick_frames = [lick_frames;(all_licks{curr_trial} + xsg_bhv_offset(curr_trial))*framerate];
    end
    
    %%%%
    %%%% Second: deconv spikes
    %%%%
    % Do AK baseline on concat traces
    [roi_trace_df spike_estimates] = AK_baselineEstimation(roi_trace_long,framerate,1);
    
    % don't think it's necessary to ignore junctions
    
    spike_train = zeros(size(roi_trace_df));
    for curr_cell = 1:size(roi_trace_df,1)
        F = roi_trace_df(curr_cell,:);
        V = struct;
        V.smc_iter_max = 0;
        V.dt = 1/framerate;
        
        P = struct;
        % estimate number of spikes with Alex's method
        P.k = sum(spike_estimates(curr_cell,:));
        
        [M P V] = smc_oopsi(F,V,P);
        [n_best P_best V_est]=fast_oopsi(F,V);
        spike_train(curr_cell,:) = n_best;
        disp(num2str(curr_cell/size(roi_trace_df,1)));
    end
    
    % get rid of spikes 20 frames+ junction points (assume 8 files/loop)
    file_frames = cellfun(@length,roi_trace_long_split);
    file_frames_start = [1;cumsum(file_frames)+1];
    loop_frames_start = file_frames_start(1:8:end);
    for i = 1:length(loop_frames_start);
        if i == 1
            spike_train(:,loop_frames_start(i):loop_frames_start(i)+20) = 0;
        else
            spike_train(:,loop_frames_start(i)-20:loop_frames_start(i)+20) = 0;
        end
    end
    
    %%%%
    %%%% Third: weighted avg
    %%%%
    % get nonzero values of spike train
    time_back = 2000; %ms
    time_forward = 2000; %ms
    
    samples_back = round((time_back/1000)*(xsg_sample_rate));
    samples_forward = round((time_forward/1000)*(xsg_sample_rate));
    weighted_force = zeros(samples_back+samples_forward+1,size(roi_trace_df,1));
    
    
    lever_force_smooth = smooth(lever_force,500);
    lever_force_diff = [diff(lever_force_smooth);0];
    
    for curr_cell = 1:size(roi_trace_df,1)
        % % do oopsi to guess at spikeish regions
        % n_best = {};
        % P_best = {};
        % V = {};
        % V.dt = 1/30;
        % % for i = 1:size(roi_trace_df,1)
        % %     [n_best{i} P_best{i} V_est]=fast_oopsi(roi_trace_df(i,:),V);
        % %     disp(num2str(i/size(roi_trace_df,1)));
        % % end
        % % n_best = [n_best{:}]';
        %
        % [n_best P_best V_est]=fast_oopsi(roi_trace_df(curr_cell,:),V);
        
        n_best = spike_train(curr_cell,:);
        
        spike_indx = find(n_best > 0.5);
        %surround_force = nan(length(spike_indx),samples_back+samples_forward+1);
        
        for i = 1:length(spike_indx)
            curr_spike = spike_indx(i);
            curr_spike_sample = round((curr_spike/framerate)*xsg_sample_rate);
            % only use if within bounds
            if curr_spike_sample-samples_back > 0 && ...
                    curr_spike_sample + samples_forward < length(lever_force)
                curr_force = lever_force_diff(curr_spike_sample-samples_back...
                    :curr_spike_sample+samples_forward);
                %surround_force(i,:) = curr_force';
                if ~all(isnan(curr_force))
                    weighted_force(:,curr_cell) = weighted_force(:,curr_cell) +...
                        curr_force*(n_best(curr_spike)/sum(n_best(spike_indx)));
                end
            end
            %i/length(spike_indx)
        end
        % % get rid of nan rows
        % badrows = all(isnan(surround_force),2);
        % surround_force(badrows,:) = [];
        % spike_indx(badrows) = [];
        %
        % weighted_force(:,curr_cell) = surround_force'*n_best(spike_indx)';
        disp(['Finished cell ' num2str(curr_cell)])
    end
    
    %%%%
    %%%% Fourth: save
    %%%%
    save([curr_loop_folder filesep animal '_' listing_cell(trace_loop_file,:) '_weighted_avg']);
end

%% Get spike-triggered average over days
% UPDATE THIS! new code under bscope .m

clear all
animal = 'AP60';
listing_cell = [120330 120331 120401 120402 120403 120404 120409 ...
    120410 120411 120412 120413 120414 120415 120416 120417 120418];

listing_cell = num2str((listing_cell'));

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal
    curr_loop_folder = ['/usr/local/lab/People/Andy/Data/' animal filesep listing_cell(trace_loop_file,:)];
    data_filename = 'AP60_template_roi.mat';
    
    data_struct = load([curr_loop_folder filesep data_filename]);
    unstruct(data_struct);
    
    %%%%
    % Step 1: get behavior
    %%%%
    
    % Get trace and behavior times
    
    animal;
    day = listing_cell(trace_loop_file,:);
    
    data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
    
    % Get behavior
    % [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    
    warning off
    bhv = load([data_path filesep bhv_filename],'-MAT');
    warning on
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_licks = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get lever force
    
    % find and load the XSG files
    acq_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_acq = dir(acq_path);
    dir_acq_filenames = {dir_currfolder_acq.name};
    acq_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_acq_filenames);
    
    
    acq_filename_all = dir_acq_filenames(acq_filename_indx);
    % make sure they're in order
    acq_filename_all = sort(acq_filename_all);
    
    lever_force_split = cell(length(acq_filename_all),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    disp('ASSUMING: 8 FULL FILES FOR EACH LOOP')
    
    % Get framerate
    
    % Woah woah WOAH woah woah - this isn't always correct! This depends on
    % hitting the 'measure framerate' button at the beginning, and even then
    % it may not be right. Here's what we'll do instead: above, get the
    % framerate by dividing the number of frames by the length of the xsg in
    % seconds. Note that this is not robust to dropped frames.
    % Edit: at the moment, I know how many frames there should be, so I'll just
    % hardcode it.
    % 2nd Edit: This doesn't make sense - the XSG is slightly longer than the
    % acq, the ONLY information is 'framerate' in the header, which is really
    % bullshit that it's not. For now, assume it's right, but fix for AP60
    % which I know didn't have measured framerate for the first couple days
    % disp('ASSUMING 4000 FRAMES per XSG!');
    % framerate = 4000/(mean(cellfun(@length,lever_force_split))/xsg_sample_rate);
    
    % [img_filename img_path] = uigetfile('*.tif','Pick sample file','Multiselect','off');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    
    img_filename = dir_filenames(tiff_filename_indx);
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    
    if strcmp(img_filename{1}(8:11),'AP60')
        framerate = 28.7224;
        disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
    end
    
    for i = 1:length(acq_filename_all);
        xsg = load([acq_path filesep acq_filename_all{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([acq_path filesep acq_filename_all{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        
        % Cut out pieces of the lever trace that weren't imaged (dropped
        % frames) - Assume they were clipped at the end
        curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
        curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
        curr_lever_trace = xsg.data.acquirer.trace_2(1:curr_imaged_lever_samples);
        
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; curr_lever_trace];
        trial_stitch = [trial_stitch trial_list(end,2)];
    end
    
    % resample lever force to be 1:1 with numframes
    numframes = size(roi_trace_long,2);
    % to make resampling work: numframes must be even number
    % if it's not, cut out one frame at the end
    if mod(numframes,2) ~= 0
        roi_trace_long(:,end) = [];
        numframes = size(roi_trace_long,2);
    end
    % resample lever force to so # samples = # frames
    [n d] = rat(numframes/length(lever_force));
    lever_force_resample = resample(lever_force,n,d);
    % add or delete last sample if necessary
    if length(lever_force_resample) < numframes
        lever_force_resample(end+1:numframes) = 0;
    elseif length(lever_force_resample) > numframes
        lever_force_resample = lever_force_resample(1:numframes);
    end
    
    % Get trial offsets (seconds)
    xsg_bhv_offset = nan(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            continue;
        end
    end
    
    % Get dispatcher times in frames
    cue_frames = [];
    reward_frames = [];
    lick_frames = [];
    for curr_trial = 1:length(xsg_bhv_offset);
        if isnan(xsg_bhv_offset(curr_trial))
            continue
        end
        
        % relative cue time
        cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        % relative reward time
        if ~isempty(all_reward{curr_trial})
            reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        end
        % relative lick time
        lick_frames = [lick_frames;(all_licks{curr_trial} + xsg_bhv_offset(curr_trial))*framerate];
    end
    
    %%%%
    % Step 2: baseline correct traces
    %%%%
    
    roi_trace_df = AK_baselineEstimation(roi_trace_long,framerate);
    
    %%%%
    % Step 3: get large calcium transients (the dumb way)
    %%%%
    
    % Get calcium transients (the dumb way)
    
    %peakstd = [1.5*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
    peakstd = [4 4 3 3 2 2 1 1]; % definitions for max trigger, in stds
    peak_back = 20;
    local_max = 30;
    drop_percentage = .8;
    spike_frames = {};
    spike_amplitudes = {};
    
    num_rois = size(roi_trace_long,1);
    
    % this is to catch a dumb bug: if ROI is offscreen, just make noise
    for i = 1:num_rois
        if all(roi_trace_long(i,:) == roi_trace_long(i,1))
            roi_trace_long(i,:) = rand(size(roi_trace_long(i,:)));
        end
    end
    
    % Estimate noise through smoothing
    size_smooth = 30;
    smooth_type = 'loess';
    % for now - just smooth whole trace
    
    % size_file = 4000;
    % roi_trace_df_smooth = zeros(size(roi_trace_df));
    % for i = 1:size(roi_trace_df,1);
    %     for j = 1:size(roi_trace_df,2)/size_file;
    %         roi_trace_df_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_df(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    %     end
    % end
    for i = 1:size(roi_trace_df,1)
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),size_smooth,smooth_type);
    end
    
    smooth_std = zeros(size(roi_trace_df,1),1);
    smooth_std = sqrt(sum((abs(roi_trace_df - roi_trace_df_smooth).^2)/size(roi_trace_df_smooth,2),2));
    
    for i = 1:num_rois
        % Estimate noise through smoothing
        curr_trace_std = smooth_std(i);
        minpeak = curr_trace_std*peakstd;
        [spike_frames{i} spike_amplitudes{i}] = ap_schmittDetect(roi_trace_df_smooth(i,:),minpeak,peak_back,local_max,drop_percentage);
        disp([num2str(i) '/' num2str(num_rois)])
    end
    
    % make spike matrix
    spike_matrix = zeros(size(roi_trace_df));
    for i = 1:num_rois
        spike_matrix(i,spike_frames{i}) = 1;
    end
    
    
    
    %%%%
    % Step 4: make amplitude-weighted spike-triggered average
    %%%%
    
    % Make spike triggered avg weighted by spike amplitude
    
    time_back = 2000; %ms
    time_forward = 2000; %ms
    
    samples_back = round((time_back/1000)*(xsg_sample_rate));
    samples_forward = round((time_forward/1000)*(xsg_sample_rate));
    
    frames_back = round((time_back/1000)*framerate);
    frames_forward = round((time_forward/1000)*framerate);
    
    lever_spike_all = cell(size(roi_trace_df,1),1);
    trace_spike_all = cell(size(roi_trace_df,1),1);
    
    lever_force_smooth = smooth(lever_force,500);
    lever_force_diff = [0;diff(lever_force_smooth)];
    
    for curr_cell = 1:size(roi_trace_df,1);
        lever_spike = nan(length(spike_frames{curr_cell}),samples_back+samples_forward+1);
        lever_diff_spike = nan(length(spike_frames{curr_cell}),samples_back+samples_forward+1);
        trace_spike = nan(length(spike_frames{curr_cell}),frames_back+frames_forward+1);
        for i = 1:length(spike_frames{curr_cell})
            
            curr_spike = spike_frames{curr_cell}(i);
            curr_spike_sample = round((curr_spike/framerate)*xsg_sample_rate);
            % only use if within bounds
            if curr_spike_sample-samples_back > 0 && ...
                    curr_spike_sample + samples_forward < length(lever_force)
                curr_force = lever_force_diff(curr_spike_sample-samples_back...
                    :curr_spike_sample+samples_forward);
                %surround_force(i,:) = (curr_force*(n_best(curr_spike)/sum(n_best(spike_indx))))';
                if ~all(isnan(curr_force))
                    lever_diff_spike(i,:) = curr_force/spike_amplitudes{curr_cell}(i);
                    trace_spike(i,:) = roi_trace_df(curr_spike-frames_back: ...
                        curr_spike + frames_forward);
                end
            end
        end
        lever_diff_spike_all{curr_cell} = lever_diff_spike;
        trace_spike_all{curr_cell} = trace_spike;
        curr_cell
    end
    lever_diff_spike_all_mean = cellfun(@(x) nanmean(x,1)',lever_diff_spike_all,'UniformOutput',false);
    lever_diff_spike_all_mean = [lever_diff_spike_all_mean{:}]';
    lever_diff_spike_all_mean_norm = bsxfun(@times,lever_diff_spike_all_mean,1./sqrt(sum(lever_diff_spike_all_mean.^2,2)));
    
    trace_spike_all_mean = cellfun(@(x) nanmean(x,1)',trace_spike_all,'UniformOutput',false);
    trace_spike_all_mean = [trace_spike_all_mean{:}]';
    
    
    %%%%
    % Step 5: align spike-triggered averages
    %%%%
    
    % Align spike-triggered averages
    spike_error = 2000;
    lever_diff_spike_all_aligned = lever_diff_spike_all;
    
    shift = spike_error;
    center = ceil(size(lever_diff_spike_all{1},2)/2);
    
    spike_shift = spike_frames;
    
    for i = 1:length(lever_diff_spike_all);
        if size(lever_diff_spike_all{i},1) > 10;
            temp_smooth_force = smooth(lever_diff_spike_all_mean(i,:),2000);
            for j = 1:size(lever_diff_spike_all{i},1)
                temp_xcov = ...
                    xcov(lever_diff_spike_all{i}(j,:), ...
                    temp_smooth_force,spike_error);
                
                %                 xcov(lever_diff_spike_all{i} ...
                %                 (j,center-shift:center+shift), ...
                %                 mean(lever_diff_spike_all{i}(:,center-shift:center+shift)));
                xsg_shift = [];
                xsg_shift = spike_error+1-find(temp_xcov == max(temp_xcov));
                lever_diff_spike_all_aligned{i}(j,:) = ...
                    circshift(lever_diff_spike_all_aligned{i}(j,:),[0 xsg_shift]);
                % shift spikes accordingly
                if ~isempty(xsg_shift)
                    spike_shift{i}(j) = spike_shift{i}(j) + round(framerate*(xsg_shift/xsg_sample_rate));
                end
            end
        else lever_diff_spike_all_aligned{i} = nan(size(lever_diff_spike_all_aligned{i}));
        end
        i
    end
    lever_diff_spike_all_aligned_mean = [];
    lever_diff_spike_all_aligned_mean = cellfun(@(x) nanmean(x,1)',lever_diff_spike_all_aligned,'UniformOutput',false);
    lever_diff_spike_all_aligned_mean = [lever_diff_spike_all_aligned_mean{:}]';
    % L2 normalize mean
    lever_diff_spike_all_aligned_mean_norm = [];
    lever_diff_spike_all_aligned_mean_norm = bsxfun(@times,lever_diff_spike_all_aligned_mean,1./sqrt(sum(lever_diff_spike_all_aligned_mean.^2,2)));
    
    lever_diff_spike_all_aligned_median = [];
    lever_diff_spike_all_aligned_median = cellfun(@(x) nanmedian(x,1)',lever_diff_spike_all_aligned,'UniformOutput',false);
    lever_diff_spike_all_aligned_median = [lever_diff_spike_all_aligned_median{:}]';
    % L2 normalize mean
    lever_diff_spike_all_aligned_median_norm = [];
    lever_diff_spike_all_aligned_median_norm = bsxfun(@times,lever_diff_spike_all_aligned_median,1./sqrt(sum(lever_diff_spike_all_aligned_median.^2,2)));
    
    lever_diff_spike_all_aligned_std = [];
    lever_diff_spike_all_aligned_std = cellfun(@(x) nanstd(x,[],1)',lever_diff_spike_all_aligned,'UniformOutput',false);
    lever_diff_spike_all_aligned_std = [lever_diff_spike_all_aligned_std{:}]';
    
    save([curr_loop_folder filesep 'STA_analysis']);
    
end

%% Given STA analysis, go through and find significantly modulated

figure; hold on;
plot(lever_diff_spike_all_aligned_mean(16,:))
curr_std = lever_diff_spike_all_aligned_std(curr_cell,:);
plot(lever_diff_spike_all_aligned_mean(16,:)+curr_std,'r')
plot(lever_diff_spike_all_aligned_mean(16,:)-curr_std,'r')






%% ~~~~~~~~!! AFTER THIS = WORKING SCRIPTS !!~~~~~~~~~~~

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% From here on, scripts are new
%%%%%%%%%%%%%%

%% jPCA batch
clc
clear all

animal = 'AP61';
days = {120417 120419 120421 120422 120423 120424 120425 120426 ...
    120427 120428 120429 120430 120501};

for curr_day = 1:length(days)
    clearvars -except animal days curr_day
    
    day = num2str(days{curr_day});
    
    data_filename = [day '_' animal '_summed_50_analysis.mat'];
    
    data_struct = load([data_filename]);
    unstruct(data_struct);
    
    % using the top 6 PCs
    
    reward_frames = bhv.reward_frames;
    % if you want just random frames, for testing
    %reward_frames = randi([round(min(reward_frames)), ...
    %    round(max(reward_frames))], length(reward_frames));
    
    frames_span = 20;
    frames_surround = frames_span*2+1;
    
    % pop out data around the time of rewarded lever presses
    total_reward_frames = frames_surround*length(reward_frames);
    num_rois = size(roi_trace_df,1);
    lever_reward_trace = zeros(num_rois,length(reward_frames)*frames_surround);
    
    smooth_size = 30;
    roi_trace_df_smooth = nan(size(roi_trace_df));
    for i = 1:size(roi_trace_df,1);
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
    end
    
    % Create a thresholded trace
    roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
    for curr_cell = 1:size(roi_trace_df_smooth,1)
        
        % call the minimum of the raw the low cutoff
        % estimate the noise from smoothed - raw
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 3*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    
    for i = 1:length(reward_frames);
        curr_frame = round(reward_frames(i));
        curr_reward_frames = curr_frame-frames_span:curr_frame+frames_span;
        curr_reward_trace = roi_trace_df_smooth_cutoff(:,curr_reward_frames);
        lever_reward_trace(:,((i-1)*frames_surround+1): ...
            ((i-1)*frames_surround)+frames_surround) = curr_reward_trace;
    end
    
    % % subtract cell all-trial mean, they use this to separate trials
    % for i = 1:size(roi_trace_df,1)
    %     trial_reshape = reshape(lever_reward_trace(i,:), ...
    %         frames_surround,length(reward_frames))';
    %     trial_mean = mean(trial_reshape);
    %
    %     trial_meansub = (trial_reshape - ...
    %         repmat(trial_mean,length(reward_frames),1))';
    %     lever_reward_trace(i,:) = trial_meansub(:);
    % end
    
    % perform PCA on concatenated lever rewarded traces
    [coeff score latent] = princomp(lever_reward_trace');
    pcs = lever_reward_trace'*coeff(:,1:6);
    
    % dynamics: find transformation matrix from PCs to derivative of PCs
    pcs_diff = diff(pcs);
    pcs_compare = pcs(1:end-1,:);
    
    X = pcs_compare;
    X_dot = pcs_diff;
    M = X\X_dot;
    
    X_tilde = blkdiag(X,X,X,X,X,X);
    
    % find the vector k which ultimately gives M_skew
    k_start = ones(15,1)*0.001;
    options = optimset('MaxFunEvals',300000);
    k = fminsearch(@(x) Churchland2012_jPCA_fminsearch(x,X_tilde,X_dot),k_start,options);
    
    % construct final M_skew
    idx = ones(6);
    tril_idx = logical(tril(idx,-1));
    triu_idx = logical(triu(idx,1));
    M_skew = zeros(6);
    M_skew(tril_idx) = k;
    M_skew_t = M_skew';
    M_skew_t(tril_idx) = -k;
    M_skew = M_skew_t';
    
    % find the largest eigenvalues and pull out the top 2 (complex conj pair)
    clear j
    [V D] = eig(M_skew);
    D = diag(D);
    jPCA_idx = find(abs(D) == max(abs(D)));
    jPCA_1 = V(:,jPCA_idx(1)) + V(:,jPCA_idx(2));
    jPCA_2 = j*(V(:,jPCA_idx(1)) - V(:,jPCA_idx(2)));
    
    % get the final jPCs, project data onto jPC space
    pop_jPCs = X*[jPCA_1 jPCA_2];
    jpc_fig = figure;hold on;
    title(['jPCs:' day])
    for i = 1:length(reward_frames)
        curr_frames = (i-1)*frames_surround+1:(i)*frames_surround-1;
        z = zeros(size(curr_frames));
        x = pop_jPCs(curr_frames,1)';
        y = pop_jPCs(curr_frames,2)';
        col = curr_frames - curr_frames(1);
        surface([x;x],[y;y],[z;z],[col;col],'facecol','no', ...
            'edgecol','interp','linew',2);
        arrow([x(end-1) y(end-1)],[x(end) y(end)],'length',10)
    end
    
    print(jpc_fig,['jPC_analysis/' day '_jPCs'],'-dpng');
    close all
    % % plot jPC coefficients
    % figure; title('jPC coefficients');
    % plot(lever_reward_trace(:,1:end-1)'\pop_jPCs);
    disp(['Finished day ' day]);
end

%% Get all bhv-aligned activity

animal = 'AP74';

% find days
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

for curr_day = 1:length(days)
    
    clearvars -except animal days curr_day activity_aligned ...
        roi_image activity rewarded_move_start
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_name = [day '_' animal '_bgROI_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    if curr_day == 1;
        activity_aligned = cell(size(im.roi_trace_df,1),1);
        roi_image.max = cell(size(im.roi_trace_df,1),1);
        roi_image.mean = cell(size(im.roi_trace_df,1),1);
        activity.peaks = cell(length(days),1);
        activity.tauish = cell(length(days),1);
    end
    
    % get lever movement
    % get xsg files
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    % get binary traces for movement, quiescence
    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    loop_size_split = cellfun(@length,im.roi_trace_split);
    loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
    loop_size_sum = cumsum(loop_size_split);
    loop_frames = [0;loop_size_sum(8:8:end-7)];
    
    lever_active_full = cell(length(xsg_filenames),1);
    lever_velocity_full = cell(length(xsg_filenames),1);
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        
        [lever_active lever_force_smooth lever_force_velocity] = ...
            AP_parseLeverMovement(xsg_data);
        
        %save the active lever portions in frames
        curr_lever_active_frames = round((find(lever_active)/1000) ...
            *bhv.framerate);
        lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
            curr_lever_active_frames <= loop_size(curr_xsg) & ...
            curr_lever_active_frames > 0)) + loop_frames(curr_xsg);
        
        [n d] = rat(loop_size(curr_xsg)/length(lever_force_velocity));
        lever_force_velocity_resample = resample(lever_force_velocity,n,d);
        lever_velocity_full{curr_xsg} = ...
            lever_force_velocity_resample(1:loop_size(curr_xsg));
    end
    
    % get movement epochs
    lever_movement = zeros(size(im.roi_trace_df,2),1);
    lever_movement(vertcat(lever_active_full{:})) = 1;
    lever_stopstart = diff([0;lever_movement;0]);
    
    lever_start = num2cell(find(lever_stopstart == 1));
    lever_stop = num2cell(find(lever_stopstart == -1));
    
    lever_movement_epochs = cellfun(@(x,y) [x:y],lever_start,lever_stop,'uni',false);
    lever_quiescent_epochs = cellfun(@(x,y) [x:y],lever_stop(1:end-1),lever_start(2:end),'uni',false);
    
    for curr_cell = 1:size(im.roi_trace_df,1)
        
        % grab reward-aligned traces
        frames_surround = 200;
        roi_reward = nan(length(bhv.reward_frames),frames_surround*2+1);
        for i = 1:length(bhv.reward_frames);
            curr_reward = round(bhv.reward_frames(i));
            roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
                continue
            end
            roi_reward(i,:) = im.roi_trace_df(curr_cell,roi_reward_frames);
        end
        
        % grab cue-aligned traces
        % if rewared move onsets provided, use them (from
        % AP_leverPaperPrep)
        if exist('rewarded_move_start','var')
            use_frames_cue = rewarded_move_start{curr_day};
        else
           use_frames_cue = bhv.cue_frames; 
        end
        frames_surround = 200;
        cue_reward = nan(length(use_frames_cue),frames_surround*2+1);
        for i = 1:length(use_frames_cue);
            curr_reward = round(use_frames_cue(i));
            roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
                continue
            end
            cue_reward(i,:) = im.roi_trace_df(curr_cell,roi_reward_frames);
        end
        
        %         % plot around movement, movement start-aligned
        %         frames_surround = 200;
        %         move_aligned = nan(length(move_start_frames),frames_surround*2+1);
        %         for i = 1:length(move_start_frames);
        %             curr_reward = round(move_start_frames(i));
        %             roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        %             if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
        %                 continue
        %             end
        %             move_aligned(i,:) = im.roi_trace_df(curr_cell,roi_reward_frames);
        %         end
        
        % plot during movement epochs
        movement_activity = cell(size(lever_movement_epochs));
        for i = 1:length(lever_movement_epochs)
            movement_activity{i} = im.roi_trace_df(curr_cell,lever_movement_epochs{i});          
        end
        
        % plot during quiescent epochs
        quiescent_activity = cell(size(lever_quiescent_epochs));
        for i = 1:length(lever_quiescent_epochs)
            quiescent_activity{i} = im.roi_trace_df(curr_cell,lever_quiescent_epochs{i});          
        end
        
        activity_aligned{curr_cell}.reward{curr_day} = roi_reward;
        activity_aligned{curr_cell}.cue{curr_day} = cue_reward;
        %activity_aligned{curr_cell}.move_start{curr_day} = move_aligned;
        activity_aligned{curr_cell}.movement_activity{curr_day} = movement_activity;
        activity_aligned{curr_cell}.quiescent_activity{curr_day} = quiescent_activity;
    end
    
%     % get number of peaks and tau-ish estimate
%     [peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);
%     
%     % threshold based on noise and get frames above
%     roi_trace_df_thresh = zeros(size(im.roi_trace_df));
%     tauish = zeros(size(im.roi_trace_df,1),1);
%     for curr_cell = 1:size(im.roi_trace_df,1);
%         if isempty(peak_frames{curr_cell});
%             continue
%         end
%         % threshold trace at 80% smallest amplitude peak
%         trace_thresh = 0.8*min(peak_amplitudes{curr_cell});
%         over_thresh = im.roi_trace_df(curr_cell,:) > trace_thresh;
%         tauish(curr_cell) = ...
%             sum(over_thresh) / ...
%             sum(peak_amplitudes{curr_cell}-trace_thresh);
%     end
%     
%     activity.peaks{curr_day} = cellfun(@length,peak_frames);
%     activity.tauish{curr_day} = tauish;
    
    disp(['Grabbed activity of day ' day ', getting roi images (%)...'])
    
    % pull out an max image of the ROI too
    
    % get the summed movie
    im_path = ['/usr/local/lab/People/Andy/Data/' animal filesep day];
    im_file = [day '_' animal '_summed_50.tif'];
    
    curr_imfinfo = imfinfo([im_path filesep im_file]);
    curr_im = zeros(512,512,length(curr_imfinfo),'uint16');
    for curr_frame = 1:length(curr_imfinfo);
        curr_im(:,:,curr_frame) = ...
            imread([im_path filesep im_file],'tiff',curr_frame);
    end
    
    % get the ROIs 
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
        animal '_roi_template'];
    roi_file = [day '_' animal '_summed_50_template.roi'];
    rois = load([roi_path filesep roi_file],'-MAT');
    
    fprintf('%2d',0);
    % loop through ROIs, get max projection
    for curr_cell = 1:length(rois.polygon.ROI)
        % get center of polygon
        [parea cx cy] = polycenter(rois.polygon.ROI{curr_cell}(:,1), ...
            rois.polygon.ROI{curr_cell}(:,2));
        cx = round(cx);
        cy = round(cy);
        %max_size = ceil(max(max(rois.polygon.ROI{curr_cell}) - ...
        %    [cx cy]))+10;
        max_size = 20;
        
        cell_pixels = [cx-max_size:cx+max_size;cy-max_size:cy+max_size]';
        x_oob_min = find(cell_pixels(:,1) < 1,1,'last');
        x_oob_max = find(cell_pixels(:,1) > size(curr_im,1),1,'first');
        y_oob_min = find(cell_pixels(:,1) < 1,1,'last');
        y_oob_max = find(cell_pixels(:,1) > size(curr_im,1),1,'first');
        
        cell_pixels(cell_pixels < 1) = 1;
        cell_pixels(cell_pixels > size(curr_im,1)) = size(curr_im,1);
        
        roi_image.max{curr_cell}{curr_day} = max(curr_im( ...
            cell_pixels(:,2),cell_pixels(:,1),:),[],3);
        roi_image.mean{curr_cell}{curr_day} = mean(curr_im( ...
            cell_pixels(:,2),cell_pixels(:,1),:),3);
        
        % if pixels out of bounds, make them black
        roi_image.max{curr_cell}{curr_day}(:, ...
            [1:x_oob_min x_oob_max:end]) = 0;
        roi_image.max{curr_cell}{curr_day}(:, ...
            [1:y_oob_min y_oob_max:end]) = 0;
        
        roi_image.mean{curr_cell}{curr_day}(:, ...
            [1:x_oob_min x_oob_max:end]) = 0;
        roi_image.mean{curr_cell}{curr_day}(:, ...
            [1:y_oob_min y_oob_max:end]) = 0;
        
        
        fprintf('%c%c%2d',8,8,round(100*(curr_cell/length(rois.polygon.ROI))));
    end
    
    disp(['Finished day ' day])
end

disp('Finished all')

%% Plot aligned activity from above
curr_cell = 64;

% plot all roi images
figure;
set(gcf,'Name',[num2str(curr_cell) ' ROI Image']);
for i = 1:length(activity_aligned{curr_cell}.cue)
    subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
    imagesc(roi_image.mean{curr_cell}{i});
    colormap(gray);
end
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1)+fig_size(3) fig_size(2) ...
    430 screen_size(4)])

% plot all cue-aligned
figure;
set(gcf,'Name',[num2str(curr_cell) ' Cue-aligned']);
for i = 1:length(activity_aligned{curr_cell}.cue)
    subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
    imagesc(activity_aligned{curr_cell}.cue{i});
    colormap(gray);
    caxis([0 2]);
    %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',1,'color','w',...
    %    'LineStyle','--');
end
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1)-fig_size(3) fig_size(2) ...
    fig_size(3) screen_size(4)])

% plot all reward-aligned
figure;
set(gcf,'Name',[num2str(curr_cell) ' Reward-aligned']);
for i = 1:length(activity_aligned{curr_cell}.reward)
    subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
    imagesc(activity_aligned{curr_cell}.reward{i});
    colormap(gray);
    caxis([0 2]);
    %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',1,'color','w',...
    %    'LineStyle','--');
end
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1:3) screen_size(4)])

% % plot activity: number of peaks, tauish
% figure; hold on
% set(gcf,'Name',[num2str(curr_cell) ' Activity']);
% num_days = 1:length(activity.peaks);
% curr_peaks = cellfun(@(x) x(curr_cell),activity.peaks);
% curr_tauish = cellfun(@(x) x(curr_cell),activity.tauish);
% [haxes] = plotyy(num_days,curr_peaks,num_days,curr_tauish);
% axes(haxes(1)); ylabel('Number of peaks');
% axes(haxes(2)); ylabel('Tau-ish');
% xlabel('Session');
% 
% % plot all movement-aligned
% figure;
% set(gcf,'Name',[num2str(curr_cell) ' Movement-aligned']);
% for i = 1:length(activity_aligned{curr_cell}.movement_activity)
%     subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
%     max_length_activity = max(cellfun(@length,...
%         activity_aligned{curr_cell}.movement_activity{i}));
%     curr_activity_pad = cellfun(@(x) padarray(x,[0 max_length_activity - ...
%         length(x)],'post'),activity_aligned{curr_cell}.movement_activity{i}, ...
%         'uni',false);
%     imagesc(cell2mat(curr_activity_pad));
%     colormap(gray);
%     caxis([0 1]);
%     %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',1,'color','w',...
%     %    'LineStyle','--');
% end
% screen_size = get(0,'ScreenSize');
% fig_size = get(gcf,'Position');
% set(gcf,'Position',[fig_size(1)-fig_size(3) fig_size(2) ...
%     fig_size(3) screen_size(4)])
% 
% % plot all quiescent-aligned
% figure;
% set(gcf,'Name',[num2str(curr_cell) ' Quiescent-aligned']);
% for i = 1:length(activity_aligned{curr_cell}.quiescent_activity)
%     subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
%     max_length_activity = max(cellfun(@length,...
%         activity_aligned{curr_cell}.quiescent_activity{i}));
%     curr_activity_pad = cellfun(@(x) padarray(x,[0 max_length_activity - ...
%         length(x)],'post'),activity_aligned{curr_cell}.quiescent_activity{i}, ...
%         'uni',false);
%     imagesc(cell2mat(curr_activity_pad));
%     colormap(gray);
%     caxis([0 1]);
%     %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',1,'color','w',...
%     %    'LineStyle','--');
% end
% screen_size = get(0,'ScreenSize');
% fig_size = get(gcf,'Position');
% set(gcf,'Position',[fig_size(1:3) screen_size(4)])

% % plot all move-start-aligned
% figure;
% set(gcf,'Name',[num2str(curr_cell) ' Move start-aligned']);
% for i = 1:length(activity_aligned{curr_cell}.move_start)
%     subplot(7,2,i)
%     imagesc(activity_aligned{curr_cell}.move_start{i});
%     caxis([0 3]);
%     %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',1,'color','w',...
%     %    'LineStyle','--');
% end
% screen_size = get(0,'ScreenSize');
% fig_size = get(gcf,'Position');
% set(gcf,'Position',[fig_size(1:3) screen_size(4)])
% 
figure;
set(gcf,'Name',[num2str(curr_cell) ' Cue-aligned']);
for i = 1:length(activity_aligned{curr_cell}.cue)
    subplot(1,length(activity_aligned{curr_cell}.cue),i)

    cue_mean = nanmean(activity_aligned{curr_cell}.cue{i});
    cue_sem = nanstd(activity_aligned{curr_cell}.cue{i})./ ...
        sqrt(size(activity_aligned{curr_cell}.cue{i},1));
    plot(cue_mean,'k');
    jbfill([1:length(cue_mean)],cue_mean+cue_sem, ...
    cue_mean-cue_sem,[0.5 0.5 0.5]);
    ylim([-0.2 2])   
end
% 
% figure;
% set(gcf,'Name',[num2str(curr_cell) ' Reward-aligned']);
% for i = 1:length(activity_aligned{curr_cell}.reward)
%     subplot(1,length(activity_aligned{curr_cell}.reward),i)
% 
%     reward_mean = nanmean(activity_aligned{curr_cell}.reward{i});
%     reward_sem = nanstd(activity_aligned{curr_cell}.reward{i})./ ...
%         sqrt(size(activity_aligned{curr_cell}.reward{i},1));
%     plot(reward_mean,'k');
%     jbfill([1:length(reward_mean)],reward_mean+reward_sem, ...
%     reward_mean-reward_sem,[0.5 0.5 0.5]);
%    ylim([-0.2 2]) 
% end

% figure;
% set(gcf,'Name',[num2str(curr_cell) ' Move start-aligned']);
% for i = 1:length(activity_aligned{curr_cell}.move_start)
%     subplot(1,length(activity_aligned{curr_cell}.move_start),i)
% 
%     move_start_mean = nanmean(activity_aligned{curr_cell}.move_start{i});
%     move_start_sem = nanstd(activity_aligned{curr_cell}.move_start{i})./ ...
%         sqrt(size(activity_aligned{curr_cell}.move_start{i},1));
%     plot(move_start_mean,'k');
%     jbfill([1:length(move_start_mean)],move_start_mean+move_start_sem, ...
%     move_start_mean-move_start_sem,[0.5 0.5 0.5]);
%    ylim([-0.2 2]) 
% end

%% Plot above (all cells of a class)

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
curr_animal = 5;
animal = animals{curr_animal};

% load ROI labels and identify cell type
analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
roilabel_name = [animal '_roilabels.roilabel'];
load([analysis_path filesep roilabel_name],'-MAT')
cells = 1:length(roi_labels);
gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
pyr_cells = ~gad_cells;
filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';

pyr_unfilled_idx = find(pyr_cells & ~filled_cells);
gad_unfilled_idx = find(gad_cells & ~filled_cells);

img_max_fig = figure;
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1)-fig_size(3) fig_size(2) 430 screen_size(4)])

img_mean_fig = figure;
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1:2) 430 screen_size(4)])

activity_fig = figure;
screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[fig_size(1)+fig_size(3) fig_size(2) ...
    430 screen_size(4)])

%for curr_cell = gad_unfilled_idx
for curr_cell = [8 13 14 75 147 148]

    % plot all roi max images
    figure(img_max_fig);
    set(gcf,'Name',[num2str(curr_cell) ' ROI Image']);
    for i = 1:length(activity_aligned{curr_cell}.cue)
        subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
        imagesc(roi_image.max{curr_cell}{i});
        colormap(gray);
    end
    
    % plot all roi mean images
    figure(img_mean_fig);
    set(gcf,'Name',[num2str(curr_cell) ' ROI Image']);
    for i = 1:length(activity_aligned{curr_cell}.cue)
        subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
        imagesc(roi_image.mean{curr_cell}{i});
        colormap(gray);
    end
    saveas(img_mean_fig,[num2str(curr_animal) '_' num2str(curr_cell) '_mean']);
    
    figure(activity_fig);
    set(gcf,'Name',[num2str(curr_cell) ' Cue-aligned']);
    for i = 1:length(activity_aligned{curr_cell}.cue)
        subplot(ceil(length(activity_aligned{curr_cell}.cue)/2),2,i)
        
        cue_mean = nanmean(activity_aligned{curr_cell}.cue{i});
        cue_sem = nanstd(activity_aligned{curr_cell}.cue{i})./ ...
            sqrt(size(activity_aligned{curr_cell}.cue{i},1));
        plot(cue_mean,'k');
        jbfill([1:length(cue_mean)],cue_mean+cue_sem, ...
            cue_mean-cue_sem,[0.5 0.5 0.5]);
        ylim([-0.2 0.5])
    end
    saveas(activity_fig,[num2str(curr_animal) '_' num2str(curr_cell) '_act']);
    
    %waitforbuttonpress
end

%% Batch cell-cell correlations (based on thresholded traces)

animal = 'AP61';

days = {120417 120418 120419 120421 120422 120423 120424 120425 ...
    120426 120427 120428 120429 120430 120501};

cell_corr = cell(length(days),1);
for curr_day = 1:length(days)
    clearvars -except animal days curr_day cell_corr
    
    day = num2str(days{curr_day});
    
    load([day '_' animal '_summed_50_analysis.mat']);
    
    roi_trace_df = AP_baselineEstimation(roi_trace_long,28);
    
    size_smooth = 30;
    smooth_type = 'loess';
    for i = 1:size(roi_trace_df,1)
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),size_smooth,smooth_type);
    end
    
    % Create a thresholded trace
    roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
    for curr_cell = 1:size(roi_trace_df_smooth,1)
        % estimate the noise from smoothed - raw
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 2*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    
    cell_corr{curr_day} = corrcoef(roi_trace_df_smooth_cutoff');
    
    disp(['Finished correlations: day ' day]);
end

% this mostly just turns into a big pile of mud when you plot it

% % plot the correlations
% for curr_day = 1:length(days)
%     figure; hold on
%     title(num2str(days{curr_day}))
%     spacing = (2*pi)/length(cell_corr{end});
%     cell_dots = [0:spacing:2*pi]';
%     [cell_cart_x cell_cart_y] = pol2cart(cell_dots,ones(length(cell_dots),1));
%     % plot dots for cells
%     plot(cell_cart_x,cell_cart_y,'.r','MarkerSize',10);
%     % plot all correlations as lines
%     for i = 1:length(cell_corr{curr_day})
%         for j = 1:length(cell_corr{curr_day});
%             if i < j
%                 continue
%             end
%             curr_corr = cell_corr{curr_day}(i,j);
%             if curr_corr > 0
%                 r = patch([cell_cart_x(i) cell_cart_x(j)], ...
%                     [cell_cart_y(i) cell_cart_y(j)],'black');
%                 set(r,'facealpha',curr_corr);
%                 set(r,'edgealpha',curr_corr);
%             elseif curr_corr < 0
%                 curr_corr = -curr_corr;
%                 r = patch([cell_cart_x(i) cell_cart_x(j)], ...
%                     [cell_cart_y(i) cell_cart_y(j)],'red');
%                 set(r,'facealpha',curr_corr);
%                 set(r,'edgealpha',curr_corr);
%             end
%         end
%     end
% end



%% Batch cell-movement correlations (based on thresholded traces)

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 120613 ...
    120614 120615 120616 120617 120618 120619};

% AP61 20418 is fucked up (2bhv files)

movecorr_all = cell(length(days),1);
movecorr_value_all = cell(length(days),1);
for curr_day = 1:length(days)
    
    %%% INITIALIZE THE DATA
    
    clearvars -except animal days curr_day movecorr_all movecorr_value_all
    
    day = num2str(days{curr_day});
    
    load([day '_' animal '_summed_50_analysis.mat']);
    
    size_smooth = 30;
    smooth_type = 'loess';
    for i = 1:size(roi_trace_df,1)
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),size_smooth,smooth_type);
    end
    %%%
    
    disp('Finding movement-correlated cells')
    %%% PROCESS THE DATA
    
    % take velocity of lever to correlate
    lever_movement = smooth(abs(diff(bhv.lever_force_resample)),2);
    % make the first value 0 to keep same size
    lever_movement = [0;lever_movement];
    
    movecorr = zeros(size(roi_trace_df,1),1);
    for curr_cell = 1:size(roi_trace_df,1)
        
        % call the minimum of the raw the low cutoff
        % estimate the noise from smoothed - raw
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 4*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        curr_trace_cutoff(over_cutoff_indx) = roi_trace_df_smooth(curr_cell,over_cutoff_indx);
        
        % if less than 10 seconds is above cutoff, then continue
        if sum(over_cutoff_indx) < bhv.framerate*10
            continue
        end
        
        active_lever = mean(lever_movement(over_cutoff_indx));
        inactive_lever = mean(lever_movement(~over_cutoff_indx));
        
        tic
        active_lever_bootstrp = bootstrp(1000,@mean,lever_movement(over_cutoff_indx));
        active_lever_bootstrp_ci = prctile(active_lever_bootstrp,[0.5 99.5]);
        active_lever_bootstrp_ci_range = abs(active_lever - active_lever_bootstrp_ci);
        toc
        inactive_lever_bootstrp = bootstrp(1000,@mean,lever_movement(~over_cutoff_indx));
        inactive_lever_bootstrp_ci = prctile(inactive_lever_bootstrp,[0.5 99.5]);
        inactive_lever_bootstrp_ci_range = abs(inactive_lever - inactive_lever_bootstrp_ci);
        toc
        
        % Shuffle lever movement and check for significance, but only
        % bother shuffling if they're different off the bat
        
        % criterion for movement correlated
        if active_lever_bootstrp_ci(1) > inactive_lever_bootstrp_ci(2) || ...
                active_lever_bootstrp_ci(2) < inactive_lever_bootstrp_ci(1)
            
            % shuffle the lever_shuffle movement
            shuffle_indx = randperm(length(lever_movement));
            lever_shuffle_movement = lever_movement(shuffle_indx);
            
            active_lever_shuffle = mean(lever_shuffle_movement(over_cutoff_indx));
            inactive_lever_shuffle = mean(lever_shuffle_movement(~over_cutoff_indx));
            
            tic
            active_lever_shuffle_bootstrp = bootstrp(1000,@mean,lever_shuffle_movement(over_cutoff_indx));
            active_lever_shuffle_bootstrp_ci = prctile(active_lever_shuffle_bootstrp,[0.5 99.5]);
            active_lever_shuffle_bootstrp_ci_range = abs(active_lever_shuffle - active_lever_shuffle_bootstrp_ci);
            toc
            inactive_lever_shuffle_bootstrp = bootstrp(1000,@mean,lever_shuffle_movement(~over_cutoff_indx));
            inactive_lever_shuffle_bootstrp_ci = prctile(inactive_lever_shuffle_bootstrp,[0.5 99.5]);
            inactive_lever_shuffle_bootstrp_ci_range = abs(inactive_lever_shuffle - inactive_lever_shuffle_bootstrp_ci);
            toc
            
            % criterion for movement correlated
            if active_lever_bootstrp_ci(1) > active_lever_shuffle_bootstrp_ci(2)
                movecorr(curr_cell) = 1;
            elseif active_lever_bootstrp_ci(2) < active_lever_shuffle_bootstrp_ci(1)
                movecorr(curr_cell) = -1;
            end
            
        end
        
        curr_cell/size(roi_trace_df,1)
        
    end
    
    % take raw corelation for strength of relationship
    movecorr_value = zeros(size(movecorr));
    for curr_cell = 1:length(movecorr)
        if movecorr(curr_cell) == 0
            continue
        end
        
        % call the minimum of the raw the low cutoff
        % estimate the noise from smoothed - raw
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 3*mean(noise_est);
        curr_trace_cutoff = zeros(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        curr_trace_cutoff(over_cutoff_indx) = roi_trace_df_smooth(curr_cell,over_cutoff_indx);
        
        temp_corr = corrcoef(curr_trace_cutoff,lever_movement);
        movecorr_value(curr_cell) = temp_corr(2);
    end
    
    %%%
    
    %%% KEEP PROCESSED DATA
    movecorr_all{curr_day} = movecorr;
    movecorr_value_all{curr_day} = movecorr_value;
    
    disp(['Finished movement correlations: day ' day]);
end

% Get rid of all excess variables
clearvars -except animal days curr_day movecorr_all movecorr_value_all


%% Lever force PCA over all days

clc
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 120613 ...
    120614 120615 120616 120617 120618 120619};

%days = {120619};

lever_reward_all = cell(size(days));

for curr_day = 1:length(days)
    
    clearvars -except curr_cell animal days lever_reward_all curr_day
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    % plot around reward
    frames_surround = 30;
    xsg_surround = 10714;
    lever_reward = nan(length(bhv.reward_frames),xsg_surround*2+1);
    for i = 1:length(bhv.reward_frames);
        curr_reward = round(bhv.reward_frames(i));
        curr_reward_xsg = round((curr_reward/bhv.framerate)*10000);
        lever_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        lever_reward_xsg = ...
            curr_reward_xsg-xsg_surround:curr_reward_xsg+xsg_surround;
        if any(lever_reward_frames<=0) || any(lever_reward_frames>size(roi_trace_df,2));
            continue
        end
        %lever_reward(i,:) = bhv.lever_force_resample(lever_reward_frames);
        lever_reward(i,:) = bhv.lever_force(lever_reward_xsg);
    end
    
    lever_reward_all{curr_day} = lever_reward;
    
    disp(['Grabbed day ' day])
end


% plot all mean +/- std
figure;
set(gcf,'Name','Rewarded lever force');
for i = 1:length(lever_reward_all)
    hold on;
    curr_mean = mean(lever_reward_all{i});
    curr_std = std(lever_reward_all{i});
    plot(curr_mean,'color',[0 i/length(lever_reward_all) 0]);
    plot(curr_mean+curr_std,'color',[i/length(lever_reward_all) 0 0])
    plot(curr_mean-curr_std,'color',[i/length(lever_reward_all) 0 0]);
end

% plot all mean +/- std
figure;
set(gcf,'Name','Rewarded lever force');
for i = 1:length(lever_reward_all)
    subplot(7,2,i)
    hold on;
    curr_mean = nanmean(lever_reward_all{i});
    curr_std = nanstd(lever_reward_all{i});
    plot(curr_mean,'k');
    plot(curr_mean+curr_std,'r')
    plot(curr_mean-curr_std,'r');
    line([frames_surround+1 frames_surround+1],ylim,'LineWidth',2,'color','w',...
        'LineStyle','--');
end

% plot all rewarded presses
figure;
set(gcf,'Name','Rewarded lever force');
for i = 1:length(lever_reward_all)
    subplot(7,2,i)
    imagesc(lever_reward_all{i});
    line([frames_surround+1 frames_surround+1],ylim,'LineWidth',2,'color','w',...
        'LineStyle','--');
end

a = vertcat(lever_reward_all{:});
for i = 1:size(a,1)
    a(i,:) = smooth(a(i,:),300);
end
a = diff(a,1,2);

% a = [];
% for i = 1:length(lever_reward_all)
%     a = [a;lever_reward_all{i}-...
%         repmat(mean(lever_reward_all{i}),size(lever_reward_all{i},1),1)];
% end

arange = [9000:12000];
a = a(:,arange);

[coeff score latent] = princomp(a');
pcs = a'*coeff(:,1:4);

% get number of presses in each day
lever_num = cellfun(@(x) size(x,1),lever_reward_all);

% plot them in PCA score space
figure;hold on
for i = 1:length(lever_reward_all)
    pc_scores = lever_reward_all{i}(:,arange)*pcs;
    plot3(pc_scores(:,1),pc_scores(:,2),pc_scores(:,3),'.', ...
        'color',[0 0 i/size(lever_reward_all,1)]);
end

%% Interneuron correlation over days

animal = 'AP71';

% find days
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
    
grab_var = cell(size(days));

for curr_day = 1:length(days)
    
    clearvars -except curr_cell animal days curr_day grab_var kmeans_idx
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    
    smooth_size = 30;
    im.roi_trace_df_smooth = nan(size(im.roi_trace_df));
    for i = 1:size(im.roi_trace_df,1);
        im.roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess');
    end
    
    %     % try a cutoff
    %     % Create a thresholded trace
    %     im.roi_trace_df_smooth_cutoff = zeros(size(im.roi_trace_df_smooth));
    %     for curr_cell = 1:size(im.roi_trace_df_smooth,1)
    %
    %         % call the minimum of the raw the low cutoff
    %         % estimate the noise from smoothed - raw
    %         noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
    %         cutoff = 4*mean(noise_est);
    %         curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
    %         over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
    %         im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
    %             im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    %     end
    
    
    cells = [1:length(roi_labels)];

    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';
    
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    gad_cells(ismember(gad_cells,contaminated)) = [];

    %
    %     % sort by group
    %     roi_trace_kmeans = [im.roi_trace_df_smooth_cutoff(gad_cells,:) kmeans_idx];
    %     roi_trace_kmeans = sortrows(roi_trace_kmeans,size(roi_trace_kmeans,2));
    %     roi_trace_kmeans = roi_trace_kmeans(:,1:end-1);
    %
    %     curr_cc = corrcoef(roi_trace_kmeans');
    curr_cc = corrcoef(im.roi_trace_df_smooth(gad_cells,:)');
    grab_var{curr_day} = curr_cc;
    
    disp(['Finished day: ' day]);
    
end

figure;hold on;
for i = 1:length(grab_var)
    subplot(7,2,i)
    imagesc(grab_var{i});
end


%% Quantify lever activity over days
clear all

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days
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
    
    grab_var.pre_cue_movement= cell(size(days));
    grab_var.response_movement = cell(size(days));
    grab_var.reward_reaction = cell(size(days));
    grab_var.reaction_time = cell(size(days));
    grab_var.cued_movement = cell(size(days));
    grab_var.iti_movement = cell(size(days));
    grab_var.rewarded_lever_force = cell(size(days));
    grab_var.move_reward_time = cell(size(days));
    grab_var.reward_move_ratio = cell(size(days));
    grab_var.reward_move_postcue_time = cell(size(days));
    grab_var.reward_move_precue_time = cell(size(days));
    grab_var.reward_move_postcue = cell(size(days));
    grab_var.iti_time = cell(size(days));
    grab_var.early_movement_ratio = cell(size(days));
    grab_var.rewarded_lever_force = cell(size(days));
    
    for curr_day = 1:length(days)
        % Initialize
        clearvars -except animals curr_animal animal days curr_day grab_var data_path
        
        day = num2str(days{curr_day});
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        grab_var.pre_cue_movement{curr_day} = nan(num_trials,1);
        grab_var.response_movement{curr_day} = nan(num_trials,1);
        grab_var.reward_reaction{curr_day} = nan(num_trials,1);
        grab_var.reaction_time{curr_day} = nan(num_trials,1);
        grab_var.cued_movement{curr_day} = false(num_trials,1);
        grab_var.iti_movement{curr_day} = nan(num_trials,1);
        grab_var.rewarded_lever_force{curr_day} = cell(num_trials,1);
        grab_var.rewarded_lever_force_fixtime{curr_day} = cell(num_trials,1);
        grab_var.move_reward_time{curr_day} = nan(num_trials,1);
        grab_var.reward_move_ratio{curr_day} = nan(num_trials,1);
        grab_var.reward_move_postcue_time{curr_day} = nan(num_trials,1);
        grab_var.reward_move_precue_time{curr_day} = nan(num_trials,1);
        grab_var.reward_move_postcue{curr_day} = nan(num_trials,1);
        grab_var.iti_time{curr_day} = nan(num_trials,1);
        grab_var.early_movement_ratio{curr_day} = nan(num_trials,1);
        % Process
        
        % Only look at rewarded trials
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));
        
        time_surround = 1; % in seconds
        sample_surround = time_surround*xsg_samplerate;
        rewarded_lever_trace = nan(num_trials,sample_surround*2+1);
        rewarded_lever_trace_vel = nan(num_trials,sample_surround*2+1);
        
        % go through all xsg files, grab movement properties
        plots_on = 0;
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            % create movement velocity, active times
            % quantify movement before
            
            [lever_active lever_force_smooth lever_velocity_envelope_smooth]...
                = AP_parseLeverMovement(xsg_data);
            
            lever_stopstart = [0;diff(lever_active)];
            lever_stopstart = [0;diff(lever_active)];
            
            %             move_start_frames = find(lever_stopstart == 1);
            %             move_stop_frames = find(lever_stopstart == -1);
            %             % get rid of unpaired movement bouts at start and end
            %             if move_stop_frames(1) < move_start_frames(1)
            %                 move_stop_frames(1) = [];
            %             end
            %
            %             if move_start_frames(end) > move_stop_frames(end)
            %                 move_start_frames(end) = [];
            %             end
%             if plots_on
%                 figure;
%                 h1 = subplot(2,1,1);
%                 plot(lever_force_smooth,'k')
%                 temp = lever_force_smooth;
%                 temp(~lever_active) = NaN;
%                 hold on; plot(temp,'r');
%                 title('Position')
%                 h2 = subplot(2,1,2); hold on;
%                 plot(lever_velocity_resample_smooth,'m')
%                 plot(lever_velocity_envelope_smooth,'k')
%                 temp = lever_velocity_envelope_smooth;
%                 temp(~lever_active) = NaN;
%                 hold on; plot(temp,'r');
%                 hold on; plot(temp,'r');
%                 line(xlim,[movethresh movethresh],'color','m');
%                 title('Velocity')
%                 linkaxes([h1 h2],'x');
%             end
            
            % get cue and reward times in resampled xsg
            cue_sample = nan(length(curr_trial_list(:,2)),1);
            reward_start_sample = nan(length(curr_trial_list(:,2)),1);
            reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
            for curr_trial_idx = 1:length(curr_trial_list(:,2));
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials) 
                    continue
                end

                % get trial start in dispatcher
                curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
                % get cue and reward sample in xsg (downsampled to ms)
                curr_bhv_cue_sample_rel =  (bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                cue_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_cue_sample_rel);
                
                curr_bhv_reward_start_sample_rel =  (bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_start_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_start_sample_rel);     
                
                curr_bhv_reward_stop_sample_rel =  (bhv{curr_trial}.states.reward(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_stop_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_stop_sample_rel); 
            end
            
            for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
                
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials) || ...
                    length(lever_force_smooth) < cue_sample(curr_trial_idx+1)
                    continue
                end
               
                % ignore if out of bounds
                if cue_sample(curr_trial_idx+1) > length( ...
                        lever_force_smooth)
                    continue
                end
                
                % classify trial cue/movement overlap
                % rewarded-movement/cue overlap can't be > 50ms
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                cued_movement = rewarded_move_start > cue_sample(curr_trial_idx)-50;
%                 cued_movement = ~any(lever_active( ...
%                     cue_sample(curr_trial_idx)-5:cue_sample(curr_trial_idx)));
                grab_var.cued_movement{curr_day}(curr_trial) = cued_movement;
                
                % find movement parameters
                
                % get movement before cues
                pre_cue_sec = 1;
                pre_cue_samples = round(pre_cue_sec*(xsg_samplerate/10));
                if cue_sample(curr_trial_idx)-pre_cue_samples > 0
                    grab_var.pre_cue_movement{curr_day}(curr_trial) = ...
                        sum(lever_velocity_envelope_smooth( ...
                        cue_sample(curr_trial_idx)-pre_cue_samples: ...
                        cue_sample(curr_trial_idx)));
                end
                
                % get movement between cue and reward
                grab_var.response_movement{curr_day}(curr_trial) = ...
                    sum(lever_velocity_envelope_smooth( ...
                    cue_sample(curr_trial_idx):reward_stop_sample(curr_trial_idx))); %...
                    %/(reward_sample(curr_trial_idx)-cue_sample(curr_trial_idx));
                
                % get time between cue and reward
                grab_var.reward_reaction{curr_day}(curr_trial) = ...
                    reward_start_sample(curr_trial_idx) - cue_sample(curr_trial_idx);
                
                % get time to first movement after cue
%                 move_thresh = lever_velocity_envelope_smooth ...
%                     (cue_sample(curr_trial_idx):reward_start_sample(curr_trial_idx)) ...
%                     > movethresh;
%                 if any(move_thresh)
%                     grab_var.reaction_time{curr_day}(curr_trial) = find(move_thresh,1);
%                 end
                grab_var.reaction_time{curr_day}(curr_trial) = ...
                    find(lever_active(cue_sample(curr_trial_idx):...
                    reward_start_sample(curr_trial_idx)),1);
                
                % get movement during the iti
                if curr_trial_idx ~= length(curr_trial_list(:,2)) && ...
                        cue_sample(curr_trial_idx+1) < length( ...
                        lever_velocity_envelope_smooth);
                    grab_var.iti_movement{curr_day}(curr_trial) = ...
                        sum(lever_velocity_envelope_smooth( ...
                        reward_stop_sample(curr_trial_idx): ...
                        cue_sample(curr_trial_idx+1))); % ...
                    % /(cue_sample(curr_trial_idx+1) - ...
                    % reward_sample(curr_trial_idx));
                    
                    % get iti time
                    grab_var.iti_time{curr_day}(curr_trial) = ...
                        cue_sample(curr_trial_idx+1) - ...
                        reward_stop_sample(curr_trial_idx);
                end
                    
                % extract the rewarded movement
                % get start of rewarded movement bout
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                rewarded_move_end = reward_start_sample(curr_trial_idx) + ...
                    find(lever_active(...
                    reward_start_sample(curr_trial_idx):end) == 0,1);
                % if the movement goes until the end of the trace, just
                % take what's available
                if isempty(rewarded_move_end)
                    rewarded_move_end = length(lever_force_smooth);
                end
                grab_var.rewarded_lever_force{curr_day}{curr_trial} =  ...
                    lever_force_smooth(rewarded_move_start: ...
                    rewarded_move_end);
                
                % get lever press for fixed time around reward
                if reward_start_sample(curr_trial_idx)-300>0 && ...
                        reward_start_sample(curr_trial_idx)+300<length(lever_force_smooth);
                    grab_var.rewarded_lever_force_fixtime{curr_day}{curr_trial} = ...
                        lever_force_smooth(reward_start_sample(curr_trial_idx)-300: ...
                        reward_start_sample(curr_trial_idx)+300);
                end
                
                grab_var.move_reward_time{curr_day}(curr_trial) = ...
                    reward_start_sample(curr_trial_idx) - rewarded_move_start;
                
                % get the ratio of post-cue rewarded movements/time to all
                % other movement/time
                if rewarded_move_start > cue_sample(curr_trial_idx)
                    temp_start = rewarded_move_start;
                else
                    temp_start = cue_sample(curr_trial_idx);                   
                end
                if curr_trial_idx < length(cue_sample)
                    if ~isnan(cue_sample(curr_trial_idx+1))
                        temp_reward_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (temp_start:...
                            rewarded_move_end));%/(temp_start-cue_sample(curr_trial_idx));
                        temp_all_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx):...
                            cue_sample(curr_trial_idx+1))); %/(cue_sample(curr_trial_idx+1)- ...
                        %cue_sample(curr_trial_idx));
                        grab_var.reward_move_ratio{curr_day}(curr_trial) = ...
                            temp_reward_move/temp_all_move;
                    end
                end
                
                
                % get movement ratio from fixed early time to all time
                % (median response time of all animals, 490ms, +400 ms water
                % time)
                
                pre_cue_time = 0;
                post_cue_time = 2000; % was previously 890
                if curr_trial_idx > 1 && curr_trial_idx < length(cue_sample)
                    if ~isnan(cue_sample(curr_trial_idx+1))
                        temp_early_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx)-pre_cue_time: ...
                            cue_sample(curr_trial_idx)+post_cue_time));
                        temp_all_move = ...
                            sum(lever_velocity_envelope_smooth ...
                            (cue_sample(curr_trial_idx): ...
                            cue_sample(curr_trial_idx+1)));
                        
                        grab_var.early_movement_ratio{curr_day}(curr_trial) = ...
                            temp_early_move/temp_all_move;
                    end
                end
                
                % If this is needed in the future: decide reward start/end
%                 % get amount of movement before cue/after reward
%                 grab_var.reward_move_postcue{curr_day}(curr_trial) = ...
%                     sum(lever_velocity_envelope_smooth( ...
%                     rewarded_move_end - reward_sample(curr_trial_idx)));
%                 
%                 grab_var.reward_move_postcue_time{curr_day}(curr_trial) = ...
%                     rewarded_move_end - reward_sample(curr_trial_idx);
%                 grab_var.reward_move_precue_time{curr_day}(curr_trial) = ...
%                     cue_sample(curr_trial_idx) - rewarded_move_start;
                          
                
                if plots_on
                    subplot(h1)
                    line([reward_sample_start(curr_trial_idx) reward_sample_start(curr_trial_idx)],ylim);
                    line([reward_sample_stop(curr_trial_idx) reward_sample_stop(curr_trial_idx)],ylim);
                    if cued_movement
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
                    else
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
                    end
                    subplot(h2)
                    line([reward_sample_start(curr_trial_idx) reward_sample_start(curr_trial_idx)],ylim);
                    line([reward_sample_stop(curr_trial_idx) reward_sample_stop(curr_trial_idx)],ylim);
                    if cued_movement
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','g');
                    else
                        line([cue_sample(curr_trial_idx) cue_sample(curr_trial_idx)],ylim,'color','r');
                    end
                end
            end
        end
        disp(['Finished day: ' day]);
    end
    disp('Finished all')
save([animal 'lever_params.mat'],'grab_var');
end
disp('Finished all animals')
 

%% plot things from above

clear all
animals = {'AP76'};
curr_path = '/usr/local/lab/People/Andy/Data/lever_data';
animal_data = cell(length(animals),1);
animal_data_cued_movement = cell(length(animals),1);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    animal_data{curr_animal} = cell(7,1);
    load([curr_path filesep animal 'lever_params.mat']);
    animal_data_cued_movement{curr_animal} = grab_var.cued_movement;
    days = length(grab_var.pre_cue_movement);
    
    figure; set(gcf,'Name',animal);
   
    % plot response movement / all movement
    %     movement_ratio = cellfun(@(x,y) x./(x+y),grab_var.response_movement, ...
    %         grab_var.iti_movement,'UniformOutput',false);
    %     movement_ratio_concat = vertcat(movement_ratio{:});
    %     movement_ratio_smoothed = cellfun(@(x) smooth(x,30),movement_ratio, ...
    %         'UniformOutput',false);
    %     movement_ratio_smoothed_concat = vertcat(movement_ratio_smoothed{:});
    
    movement_ratio_concat = vertcat(grab_var.early_movement_ratio{:});
    movement_ratio_smoothed = cellfun(@(x) smooth(x,30), ...
        grab_var.early_movement_ratio,'UniformOutput',false);
    movement_ratio_smoothed_concat = vertcat(movement_ratio_smoothed{:});
    
    prop_sum = cumsum(cellfun(@length,grab_var.early_movement_ratio))+1;
    movement_ratio_smoothed_concat(prop_sum) = NaN;
    
    subplot(2,1,1); hold on;
     
    plot(movement_ratio_concat,'.k')
    plot(movement_ratio_smoothed_concat,'r','linewidth',2)
    xlim([0 prop_sum(end)]);
    title('movement ratio (fixed time)')
    for i = 1:length(prop_sum)
        line([prop_sum(i) prop_sum(i)],ylim,'color','b')
    end 
    
    % plot response time
    reward_reaction_concat = vertcat(grab_var.reward_reaction{:});
    reward_reaction_smooth = cellfun(@(x) smooth(x,30),...
        grab_var.reward_reaction, 'UniformOutput',false);
    reward_reaction_smooth_concat = vertcat(reward_reaction_smooth{:});
    prop_sum = cumsum(cellfun(@length,grab_var.reward_reaction))+1;
    reward_reaction_smooth_concat(prop_sum) = NaN;
    
    subplot(2,1,2); hold on;
    
    plot(reward_reaction_concat,'.k')
    plot(reward_reaction_smooth_concat,'r','linewidth',2)
    xlim([0 prop_sum(end)]);
    title('response time')
    for i = 1:length(prop_sum)
        line([prop_sum(i) prop_sum(i)],ylim,'color','b')
    end  
    ylim([-50 5000]);
    %saveas(gcf,[curr_path filesep animal '_move_ratio.png'],'png');
    
%     % this plot doesn't make sense, but it looks nice (sometimes)
%     figure; hold on;
%     response_vel = cellfun(@(x,y) x./y,grab_var.response_movement,grab_var.reward_reaction,'UniformOutput',false);
%     iti_vel = cellfun(@(x,y) x./y,grab_var.iti_movement,grab_var.iti_time,'UniformOutput',false);   
%     vel_ratio = cellfun(@(x,y) x./(x+y),response_vel,iti_vel,'UniformOutput',false);
%     vel_ratio_concat = vertcat(vel_ratio{:});
%     vel_ratio_smooth = cellfun(@(x) smooth(x,20),vel_ratio,'UniformOutput',false);
%     vel_ratio_smooth_concat = vertcat(vel_ratio_smooth{:});
%     prop_sum = cumsum(cellfun(@length,vel_ratio))+1;
%     vel_ratio_smooth_concat(prop_sum) = NaN;
%     
%     
%     plot(vel_ratio_concat,'.k')
%     plot(vel_ratio_smooth_concat,'r','linewidth',2)    
%     title('velocity ratio')
%     for i = 1:length(prop_sum)
%         line([prop_sum(i) prop_sum(i)],ylim,'color','b')
%     end  
    
end
% 
% % plot the animals
% clear all
% animals = {'AP61' 'AP62' 'AP63' 'AP64' 'AP65' 'AP66' ...
%     'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
% 
% animal_data = cell(length(animals),1);
% animal_data_cued_movement = cell(length(animals),1);
% for curr_animal = 1:length(animals)
%     animal = animals{curr_animal};
%     animal_data{curr_animal} = cell(7,1);
%     load([animal 'lever_params.mat']);
%     animal_data_cued_movement{curr_animal} = grab_var.cued_movement;
%     days = length(grab_var.pre_cue_movement);
%     figure;
%     for plots = 1:7
%         if plots ~= 7
%             switch plots
%                 case 1
%                     plotprop = grab_var.pre_cue_movement;
%                     name = 'pre cue movement';
%                 case 2
%                     plotprop = grab_var.response_movement;
%                     name = 'response movement';
%                 case 3
%                     plotprop = grab_var.reward_reaction;
%                     name = 'reward reaction';
%                 case 4
%                     plotprop = grab_var.reaction_time;
%                     name = 'reaction time';
%                 case 5
%                     plotprop = grab_var.iti_movement;
%                     name = 'iti movement';
%                 case 6
%                     plotprop = grab_var.move_reward_time;
%                     name = 'move reward time';
%                     
%             end
%             subplot(6,2,plots+((ceil(plots/2)-1)*2))
%             plot_color = [];
%             plot_color = vertcat(grab_var.cued_movement{:});
%             plot_color = [plot_color zeros(length(plot_color),2)];
%             plotdata = vertcat(plotprop{:});
%             scatter(1:length(plotdata),plotdata,[],plot_color,'.')
%             corr_sum = cumsum(cellfun(@length,plotprop))+0.5;
%             for i = 1:length(corr_sum)
%                 line([corr_sum(i) corr_sum(i)],ylim,'color','b')
%             end
%             xlim([0-20 corr_sum(end)+20])
%             ylim([min(vertcat(plotprop{:}))-0.5*nanstd(vertcat(plotprop{:}))...
%                 nanmean(vertcat(plotprop{:}))+4*nanstd(vertcat(plotprop{:}))])
%             title(name);
%             subplot(6,2,plots+((ceil(plots/2))*2))
%             hold on
%             errorbar(cellfun(@(x,y) nanmean(x(~y)), ...
%                 plotprop,grab_var.cued_movement), ...
%                 cellfun(@(x,y) nanstd(x(~y)), ...
%                 plotprop,grab_var.cued_movement),'k');
%             errorbar(cellfun(@(x,y) nanmean(x(y)), ...
%                 plotprop,grab_var.cued_movement), ...
%                 cellfun(@(x,y) nanstd(x(y)), ...
%                 plotprop,grab_var.cued_movement),'r');
%             
%             animal_data{curr_animal}{plots} = plotprop;
%             
%             xlim([0 days+1])
%             alldata = [[cellfun(@(x,y) nanmean(x(y)), ...
%                 plotprop,grab_var.cued_movement) ...
%                 cellfun(@(x,y) nanmean(x(~y)), ...
%                 plotprop,grab_var.cued_movement)] + ...
%                 [cellfun(@(x,y) nanstd(x(y)), ...
%                 plotprop,grab_var.cued_movement) ...
%                 cellfun(@(x,y) nanstd(x(~y)), ...
%                 plotprop,grab_var.cued_movement)] ...
%                 [cellfun(@(x,y) nanmean(x(y)), ...
%                 plotprop,grab_var.cued_movement) ...
%                 cellfun(@(x,y) nanmean(x(~y)), ...
%                 plotprop,grab_var.cued_movement)] - ...
%                 [cellfun(@(x,y) nanstd(x(y)), ...
%                 plotprop,grab_var.cued_movement) ...
%                 cellfun(@(x,y) nanstd(x(~y)), ...
%                 plotprop,grab_var.cued_movement)]];
%             
%             ylim([min(alldata) max(alldata)])
%             
%         else
%             
%             lever_force_xcorr = cell(length(grab_var.rewarded_lever_force),1);
%             for curr_day = 1:length(grab_var.rewarded_lever_force)
%                 % skip if missing data/no data
%                 if sum(grab_var.cued_movement{curr_day}) == 0
%                     continue
%                 end
%                 % pad movements with zeros
%                 move_time = [];
%                 move_time = cellfun(@length,grab_var.rewarded_lever_force{curr_day});
%                 rewarded_lever_force_pad = cellfun(@(x) padarray(x, ...
%                     [max(move_time)+1-length(x) 0],'post'), ...
%                     grab_var.rewarded_lever_force{curr_day},'UniformOutput',false);
%                 % get cross-correlation
%                 max_xcorr = max(xcorr([rewarded_lever_force_pad{...
%                     grab_var.cued_movement{curr_day}}],'coeff'));
%                 lever_force_xcorr{curr_day} = max_xcorr';
%                 disp(['Xcorr''d ' num2str(curr_day)]);
%             end
%             
%             subplot(6,2,plots+((ceil(plots/2)-1)*2))
%             errorbar(cellfun(@(x) mean(x(x~=1)),lever_force_xcorr), ...
%                 cellfun(@std,lever_force_xcorr),'k')
%             xlim([0 days+1])
%             title('Cued-rewarded press xcorr')
%             subplot(6,2,plots+((ceil(plots/2))*2))
%             plot(cellfun(@sum,grab_var.cued_movement),'k')
%             xlim([0 days+1])
%             title('Cued-rewarded press number')
%             
%             if sum(grab_var.cued_movement{curr_day}) > 0
%                 animal_data{curr_animal}{6} = lever_force_xcorr';
%                 animal_data{curr_animal}{7} = cellfun(@sum,grab_var.cued_movement);
%             end
%         end
%     end
%     
%     %saveas(gcf,[animal '_rewarded_cue_movements.png'],'png');
%     saveas(gcf,[animal '_move_rewarded_time.png'],'png');
%     close all
% end
% save('all_animal_plots','animal_data', 'animal_data_cued_movement');
% 
% % Combine days from all animals
% % get the maximum number of days
% max_days = max(cellfun(@(x) length(x),animal_data_cued_movement));
% animal_data_combined_mean = zeros(7,max_days);
% animal_data_combined_std = zeros(7,max_days);
% % get means of all for cued/uncued
% % first unpack all data
% animal_data_concat = cellfun(@(x) cellfun(@(y) ...
%     vertcat(y{:}),x(1:6),'UniformOutput',false), ...
%     animal_data,'UniformOutput',false);
% % next go through and extract cued/uncued (couldn't do cellfun)
% animal_data_cued_norm = cell(size(animal_data_concat));
% animal_data_uncued_norm = cell(size(animal_data_concat));
% for curr_mouse = 1:length(animal_data_concat)
%     cue_indx = vertcat(animal_data_cued_movement{curr_mouse}{:})>0;
%     for curr_condition = 1:5
%         
%        temp_cued_mean = nanmean( ...
%             animal_data_concat{curr_mouse}{curr_condition});
%         temp_uncued_mean = nanmean( ...
%             animal_data_concat{curr_mouse}{curr_condition});
%         
%         animal_data_cued_norm{curr_mouse}{curr_condition} = ...
%             cellfun(@(x,y) x(y)/temp_cued_mean, ...
%             animal_data{curr_mouse}{curr_condition}, ...
%             animal_data_cued_movement{curr_mouse},'UniformOutput',false);
%         animal_data_uncued_norm{curr_mouse}{curr_condition} = ...
%             cellfun(@(x,y) x(~y)/temp_uncued_mean, ...
%             animal_data{curr_mouse}{curr_condition}, ...
%             animal_data_cued_movement{curr_mouse},'UniformOutput',false);
%         
%         % set up dummy cell so they're all the same length
%         animal_data_cued_norm{curr_mouse}{curr_condition}{max_days+1} = [];
%         animal_data_uncued_norm{curr_mouse}{curr_condition}{max_days+1} = [];
%     end
% end
% 
% % collect all days
% day_data_cued_concat = cell(max_days,1);
% day_data_uncued_concat = cell(max_days,1);
% for curr_day = 1:max_days
%     day_data_cued_concat{curr_day} = cell(5,1);
%     day_data_uncued_concat{curr_day} = cell(5,1);
%     for curr_condition = 1:5
%        temp_combined_cued = cellfun(@(x) x{curr_condition}{curr_day}, ...
%            animal_data_cued_norm,'UniformOutput',false);
%        day_data_cued_concat{curr_day}{curr_condition} = vertcat(temp_combined_cued{:});
%         
%        temp_combined_uncued = cellfun(@(x) x{curr_condition}{curr_day}, ...
%            animal_data_uncued_norm,'UniformOutput',false);
%        day_data_uncued_concat{curr_day}{curr_condition} = vertcat(temp_combined_uncued{:});
%     end
% end
% 
% % get mean/std of all days
% all_cued_mean = nan(max_days,5);
% all_cued_std = nan(max_days,5);
% all_uncued_mean = nan(max_days,5);
% all_uncued_std = nan(max_days,5);
% for curr_day = 1:max_days
%     for curr_condition = 1:5
%         all_cued_mean(curr_day,curr_condition) = nanmean( ...
%             day_data_cued_concat{curr_day}{curr_condition});
%         all_cued_std(curr_day,curr_condition) = nanstd( ...
%             day_data_cued_concat{curr_day}{curr_condition});
%         all_uncued_mean(curr_day,curr_condition) = nanmean( ...
%             day_data_uncued_concat{curr_day}{curr_condition});
%         all_uncued_std(curr_day,curr_condition) = nanstd( ...
%             day_data_uncued_concat{curr_day}{curr_condition});
%     end
% end
% 
% figure
% for i = 1:5
%     subplot(5,1,i)
%     hold on
%     errorbar(all_uncued_mean(:,i),all_uncued_std(:,i),'k')
%     errorbar(all_cued_mean(:,i),all_cued_std(:,i),'r')
%     
%     switch i
%                 case 1
%                     name = 'pre cue movement';
%                 case 2
%                     name = 'response movement';
%                 case 3
%                     name = 'reward reaction';
%                 case 4
%                     name = 'reaction time';
%                 case 5
%                     name = 'iti movement';
%     end
%     title(name);
% end

%% Get relative activity of pyramidal cells and interneurons
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    % Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    % Process
    smooth_size = 30;
    roi_trace_df_smooth = nan(size(roi_trace_df));
    for i = 1:size(roi_trace_df,1);
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
    end
    
    % Create a thresholded trace
    roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
    for curr_cell = 1:size(roi_trace_df_smooth,1)
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 4*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    
    % Only look at rewarded trials
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
    good_trials = find(rewarded_trials & bhv.imaged_trials);
    reward_frames_good = cellfun(@(x) x.states.reward(1), ...
        bhv.bhv_frames(good_trials));
    cue_frames_good = cellfun(@(x) x.states.cue(1), ...
        bhv.bhv_frames(good_trials));
    
    % compensate for loops
    reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
    cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
    
    % Initialize vars to save
    cue_reward_cross = nan(size(roi_trace_df,1),length(cue_frames_good));
    cue_reward_max = nan(size(roi_trace_df,1),length(cue_frames_good));
    for curr_cell = 1:size(roi_trace_df,1);
        frames_back = 5;
        frames_forward = 40;
        frames_surround = frames_back+frames_forward+1;
        
        % plot around cue
        cue_reward = nan(length(cue_frames_good),frames_surround);
        for i = 1:length(cue_frames_good);
            curr_reward = round(bhv.cue_frames(i));
            roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
                continue
            end
            cue_reward(i,:) = roi_trace_df_smooth_cutoff(curr_cell,roi_reward_frames);
            if any(cue_reward(i,:) > 0);
                cue_reward_cross(curr_cell,i) = find(cue_reward(i,:) > 0,1);
                cue_reward_max(curr_cell,i) = max(cue_reward(i,:));
            end
        end
    end
    
    % get start times relative to rewards
    lever_reward_cross = cue_reward_cross - ...
        repmat([reward_frames_good-cue_frames_good]',size(cue_reward_cross,1),1);
    
    % pick cell alignment by which has smaller start sem
    cue_reward_sem = nanstd(cue_reward_cross,[],2)./ ...
        sqrt(sum(cue_reward_cross > 0,2));
    lever_reward_sem = nanstd(lever_reward_cross,[],2)./ ...
        sqrt(sum(lever_reward_cross > 0,2));
    % take out cells that don't fire on enough trials
    sparse_cells = sum(cue_reward_cross > 0,2) < 10;
    cue_cells = cue_reward_sem < lever_reward_sem & ~sparse_cells;
    lever_cells = lever_reward_sem < cue_reward_sem & ~sparse_cells;
    unaligned_cells = ~cue_cells & ~lever_cells | sparse_cells;
    
    %%% For TK: activity of pyramidal cells vs. gad cells
    cells = [1:size(roi_trace_df,1)];
    gad_cells = fix this;
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    cue_reward_bin = zeros(size(cue_reward_cross));
    cue_reward_bin(cue_reward_cross > 0) = 1;
    
    pyr_num = sum(cue_reward_bin(pyr_cells,:))/length(pyr_cells);
    gad_num = sum(cue_reward_bin(gad_cells,:))/length(gad_cells);
    %%%
    
    % Save
    grab_var.pyr_num{curr_day} = pyr_num;
    grab_var.gad_num{curr_day} = gad_num;
    
    disp(['Finished day: ' day]);
    
end

figure
for i = 1:14
    subplot(7,2,i)
    plot(grab_var.pyr_num{i},grab_var.gad_num{i},'.')
end


%% Get max df/f after every trial for all cells over days
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    % Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    % Process
    smooth_size = 30;
    roi_trace_df_smooth = nan(size(roi_trace_df));
    for i = 1:size(roi_trace_df,1);
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
    end
    
    % Create a thresholded trace
    roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
    for curr_cell = 1:size(roi_trace_df_smooth,1)
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = 4*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    
    % Only look at rewarded trials
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
    good_trials = find(rewarded_trials & bhv.imaged_trials);
    reward_frames_good = cellfun(@(x) x.states.reward(1), ...
        bhv.bhv_frames(good_trials));
    cue_frames_good = cellfun(@(x) x.states.cue(1), ...
        bhv.bhv_frames(good_trials));
    
    % compensate for loops
    reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
    cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
    
    % Initialize vars to save
    cue_reward_cross = nan(size(roi_trace_df,1),length(cue_frames_good));
    cue_reward_max = zeros(size(roi_trace_df,1),length(cue_frames_good));
    for curr_cell = 1:size(roi_trace_df,1);
        frames_back = 0;
        frames_forward = 20;
        frames_surround = frames_back+frames_forward+1;
        
        % plot around cue
        cue_reward = nan(length(cue_frames_good),frames_surround);
        for i = 1:length(cue_frames_good);
            curr_reward = round(bhv.cue_frames(i));
            roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
                continue
            end
            cue_reward(i,:) = roi_trace_df_smooth_cutoff(curr_cell,roi_reward_frames);
            if any(cue_reward(i,:) > 0);
                cue_reward_cross(curr_cell,i) = find(cue_reward(i,:) > 0,1);
                cue_reward_max(curr_cell,i) = max(cue_reward(i,:));
            end
        end
    end
    
    
    % Save
    grab_var{curr_day} = cue_reward_max;
    
    disp(['Finished day: ' day]);
    
end

figure
for i = 1:14
    subplot(7,2,i)
    imagesc(corrcoef(grab_var{i}))
end



%% Pyr-Gad distance relationship
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    % Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    load([day '_' animal '_summed_50_template.roi'],'-MAT');
    
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for i  = 1:num_rois
        [area roi_center_x(i) roi_center_y(i)] = ...
            polycenter(polygon.ROI{i}(:,1),polygon.ROI{i}(:,2));
    end
    
    % create matrix of distances
    roi_dist = nan(length(roi_center_x));
    for i = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(i) - roi_center_x(j))^2 + ...
                (roi_center_y(i) - roi_center_y(j))^2);
            roi_dist(i,j) = curr_dist;
        end
    end
    
    % Process
    
    smooth_size = 30;
    roi_trace_df_smooth = nan(size(roi_trace_df));
    for i = 1:size(roi_trace_df,1);
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
    end
    
    cells = [1:size(roi_trace_df,1)];
    
    % new estimate
    gad_cells = fix this;
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    % put gad cells in order: PV (1-13), SOM (end)
    pv = [109   111    21    40    51   114    99   119  ...
        97    85    11    47   118];
    som = [29    12    49 ...
        192    96   204    76   101    55     5   134   102];
    
    % do oopsi to guess at spikeish regions
    n_best = {};
    P_best = {};
    V = {};
    V.dt = 1/bhv.framerate;
    V.fast_poiss = 1;
    for i = 1:size(roi_trace_df,1)
        [n_best{i} P_best{i} V_est]=fast_oopsi(roi_trace_df(i,:),V);
        %[n_best{i} P_best{i} V_est]=fast_oopsi(roi_trace_df(i,:),V_est);
        disp(num2str(i/size(roi_trace_df,1)));
    end
    n_best = [n_best{:}]';
    
    n_best_thresh = n_best;
    n_best_thresh(n_best_thresh < 0.2) = 0;
    n_best_thresh(:,1:5) = 0;
    
    n_best_thresh_pyr = sum(n_best_thresh(pyr_cells,:));
    n_best_thresh_gad = sum(n_best_thresh(gad_cells,:));
    
    % get corrcoef of each gad/som with pyramidal cells in different radii
    radii = 0:10:500;
    gad_rad_som = zeros(length(radii),length(som));
    for radius_idx = 1:length(radii);
        radius = radii(radius_idx);
        for curr_gad = 1:size(gad_rad_som,2)
            curr_som = som(curr_gad);
            curr_pyr_radii = find(roi_dist(curr_som,pyr_cells) <= radius);
            curr_pyr_act = sum(n_best_thresh(pyr_cells(curr_pyr_radii),:),1);
            curr_corrcoef = corrcoef(curr_pyr_act, ...
                roi_trace_df_smooth(curr_som,:));
            gad_rad_som(radius_idx,curr_gad) = curr_corrcoef(2);
        end
    end
    
    gad_rad_pv = zeros(length(radii),length(pv));
    for radius_idx = 1:length(radii);
        radius = radii(radius_idx);
        for curr_gad = 1:size(gad_rad_pv,2)
            curr_pv = pv(curr_gad);
            curr_pyr_radii = find(roi_dist(curr_pv,pyr_cells) <= radius);
            curr_pyr_act = sum(n_best_thresh(pyr_cells(curr_pyr_radii),:),1);
            curr_corrcoef = corrcoef(curr_pyr_act, ...
                roi_trace_df_smooth(curr_pv,:));
            gad_rad_pv(radius_idx,curr_gad) = curr_corrcoef(2);
        end
    end
    
    % Save
    grab_var{curr_day}.gad_rad_som = gad_rad_som;
    grab_var{curr_day}.gad_rad_pv = gad_rad_pv;
    
    disp(['Finished day: ' day]);
    
end


figure
for i = 1:14
    subplot(7,2,i); hold on;
    errorbar(nanmean(grab_var{i}.gad_rad_som,2), ...
        2*nanstd(grab_var{i}.gad_rad_som,[],2)/...
        sqrt(size(grab_var{i}.gad_rad_som,2)));
    errorbar(nanmean(grab_var{i}.gad_rad_pv,2), ...
        2*nanstd(grab_var{i}.gad_rad_pv,[],2)/...
        sqrt(size(grab_var{i}.gad_rad_pv,2)),'r');
end


%% TK length constant by activity
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    % Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    roi_name = [day '_' animal '_summed_50_template.roi'];
    load(roi_name,'-MAT');
    
    % Process
    
    cells = [1:size(roi_trace_df,1)];
    
    % define cells
    gad_cells = fix this;
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    % get rid of ones in corner, and define potential pv/som
    pv = [109   111    21    40    51   114    99   119  ...
        97    85    11    47   118];
    som = [29    12    49 ...
        192    96   204    76   101    55     5   134   102];
    gad_cells = [pv som];
    
    % loess smooth trace
    smooth_size = 30;
    roi_trace_df_smooth = nan(size(roi_trace_df));
    for i = 1:size(roi_trace_df,1);
        roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
    end
    
    % Create a thresholded trace (seperate thresh for gad/pyr)
    pyr_cutoff = 5;
    gad_cutoff = 2;
    
    roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
    for curr_pyr = 1:length(pyr_cells)
        curr_cell = pyr_cells(curr_pyr);
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = pyr_cutoff*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    for curr_gad = 1:length(gad_cells)
        curr_cell = gad_cells(curr_gad);
        noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
        cutoff = gad_cutoff*mean(noise_est);
        curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
        over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
        roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
            roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    end
    
    % Only look at rewarded trials
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
    good_trials = find(rewarded_trials & bhv.imaged_trials);
    reward_frames_good = cellfun(@(x) x.states.reward(1), ...
        bhv.bhv_frames(good_trials));
    cue_frames_good = cellfun(@(x) x.states.cue(1), ...
        bhv.bhv_frames(good_trials));
    
    % compensate for loops
    reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
    cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
    
    % find active cells around the cue
    cue_reward_cross = nan(size(roi_trace_df,1),length(cue_frames_good));
    cue_reward_max = nan(size(roi_trace_df,1),length(cue_frames_good));
    cue_reward_trace = cell(size(roi_trace_df,1),length(cue_frames_good));
    for curr_cell = 1:size(roi_trace_df,1);
        frames_back = 10;
        frames_forward = 40;
        frames_surround = frames_back+frames_forward+1;
        
        cue_reward = nan(length(cue_frames_good),frames_surround);
        for i = 1:length(cue_frames_good);
            curr_reward = round(bhv.cue_frames(i));
            roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
                continue
            end
            cue_reward(i,:) = roi_trace_df_smooth_cutoff(curr_cell,roi_reward_frames);
            if any(cue_reward(i,:) > 0);
                cue_reward_cross(curr_cell,i) = find(cue_reward(i,:) > 0,1);
                cue_reward_max(curr_cell,i) = max(cue_reward(i,:));
            end
            % save the whole trace
            cue_reward_trace{curr_cell,i} = cue_reward(i,:);
        end
    end
    
    % binarize activity around cue
    cue_reward_bin = zeros(size(cue_reward_cross));
    cue_reward_bin(cue_reward_cross > 0) = 1;
    
    pyr_num = sum(cue_reward_bin(pyr_cells,:))/length(pyr_cells);
    gad_num = sum(cue_reward_bin(gad_cells,:))/length(gad_cells);
    
    % create matrix of distances
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for i  = 1:num_rois
        [area roi_center_x(i) roi_center_y(i)] = ...
            polycenter(polygon.ROI{i}(:,1),polygon.ROI{i}(:,2));
    end
    
    roi_dist = nan(length(roi_center_x));
    for i = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(i) - roi_center_x(j))^2 + ...
                (roi_center_y(i) - roi_center_y(j))^2);
            roi_dist(i,j) = curr_dist;
        end
    end
    
    % get correlation matrix between cells
    cue_reward_bin_corr = corrcoef(cue_reward_);
    % set diagonal to NaN
    cue_reward_bin_corr(eye(length(cue_reward_bin_corr)) ~= 0) = NaN;
    
    % prepare bins
    % num_bins = 50;
    % binsize = linspace(min(roi_dist(:)),max(roi_dist(:)),num_bins);
    binsize = 0:10:200;
    num_bins = length(binsize);
    
    % pyr-pyr pairs
    pyrpyr_corr = cue_reward_bin_corr(pyr_cells,pyr_cells);
    pyrpyr_dist = roi_dist(pyr_cells,pyr_cells);
    
    [n,pyrpyr_bin] = histc(pyrpyr_dist(:),binsize);
    pyrpyr_binMean = [];
    pyrpyr_binSem = [];
    for i = 1:num_bins
        flagBinMembers = (pyrpyr_bin == i);
        binMembers = pyrpyr_corr(flagBinMembers);
        pyrpyr_binMean(i) = nanmean(binMembers);
        pyrpyr_binSem(i) = nanstd(binMembers)/sqrt(length(binMembers));
    end
    
    % gad-gad pairs
    gadgad_corr = cue_reward_bin_corr(gad_cells,gad_cells);
    gadgad_dist = roi_dist(gad_cells,gad_cells);
    
    [n,gadgad_bin] = histc(gadgad_dist(:),binsize);
    gadgad_binMean = [];
    gadgad_binSem = [];
    for i = 1:num_bins
        flagBinMembers = (gadgad_bin == i);
        binMembers = gadgad_corr(flagBinMembers);
        gadgad_binMean(i) = nanmean(binMembers);
        gadgad_binSem(i) = nanstd(binMembers)/sqrt(length(binMembers));
    end
    
    % pyr-gad pairs
    pyrgad_corr = cue_reward_bin_corr(pyr_cells,gad_cells);
    pyrgad_dist = roi_dist(pyr_cells,gad_cells);
    
    [n,pyrgad_bin] = histc(pyrgad_dist(:),binsize);
    pyrgad_binMean = [];
    pyrgad_binSem = [];
    for i = 1:num_bins
        flagBinMembers = (pyrgad_bin == i);
        binMembers = pyrgad_corr(flagBinMembers);
        pyrgad_binMean(i) = nanmean(binMembers);
        pyrgad_binSem(i) = nanstd(binMembers)/sqrt(length(binMembers));
    end
    
    
    % Save
    grab_var{curr_day}.pyrpyr_binMean = pyrpyr_binMean;
    grab_var{curr_day}.pyrgad_binMean = pyrgad_binMean;
    grab_var{curr_day}.gadgad_binMean = gadgad_binMean;
    
    grab_var{curr_day}.pyrpyr_binSem = pyrpyr_binSem;
    grab_var{curr_day}.pyrgad_binSem = pyrgad_binSem;
    grab_var{curr_day}.gadgad_binSem = gadgad_binSem;
    
    disp(['Finished day: ' day]);
    
end

% figure;hold on;
% for i = 1:14
%     plot(grab_var{i}.pyrgad_binMean,'color',[i/14 0 0]);
% end
%
% figure;hold on;
% for i = 1:14
%     plot(grab_var{i}.pyrpyr_binMean,'color',[i/14 0 0]);
% end

%% Classify cells over days: peri-event
% this isn't working properly with oopsi
clear all

animal = 'AP71';

% find days
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

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    
    %%%% Process
    
    % do twice: once with threshold, once with deconv
    for thresh_deconv = 1:2
        
        if thresh_deconv == 1
            smooth_size = 30;
            im.roi_trace_df_smooth = nan(size(im.roi_trace_df));
            for i = 1:size(im.roi_trace_df,1);
                im.roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess');
            end
            
            % Create a thresholded trace
            im.roi_trace_df_smooth_cutoff = zeros(size(im.roi_trace_df_smooth));
            for curr_cell = 1:size(im.roi_trace_df_smooth,1)
                noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
                cutoff = 4*mean(noise_est);
                curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
                over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
                im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
                    im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
            end
            % cut off beginning
            im.roi_trace_df_smooth_cutoff(:,1:10) = 0;
        elseif thresh_deconv == 2
            im.roi_trace_df_smooth_cutoff = im.oopsi;
            im.roi_trace_df_smooth_cutoff(im.oopsi < 0.2) = 0;
        end
        
        % Only look at rewarded trials
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
        good_trials = find(rewarded_trials & bhv.imaged_trials);
        reward_frames_good = cellfun(@(x) x.states.reward(1), ...
            bhv.bhv_frames(good_trials));
        cue_frames_good = cellfun(@(x) x.states.cue(1), ...
            bhv.bhv_frames(good_trials));
        
        % compensate for loops
        reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
        cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
        
        % make another category: pre-cue (like reward vs. cue)
        median_trial_time = median(reward_frames_good - cue_frames_good);
        precue_frames_good = cue_frames_good - median_trial_time;
        
        % round for indexing
        reward_frames_good_round = round(reward_frames_good);
        cue_frames_good_round = round(cue_frames_good);
        precue_frames_good_round = round(precue_frames_good);
        
        % make a matrix of trial times for purposes of searching for transients
        trial_time = reward_frames_good - cue_frames_good;
        
        % define the cells to look at
        cells = [1:size(im.roi_trace_df,1)];
        
        gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
        pyr_cells = cells(~ismember(cells,gad_cells));
        
        contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels));
        
        pyr_cells(ismember(pyr_cells,contaminated)) = [];
        gad_cells(ismember(gad_cells,contaminated)) = [];
        
        frames_back = 1;
        
        num_rois = size(im.roi_trace_df,1);
        num_frames = size(im.roi_trace_df,2);
        
        precue_active = zeros(num_rois,length(good_trials));
        cue_active = zeros(num_rois,length(good_trials));
        reward_active = zeros(num_rois,length(good_trials));
        
        cue_preactive = zeros(num_rois,length(good_trials));
        reward_preactive = zeros(num_rois,length(good_trials));
        
        cue_aligned_trials = zeros(num_rois,length(good_trials));
        reward_aligned_trials = zeros(num_rois,length(good_trials));
        
        cue_alignment_indx = zeros(size(im.roi_trace_df,1),1);
        reward_alignment_indx = zeros(size(im.roi_trace_df,1),1);
        for curr_cell = [pyr_cells gad_cells];
            
            % go through valid trials, check for activity
            for i = 1:length(good_trials)
                % find if activity after trial time
                frames_forward = trial_time(i);
                if frames_forward + reward_frames_good_round(i) > num_frames || ...
                        precue_frames_good_round(i) < 1;
                    continue
                end
                precue_frames_search = precue_frames_good_round(i) : ...
                    precue_frames_good_round(i) + frames_forward;
                cue_frames_search = cue_frames_good_round(i): ...
                    cue_frames_good_round(i) + frames_forward;
                reward_frames_search = reward_frames_good_round(i): ...
                    reward_frames_good_round(i) + frames_forward;
                
                precue_active(curr_cell,i) = ...
                    any(im.roi_trace_df_smooth_cutoff(curr_cell,precue_frames_search));
                cue_active(curr_cell,i) = ...
                    any(im.roi_trace_df_smooth_cutoff(curr_cell,cue_frames_search));
                reward_active(curr_cell,i) = ...
                    any(im.roi_trace_df_smooth_cutoff(curr_cell,reward_frames_search));
                
                % find if activity before alignment points
                frames_forward = trial_time(i);
                cue_preframes_search = cue_frames_good_round(i) - frames_back: ...
                    cue_frames_good_round(i);
                reward_preframes_search = reward_frames_good_round(i) - frames_back: ...
                    reward_frames_good_round(i);
                
                cue_preactive(curr_cell,i) = ...
                    any(im.roi_trace_df_smooth_cutoff(curr_cell,cue_preframes_search));
                reward_preactive(curr_cell,i) = ...
                    any(im.roi_trace_df_smooth_cutoff(curr_cell,reward_preframes_search));
            end
            
            % to classify: cross thresh in given time period, not cross thresh in
            % the time before, not cross thresh in the other given time period
            cue_aligned_trials(curr_cell,:) = cue_active(curr_cell,:) & ...
                ~precue_active(curr_cell,:) & ~cue_preactive(curr_cell,:);
            reward_aligned_trials(curr_cell,:) = reward_active(curr_cell,:) & ...
                ~cue_active(curr_cell,:) & ~reward_preactive(curr_cell,:);
            
            % calculate alignment index
            num_frames = size(im.roi_trace_df_smooth_cutoff,2);
            num_trials = length(cue_frames_good);
            % define chance overlap of calcium transient with surround frames
            thresh_start = diff(im.roi_trace_df_smooth_cutoff(curr_cell,:) > 0);
            if im.roi_trace_df_smooth_cutoff(curr_cell,1) > 0;
                thresh_start(find(thresh_start < 0,1)) = 0;
            end
            if im.roi_trace_df_smooth_cutoff(curr_cell,end) > 0;
                thresh_start(find(thresh_start > 0,1,'last')) = 0;
            end
            thresh_run = find(thresh_start == -1)-find(thresh_start == 1);
            mean_thresh_run = mean(thresh_run);
            thresh_freq = length(thresh_run)/num_frames;
            p_overlap_align = thresh_freq*(median_trial_time+mean_thresh_run);
            p_nonoverlap_prealign = (1-thresh_freq)^...
                ((median_trial_time > mean_thresh_run)*median_trial_time + ...
                (mean_thresh_run > median_trial_time)*mean_thresh_run);
            chance_align = p_overlap_align*p_nonoverlap_prealign*num_trials;
            
            % calculate alignment index
            reward_alignment_indx(curr_cell) = (sum(reward_aligned_trials ...
                (curr_cell,:)) - chance_align)/num_trials;
            cue_alignment_indx(curr_cell) = (sum(cue_aligned_trials ...
                (curr_cell,:)) - chance_align)/num_trials;
            
            
        end
        
        if thresh_deconv == 1
            reward_alignment_thresh_indx = reward_alignment_indx;
            cue_alignment_thresh_indx = cue_alignment_indx;
        elseif thresh_deconv == 2
            reward_alignment_oopsi_indx = reward_alignment_indx;
            cue_alignment_oopsi_indx = cue_alignment_indx;
        end
        
    end
    
    %%%% Save
    grab_var.reward_alignment_thresh_indx{curr_day} =  reward_alignment_thresh_indx;
    grab_var.cue_alignment_thresh_indx{curr_day} = cue_alignment_thresh_indx;
    
    grab_var.reward_alignment_oopsi_indx{curr_day} =  reward_alignment_oopsi_indx;
    grab_var.cue_alignment_oopsi_indx{curr_day} = cue_alignment_oopsi_indx;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);



%% Classify cells over days: deconv/xcorr
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    %%%% Process
    
    % at the moment: just assign a value for each condition to cells
    
    %%% Find movement times
    
    % define velocity threshold for active movement
    movethresh = 0.02;
    
    % get velocity envelope, smooth over 5 frames (~150ms)
    lever_velocity = [0;diff(bhv.lever_force_resample)];
    lever_velocity_hilbert = hilbert(lever_velocity);
    lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
        conj(lever_velocity_hilbert));
    lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,5);
    
    % define active parts of the lever, fill in gaps of < 1 second (30 frames)
    lever_active = lever_velocity_envelope_smooth > movethresh;
    lever_active_switch = [0;diff(lever_active)];
    lever_active_edges = find(lever_active_switch);
    lever_active_spaces = [0;diff(lever_active_switch(lever_active_edges))];
    lever_active_up = find(lever_active_spaces == 2);
    lever_active_down = lever_active_up-1;
    lever_active_space_dist = lever_active_edges(lever_active_up) - ...
        lever_active_edges(lever_active_down);
    lever_active_fill = lever_active_space_dist < 30;
    for i = find(lever_active_fill)'
        lever_active(lever_active_edges(lever_active_down(i)): ...
            lever_active_edges(lever_active_up(i))) = 1;
    end
    
    lever_stopstart = [0;diff(lever_active)];
    move_start_frames = find(lever_stopstart == 1);
    move_stop_frames = find(lever_stopstart == -1);
    
    
    % Spread out movement times for surrounding window
    spread_frames = 15;
    spread_frames_filter = ones(1,spread_frames);
    maxlags = 40;
    
    move_start_trace = zeros(1,size(roi_trace_df,2));
    move_start_trace(move_start_frames) = 1;
    move_start_trace = conv(move_start_trace,spread_frames_filter,'same');
    
    move_stop_trace = zeros(1,size(roi_trace_df,2));
    move_stop_trace(move_stop_frames) = 1;
    move_stop_trace = conv(move_stop_trace,spread_frames_filter,'same');
    
    %%% Find the maximum xcorr for each cell for movement start/stop
    
    move_start_max_xcorr = zeros(size(roi_trace_df,1),1);
    move_start_max_xcorr_lag = zeros(size(roi_trace_df,1),1);
    move_stop_max_xcorr = zeros(size(roi_trace_df,1),1);
    move_stop_max_xcorr_lag = zeros(size(roi_trace_df,1),1);
    for curr_cell = 1:size(roi_trace_df,1);
        [temp_xcorr temp_xcorr_lags] = ...
            xcorr(im.oopsi(curr_cell,:),move_start_trace,maxlags);
        move_start_max_xcorr(curr_cell) = max(temp_xcorr);
        move_start_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
            find(temp_xcorr == max(temp_xcorr)));
        
        [temp_xcorr temp_xcorr_lags] = ...
            xcorr(im.oopsi(curr_cell,:),move_stop_trace,maxlags);
        move_stop_max_xcorr(curr_cell) = max(temp_xcorr);
        move_stop_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
            find(temp_xcorr == max(temp_xcorr)));
        disp(curr_cell/size(roi_trace_df,1))
    end
    
    %%% Find the spiking within movement/non-movement times
    move_dot = im.oopsi*lever_active;
    nomove_dot = im.oopsi*~lever_active;
    
    
    %%% Find cue/reward times
    
    % Only look at rewarded trials
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
    good_trials = find(rewarded_trials & bhv.imaged_trials);
    reward_frames_good = cellfun(@(x) x.states.reward(1), ...
        bhv.bhv_frames(good_trials));
    cue_frames_good = cellfun(@(x) x.states.cue(1), ...
        bhv.bhv_frames(good_trials));
    
    % compensate for loops
    reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
    cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
    
    % round for indexing
    reward_frames_good_round = round(reward_frames_good);
    cue_frames_good_round = round(cue_frames_good);
    
    % Spread out movement times for surrounding window
    cue_trace = zeros(1,size(roi_trace_df,2));
    cue_trace(cue_frames_good_round) = 1;
    cue_trace = conv(cue_trace,spread_frames_filter,'same');
    
    reward_trace = zeros(1,size(roi_trace_df,2));
    reward_trace(reward_frames_good_round) = 1;
    reward_trace = conv(reward_trace,spread_frames_filter,'same');
    
    %%% Find the maximum xcorr for each cell for cue/reward
    
    cue_max_xcorr = zeros(size(roi_trace_df,1),1);
    cue_max_xcorr_lag = zeros(size(roi_trace_df,1),1);
    reward_max_xcorr = zeros(size(roi_trace_df,1),1);
    reward_max_xcorr_lag = zeros(size(roi_trace_df,1),1);
    for curr_cell = 1:size(roi_trace_df,1);
        [temp_xcorr temp_xcorr_lags] = ...
            xcorr(im.oopsi(curr_cell,:),cue_trace,maxlags);
        cue_max_xcorr(curr_cell) = max(temp_xcorr);
        cue_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
            find(temp_xcorr == max(temp_xcorr)));
        
        [temp_xcorr temp_xcorr_lags] = ...
            xcorr(im.oopsi(curr_cell,:),reward_trace,maxlags);
        reward_max_xcorr(curr_cell) = max(temp_xcorr);
        reward_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
            find(temp_xcorr == max(temp_xcorr)));
        disp(curr_cell/size(roi_trace_df,1))
    end
    
    %%% Normalize maximum xcorr norms
    move_start_max_xcorr_norm = move_start_max_xcorr./ ...
        (length(move_start_frames)*sum(im.oopsi,2));
    move_stop_max_xcorr_norm = move_stop_max_xcorr./ ...
        (length(move_stop_frames)*sum(im.oopsi,2));
    cue_max_xcorr_norm = cue_max_xcorr./ ...
        (length(cue_frames_good)*sum(im.oopsi,2));
    reward_max_xcorr_norm = reward_max_xcorr./ ...
        (length(reward_frames_good)*sum(im.oopsi,2));
    move_dot_norm = move_dot./ ...
        (sum(lever_active)*sum(im.oopsi,2));
    nomove_dot_norm = nomove_dot./ ...
        (sum(~lever_active)*sum(im.oopsi,2));
    
    % make TK-style index: normalized above-chance coincidence
    %%% Make tk correlation index type measure
    % move start
    num_frames = size(roi_trace_df,2);
    move_start_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = move_start_max_xcorr;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(move_start_trace,2)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(move_start_trace,2)/num_frames)');
    move_start_tk_ci = (coincident_events - chance_events)./mean_events;
    
    % move stop
    num_frames = size(roi_trace_df,2);
    move_stop_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = move_stop_max_xcorr;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(move_stop_trace,2)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(move_stop_trace,2)/num_frames)');
    move_stop_tk_ci = (coincident_events - chance_events)./mean_events;
    
    % move times
    num_frames = size(roi_trace_df,2);
    move_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = move_dot;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(lever_active)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(lever_active)/num_frames)');
    move_tk_ci = (coincident_events - chance_events)./mean_events;
    
    % nomove times
    num_frames = size(roi_trace_df,2);
    nomove_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = nomove_dot;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(~lever_active)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(~lever_active)/num_frames)');
    nomove_tk_ci = (coincident_events - chance_events)./mean_events;
    
    % cue
    num_frames = size(roi_trace_df,2);
    cue_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = cue_max_xcorr;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(cue_trace,2)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(cue_trace,2)/num_frames)');
    cue_tk_ci = (coincident_events - chance_events)./mean_events;
    
    % reward
    num_frames = size(roi_trace_df,2);
    reward_tk_ci = zeros(size(roi_trace_df,1));
    
    coincident_events = reward_max_xcorr;
    chance_events = num_frames*(sum(im.oopsi,2)/num_frames)* ...
        (sum(reward_trace,2)/num_frames)';
    mean_events =   num_frames*sqrt((sum(im.oopsi,2)/num_frames)* ...
        (sum(reward_trace,2)/num_frames)');
    reward_tk_ci = (coincident_events - chance_events)./mean_events;
    
    
    %%%% Save
    grab_var.move_start_max_xcorr_norm{curr_day} = move_start_max_xcorr_norm;
    grab_var.move_stop_max_xcorr_norm{curr_day} = move_stop_max_xcorr_norm;
    grab_var.cue_max_xcorr_norm{curr_day} = cue_max_xcorr_norm;
    grab_var.reward_max_xcorr_norm{curr_day} = reward_max_xcorr_norm;
    
    grab_var.move_start_max_xcorr_lag{curr_day} = move_start_max_xcorr_lag;
    grab_var.move_stop_max_xcorr_lag{curr_day} = move_stop_max_xcorr_lag;
    grab_var.cue_max_xcorr_lag{curr_day} = cue_max_xcorr_lag;
    grab_var.reward_max_xcorr_lag{curr_day} = reward_max_xcorr_lag;
    
    grab_var.move_start_tk_ci{curr_day} = move_start_tk_ci;
    grab_var.move_stop_tk_ci{curr_day} = move_stop_tk_ci;
    grab_var.cue_tk_ci{curr_day} = cue_tk_ci;
    grab_var.reward_tk_ci{curr_day} = reward_tk_ci;
    
    grab_var.move_tk_ci{curr_day} = move_tk_ci;
    grab_var.nomove_tk_ci{curr_day} = nomove_tk_ci;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);

% if the lags are negative, set the xcorr to zero
cue_max_xcorr_tk_ci_concat = cell2mat(grab_var.cue_tk_ci);
cue_max_xcorr_lag_concat = cell2mat(grab_var.cue_max_xcorr_lag);
cue_max_xcorr_tk_ci_concat(cue_max_xcorr_lag_concat < 0) = 0;

reward_max_xcorr_tk_ci_concat = cell2mat(grab_var.reward_tk_ci);
reward_max_xcorr_lag_concat = cell2mat(grab_var.reward_max_xcorr_lag);
reward_max_xcorr_tk_ci_concat(reward_max_xcorr_lag_concat < 0) = 0;

%% Behavior-triggered trace avg
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    %%%% Process
    %im.oopsi(im.oopsi < 0.2) = 0;
    
    % find rewarded trials, use only those cues and reward times
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
    good_trials = find(rewarded_trials & bhv.imaged_trials);
    reward_frames_good = cellfun(@(x) x.states.reward(1), ...
        bhv.bhv_frames(good_trials));
    cue_frames_good = cellfun(@(x) x.states.cue(1), ...
        bhv.bhv_frames(good_trials));
    % compensate for loops
    reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
    cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);
    
    % define velocity threshold for active movement
    movethresh = 0.02;
    
    % get velocity envelope, smooth over 5 frames (~150ms)
    lever_velocity = [0;diff(bhv.lever_force_resample)];
    lever_velocity_hilbert = hilbert(lever_velocity);
    lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
        conj(lever_velocity_hilbert));
    lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,5);
    
    % define active parts of the lever, fill in gaps of < 1/3 second
    lever_active = lever_velocity_envelope_smooth > movethresh;
    lever_active_switch = [0;diff(lever_active)];
    lever_active_edges = find(lever_active_switch);
    lever_active_spaces = [0;diff(lever_active_switch(lever_active_edges))];
    lever_active_up = find(lever_active_spaces == 2);
    lever_active_down = lever_active_up-1;
    lever_active_space_dist = lever_active_edges(lever_active_up) - ...
        lever_active_edges(lever_active_down);
    lever_active_fill = lever_active_space_dist < bhv.framerate/3;
    for i = find(lever_active_fill)'
        lever_active(lever_active_edges(lever_active_down(i)): ...
            lever_active_edges(lever_active_up(i))) = 1;
    end
    
    lever_stopstart = [0;diff(lever_active)];
    move_start_frames = find(lever_stopstart == 1);
    move_stop_frames = find(lever_stopstart == -1);
    
    % get rid of unpaired movement bouts at start and end
    if move_stop_frames(1) < move_start_frames(1)
        move_stop_frames(1) = [];
    end
    
    if move_start_frames(end) > move_stop_frames(end)
        move_start_frames(end) = [];
    end
    % get iti movement frames (movement cannot result in reward)
    % parse movements
    movement_frames = cell(size(move_start_frames),1);
    movements_p = cell(size(move_start_frames),1);
    movements_v = cell(size(move_start_frames),1);
    for i = 1:length(movement_frames)
        movement_frames{i} = move_start_frames(i):move_stop_frames(i);
        movements_p{i} = bhv.lever_force_resample(move_start_frames(i):...
            move_stop_frames(i));
        movements_v{i} = lever_velocity(move_start_frames(i):...
            move_stop_frames(i));
    end
    reward_overlap_movements = cellfun(@(x) any(ismember(x, ...
        round(reward_frames_good))), movement_frames);
    % give 5 frame leeway (sometimes close = overlaps by accident)
    cue_overlap_movements = cellfun(@(x) any(ismember(x, ...
        round(cue_frames_good-5))), movement_frames);
    
    iti_movements = ~reward_overlap_movements & ~cue_overlap_movements;
    cued_rewarded_movements = reward_overlap_movements & ~cue_overlap_movements;
    uncued_rewarded_movements = reward_overlap_movements & cue_overlap_movements;
    
    % get indicies for which cues are overlapped by movement
    cue_overlap_idx = cellfun(@(x) find(ismember(round(cue_frames_good),x)), ...
        movement_frames,'UniformOutput',false);
    
    % get movement properties
    move_time = move_stop_frames - move_start_frames;
    move_max_p = cellfun(@max,movements_p);
    move_max_v = cellfun(@max,movements_v);
    move_min_p = cellfun(@min,movements_p);
    move_min_v = cellfun(@min,movements_v);
    % pad movements with zeros, xcorr all
    movements_p_pad = cellfun(@(x) padarray(x, ...
        [max(move_time)+1-length(x) 0],'post'), ...
        movements_p,'UniformOutput',false);
    %max_xcorr = reshape(max(xcorr([movements_p_pad{:}],'coeff')),length(move_time),...
    %    length(move_time));
    
    % find sum deconv spike for each movement/cell
    movement_sum_oopsi = zeros(size(roi_trace_df,1),length(movement_frames));
    for i = 1:length(movement_frames)
    movement_sum_oopsi(:,i) = sum(im.oopsi(:,movement_frames{i}),2);
    end
    
    for spike_loop = 1:2
        
        if spike_loop == 1
            spike_frames = round(cue_frames_good);
        elseif spike_loop == 2
            spike_frames = move_start_frames;
        end
        
        lever_spike_all{spike_loop} = cell(size(roi_trace_df,1),1);
        lever_diff_spike_all{spike_loop} = cell(size(roi_trace_df,1),1);
        trace_spike_all{spike_loop} = cell(size(roi_trace_df,1),1);
        oopsi_spike_all{spike_loop} = cell(size(roi_trace_df,1),1);
        
        time_back = 2000; %ms
        time_forward = 2000; %ms
        
        xsg_sample_rate = 10000;
        
        samples_back = round((time_back/1000)*(xsg_sample_rate));
        samples_forward = round((time_forward/1000)*(xsg_sample_rate));
        
        frames_back = round((time_back/1000)*bhv.framerate);
        frames_forward = round((time_forward/1000)*bhv.framerate);
        
        %lever_force_smooth = smooth(bhv.lever_force,500);
        %lever_force_diff = [0;diff(lever_force_smooth)];
        
        % get frames of interest
        spike_frames_idx = repmat(spike_frames,1,frames_back+frames_forward+1);
        spike_frames_add = repmat([-frames_back:frames_forward],length(spike_frames),1);
        spike_frames_idx = spike_frames_idx + spike_frames_add;
        % knock out alignments where frames are out of bounds
        spike_frames_idx(any(spike_frames_idx <= 0,2) | ...
            any(spike_frames_idx > size(roi_trace_df,2),2),:) = [];
        
        for curr_cell = 1:size(roi_trace_df,1);
            clear lever_spike trace_spike oopsi_spike
            %lever_diff_spike = nan(length(spike_frames),samples_back+samples_forward+1);
            lever_spike = bhv.lever_force_resample(spike_frames_idx);
            lever_spike = reshape(lever_spike,size(spike_frames_idx));
            
            trace_spike = roi_trace_df(curr_cell,spike_frames_idx);
            trace_spike = reshape(trace_spike,size(spike_frames_idx));
            
            oopsi_spike = im.oopsi(curr_cell,spike_frames_idx);
            oopsi_spike = reshape(oopsi_spike,size(spike_frames_idx));
            
            lever_spike_all{spike_loop}{curr_cell} = lever_spike;
            %lever_diff_spike_all{spike_loop}{curr_cell} = lever_diff_spike;
            trace_spike_all{spike_loop}{curr_cell} = trace_spike;
            oopsi_spike_all{spike_loop}{curr_cell} = oopsi_spike;
        end
    end
    
    % for all movement parameters, find difference between cued/rewarded, iti
    move_params_act = {};
    for curr_param = 1:5;
        switch curr_param
            case 1
                move_params = move_max_p;
            case 2
                move_params = move_max_v;
            case 3
                move_params = move_min_p;
            case 4
                move_params = move_min_v;
            case 5
                move_params = move_time;
        end
        
        % bin movement parameters from highest 25% to lowest 75%
        cr_prctile = prctile(move_params(cued_rewarded_movements),[25 75]);
        iti_prctile = prctile(move_params(iti_movements),[25 75]);
        bin_lims = max([cr_prctile; iti_prctile]);
        bin_edges = bin_lims(1):diff(bin_lims)/5:bin_lims(2);
        
        [n,move_params_cr_bin] = histc(move_params(cued_rewarded_movements),bin_edges);
        [n,move_params_iti_bin] = histc(move_params(iti_movements),bin_edges);
        
        % save as cells 1: edges, 2: cr activity, 3: cr groups
        % 4: iti activity, 5: iti groups
        move_params_act{curr_param}{1} = bin_edges;
        move_params_act{curr_param}{2} = ...
            movement_sum_oopsi(:,cued_rewarded_movements);
        move_params_act{curr_param}{3} = move_params_cr_bin;
        move_params_act{curr_param}{4} = ...
            movement_sum_oopsi(:,iti_movements);
        move_params_act{curr_param}{5} = move_params_iti_bin;
        
        % this is for summary statistics, box plots are probably better
        %     cr_bin_mean = grpstats(movement_max_oopsi(curr_cell, ...
        %         cued_rewarded_movements),move_params_cued_bin,'mean');
        %     iti_bin_mean = grpstats(movement_max_oopsi(curr_cell, ...
        %         iti_movements),move_params_iti_bin,'mean');
        %
        %     cr_bin_sem = grpstats(movement_max_oopsi(curr_cell, ...
        %         cued_rewarded_movements),move_params_cued_bin,'sem');
        %     iti_bin_sem = grpstats(movement_max_oopsi(curr_cell, ...
        %         iti_movements),move_params_iti_bin,'sem');
    end
    %%%% Save
    grab_var.trace_spike_all{curr_day} = trace_spike_all;
    grab_var.oopsi_spike_all{curr_day} = oopsi_spike_all;
    grab_var.move_params_act{curr_day} = move_params_act;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);

% plot a cell over days
curr_cell = 7;
figure('Name',num2str(curr_cell));
for i = 1:14
   subplot(2,7,i);
   hold on;
   
   cr_bin_mean = grpstats(grab_var.move_params_act{i}{1}{2}(curr_cell,:), ...
       grab_var.move_params_act{i}{1}{3},'mean');
   iti_bin_mean = grpstats(grab_var.move_params_act{i}{1}{4}(curr_cell,:), ...
       grab_var.move_params_act{i}{1}{5},'mean');
   
   cr_bin_std = grpstats(grab_var.move_params_act{i}{1}{2}(curr_cell,:), ...
       grab_var.move_params_act{i}{1}{3},'std');
   iti_bin_std = grpstats(grab_var.move_params_act{i}{1}{4}(curr_cell,:), ...
       grab_var.move_params_act{i}{1}{5},'std');
   
   errorbar([1:length(cr_bin_mean)],cr_bin_mean, ...
       cr_bin_std,'r')
   errorbar([1:length(iti_bin_mean)],iti_bin_mean,...
       iti_bin_std,'k')

%    boxplot(grab_var.move_params_act{i}{1}{2}(curr_cell,:), ...
%        grab_var.move_params_act{i}{1}{3},'color','r','notch','on');
%    boxplot(grab_var.move_params_act{i}{1}{4}(curr_cell,:), ...
%        grab_var.move_params_act{i}{1}{5},'color','k','notch','on');
   %axis off
   title(i)
end

% plot difference over days
curr_cell = [7];
act_diff = zeros(14,1);
act_diff_sem = zeros(14,1);
figure;title(curr_cell);
for j = 1:5
    subplot(1,5,j);
    for i = 1:14
        cr_binned = grab_var.move_params_act{i}{j}{3} > 0;
        iti_binned = grab_var.move_params_act{i}{j}{5} > 0;
        
        act_diff(i) = mean(grab_var.move_params_act{i}{j}{2}(curr_cell,cr_binned)) - ...
            mean(grab_var.move_params_act{i}{j}{4}(curr_cell,iti_binned));
        act_diff_sem(i) = sqrt(...
            (std(grab_var.move_params_act{i}{j}{2}(curr_cell,cr_binned))^2/ ...
            sum(cr_binned)) + ...
            (std(grab_var.move_params_act{i}{j}{4}(curr_cell,iti_binned))^2/ ...
            sum(iti_binned)));
    end
    errorbar(act_diff,2*act_diff_sem)
    line(xlim,[0 0],'color','k');
end


%% Extract raw rewarded lever presses
clear all

animals = {'AP66' 'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
for curr_animal = 1:length(animals)

animal = animals{curr_animal};

% find days
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
    
grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var data_path animals curr_animal
    
    day = num2str(days{curr_day});
    
    % get behavior file
    dir_currfolder = dir([data_path filesep day]);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
    
    % get xsg files
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    
    % get behavior data and xsg sample rate
    try
        % ignore if something's wrong with datafile (usually, >1 of them)
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
    catch me
        continue
    end
    bhv = saved_history.ProtocolsSection_parsed_events;
    num_trials = length(bhv);
    
    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    %%%% Process
    
    % Only look at rewarded trials
    rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));
   
    time_surround = 1; % in seconds
    sample_surround = time_surround*xsg_samplerate;
    rewarded_lever_trace = nan(num_trials,sample_surround*2+1);
    rewarded_lever_trace_vel = nan(num_trials,sample_surround*2+1);
    
    % go through all xsg files, grab rewarded lever presses
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        lever_vel = [0;diff(smooth(xsg_data.data.acquirer.trace_2,500))];
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        for curr_trial = curr_trial_list(:,2)';
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % get trial/xsg sample offset
            curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
            % get reward sample in xsg
            curr_bhv_reward_sample_rel =  (bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*xsg_samplerate;
            curr_bhv_reward_sample_abs = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *xsg_samplerate + ...
                curr_bhv_reward_sample_rel);
            % ignore if out of bounds
            
            if curr_bhv_reward_sample_abs + sample_surround > length( ...
                    xsg_data.data.acquirer.trace_2) || ...
                    curr_bhv_reward_sample_abs - sample_surround < 0;
                continue
            end
            
            rewarded_lever_trace(curr_trial,:) = ...
                xsg_data.data.acquirer.trace_2( ...
                curr_bhv_reward_sample_abs - sample_surround : ...
                curr_bhv_reward_sample_abs + sample_surround);
            rewarded_lever_trace_vel(curr_trial,:) = ...
                lever_vel( ...
                curr_bhv_reward_sample_abs - sample_surround : ...
                curr_bhv_reward_sample_abs + sample_surround);           
        end
    end
           
    %%%% Save
    grab_var.rewarded_lever_trace{curr_day} =  rewarded_lever_trace;
    grab_var.rewarded_lever_trace_vel{curr_day} = rewarded_lever_trace_vel;
    
    disp(['Finished day: ' day]);
    
end
save(['/usr/local/lab/People/Andy/Data/lever_data/' animal '_rewarded_lever'], 'grab_var');
disp(['Finished all.']);
end

% mean +/- std
animals = {'AP61' 'AP62' 'AP63' 'AP64' 'AP65' 'AP66' 'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([animal '_rewarded_lever.mat'])
    figure;
    days = length(grab_var.rewarded_lever_trace);
    sqrtdays = ceil(sqrt(days));
    for i = 1:days
        subplot(sqrtdays,sqrtdays,i)
        hold on
        plot(nanmean(grab_var.rewarded_lever_trace{i}))
        plot(nanmean(grab_var.rewarded_lever_trace{i}) + ...
            nanstd(grab_var.rewarded_lever_trace{i}),'r')
        plot(nanmean(grab_var.rewarded_lever_trace{i}) - ...
            nanstd(grab_var.rewarded_lever_trace{i}),'r')
        %ylim([-0.0008 0.001])
        ylim([0.5 2.5])
        xlim([0 20001])
    end   
    saveas(gcf,[animal '_rewarded_lever.png'],'png');
    close all
end

% heatmap
animals = {'AP61' 'AP62' 'AP63' 'AP64' 'AP65' 'AP66' 'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([animal '_rewarded_lever.mat'])
    figure;
    days = length(grab_var.rewarded_lever_trace);
    sqrtdays = ceil(sqrt(days));
    for i = 1:days
        subplot(sqrtdays,sqrtdays,i)
        imagesc(grab_var.rewarded_lever_trace{i})
        caxis([0.5 2.5])
        axis off
    end   
    saveas(gcf,[animal '_rewarded_lever_heat.png'],'png');
    close all
end

% corrcoef
animals = {'AP61' 'AP62' 'AP63' 'AP64' 'AP65' 'AP66' 'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([animal '_rewarded_lever.mat'])
    comp_press = nanmean(grab_var.rewarded_lever_trace{end}(:,7000:15000))';
    days = length(grab_var.rewarded_lever_trace);
    sqrtdays = ceil(sqrt(days));
    corr_days = cell(length(days),1);
    for i = 1:days
        if isempty(grab_var.rewarded_lever_trace{i})
            continue
        end
        temp_corr = corrcoef([ ...
            grab_var.rewarded_lever_trace{i}(:,7000:15000)' ...
            comp_press]);
        corr_days{i} = temp_corr(end,:);
    end   
    figure;
    subplot(3,1,1)
    plot(cell2mat(corr_days),'.k')
    corr_sum = cumsum(cellfun(@length,corr_days))+0.5;
    for i = 1:length(corr_days)
        line([corr_sum(i) corr_sum(i)],ylim,'color','r')
    end
    xlim([-50 corr_sum(end)+50])
    subplot(3,1,2);
    errorbar(cellfun(@nanmean,corr_days),cellfun(@nanstd,corr_days),'k');
    xlim([0 days+1])
    subplot(3,1,3);
    plot(comp_press,'k')
    saveas(gcf,[animal '_rewarded_lever_corr_300end.png'],'png');
    close all
end

% movement properties
animals = {'AP61' 'AP62' 'AP63' 'AP64' 'AP65' 'AP66' 'AP67' 'AP68' 'AP69' 'AP70' 'AP71' 'AP72'};
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([animal '_rewarded_lever.mat'])
    days = length(grab_var.rewarded_lever_trace);
    
    move_max_p = cellfun(@(x) max(x,[],2),grab_var.rewarded_lever_trace,'UniformOutput',false);
    move_max_v = cellfun(@(x) max(x,[],2),grab_var.rewarded_lever_trace_vel,'UniformOutput',false);
    move_min_p = cellfun(@(x) min(x,[],2),grab_var.rewarded_lever_trace,'UniformOutput',false);
    move_min_v = cellfun(@(x) min(x,[],2),grab_var.rewarded_lever_trace_vel,'UniformOutput',false);
    
    % get movement times
    movement_times = cell(days,1);
    movethresh = 0.0002;
    
    for curr_day = 1:days
        lever_velocity = grab_var.rewarded_lever_trace_vel{curr_day}';
        lever_velocity = lever_velocity(:);

        for press = 1:size(grab_var.rewarded_lever_trace{curr_day},1);
            lever_velocity = grab_var.rewarded_lever_trace_vel{curr_day}(press,:);
            lever_velocity_resample = resample(lever_velocity,1,10);
            lever_velocity_resample_smooth = smooth(lever_velocity_resample,5);
            lever_velocity_resample_smooth(isnan(lever_velocity_resample)) = NaN;
            
            % get velocity envelope, smooth over 5 frames (~150ms)
            lever_velocity_hilbert = hilbert(lever_velocity_resample_smooth);
            lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
                conj(lever_velocity_hilbert));
            lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,20);
            
            % define active parts of the lever, fill in gaps of < 1 second (30 frames)
            lever_active = lever_velocity_envelope_smooth > movethresh;
            lever_active_switch = [0;diff(lever_active)];
            lever_active_edges = find(lever_active_switch);
            lever_active_spaces = [0;diff(lever_active_switch(lever_active_edges))];
            lever_active_up = find(lever_active_spaces == 2);
            lever_active_down = lever_active_up-1;
            lever_active_space_dist = lever_active_edges(lever_active_up) - ...
                lever_active_edges(lever_active_down);
            lever_active_fill = lever_active_space_dist < 10;
            for i = find(lever_active_fill)'
                lever_active(lever_active_edges(lever_active_down(i)): ...
                    lever_active_edges(lever_active_up(i))) = 1;
            end
            
            lever_stopstart = [0;diff(lever_active)];
           
            % get the first movement stop after push
            midpoint = ceil(length(lever_stopstart)/2);
            move_stop_sample = find(lever_stopstart(midpoint:end) == -1,1);
            move_start_sample = find(lever_stopstart(midpoint:-1:1) == 1,1);
            if ~isempty(move_stop_sample) && ~isempty(move_start_sample)
                movement_times{curr_day}(press,1) = move_start_sample+move_stop_sample;
            else
                movement_times{curr_day}(press,1) = NaN;
            end
        end
    end
    
    figure;
    for plots = 1:5
        subplot(2,3,plots)
        switch plots
            case 1
                plotprop = move_max_p;
                name = 'move max p';
            case 2
                plotprop = move_min_p;
                name = 'move min p';
            case 3
                plotprop = move_max_v;
                name = 'move max v';
            case 4
                plotprop = move_min_v;
                name = 'move min v';
            case 5
                plotprop = movement_times;
                name = 'movement times';
        end
        plot(vertcat(plotprop{:}),'.k')
        corr_sum = cumsum(cellfun(@length,plotprop))+0.5;
        for i = 1:length(corr_sum)
            line([corr_sum(i) corr_sum(i)],ylim,'color','r')
        end
        title(name);
    end
    
    saveas(gcf,[animal '_rewarded_lever_params.png'],'png');
    close all
    
end

   


%% Get mean trial max trace for cued-rewarded presses
clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    %%%% Process
    loop_size = cellfun(@length,roi_trace_long_split);
    loop_size_sum = cumsum(loop_size);
    loop_frames = [0;loop_size_sum(8:8:end-7)];
    
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    
    % get behavior file
    dir_currfolder = dir([data_path filesep day]);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
    
    % get xsg files
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    
    % get behavior data and xsg sample rate
    try
        % ignore if something's wrong with datafile (usually, >1 of them)
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
    catch me
        error('Something wrong with data file')
    end
    raw_bhv = saved_history.ProtocolsSection_parsed_events;
    num_trials = length(raw_bhv);
    
    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    % Only look at rewarded trials
    rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), raw_bhv));
    cued_movement = nan(num_trials,1);
    cue_frames = nan(num_trials,1);
    reward_frames = nan(num_trials,1);
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        
        % create movement velocity, active times
        % quantify movement before
        
        %lever_force_smooth = smooth(xsg_data.data.acquirer.trace_2,500);
        lever_force_resample = resample(xsg_data.data.acquirer.trace_2,1,10);
        % normalized cutoff freq is fraction of nyquist
        butterworth_freq = 20/500;
        [b a] = butter(4, butterworth_freq);
        % using filtfilt instead of filter allows for zero-phase shift
        lever_force_smooth = filtfilt(b,a,lever_force_resample);
        
        lever_velocity_resample = [0;diff(lever_force_smooth)];
        %lever_velocity_resample = resample(lever_velocity,1,10);
        %lever_force_smooth_resample = resample(lever_force_smooth,1,10);
        lever_velocity_resample_smooth = smooth(lever_velocity_resample,5);
        lever_velocity_resample_smooth(isnan(lever_velocity_resample)) = NaN;
        
        % get velocity envelope, smooth over 5 frames (~150ms)
        lever_velocity_hilbert = hilbert(lever_velocity_resample_smooth);
        lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
            conj(lever_velocity_hilbert));
        % at the moment - don't smooth
        lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,1);
        
        % define active parts of the lever, ignore small ones and
        % fill in gaps of < 100ms
        movethresh = 0.001;
        lever_active = lever_velocity_envelope_smooth > movethresh;
        lever_v_thresh = lever_velocity_envelope_smooth;
        lever_v_thresh(lever_v_thresh < 0.002) = 0;
        % get rid of short movements, find movement times again
        for short_move_erase = 1:3
            lever_active_switch = [0;diff(lever_active)];
            lever_active_starts = find(lever_active_switch == 1);
            lever_active_stops = find(lever_active_switch == -1);
            % get rid of unpaired starts/stops
            if lever_active_stops(1) < lever_active_starts(1);
                lever_active_stops(1) = [];
            end
            if lever_active_starts(end) > lever_active_stops(end)
                lever_active_starts(end) = [];
            end
            
            lever_active_movement_times = lever_active_stops - ...
                lever_active_starts;
            lever_active_intermovement_times = lever_active_starts(2:end) - ...
                lever_active_stops(1:end-1);
            
            if short_move_erase == 1
                % fill in short gaps between movements, but only if
                % there's a movement afterwards that's a certain time
                lever_active_fill = ...
                    lever_active_intermovement_times < 20 & ...
                    lever_active_movement_times(2:end) > 50;
                for i = find(lever_active_fill)'
                    lever_active(lever_active_stops(i): ...
                        lever_active_starts(i+1)) = 1;
                end
            else
                % erase movement that's too brief
                lever_active_erase = lever_active_movement_times < 100;
                for i = find(lever_active_erase')
                    lever_active(lever_active_starts(i): ...
                        lever_active_stops(i)) = 0;
                end
            end
        end
        
        % fill in short gaps between movements
        lever_active_intermovement_times = lever_active_starts(2:end) - ...
            lever_active_stops(1:end-1);
        lever_active_fill = lever_active_intermovement_times < 200;
        for i = find(lever_active_fill)'
            lever_active(lever_active_stops(i): ...
                lever_active_starts(i+1)) = 1;
        end
        
        lever_stopstart = [0;diff(lever_active)];
        lever_stopstart = [0;diff(lever_active)];
        
        % get cue and reward times in resampled xsg
        cue_sample = nan(length(curr_trial_list(:,2)),1);
        reward_start_sample = nan(length(curr_trial_list(:,2)),1);
        reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
        for curr_trial_idx = 1:length(curr_trial_list(:,2));
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % get trial start in dispatcher
            curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
            % get cue and reward sample in xsg (downsampled to ms)
            curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            cue_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_cue_sample_rel);
            curr_bhv_cue_frame_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(bhv.framerate);
            cue_frames(curr_trial) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(bhv.framerate) + ...
                curr_bhv_cue_frame_rel) + loop_frames(curr_xsg);
            
            curr_bhv_reward_start_sample_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_start_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_start_sample_rel);
            curr_bhv_reward_start_frame_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(bhv.framerate);
            reward_frames(curr_trial) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(bhv.framerate) + ...
                curr_bhv_reward_start_frame_rel) + loop_frames(curr_xsg);
            
            curr_bhv_reward_stop_sample_rel =  (raw_bhv{curr_trial}.states.reward(2) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_stop_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_stop_sample_rel);
        end
        
        
        for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
            
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials) || ...
                    length(lever_velocity_envelope_smooth) < cue_sample(curr_trial_idx+1)
                continue
            end
            
            % ignore if out of bounds
            if cue_sample(curr_trial_idx+1) > length( ...
                    lever_velocity_envelope_smooth)
                continue
            end
            
            % classify trial cue/movement overlap
            % rewarded-movement/cue overlap can't be > 50ms
            rewarded_move_start = find(lever_active(1: ...
                reward_start_sample(curr_trial_idx)) == 0,1,'last');
            cued_movement(curr_trial) = rewarded_move_start > cue_sample(curr_trial_idx)-50;
        end
    end
    
    roi_trace_df_smooth = zeros(size(roi_trace_df));
    for curr_cell = 1:size(roi_trace_df,1)
        
        % normalized cutoff freq is fraction of nyquist
        butterworth_freq = (0.1*bhv.framerate)/(bhv.framerate/2);
        [b a] = butter(4, butterworth_freq);
        % using filtfilt instead of filter allows for zero-phase shift
        roi_trace_df_smooth(curr_cell,:) = filtfilt(b,a,roi_trace_df(curr_cell,:));
    end
    
    % Get the maximum response of cell during each trial
    cell_trial_max = nan(size(roi_trace_df,1),num_trials);
    for i = 1:num_trials-1
        if ~any(isnan(cue_frames(i:i+1)));
            cell_trial_max(:,i) = max(roi_trace_df_smooth(:, ...
                cue_frames(i):cue_frames(i+1)),[],2);
        end
    end
    
    % get sum deconvoluted spikes
    cell_trial_oopsi = nan(size(roi_trace_df,1),num_trials);
    oopsi_thresh = im.oopsi;
    oopsi_thresh(oopsi_thresh < 0.2) = 0;
    for i = 1:num_trials-1
        if ~any(isnan(cue_frames(i:i+1)));
            cell_trial_oopsi(:,i) = sum(im.oopsi(:, ...
                cue_frames(i):cue_frames(i+1)),2);
        end
    end
    
    
    % pull out the traces for cued-rewarded trials
    trial_traces = cell(num_trials,1);
    for i = 1:num_trials-1;
        if ~isnan(cue_frames(i)) && ~isnan(cue_frames(i+1))
            trial_traces{i} = roi_trace_df(:,cue_frames(i):cue_frames(i+1));
        end
    end



    %%%% Save
    grab_var.cell_trial_max{curr_day} = cell_trial_max;
    grab_var.cell_trial_oopsi{curr_day} = cell_trial_oopsi;
    grab_var.trial_traces{curr_day} = trial_traces;
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);




%% Get cells significantly modulated by (cued movement/ all movement)
clear all

animal = 'AP74';

% find days
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

grab_var = {};

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    % temporary, while deconv hasn't caught up
    if ~isfield(im,'roi_trace_long_split')
        im.roi_trace_long_split = roi_trace_long_split;
        im.roi_trace_df = roi_trace_df;
    end
    
    %%%% Process
    
    move_trial_mean = [];
    still_trial_mean = [];
    
    % Break up baseline corrected trace into loops
    num_loops = length(im.roi_trace_long_split)/8;
    file_lengths = cellfun(@length,im.roi_trace_long_split);
    file_lengths_loops = reshape(file_lengths,8,num_loops);
    loop_lengths = sum(file_lengths_loops);
    oopsi_split = mat2cell(im.oopsi,size(im.oopsi,1), ...
        loop_lengths);
    
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    
    % get behavior file
    dir_currfolder = dir([data_path filesep day]);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
    
    % get xsg files
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    
    % get behavior data and xsg sample rate
    try
        % ignore if something's wrong with datafile (usually, >1 of them)
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
    catch me
        error('Problem with behavior file')
    end
    raw_bhv = saved_history.ProtocolsSection_parsed_events;
    num_trials = length(raw_bhv);
    
    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames));
    
    lever_active_frames = {};
    movement_epochs = {};
    % go through all xsg files, find movement times
    plots_on = 0;
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        
        % create movement velocity, active times
        % quantify movement before
        
        %lever_force_smooth = smooth(xsg_data.data.acquirer.trace_2,500);
        lever_force_resample = resample(xsg_data.data.acquirer.trace_2,1,10);
        % normalized cutoff freq is fraction of nyquist
        butterworth_freq = 20/500;
        [b a] = butter(4, butterworth_freq);
        % using filtfilt instead of filter allows for zero-phase shift
        lever_force_smooth = filtfilt(b,a,lever_force_resample);
        
        lever_velocity_resample = [0;diff(lever_force_smooth)];
        %lever_velocity_resample = resample(lever_velocity,1,10);
        %lever_force_smooth_resample = resample(lever_force_smooth,1,10);
        lever_velocity_resample_smooth = smooth(lever_velocity_resample,5);
        lever_velocity_resample_smooth(isnan(lever_velocity_resample)) = NaN;
        
        % get velocity envelope, smooth over 5 frames (~150ms)
        lever_velocity_hilbert = hilbert(lever_velocity_resample_smooth);
        lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
            conj(lever_velocity_hilbert));
        % at the moment - don't smooth
        lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,1);
        
        % define active parts of the lever, ignore small ones and
        % fill in gaps of < 100ms
        movethresh = 0.001;
        lever_active = lever_velocity_envelope_smooth > movethresh;
        lever_v_thresh = lever_velocity_envelope_smooth;
        lever_v_thresh(lever_v_thresh < 0.002) = 0;
        % get rid of short movements, find movement times again
        for short_move_erase = 1:3
            lever_active_switch = [0;diff(lever_active)];
            lever_active_starts = find(lever_active_switch == 1);
            lever_active_stops = find(lever_active_switch == -1);
            % get rid of unpaired starts/stops
            if lever_active_stops(1) < lever_active_starts(1);
                lever_active_stops(1) = [];
            end
            if lever_active_starts(end) > lever_active_stops(end)
                lever_active_starts(end) = [];
            end
            
            lever_active_movement_times = lever_active_stops - ...
                lever_active_starts;
            lever_active_intermovement_times = lever_active_starts(2:end) - ...
                lever_active_stops(1:end-1);
            
            if short_move_erase == 1
                % fill in short gaps between movements, but only if
                % there's a movement afterwards that's a certain time
                lever_active_fill = ...
                    lever_active_intermovement_times < 20 & ...
                    lever_active_movement_times(2:end) > 50;
                for i = find(lever_active_fill)'
                    lever_active(lever_active_stops(i): ...
                        lever_active_starts(i+1)) = 1;
                end
            else
                % erase movement that's too brief
                lever_active_erase = lever_active_movement_times < 100;
                for i = find(lever_active_erase')
                    lever_active(lever_active_starts(i): ...
                        lever_active_stops(i)) = 0;
                end
            end
        end
        
        % fill in short gaps between movements
        lever_active_intermovement_times = lever_active_starts(2:end) - ...
            lever_active_stops(1:end-1);
        lever_active_fill = lever_active_intermovement_times < 200;
        for i = find(lever_active_fill)'
            lever_active(lever_active_stops(i): ...
                lever_active_starts(i+1)) = 1;
        end
        
        lever_stopstart = [0;diff(lever_active)];
        
        move_frames = round(bhv.framerate*...
            (find(lever_active)./(xsg_samplerate/10)));
        still_frames = round(bhv.framerate*...
            (find(~lever_active)./(xsg_samplerate/10)));
        
        % split up the movement into movement epochs (for shuffling)
        lever_active_frames{curr_xsg} = ...
            zeros(1,size(oopsi_split{curr_xsg},2));
        lever_active_frames{curr_xsg}(intersect(1:length( ...
            lever_active_frames{curr_xsg}),move_frames)) = 1;
        % find size of movement epochs
        movement_epoch_starts = find(diff([0 lever_active_frames{curr_xsg} 0]) == 1);
        movement_epoch_ends = find(diff([0 lever_active_frames{curr_xsg} 0]) == -1);
        movement_epoch_sizes = movement_epoch_ends - movement_epoch_starts;
        % create cell array with grouped movement sizes
        movement_epoch_cell = false(size(lever_active_frames{curr_xsg}));
        movement_epoch_cell(1:sum(movement_epoch_sizes)) = true;
        movement_epoch_cell = mat2cell(movement_epoch_cell,...
            1, [movement_epoch_sizes ones(1,sum(~lever_active_frames{curr_xsg}))]);
        movement_epochs{curr_xsg} = movement_epoch_cell;
                
        % get cue/reward times in current loop
        cue_sample = nan(length(curr_trial_list(:,2)),1);
        reward_start_sample = nan(length(curr_trial_list(:,2)),1);
        reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
        for curr_trial_idx = 1:length(curr_trial_list(:,2));
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % get trial start in dispatcher
            curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
            % get cue and reward sample in xsg (downsampled to ms)
            curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            cue_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_cue_sample_rel);
            
            curr_bhv_reward_start_sample_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_start_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_start_sample_rel);
            
            curr_bhv_reward_stop_sample_rel =  (raw_bhv{curr_trial}.states.reward(2) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_stop_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_stop_sample_rel);
        end
  
        % find modulated cells based on responses to cued-rewarded movements
        for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;      
            curr_trial = curr_trial_list(curr_trial_idx,2);
  
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % classify trial cue/movement overlap
            % rewarded-movement/cue overlap can't be > 50ms
            rewarded_move_start_curr = find(lever_active(1: ...
                reward_start_sample(curr_trial_idx)) == 0,1,'last');
            curr_cued_movement = rewarded_move_start_curr > cue_sample(curr_trial_idx)-50;
            
            % save the movement start time and classification
            cued_rewarded_move_start_frame(curr_trial) = ...
                round(rewarded_move_start_curr/(xsg_samplerate/10) ...
                *bhv.framerate+bhv.add_frames(curr_trial));
            cued_movement(curr_trial) = curr_cued_movement;
        end     
    end   
        
    %%%% Prepare the activity measure %%%%
    
    % set activity around loop stitches to zero
    activity_thresh = im.roi_trace_df;
    activity_thresh(activity_thresh < 0.5) = 0;
    
    % Break up baseline corrected trace into loops
    roi_trace_df_split = mat2cell(im.roi_trace_df,size(im.roi_trace_df,1), ...
        loop_lengths);
    split_indx = mat2cell(1:size(im.roi_trace_df,2),1, ...
        loop_lengths);
    
    inactive_thresh = 0.1;
    
    stitching_indx_all = cell(size(im.roi_trace_df,1),1);
    for curr_cell = 1:size(im.roi_trace_df,1)
        
        % Stitch together loops at inactive parts, get index of stitched frames
        loop_start_stitches = cellfun(@(x) find(x(curr_cell,:) ...
            < inactive_thresh & x(curr_cell,:) > 0,1), ...
            roi_trace_df_split,'UniformOutput',false);
        loop_stop_stitches = cellfun(@(x) find(x(curr_cell,:) ...
            < inactive_thresh & x(curr_cell,:) > 0,1,'last'), ...
            roi_trace_df_split,'UniformOutput',false);
        F_stitched = cellfun(@(x,y,z) x(curr_cell,y:z),roi_trace_df_split,...
            loop_start_stitches,loop_stop_stitches,'UniformOutput',false);
        stitching_indx = cellfun(@(x,y,z) x(y:z), split_indx, ...
            loop_start_stitches, loop_stop_stitches,'UniformOutput',false);
        
        stitching_indx_all{curr_cell} = [stitching_indx{:}];
        stitch_out = true(size(im.roi_trace_df,2),1);
        stitch_out(stitching_indx_all{curr_cell}) = false;
        activity_thresh(curr_cell,stitch_out) = 0;        
    end   
    
    num_rois = size(activity_thresh,1);
    
    %%%% Find significantly modulated cells via shuffling %%%%
    % this is to find modulated cells from cued-rewarded movements
    % time to look post cue for activity (in seconds)
    postcue_time = 1;
    postcue_frames = round(postcue_time*bhv.framerate);
    
    % get the real activity/behavior measure
    cue_cued_movement = num2cell(cellfun(@(x) round(x.states.cue(1)), ...
        bhv.bhv_frames(cued_movement))+bhv.add_frames(cued_movement));
    reward_cued_movement = num2cell(cellfun(@(x) round(x.states.reward(1)), ...
        bhv.bhv_frames(cued_movement))+bhv.add_frames(cued_movement));
    
    cued_movement_all = zeros(1,size(activity_thresh,2));
   
    cued_movement_start = num2cell(cued_rewarded_move_start_frame(cued_movement));
    
    cued_movement_frames = cellfun(@(x,y) x:y+postcue_frames,...
        cued_movement_start',reward_cued_movement,'UniformOutput',false);
    cued_movement_all([cued_movement_frames{:}]) = 1;
    activity_behavior_real = activity_thresh*cued_movement_all';

    % prepare and do shuffling
    num_cued_movements = sum(cued_movement);
    cued_movement_epochs = cellfun(@(x) ones(1,length(x)), ...
        cued_movement_frames,'UniformOutput',false);
    cued_movement_zeros = num2cell(zeros(1,size( ...
        activity_thresh,2) - sum(cellfun(@length, cued_movement_frames))));
    cued_movement_epochs = [cued_movement_epochs' cued_movement_zeros];
    
    num_rep = 1000;
    activity_behavior_shuffle = zeros(num_rois,num_rep);
    for i = 1:num_rep
        curr_perm = randperm(length(cued_movement_epochs));
        activity_behavior_shuffle(:,i) = activity_thresh*[cued_movement_epochs{curr_perm}]';
    end
    cued_movement_shuffle_prctile = prctile(activity_behavior_shuffle',[0.5 99.5]);
    move_cells = find(activity_behavior_real' > cued_movement_shuffle_prctile(2,:));
    still_cells = find(activity_behavior_real' < cued_movement_shuffle_prctile(1,:));
    
    %%%%%%%%%%%%%%%%%
    % this is thrown in here for now: get max amplitude of each cell for
    % each cued movement period to look at robustness
    cued_movement_activity_amplitude = cell2mat(cellfun(@(x) ...
        max(activity_thresh(:,x),[],2), cued_movement_frames,'UniformOutput',false)');
    %%%%%%%%%%%%%%%%%
    
    % this is to find modulated cells for all movement
%     num_rois = size(im.roi_trace_df,1);
%     movement_epochs_full = [movement_epochs{:}];
%     move_deconv = activity_thresh*[lever_active_frames{:}]';
%     
%     num_rep = 1000;
%     move_deconv_shuffle = zeros(num_rois,num_rep);
%     for i = 1:num_rep
%         curr_perm = randperm(length(movement_epochs_full));
%         move_deconv_shuffle(:,i) = activity_thresh*[movement_epochs_full{curr_perm}]';
%     end
%     move_deconv_prctile = prctile(move_deconv_shuffle',[0.5 99.5]);
%     move_cells = find(move_deconv' > move_deconv_prctile(2,:));
%     still_cells = find(move_deconv' < move_deconv_prctile(1,:));
    
    % get the std away from the mean for each cell's activity
    activity_stds = (activity_behavior_real - mean(activity_behavior_shuffle,2))./ ...
        std(activity_behavior_shuffle,[],2);
    
    
    %%%% Save
    grab_var.move_cells{curr_day} =  move_cells;
    grab_var.still_cells{curr_day} = still_cells;
    grab_var.activity_stds{curr_day} = activity_stds;
    grab_var.cued_movement_activity_amplitude{curr_day} = cued_movement_activity_amplitude;
    grab_var.cued_movement_frames{curr_day} = cued_movement_frames;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);


%% Align cellular activity by cued movement (requires last cell)

% animal is taken from before

% find days
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

aligned_activity = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var aligned_activity
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    %%%% Process
    curr_aligned_activity = cell(size(im.roi_trace_df,1),1);
    for curr_cell = 1:size(im.roi_trace_df,1);
        curr_aligned_activity{curr_cell} = cellfun(@(x) ...
            im.roi_trace_df(curr_cell,x), grab_var.cued_movement_frames{curr_day}, ...
            'UniformOutput',false);
    end
    
    %%%% Save
    aligned_activity{curr_day} =  curr_aligned_activity;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);

% plot a cell
curr_cell = 209;
figure;
for i = 1:14
    subplot(7,2,i)
    % pad activity until all traces are same length
    aligned_activity_curr = aligned_activity{i}{curr_cell};
    max_frames = max(cellfun(@length,aligned_activity_curr));
    aligned_activity_padded = cell2mat(cellfun(@(x) ...
        padarray(x,[0 max_frames-length(x)],'post'),aligned_activity_curr,...
        'UniformOutput',false));
    imagesc(aligned_activity_padded);
    caxis([0.5 2])
end

%% Get length correlations over days

animal = 'AP74';

% find days
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

length_corr = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var length_corr
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    % load roi file
    roi_name = [day '_' animal '_summed_50_template.roi'];
    load([analysis_path filesep roi_name],'-MAT')
    
    % load labels file
    labels_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep labels_name],'-MAT');
    
    % get cell labels
    cells = 1:size(im.roi_trace_df,1);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels));
    
    % exclude contaminated cells
    gad_cells(ismember(gad_cells,contaminated)) = [];
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    
    %%%% Process
    
    % Get centers/distances of all ROIs
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for i  = 1:num_rois
        [area roi_center_x(i) roi_center_y(i)] = ...
            polycenter(polygon.ROI{i}(:,1),polygon.ROI{i}(:,2));
    end
    
    % create matrix of distances
    roi_dist = nan(length(roi_center_x));
    for i = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(i) - roi_center_x(j))^2 + ...
                (roi_center_y(i) - roi_center_y(j))^2);
            roi_dist(i,j) = curr_dist;
        end
    end
    
    % get thresholded trace
    thresholded_trace = im.roi_trace_df;
    thresholded_trace(thresholded_trace < 0.5) = 0;
    
    % get rid of stitching artifacts: cut first 10 frames of all loops
    num_loops = length(im.roi_trace_long_split)/8;
    file_lengths = cellfun(@length,im.roi_trace_long_split);
    file_lengths_loops = reshape(file_lengths,8,num_loops);
    loop_lengths = cumsum(sum(file_lengths_loops));
    loop_starts = [1 loop_lengths(1:end-1) + 1];
    loop_start_cut = repmat(loop_starts,10,1) + repmat([0:9]',1,num_loops);
    thresholded_trace(:,loop_start_cut) = [];
    
    % get all correlations
    trace_correlation = corrcoef(thresholded_trace');
    
    % prepare bins
    binsize = [1 801];
    num_bins = length(binsize);
    
    mod_cells = grab_var.move_cells{curr_day};
    %mod_cells = grab_var.mod_cells;
    unmod_cells = setdiff(1:size(im.roi_trace_df,1),mod_cells);
    unmod_cells_active = unmod_cells(any(thresholded_trace(unmod_cells,:),2));
    
    % unmod pairs
    unmod_corr = trace_correlation(unmod_cells_active,unmod_cells_active);
    unmod_dist = roi_dist(unmod_cells_active,unmod_cells_active);
    
    % get rid of double pairs and self-comparisons
    unmod_corr = tril(unmod_corr,-1);
    unmod_dist = tril(unmod_dist,-1);
    
    [n,unmod_bin] = histc(unmod_dist(:),binsize);
    unmod_means = grpstats(unmod_corr(:),unmod_bin);
    unmod_sems = grpstats(unmod_corr(:),unmod_bin,'sem');
    
    % mod pairs
    mod_corr = trace_correlation(mod_cells,mod_cells);
    mod_dist = roi_dist(mod_cells,mod_cells);
    
    % get rid of double pairs and self-comparisons
    mod_corr = tril(mod_corr,-1);
    mod_dist = tril(mod_dist,-1);
    
    [n,mod_bin] = histc(mod_dist(:),binsize);
    mod_means = grpstats(mod_corr(:),mod_bin);
    mod_sems = grpstats(mod_corr(:),mod_bin,'sem');
    
    
    %%%% Save
    length_corr.unmod_means{curr_day} =  unmod_means;
    length_corr.unmod_sems{curr_day} =  unmod_sems;
    length_corr.mod_means{curr_day} =  mod_means;
    length_corr.mod_sems{curr_day} =  mod_sems;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);
% 
% % plot them over days
% figure; hold on
% for i = 1:length(days)
%    errorbar(length_corr.mod_means{i}(2:end),length_corr.mod_sems{i}(2:end) ...
%        ,'color',[i/length(days) 0 0]); 
%    errorbar(length_corr.unmod_means{i}(2:end),length_corr.unmod_sems{i}(2:end) ...
%        ,'color',[0 i/length(days) 0]);
% end

% plot them over days
for i = 1:length(days)
    figure;
   errorbar(length_corr.mod_means{i}(2:end),length_corr.mod_sems{i}(2:end)); 
   ylim([0 0.3])
end

%% Get number of active pyr/gad cells for each trial

clear all

animal = 'AP71';

% find days
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

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    %%%% Process
    loop_size = cellfun(@length,im.roi_trace_long_split);
    loop_size_sum = cumsum(loop_size);
    loop_frames = [0;loop_size_sum(8:8:end-7)];
    
    % get behavior file
    dir_currfolder = dir([data_path filesep day]);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
    
    % get xsg files
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    
    % get behavior data and xsg sample rate
    try
        % ignore if something's wrong with datafile (usually, >1 of them)
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
    catch me
        error('Something wrong with data file')
    end
    raw_bhv = saved_history.ProtocolsSection_parsed_events;
    num_trials = length(raw_bhv);
    
    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    % Only look at rewarded trials
    rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), raw_bhv));
    cued_movement = nan(num_trials,1);
    cue_frames = nan(num_trials,1);
    reward_frames = nan(num_trials,1);
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        
        % create movement velocity, active times
        % quantify movement before
        
        %lever_force_smooth = smooth(xsg_data.data.acquirer.trace_2,500);
        lever_force_resample = resample(xsg_data.data.acquirer.trace_2,1,10);
        % normalized cutoff freq is fraction of nyquist
        butterworth_freq = 20/500;
        [b a] = butter(4, butterworth_freq);
        % using filtfilt instead of filter allows for zero-phase shift
        lever_force_smooth = filtfilt(b,a,lever_force_resample);
        
        lever_velocity_resample = [0;diff(lever_force_smooth)];
        %lever_velocity_resample = resample(lever_velocity,1,10);
        %lever_force_smooth_resample = resample(lever_force_smooth,1,10);
        lever_velocity_resample_smooth = smooth(lever_velocity_resample,5);
        lever_velocity_resample_smooth(isnan(lever_velocity_resample)) = NaN;
        
        % get velocity envelope, smooth over 5 frames (~150ms)
        lever_velocity_hilbert = hilbert(lever_velocity_resample_smooth);
        lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
            conj(lever_velocity_hilbert));
        % at the moment - don't smooth
        lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,1);
        
        % define active parts of the lever, ignore small ones and
        % fill in gaps of < 100ms
        movethresh = 0.001;
        lever_active = lever_velocity_envelope_smooth > movethresh;
        lever_v_thresh = lever_velocity_envelope_smooth;
        lever_v_thresh(lever_v_thresh < 0.002) = 0;
        % get rid of short movements, find movement times again
        for short_move_erase = 1:3
            lever_active_switch = [0;diff(lever_active)];
            lever_active_starts = find(lever_active_switch == 1);
            lever_active_stops = find(lever_active_switch == -1);
            % get rid of unpaired starts/stops
            if lever_active_stops(1) < lever_active_starts(1);
                lever_active_stops(1) = [];
            end
            if lever_active_starts(end) > lever_active_stops(end)
                lever_active_starts(end) = [];
            end
            
            lever_active_movement_times = lever_active_stops - ...
                lever_active_starts;
            lever_active_intermovement_times = lever_active_starts(2:end) - ...
                lever_active_stops(1:end-1);
            
            if short_move_erase == 1
                % fill in short gaps between movements, but only if
                % there's a movement afterwards that's a certain time
                lever_active_fill = ...
                    lever_active_intermovement_times < 20 & ...
                    lever_active_movement_times(2:end) > 50;
                for i = find(lever_active_fill)'
                    lever_active(lever_active_stops(i): ...
                        lever_active_starts(i+1)) = 1;
                end
            else
                % erase movement that's too brief
                lever_active_erase = lever_active_movement_times < 100;
                for i = find(lever_active_erase')
                    lever_active(lever_active_starts(i): ...
                        lever_active_stops(i)) = 0;
                end
            end
        end
        
        % fill in short gaps between movements
        lever_active_intermovement_times = lever_active_starts(2:end) - ...
            lever_active_stops(1:end-1);
        lever_active_fill = lever_active_intermovement_times < 200;
        for i = find(lever_active_fill)'
            lever_active(lever_active_stops(i): ...
                lever_active_starts(i+1)) = 1;
        end
        
        lever_stopstart = [0;diff(lever_active)];
        lever_stopstart = [0;diff(lever_active)];
        
        % get cue and reward times in resampled xsg
        cue_sample = nan(length(curr_trial_list(:,2)),1);
        reward_start_sample = nan(length(curr_trial_list(:,2)),1);
        reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
        for curr_trial_idx = 1:length(curr_trial_list(:,2));
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % get trial start in dispatcher
            curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
            % get cue and reward sample in xsg (downsampled to ms)
            curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            cue_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_cue_sample_rel);
            curr_bhv_cue_frame_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(bhv.framerate);
            cue_frames(curr_trial) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(bhv.framerate) + ...
                curr_bhv_cue_frame_rel) + loop_frames(curr_xsg);
            
            curr_bhv_reward_start_sample_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_start_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_start_sample_rel);
            curr_bhv_reward_start_frame_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(bhv.framerate);
            reward_frames(curr_trial) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(bhv.framerate) + ...
                curr_bhv_reward_start_frame_rel) + loop_frames(curr_xsg);
            
            curr_bhv_reward_stop_sample_rel =  (raw_bhv{curr_trial}.states.reward(2) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_stop_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_stop_sample_rel);
        end
        
        
        for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
            
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials) || ...
                    length(lever_velocity_envelope_smooth) < cue_sample(curr_trial_idx+1)
                continue
            end
            
            % ignore if out of bounds
            if cue_sample(curr_trial_idx+1) > length( ...
                    lever_velocity_envelope_smooth)
                continue
            end
            
            % classify trial cue/movement overlap
            % rewarded-movement/cue overlap can't be > 50ms
            rewarded_move_start = find(lever_active(1: ...
                reward_start_sample(curr_trial_idx)) == 0,1,'last');
            cued_movement(curr_trial) = rewarded_move_start > cue_sample(curr_trial_idx)-50;
        end
    end
    
    roi_trace_df_smooth = zeros(size(im.roi_trace_df));
    for curr_cell = 1:size(im.roi_trace_df,1)
        
        % normalized cutoff freq is fraction of nyquist
        butterworth_freq = (0.1*bhv.framerate)/(bhv.framerate/2);
        [b a] = butter(4, butterworth_freq);
        % using filtfilt instead of filter allows for zero-phase shift
        roi_trace_df_smooth(curr_cell,:) = filtfilt(b,a,im.roi_trace_df(curr_cell,:));
    end
    
    % Get the maximum response of cell during each trial
    cell_trial_max = nan(size(im.roi_trace_df,1),num_trials);
    for i = 1:num_trials-1
        if ~any(isnan(cue_frames(i:i+1)));
            cell_trial_max(:,i) = max(roi_trace_df_smooth(:, ...
                cue_frames(i):cue_frames(i+1)),[],2);
        end
    end
    
    
    roilabel_filename = [data_path filesep animal '_roi_template' filesep...
        animal '_roilabels.roilabel'];
    load(roilabel_filename,'-MAT')
    cells = [1:length(roi_labels)];
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    gad_thresh = 0.3;
    pyr_thresh = 0.5;
    
    bin_thresh = zeros(length(roi_labels),1);
    bin_thresh(gad_cells) = gad_thresh;
    bin_thresh(pyr_cells) = pyr_thresh;
    bin_thresh = repmat(bin_thresh,[1 size(cell_trial_max,2)]);
    bin_activity = cell_trial_max > bin_thresh;
    
    % get numbers of active pyr/gad cells
    gad_num_active = sum(bin_activity(gad_cells,:));
    pyr_num_active = sum(bin_activity(pyr_cells,:));
    
    %%%% Save
    grab_var.gad_num_active{curr_day} =  gad_num_active;
    grab_var.pyr_num_active{curr_day} =  pyr_num_active;
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);


% plot all
figure; hold on;
for i = 1:length(grab_var.gad_num_active)
    plot(grab_var.pyr_num_active{i}, grab_var.gad_num_active{i}, ...
        'color',[i/length(grab_var.gad_num_active) 0 0],'linestyle','.');
end
xlabel('Number of pyramidal cells active');
ylabel('Number of gad cells active');


%% Get population activity over days
clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
peak_matrix_full = cell(size(animals));
concat_peaks_10000 = cell(size(animals));

for curr_animal = 1:length(animals);
    clearvars -except animals peak_matrix_full concat_peaks_10000 curr_animal
    
    animal = animals{curr_animal};
    
    % find days
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
    
    grab_var = cell(size(days));
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    for curr_day = 1:length(days)
        %%%% Initialize
        
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        %%%% Process
                
        % get number of peaks for pyramidal cells
        [peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);
        num_peaks = [];
        num_peaks = cellfun(@length,peak_frames);
        peak_matrix = zeros(size(im.roi_trace_df));
        for i = pyr_cells
            peak_matrix(i,peak_frames{i}) = 1;
        end
        
        % for gad cells, use threshold
        for i = gad_cells
            temp_smooth = smooth(im.roi_trace_df(i,:),30,'loess')';
            noise_est = mean(abs(im.roi_trace_df(i,:) - temp_smooth));
            noise_thresh = noise_est*3;
            peak_matrix(i,:) = +(im.roi_trace_df(i,:) ...
                > noise_thresh);
        end
        
        %%%% Save
        num_frames = size(im.roi_trace_df,2);
        %grab_var.over_thresh{curr_day} =  roi_trace_df_thresh;
        grab_var.num_peaks{curr_day} =  peak_matrix;
        
        disp(['Finished day: ' day]);
        
    end
    disp(['Finished all.']);
    
    % save peak matricies
    peak_matrix_full{curr_animal} = grab_var.num_peaks;
    
    % create smoothed version for long term changes
    concat_peaks = [peak_matrix_full{curr_animal}{:}];
    box_filter = ones(1,10000);
    concat_peaks_10000{curr_animal} = conv2(concat_peaks,box_filter,'same');
    
    
    disp(['Finished animal ' animal])
    
end
save('/usr/local/lab/People/Andy/Data/figures/population_activity/peak_matrix_full.mat', ...
    'peak_matrix_full','concat_peaks_10000','animals','-v7.3')


%% Do PCA and plot above analysis

% generate random correlated data
x = randn(100,1000);
for i = 2:100 
    curr_corr = 1-(1/i);%1/100*i;%1-(1/i);%0.8
    x(i,:) = x(i-1,:)*curr_corr+x(i,:)*sqrt(1-curr_corr^2);
end

% get correlation of average activity between days/across animals
a = nan(15,15,8);
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
cells_all = [];
for i = 1:8
    animal = animals{i};

    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
   
    curr_cells = pyr_cells;
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2),peak_matrix_full{i},'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp_norm = bsxfun(@times,temp,1./max(temp,[],2));
    cells_all = [cells_all;temp(:,1:11)];
    temp_cc = corrcoef(temp);
    a(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
end
% Don't use day 1-2 of AP71
a(1:2,:,1) = NaN;
a(:,1:2,1) = NaN;

% get rid of all cells consecutively - find biggest contributors to CC
cc_cell_matrix = cell(length(animals),1);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};

    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
   
    curr_cells = pyr_cells;
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2),peak_matrix_full{curr_animal},'UniformOutput',false));
    temp(isnan(temp)) = 0;
    
    cc_cell_matrix{curr_animal} = zeros(size(temp,2),size(temp,2),size(temp,1));
    
    for cc1 = 1:size(temp,2)
        for cc2 = 1:size(temp,2)
            cc = corrcoef(temp(:,cc1),temp(:,cc2));
            cc_diff = zeros(size(temp,1),1);
            for i = 1:size(temp,1)
                curr_temp = temp;
                curr_temp(i,:) = [];
                curr_cc = corrcoef(curr_temp(:,cc1),curr_temp(:,cc2));
                
                cc_diff(i) = curr_cc(2);
            end
            cc_diff2 = cc_diff - cc(2);
            cc_cell_matrix{curr_animal}(cc1,cc2,:) = cc_diff2;
        end
    end

end


% get correlation of average activity between days/across animals in 1/4th
% of days
quarter_activity = cell(8,1);
a = nan(15,15,8);
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for i = 1:8
    animal = animals{i};

    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    cell_type = gad_cells;
    
    quarter_indx = cell(length(peak_matrix_full{i}),1);
    for j = 1:length(peak_matrix_full{i});
%         quarter_length = cellfun(@(x) reshape(1:size(x,2),ceil(size(x,2)/4),4),...
%             peak_matrix_full{i},'UniformOutput',false);
%         quarter_indx = cellfun(@(x) num2cell(x,1),quarter_length,'UniformOutput',false);
        quarter_length = 1:ceil(size(peak_matrix_full{i}{j},2)/4): ...
            size(peak_matrix_full{i}{j},2);
        quarter_indx{j}{1} = quarter_length(1):quarter_length(2);
        quarter_indx{j}{2} = quarter_length(2):quarter_length(3);
        quarter_indx{j}{3} = quarter_length(3):quarter_length(4);
        quarter_indx{j}{4} = quarter_length(4):size(peak_matrix_full{i}{j},2);
    end
    
    for j = 1:length(peak_matrix_full{i})
        quarter_activity{i}{j} = cellfun(@(x) ...
            nanmean(peak_matrix_full{i}{j}(cell_type,x),2),quarter_indx{j},'UniformOutput',false);
    end

%     temp = cell2mat(cellfun(@(x) nanmean(x(pyr_cells,:),2),peak_matrix_full{i},'UniformOutput',false));
%     temp(isnan(temp)) = 0;
%     temp_cc = corrcoef(temp);
%     a(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
end
% from the last stuff
pop_corr = nan(15*4,15*4,8);
for i = 1:8
    temp = [quarter_activity{i}{:}];
    temp = [temp{:}];
    temp_cc = corrcoef(temp);
    pop_corr(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
end
% Don't use day 1-2 of AP71
pop_corr(1:8,:,1) = NaN;
pop_corr(:,1:8,1) = NaN;

% compare different number of days surrounding
figure; hold on;
pop_corr_2wk = nanmean(pop_corr,3);
pop_corr_2wk = pop_corr_2wk(1:56,1:56);
pop_corr_2wk(pop_corr_2wk == 1) = NaN;
for i = 1:length(pop_corr_2wk)
   curr_corr = zeros(size(pop_corr_2wk));
   for j = 1:i+1
      %curr_corr = curr_corr + diag(diag(pop_corr_2wk,j-1),j-1); 
      curr_corr = curr_corr + diag(diag(pop_corr_2wk,-j+1),-j+1);
   end
   curr_corr(curr_corr == 0) = NaN;
   plot(nanmean(curr_corr),'color',[0 0 i/length(pop_corr_2wk)]);     
end

% compare different number of days surrounding (each number seperate val)
pop_corr_2wk = pop_corr(1:56,1:56,:);
pop_corr_2wk(pop_corr_2wk == 1) = NaN;
i = 4*4;
curr_corr_all = [];
for k = 1:size(pop_corr_2wk,3);
    curr_corr = zeros(56,56,1);
    for j = 1:i+1
        % past
        %curr_corr = curr_corr + diag(diag(pop_corr_2wk,j-1),j-1);
        % future
        curr_corr = curr_corr + diag(diag(pop_corr_2wk(:,:,k),-j+1),-j+1);
    end
    curr_corr(curr_corr == 0) = NaN;
    curr_corr_all = [curr_corr_all;curr_corr(:,1:44)];
end


%for curr_animal = 1:length(peak_matrix_full)
    
    animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
    animal = animals{curr_animal};
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    %concat_peaks = [peak_matrix_full{curr_animal}{:}];
    %box_filter = ones(1,10000);
    %concat_peaks_spread = conv2(concat_peaks,box_filter,'same');
    concat_peaks_spread = concat_peaks_10000{curr_animal};
    
    % do all the following for pyr and gad
    for curr_cell_type = 1:1;
        
        if curr_cell_type == 1
            cell_type = 'pyr';
            curr_cells = pyr_cells;
        elseif curr_cell_type == 2;
            cell_type = 'gad';
            curr_cells = gad_cells;
        end
        
        concat_peaks_downsamp = concat_peaks_spread(curr_cells,1:1000:end);
        
        % normalize each row by the max
        concat_peaks_downsamp_norm = bsxfun(@times,concat_peaks_downsamp, ...
            1./max(concat_peaks_downsamp,[],2));
        concat_peaks_downsamp_norm(isnan(concat_peaks_downsamp_norm)) = 0;
        
%         concat_peaks_downsamp_norm = bsxfun(@times,concat_peaks_downsamp, ...
%             1./sqrt(sum(concat_peaks_downsamp.^2,2)));
%         concat_peaks_downsamp_norm(isnan(concat_peaks_downsamp_norm)) = 0;
        
        concat_peaks_downsamp = concat_peaks_downsamp_norm;
        
        
        % kmeans sort traces
        kidx = kmeans(concat_peaks_downsamp,20);
        [temp sort_idx] = sort(kidx);
        figure;imagesc(concat_peaks_downsamp(sort_idx,:));colormap(gray);
        a = cellfun(@(x) size(x,2)/1000, peak_matrix_full{curr_animal});
        b = cumsum(a);
        for i = 1:length(b);
            line([b(i) b(i)],ylim,'color','b','linewidth',2);
        end
        title('kmeans sorted')
        % get the max of the mean activity in each sorted group
        [temp kmeans_peak_activity] = max(grpstats(concat_peaks_downsamp,kidx),[],2);
        numel_kmeans = grpstats(ones(length(kidx),1),kidx,'numel');
        
        
        [coeff score latents] = princomp(concat_peaks_downsamp);
        
        figure; hold on
        plot(coeff(:,1),'k');
        plot(coeff(:,2),'r');
        a = cellfun(@(x) size(x,2)/1000, peak_matrix_full{curr_animal});
        b = cumsum(a);
        for i = 1:length(b);
            line([b(i) b(i)],ylim,'color','b','linewidth',2);
        end
        title([animals{curr_animal} ': ' cell_type]);
        
%         % sort cells by PC coefficient difference
%         [pc_diff pc_idx] = sort(coeff(:,1) - coeff(:,2));
%         figure;imagesc(concat_peaks_downsamp(pc_idx,:));colormap(gray)
%         title([animals{curr_animal} ': ' cell_type ', PC1-PC2'])
        
        % sort cells by PC2
        [temp coeff2_idx] = sort(score(:,2));
        figure;imagesc(concat_peaks_downsamp(coeff2_idx,:));colormap(gray)
        title([animals{curr_animal} ': ' cell_type ', PC2'])
        
        % plot map of ROIs colored by PC2 score
        roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
            animal '_roi_template'];
        temp_dir = dir(roi_path);
        temp_day = temp_dir(3).name(1:6);
        roi_file = [temp_day '_' animal '_summed_50_template.roi'];
        rois = load([roi_path filesep roi_file],'-MAT');
        figure;
        roi_handle = [];
        for i = 1:length(curr_cells)
            curr_cell = curr_cells(i);
            roi_handle(i) = patch(rois.polygon.ROI{curr_cell}(:,1), ...
                -rois.polygon.ROI{curr_cell}(:,2), ...
                [score(i,2)]);
        end
        ylim([-512 0]);
        xlim([0 512]);
        caxis([-max(caxis) max(caxis)]);
        title('ROIs by PC2')
        
        % sort cells by peak of day activity
        smooth_filter = ones(1,50);
        concat_peaks_downsamp_smooth = conv2(concat_peaks_downsamp,smooth_filter,'same');
        [coeff2 score2 latents2] = princomp(concat_peaks_downsamp_smooth');
        [max_peaks max_idx] = max(concat_peaks_downsamp_smooth,[],2);
        [temp max_idx_sort] = sort(max_idx);
        figure; subplot(2,1,1);
        imagesc(concat_peaks_downsamp(max_idx_sort,:));colormap(gray)
        title([animals{curr_animal} ': ' cell_type ', Peak activity'])
        subplot(2,1,2);plot(sort(max_idx),1:length(max_idx),'k')
        for i = 1:length(b);
            line([b(i) b(i)],ylim,'color','r','linewidth',1,'linestyle','--');
        end
        set(gca,'YDir','reverse');
        xlim([1 size(concat_peaks_downsamp,2)]);
        
        % attempt to sort cells by start and length of activity
        activity_thresh = prctile(concat_peaks_downsamp_smooth(:),75);
        concat_peaks_thresh = concat_peaks_downsamp_smooth > activity_thresh;
        activity_starts = zeros(size(concat_peaks_thresh,1),1);
        activity_stops = zeros(size(concat_peaks_thresh,1),1);
        activity_lengths = zeros(size(concat_peaks_thresh,1),1);
        % get longest run - why is this always difficult?
        for i = 1:size(concat_peaks_thresh,1)
            concat_peaks_thresh_start = find(diff([0 ...
                concat_peaks_thresh(i,:) 0]) == 1);
            concat_peaks_thresh_stop = find(diff([0 ...
                concat_peaks_thresh(i,:) 0]) == -1);
            activity_run = concat_peaks_thresh_stop -concat_peaks_thresh_start;
            if ~isempty(activity_run)
                [activity_lengths(i) run_idx] = max(activity_run);
                activity_starts(i) = concat_peaks_thresh_start(run_idx); 
                activity_stops(i) = concat_peaks_thresh_stop(run_idx);   
            end
        end
        [temp run_sort_idx] = sort(activity_starts);
        figure;imagesc(concat_peaks_downsamp(run_sort_idx,:));colormap(gray);
        title([animals{curr_animal} ': ' cell_type ', Start activity'])
        [temp run_sort_idx] = sort(activity_stops);
        figure;imagesc(concat_peaks_downsamp(run_sort_idx,:));colormap(gray);
        title([animals{curr_animal} ': ' cell_type ', Stop activity'])
        
        % normalize each row by the max
        concat_peaks_downsamp_norm = bsxfun(@times,concat_peaks_downsamp, ...
           1./max(concat_peaks_downsamp,[],2)); 
        
        
    end
%     figure; hold on
%     plot(coeff(:,1),'.k','MarkerSize',20);
%     plot(coeff(:,2),'.r','MarkerSize',5);
%     for i = 1:length(coeff(:,1));
%         line([i i], [coeff(i,1) coeff(i,2)])
%     end
%end

%% Trial-by-trial analysis

clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
trial_activity = struct;
trial_type = struct;
for curr_animal = 1:length(animals);
    clearvars -except animals trial_type trial_activity curr_animal
    
    animal = animals{curr_animal};
      
    % find days
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
    
    grab_var = cell(size(days));
    
    for curr_day = 1:length(days)
        %%%% Initialize
        clearvars -except animal days curr_day grab_var ...
            animals trial_type trial_activity curr_animal
        
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        % load ROI labels and identify cell type
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
        roilabel_name = [animal '_roilabels.roilabel'];
        load([analysis_path filesep roilabel_name],'-MAT')
        cells = 1:length(roi_labels);
        gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
        pyr_cells = cells(~ismember(cells,gad_cells));
        filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
        
        %%%% Process
        
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            error('Something wrong with data file')
        end
        raw_bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(raw_bhv);
        
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        % Only look at rewarded trials
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), raw_bhv);
        cued_movement = false(num_trials,1);
        cue_frames = nan(num_trials,1);
        reward_frames = nan(num_trials,1);
        lever_active_full = cell(length(xsg_filenames),1);
        rewarded_move_starts = nan(num_trials,1);
        rewarded_move_stops = nan(num_trials,1);
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            % create movement velocity, active times
            % quantify movement before
            
            [lever_active lever_force_smooth] = AP_parseLeverMovement(xsg_data);
            
            % save the active lever portions in frames
            curr_lever_active_frames = round((find(lever_active)/1000) ...
              *bhv.framerate);
            lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
                curr_lever_active_frames <= loop_size(curr_xsg))) + loop_frames(curr_xsg);
            
            % get cue and reward times in resampled xsg
            cue_sample = nan(length(curr_trial_list(:,2)),1);
            reward_start_sample = nan(length(curr_trial_list(:,2)),1);
            reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
            for curr_trial_idx = 1:length(curr_trial_list(:,2));
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if not imaged trial
                if curr_trial > length(bhv.imaged_trials) || ...
                    bhv.imaged_trials(curr_trial) == 0
                    continue
                end
                
                % get trial start in dispatcher
                curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
                % get cue and reward sample in xsg (downsampled to ms)
                curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                cue_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_cue_sample_rel);
                curr_bhv_cue_frame_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(bhv.framerate);
                cue_frames(curr_trial) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(bhv.framerate) + ...
                    curr_bhv_cue_frame_rel) + loop_frames(curr_xsg);
                
                % find reward times if rewarded trial
                if ~rewarded_trials(curr_trial)
                    continue
                end
                
                curr_bhv_reward_start_sample_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_start_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_start_sample_rel);
                curr_bhv_reward_start_frame_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(bhv.framerate);
                reward_frames(curr_trial) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(bhv.framerate) + ...
                    curr_bhv_reward_start_frame_rel) + loop_frames(curr_xsg);
                
                curr_bhv_reward_stop_sample_rel =  (raw_bhv{curr_trial}.states.reward(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_stop_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_stop_sample_rel);
            end
            
            % find cued-movement rewarded trials, and rewarded movement
            % start/stop times
            for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
                
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial or out of bounds
                if ~rewarded_trials(curr_trial) || ...
                        bhv.imaged_trials(curr_trial) == 0
                    continue
                end
                
                % classify trial cue/movement overlap
                % rewarded-movement/cue overlap can't be > 150ms
                rewarded_move_start_rel = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');              
                cued_movement(curr_trial) = rewarded_move_start_rel > ...
                    cue_sample(curr_trial_idx)-150;
                
                % get absolute rewarded move time
                rewarded_move_starts(curr_trial) = round(find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last') ...
                    * bhv.framerate/1000) + loop_frames(curr_xsg);
                
                rewarded_move_stops(curr_trial) = round((find(lever_active(...
                    reward_start_sample(curr_trial_idx):end) == 0,1,'first') + ...
                    reward_start_sample(curr_trial_idx)) * bhv.framerate/1000) ...
                    + loop_frames(curr_xsg);

            end
            
        end
        
        roi_trace_df_smooth = zeros(size(im.roi_trace_df));
        for curr_cell = 1:size(im.roi_trace_df,1)           
            % normalized cutoff freq is fraction of nyquist
            butterworth_freq = 5/(bhv.framerate/2);
            [b a] = butter(4, butterworth_freq);
            % using filtfilt instead of filter allows for zero-phase shift
            roi_trace_df_smooth(curr_cell,:) = filtfilt(b,a,im.roi_trace_df(curr_cell,:));
        end
        
        lever_moving = vertcat(lever_active_full{:});
        lever_quiescent = setdiff(1:size(im.roi_trace_df,2),lever_moving);
        
        % Get the maximum response of cell during each trial
        cell_trial_max_moving = nan(size(im.roi_trace_df,1),num_trials);
        cell_trial_max_quiescent = nan(size(im.roi_trace_df,1),num_trials);
        for i = 1:num_trials-1
            if ~any(isnan(cue_frames(i:i+1))) && ...
                    rewarded_trials(i);
                curr_frames = cue_frames(i):cue_frames(i+1);
                %curr_moving_frames = intersect(rewarded_move_starts(i): ...
                %    rewarded_move_stops(i),curr_frames);
                curr_moving_frames = rewarded_move_starts(i): ...
                    rewarded_move_stops(i);
                curr_quiescent_frames = intersect(lever_quiescent,curr_frames);
                
                if ~isempty(curr_moving_frames)
                    cell_trial_max_moving(:,i) = max(im.roi_concat_oopsi(:, ...
                        curr_moving_frames),[],2);
                end
                
                if ~isempty(curr_quiescent_frames)
                    cell_trial_max_quiescent(:,i) = max(im.roi_concat_oopsi(:, ...
                        curr_quiescent_frames),[],2);
                end
            end
        end
        
        % get presence or absence of peak during trial - DO DIFFERENTLY FOR
        % PYR / GAD, pyr = peak, gad = >thresh
        peak_matrix = AP_caEvents(im.roi_trace_df,pyr_cells,gad_cells);
%         [peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);
%         peak_matrix = zeros(size(im.roi_trace_df));
%         for curr_cell = pyr_cells;
%            peak_matrix(curr_cell,peak_frames{curr_cell}) = 1; 
%         end
%         for curr_cell = gad_cells
%             temp_smooth = smooth(im.roi_trace_df(curr_cell,:),30,'loess')';
%             noise_est = mean(abs(im.roi_trace_df(curr_cell,:) - temp_smooth));
%             noise_thresh = noise_est*4;
%             peak_matrix(curr_cell,(im.roi_trace_df(curr_cell,:) ...
%                 > noise_thresh)) = 1;
%         end
        
        cell_trial_peak_moving = nan(size(im.roi_trace_df,1),num_trials);
        cell_trial_peak_quiescent = nan(size(im.roi_trace_df,1),num_trials);
        for i = 1:num_trials-1
            if ~any(isnan(cue_frames(i:i+1))) && ...
                    rewarded_trials(i);
                curr_frames = cue_frames(i):cue_frames(i+1);
                %curr_moving_frames = intersect(rewarded_move_starts(i): ...
                %    rewarded_move_stops(i),curr_frames);
                curr_moving_frames = rewarded_move_starts(i): ...
                    rewarded_move_stops(i);
                curr_quiescent_frames = intersect(lever_quiescent,curr_frames);
                
                if ~isempty(curr_moving_frames)
                    cell_trial_peak_moving(:,i) = max(peak_matrix(:, ...
                        curr_moving_frames),[],2);
                end
                
                if ~isempty(curr_quiescent_frames)
                    cell_trial_peak_quiescent(:,i) = max(peak_matrix(:, ...
                        curr_quiescent_frames),[],2);
                end
            end
        end
        
        cell_trial_df_moving = nan(size(im.roi_trace_df,1),num_trials);
        cell_trial_df_quiescent = nan(size(im.roi_trace_df,1),num_trials);
        for i = 1:num_trials-1
            if ~any(isnan(cue_frames(i:i+1))) && ...
                    rewarded_trials(i);
                curr_frames = cue_frames(i):cue_frames(i+1);
                %curr_moving_frames = intersect(rewarded_move_starts(i): ...
                %    rewarded_move_stops(i),curr_frames);
                curr_moving_frames = rewarded_move_starts(i): ...
                    rewarded_move_stops(i);
                curr_quiescent_frames = intersect(lever_quiescent,curr_frames);
                
                if ~isempty(curr_moving_frames)
                    cell_trial_df_moving(:,i) = max(im.roi_trace_df(:, ...
                        curr_moving_frames),[],2);
                end
                
                if ~isempty(curr_quiescent_frames)
                    cell_trial_df_quiescent(:,i) = max(im.roi_trace_df(:, ...
                        curr_quiescent_frames),[],2);
                end
            end
        end
        
        % ignore catch trials (or make new matrix of zeros, if no catch
        % trials)
        if ~isfield(bhv,'catch_trials');
            bhv.catch_trials = zeros(length(cued_movement),1);
        end
        
        
        %%%% Save
        % save oopsi
        grab_var.cell_trial_max.cued.moving{curr_day} = ...
            cell_trial_max_moving(:,cued_movement & ~bhv.catch_trials);
        grab_var.cell_trial_max.cued.quiescent{curr_day} = ...
            cell_trial_max_quiescent(:,cued_movement & ~bhv.catch_trials);
        
        grab_var.cell_trial_max.uncued.moving{curr_day} = ...
            cell_trial_max_moving(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        grab_var.cell_trial_max.uncued.quiescent{curr_day} = ...
            cell_trial_max_quiescent(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        
        % save peak
        grab_var.cell_trial_peak.cued.moving{curr_day} = ...
            cell_trial_peak_moving(:,cued_movement & ~bhv.catch_trials);
        grab_var.cell_trial_peak.cued.quiescent{curr_day} = ...
            cell_trial_peak_quiescent(:,cued_movement & ~bhv.catch_trials);
        
        grab_var.cell_trial_peak.uncued.moving{curr_day} = ...
            cell_trial_peak_moving(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        grab_var.cell_trial_peak.uncued.quiescent{curr_day} = ...
            cell_trial_peak_quiescent(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        
        % save df
        grab_var.cell_trial_df.cued.moving{curr_day} = ...
            cell_trial_df_moving(:,cued_movement & ~bhv.catch_trials);
        grab_var.cell_trial_df.cued.quiescent{curr_day} = ...
            cell_trial_df_quiescent(:,cued_movement & ~bhv.catch_trials);
        
        grab_var.cell_trial_df.uncued.moving{curr_day} = ...
            cell_trial_df_moving(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        grab_var.cell_trial_df.uncued.quiescent{curr_day} = ...
            cell_trial_df_quiescent(:,~cued_movement & ~bhv.catch_trials & rewarded_trials);
        
        % save trial numbers
        grab_var.cued_movement{curr_day} = cued_movement;
        grab_var.rewarded_trials{curr_day} = rewarded_trials;
        grab_var.catch_trials{curr_day} = bhv.catch_trials;
        
        disp(['Finished day: ' day]);
        
    end
    disp(['Finished all ' animal]);
    
    trial_activity.oopsi{curr_animal} = grab_var.cell_trial_max;
    trial_activity.peak{curr_animal} = grab_var.cell_trial_peak;
    trial_activity.df{curr_animal} = grab_var.cell_trial_df;
    trial_type.cued_movement{curr_animal} = grab_var.cued_movement;
    trial_type.rewarded_trials{curr_animal} = grab_var.rewarded_trials;
    trial_type.catch_trials{curr_animal} = grab_var.catch_trials;
end
save('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/trial_activity_caEvents', ...
    'trial_type','trial_activity');

%% plot above

% correlate mean of day with mean of other days
cued_moving = nan(15,15,8);
uncued_moving = nan(15,15,8);
cued_quiescent = nan(15,15,8);
uncued_quiescent = nan(15,15,8);
moving_quiescent = nan(15,15,8);
num_moving = {};
num_quiescent = {};
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for i = 1:8
    animal = animals{i};
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat']);
    
    % pick max/peak as measure of activity
    curr_measure = cell_trial_peak_full{i};
    curr_cells = setdiff(pyr_cells,filled_cells);
    
    % exclude non-move cells
    move_cells = vertcat(classified_cells.move_cells_oopsi{:});
    exclude_cells = [];%~move_cells(:,curr_cells)';
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2), ...
        curr_measure.cued.moving,'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp(exclude_cells) = 0;
    temp_cc = corrcoef(temp);
    cued_moving(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2), ...
        curr_measure.uncued.moving,'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp_cc = corrcoef(temp);
    uncued_moving(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2), ...
        curr_measure.cued.quiescent,'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp(exclude_cells) = 0;
    temp_cc = corrcoef(temp);
    cued_quiescent(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
    
    temp = cell2mat(cellfun(@(x) nanmean(x(curr_cells,:),2), ...
        curr_measure.uncued.quiescent,'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp(exclude_cells) = 0;
    temp_cc = corrcoef(temp);
    uncued_quiescent(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
    
    temp = cell2mat(cellfun(@(x,y) ...
        nanmean([x(curr_cells,:) y(curr_cells,:)],2), ...
        curr_measure.cued.moving, ...
        curr_measure.uncued.moving,'UniformOutput',false));
    temp(isnan(temp)) = 0;
    temp(exclude_cells) = 0;
    temp2 = cell2mat(cellfun(@(x,y) ...
        nanmean([x(curr_cells,:) y(curr_cells,:)],2), ...
        curr_measure.cued.quiescent, ...
        curr_measure.uncued.quiescent,'UniformOutput',false));
    temp2(isnan(temp)) = 0;
    temp(exclude_cells) = 0;
    temp_cc = corrcoef([temp temp2]);
    moving_quiescent(1:length(temp_cc)/2,1:length(temp_cc)/2,i) = ...
        temp_cc(length(temp_cc)/2+1:end,1:length(temp_cc)/2);
    
    num_moving{i} = nan(1,size(temp,2));
    num_quiescent{i} = nan(1,size(temp2,2));
    
    num_moving{i}(1:size(temp,2)) = nanmean(temp);
    num_quiescent{i}(1:size(temp2,2)) = nanmean(temp2);
end
% for goddamn AP71 bad computer
cued_moving(:,1:2,1) = NaN;
uncued_moving(:,1:2,1) = NaN;
cued_quiescent(:,1:2,1) = NaN;
uncued_quiescent(:,1:2,1) = NaN;
moving_quiescent(:,1:2,1) = NaN;
cued_moving(1:2,:,1) = NaN;
uncued_moving(1:2,:,1) = NaN;
cued_quiescent(1:2,:,1) = NaN;
uncued_quiescent(1:2,:,1) = NaN;
moving_quiescent(1:2,:,1) = NaN;
num_moving{1}(1:2) = NaN;
num_quiescent{1}(1:2) = NaN;


% get mean correlation of all trials with all other trials
a = nan(15,15,8);
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for i = 1:8
    animal = animals{i};
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));   
    temp = cell2mat(cellfun(@(x) nanmean(x(pyr_cells,:),2),...
        cell_trial_max_full{i}.cued.moving,'uni',false));
    temp(isnan(temp)) = 0;
    temp_cc = corrcoef(temp);
    a(1:length(temp_cc),1:length(temp_cc),i) = temp_cc;
end


%% Batch template
clear all

animal = 'AP76';

% find days
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

grab_var = cell(size(days));

for curr_day = 1:length(days)
    %%%% Initialize
    clearvars -except animal days curr_day grab_var
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_name = [day '_' animal '_bgROI_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    %%%% Process
    
    %%%% Save
    grab_var{curr_day} =  [];
    
    disp(['Finished day: ' day]);
    
end
disp(['Finished all.']);
