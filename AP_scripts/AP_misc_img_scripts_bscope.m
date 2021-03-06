%% Get rid of loop-start fluorscence spike, for correlations
% NOTE! The behavior assumes that these frames are here!

% Go through roi trace split, assume 8 files per loop and cut out first 10
% frames (should be ~ 5 that are artifact)
roi_trace_loopcut = [];
for curr_file = 1:length(roi_trace_long_split)
    curr_trace = [];
    curr_trace = roi_trace_long_split{curr_file};
    
    if mod(curr_file-1,8) == 0
        curr_trace = curr_trace(:,10:end);
    end
    
    roi_trace_loopcut = [roi_trace_loopcut curr_trace];
end

% get df/f normalization
roi_trace_loopcut_df = AP_baselineEstimation(roi_trace_loopcut,bhv.framerate);

% loess smooth 
smooth_size = 30;
roi_trace_loopcut_df_smooth = nan(size(roi_trace_loopcut_df));
for i = 1:size(roi_trace_loopcut_df,1);
    roi_trace_loopcut_df_smooth(i,:) = smooth(roi_trace_loopcut_df(i,:),smooth_size,'loess');
end

% Create a thresholded trace
roi_trace_loopcut_df_smooth_cutoff = zeros(size(roi_trace_loopcut_df_smooth));
for curr_cell = 1:size(roi_trace_loopcut_df_smooth,1)
    
    % call the minimum of the raw the low cutoff
    % estimate the noise from smoothed - raw
    noise_est = abs(roi_trace_loopcut_df(curr_cell,:) - roi_trace_loopcut_df_smooth(curr_cell,:));
    cutoff = 4*mean(noise_est);
    curr_trace_cutoff = nan(size(roi_trace_loopcut_df_smooth(curr_cell,:)));
    over_cutoff_indx = roi_trace_loopcut_df_smooth(curr_cell,:) > cutoff;
    roi_trace_loopcut_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
        roi_trace_loopcut_df_smooth(curr_cell,over_cutoff_indx);
end

%% Go through all saved concat files and make behavior, df/f traces, loopcut, oopsi

% update this! make this also produce a loop-stitching corrected trace

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

for curr_day = 1:length(days)
        
    clearvars -except animal days curr_day
    
    day = days{curr_day};
    
    % load raw trace
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal ...
        filesep animal '_roi_template'];
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
     % save imaging variables into image structure
     if exist('roi_trace_long')
         im.roi_trace_long = roi_trace_long;
         im.roi_trace_long_split = roi_trace_long_split;
     end
    
    % get behavior
    bhv = AP_getBehavior_leverBscope(animal,day,im.roi_trace_long);

    % baseline the trace
    im.roi_trace_df = AP_baselineEstimation(im.roi_trace_long,bhv.framerate);
     
     % this was an old,unused variable, don't save this
     if isfield(im,'roi_trace_df_loopcut')
         im = rmfield(im,'roi_trace_df_loopcut');
     end
     
% %%%%%
% %%%%% do deconvolution
% %%%%%
% 
% % load labels
% label_filename = [animal '_roilabels.roilabel'];
% load([analysis_path filesep label_filename],'-MAT');
% gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
% 
% % set options/parameters
% clear V
% V.dt = 1/bhv.framerate;
% V.est_sig = 1;
% V.est_lam = 1;
% V.est_gam = 1;
% V.est_b = 1;
% V.est_a = 1;
% 
% clear P
% tau = 1;
% P.gam = 1-(V.dt/tau);
% 
% roi_oopsi = zeros(size(roi_trace_df));
% 
% inactive_thresh = 0.1; % df threshold for stitching
% num_loops = length(roi_trace_long_split)/8;
% file_lengths = cellfun(@length,roi_trace_long_split);
% file_lengths_loops = reshape(file_lengths,8,num_loops);
% loop_lengths = sum(file_lengths_loops);
% 
% % Break up baseline corrected trace into loops
% roi_trace_df_split = mat2cell(im.roi_trace_df,size(im.roi_trace_df,1), ...
%     loop_lengths);
% split_indx = mat2cell(1:size(im.roi_trace_df,2),1, ...
%     loop_lengths);
% 
% disp('Deconvolving traces...')
% fprintf('%2d',0);
% for curr_cell = 1:size(roi_trace_df,1)
%     
%     % Stitch together loops at inactive parts, get index of stitched frames
%  % Skip the first 10 frames in case of weird things
%         loop_start_stitches = cellfun(@(x) find(x(curr_cell,11:end) ...
%             < inactive_thresh & x(curr_cell,:) > 0,1)+10, ...
%             roi_trace_df_split,'UniformOutput',false);
%     loop_stop_stitches = cellfun(@(x) find(x(curr_cell,:) ...
%         < inactive_thresh & x(curr_cell,:) > 0,1,'last'), ...
%         roi_trace_df_split,'UniformOutput',false);
%     F_stitched = cellfun(@(x,y,z) x(curr_cell,y:z),roi_trace_df_split,...
%         loop_start_stitches,loop_stop_stitches,'UniformOutput',false);
%     stitching_indx = cellfun(@(x,y,z) x(y:z), split_indx, ...
%         loop_start_stitches, loop_stop_stitches,'UniformOutput',false);
%     
%     F = [F_stitched{:}];
%     stitch_indx_long = [stitching_indx{:}];
%     
%     if any(isnan(roi_trace_df(curr_cell,:)))
%         continue
%     end
%     
%     % jovo added this
%     F = F - min(F);
% 
%     if ismember(curr_cell,gad_cells)
%         P.lam = 0.5;
%     else
%         P.lam = 0.5;
%     end
% 
%     % deconvolve filtered trace
%     V = struct;
%     V.dt = 1/bhv.framerate;
%     [n_est P_est V_est] = fast_oopsi(F,V,P);
%     roi_oopsi(curr_cell,stitch_indx_long) = n_est/max(n_est);
%     fprintf('%c%c%2d',8,8,round(100*curr_cell/size(roi_trace_df,1)));
% end
% 
% im.oopsi_lambda_0_5 = roi_oopsi;

save([analysis_path filesep analysis_name], ...
    'bhv', ...
    'im');

disp(['Finished day ' day])

end
disp('Done all')

%% Go through all df traces, do fast oopsi

clear all

animal = 'AP71';

days = {120606 120607 120608 120609 120610 120611 120612 ...
    120613 120614 120615 120616 120617 120618 120619};

for curr_day = 1:length(days)
        
    clearvars -except animal days curr_day
    
    day = num2str(days{curr_day});
    
    % load raw trace
    analysis_name = [day '_' animal '_summed_50_template_analysis.mat'];
    load(analysis_name);
    
    roi_oopsi = zeros(size(roi_trace_df));
    filt_trace = zeros(size(roi_trace_df));
    for curr_cell = 1:size(roi_trace_df,1)
        if all(isnan(roi_trace_df(curr_cell,:)));
            continue
        end
%         % old: 0.1
%         % normalized cutoff freq is fraction of nyquist
%         butterworth_freq = (0.05*bhv.framerate)/(bhv.framerate/2);
%         [b a] = butter(4, butterworth_freq);
%         % using filtfilt instead of filter allows for zero-phase shift
%         curr_filt_trace = filtfilt(b,a,roi_trace_df(curr_cell,:));
%         filt_trace(curr_cell,:) = curr_filt_trace;
        % deconvolve filtered trace
        V = struct;
        V.dt = 1/bhv.framerate;
        [n_best P_best V_est] = fast_oopsi(roi_trace_df(curr_cell,:),V);
        roi_oopsi(curr_cell,:) = n_best;
        disp(curr_cell/size(roi_trace_df,1))
    end
    
    % cut out frames around loop junctions (high fluor at new loop)
    total_frames = cumsum(cellfun(@length, roi_trace_long_split));
    loop_frames = [1;total_frames(8:8:end)];
    loopcut_frames = repmat(loop_frames,1,15);
    loopcut_addon = repmat([-5:9],size(loopcut_frames,1),1);
    loopcut_frames = [loopcut_frames+loopcut_addon]';
    loopcut_frames = loopcut_frames(:);
    loopcut_frames(loopcut_frames < 1) = [];
    loopcut_frames(loopcut_frames > size(roi_trace_long,2)) = [];
    roi_oopsi(:,loopcut_frames) = 0;
    im.fast_oopsi = roi_oopsi;
    
    save(analysis_name, ...
        'roi_trace_long', ...
        'roi_trace_long_split', ...
        'roi_trace_df', ...
        'bhv', ...
        'im');
    
    disp(['Finished day ' day])
    
end

%% Get centers/distances of all ROIs

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

%% Get trace and behavior times

animal  = 'AP61';
day = '120417';

data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];

disp('Getting behavior')

disp('Getting filenames')

% Get behavior filename
dir_currfolder = dir(data_path);
dir_filenames = {dir_currfolder.name};
bhv_file = strncmp('data_@',dir_filenames,6);
bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
bhv_fullfilename = [data_path filesep bhv_filename];

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

disp('Getting framerate from first file')
% Get framerate (pick first image file, read header)
dir_currfolder = dir(data_path);
dir_filenames = {dir_currfolder.name};
tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
img_filename = dir_filenames(tiff_filename_indx);
img_info = imfinfo([data_path filesep img_filename{1}]);
img_info = img_info(1).ImageDescription;
[img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
framerate = str2num(img_value{framerate_indx});
% AP60 is an exception
if strcmp(img_filename{1}(8:11),'AP60')
    framerate = 28.7224;
    disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
end

% Convert behavior to frames
disp('Converting behavior times to image frames')
[bhv_frames imaged_trials xsg_trials] = AP_getBehavior_dispatcher(framerate,bhv_fullfilename,xsg_fullfilenames);

num_trials = length(bhv_frames);

% Compensate for loops to have bhv_frames be in continuous time

% At the moment, concatenate all and pretend no interruption
if mod(size(img_filename,2),8) ~= 0
    error('Number of files not divisible by 8!')
end

% Check all last frames to watch for dropped frames
disp('Double checking frame number and logging dropped frames at ends')
last_frames = zeros(size(img_filename,2)/8,1);
for i = 1:size(img_filename,2)/8
    last_frames(i) = length(imfinfo([data_path filesep img_filename{8*i}]));
end

% Sanity check that the number of frames matches calculated number
if exist('roi_trace_long','var')
   total_est_frames = 500*7*size(img_filename,2)/8+sum(last_frames);
   if size(roi_trace_long,2) ~= total_est_frames
       error('Incorrect number of total estimated frames')
   end
end

% Prepare vector for adding frames based on loop number
add_frames = zeros(num_trials,1);
for i = find(imaged_trials)'
    if xsg_trials(i) <= 1
        add_frames(i) = 0;
    else
        add_frames(i) = (xsg_trials(i)-1)*500*7+sum(last_frames(1:xsg_trials(i)-1));
    end
end

% get relevant behaviors, add appropriate number of frames for concat
disp('Pulling out relevant behaviors')
all_reward = cellfun(@(x) x.states.reward, bhv_frames,'UniformOutput',0);
all_cue = cellfun(@(x) x.states.cue, bhv_frames,'UniformOutput',0);
all_iti = cellfun(@(x) x.states.iti, bhv_frames,'UniformOutput',0);
all_licks = cellfun(@(x) x.pokes.C(:,1), bhv_frames,'UniformOutput',0);

% all_reward = all_reward + add_frames;
% all_cue = all_cue + add_frames;
% all_iti = all_iti + add_frames;
% all_licks = all_licks + add_frames;

% get trial starts and ends, add appropriate number of frames for concat
all_trial_start = cellfun(@(x) x.states.state_0(1,2), bhv_frames,'UniformOutput',0);
all_trial_end = cellfun(@(x) x.states.state_0(2,1), bhv_frames,'UniformOutput',0);

% all_trial_start = all_trial_start + add_frames;
% all_trial_end = all_trial_end + add_frames;

reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];

% NOTE ON SETUP: C = lickport, L = lever
% Index if/which lever press ended the cue state
lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), bhv_frames,'UniformOutput',0);

% get lever data, downsample
lever_force_split = cell(length(xsg_filename_all),1);

lever_force = [];
trial_list = [];
lever_force = [];
bitcode_trace = [];
trial_stitch = [];

disp('Getting lever traces')
for i = 1:length(xsg_filename_all);
    xsg = load([xsg_path filesep xsg_filename_all{i}],'-MAT');
    
    % Get all trial offsets in samples
    xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;

    % Cut out pieces of the lever trace that weren't imaged (dropped
    % frames) - Assume they were clipped at the end, assume that there are
    % 500 frames total normally too
    curr_loop_frames = 500*7+last_frames(i);
    %curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
    curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
    curr_lever_trace = xsg.data.acquirer.trace_2(1:curr_imaged_lever_samples);
    
    lever_force_split{i} = xsg.data.acquirer.trace_2;
    lever_force = [lever_force; curr_lever_trace];
    
    % could have getBehavior_dispatcher do this, and not nec. anymore
%     trial_list = AP_ReadBitCode([xsg_path filesep xsg_filename_all{i}]);
%     
%     if ~isempty(trial_list)
%         trial_stitch = [trial_stitch;trial_list];
%     end
end

% resample lever force to be 1:1 with numframes
disp('Resampling lever traces')
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

% Don't need this anymore with getBehavior_dispatcher - already in frames
% % Get trial offsets (seconds)
% xsg_bhv_offset = nan(num_trials,1);
% for curr_trial = 1:num_trials
%     xsg_trial_indx = find(trial_stitch(:,2) == curr_trial,1);
%     if ~isempty(xsg_trial_indx)
%         xsg_bhv_offset(curr_trial) = ...
%             trial_trial_stitch(xsg_trial_indx,1) - ...
%             all_trial_start{curr_trial}(1);
%     else
%         continue;
%     end
% end

% Prepare behavior frame times for concatenation
disp('Preparing behavior for concatenation')
cue_frames = [];
reward_frames = [];
lick_frames = [];
for curr_trial = find(imaged_trials)';    
    % relative cue time
    cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + add_frames(curr_trial))];
    % relative reward time
    if ~isempty(all_reward{curr_trial})
        reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + add_frames(curr_trial))];
    end
    % relative lick time
    lick_frames = [lick_frames;(all_licks{curr_trial} + add_frames(curr_trial))];
end

% Prepare bhv_frames by adding previous xsg frames for concatenation

disp('Finished getting behavior.')

%% Plot behavior
% assuming last cell run
% Reaction times, licking, press stereotypy
% most of this is now obsolete

% make lick raster
figure; 
surround_back = 20;
surround_forward = 180;

catch_trials = zeros(size(bhv.cue_frames));
if isfield(bhv,'catch_trials')
    catch_trials = bhv.catch_trials(bhv.imaged_trials);
end

reward_frames = bhv.reward_frames(~catch_trials);
catch_frames = bhv.reward_frames(catch_trials);

subplot(2,1,1); hold on;
for i = 1:length(bhv.reward_frames)
    temp_lick = bhv.lick_frames - bhv.reward_frames(i);
    plot(temp_lick(temp_lick > -surround_back & ...
        temp_lick < surround_forward),i,'.k')
end
title('Lick raster - rewarded trials')

subplot(2,1,2); hold on;


catch_trials = zeros(size(bhv.cue_frames));
if isfield(bhv,'catch_trials')
    catch_trials = bhv.catch_trials(bhv.imaged_trials);
end
title('Lick raster - catch trials')

% reaction times: reward - cue
rxn_time = zeros(length(all_cue),1);
for i = 1:length(all_cue)
    if ~isempty(all_reward{i})
        rxn_time(i) = all_reward{i}(1) - all_cue{i}(1);
    end
end

% relative lick times to reward
lick_times = cellfun(@(x) x.pokes.C(:,1), bhv_frames,'UniformOutput',0);
for i = 1:length(lick_times)
    if ~isempty(lick_times{i}) && reward_trials(i) == 1;
        lick_times{i} = lick_times{i} - all_reward{i}(1);
    end
end

figure; hold on
for i = find(imaged_trials)';%1:length(lick_times)
    if ~isempty(lick_times{i}) && reward_trials(i) == 1;
        plot(lick_times{i},i,'.k');
        plot(all_trial_start{i} - all_reward{i}(1),i,'.g');
        plot(all_trial_end{i} - all_reward{i}(1),i,'.r');
        plot(all_cue{i}(1) - all_reward{i}(1),i,'.b');
    end
end
title('Lick Raster')

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
lever_surround_sec = 0.3;
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

figure; hold on;
plot(lever_force_reward');
line([round(lever_surround_sample/2) round(lever_surround_sample/2)],ylim,'color','k','linestyle','--');
title('Rewarded lever force')

avg_lever_force_reward = nanmean(lever_force_reward)';
lever_force_reward_norm = bsxfun(@times,lever_force_reward,1./sqrt(sum(lever_force_reward.^2,2)));
lever_force_reward_dot = lever_force_reward_norm*avg_lever_force_reward;

% TO DO: smooth by 50 ms, take L2 norm and dot product for stereotypy



%% Get calcium transients (the dumb way)

%peakstd = [1.5*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
peakstd = [5 5]; % definitions for max trigger, in stds
peak_back = 20;
local_max = 30;
drop_percentage = .8;
spike_frames = {};
spike_amplitudes = {};

num_rois = size(im.roi_trace_df,1);

% this is to catch a dumb bug: if ROI is offscreen, just make noise
for i = 1:num_rois
    if all(im.roi_trace_df(i,:) == im.roi_trace_df(i,1))
        im.roi_trace_df(i,:) = rand(size(im.roi_trace_df(i,:)));
    end
end

% Estimate noise through smoothing
size_smooth = 30;
smooth_type = 'loess';
    % for now - just smooth whole trace

for i = 1:size(im.roi_trace_df,1)
    roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),size_smooth,smooth_type);
end

smooth_std = zeros(size(im.roi_trace_df,1),1);
smooth_std = sqrt(sum((abs(im.roi_trace_df - roi_trace_df_smooth).^2)/size(roi_trace_df_smooth,2),2));

for i = 1:num_rois
    % Estimate noise through smoothing
    curr_trace_std = smooth_std(i);
    minpeak = curr_trace_std*peakstd;
    minpeak = 0.5;
    [spike_frames{i} spike_amplitudes{i}] = ap_schmittDetect(roi_trace_df_smooth(i,:),minpeak,peak_back,local_max,drop_percentage);
    disp([num2str(i) '/' num2str(num_rois)])
end

% make spike matrix
spike_matrix = zeros(size(roi_trace_df));
for i = 1:num_rois
    spike_matrix(i,spike_frames{i}) = 1;
end

% % eliminate peaks at loop junctions
% for curr_roi = 1:num_rois;
%     file_frames = cellfun(@length,roi_trace_long_split);
%     file_frames_start = [1;cumsum(file_frames)+1];
%     loop_frames_start = file_frames_start(1:8:end);
%     for i = 1:length(loop_frames_start);
%         if i == 1
%             spike_train(:,loop_frames_start(i):loop_frames_start(i)+20) = 0;
%         else
%             spike_train(:,loop_frames_start(i)-20:loop_frames_start(i)+20) = 0;
%         end
%     end
% end

% % eliminate peaks at file junctions because they're mistimed
% for i = 1:num_rois
%     if ~isempty(spike_frames{i})
%         file_stitch = [];
%         file_stitch = find(mod(spike_frames{i}-1,movie_frames) == 0);
%         peak_frames(file_stitch,:) = [];
%     end
% end

%%%% what's been kind of working
%[peak_frames] = ap_schmittDetect(roi_trace_long(24,:),minpeak,5,20,0.5);
%%%%

%% Plot stack of all traces
% roi_trace_df_smooth = zeros(size(im.roi_trace_df));
% for curr_cell = 1:size(im.roi_trace_df,1)
%     
%     % normalized cutoff freq is fraction of nyquist
%     butterworth_freq = (0.1*bhv.framerate)/(bhv.framerate/2);
%     [b a] = butter(4, butterworth_freq);
%     % using filtfilt instead of filter allows for zero-phase shift
%     roi_trace_df_smooth(curr_cell,:) = filtfilt(b,a,im.roi_trace_df(curr_cell,:));
% end
% %plot_color = zeros(length(gad_cells),3);
% gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';

curr_traces = [im.roi_trace_df(curr_move(sort_idx),:)]; 
% color_idx = kidx;
% plot_color(color_idx == 1,1) = 0.5;
% plot_color(color_idx == 2,2) = 0.5;
% plot_color(color_idx == 3,3) = 0.5;

figure;
hold on;
spacing = 2;
for i = 1:size(curr_traces,1)
    if exist('plot_color','var')
        plot([1:size(curr_traces,2)], ... %/bhv.framerate, ...
            curr_traces(i,:) + spacing*i,'-k','linewidth',2,'color',plot_color(i,:));
    else
        plot([1:size(curr_traces,2)], ...% /bhv.framerate, ...
            curr_traces(i,:) + spacing*i,'-k','linewidth',2);
    end
end

% % plot cues as magenta lines
% for i = 1:length(bhv.cue_frames)
%     line([bhv.cue_frames(i)/bhv.framerate bhv.cue_frames(i)/bhv.framerate],ylim,'color','m');
% end
% % plot rewards as green lines
% for i = 1:length(bhv.reward_frames)
%     line([bhv.reward_frames(i)/bhv.framerate bhv.reward_frames(i)/bhv.framerate],ylim,'color','g');
% end


%% TEMP


s = 37000; % starting time
l = 1200; % length of time
a = bhv.cue_frames(bhv.cue_frames > s & bhv.cue_frames < s+l);
b = bhv.reward_frames(bhv.reward_frames > s & bhv.reward_frames < s+l);
c = bhv.lick_frames(bhv.lick_frames > s & bhv.lick_frames < s+l);

figure;

p1 = subplot(8,1,1); hold on;
for i = 1:length(c);
    line([c(i)-s c(i)-s],[0 1],'color','k')  
end


p2 = subplot(8,1,2:8); hold on;
for i = 1:length(a)
    patch([a(i)-s a(i)-s b(i)-s b(i)-s],[0 3 3 0],[1 0 0],'EdgeColor','none');
end
plot(bhv.lever_force_resample(s:s+l),'k','linewidth',2)

linkaxes([p1 p2],'x');


%% Plot trace / behavior / cell
% requires: roi_trace_df, last cell to be run

curr_cell = 103;


if ~exist('im','var')
    im.roi_trace_df = roi_trace_df; 
end

figure;hold on;
plot([1:size(im.roi_trace_df,2)],im.roi_trace_df(curr_cell,:),'k');
%plot([1:size(im.roi_trace_df,2)],im.oopsi(curr_cell,:),'r');
plot([1:length(bhv.lever_force_resample)],bhv.lever_force_resample,'r');

catch_trials = zeros(size(bhv.cue_frames));
if isfield(bhv,'catch_trials')
    catch_trials = bhv.catch_trials(bhv.imaged_trials);
end

% if you want to plot the raw lever_force, but not really necessary
%plot([1:length(bhv.lever_force)]/(length(bhv.lever_force)/(size(im.roi_trace_df,2))),bhv.lever_force,'b');
if exist('roi_trace_df_smooth')
    plot([1:size(roi_trace_df_smooth,2)],roi_trace_df_smooth(curr_cell,:),'--b','LineWidth',2);
end

if exist('peak_frames','var') && ~isempty(peak_frames{curr_cell})
    plot(peak_frames{curr_cell},0,'.m','MarkerSize',15);
end

if isfield(im,'roi_concat_oopsi');
    plot(im.roi_concat_oopsi(curr_cell,:),'r');
end

if exist('peak_frames')
     plot(peak_frames{curr_cell}/bhv.framerate,0,'.k','MarkerSize',30);
end
 
% plot cues as magenta lines
for i = 1:length(bhv.cue_frames)
    temph = line([bhv.cue_frames(i) bhv.cue_frames(i)],ylim,'color','m');
    if catch_trials(i)
       set(temph,'color','b') 
    end
end
% plot rewards as green lines
for i = 1:length(bhv.reward_frames)
    line([bhv.reward_frames(i) bhv.reward_frames(i)],ylim,'color','g');
end
% plot licks as black dots
plot(bhv.lick_frames,2,'.r');
% plot lick rate
%somehow...

% show cue/reward aligned activity
if exist('reward_aligned_active','var')
    figure; hold on;
    plot(bhv.lever_force_resample,'m');
    plot(bhv.lick_frames,2,'.k');
    plot(im.roi_trace_df(curr_cell,:),'k')
    plot(im.roi_trace_df_smooth_cutoff(curr_cell,:),'b','linewidth',2)
    for i = 1:length(cue_frames_good)
        line([cue_frames_good(i) cue_frames_good(i)],ylim,'color',[0.5 0 0]);
    end
    for i = 1:length(reward_frames_good)
        line([reward_frames_good(i) reward_frames_good(i)],ylim,'color',[1 0 0]);
    end
    
    for i = find(cue_aligned_active(curr_cell,:))
        line([cue_frames_good(i) cue_frames_good(i)],ylim,'color',[0 0.5 0]);
    end
    for i = find(reward_aligned_active(curr_cell,:))
        line([reward_frames_good(i) reward_frames_good(i)],ylim,'color',[0 1 0]);
    end
end

% plot around reward
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
figure;set(gcf,'name',num2str(curr_cell));imagesc(roi_reward);
colormap(gray);
%line([frames_surround+1 frames_surround+1],ylim,'LineWidth',2,'color','w',...
%    'LineStyle','--');
title('Reward-aligned');

% plot around cue
cue_reward = nan(length(bhv.cue_frames),frames_surround*2+1);
for i = 1:length(bhv.cue_frames);
    curr_reward = round(bhv.cue_frames(i));
    roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
    if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
        continue
    end
    cue_reward(i,:) = im.roi_trace_df(curr_cell,roi_reward_frames);
end
figure;imagesc(cue_reward);
colormap(gray);
%line([frames_surround+1 frames_surround+1],ylim,'LineWidth',2,'color','w',...
%    'LineStyle','--');
title('Cue-aligned');


% if available, plot around movement starts
if exist('move_starts','var');
    move_align = nan(length(move_starts),frames_surround*2+1);
    for i = 1:length(move_starts);
        curr_reward = round(move_starts(i));
        roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
            continue
        end
        move_align(i,:) = im.roi_trace_df(curr_cell,roi_reward_frames);
    end
    figure;imagesc(move_align);
    colormap(gray);
    %line([frames_surround+1 frames_surround+1],ylim,'LineWidth',2,'color','w',...
    %    'LineStyle','--');
    title('Move-start aligned');
end


% plot mean of cue-aligned and reward_aligned with std
figure; hold on
roi_reward_mean = nanmean(roi_reward);
roi_reward_sem = nanstd(roi_reward)./sqrt(size(roi_reward,1));
cue_reward_mean = nanmean(cue_reward);
cue_reward_sem = nanstd(cue_reward)./sqrt(size(cue_reward,1));

%plot(roi_reward_mean,'k');
plot(cue_reward_mean,'k');
%jbfill([1:length(roi_reward_mean)],roi_reward_mean+roi_reward_sem, ...
%    roi_reward_mean-roi_reward_sem,[0.5 0.5 0.5]);
jbfill([1:length(cue_reward_mean)],cue_reward_mean+cue_reward_sem, ...
    cue_reward_mean-cue_reward_sem,[0.5 0.5 0.5]);
%legend({'Reward', 'Cue'});
title(curr_cell)
figure;
plot(roi_reward_mean,'k');
%plot(cue_reward_mean,'k');
jbfill([1:length(roi_reward_mean)],roi_reward_mean+roi_reward_sem, ...
    roi_reward_mean-roi_reward_sem,[0.5 0.5 0.5]);
%jbfill([1:length(cue_reward_mean)],cue_reward_mean+cue_reward_sem, ...
%    cue_reward_mean-cue_reward_sem,[0.5 0.5 0.5]);
%legend({'Reward', 'Cue'});
title(curr_cell)

%% Show ROIs / color

% define the cells to look at
cells = [1:length(roi_labels)];

if exist('roi_labels','var')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels);
    pyr_cells = ~gad_cells;
    
    filled = cellfun(@(x) any(strcmp('filled',x)),roi_labels);
    
else
    pyr_cells = cells;
    gad_cells = [];
end

pyr_unfilled_idx = find(pyr_cells & ~filled);

figure;

cell_color = ones(length(pyr_unfilled_idx),3)*0.5;

% a = classified_cells.move_cells_peak{7};
% b = classified_cells.move_cells_oopsi{7};
% cell_color(:,1) = a;
% cell_color(:,2) = b;

roi_handle = [];
for i = 1:length(pyr_unfilled_idx)
    curr_cell = pyr_unfilled_idx(i);
%for i = setdiff(pyr_cells,filled);
    roi_handle(i) = patch(polygon.ROI{curr_cell}(:,1),-polygon.ROI{curr_cell}(:,2), ...
        cell_color(i,:));
        %[max_idx(i)/max(max_idx) max_peaks(i)/max(max_peaks) 0]);
end
ylim([-512 0]);
xlim([0 512]);


% cell_group = {};
% cell_group{1} = [a]; cell_color{1} = [1 0 0];
% cell_group{2} = [b]; cell_color{2} = [0 1 0];
% cell_group{2} = [intersect(a,b)]; cell_color{2} = [1 1 0];
% 
% 
% if ~isempty(cell_group)
%     for curr_group = 1:length(cell_group)
%         for i = 1:length(cell_group{curr_group})
%             set(roi_handle(cell_group{curr_group}(i)),'FaceColor', ...
%                 cell_color{curr_group})
%         end
%     end
% end



% % make legend
% legend_handles = [];
% for i = 1:length(cell_group)
%    legend_handles(i) = patch([0 1],[0 1],cell_color{i});
% end
% legend([legend_handles(:)],{'cue' 'reward' 'movestart' 'movestop'});


% display ROI numbers
% write numbers in the centers for number of spikes
% for i = 1:length(polygon.ROI);
%     roi_numbers_h(i) = text(roi_center_x(i),-roi_center_y(i),num2str(i),'color','r');
% end

%title(trial);
%print(gcf,['AP71_120619_PyrGadActivity_' num2str(trial) '.png'],'-dpng')

%     paperUnits = get(gcf,'PaperUnits');
%     set(gcf,'PaperUnits','inches');
%     paperSize = get(gcf,'PaperSize');
%     paperPosition = [0.5 0.5 paperSize - 0.5];
%     set(gcf,'PaperPosition',paperPosition);
%     set(gcf,'PaperUnits',paperUnits);
%     print(gcf,['AP71_120619_PyrGadActivity_' num2str(page) '.png'],'-dpng');
%     close all



%% Plot all ROIs with lever and cue/presses

figure;
a = subplot(10,1,1:3); hold on;
plot(bhv.lever_force_resample,'k')
% plot cues as magenta lines
for i = 1:length(bhv.cue_frames)
    line([bhv.cue_frames(i) bhv.cue_frames(i)],ylim,'color','m');
end
% plot rewards as green lines
for i = 1:length(bhv.reward_frames)
    line([bhv.reward_frames(i) bhv.reward_frames(i)],ylim,'color','g');
end
% plot licks as black dots
plot(bhv.lick_frames,2,'.k');

b = subplot(10,1,4:10);
imagesc(roi_trace_df)
% plot cues as magenta lines
for i = 1:length(bhv.cue_frames)
    line([bhv.cue_frames(i) bhv.cue_frames(i)],ylim,'color','m');
end
% plot rewards as green lines
for i = 1:length(bhv.reward_frames)
    line([bhv.reward_frames(i) bhv.reward_frames(i)],ylim,'color','g');
end

linkaxes([a b],'x');

%% Make spike triggered avg

% im.oopsi(im.oopsi < 0.2) = 0;
% 
% % Try using peaks of deconv trace
% spike_frames = cell(size(roi_trace_df,1),1);
% for curr_cell = 1:size(roi_trace_df,1)
%     [spike_frames{curr_cell} mintab] = peakdet(im.oopsi(curr_cell,:),0.3);    
% end

% smooth_size = 30;
% roi_trace_df_smooth = nan(size(im.roi_trace_df));
% for i = 1:size(im.roi_trace_df,1);
%     roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess');
% end

caEvent_matrix = AP_caEvents(im.roi_trace_df,pyr_unfilled,[]);
event_onsets = [false(size(caEvent_matrix,1),1) ...
    diff(caEvent_matrix>0,[],2) == 1];
pyr_activity_onset = zeros(size(caEvent_matrix));
pyr_activity_onset(event_onsets) = caEvent_matrix(event_onsets);
spike_frames = arrayfun(@(x) find(pyr_activity_onset(x,:)), ...
    1:size(pyr_activity_onset,1),'uni',false);

time_back = 10000; %ms
time_forward = 10000; %ms

xsg_sample_rate = 10000;

samples_back = round((time_back/1000)*(xsg_sample_rate));
samples_forward = round((time_forward/1000)*(xsg_sample_rate));

frames_back = round((time_back/1000)*bhv.framerate);
frames_forward = round((time_forward/1000)*bhv.framerate);

lever_spike_all = cell(size(im.roi_trace_df,1),1);
lever_diff_spike_all = cell(size(im.roi_trace_df,1),1);
trace_spike_all = cell(size(im.roi_trace_df,1),1);

lever_force_smooth = smooth(bhv.lever_force,500);
lever_force_diff = [0;diff(lever_force_smooth)];

for curr_cell = 1:size(im.roi_trace_df,1);
    lever_spike = nan(length(spike_frames{curr_cell}),samples_back+samples_forward+1);
    lever_diff_spike = nan(length(spike_frames{curr_cell}),samples_back+samples_forward+1);
    trace_spike = nan(length(spike_frames{curr_cell}),frames_back+frames_forward+1);
    for i = 1:length(spike_frames{curr_cell})
        
        curr_spike = spike_frames{curr_cell}(i);
        curr_spike_sample = round((curr_spike/bhv.framerate)*xsg_sample_rate);
        % only use if within bounds
        if curr_spike_sample-samples_back > 0 && ...
                curr_spike_sample + samples_forward < length(bhv.lever_force)
            curr_force_diff = lever_force_diff(curr_spike_sample-samples_back...
                :curr_spike_sample+samples_forward);
            curr_force = bhv.lever_force(curr_spike_sample-samples_back...
                :curr_spike_sample+samples_forward);
            %surround_force(i,:) = (curr_force*(n_best(curr_spike)/sum(n_best(spike_indx))))';
            if ~all(isnan(curr_force))
                lever_spike(i,:) = curr_force;
                lever_diff_spike(i,:) = curr_force_diff;%/spike_amplitudes{curr_cell}(i);
                trace_spike(i,:) = im.roi_trace_df(curr_cell,curr_spike-frames_back: ...
                    curr_spike + frames_forward);
            end
        end
    end
    lever_spike_all{curr_cell} = lever_spike;
    lever_diff_spike_all{curr_cell} = lever_diff_spike;
    trace_spike_all{curr_cell} = trace_spike;
    curr_cell
end
lever_diff_spike_all_mean = cellfun(@(x) nanmean(x,1)',lever_diff_spike_all,'UniformOutput',false);
lever_diff_spike_all_mean = [lever_diff_spike_all_mean{:}]';
lever_diff_spike_all_mean_norm = bsxfun(@times,lever_diff_spike_all_mean,1./sqrt(sum(lever_diff_spike_all_mean.^2,2)));

lever_spike_all_mean = cellfun(@(x) nanmean(x,1)',lever_spike_all,'UniformOutput',false);
lever_spike_all_mean = [lever_spike_all_mean{:}];
trace_spike_all_mean = cellfun(@(x) nanmean(x,1)',trace_spike_all,'UniformOutput',false);
trace_spike_all_mean = [trace_spike_all_mean{:}]';


% get gaussian-assumed confidence intervals 
%%%% NOT STARTED YET
% lever_diff_spike_all_ci = cell(size(roi_trace_df,1),1);
% for curr_cell = 1:size(roi_trace_df,1);
%     lever_diff_spike_all_ci{curr_cell} = nan(2,size(lever_diff_spike_all{curr_cell},2));
%     if size(lever_diff_spike_all{curr_cell},1) > 2;
%         for i = 1:size(lever_diff_spike_all{curr_cell},2)
%             lever_diff_spike_all_ci{curr_cell}(:,i) = bootci(100,@nanmean,lever_diff_spike_all{curr_cell}(:,i));
%         end
%     end
%     disp(['Finished conf intervals for cell ' num2str(curr_cell)])
% end

% % get bootstrap confidence intervals
% num_rep = 100;
% lever_diff_spike_all_ci = cell(size(roi_trace_df,1),1);
% for curr_cell = 1:size(roi_trace_df,1);
%     lever_diff_spike_all_ci{curr_cell} = nan(2,size(lever_diff_spike_all{curr_cell},2));
%     if size(lever_diff_spike_all{curr_cell},1) > 2;
%         for i = 1:size(lever_diff_spike_all{curr_cell},2)
%             lever_diff_spike_all_ci{curr_cell}(:,i) = bootci(100,@nanmean,lever_diff_spike_all{curr_cell}(:,i));
%         end
%     end
%     disp(['Finished conf intervals for cell ' num2str(curr_cell)])
% end


% figure;
% imagesc(lever_spike);
% figure;
% imagesc(trace_spike);
% figure; hold on
% plot(nanmean(lever_spike));
% plot(nanmean(trace_spike),'r');

%% Make behavior-triggered df/deconv average

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

% get indicies for which cues are overlapped by movement or not
cue_overlap_idx = cellfun(@(x) 5+find(ismember(round(cue_frames_good-5),x)), ...
    movement_frames,'UniformOutput',false);
reward_overlap_idx = cellfun(@(x) find(ismember(round(reward_frames_good),x)), ...
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
max_xcorr = reshape(max(xcorr([movements_p_pad{:}])),length(move_time),...
    length(move_time));
max_xcorr_norm = reshape(max(xcorr([movements_p_pad{:}],'coeff')),length(move_time),...
    length(move_time));

% find sum deconv spike for each movement/cell
movement_sum_oopsi = zeros(size(roi_trace_df,1),length(movement_frames));
for i = 1:length(movement_frames)
    movement_sum_oopsi(:,i) = sum(im.oopsi(:,movement_frames{i}),2);
end

for spike_loop = 1:3
    
    switch spike_loop
        case 1
            spike_frames = round(cue_frames_good);
        case 2
            spike_frames = move_start_frames;
        case 3
            spike_frames = round(reward_frames_good);
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
    spike_frames_idx(any(spike_frames_idx < 0,2) | ...
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

% plot the first measurement, over all bins
curr_cell = 7;
cr_bin_mean = grpstats(move_params_act{1}{2}(curr_cell,:), ...
    move_params_act{1}{3},'mean');
iti_bin_mean = grpstats(move_params_act{1}{4}(curr_cell,:), ...
    move_params_act{1}{5},'mean');

cr_bin_std = grpstats(move_params_act{1}{2}(curr_cell,:), ...
    move_params_act{1}{3},'std');
iti_bin_std = grpstats(move_params_act{1}{4}(curr_cell,:), ...
    move_params_act{1}{5},'std');
figure;hold on;
errorbar([1:length(cr_bin_mean)],cr_bin_mean, ...
    cr_bin_std,'r');
errorbar([1:length(iti_bin_mean)],iti_bin_mean,...
    iti_bin_std,'k');

%% Align spike-triggered averages
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
            % get lag from activity to movement
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

% get bootstrap confidence intervals (using only 1/100 of the force trace)
num_rep = 100;
lever_diff_spike_all_aligned_ci = cell(size(roi_trace_df,1),1);
for curr_cell = 1:size(roi_trace_df,1);
    lever_diff_spike_all_ci{curr_cell} = nan(2,size(lever_diff_spike_all{curr_cell},2));
    if size(lever_diff_spike_all{curr_cell},1) > 2;
        for i = 1:100:size(lever_diff_spike_all{curr_cell},2)
            lever_diff_spike_all_aligned_ci{curr_cell}(:,i) = ...
                bootci(100,@nanmedian,lever_diff_spike_all_aligned{curr_cell}(:,i));
        end
    end
    disp(['Finished bootstrap conf intervals for cell ' num2str(curr_cell)])
end

% % get modulated cells
% % Do this
% figure; hold on;
% plot(lever_diff_spike_all_aligned_ci{29}(:,1:100:end)','-r')
% indx_lo = lever_diff_spike_all_aligned_ci{29}(1,1:100:end)' > 0;
% indx_hi = lever_diff_spike_all_aligned_ci{29}(2,1:100:end)' < 0;
% temp = lever_diff_spike_all_aligned_ci{29}(:,1:100:end);
% plot(find(indx_lo|indx_hi),temp(:,indx_lo|indx_hi)','.g')




%% Get lag between spikes and smoothed convolution of lever velocity and mean
% nope doesn't really work/need to be done - was already done in alignment!

movement_lag = nan(size(roi_trace_df,1),1);
for i = 1:size(roi_trace_df,1);
    % convolute aligned mean veloctiy with all velocity
    temp_conv = [];
    temp_conv = conv(lever_force_diff,nanmedian(lever_diff_spike_all_aligned{i}),'same');
    % resample to get length lever force = frames
    [n d] = rat(size(roi_trace_df,2)/length(lever_force));
    temp_conv_resample = resample(temp_conv,n,d);
    % threshold to get only positives
    temp_conv_resample_thresh = temp_conv_resample;
    temp_conv_resample_thresh(temp_conv_resample_thresh < 0) = 0;
    % smooth (value by looking at one cell by example)
    temp_conv_resample_thresh_smooth = smooth(temp_conv_resample_thresh,1);
    
    % get max covariance between spikes and smoothed conv, limit 100 frames
    % use the shifted spikes from initial alignment
    temp_spike_train = zeros(1,size(roi_trace_df,2));
    temp_spike_train(spike_shift{i}) = 1;
    
    temp_xcov = [];
    temp_xcov = xcov(temp_spike_train,temp_conv_resample_thresh_smooth,200);
    maxlag = find(temp_xcov == max(temp_xcov));
    if ~isempty(maxlag)
        movement_lag(i) = 101 - maxlag; 
    end
    disp(['Got movement lag: cell ' num2str(i)])
end

%% Make reward-triggered average
% plot around reward
frames_surround = 200;
roi_reward_cell = cell(size(roi_trace_df,1),1);
for curr_cell = 1:length(roi_reward_cell)
    roi_reward_cell{curr_cell} = nan(length(reward_frames),frames_surround*2+1);
    for i = 1:length(reward_frames);
        curr_reward = round(reward_frames(i));
        roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
            continue
        end
        roi_reward_cell{curr_cell}(i,:) = roi_trace_df(curr_cell,roi_reward_frames);
    end
end
roi_reward_all = cellfun(@(x) nanmean(x)',roi_reward_cell,'UniformOutput',false);
roi_reward_all = [roi_reward_all{:}]';

%% PCA of the lever presses
% requires all of the above
lever_spike_nonan = lever_spike;
lever_spike_nonan(all(isnan(lever_spike),2),:) = [];
[coeff,score,latent] = princomp(lever_spike_nonan);


%% Regress lever trace from cell traces
% TO DO: multiple regression: past (/future?) time points from lever and
% traces, lick rate, threshold by standard deviation, make matrix of
% constants instead of guessing the 1.3

% smooth, if that helps

% try a bunch of smoothing
for i = 5
    smooth_factor = i;
    
    smooth_kernel = ones(1,smooth_factor)/smooth_factor;
    roi_trace_df_smooth = [];
    roi_trace_df_smooth = filter2(smooth_kernel,roi_trace_df);
    roi_trace_df_decimated = [];
    roi_trace_df_decimated = roi_trace_df_smooth(:,1:smooth_factor:end);
    
    lever_force_resample_decimated = [];
    lever_force_resample_decimated = smooth(lever_force_resample,smooth_factor);
    lever_force_resample_decimated = lever_force_resample_decimated(1:smooth_factor:end)';
    lever_force_resample_decimated_diff = [diff(lever_force_resample_decimated) 0];
    
    % this is to test the shifting
    %shift_values = -round(framerate)*10:round(framerate/2):round(framerate)*10;
    shift_values = 0;
    n_bases = 1;
    
    % weight_thresh = cell(length(shift_values),size(roi_trace_df_decimated,1)+1);
    % weight_thresh_diff = cell(length(shift_values),size(roi_trace_df_decimated,1)+1);
    %
    % predicted_traces = cell(length(shift_values),size(roi_trace_df_decimated,2));
    % predicted_traces_diff = cell(length(shift_values),size(roi_trace_df_decimated,2));
    %
    % predicted_corrcoef = cell(length(shift_values),1);
    % predicted_diff_corrcoef = cell(length(shift_values),1);
    
    weight_thresh = cell(shift_values,1);
    weight_thresh_diff = cell(shift_values,1);
    
    predicted_traces = cell(shift_values,1);
    predicted_traces_diff = cell(shift_values,1);
    
    predicted_corrcoef = cell(length(shift_values),1);
    predicted_diff_corrcoef = cell(length(shift_values),1);
    
    
    % set up basis sets: random orderings of indicies
    basis_indx = zeros(n_bases,size(roi_trace_df_decimated,2));
    for i = 1:n_bases
        basis_indx(i,:) = randperm(size(roi_trace_df_decimated,2),size(roi_trace_df_decimated,2));
    end
    
    for i = 1:length(shift_values);
        %     temp_thresh = [];
        %     temp_thresh = roi_trace_df_decimated;
        %     temp_thresh(temp_thresh < thresh_values(i)) = 0;
        
        % threshold based on std and minumum
        roi_trace_thresh = [];
        roi_trace_thresh = roi_trace_df_decimated;
        for curr_roi = 1:size(roi_trace_df_decimated,1);
            trace_median = median(roi_trace_df_decimated(curr_roi,:));
            trace_std = std(roi_trace_df_decimated(curr_roi,:));
            roi_trace_thresh(curr_roi,:) = roi_trace_thresh(curr_roi,:) - trace_median;
            
            trace_min = min(roi_trace_thresh(curr_roi,:));
            trace_min_std = abs(trace_min/trace_std)*2;
            
            thresh = trace_min_std*trace_std;
            % subtract the bad new estimate for baseline
            
            roi_trace_thresh(curr_roi,roi_trace_thresh(curr_roi,:)<thresh) = 0;
        end
        
        % TESTING
        roi_trace_thresh = circshift(roi_trace_thresh,[0 shift_values(i)]);
        % END TESTING
        
        roi_trace_thresh_constantoffset = [roi_trace_thresh;ones(1,size(roi_trace_thresh,2))];
        
        % fit 80%, predict 20%, n_bases times
        for basis = 1:n_bases
            for frac = 1:5 % because 80/20
                numframes = size(basis_indx,2);
                basis_frames_indx = true(1,numframes);
                basis_frames_indx(1:ceil(numframes*0.2)) = false;
                basis_frames_indx = circshift(basis_frames_indx, ...
                    [0 ceil(numframes*0.2)*frac]);
                basis_frames = basis_indx(basis,basis_frames_indx);
                predicted_frames = ~basis_frames_indx;
                
                warning off % it's rank deficient, don't get it
                weight_thresh{i}(basis,:) = lever_force_resample_decimated(basis_frames)/...
                    roi_trace_thresh_constantoffset(:,basis_frames);
                weight_thresh_diff{i}(basis,:) = lever_force_resample_decimated_diff(basis_frames)/ ...
                    roi_trace_thresh_constantoffset(:,basis_frames);
                warning on
                
                predicted_traces{i}(basis,predicted_frames) = ...
                    weight_thresh{i}(basis,:)*roi_trace_thresh_constantoffset(:,predicted_frames);
                predicted_traces_diff{i}(basis,predicted_frames) = ...
                    weight_thresh_diff{i}(basis,:)*roi_trace_thresh_constantoffset(:,predicted_frames);
            end
            
            temp_corrcoef = corrcoef(lever_force_resample_decimated,predicted_traces{i}(basis,:));
            temp_diff_corrcoef = corrcoef(lever_force_resample_decimated_diff,predicted_traces_diff{i}(basis,:));
            
            predicted_corrcoef{i}(basis) = temp_corrcoef(2);
            predicted_diff_corrcoef{i}(basis) = temp_diff_corrcoef(2);
            
            disp([num2str(100*basis/n_bases) '% bases sets,' num2str(100*i/length(shift_values)) '% shift values']);
        end
    end
    
    %%%% WARNING! This next code is dumb. could do the same with regular
    %%%% corrcoef, since weighting doesn't change anything. whoops.
    
    % get correlations given by single cells
    cell_corr = cell(length(shift_values),1);
    for i = 1:length(shift_values)
        for curr_cell = 1:size(roi_trace_df,1);
            for basis = 1:n_bases;
                temp_pred = [];
                for frac = 1:5 % because 80/20
                    numframes = size(basis_indx,2);
                    basis_frames_indx = true(1,numframes);
                    basis_frames_indx(1:ceil(numframes*0.2)) = false;
                    basis_frames_indx = circshift(basis_frames_indx, ...
                        [0 ceil(numframes*0.2)*frac]);
                    basis_frames = basis_indx(basis,basis_frames_indx);
                    predicted_frames = ~basis_frames_indx;
                    
                    temp_pred(predicted_frames) = ...
                        weight_thresh{i}(basis,curr_cell)*...
                        roi_trace_thresh_constantoffset(curr_cell,predicted_frames);
                end
                temp_corr = [];
                temp_corr = corrcoef(lever_force_resample_decimated,temp_pred);
                cell_corr{i}(curr_cell,basis) = temp_corr(2);
            end
            disp([num2str(100*i/length(shift_values)) '% shift values, cell ' num2str(curr_cell)])
        end
    end
    
    %%%% end dumb code
    %%%% end dumb code
    
    % get correlations between pieces of cells which are active and trace
    cell_thresh_corr = zeros(size(roi_trace_thresh,1),1);
    for i = 1:size(roi_trace_thresh,1)
        active_indx = roi_trace_thresh(i,:) ~= 0;
        temp_corr = [];
        temp_corr = corrcoef(roi_trace_thresh(i,active_indx), ...
            lever_force_resample_decimated(active_indx));
        if sum(active_indx) > 2
            cell_thresh_corr(i) = temp_corr(2)* ...
                (sum(active_indx)/length(lever_force_resample_decimated));
        end
        disp(num2str(i/size(roi_trace_thresh,1)));
    end
    
    
    % prediction_error = [];
    % prediction_error = abs(repmat(lever_force_resample_decimated,...
    %     length(thresh_values),1) - predicted_traces).^2;
    %
    % figure; hold on;
    % plot(lever_force_resample_decimated,'k')
    % plot(predicted_traces,'--r')
    % sq_error = sum(abs(lever_force_resample_decimated - predicted_traces).^2)
    
    % plot(mean(prediction_error,2))
    % set(gca,'XTick',[1:length(thresh_values)])
    % set(gca,'XTickLabel',thresh_values)
    
    
end

%% Threshold trace like Tank lab ish
% just for checking it out, don't do any corrections, just raw med/std

roi_trace_thresh = roi_trace_df;
for i = 1:size(roi_trace_df,1);
    trace_median = median(roi_trace_df(i,:));
    trace_std = std(roi_trace_df(i,:));
    trace_min = min(roi_trace_df(i,:));
    trace_min_std = abs(trace_min/trace_std)+0.5;
    
    thresh = trace_min_std*trace_std-trace_median;
    % subtract the bad new estimate for baseline
    roi_trace_thresh(i,:) = roi_trace_thresh(i,:) - trace_median;
    roi_trace_thresh(i,roi_trace_thresh(i,:)<thresh) = 0;
end


%% Trying Hatsopoulos movement vector preference

% % make long df trace with first 10 loop frames as 0
% % identify the bad frames (first 10 of every loop)
% badtrace_mask = roi_trace_long_split;
% badtrace_mask = cellfun(@(x) false(1,length(x)), badtrace_mask,'UniformOutput',false);
% badtrace_mask(1:8:end) = cellfun(@(x) [true(1,10) ...
%     false(1,length(x)-10)], badtrace_mask(1:8:end),'UniformOutput',false);
% badtrace_mask = [badtrace_mask{:}];
%
% roi_trace_df(badtrace_mask) = 0;

% do oopsi to guess at spikeish regions
n_best = {};
P_best = {};
V = {};
V.dt = 1/bhv.framerate;
for i = 1:size(roi_trace_df,1)
    [n_best{i} P_best{i} V_est]=fast_oopsi(roi_trace_df(i,:),V);
    disp(num2str(i/size(roi_trace_df,1)));
end
n_best = [n_best{:}]';

% normalize n_best
n_best = n_best/max(n_best);


% at the moment, downsampling by 10

segment_back = 100; % in ms
segment_forward = 300; % in ms
segment_separate_frames = 2; % in frames

segment_separate_samples = round((segment_separate_frames/framerate)...
    *xsg_sample_rate);
segment_back_samples = round(segment_back*(xsg_sample_rate/1000));
segment_forward_samples = round(segment_forward*(xsg_sample_rate/1000));

segment_back_frames = round((segment_back/1000)*framerate);
segment_forward_frames = round((segment_forward/1000)*framerate);

% smooth lever force by 50ms
lever_force_smooth = smooth(lever_force,500);
% take out segment_length fragments every segment_separate ms
segment_centers = segment_back_samples+1:segment_separate_samples:...
    length(lever_force_smooth)-(segment_forward_samples+1);
lever_force_segments = zeros(length(segment_centers),...
    segment_back_samples+segment_forward_samples+1);
for j = 1:length(segment_centers)
    i = segment_centers(j);
    lever_force_segments(j,:) = lever_force_smooth(i-segment_back_samples:...
        i+segment_forward_samples)';
    j/length(segment_centers)
end

% get velocity vector
segment_diff = diff(lever_force_segments,1,2);

% L2 normalize velocity vector
segment_diff_norm = bsxfun(@times,segment_diff,1./sqrt(sum(segment_diff.^2,2)));

% get frame start/end times for each segment
segment_frames = zeros(length(segment_centers),2);
segment_frames(:,1) = round(((segment_centers/(xsg_sample_rate)))*framerate)-1;
segment_frames(:,2) = round(((segment_centers/(xsg_sample_rate)))*framerate)+1;


% % check for significant values during segments
% roi_trace_thresh_bin = roi_trace_thresh;
% roi_trace_thresh_bin = find(roi_trace_thresh_bin > 0);
% segment_sig = zeros(size(segment_frames,1),1);
% for i = 1:size(segment_frames,1);
%    if  any(ismember(roi_trace_thresh_bin,...
%            [segment_frames(i,1):segment_frames(i,2)]))
%        segment_sig(i) = 1;
%    end
% end

% get rid of the first for everything
% DO THIS HERE?%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%5

% do PCA on segment velocities
[coeff score latent] = princomp(segment_diff_norm');
PCs_all = segment_diff_norm'*coeff(:,1:10);
PCs = PCs_all;

% try mean of spiking
spike_sum = zeros(length(segment_centers),1);
for i = 1:length(segment_centers)
    spike_sum(i) = mean(n_best(segment_frames(i,1):segment_frames(i,2)));
end

% trying right on center?
spike_sum = zeros(length(segment_centers),1);
for i = 1:length(segment_centers)
    spike_sum(i) = n_best(ceil((segment_centers(i)/xsg_sample_rate)*framerate-2));
end

norm_factor = segment_diff_norm*segment_diff_norm';
norm_factor = diag(norm_factor);

starting = ones(11,1);
starting(2:end) = starting(2:end)*(1/10000);
% pathlet_coeffs = fminsearch(@(x) ...fi
%     sum((spike_sum - exp(x(1)+segment_diff_norm*(PCs*x(2:end)))).^2), starting);
pathlet_coeffs = fminsearch(@(x) ...
    sum((spike_sum - (x(1)*segment_diff_norm*(PCs*x(2:end)))).^2), starting);

% test by circshifting spikes
pathlet_coeffs = fminsearch(@(x) ...
    sum((circshift(spike_sum,[-100 0]) - (x(1)+segment_diff_norm*(PCs*x(2:end)))).^2), starting);

%% Deconvolve spikes based on oopsi
% single exponential up and down

starting = [0.1 2 10];

event_params = fminsearch(@(x) AP_deconvSpike(x,n_best,roi_trace_df), starting);

% this info, tau_down at least, is under P.gam, where tau = dt/(1-gam)

%% Get weighted average of surrounding lever press
% requires oopsi

% get nonzero values of spike train
time_back = 2000; %ms
time_forward = 2000; %ms

samples_back = round((time_back/1000)*(xsg_sample_rate));
samples_forward = round((time_forward/1000)*(xsg_sample_rate));
weighted_force = zeros(samples_back+samples_forward+1,size(roi_trace_df,1));


lever_force_smooth = smooth(lever_force,500);
lever_force_diff = [0;diff(lever_force_smooth)];

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
        %surround_force(i,:) = (curr_force*(n_best(curr_spike)/sum(n_best(spike_indx))))';
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



%% Fast/SMC Oopsi

tic
roi_df_oopsi = zeros(size(im.roi_trace_df));

inactive_thresh = 0.1; % threshold for stitching
num_loops = length(im.roi_trace_split)/8;
file_lengths = cellfun(@length,im.roi_trace_split);
file_lengths_loops = reshape(file_lengths,8,num_loops);
loop_lengths = sum(file_lengths_loops);

% Break up baseline corrected trace into loops
roi_trace_df_split = mat2cell(im.roi_trace_df,[size(im.roi_trace_df,1)], ...
    loop_lengths);

split_indx = mat2cell(1:size(im.roi_trace_df,2),1, ...
    loop_lengths);

for curr_cell = 1:size(im.roi_trace_df,1)
    
    % Stitch together loops at inactive parts
    % Skip the first 10 frames in case of weird things
    loop_start_stitches = cellfun(@(x) find(x(curr_cell,11:end) ...
            < inactive_thresh & x(curr_cell,11:end) > 0,1)+10, ...
            roi_trace_df_split,'UniformOutput',false);
    loop_stop_stitches = cellfun(@(x) find(x(curr_cell,:) ...
        < inactive_thresh & x(curr_cell,:) > 0,1,'last'), ...
        roi_trace_df_split,'UniformOutput',false);
    F_stitched = cellfun(@(x,y,z) x(curr_cell,y:z),roi_trace_df_split,...
        loop_start_stitches,loop_stop_stitches,'UniformOutput',false);
    stitching_indx = cellfun(@(x,y,z) x(y:z), split_indx, ...
        loop_start_stitches, loop_stop_stitches,'UniformOutput',false);
    
    F = [F_stitched{:}];
    
    % set options/parameters
    clear V
    V.dt = 1/bhv.framerate;
    V.est_sig = 1;
    V.est_lam = 1;
    V.est_gam = 1;
    V.est_b = 1;
    V.est_a = 1;
    %V.fast_thr = 1;
    
    clear P
 
    %F = F - min(F);
    tau = 1.5; %2.8;
    P.gam = 1-(V.dt/tau);
    P.lam = 0.1; 
%     P.a = 1;
%     P.b = 0;
    
    [n_out P_out V_out] = fast_oopsi(F,V,P);
    temp = n_out/max(n_out);
    roi_df_oopsi(curr_cell,[stitching_indx{:}]) = temp;
    
    disp(['Finished oopsi: ' num2str(curr_cell/size(roi_trace_df,1))]);
end
toc


%% STA (test/development)

cell_num = 32;

list = [120330 120331 120401 120402 120403 120409 120410 120412 ...
    120413 120414 120415];
list = num2str(list');

figure; hold on;
for i = 1:size(list,1)
    cd(list(i,:))
    clearvars -except list i cell_num
    load('STA_analysis','lever_diff_spike_all_aligned')
    plot(nanmedian(lever_diff_spike_all_aligned{cell_num}),'color',[0 i/16 0]);
    cd ..
end

%% Find movement-correlated cells

lever_force_resample = bhv.lever_force_resample;
framerate = bhv.framerate;

plot_movecorr = 0;

% smooth trace
im.roi_trace_df_smooth = zeros(size(im.roi_trace_df));
size_smooth = 10;
smooth_type = 'loess';
% for now - just smooth whole trace
for i = 1:size(im.roi_trace_df,1)
    im.roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),size_smooth,smooth_type);
end


movecorr = zeros(size(im.roi_trace_df,1),1);
for curr_cell = 1:size(im.roi_trace_df,1)
      
    % take velocity of lever to correlate
    lever_movement = smooth(abs(diff(lever_force_resample)),2);
    % make the first value 0 to keep same size
    lever_movement = [0;lever_movement];
    
    if plot_movecorr == 1;
        figure;hold on;
        plot([1:size(im.roi_trace_df,2)]/framerate,im.roi_trace_df(curr_cell,:),'k');
        plot([1:length(lever_force_resample)]/framerate,lever_force_resample,'k');
        plot([1:length(lever_force_resample)]/framerate,lever_movement,'b');
        
        % if you want to plot the raw lever_force, but not really necessary
        %plot([1:length(lever_force)]/xsg_sample_rate,lever_force,'--r');
        if exist('im.roi_trace_df_smooth')
            plot([1:size(im.roi_trace_df_smooth,2)]/framerate,im.roi_trace_df_smooth(curr_cell,:),'--r','LineWidth',2);
        end
    end
    
    % call the minimum of the raw the low cutoff
    % estimate the noise from smoothed - raw
    noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
    %cutoff = 4*mean(noise_est);
    cutoff = 0.5;
    curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
    over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
    curr_trace_cutoff(over_cutoff_indx) = im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
    
    % if less than 10 seconds is above cutoff, then continue
    if sum(over_cutoff_indx) < framerate*10
        continue
    end
    
    if plot_movecorr == 1;
        plot([1:size(curr_trace_cutoff,2)]/framerate,curr_trace_cutoff,'--g','LineWidth',2);
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
%     
%     ACTIVATE THIS WHEN YOU WANT TO SHUFFLE
%     
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
    
    if plot_movecorr == 1;
        figure; hold on
        errorbar([1 2],[active_lever inactive_lever],...
            [active_lever_bootstrp_ci_range(1) inactive_lever_bootstrp_ci_range(1)], ...
            [active_lever_bootstrp_ci_range(2) inactive_lever_bootstrp_ci_range(2)],'r');
        errorbar([1 2],[active_lever_shuffle inactive_lever_shuffle],...
            [active_lever_shuffle_bootstrp_ci_range(1) inactive_lever_shuffle_bootstrp_ci_range(1)], ...
            [active_lever_shuffle_bootstrp_ci_range(2) inactive_lever_shuffle_bootstrp_ci_range(2)],'k');
    end
    
    
    % criterion for movement correlated
    if active_lever_bootstrp_ci(1) > inactive_lever_bootstrp_ci(2) && ...
        active_lever_bootstrp_ci(1) > active_lever_shuffle_bootstrp_ci(2)
        movecorr(curr_cell) = 1;
    elseif active_lever_bootstrp_ci(2) < inactive_lever_bootstrp_ci(1) && ...
        active_lever_bootstrp_ci(2) < active_lever_shuffle_bootstrp_ci(1)
        movecorr(curr_cell) = -1;
    end

    curr_cell/size(im.roi_trace_df,1)
    
end

% take raw corelation for strength of relationship
movecorr_value = zeros(size(movecorr));
for curr_cell = 1:length(movecorr)
   if movecorr(curr_cell) == 0
       continue
   end
   
   % call the minimum of the raw the low cutoff
   % estimate the noise from smoothed - raw
   noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
   cutoff = 3*mean(noise_est);
   curr_trace_cutoff = zeros(size(im.roi_trace_df_smooth(curr_cell,:)));
   over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
   curr_trace_cutoff(over_cutoff_indx) = im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
   
   temp_corr = corrcoef(curr_trace_cutoff,lever_movement);
   movecorr_value(curr_cell) = temp_corr(2);
end

% Create a thresholded trace
im.roi_trace_df_smooth_cutoff = zeros(size(im.roi_trace_df_smooth));
for curr_cell = 1:size(im.roi_trace_df_smooth,1)

    % call the minimum of the raw the low cutoff
    % estimate the noise from smoothed - raw
    noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
    cutoff = 2*mean(noise_est);
    curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
    over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
    im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
        im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
end

move_cells = find(movecorr == 1);
antimove_cells = find(movecorr == -1);
uncorr_cells = find(movecorr == 0);

%cell_order_indx = [move_cells;antimove_cells;uncorr_cells];

% order by correlation value instead of just binary
cell_order_indx = [movecorr_value (1:length(movecorr_value))'];
cell_order_indx = sortrows(cell_order_indx,1);
cell_order_indx = flipud(cell_order_indx(:,2));

trace_sort = [];
trace_sort = im.roi_trace_df_smooth_cutoff(cell_order_indx,:);

% trace_sort = [];
% trace_sort = [trace_sort;im.roi_trace_df_smooth_cutoff(move_cells,:)];
% trace_sort = [trace_sort;im.roi_trace_df_smooth_cutoff(antimove_cells,:)];
% trace_sort = [trace_sort;im.roi_trace_df_smooth_cutoff(uncorr_cells,:)];

plot_color = zeros(size(im.roi_trace_df_smooth_cutoff,1),3);
plot_color(1:length(move_cells),1) = 1;
plot_color(length(move_cells)+1:length(move_cells)+length(antimove_cells),2) = 1;



%% Find correlation between unprocessed thresholded traces, kmeans cluster

smooth_size = 10;
roi_trace_df_smooth = nan(size(im.roi_trace_df));
for i = 1:size(im.roi_trace_df,1);
   roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess'); 
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

% actually... just do kmeans on the regular smoothed trace
roi_trace_df_smooth_cutoff = roi_trace_df_smooth;

% for kmeans: get rid of traces with no activity
active_cells = find(sum(roi_trace_df_smooth_cutoff,2) ~= 0);

% group by k-means via correlation
kmeans_idx = kmeans( ...
    roi_trace_df_smooth_cutoff(active_cells,:),3,'Distance','correlation');

% sort by group
roi_trace_kmeans = [roi_trace_df_smooth_cutoff(active_cells,:) kmeans_idx];
roi_trace_kmeans = sortrows(roi_trace_kmeans,size(roi_trace_kmeans,2));
roi_trace_kmeans = roi_trace_kmeans(:,1:end-1);


%% jPCA/rotation analysis a la Churchland & Shenoy
% using the top 6 PCs

reward_frames = bhv.reward_frames;
% if you want just random frames, for testing
%reward_frames = randi([round(min(reward_frames)), ...
%    round(max(reward_frames))], length(reward_frames));

frames_span = 20;
frames_surround = frames_span*2+2;

% pop out data around the time of rewarded lever presses
total_reward_frames = frames_surround*length(reward_frames);
num_rois = size(im.roi_trace_df,1);
lever_reward_trace = zeros(num_rois,length(reward_frames)*frames_surround);

smooth_size = 30;
im.roi_trace_df_smooth = nan(size(im.roi_trace_df));
for i = 1:size(im.roi_trace_df,1);
   im.roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess'); 
end
    
% Create a thresholded trace
im.roi_trace_df_smooth_cutoff = zeros(size(im.roi_trace_df_smooth));
for curr_cell = 1:size(im.roi_trace_df_smooth,1)

    % call the minimum of the raw the low cutoff
    % estimate the noise from smoothed - raw
    noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
    cutoff = 3*mean(noise_est);
    curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
    over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
    im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
        im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
end

% % Use positive diff of im.roi_trace_df_smooth as PCA
% im.roi_trace_df_smooth_diff = diff(im.roi_trace_df_smooth,1,2);
% im.roi_trace_df_smooth_diff(im.roi_trace_df_smooth_diff < 0) = 0;

%%%%%%%%%%%% DEFINE WHAT THE jPCA IS DONE ON

%%%%%%%%%%%% THIS IS FOR SNIPS AROUND CUE, THEN REWARD, TIME

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%   TRY THIS AGAIN! DIDN'T WORK FIRST TIME BECAUSE XSG OFFSET NOT DONE
%%%%%%%%%%%%%%%%%%%%%

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

lever_reward_trace = zeros(size(im.roi_trace_df,1), ...
    frames_surround*length(good_trials));

for i = 1:length(reward_frames_good);
   curr_reward_frame = round(reward_frames_good(i)); 
   curr_reward_frames = curr_reward_frame:curr_reward_frame+frames_span;
   curr_reward_trace = im.roi_trace_df_smooth_cutoff(:,curr_reward_frames);
   
   curr_cue_frame = round(cue_frames_good(i));
   curr_cue_frames = curr_cue_frame:curr_cue_frame+frames_span;
   curr_cue_trace = im.roi_trace_df_smooth_cutoff(:,curr_cue_frames);
   
   lever_reward_trace(:,((i-1)*frames_surround+1): ...
        ((i-1)*frames_surround)+frames_surround) = ...
        [curr_reward_trace curr_cue_trace];
end


%%%%%%%%%%%% THIS IS FOR SNIPS AROUND THE REWARD TIME
% for i = 1:length(bhv.reward_frames);
%    curr_frame = round(bhv.reward_frames(i)); 
%    curr_reward_frames = curr_frame-frames_span:curr_frame+frames_span;
%    curr_reward_trace = im.roi_trace_df_smooth_cutoff(:,curr_reward_frames);
%    lever_reward_trace(:,((i-1)*frames_surround+1): ...
%         ((i-1)*frames_surround)+frames_surround) = curr_reward_trace;
% end

% %%%%%%%%%%%% THIS IS FOR SNIPS AROUND THE CUE TIME
% total_cue_frames = frames_surround*length(bhv.cue_frames);
% num_rois = size(im.roi_trace_df,1);
% lever_reward_trace = zeros(num_rois,length(bhv.cue_frames)*frames_surround);
% for i = 1:length(bhv.cue_frames);
%    curr_frame = round(bhv.cue_frames(i)); 
%    curr_cue_frames = curr_frame:curr_frame+frames_span-1;
%    curr_cue_trace = im.roi_trace_df_smooth_cutoff(:,curr_cue_frames);
%    lever_reward_trace(:,((i-1)*frames_surround+1): ...
%         ((i-1)*frames_surround)+frames_surround) = curr_cue_trace;
% end

%%%%%%%%%%%%

% % subtract cell all-trial mean, they use this to separate trials
% for i = 1:size(im.roi_trace_df,1)
%     trial_reshape = reshape(lever_reward_trace(i,:), ...
%         frames_surround,length(reward_frames))';
%     trial_mean = mean(trial_reshape);
%     
%     trial_meansub = (trial_reshape - ...
%         repmat(trial_mean,length(reward_frames),1))';
%     lever_reward_trace(i,:) = trial_meansub(:);
% end


pca_split = 'rotation';

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

switch pca_split
    case 'rotation'
        k = fminsearch(@(x) Churchland2012_jPCA_rotation_fminsearch ...
            (x,X_tilde,X_dot),k_start,options);
    case 'expansion'
        k = fminsearch(@(x) Churchland2012_jPCA_expansion_fminsearch ...
            (x,X_tilde,X_dot),k_start,options);
end

% construct final M_skew
idx = ones(6);
tril_idx = logical(tril(idx,-1));
triu_idx = logical(triu(idx,1));
M_skew = zeros(6);
M_skew(tril_idx) = k;
M_skew_t = M_skew';
switch pca_split
    case 'rotation'
        M_skew_t(tril_idx) = -k;
    case 'expansion'
        M_skew_t(tril_idx) = k;
end
M_skew = M_skew_t';

% find the largest eigenvalues and pull out the top 2 (complex conj pair)
clear j
[V D] = eig(M_skew);
D = diag(D);
switch pca_split
    case 'rotation'
        jPCA_idx = find(abs(D) == max(abs(D)));
        jPCA_1 = V(:,jPCA_idx(1)) + V(:,jPCA_idx(2)); 
        jPCA_2 = j*(V(:,jPCA_idx(1)) - V(:,jPCA_idx(2))); 
    case 'expansion'
        jPCA_idx = [D [1:length(D)]'];
        jPCA_idx = sortrows(jPCA_idx,1);
        jPCA_idx = jPCA_idx(end:-1:end-1,2);
        jPCA_1 = V(:,jPCA_idx(1)) + V(:,jPCA_idx(2)); 
        jPCA_2 = V(:,jPCA_idx(1)) - V(:,jPCA_idx(2)); 
end

% get the final jPCs, project data onto jPC space
pop_jPCs = X*[jPCA_1 jPCA_2];
figure;hold on;
title('jPCs')
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

% plot jPC coefficients
jPC_coeffs = lever_reward_trace(:,1:end-1)'\pop_jPCs;
figure; title('jPC coefficients');
plot(jPC_coeffs(:,1) - jPC_coeffs(:,2),'.');

% My attempt at error bars on this process: don't think it makes sense
% % NOPE this doesn't really make sense, instead do bootstrapped errorbars
% % get coefficients for shuffled traces
% num_shuffle = 30;
% jPC_coeffs_shuffle = nan(size(lever_reward_trace,1),num_shuffle);
% for curr_shuffle= 1:num_shuffle;
%     curr_perm = randperm(size(lever_reward_trace,2)-1);
%     curr_coeffs = [];
%     % turn warning off, since usually rank deficient
%     warning off
%     curr_coeffs = lever_reward_trace(:,curr_perm)'\pop_jPCs;
%     warning on
%     jPC_coeffs_shuffle(:,curr_shuffle) = -diff(curr_coeffs,1,2);
%     curr_shuffle
% end
% % get percentiles to errorbar shuffled data correlations
% jPC_error = prctile(jPC_coeffs_shuffle',[0.5 99.5]);
% 
% %%%% try here bootstrapped instead
% warning off
% jPC_error = bootstrp(5,@(x) -diff(x'\pop_jPCs,1,2),lever_reward_trace(:,1:end-1));
% warning on 

% plot the jPCs for each rewarded press
figure; 
for jPC_num = 1:2
    % tack on an extra point? lost somewhere...
    curr_jPC = [pop_jPCs(:,jPC_num);pop_jPCs(end,jPC_num)];
    jPC_trial_reshape = reshape(curr_jPC, ...
        frames_surround,length(bhv.cue_frames))';
    subplot(2,1,jPC_num)
    imagesc(jPC_trial_reshape);
end
title('Trial-by-trial jPCs');

% get pc scores
jPC_scores = lever_reward_trace(:,1:end-1)*[pop_jPCs(:,1) -pop_jPCs(:,2)];
% get pc scores above zero
pop_jPCs_pos = pop_jPCs;
pop_jPCs_pos(pop_jPCs_pos < 0) = 0;
jPC_scores_pos = lever_reward_trace(:,1:end-1)* ...
    [pop_jPCs_pos(:,1) -pop_jPCs_pos(:,2)];


figure; plot3([1:size(lever_reward_trace,1)],jPC_scores(:,1), ...
    jPC_scores(:,2),'.');

figure;plot(jPC_scores(:,1) - jPC_scores(:,2),'.');

% get straight correlation of lever-surrounding traces, sort kmeans
active_cells = find(sum(lever_reward_trace,2) ~= 0);
lever_reward_trace_corr = corrcoef(lever_reward_trace');
kmeans_idx = kmeans(lever_reward_trace(active_cells,:), ...
    10,'Distance','correlation');
% sort by group
roi_trace_kmeans = [lever_reward_trace(active_cells,:) kmeans_idx];
roi_trace_kmeans = sortrows(roi_trace_kmeans,size(roi_trace_kmeans,2));
roi_trace_kmeans = roi_trace_kmeans(:,1:end-1);


% 
% % get correlations between cells and jPCs
% jPC_long = im.roi_trace_df_smooth(:,1:end-1)'*jPC_coeffs;
% jPC_corr = zeros(size(im.roi_trace_df_smooth,1),2);
% for i = 1:size(lever_reward_trace,1)
%    temp_corr1 = corrcoef(jPC_long(:,1),im.roi_trace_df_smooth(i,1:end-1));
%    temp_corr2 = corrcoef(jPC_long(:,2),im.roi_trace_df_smooth(i,1:end-1));
%    jPC_corr(i,1) = temp_corr1(2);
%    jPC_corr(i,2) = temp_corr2(2);
% end



% subplot(2,1,2)
% title('PCs')
% for i = 1:length(reward_frames)
%     curr_frames = (i-1)*frames_surround+1:(i)*frames_surround-1;
%     z = zeros(size(curr_frames));
%     x = pcs(curr_frames,1)';
%     y = pcs(curr_frames,2)';
%     z = pcs(curr_frames,3)';
%     col = curr_frames - curr_frames(1);
%     surface([x;x],[y;y],[z;z],[col;col],'facecol','no', ...
%         'edgecol','interp','linew',2);
%     arrow([x(end-1) y(end-1)],[x(end) y(end)],'length',10)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% testing stuff



% % for kmeans: get rid of traces with no activity
% active_cells = find(sum(lever_reward_trace,2) ~= 0);
% 
% % group by k-means via correlation
% kmeans_idx = kmeans( ...
%     lever_reward_trace(active_cells,:),5,'Distance','correlation');
% 
% % sort by group
% roi_trace_kmeans = [lever_reward_trace(active_cells,:) kmeans_idx];
% roi_trace_kmeans = sortrows(roi_trace_kmeans,size(roi_trace_kmeans,2));
% roi_trace_kmeans = roi_trace_kmeans(:,1:end-1);
% 
% pc_score = lever_reward_trace*pcs;
% figure;plot(pc_score(:,1),pc_score(:,2),'.');

%% Correlate cell with other cells

curr_cell = 75;

smooth_size = 30;
roi_trace_df_smooth = nan(size(im.roi_trace_df));
for i = 1:size(im.roi_trace_df,1);
    roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess');
end

% try a cutoff
% Create a thresholded trace
roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
for i = 1:size(roi_trace_df_smooth,1)
    
    % call the minimum of the raw the low cutoff
    % estimate the noise from smoothed - raw
    noise_est = abs(im.roi_trace_df(i,:) - roi_trace_df_smooth(i,:));
    cutoff = 2*mean(noise_est);
    curr_trace_cutoff = nan(size(roi_trace_df_smooth(i,:)));
    over_cutoff_indx = roi_trace_df_smooth(i,:) > cutoff;
    roi_trace_df_smooth_cutoff(i,over_cutoff_indx) = ...
        roi_trace_df_smooth(i,over_cutoff_indx);
end


% try a diff, set negs to zero
roi_trace_df_smooth_diff = diff(roi_trace_df_smooth);
roi_trace_df_smooth_diff(roi_trace_df_smooth_diff < 0) = 0;

all_cc = corrcoef(roi_trace_df_smooth_cutoff');
curr_cc = all_cc(curr_cell,:);
curr_cc(curr_cell) = NaN;

figure; hold on;
plot(pyr_cells,curr_cc(pyr_cells),'.')
plot(gad_cells,curr_cc(gad_cells),'.r')

% Find the maximum deconv xcorr
curr_max_xcorr = zeros(size(im.roi_trace_df,1),1);
for i = 1:size(im.roi_trace_df,1);
   temp_xcorr = xcorr(im.oopsi(curr_cell,:),im.roi_concat_oopsi(i,:));
   curr_max_xcorr(i) = max(temp_xcorr)/(sum(im.roi_concat_oopsi(curr_cell,:))* ...
       sum(im.oopsi(i,:)));    
end
curr_max_xcorr(curr_cell) = 0;



%% Get smoothed trace around reward, do PCA
smooth_size = 30;
roi_trace_df_smooth = nan(size(roi_trace_df));
for i = 1:size(roi_trace_df,1);
    roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
end

roi_reward_all = cell(size(roi_trace_df,1),1);
fprintf('%2d',0)
for curr_cell = 1:size(roi_trace_df,1)
    % plot around reward
    frames_surround = 30;
    roi_reward = nan(length(bhv.reward_frames),frames_surround*2+1);
    for i = 1:length(bhv.reward_frames);
        curr_reward = round(bhv.reward_frames(i));
        roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
            continue
        end
        roi_reward(i,:) = roi_trace_df_smooth(curr_cell,roi_reward_frames);
    end
    roi_reward_all{curr_cell} = roi_reward;
   fprintf('%c%c%2d',8,8,100*floor(curr_cell/size(roi_trace_df,1))); 
end

roi_reward_all_concat = horzcat(roi_reward_all{:})';
% make NaNs zeros
roi_reward_all_concat(isnan(roi_reward_all_concat)) = 0;
% subtract the mean
roi_reward_all_concat_meansub = roi_reward_all_concat - ...
    repmat(mean(roi_reward_all_concat,2),1,size(roi_reward_all_concat,2));
coeffs = princomp(roi_reward_all_concat_meansub);
pcs = roi_reward_all_concat_meansub*coeffs(:,1:3);
pc_score = roi_reward_all_concat_meansub'*pcs;



roi_reward_all_long = cellfun(@(x) x',roi_reward_all,'UniformOutput',false);
roi_reward_all_long = cellfun(@(x) x(:),roi_reward_all,'UniformOutput',false);
roi_reward_all_long = [roi_reward_all_long{:}]';
% replace nans with zeros
roi_reward_all_long(isnan(roi_reward_all_long)) = 0;


%% Classify behavior/movement aligned cells: thresholded trace peri-event
roi_trace_df = im.roi_trace_df;
% smooth_size = 10;
% roi_trace_df_smooth = nan(size(roi_trace_df));
% for i = 1:size(roi_trace_df,1);
%     roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
% end
% 
% % Create a thresholded trace
% roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
% for curr_cell = 1:size(roi_trace_df_smooth,1)
%     noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
%     cutoff = 4*mean(noise_est);
%     curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
%     over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
%     roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
%         roi_trace_df_smooth(curr_cell,over_cutoff_indx);
% end

thresh = 0.2;
roi_trace_df_smooth_cutoff = im.roi_trace_df;
roi_trace_df_smooth_cutoff(roi_trace_df_smooth_cutoff < thresh) = 0;

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
cells = [1:size(roi_trace_df,1)];    

if exist('roi_labels','var')
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';
    
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    gad_cells(ismember(gad_cells,contaminated)) = [];
else
    pyr_cells = cells;
    gad_cells = [];
end
frames_back = 1;

num_rois = size(roi_trace_df,1);
num_frames = size(roi_trace_df,2);

precue_active = zeros(num_rois,length(good_trials));
cue_active = zeros(num_rois,length(good_trials));
reward_active = zeros(num_rois,length(good_trials));

cue_preactive = zeros(num_rois,length(good_trials));
reward_preactive = zeros(num_rois,length(good_trials));

cue_aligned_trials = zeros(num_rois,length(good_trials));
reward_aligned_trials = zeros(num_rois,length(good_trials));

cue_alignment_indx = zeros(size(roi_trace_df,1),1);
reward_alignment_indx = zeros(size(roi_trace_df,1),1);
for curr_cell = [pyr_cells gad_cells];
    
    % EVENTS: go through valid trials, check for activity 
    for i = 1:length(good_trials)
        % find if activity after trial time
        frames_forward = trial_time(i);
        % check that this isn't out of bounds
        if frames_forward + reward_frames_good_round(i) > num_frames
            continue
        end
        precue_frames_search = precue_frames_good_round(i) : ...
            precue_frames_good_round(i) + frames_forward;
        cue_frames_search = cue_frames_good_round(i): ...
            cue_frames_good_round(i) + frames_forward;
        reward_frames_search = reward_frames_good_round(i): ...
            reward_frames_good_round(i) + frames_forward;
        
        precue_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,precue_frames_search));
        cue_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,cue_frames_search));
        reward_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,reward_frames_search));
        
        % find if activity before alignment points
        frames_forward = trial_time(i);
        cue_preframes_search = cue_frames_good_round(i) - frames_back: ...
            cue_frames_good_round(i);
        reward_preframes_search = reward_frames_good_round(i) - frames_back: ...
            reward_frames_good_round(i);

        cue_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,cue_preframes_search));
        reward_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,reward_preframes_search));      
    end
     
    % to classify: cross thresh in given time period, not cross thresh in 
    % the time before, not cross thresh in the other given time period
    cue_aligned_trials(curr_cell,:) = cue_active(curr_cell,:) & ...
        ~precue_active(curr_cell,:) & ~cue_preactive(curr_cell,:);
    reward_aligned_trials(curr_cell,:) = reward_active(curr_cell,:) & ...
        ~cue_active(curr_cell,:) & ~reward_preactive(curr_cell,:);
    
    % calculate alignment index
    num_frames = size(roi_trace_df_smooth_cutoff,2);
    num_trials = length(cue_frames_good);
    % define chance overlap of calcium transient with surround frames
    thresh_start = diff(roi_trace_df_smooth_cutoff(curr_cell,:) > 0);
    if roi_trace_df_smooth_cutoff(curr_cell,1) > 0;
        thresh_start(find(thresh_start < 0,1)) = 0;
    end
    if roi_trace_df_smooth_cutoff(curr_cell,end) > 0;
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

%% Classify behavior/movement aligned cells: peri-event w/ deconvolution

% Use the thresholded oopsi trace
roi_trace_df_smooth_cutoff = im.oopsi;
roi_trace_df_smooth_cutoff(im.oopsi < 0.2) = 0;

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
cells = [1:size(roi_trace_df,1)];    

if exist('roi_labels','var')
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';
    
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    gad_cells(ismember(gad_cells,contaminated)) = [];
else
    pyr_cells = cells;
    gad_cells = [];
end
frames_back = 1;

num_rois = size(roi_trace_df,1);
num_frames = size(roi_trace_df,2);

precue_active = zeros(num_rois,length(good_trials));
cue_active = zeros(num_rois,length(good_trials));
reward_active = zeros(num_rois,length(good_trials));

cue_preactive = zeros(num_rois,length(good_trials));
reward_preactive = zeros(num_rois,length(good_trials));

cue_aligned_trials = zeros(num_rois,length(good_trials));
reward_aligned_trials = zeros(num_rois,length(good_trials));

cue_alignment_indx = zeros(size(roi_trace_df,1),1);
reward_alignment_indx = zeros(size(roi_trace_df,1),1);
for curr_cell = [pyr_cells gad_cells];
    
    % EVENTS: go through valid trials, check for activity 
    for i = 1:length(good_trials)
        % find if activity after trial time
        frames_forward = trial_time(i);
        % check that this isn't out of bounds
        if frames_forward + reward_frames_good_round(i) > num_frames
            continue
        end
        precue_frames_search = precue_frames_good_round(i) : ...
            precue_frames_good_round(i) + frames_forward;
        cue_frames_search = cue_frames_good_round(i): ...
            cue_frames_good_round(i) + frames_forward;
        reward_frames_search = reward_frames_good_round(i): ...
            reward_frames_good_round(i) + frames_forward;
        
        precue_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,precue_frames_search));
        cue_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,cue_frames_search));
        reward_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,reward_frames_search));
        
        % find if activity before alignment points
        frames_forward = trial_time(i);
        cue_preframes_search = cue_frames_good_round(i) - frames_back: ...
            cue_frames_good_round(i);
        reward_preframes_search = reward_frames_good_round(i) - frames_back: ...
            reward_frames_good_round(i);

        cue_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,cue_preframes_search));
        reward_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,reward_preframes_search));      
    end
     
    % to classify: cross thresh in given time period, not cross thresh in 
    % the time before, not cross thresh in the other given time period
    cue_aligned_trials(curr_cell,:) = cue_active(curr_cell,:) & ...
        ~precue_active(curr_cell,:) & ~cue_preactive(curr_cell,:);
    reward_aligned_trials(curr_cell,:) = reward_active(curr_cell,:) & ...
        ~cue_active(curr_cell,:) & ~reward_preactive(curr_cell,:);
    
    % calculate alignment index
    num_frames = size(roi_trace_df_smooth_cutoff,2);
    num_trials = length(cue_frames_good);
    % define chance overlap of calcium transient with surround frames
    thresh_start = diff(roi_trace_df_smooth_cutoff(curr_cell,:) > 0);
    if roi_trace_df_smooth_cutoff(curr_cell,1) > 0;
        thresh_start(find(thresh_start < 0,1)) = 0;
    end
    if roi_trace_df_smooth_cutoff(curr_cell,end) > 0;
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

%% Get distance-activity relationships

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

cells = [1:size(im.roi_trace_df,1)];    

gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));
contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';

pyr_cells(ismember(pyr_cells,contaminated)) = [];
gad_cells(ismember(gad_cells,contaminated)) = [];

%gad_cells = pv;

% % these cells are particularly strong
% contaminated = [pyr_cells([102 94 93 78 52 33 11 9 5 3])];
% 
% pyr_cells(ismember(pyr_cells,contaminated)) = [];
% gad_cells(ismember(gad_cells,contaminated)) = [];

% this excludes the ones in the lower right corner
%gad_cells = [29 12 49 109 192 111 21 96 40 51 114 99 204 ...
%    76 101 119 97 85 55 11 47 5 134 118 102];



%gad_cells = [pv som];

% this includes only the group from kmeans
%gad_cells = gad_cells(active_cells(kmeans_idx == 2));
%gad_cells =  [109    32    22    21    51   114    99   251   119  ...
%                    6    85   1     4      118];
%gad_cells = [111 11 85 192 204 114 99];
%gad_cells = [85 114 99];

smooth_size = 30;
im.roi_trace_df_smooth = nan(size(im.roi_trace_df));
for i = 1:size(im.roi_trace_df,1);
    im.roi_trace_df_smooth(i,:) = smooth(im.roi_trace_df(i,:),smooth_size,'loess');
end

% Create a thresholded trace (seperate thresh for gad/pyr)
im.roi_trace_df_smooth_cutoff = zeros(size(im.roi_trace_df_smooth));
im.roi_trace_df_cutoff = zeros(size(im.roi_trace_df_smooth));

im.roi_trace_df_smooth_cutoff =im.roi_trace_df_smooth;
im.roi_trace_df_cutoff = im.roi_trace_df_smooth;

% % get rid of fluorescence spikes at loop starts (first 10 frames)
% num_loops = length(im.roi_trace_long_split)/8;
% file_lengths = cellfun(@length,im.roi_trace_long_split);
% file_lengths_loops = reshape(file_lengths,8,num_loops);
% loop_lengths = cumsum(sum(file_lengths_loops));
% loop_starts = [1 loop_lengths(1:end-1) + 1];
% loop_start_cut = repmat(loop_starts,10,1) + repmat([0:9]',1,num_loops);
% im.roi_trace_df_cutoff(:,loop_start_cut) = [];

cutoff = 0.5;
im.roi_trace_df_cutoff(im.roi_trace_df_cutoff < cutoff) = 0;
    
% % this was a really dumb way to do this.
% for curr_pyr = 1:length(pyr_cells)
%     curr_cell = pyr_cells(curr_pyr);
%     noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
%     cutoff = 5*mean(noise_est);
%     curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
%     over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
%     im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
%         im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
%     im.roi_trace_df_cutoff(curr_cell,over_cutoff_indx) = ...
%         im.roi_trace_df(curr_cell,over_cutoff_indx);
% end
% for curr_gad = 1:length(gad_cells)
%     curr_cell = gad_cells(curr_gad);
%     noise_est = abs(im.roi_trace_df(curr_cell,:) - im.roi_trace_df_smooth(curr_cell,:));
%     cutoff = 4*mean(noise_est);
%     curr_trace_cutoff = nan(size(im.roi_trace_df_smooth(curr_cell,:)));
%     over_cutoff_indx = im.roi_trace_df_smooth(curr_cell,:) > cutoff;
%     im.roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
%         im.roi_trace_df_smooth(curr_cell,over_cutoff_indx);
%     im.roi_trace_df_cutoff(curr_cell,over_cutoff_indx) = ...
%         im.roi_trace_df(curr_cell,over_cutoff_indx);
% end

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
cue_reward_cross = nan(size(im.roi_trace_df,1),length(cue_frames_good));
cue_reward_max = nan(size(im.roi_trace_df,1),length(cue_frames_good));
cue_reward_trace = cell(size(im.roi_trace_df,1),length(cue_frames_good));
for curr_cell = 1:size(im.roi_trace_df,1);
    frames_back = 10;
    frames_forward = 40;
    frames_surround = frames_back+frames_forward+1;

    % plot around cue
    cue_reward = nan(length(cue_frames_good),frames_surround);
    for i = 1:length(cue_frames_good);
        curr_reward = round(bhv.cue_frames(i));
        roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(im.roi_trace_df,2));
            continue
        end
        cue_reward(i,:) = im.roi_trace_df_smooth_cutoff(curr_cell,roi_reward_frames);
        if any(cue_reward(i,:) > 0);
            cue_reward_cross(curr_cell,i) = find(cue_reward(i,:) > 0,1);
            cue_reward_max(curr_cell,i) = max(cue_reward(i,:));
        end
        % save the whole trace
        cue_reward_trace{curr_cell,i} = cue_reward(i,:);
    end
end

% get start times relative to rewards
lever_reward_cross = cue_reward_cross - ...
    repmat([reward_frames_good-cue_frames_good]',size(cue_reward_cross,1),1);

% plot the start times from the cue
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

% figure;hold on;
% errorb(find(cue_cells), nanmean(cue_reward_cross(cue_cells,:),2), ...
%     cue_reward_sem(cue_cells),'color','k');
% errorb(find(lever_cells), nanmean(lever_reward_cross(lever_cells,:),2), ...
%     lever_reward_sem(lever_cells),'color','r');

%%% For TK: activity of pyramidal cells vs. gad cells

cue_reward_bin = zeros(size(cue_reward_cross));
cue_reward_bin(cue_reward_cross > 0) = 1;

pyr_num = sum(cue_reward_bin(pyr_cells,:))/length(pyr_cells);
gad_num = sum(cue_reward_bin(gad_cells,:))/length(gad_cells);

cue_reward_max_zeros = cue_reward_max;
cue_reward_max_zeros(isnan(cue_reward_max)) = 0;
pyr_max = sum(cue_reward_max_zeros(pyr_cells,:));
gad_max = sum(cue_reward_max_zeros(gad_cells,:));

%%% Also for TK: correlated behavior-related activity by cell type
% replaces nans with zeros
cue_reward_max_zeros = cue_reward_max;
cue_reward_max_zeros(isnan(cue_reward_max_zeros)) = 0;
%cue_reward_max_corr = corrcoef(cue_reward_max_zeros');
cue_reward_max_corr = corrcoef(im.roi_trace_df_cutoff');

% create matrix of distances
roi_dist = nan(length(roi_center_x));
for i = 1:length(roi_dist)
    for j = 1:length(roi_dist)
    curr_dist = sqrt((roi_center_x(i) - roi_center_x(j))^2 + ...
        (roi_center_y(i) - roi_center_y(j))^2);
    roi_dist(i,j) = curr_dist;        
    end
end

% set all self-correlations to NaN
cue_reward_max_corr(eye(length(cue_reward_max_corr)) ~= 0) = NaN;

% prepare bins
% num_bins = 50;
% binsize = linspace(min(roi_dist(:)),max(roi_dist(:)),num_bins);
binsize = 1:20:401;
num_bins = length(binsize);

% pyr-pyr pairs
pyrpyr_corr = tril(cue_reward_max_corr(pyr_cells,pyr_cells),-1);
pyrpyr_dist = tril(roi_dist(pyr_cells,pyr_cells),-1);

[n,pyrpyr_bin] = histc(pyrpyr_dist(:),binsize);
pyrpyr_binMean = grpstats(pyrpyr_corr(:),pyrpyr_bin);
pyrpyr_binSem = grpstats(pyrpyr_corr(:),pyrpyr_bin,'sem');
    
% gad-gad pairs
gadgad_corr = tril(cue_reward_max_corr(gad_cells,gad_cells),-1);
gadgad_dist = tril(roi_dist(gad_cells,gad_cells),-1);

[n,gadgad_bin] = histc(gadgad_dist(:),binsize);
gadgad_binMean = grpstats(gadgad_corr(:),gadgad_bin);
gadgad_binSem = grpstats(gadgad_corr(:),gadgad_bin,'sem');

% pyr-gad pairs
pyrgad_corr = tril(cue_reward_max_corr(pyr_cells,gad_cells),-1);
pyrgad_dist = tril(roi_dist(pyr_cells,gad_cells),-1);

[n,pyrgad_bin] = histc(pyrgad_dist(:),binsize);
pyrgad_binMean = grpstats(pyrgad_corr(:),pyrgad_bin);
pyrgad_binSem = grpstats(pyrgad_corr(:),pyrgad_bin,'sem');

figure; hold on
errorbar(pyrpyr_binMean(2:end),pyrpyr_binSem(2:end),'color','k','linewidth',2)
errorbar(gadgad_binMean(2:end),gadgad_binSem(2:end),'color','r','linewidth',2)
errorbar(pyrgad_binMean(2:end),pyrgad_binSem(2:end),'color','b','linewidth',2)

set(gca,'XTick',[1:1:num_bins])
set(gca,'XTickLabel',round(binsize(end)/num_bins)*get(gca,'XTick'))
xlabel('Distance (\mum)')
ylabel('Correlation coefficient')
legend({'Pyr-Pyr','Gad-Gad','Pyr-Gad'})
    

% get TK coincidence index (optional)
tkidx = 0;

if tkidx == 1;
    
    coin_data = im.oopsi;
    smooth_kernel = ones(1,10)./10;
    coin_data_smooth = conv2(coin_data,smooth_kernel);
    
    % pyr-pyr pairs
    pyrpyr_event_dot = coin_data_smooth(pyr_cells,:)* ...
        coin_data_smooth(pyr_cells,:)';
    pyrpyr_event_mean = (sum(coin_data_smooth(pyr_cells,:),2)/size(coin_data_smooth,2)* ...
        sum(coin_data_smooth(pyr_cells,:),2)'/size(coin_data_smooth,2));
    pyrpyr_event_coinChance = pyrpyr_event_mean* ...
        size(coin_data_smooth,2);
    pyrpyr_event_coinChanceGeo = sqrt(pyrpyr_event_mean)* ...
        size(coin_data_smooth,2);
    pyrpyr_tkCorr = (pyrpyr_event_dot-pyrpyr_event_coinChance)./ ...
        pyrpyr_event_coinChanceGeo;
    % get rid of self-correlations
    pyrpyr_tkCorr(eye(size(pyrpyr_tkCorr))>0) = NaN;
    for i = 1:num_bins
        flagBinMembers = (pyrpyr_bin == i);
        binMembers = pyrpyr_tkCorr(flagBinMembers);
        pyrpyr_binMean(i) = nanmean(binMembers);
        pyrpyr_binSem(i) = nanstd(binMembers)/(sqrt(sum(~isnan(binMembers)))/2);
    end
    
    % gad-gad pairs
    gadgad_event_dot = coin_data_smooth(gad_cells,:)* ...
        coin_data_smooth(gad_cells,:)';
    gadgad_event_mean = (sum(coin_data_smooth(gad_cells,:),2)/size(coin_data_smooth,2)* ...
        sum(coin_data_smooth(gad_cells,:),2)'/size(coin_data_smooth,2));
    gadgad_event_coinChance = gadgad_event_mean* ...
        size(coin_data_smooth,2);
    gadgad_event_coinChanceGeo = sqrt(gadgad_event_mean)* ...
        size(coin_data_smooth,2);
    gadgad_tkCorr = (gadgad_event_dot-gadgad_event_coinChance)./ ...
        gadgad_event_coinChanceGeo;
    % get rid of self-correlations
    gadgad_tkCorr(eye(size(gadgad_tkCorr))>0) = NaN;
    for i = 1:num_bins
        flagBinMembers = (gadgad_bin == i);
        binMembers = gadgad_tkCorr(flagBinMembers);
        gadgad_binMean(i) = nanmean(binMembers);
        gadgad_binSem(i) = nanstd(binMembers)/(sqrt(sum(~isnan(binMembers)))/2);
    end
    
    % pyr-gad pairs
    pyrgad_event_dot = coin_data_smooth(pyr_cells,:)* ...
        coin_data_smooth(gad_cells,:)';
    pyrgad_event_mean = (sum(coin_data_smooth(pyr_cells,:),2)/size(coin_data_smooth,2)* ...
        sum(coin_data_smooth(gad_cells,:),2)'/size(coin_data_smooth,2));
    pyrgad_event_coinChance = pyrgad_event_mean* ...
        size(coin_data_smooth,2);
    pyrgad_event_coinChanceGeo = sqrt(pyrgad_event_mean)* ...
        size(coin_data_smooth,2);
    pyrgad_tkCorr = (pyrgad_event_dot-pyrgad_event_coinChance)./ ...
        pyrgad_event_coinChanceGeo;
    % get rid of self-correlations
    pyrgad_tkCorr(eye(size(pyrgad_tkCorr))>0) = NaN;
    for i = 1:num_bins
        flagBinMembers = (pyrgad_bin == i);
        binMembers = pyrgad_tkCorr(flagBinMembers);
        pyrgad_binMean(i) = nanmean(binMembers);
        pyrgad_binSem(i) = nanstd(binMembers)/(sqrt(sum(~isnan(binMembers)))/2);
    end

end

%%%


%% Quantify lever presses

% quantify movement before
lever_force_resample_vel = smooth(abs(diff(bhv.lever_force_resample)),2);

rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames);
good_trials = find(rewarded_trials & bhv.imaged_trials);
reward_frames_good = cellfun(@(x) x.states.reward(1), ...
    bhv.bhv_frames(good_trials));
cue_frames_good = cellfun(@(x) x.states.cue(1), ...
    bhv.bhv_frames(good_trials));

% compensate for loops
reward_frames_good = reward_frames_good + bhv.add_frames(good_trials);
cue_frames_good = cue_frames_good + bhv.add_frames(good_trials);

% get movement right before cues
pre_cue_sec = 1;
pre_cue_frames = round(pre_cue_sec*bhv.framerate);
pre_cue_movement = nan(length(cue_frames_good),1);
for i = 1:length(cue_frames_good)
    pre_cue_movement(i) = sum(lever_force_resample_vel( ...
        round(cue_frames_good(i))-pre_cue_frames: ...
        round(cue_frames_good(i))-1));
end

% get movement between cue and reward
response_movement = nan(length(cue_frames_good),1);
for i = 1:length(reward_frames_good)
    response_frames = round(cue_frames_good(i)):round(reward_frames_good(i));
    response_movement(i) = sum(lever_force_resample_vel(response_frames));
end

% get time between cue and reward
reward_reaction_frames = reward_frames_good - cue_frames_good;

% get time to first movement after cue
reaction_frames = nan(length(cue_frames_good),1);
for i = 1:length(reward_frames_good)
    response_frames = round(cue_frames_good(i)):round(reward_frames_good(i));
    move_thresh = lever_force_resample_vel(response_frames) > 0.04;
    if any(move_thresh)
        reaction_frames(i) = find(move_thresh,1);
    end
end

% get movement during iti before cue
iti_movement = nan(length(cue_frames_good),1);
for i = 2:length(cue_frames_good)
    iti_frames = round(reward_frames_good(i-1)):round(cue_frames_good(i));
    iti_movement(i) = sum(lever_force_resample_vel(iti_frames));
end



%% Get pyr/gad distance-activity radius

figure;hold on
cue_reward_bin_pyr = cue_reward_bin(pyr_cells,:);
cue_reward_bin_gad = cue_reward_bin(gad_cells,:);

radii = 0:0.5:200;
radius_corr = zeros(length(radii),1);

gad_rad = zeros(length(radii),size(cue_reward_bin,2));
for radius = radii;
    for trial = 1:size(cue_reward_bin,2);
        curr_gad_rad_indx = ...
            roi_dist(pyr_cells(find(cue_reward_bin_pyr(:,trial))), ...
            gad_cells(find(cue_reward_bin_gad(:,trial)))) <= radius;
        curr_gad_rad = sum(sum(curr_gad_rad_indx) > 0);
        gad_rad(find(radii == radius),trial) = curr_gad_rad;
    end
end

plot(sum(gad_rad,2))

for test = 1:50
a = randperm(length(gad_cells));
gad_cells_rand = gad_cells(a);

cue_reward_bin_pyr = cue_reward_bin(pyr_cells,:);
cue_reward_bin_gad = cue_reward_bin(gad_cells_rand,:);

radii = 0:0.5:200;
radius_corr = zeros(length(radii),1);

gad_rad = zeros(length(radii),size(cue_reward_bin,2));
for radius = radii;
    for trial = 1:size(cue_reward_bin,2);
        curr_gad_rad_indx = ...
            roi_dist(pyr_cells(find(cue_reward_bin_pyr(:,trial))), ...
            gad_cells(find(cue_reward_bin_gad(:,trial)))) <= radius;
        curr_gad_rad = sum(sum(curr_gad_rad_indx) > 0);
        gad_rad(find(radii == radius),trial) = curr_gad_rad;
    end
end

plot(sum(gad_rad,2),'r')
end



%% Get interneuron correlation with population/radius

cells = [1:size(im.roi_trace_df,1)];    

gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));

contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';

pyr_cells(ismember(pyr_cells,contaminated)) = [];
gad_cells(ismember(gad_cells,contaminated)) = [];

n_best_thresh_pyr = sum(im.oopsi(pyr_cells,:));
n_best_thresh_gad = sum(im.oopsi(gad_cells,:));

% get corrcoef of each gad/som with pyramidal cells in different radii
radii = 1:20:300;
radius_corr = zeros(length(radii),1);

gad_rad = zeros(length(radii),length(gad_cells));
for radius_idx = 1:length(radii);
    radius = radii(radius_idx);
    for curr_gad = 1:size(gad_rad,2)
        curr_gad_cell = gad_cells(curr_gad);
        curr_pyr_radii = find(roi_dist(curr_gad_cell,pyr_cells) <= radius);
        curr_pyr_act = sum(im.oopsi(pyr_cells(curr_pyr_radii),:),1);
        curr_corrcoef = corrcoef(curr_pyr_act, ...
            im.roi_trace_df(curr_gad_cell,:));
        gad_rad(radius_idx,curr_gad) = curr_corrcoef(2);
    end
    disp(num2str(radius_idx/length(radii)));
end

%% Look at interneuron correlations ~110 microns away
radius = 110;
for curr_gad = 1:length(gad_cells)
    curr_som = gad_cells(curr_gad);
    curr_pyr_radii = find(roi_dist(curr_som,pyr_cells) <= radius);
    curr_pyr_act = sum(im.oopsi(pyr_cells(curr_pyr_radii),:),1);
    curr_corrcoef = corrcoef(curr_pyr_act, ...
        roi_trace_df_smooth(curr_som,:));
    gad_rad(radius_idx,curr_gad) = curr_corrcoef(2);
end


%% Classify movement onset/offset related: like cue/reward

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


%%%%% Classify cells as movement onset/offset related

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

% get rid of unpaired movement bouts at start and end
if move_stop_frames(1) < move_start_frames(1) 
    move_stop_frames(1) = [];
end
    
if move_start_frames(end) > move_stop_frames(end)
    move_start_frames(end) = [];
end

% make another category: pre-cue (like reward vs. cue)
median_move_time = median(move_stop_frames - move_start_frames);
premove_frames = move_start_frames - median_move_time;

% make a matrix of trial times for purposes of searching for transients
move_time = move_stop_frames - move_start_frames;

% define the cells to look at
cells = [1:size(roi_trace_df,1)];    

gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));

contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';

pyr_cells(ismember(pyr_cells,contaminated)) = [];
gad_cells(ismember(gad_cells,contaminated)) = [];

frames_back = 1;

num_rois = size(roi_trace_df,1);
num_frames = size(roi_trace_df,2);

premove_active = zeros(num_rois,length(move_start_frames));
move_start_active = zeros(num_rois,length(move_start_frames));
move_stop_active = zeros(num_rois,length(move_start_frames));

move_start_preactive = zeros(num_rois,length(move_start_frames));
move_stop_preactive = zeros(num_rois,length(move_start_frames));

move_start_aligned_trials = zeros(num_rois,length(move_start_frames));
move_stop_aligned_trials = zeros(num_rois,length(move_start_frames));

move_start_alignment_indx = zeros(size(roi_trace_df,1),1);
move_stop_alignment_indx = zeros(size(roi_trace_df,1),1);
for curr_cell = [pyr_cells gad_cells];
    
    % EVENTS: go through valid trials, check for activity 
    for i = 1:length(move_start_frames)
        % find if activity after trial time
        frames_forward = move_time(i);
        % check that this isn't out of bounds
        if frames_forward + move_stop_frames(i) > num_frames
            continue
        end
        premove_frames_search = premove_frames(i) : ...
            premove_frames(i) + frames_forward;
        move_start_frames_search = move_start_frames(i): ...
            move_start_frames(i) + frames_forward;
        move_stop_frames_search = move_stop_frames(i): ...
            move_stop_frames(i) + frames_forward;
        
        premove_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,premove_frames_search));
        move_start_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_start_frames_search));
        move_start_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_stop_frames_search));
        
        % find if activity before alignment points
        frames_forward = move_time(i);
        move_start_preframes_search = move_start_frames(i) - frames_back: ...
            move_start_frames(i);
        move_stop_preframes_search = move_stop_frames(i) - frames_back: ...
            move_stop_frames(i);

        move_start_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_start_preframes_search));
        move_stop_preactive(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_stop_preframes_search));      
    end
     
    % to classify: cross thresh in given time period, not cross thresh in 
    % the time before, not cross thresh in the previous time period
    move_start_aligned_trials(curr_cell,:) = move_start_active(curr_cell,:) & ...
        ~premove_active(curr_cell,:) & ~move_start_preactive(curr_cell,:);
    move_stop_aligned_trials(curr_cell,:) = move_stop_active(curr_cell,:) & ...
        ~move_start_active(curr_cell,:) & ~move_stop_preactive(curr_cell,:);
    
    % calculate alignment index
    num_frames = size(roi_trace_df_smooth_cutoff,2);
    num_trials = length(move_start_frames);
    % define chance overlap of calcium transient with surround frames
    thresh_start = diff(roi_trace_df_smooth_cutoff(curr_cell,:) > 0);
    if roi_trace_df_smooth_cutoff(curr_cell,1) > 0;
        thresh_start(find(thresh_start < 0,1)) = 0;
    end
    if roi_trace_df_smooth_cutoff(curr_cell,end) > 0;
        thresh_start(find(thresh_start > 0,1,'last')) = 0;
    end
    thresh_run = find(thresh_start == -1)-find(thresh_start == 1);
    mean_thresh_run = mean(thresh_run);
    thresh_freq = length(thresh_run)/num_frames;
    p_overlap_align = thresh_freq*(median_move_time+mean_thresh_run);
    p_nonoverlap_prealign = (1-thresh_freq)^...
        ((median_move_time > mean_thresh_run)*median_move_time + ...
        (mean_thresh_run > median_move_time)*mean_thresh_run);
    chance_align = p_overlap_align*p_nonoverlap_prealign*num_trials;
    
    % calculate alignment index
    move_start_alignment_indx(curr_cell) = (sum(move_start_aligned_trials ...
        (curr_cell,:)) - chance_align)/num_trials;
    move_stop_alignment_indx(curr_cell) = (sum(move_stop_aligned_trials ...
        (curr_cell,:)) - chance_align)/num_trials;
 
end


%% Classify movement onset/offset related: starting in surrounding time
% test: the calcium transient just has to start and stop around the time

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


%%%%% Classify cells as movement onset/offset related

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

% get rid of unpaired movement bouts at start and end
if move_stop_frames(1) < move_start_frames(1) 
    move_stop_frames(1) = [];
end
    
if move_start_frames(end) > move_stop_frames(end)
    move_start_frames(end) = [];
end

median_move_time = median(move_stop_frames - move_start_frames);

% make a matrix of trial times for purposes of searching for transients
move_time = move_stop_frames - move_start_frames;

% define the cells to look at
cells = [1:size(roi_trace_df,1)];    

gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));

contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';

pyr_cells(ismember(pyr_cells,contaminated)) = [];
gad_cells(ismember(gad_cells,contaminated)) = [];

num_rois = size(roi_trace_df,1);
num_frames = size(roi_trace_df,2);

premove_active = zeros(num_rois,length(move_start_frames));
move_start_active = zeros(num_rois,length(move_start_frames));
move_stop_active = zeros(num_rois,length(move_start_frames));

move_start_preactive = zeros(num_rois,length(move_start_frames));
move_stop_preactive = zeros(num_rois,length(move_start_frames));

move_start_aligned_trials = zeros(num_rois,length(move_start_frames));
move_stop_aligned_trials = zeros(num_rois,length(move_start_frames));

move_start_alignment_indx = zeros(size(roi_trace_df,1),1);
move_stop_alignment_indx = zeros(size(roi_trace_df,1),1);

frames_surround = 15;

for curr_cell = [pyr_cells gad_cells];
    
    % EVENTS: go through valid trials, check for activity 
    for i = 1:length(move_start_frames)

        % check that this isn't out of bounds
        if frames_surround + move_stop_frames(i) > num_frames
            continue
        end
        
        move_start_frames_search = move_start_frames(i): ...
            move_start_frames(i) + frames_surround;
        move_stop_frames_search = move_stop_frames(i): ...
            move_stop_frames(i) + frames_surround;
        
        % the trace has to be both high and low in the bounds
        move_start_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_start_frames_search)) && ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_start_frames_search) == 0);
        move_start_active(curr_cell,i) = ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_stop_frames_search))  && ...
            any(roi_trace_df_smooth_cutoff(curr_cell,move_stop_frames_search) == 0);
   
    end
     
    % to classify: cross thresh in given time period, not cross thresh in 
    % the time before, not cross thresh in the previous time period
    move_start_aligned_trials(curr_cell,:) = move_start_active(curr_cell,:) & ...
        ~move_stop_active(curr_cell,:);
    move_stop_aligned_trials(curr_cell,:) = move_stop_active(curr_cell,:) & ...
        ~move_start_active(curr_cell,:);
    
    % calculate alignment index
    num_frames = size(roi_trace_df_smooth_cutoff,2);
    num_trials = length(move_start_frames);
    % define chance overlap of calcium transient with surround frames
    thresh_start = diff(roi_trace_df_smooth_cutoff(curr_cell,:) > 0);
    if roi_trace_df_smooth_cutoff(curr_cell,1) > 0;
        thresh_start(find(thresh_start < 0,1)) = 0;
    end
    if roi_trace_df_smooth_cutoff(curr_cell,end) > 0;
        thresh_start(find(thresh_start > 0,1,'last')) = 0;
    end
    thresh_run = find(thresh_start == -1)-find(thresh_start == 1);
    mean_thresh_run = mean(thresh_run);
    thresh_freq = length(thresh_run)/num_frames;
    p_overlap_align = thresh_freq*(frames_surround);
    p_nonoverlap_prealign = (1-thresh_freq)^frames_surround;
    chance_align = p_overlap_align*p_nonoverlap_prealign*num_trials;
    
    % calculate alignment index
    move_start_alignment_indx(curr_cell) = (sum(move_start_aligned_trials ...
        (curr_cell,:)) - chance_align)/num_trials;
    move_stop_alignment_indx(curr_cell) = (sum(move_stop_aligned_trials ...
        (curr_cell,:)) - chance_align)/num_trials;
 
end

%% Classify behavior related: oopsi xcorr/dot

im.oopsi(im.oopsi < 0.1) = 0;

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
        find(temp_xcorr == max(temp_xcorr),1));
    
    [temp_xcorr temp_xcorr_lags] = ...
        xcorr(im.oopsi(curr_cell,:),move_stop_trace,maxlags);
    move_stop_max_xcorr(curr_cell) = max(temp_xcorr);
    move_stop_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
        find(temp_xcorr == max(temp_xcorr),1));
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
        find(temp_xcorr == max(temp_xcorr),1));
    
    [temp_xcorr temp_xcorr_lags] = ...
        xcorr(im.oopsi(curr_cell,:),reward_trace,maxlags);
    reward_max_xcorr(curr_cell) = max(temp_xcorr);
    reward_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
        find(temp_xcorr == max(temp_xcorr),1));
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

% % make things zero based on lag
% reward_tk_ci(reward_max_xcorr_lag < 0) = 0;
% cue_tk_ci(cue_max_xcorr_lag < 0) = 0;

%% Classify behavior related: peakdet xcorr/dot
[peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);

peak_matrix = zeros(size(im.roi_trace_df));
for i = 1:length(peak_frames)
    peak_matrix(i,peak_frames{i}) = 1;
end

spread_frames = 10;
spread_frames_filter = ones(1,spread_frames);

cue_trace = zeros(1,size(im.roi_trace_df,2));
cue_trace(round(bhv.cue_frames)) = 1;
cue_trace = conv(cue_trace,spread_frames_filter,'same');

reward_trace = zeros(1,size(im.roi_trace_df,2));
reward_trace(round(bhv.reward_frames)) = 1;
reward_trace = conv(reward_trace,spread_frames_filter,'same');

% load ROI labels and identify cell type
analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
roilabel_name = [animal '_roilabels.roilabel'];
load([analysis_path filesep roilabel_name],'-MAT')
cells = 1:length(roi_labels);
gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));

maxlags = 20;
cue_max_xcorr = zeros(size(im.roi_trace_df,1),1);
cue_max_xcorr_lag = zeros(size(im.roi_trace_df,1),1);
reward_max_xcorr = zeros(size(im.roi_trace_df,1),1);
reward_max_xcorr_lag = zeros(size(im.roi_trace_df,1),1);
for curr_cell = 1:size(im.roi_trace_df,1);
    [temp_xcorr temp_xcorr_lags] = ...
        xcorr(peak_matrix(curr_cell,:),cue_trace,maxlags);
    temp_xcorr = temp_xcorr(temp_xcorr_lags >= 0);
    temp_xcorr_lags = temp_xcorr_lags(temp_xcorr_lags >= 0);
    cue_max_xcorr(curr_cell) = max(temp_xcorr);
    cue_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
        find(temp_xcorr == max(temp_xcorr),1));
    
    [temp_xcorr temp_xcorr_lags] = ...
        xcorr(peak_matrix(curr_cell,:),reward_trace,maxlags);
    temp_xcorr = temp_xcorr(temp_xcorr_lags >= 0);
    temp_xcorr_lags = temp_xcorr_lags(temp_xcorr_lags >= 0);
    reward_max_xcorr(curr_cell) = max(temp_xcorr);
    reward_max_xcorr_lag(curr_cell) = temp_xcorr_lags( ...
        find(temp_xcorr == max(temp_xcorr),1));
    disp(curr_cell/size(im.roi_trace_df,1))
end

% get significance
num_shuff = 1000;
cue_shuffle = zeros(num_shuff,size(im.roi_trace_df,1));
for i = 1:num_shuff;
    curr_shuff = randperm(length(cue_trace));
    cue_shuffle(i,:) = peak_matrix*cue_trace(curr_shuff)';
    i
end
cue_errorbars = prctile(cue_shuffle,[0.05 99.95]);
sig_cue_cells = find(cue_max_xcorr > cue_errorbars(2,:)');

reward_shuffle = zeros(num_shuff,size(im.roi_trace_df,1));
for i = 1:num_shuff;
    curr_shuff = randperm(length(reward_trace));
    reward_shuffle(i,:) = peak_matrix*reward_trace(curr_shuff)';
    i
end
reward_errorbars = prctile(reward_shuffle,[0.05 99.95]);
sig_reward_cells = find(reward_max_xcorr > reward_errorbars(2,:)');

%% Threshold correlate

%%% Smooth, threshold trace
% smooth trace
smooth_size = 30;
roi_trace_df_smooth = nan(size(roi_trace_df));
for i = 1:size(roi_trace_df,1);
    roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),smooth_size,'loess');
end
% threshold trace
roi_trace_df_smooth_cutoff = zeros(size(roi_trace_df_smooth));
for curr_cell = 1:size(roi_trace_df_smooth,1)
    noise_est = abs(roi_trace_df(curr_cell,:) - roi_trace_df_smooth(curr_cell,:));
    cutoff = 4*mean(noise_est);
    curr_trace_cutoff = nan(size(roi_trace_df_smooth(curr_cell,:)));
    over_cutoff_indx = roi_trace_df_smooth(curr_cell,:) > cutoff;
    roi_trace_df_smooth_cutoff(curr_cell,over_cutoff_indx) = ...
        roi_trace_df_smooth(curr_cell,over_cutoff_indx);
end

%%% Find times of movement
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

movecorr = zeros(size(roi_trace_df_smooth),1);
for curr_cell = 1:size(roi_trace_df_smooth,1)
    tempcorr = corrcoef(roi_trace_df_smooth(curr_cell,:), ...
        +lever_active');
    movecorr(curr_cell) = tempcorr(2);
end

%% Try Friedrich deconvolution

roi_oopsi = zeros(size(roi_trace_df));
for curr_cell = 1:size(roi_trace_df,1)
    
    % normalized cutoff freq is fraction of nyquist
    butterworth_freq = (0.1*bhv.framerate)/(bhv.framerate/2);
    [b a] = butter(4, butterworth_freq);
    % using filtfilt instead of filter allows for zero-phase shift
    filt_trace = filtfilt(b,a,roi_trace_df(curr_cell,:));
    
    % deconvolve filtered trace
    V = struct;
    V.dt = 1/bhv.framerate;
    [n_best P_best V_est] = fast_oopsi(filt_trace,V);
    roi_oopsi(curr_cell,:) = n_best;
    disp(curr_cell/size(roi_trace_df,1))
end

% cut out frames around loop junctions (high fluor at new loop)
% cut out first 10 frames of every loop/8 files (fluorescence spike)
total_frames = cumsum(cellfun(@length, roi_trace_long_split));
loop_frames = [1;total_frames(8:8:end)];
loopcut_frames = repmat(loop_frames,1,15);
loopcut_addon = repmat([-5:9],size(loopcut_frames,1),1);
loopcut_frames = [loopcut_frames+loopcut_addon]';
loopcut_frames = loopcut_frames(:);
loopcut_frames(loopcut_frames < 1) = [];
loopcut_frames(loopcut_frames > size(roi_trace_long,2)) = [];
roi_oopsi(:,loopcut_frames) = 0;

% threshold oopsi
roi_oopsi_thresh = roi_oopsi;
roi_oopsi_thresh(roi_oopsi_thresh < 0.2) = 0;

% get TK coincidence index (digital)
num_frames = size(roi_trace_df,2);
roi_oopsi_bin = roi_oopsi_thresh;
roi_oopsi_bin(roi_oopsi_bin > 0) = 1;
tk_ci = zeros(size(roi_trace_df,1));

coincident_events = roi_oopsi_bin*roi_oopsi_bin';
chance_events = num_frames*(sum(roi_oopsi_bin,2)/num_frames)* ...
            (sum(roi_oopsi_bin,2)/num_frames)';
mean_events =   num_frames*sqrt((sum(roi_oopsi_bin,2)/num_frames)* ...
            (sum(roi_oopsi_bin,2)/num_frames)');  
tk_ci = (coincident_events - chance_events)./mean_events;        
        
% get TK coincidence index (analog)
num_frames = size(roi_trace_df,2);
tk_ci = zeros(size(roi_trace_df,1));

coincident_events = roi_oopsi*roi_oopsi';
chance_events = num_frames*(sum(roi_oopsi,2)/num_frames)* ...
            (sum(roi_oopsi,2)/num_frames)';
mean_events =   num_frames*sqrt((sum(roi_oopsi,2)/num_frames)* ...
            (sum(roi_oopsi,2)/num_frames)');  
tk_ci = (coincident_events - chance_events)./mean_events;        
                
        
        
%% XCorr response profile PCA

% define the cells to look at
cells = [1:length(grab_var.cue_tk_ci{1})];

if exist('roi_labels','var')
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';
    
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    gad_cells(ismember(gad_cells,contaminated)) = [];
else
    pyr_cells = cells;
    gad_cells = [];
end

cue_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.cue_tk_ci,'UniformOutput',false));
reward_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.reward_tk_ci,'UniformOutput',false));
move_start_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.move_start_tk_ci,'UniformOutput',false));
move_stop_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.move_stop_tk_ci,'UniformOutput',false));
move_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.move_tk_ci,'UniformOutput',false));
nomove_concat = cell2mat(cellfun(@(x) x(pyr_cells'),grab_var.nomove_tk_ci,'UniformOutput',false));

[coeffs_onset score_onset latent_onset] = princomp([cue_concat move_start_concat]);
[coeffs_offset score_offset latent_offset] = princomp([reward_concat move_stop_concat]);

[coeffs_offset score_offset latent_offset] = princomp([nomove_concat]);

[coeffs_all score_all latent_all] = princomp([cue_concat move_start_concat ...
    reward_concat move_stop_concat]);

% try kmeans correlation
all_concat = [cue_concat move_start_concat reward_concat ...
    move_stop_concat move_concat nomove_concat];
kmeans_concat = [cue_concat move_start_concat];
kmeans_idx = kmeans(kmeans_concat,8,'Distance','Correlation');
[temp_order idx] = sort(kmeans_idx);
temp = kmeans_concat(idx,:);
figure;imagesc(corrcoef(temp'));

% try normalizing, then doing stuff
curr_concat = [cue_concat move_start_concat ...
    reward_concat move_stop_concat];
curr_concat_norm = ...
     bsxfun(@times,curr_concat,1./sqrt(sum(curr_concat.^2,2)));
[coeffs_norm scores_norm latent_norm] = princomp(curr_concat_norm);




%% TTest behavior-related cells

% use thresholded oopsi trace
oopsi_cutoff = im.oopsi;
oopsi_cutoff(oopsi_cutoff < 0.2) = 0;

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

% classify movement/inactive times
% define velocity threshold for active movement
movethresh = 0.02;

% get velocity envelope, smooth over 5 frames (~150ms)
lever_force_resample_smooth = smooth(bhv.lever_force_resample,3);
lever_velocity = [0;diff(lever_force_resample_smooth)];
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

% define the cells to look at
cells = [1:size(roi_trace_df,1)];    

if exist('roi_labels','var')
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels))';
    
    pyr_cells(ismember(pyr_cells,contaminated)) = [];
    gad_cells(ismember(gad_cells,contaminated)) = [];
else
    disp('No ROI labels loaded')
    pyr_cells = cells;
    gad_cells = [];
end

frames_back = 0;
frames_forward = 5;
num_rois = size(roi_trace_df,1);
num_frames = size(roi_trace_df,2);
num_trials = length(good_trials);

bhv_activity_max = cell(2,1);
bhv_activity_mean = cell(2,1);
iti_activity_max = cell(2,1);
iti_activity_mean = cell(2,1);
bhv_activity_sig = cell(2,1);
% get mean/max activity for each trial
for behavior = 1:2
    switch behavior
        case 1
            bhv_frames = round(cue_frames_good);
        case 2
            bhv_frames = round(reward_frames_good);
    end
    bhv_activity_max{behavior} = zeros(num_rois,1);
    bhv_activity_mean{behavior} = zeros(num_rois,1);
    
    iti_activity_max{behavior} = zeros(num_rois,1);
    iti_activity_mean{behavior} = zeros(num_rois,1);
    
    movement_activity_max{behavior} = zeros(num_rois,1);
    movement_activity_mean{behavior} = zeros(num_rois,1);
    for i = 1:num_trials
        frames_search = bhv_frames(i):bhv_frames(i) + frames_forward;
        if max(frames_search) > size(roi_trace_df,2)
            continue
        end
        bhv_activity_max{behavior}(:,i) = ...
            max(oopsi_cutoff(:,frames_search),[],2);
        bhv_activity_mean{behavior}(:,i) = ...
            mean(oopsi_cutoff(:,frames_search),2);
        
        curr_trial_end = bhv.cue_frames(find(...
            bhv.cue_frames > bhv.cue_frames(i),1));
        curr_trial_unmoving = intersect(find(~lever_active), ...
            bhv_frames(i):curr_trial_end);
        curr_trial_moving = intersect(find(lever_active), ...
            bhv_frames(i):curr_trial_end);
        if isempty(curr_trial_unmoving) || isempty(curr_trial_moving)
            iti_activity_max{behavior}(:,i) = NaN;
            iti_activity_mean{behavior}(:,i) = NaN;
            
            movement_activity_max{behavior}(:,i) = NaN;
            movement_activity_mean{behavior}(:,i) = NaN;
            continue
        end
        iti_activity_max{behavior}(:,i) = ...
            max(oopsi_cutoff(:,curr_trial_unmoving),[],2);
        iti_activity_mean{behavior}(:,i) = ...
            mean(oopsi_cutoff(:,curr_trial_unmoving),2);
        
        movement_activity_max{behavior}(:,i) = ...
            max(oopsi_cutoff(:,curr_trial_moving),[],2);
        movement_activity_mean{behavior}(:,i) = ...
            mean(oopsi_cutoff(:,curr_trial_moving),2);
    end
    
    % find significantly modulated cells via bootstrap
    % max
    iti_prctile_max = ...
        prctile(bootstrp(1000,@nanmean,iti_activity_max{behavior}'),[0.05 99.95]);
    bhv_prctile_max = ...
        prctile(bootstrp(1000,@nanmean,bhv_activity_max{behavior}'),[0.05 99.95]);
    movement_prctile_max = ...
        prctile(bootstrp(1000,@nanmean,movement_activity_max{behavior}'),[0.05 99.95]);
    
    sig_up_cells_max{behavior} = find(bhv_prctile_max(1,:) > iti_prctile_max(2,:));
    sig_down_cells_max{behavior} = find(bhv_prctile_max(2,:) < iti_prctile_max(1,:));
    % mean
    iti_prctile_mean = ...
        prctile(bootstrp(1000,@nanmean,iti_activity_mean{behavior}'),[0.05 99.95]);
    bhv_prctile_mean = ...
        prctile(bootstrp(1000,@nanmean,bhv_activity_mean{behavior}'),[0.05 99.95]);
    movement_prctile_mean = ...
        prctile(bootstrp(1000,@nanmean,movement_activity_mean{behavior}'),[0.05 99.95]);
    
    sig_up_cells_mean{behavior} = find(bhv_prctile_mean(1,:) > iti_prctile_mean(2,:));
    sig_down_cells_mean{behavior} = find(bhv_prctile_mean(2,:) < iti_prctile_mean(1,:));
    
    if behavior == 2
        move_up_cells_max = find(movement_prctile_max(1,:) > iti_prctile_max(2,:));
        move_down_cells_max = find(movement_prctile_max(2,:) < iti_prctile_max(1,:));

        move_up_cells_mean = find(movement_prctile_mean(1,:) > iti_prctile_mean(2,:));
        move_down_cells_mean = find(movement_prctile_mean(2,:) < iti_prctile_mean(1,:));
    end
    
end

frames_surround = 200;
% plot the cue-related cells
cue_related_activity = cell(size(sig_up_cells_mean{1}));
for cell_idx = 1:length(cue_related_activity)
    curr_cell = sig_up_cells_mean{1}(cell_idx);
    % plot around cue
    cue_reward = nan(length(bhv.cue_frames),frames_surround*2+1);
    for i = 1:length(bhv.cue_frames);
        curr_reward = round(bhv.cue_frames(i));
        roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
            continue
        end
        cue_reward(i,:) = roi_trace_df(curr_cell,roi_reward_frames);
    end
    cue_related_activity{cell_idx} = cue_reward;
end
figure;imagesc([cue_related_activity{:}]);
title('cue-related')

% plot the reward-related cells
reward_related_activity = cell(size(sig_up_cells_mean{2}));
for cell_idx = 1:length(reward_related_activity)
    curr_cell = sig_up_cells_mean{2}(cell_idx);
    % plot around reward
    reward_reward = nan(length(bhv.reward_frames),frames_surround*2+1);
    for i = 1:length(bhv.reward_frames);
        curr_reward = round(bhv.reward_frames(i));
        roi_reward_frames = curr_reward-frames_surround:curr_reward+frames_surround;
        if any(roi_reward_frames<=0) || any(roi_reward_frames>size(roi_trace_df,2));
            continue
        end
        reward_reward(i,:) = roi_trace_df(curr_cell,roi_reward_frames);
    end
    reward_related_activity{cell_idx} = reward_reward;
end
figure;imagesc([reward_related_activity{:}]);
title('reward-related')




%% P-value for movement relatedness

animal = 'AP76';
day = '120903';

move_trial_mean = [];
still_trial_mean = [];

% Break up baseline corrected trace into loops
num_loops = length(im.roi_trace_split)/8;
file_lengths = cellfun(@length,im.roi_trace_split);
file_lengths_loops = reshape(file_lengths,8,num_loops);
loop_lengths = sum(file_lengths_loops);
%oopsi_split = mat2cell(im.oopsi,size(im.oopsi,1), ...
%    loop_lengths);

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
        zeros(1,loop_lengths(curr_xsg));
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
    
    % apparently the idea for trial analysis was a miscommunication
    
%     % get cue times in current loop
%     clear cue_sample cue_frames
%     for curr_trial_idx = 1:length(curr_trial_list(:,2));
%         curr_trial = curr_trial_list(curr_trial_idx,2);
% 
%         % get trial start in dispatcher
%         curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
%         % get cue and reward sample in xsg (downsampled to ms)
%         curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
%             curr_bhv_start)*(xsg_samplerate/10);
%         cue_sample(curr_trial_idx) = ...
%             round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
%             *(xsg_samplerate/10) + ...
%             curr_bhv_cue_sample_rel);
%         curr_bhv_cue_frame_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
%             curr_bhv_start)*(bhv.framerate);
%         cue_frames(curr_trial_idx) = ...
%             round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
%             *(bhv.framerate) + ...
%             curr_bhv_cue_frame_rel);
%     end
%     
%     for curr_trial_idx = 1:length(curr_trial_list(:,2))-1;
%         curr_trial = curr_trial_list(curr_trial_idx,2);
%         curr_trial_frames  = cue_frames(curr_trial_idx): ...
%             cue_frames(curr_trial_idx+1);
%         curr_trial_move_frames = intersect(curr_trial_frames,move_frames);
%         curr_trial_still_frames = intersect(curr_trial_frames,still_frames);
%         
%         move_trial_mean(:,curr_trial) = sum(oopsi_split{curr_xsg} ...
%             (:,curr_trial_move_frames),2);
%         still_trial_mean(:,curr_trial) = sum(oopsi_split{curr_xsg} ...
%             (:,curr_trial_still_frames),2);
%     end
end

oopsi_thresh = im.roi_trace_df;
oopsi_thresh(oopsi_thresh < 0.0) = 0;

% find significantly modulated cells via shuffling
num_rois = size(im.roi_trace_df,1);
movement_epochs_full = [movement_epochs{:}];
move_deconv = oopsi_thresh*[lever_active_frames{:}]';

num_rep = 1000;
move_deconv_shuffle = zeros(num_rois,num_rep);
for i = 1:num_rep
    curr_perm = randperm(length(movement_epochs_full));
    move_deconv_shuffle(:,i) = oopsi_thresh*[movement_epochs_full{curr_perm}]';
    i
end
move_deconv_prctile = prctile(move_deconv_shuffle',[0.5 99.5]);
move_cells = find(move_deconv' > move_deconv_prctile(2,:));
still_cells = find(move_deconv' < move_deconv_prctile(1,:));

% still_prctile_mean = ...
%     prctile(bootstrp(1000,@nanmean,still_trial_mean'),[0.05 99.95]);
% move_prctile_mean = ...
%     prctile(bootstrp(1000,@nanmean,move_trial_mean'),[2.5 97.5]);
% 
% sig_move_cells = find(move_prctile_mean(1,:) > still_prctile_mean(2,:));
% sig_still_cells = find(move_prctile_mean(2,:) < still_prctile_mean(1,:));















%% Do trial-by-trial analysis

% split up roi_trace_df into roi_trace_long_split*8 sizes
loop_size = cellfun(@length,im.roi_trace_split);
loop_size_sum = cumsum(loop_size);
loop_frames = [0;loop_size_sum(8:8:end-7)];
% roi_trace_df_loops = mat2cell(roi_trace_df,1,loop_frames);

% for given animal, figure out which trials are cued/uncued rewarded, which
% trials are fully imaged and rewarded
animal = 'AP79';
curr_day = 5;

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

day = num2str(days{5});

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
    
    [lever_active lever_force_smooth] = AP_parseLeverMovement(xsg_data);
    
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
                length(lever_active) < cue_sample(curr_trial_idx+1)
            continue
        end
        
        % ignore if out of bounds
        if cue_sample(curr_trial_idx+1) > length( ...
                lever_active)
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
    butterworth_freq = 5/(bhv.framerate/2);
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

% 
% % get correlations between gad and pyramidal cells
% max_corrcoef = corrcoef(cell_trial_max');
% 
% % pull out the traces for cued-rewarded trials
% trial_traces = cell(num_trials,1);
% for i = 1:num_trials-1;
%     if ~isnan(cue_frames(i)) && ~isnan(cue_frames(i+1))
%         trial_traces{i} = roi_trace_df(:,cue_frames(i):cue_frames(i+1));
%     end
% end



%% Get length constant for modulated/unmodulated cells

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
binsize = 1:50:501;
num_bins = length(binsize);

mod_cells = grab_var.move_cells{curr_day};
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







%% Playground

spread_frames = 30;
spread_frames_filter = ones(1,spread_frames);
peak_matrix_conv = conv2(peak_matrix,spread_frames_filter,'same');

peak_active = sum(peak_matrix_conv(pyr_cells,:)) > 0;

active_starts = find(diff([0 peak_active 0]) == 1);
active_stops = find(diff([0 peak_active 0]) == -1);

peak_cells = cell(length(active_starts),1);
cell_active_matrix = zeros(length(active_starts),size(im.roi_trace_df,1));
for i = 1:length(active_starts);
    peak_cells{i} = ...
        find(any(peak_matrix_conv(:,active_starts(i):active_stops(i)),2)); 
    cell_active_matrix(i,peak_cells{i}) = 1;
end








