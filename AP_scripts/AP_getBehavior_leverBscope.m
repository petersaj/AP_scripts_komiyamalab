function bhv = AP_getBehavior_leverBscope(animal,day,roi_trace_long)
% bhv = AP_getBehavior_leverBscope(animal,day_roi_trace_long)


%% Get trace and behavior times

data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];

disp('Getting behavior:')
disp('Specific to Andy')

disp('Getting filenames')

% Get behavior filename
dir_currfolder = dir(data_path);
dir_filenames = {dir_currfolder.name};
bhv_file = strncmp('data_@',dir_filenames,6);
bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
bhv_fullfilename = [data_path filesep bhv_filename];

% Get XSG filenames
xsg_path = [data_path filesep 'AP' repmat('0',1,4-length(animal(3:end))) animal(3:end)];
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

% save framerate for output
bhv.framerate = framerate;

if ~isempty(bhv_filename)
    
    % Convert behavior to frames
    disp('Converting behavior times to image frames')
    [bhv_frames imaged_trials xsg_trials] = AP_getBehavior_dispatcher(framerate,bhv_fullfilename,xsg_fullfilenames);
    
    % if there are catch trials, get them
    warning off
    load(bhv_fullfilename,'-MAT');
    warning on
    
    if isfield(saved_history,'CatchSection_catch_trial')
        catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
        % there's one or more uncompleted trials at the end
        bhv.catch_trials = catch_trials(1:length(bhv_frames));
    end
    
    % save frame bhv for output
    bhv.bhv_frames = bhv_frames;
    bhv.imaged_trials = imaged_trials;
    bhv.xsg_trials = xsg_trials;
    
    num_trials = length(bhv_frames);
    
else
    % if no behavior file
    disp([animal ' day ' day ' has no bhv file']);
    num_trials = 0;
    imaged_trials = [];
    xsg_trials = [];
    bhv_frames = [];
end



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

if strmatch(day,'121001')
    last_frames(:) = 500;
end

% Sanity check that the number of frames matches calculated number
if exist('roi_trace_long','var')
    total_est_frames = 500*7*size(img_filename,2)/8+sum(last_frames);
    if size(roi_trace_long,2) ~= total_est_frames
        display('Incorrect number of total estimated frames')
        keyboard
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

% save add_frames for output
bhv.add_frames = add_frames;

% get relevant behaviors, add appropriate number of frames for concat
if ~isempty(bhv_filename) && isfield('reward',bhv_frames{2}.states)
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
end

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

% save lever force for output
bhv.lever_force = lever_force;

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

% save resampled lever force for output
bhv.lever_force_resample = lever_force_resample;

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

if ~isempty(bhv_filename) && isfield('reward',bhv_frames{2}.states)
    
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
    
    % save these behaviors for output
    bhv.cue_frames = cue_frames;
    bhv.reward_frames = reward_frames;
    bhv.lick_frames = lick_frames;
    
end
% Prepare bhv_frames by adding previous xsg frames for concatenation

disp('Finished getting behavior.')