function bhv = AP_getLeverBehavior_continuous(animal,data_path)
% bhv = AP_getLeverBehavior_continuous(animal,data_path)
%
% XSG folder and behavior file must be in data_path
% Animal must be named by first two initials then number (i.e. AP101)

%% Get imaging framerate

% Get framerate (pick first image file, read header)
dir_currfolder = dir(data_path);
dir_filenames = {dir_currfolder.name};
tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
img_filename = dir_filenames(tiff_filename_indx);

% If there are no image files, get out only behavior data
if isempty(img_filename)
    bhv_only = true;
elseif ~isempty(img_filename)
    bhv_only = false;
end

if ~bhv_only
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    
    % save framerate for output
    bhv.framerate = framerate;
    
end

%% Get XSG data

% Get XSG filenames, load XSG data
xsg_path = [data_path filesep animal(1:2) repmat('0',1,4-length(animal(3:end))) animal(3:end)];
xsg_data = AP_load_xsg_continuous(xsg_path);

% save xsg sample rate (unfortunately not saved in continuous)
bhv.xsg_sample_rate = 10000;

if ~bhv_only
    
    % get frame times (in seconds) from frame trigger trace
    frame_channel = cellfun(@(x) strcmp(x,'Frame'),xsg_data.channel_names);
    frame_trace = xsg_data.channels(:,frame_channel);
    frame_times = (find(frame_trace(2:end) > 2.5 & ...
        frame_trace(1:end-1) < 2.5) + 1)/bhv.xsg_sample_rate;
    
    numframes = length(frame_times);
    
    bhv.frame_times = frame_times;
    
    % get lever duing imaging, downsample, store in structure
    median_frame_time = median(diff(bhv.frame_times));
    
end

% if touch sensor was recorded, save
touch_channel = cellfun(@(x) strcmp(x,'Lever_touch'),xsg_data.channel_names);
if any(touch_channel)
   bhv.lever_touch = xsg_data.channels(:,touch_channel); 
end

% if optogenetic TTL was recorded, save
opto_channel = cellfun(@(x) strcmp(x,'Opto'),xsg_data.channel_names);
if any(opto_channel)
   bhv.opto = xsg_data.channels(:,opto_channel); 
end

% --> the round function shouldn't be necessary here, but sometimes an
% interger is not recorded as an interger for unknown reason

lever_channel = cellfun(@(x) strcmp(x,'Lever'),xsg_data.channel_names);
lever_force_raw = xsg_data.channels(:,lever_channel);

bhv.lever_force = lever_force_raw;

% If number of frames provided:
% - resample lever to be in kHz
% - resample lever force to be 1:1 with numframes for plotting
if ~bhv_only

    lever_samples_imaged = round([bhv.frame_times(1)*bhv.xsg_sample_rate: ...
        (bhv.frame_times(numframes) + median_frame_time)*bhv.xsg_sample_rate]);
    
    lever_force_imaged = lever_force_raw(lever_samples_imaged);
    
    lever_force_downsample = resample(lever_force_imaged,1,10);
    
    bhv.imaged_downsampled_lever_force = lever_force_downsample;
    
    lever_numframes = numframes;
    if mod(lever_numframes,2) ~= 0
        lever_numframes = lever_numframes - 1;
    end
    % resample lever force to so # samples = # frames
    [n d] = rat(lever_numframes/length(lever_samples_imaged));
    lever_force_resample = resample(lever_force_imaged,n,d);
    % add or delete last sample if necessary
    if length(lever_force_resample) < lever_numframes
        lever_force_resample(end+1:lever_numframes) = 0;
    elseif length(lever_force_resample) > lever_numframes
        lever_force_resample = lever_force_resample(1:lever_numframes);
    end
    
    % save resampled lever force for output
    bhv.lever_force_plot = lever_force_resample;
    
end

%% Get Dispatcher data

% Get behavior filename
dir_currfolder = dir(data_path);
dir_filenames = {dir_currfolder.name};
bhv_file = strncmp('data_@',dir_filenames,6);
bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
bhv_fullfilename = [data_path filesep bhv_filename];

% Load behavior
warning off
load(bhv_fullfilename,'-MAT');
warning on

% Store raw behavior
bhv.bhv_times = saved_history.ProtocolsSection_parsed_events;

% Get trial times and numbers
trial_channel = cellfun(@(x) strcmp(x,'Trial_number'),xsg_data.channel_names);
trial_list = AP_ReadBitCode(xsg_data.channels(:,trial_channel));
bhv.trial_list = trial_list;

% Convert behavior events into frames, pull out frame/imaged trial info
% if no dispatcher file (because I was an idiot and didn't save) then skip
if ~isempty(bhv_filename) && ~bhv_only
    
    % Convert behavior to frames
    [bhv_frames,imaged_trials,frame_times] = ...
        AP_dispatcher2frames_continuous(bhv_fullfilename,xsg_data);
    
    % if there are catch trials, get them   
    if isfield(saved_history,'CatchSection_catch_trial')
        catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
        % there's one or more uncompleted trials at the end
        bhv.catch_trials = logical(catch_trials(1:length(bhv_frames)));
    end
    
    % save frame bhv for output
    bhv.bhv_frames = bhv_frames;
    bhv.imaged_trials = imaged_trials;
    bhv.frame_times = frame_times;
    
end
