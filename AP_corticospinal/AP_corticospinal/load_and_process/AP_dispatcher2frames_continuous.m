function [bhv_frames,imaged_trials,frame_times] = AP_dispatcher2frames_continuous(bhv_fullfilename,xsg_data)
% Convert all dispatcher behavior into frames
% [bhv_frames imaged_trials xsg_trials] = AP_dispatcher2frames_continuous(framerate,bhv_fullfilename,xsg_fullfilenames);
%
% INPUT
% -framerate (in hz)
% 
% - bhv is filenames of behavior file, will be prompted to select
% if not entered
%
% - xsg_data is raw continuous xsg data as read by AP_load_xsg_continuous
%
% OUTPUT
% - bhv_frames is the behavior file from dispatcher, with all times converted
% into frames since beginning of that acquisition
%
% - imaged_trials is a logical index of trials which were imaged
%
% - frame_times: the times (in seconds) of each frame

if ~exist('bhv_fullfilename','var')
    [bhv_filename bhv_path] = uigetfile('*.mat','Pick behavior file','Multiselect','off');
    bhv_fullfilename = [bhv_path bhv_filename];
end
% Turn warning off for loading, otherwise get message about
% StateMachineAssembler class
warning off
load(bhv_fullfilename,'-MAT');
warning on
bhv = saved_history.ProtocolsSection_parsed_events;
bhv_frames = bhv;

% Go through all imaged trials, convert behavior times to frames
imaged_trials = false(length(bhv),1);

% Get all trial offsets in samples - has to be hardcoded because not
% recorded in continuous XSG raw files
xsg_sample_rate = 10000;

% Get frame times (in seconds) from frame trigger trace
frame_channel = cellfun(@(x) strcmp(x,'Frame'),xsg_data.channel_names);
frame_trace = xsg_data.channels(:,frame_channel);
frame_times = (find(frame_trace(2:end) > 2.5 & ...
    frame_trace(1:end-1) < 2.5) + 1)/xsg_sample_rate;

% Get trials in raw samples since started
trial_channel = cellfun(@(x) strcmp(x,'Trial_number'),xsg_data.channel_names);
curr_trial_list = AP_ReadBitCode(xsg_data.channels(:,trial_channel));

% loop through those trials, find the offsets
for curr_trial_indx = 1:size(curr_trial_list,1);
    
    curr_trial = curr_trial_list(curr_trial_indx,2);
    
    % skip if it's the last trial and not completed in behavior
    if curr_trial > length(bhv) || curr_trial < 1
        continue
    end
    
    imaged_trials(curr_trial) = true;
    
    % the start time is the rise of the first bitcode
    curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
    curr_xsg_bhv_offset = curr_bhv_start - curr_trial_list(curr_trial_indx,1);
    % apply the offset to all numbers within the trial
    % Find all fields in overall structure of trial
    curr_fieldnames = fieldnames(bhv_frames{curr_trial});
    
    for curr_field = 1:length(curr_fieldnames)
        % a) get subfields
        curr_subfields = fieldnames(bhv_frames{curr_trial}.(curr_fieldnames{curr_field}));
        % b) find which subfields are numeric
        curr_numeric_subfields = structfun(@isnumeric,bhv_frames{curr_trial}.(curr_fieldnames{curr_field}));
        % c) subtract offset from numeric fields and convert to frames
        for subfield_offset_fix = find(curr_numeric_subfields)';
            % get current values of the subfield
            curr_bhv_times = bhv_frames{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix});
            % compensate for offset
            curr_bhv_times = curr_bhv_times - curr_xsg_bhv_offset;
            % convert to frames (get the closest frame from frame times)
            curr_bhv_frames = nan(size(curr_bhv_times));
            for curr_el = 1:numel(curr_bhv_times)
                [~,curr_bhv_frames(curr_el)] = min(abs(frame_times-curr_bhv_times(curr_el)));
            end
            bhv_frames{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix}) = ...
                curr_bhv_frames;
        end
    end
end











