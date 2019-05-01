function [bhv_frames imaged_trials xsg_trials] = AP_getBehavior_dispatcher(framerate,bhv_fullfilename,xsg_fullfilenames)
% Convert all dispatcher behavior into frames
% [bhv_frames imaged_trials xsg_trials] = AP_getBehavior_dispatcher(framerate,bhv_fullfilename,xsg_fullfilenames);
%
% INPUT
% -framerate (in hz)
% 
% - bhv is optional filenames of behavior file, will be prompted to select
% if not entered
%
% - xsg_files is filenames (in a cell) of xsg files, will be prompted to select
% if not entered
%
% OUTPUT
% - bhv_frames is the behavior file from dispatcher, with all times converted
% into frames since beginning of that acquisition
%
% - imaged_trials is a logical index of trials which were imaged
%
% - acq_trials is the xsg file which that trial was a part of (from the order
% in which they were selected)

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

if ~exist('xsg_fullfilenames','var')
    [xsg_filename_all xsg_path] = uigetfile('*.xsg','Pick xsg files','Multiselect','on');
    xsg_filename_all = sort(xsg_filename_all);
    xsg_fullfilenames = xsg_filename_all;
    xsg_fullfilenames = cellfun(@(x) [xsg_path x],xsg_fullfilenames,'UniformOutput',false);
end

% Initialize matricies for imaged trials and acquisition trials
imaged_trials = false(length(bhv),1);
xsg_trials = zeros(length(bhv),1);

% Go through all imaged trials, convert behavior times to frames
for i = 1:length(xsg_fullfilenames);
    xsg = load(xsg_fullfilenames{i},'-MAT');
    
    % Get all trial offsets in samples
    xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
    
    % Get trials in raw samples since started
    curr_trial_list = [];
    curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{i});
    
    % loop through those trials, find the offsets
    for curr_trial_indx = 1:size(curr_trial_list,1);    
        
        curr_trial = curr_trial_list(curr_trial_indx,2);
        % skip if it's the last trial and not completed in behavior
        if curr_trial > length(bhv) || curr_trial < 1
            continue
        end
        imaged_trials(curr_trial) = true;
        xsg_trials(curr_trial) = i;
        % the start time is the rise of the first bitcode
        curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
        curr_xsg_bhv_offset = curr_bhv_start - curr_trial_list(curr_trial_indx,1);
        % apply the offset to all numbers within the trial 
        % Find all fields in overall structure of trial
        curr_fieldnames = [];
        curr_fieldnames = fieldnames(bhv_frames{curr_trial});
        
        % make sure that bhv_frames is identical for this trial first
        bhv_frames{curr_trial} = bhv{curr_trial};
        
        for curr_field = 1:length(curr_fieldnames)
            % a) get subfields
            curr_subfields = [];
            curr_subfields = fieldnames(bhv_frames{curr_trial}.(curr_fieldnames{curr_field}));
            % b) find which subfields are numeric
            curr_numeric_subfields = [];
            curr_numeric_subfields = structfun(@isnumeric,bhv_frames{curr_trial}.(curr_fieldnames{curr_field}));
            % c) subtract offset from numeric fields and convert to frames
            for subfield_offset_fix = find(curr_numeric_subfields)';
                % get current values of the subfield
                curr_bhv_times = [];
                curr_bhv_times = bhv_frames{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix});
                % compensate for offset
                curr_bhv_times = curr_bhv_times - curr_xsg_bhv_offset;
                % convert to frames (convert xsg sample rate to ms)
                curr_bhv_frames = [];
                % this was a mix up, converting seconds -> frames so fixed
                curr_bhv_frames = framerate*curr_bhv_times;%(curr_bhv_times/(xsg_sample_rate/1000));
                bhv_frames{curr_trial}.(curr_fieldnames{curr_field}).(curr_subfields{subfield_offset_fix}) = ...
                    curr_bhv_frames;
            end
        end           
    end   
end











