function AP_corticospinal_process_raw(animal,bg_flag,overwrite)
% AP_corticospinal_process_raw(animal,bg_flag,overwrite)
%
% From the ground up: create fixed background ROI files, create
% concatenated traces, create analysis files with (corrected) df/f traces
%
% bg_flag = background subtraction
%   CAN BE EITHER: 
% 0/false - no background
% 1 - annulus ROI background subtraction
% 2 - interpolated background subtraction
%
% overwrite = overwrite existing analysis files (true by default)

%% Check if this is being done locally (PC) or on compute (Linux)
    
local_comp = ispc;
    
%% Set default overwrite
if nargin < 3
   overwrite = true; 
end

%% Set paths
if ~local_comp
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
else
    data_path = ['Z:\People\Andy\Data\' animal];
end

roi_path = uigetdir(data_path,'Choose ROI path');

%roi_path = [data_path filesep animal '_batch_thresh_roi'];

if bg_flag == 1
    bg_roi_path = [roi_path '_bg'];
    use_roi_path = bg_roi_path;
else
    use_roi_path = roi_path;
end


%% Copy over XSG and behavior files from raw to processed data folder

if ~local_comp
    raw_path = '/usr/local/lab/Data/ImagingRig3';
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
else
    raw_path = 'Z:\Data\ImagingRig3';
    data_path = ['Z:\People\Andy\Data\' animal];
end


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
    
    day = num2str(days{curr_day});
    curr_source = [raw_path filesep day filesep animal];
    curr_destination = [data_days{curr_day}];
    
    % Identify non-tiff files / subfolders
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    nontiff_dir = cellfun(@(x) isempty(strfind(x,'.tif')),dir_filenames);
    nontiff_names = dir_filenames(nontiff_dir);
    
    % copy all the non-tiffs over (first two = ., .., start at 3)
    for i = 3:length(nontiff_names)
        if ~exist([curr_destination filesep nontiff_names{i}])
            copyfile([curr_source filesep nontiff_names{i}], ...
                [curr_destination filesep nontiff_names{i}]);
        end
    end
    disp(['Copied extras: ' num2str(curr_day) '/' num2str(length(days))])
    
end


%% Create annulus background ROIs (if flagged)

if bg_flag == 1    
    AP_make_bgroi(data_path,roi_path);    
end

%% Get concatenated traces, make analysis files

% Get ROI filenames
roi_dir = dir([use_roi_path filesep '*.roi']);
roi_files = sort({roi_dir.name});
roi_filenames = cellfun(@(x) [use_roi_path filesep x],roi_files,'uni',false);

for curr_session = 1:length(roi_filenames)
    
    analysis_filename = [use_roi_path filesep ...
        roi_files{curr_session}(1:end-4) '_analysis.mat'];
    
    % Don't process if already written & overwrite turned off
    if ~overwrite && exist(analysis_filename,'file')
        continue
    end
       
    % Get day and corresponding TIFF path
    day = roi_files{curr_session}(1:6);
    tiff_path = [data_path filesep day];
    
    % Get concatenated trace
    [roi_trace,roi_trace_bg] = ...
        AP_getConcatTrace_continuous_batch(roi_filenames{curr_session},tiff_path,bg_flag == 2,local_comp);
    
    % Get behavior
    bhv = AP_getLeverBehavior_continuous(animal,tiff_path);
    
    % Make df/f trace
    if bg_flag
        [~,roi_trace_bg_baseline] = ...
            AP_baselineEstimation(roi_trace_bg,bhv.framerate);
        roi_trace_bgsubtract = roi_trace - ...
            (roi_trace_bg - roi_trace_bg_baseline);
              
        [~,roi_trace_baseline] = ...
            AP_baselineEstimation(roi_trace,bhv.framerate);
        
%         % To prevent over-subtracting, rectify at baseline
%         roi_trace_bgsubtract = roi_trace - ...
%             min((roi_trace_bg - roi_trace_bg_baseline), ...
%             (roi_trace - roi_trace_baseline));
        
        roi_trace_df = ...
            (roi_trace_bgsubtract - roi_trace_baseline)./roi_trace_baseline;
    else
        roi_trace_df = ...
            AP_baselineEstimation(roi_trace,bhv.framerate);
    end
    
    im.roi_trace_df = roi_trace_df;
    im.roi_trace = roi_trace;
    if bg_flag
        im.roi_trace_bg = roi_trace_bg;
    end
    
    save(analysis_filename,'bhv','im');
    
    disp(['Finished analyzing: ' animal ' day ' day])
    
end

disp(['Finished processing ' animal]);







