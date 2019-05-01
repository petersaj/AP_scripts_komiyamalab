function [roi_trace,roi_trace_bg] = AP_getConcatTrace_continuous_batch_fastz(roi_filenames,roi_slices,tiff_path,interp_sub,local_comp)
%[roi_trace,roi_trace_bg] = AP_getConcatTrace_continuous_batch_fastz(roi_filenames,roi_slices,tiff_path,interp_sub,local_comp)
%
% Get fluorescence traces of ROIs (and background ROIs, if available) of
% movies which were taken continuously. Assumes new TIFF header format
% where frame number is logged.
%
% This is specifically for fast-z
%
% roi_filenames - ROI files for slices (any number)
%
% roi_slices - the slices which the ROI files correspond to (in order)
%
% interp_sub - get fluorescence interpolated across ROIs as background
% (false by default)
%
% local_comp - if being done locally, copy file to local drive temporarily
% and load/process (true by default if on PC)

% If ROI filenames are not specified, query
if ~exist('roi_filenames','var') || isempty(roi_filenames)
    
    [roi_files,roi_path] = uigetfile('*.roi','Choose ROI files','multiselect','on');    
    roi_files = sort(roi_files);
    roi_filenames = cellfun(@(x) [roi_path x],roi_files,'uni',false);
    
    roi_slices = str2num(cell2mat((inputdlg('Slices which ROIs correspond to:'))));
    
end

% If TIFF path is not specified, query
if ~exist('tiff_path','var') || isempty(tiff_path)
    [tiff_path] = uigetdir('Choose TIFF path');    
end

% Get TIFF filenames
tiff_dir = dir([tiff_path filesep '*.tif']);
tiff_files = sort({tiff_dir.name});

% append directory to tiff filename
tiff_filenames = cellfun(@(x) [tiff_path filesep x], ...
    tiff_files,'UniformOutput',false);

% if on local computer, create directory for temporary files
if ~exist('local_comp','var') && ispc 
    local_comp = true;
elseif ~exist('local_comp','var') && ~ispc
    local_comp = false;
end

if local_comp
    local_dir = 'C:\temp_files';
    local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
    if ~isdir(local_dir)
        mkdir(local_dir)
    end   
    % delete local temp file if it exists
    if exist(local_filename)
        delete(local_filename)
    end
end

if ~exist('interp_sub','var');
    interp_sub = false;
end

% load ROIs
rois = cell(length(roi_filenames),1);
for curr_roi_file = 1:length(roi_filenames)
    rois{curr_roi_file} = load(roi_filenames{curr_roi_file},'-MAT');
end

% Get info from first file
imageinfo = imfinfo(tiff_filenames{1},'tiff');
[img_parameter img_value] = strread( ...
    imageinfo(1).ImageDescription,'%s %s', 'delimiter','=\n');

% Dimensions of movie
M = imageinfo(1).Height;
N = imageinfo(1).Width;

% Find number of volumes
numvolumes_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.fastZNumVolumes')),img_parameter,'UniformOutput',0));
numvolumes = str2num(img_value{numvolumes_indx});

% get number of slices
numslices_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.stackNumSlices')),img_parameter,'UniformOutput',0));
numslices = str2num(img_value{numslices_indx});

% use frame tag to check if slices discarded in motion correction: get the
% frame tags from the first numslices frames and use the unique numbers
% equal or less than numslices
numslices_frame_tag = nan(numslices,1);
for curr_frame_tag = 1:numslices
    [img_parameter img_value] = strread( ...
        imageinfo(curr_frame_tag).ImageDescription,'%s %s', 'delimiter','=\n');
    frame_tag_idx = cellfun(@(x) ~isempty(strfind(x,'Frame Tag')),img_parameter);
    numslices_frame_tag(curr_frame_tag) = str2num(cell2mat(img_value(frame_tag_idx)));
end
saved_slices = numslices_frame_tag(numslices_frame_tag <= numslices);
num_saved_slices = length(saved_slices);

% create masks from ROIs
cellMask = cell(length(roi_filenames),1);
for curr_roi_file = 1:length(roi_filenames)
    cellMask{curr_roi_file} = false(N*M,length(rois{curr_roi_file}.polygon.ROI));
    for n_polygon = 1:length(rois{curr_roi_file}.polygon.ROI);
        if ~isfield(rois{curr_roi_file}.polygon,'autosort');
            temp_mask = poly2mask(rois{curr_roi_file}.polygon.ROI{n_polygon}(:,1)', ...
                rois{curr_roi_file}.polygon.ROI{n_polygon}(:,2)',N,M);
            cellMask{curr_roi_file}(:,n_polygon) = temp_mask(:);
        elseif isfield(rois{curr_roi_file}.polygon,'autosort');
            temp_mask = rois{curr_roi_file}.polygon.autosort(:,:,n_polygon);
            cellMask{curr_roi_file}(:,n_polygon) = temp_mask(:);
        end
    end
    % Normalize mask values to prepare for matrix multiplation
    cellMask_norm{curr_roi_file} = bsxfun(@times,cellMask{curr_roi_file},1./sum(cellMask{curr_roi_file},1));
    
    % Grab mask for background ROIs, if available
    if isfield(rois{curr_roi_file}.polygon,'bgROI')
        for n_polygon = 1:length(rois{curr_roi_file}.polygon.bgROI);
            cellMask_bg{curr_roi_file}(:,n_polygon) = rois{curr_roi_file}.polygon.bgROI_mask{n_polygon}(:);
        end
        % Normalize mask values to prepare for matrix multiplation
        cellMask_bg_norm{curr_roi_file} = bsxfun(@times,cellMask_bg{curr_roi_file},1./sum(cellMask_bg{curr_roi_file},1));
    end
end

% get the traces for the ROIs
disp('Getting and concatenating traces');
roi_trace = cellfun(@(x) nan(length(x.polygon.ROI),numvolumes),rois,'uni',false);
roi_trace_bg = cellfun(@(x) nan(length(x.polygon.ROI),numvolumes),rois,'uni',false);

for i = 1:length(tiff_filenames);
        
    % if on local computer, copy image to temporary file, otherwise load
    if ~local_comp
        img_filename = [tiff_filenames{i}];
    else
        disp('Copying movie to local drive...')
        tic
        copyfile(tiff_filenames{i},local_filename);
        toc
        disp('Done.')
        img_filename = local_filename;
    end
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % get number of each frame from header
    [img_parameter img_value] = arrayfun(@(x) strread( ...
        imageinfo(x).ImageDescription,'%s %s', 'delimiter','=\n'),1:length(imageinfo),'uni',false);
    frame_tag_idx = cellfun(@(x) cellfun(@(x) ...
        ~isempty(strfind(x,'Frame Tag')),x),img_parameter,'uni',false);
    frame_tag = cellfun(@(value,idx) str2num(value{idx}),img_value,frame_tag_idx);
    
    % correct frame tag by saved slices
    frame_tag_slice = ceil(frame_tag./numslices);
        
    % Load and get activity frame-by-frame to allow for any computer to do
    % many of these simultanously (as opposed to loading in whole movie,
    % which is probably faster)
    tic
    disp(['Loading file, getting activity (' img_filename ')']);
    if interp_sub
        for curr_roi_file = 1:length(roi_filenames)
            allROI_mask{curr_roi_file} = reshape(any(cellMask{curr_roi_file},2),M,N);
        end
    end
    
    load_slice_frames = cell(length(roi_filenames),1);
    for curr_roi_file = 1:length(roi_filenames);
        load_slice_frames{curr_roi_file} = roi_slices(curr_roi_file):num_saved_slices:numframes;
    end
    
    for curr_roi_file = 1:length(roi_filenames);
        for loadframe_idx = 1:length(load_slice_frames{curr_roi_file});
            loadframe = load_slice_frames{curr_roi_file}(loadframe_idx);
            curr_volume = frame_tag_slice(loadframe);
            
            % Load current frame
            curr_frame = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
            
            % Get activity from current frame
            roi_trace{curr_roi_file}(:,curr_volume) = transpose(double(curr_frame(:)')*cellMask_norm{curr_roi_file});
            
            % Get background activity from current frame if selected
            if interp_sub
                curr_frame_interp = double(curr_frame);
                curr_frame_interp(allROI_mask{curr_roi_file}) = NaN;
                curr_frame_interp = inpaint_nans(curr_frame_interp,4);
                roi_trace_bg{curr_roi_file}(:,curr_volume) = ...
                    transpose(curr_frame_interp(:)'*cellMask_norm{curr_roi_file});
            elseif isfield(rois{curr_roi_file}.polygon,'bgROI')
                roi_trace_bg{curr_roi_file}(:,curr_volume) = ...
                    transpose(double(curr_frame(:)')*cellMask_bg_norm{curr_roi_file});
            end
        end
    end
    disp('Done')
    toc
    
    % If on local computer, get rid of temporary file
    if local_comp
        delete(local_filename);
    end
end

disp('Finished getting ROI traces');

