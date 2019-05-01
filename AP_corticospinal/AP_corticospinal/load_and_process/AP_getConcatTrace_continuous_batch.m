function [roi_trace,roi_trace_bg] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp)
%[roi_trace,roi_trace_bg] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp)
%
% Get fluorescence traces of ROIs (and background ROIs, if available) of
% movies which were taken continuously. Assumes new TIFF header format
% where frame number is logged.
%
% interp_sub - get fluorescence interpolated across ROIs as background
%
% local_comp - if being done locally, copy file to local drive temporarily
% and load/process

% Get TIFF filenames
tiff_dir = dir([tiff_path filesep '*.tif']);
tiff_files = sort({tiff_dir.name});

% append directory to tiff filename
tiff_filenames = cellfun(@(x) [tiff_path filesep x], ...
    tiff_files,'UniformOutput',false);

% if on local computer, create directory for temporary files
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

% load ROIs
load(roi_filename,'-MAT')

% get size of frame from first frame of first file
imageinfo=imfinfo(tiff_filenames{1},'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% get total number of frames from the last frame of the last file
imageinfo = imfinfo(tiff_filenames{end},'tiff');
[img_parameter img_value] = strread( ...
    imageinfo(end).ImageDescription,'%s %s', 'delimiter','=\n');
frame_tag_idx = cellfun(@(x) ~isempty(strfind(x,'Frame Tag')),img_parameter);
total_frames = str2num(img_value{frame_tag_idx});

% create masks from ROIs
cellMask = false(N*M,length(polygon.ROI));
for n_polygon = 1:length(polygon.ROI);
    if ~isfield(polygon,'autosort');
        temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
        cellMask(:,n_polygon) = temp_mask(:);
    elseif isfield(polygon,'autosort');
        temp_mask = polygon.autosort(:,:,n_polygon);
        cellMask(:,n_polygon) = temp_mask(:);
    end
end
% Normalize mask values to prepare for matrix multiplation
cellMask_norm = bsxfun(@times,cellMask,1./sum(cellMask,1));

% Grab mask for background ROIs, if available
if isfield(polygon,'bgROI')
    for n_polygon = 1:length(polygon.bgROI);
        cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
    end
    % Normalize mask values to prepare for matrix multiplation
    cellMask_bg_norm = bsxfun(@times,cellMask_bg,1./sum(cellMask_bg,1));
end

% get the traces for the ROIs
disp('Getting and concatenating traces');
roi_trace = nan(length(polygon.ROI),total_frames);
roi_trace_bg = nan(length(polygon.ROI),total_frames);

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
        
    % Load and get activity frame-by-frame to allow for any computer to do
    % many of these simultanously (as opposed to loading in whole movie,
    % which is probably faster)
    tic
    disp(['Loading file, getting activity (' img_filename ')']);
    if interp_sub
       allROI_mask = reshape(any(cellMask,2),M,N);
    end
    for loadframe = 1:numframes
        % Load current frame
        curr_frame = [];
        curr_frame = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
        
        % Get activity from current frame
        roi_trace(:,frame_tag(loadframe)) = transpose(double(curr_frame(:)')*cellMask_norm);
        
        % Get background activity from current frame if selected
        if interp_sub
            curr_frame_interp = double(curr_frame);
            curr_frame_interp(allROI_mask) = NaN;
            curr_frame_interp = inpaint_nans(curr_frame_interp,4);
            roi_trace_bg(:,frame_tag(loadframe)) = transpose(curr_frame_interp(:)'*cellMask_norm);
        elseif isfield(polygon,'bgROI')
            roi_trace_bg(:,frame_tag(loadframe)) = transpose(double(curr_frame(:)')*cellMask_bg_norm);
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

