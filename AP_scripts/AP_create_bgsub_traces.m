function roi_traces = AP_create_bgsub_traces(framerate,roi_filename,summed_movie_filename,tiff_path)
% AP_create_bgsub_traces(roi,tiff_path) 
%
% Inputs//
%
% framerate: the imaging framerate.
% NOTE -- this is used for the df/f normalization step
% (AP_baselineEstimation). 
%
% roi: file containing previously drawn ROIs (if not included, prompted)
% NOTE -- the format of rois has to be the following:
% a structure called "polygon" containing a field called "ROI" which is a 
% cell array containing the verticies (x,y) in clockwise or 
% counterclockwise order for each ROI 
% (this is the format returned by matlab's roipoly function)
%
% summed_movie_filename: the filename for the summed movie (if not
% included, prompted)
% NOTE -- this is the filename for a TIFF file which is a downsampled 
% (by summing adjacent frames) version of the complete original movie. 
% I can provide code for this if it's wanted (AP_sumMovie). The reason for
% this is to have a less noisy estimate of pixel values in order to more
% efficiently exclude unwanted background ROI pixels. This step is not
% necessarily important: the downsampled movie can be created within this
% function, or maybe depending on microscope the raw movie can be used.
%
% tiff_path: the directory containing the TIFF data
% NOTE -- the only TIFF files in this directory are the data used to
% make the complete ROI traces. They have to be named in a sortable order
% (i.e. numbered)
%
% Output//
%
% roi_traces: structure with final traces
% This structure contains 3 fields:
% roi_trace_bgsub_df_split: the df/f-normalized background subtraced traces
% roi_trace_split: the non-normalized, un-subtracted ROI trace
% roi_trace_bg_split: the non-normalized background ROI trace

if nargin < 1 || isempty(framerate)
    error('No framerate provided')
end

% If no tiff directory is inputted, prompt a choice
if nargin < 4 || isempty(tiff_path)
    tiff_path = uigetdir([],'Choose TIFF directory');
end


%% Create bgROIs

disp(['Creating background ROIs...'])

%%%% Load in files

% Load summed movie (if not given, prompt a choice)
if nargin < 3 || isempty(summed_movie_filename)
    [summed_movie_file summed_movie_path] = uigetfile('*.tif','Choose summed movie file');
    summed_movie_filename = [summed_movie_path summed_movie_file];
end

% load ROI file (if not given, prompt a choice)
if nargin < 2 || isempty(roi_filename)
    [roi_file roi_path] = uigetfile('*','Choose ROI file');
    roi_filename = [roi_path roi_file];
end

imageinfo=imfinfo(summed_movie_filename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

im = zeros(N*M,numframes,'double');
for loadframe = 1:numframes
    curr_frame = imread(summed_movie_filename,'tiff',loadframe,'Info',imageinfo);
    im(:,loadframe) = curr_frame(:);
end

load(roi_filename,'-MAT');

%%%% Create background ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: SET THESE TO DESIRED VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aspect ratio of the pixels in tiff files (width / height)
aspect_ratio = 4;
% Inner radius of background ROIs to be drawn (distance from ROI to inner
% edge of background ROI
bg_roi_inner = 0;
% Outer radius of background ROIs to be drawn (distance from inner
% background ROI edge to outer ROI edge)
bg_roi_outer = 6;
% The minimum number of pixels to allow for a background ROI after
% eliminating outliers (if under minimum, no outliers are removed)
min_bgpixels = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bg_roi_outer_full = bg_roi_outer + bg_roi_inner;
n_polygon = length(polygon.ROI);

% Loop through ROIs, draw background ROIs with given inner/outer radius
for curr_roi = 1:n_polygon
    
    % get center of current polygon
    [area roi_center_x roi_center_y] = ...
        polycenter(polygon.ROI{curr_roi}(:,1), ...
        polygon.ROI{curr_roi}(:,2));
    
    % center current polygon at zero
    temp_roi = [polygon.ROI{curr_roi}(:,1) - roi_center_x ...
        polygon.ROI{curr_roi}(:,2) - roi_center_y];
    
    % expand polygon by inner and outer diameters
    temp_roi_inner = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_inner ...
        temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_inner/aspect_ratio];
    temp_roi_outer = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_outer_full ...
        temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_outer_full/aspect_ratio];
    
    % connect the inner and outer sections
    bg_roi = [temp_roi_outer;temp_roi_inner;temp_roi_outer(1,:)];
    bg_roi_centered = [bg_roi(:,1) + roi_center_x ...
        bg_roi(:,2) + roi_center_y];
    
    polygon.bgROI{curr_roi} = bg_roi_centered;
    
    curr_bgroi_mask = poly2mask(polygon.bgROI{curr_roi}(:,1)',...
        polygon.bgROI{curr_roi}(:,2)',N,M);
    
    polygon.bgROI{curr_roi} = bg_roi_centered;
    polygon.bgROI_mask{curr_roi} = curr_bgroi_mask;
end

%%%% Fix bgROIs by cutting out pixels which are abnormally brighter

for curr_roi = 1:n_polygon
    
    % get ROI trace
    curr_roi_mask = poly2mask(polygon.ROI{curr_roi}(:,1)',...
        polygon.ROI{curr_roi}(:,2)',N,M);
    curr_roi_trace = mean(im(curr_roi_mask(:),:));
    
    % get background trace
    curr_bgroi_pixel_trace = im(polygon.bgROI_mask{curr_roi}(:),:);
    
    bg_diff = ...
        bsxfun(@minus,double(curr_bgroi_pixel_trace),double(curr_roi_trace));
    
    % define outliers by the distribution of differences. There is a
    % gaussian (sometimes a little buried) - so just estimate the mean with
    % mode and the std
    
    % get mode rounded to nearest integer
    bg_diff_mode = mode(round(bg_diff(:)*1))/1;
    bg_diff_std = std(bg_diff(:));
    bg_cutoff = bg_diff_mode + 2*bg_diff_std;
    
    bg_trace_outliers = ...
        bsxfun(@gt,curr_bgroi_pixel_trace,curr_roi_trace+bg_cutoff);
    bg_mask_indx = find(polygon.bgROI_mask{curr_roi});
    bg_pixel_outliers = bg_mask_indx(any(bg_trace_outliers,2));
    
    polygon.bgROI_mask{curr_roi}(bg_pixel_outliers) = false;
    
    % if there are < min_bgpixels pixels, use whole bgroi
    if sum(polygon.bgROI_mask{curr_roi}(:)) < min_bgpixels
        polygon.bgROI_mask{curr_roi}(bg_pixel_outliers) = true;
    end
end

% Clear the summed movie from memory
clear im


%% Get concatenated traces

disp('Finished creating ROIs, getting concatenated traces');

tiff_dir = dir([tiff_path filesep '*.tif']);
tiff_filename = sort({tiff_dir(:).name});

disp('Getting and concatenating traces');
if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end

roi_trace_long_split = cell(1,length(tiff_filename));
roi_trace_long_split_bg = cell(1,length(tiff_filename));
for i = 1:length(tiff_filename);
    img_filename = [tiff_path filesep tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % Create mask for ROIs
    for n_polygon = 1:length(polygon.ROI);
        temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
        cellMask(:,n_polygon) = temp_mask(:);
    end
    
    % Grab masks for background ROIs
    for n_polygon = 1:length(polygon.bgROI);
        cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
    end
    
    im = [];
    im = zeros(N*M,numframes);
    
    tic
    disp('Loading file....')
    for loadframe = 1:numframes
        im_temp = [];
        im_temp = imread(img_filename,'tiff',loadframe);
        % image: columns are frames
        im(:,loadframe) = im_temp(:);
    end
    disp('Done')
    toc
    
    disp(['Getting activity from file ' num2str(i) '/' num2str(length(tiff_filename))]);
 
    roi_trace = zeros(size(cellMask,2),numframes);
    traceVal = zeros(size(cellMask,2),numframes);
    
    roi_trace_bg = zeros(size(cellMask,2),numframes);
    traceVal_bg = zeros(size(cellMask,2),numframes);
    
    for n_polygon = 1:size(cellMask,2)
        % sum of pixel brighness / number of pixels, avg fluor in ROI
        traceVal(n_polygon,:) = nanmean(im(cellMask(:,n_polygon),:));
        traceVal_bg(n_polygon,:) = nanmean(im(cellMask_bg(:,n_polygon),:));
    end
    
    roi_trace = traceVal;
    roi_trace_bg = traceVal_bg;
    
    % store traces
    roi_trace_long_split{i} = roi_trace;    
    roi_trace_long_bg_split{i} = roi_trace_bg;
    
    clear roi_trace roi_trace_bg im
end

%% Make final df/f background-subtracted traces

roi_trace_long_cat = horzcat(roi_trace_long_split{:});
roi_trace_long_bg_cat = horzcat(roi_trace_long_bg_split{:});

% make df/f trace
[roi_trace_bg_df roi_trace_bg_baseline] = ...
    AP_baselineEstimation(roi_trace_long_bg_cat,framerate);
roi_trace_long_bgsubtract = roi_trace_long_cat - ...
    (roi_trace_long_bg_cat - roi_trace_bg_baseline);

roi_trace_bgsub_df = ...
    AP_baselineEstimation(roi_trace_long_bgsubtract,framerate);

frames_per_file = cellfun(@(x) size(x,2),roi_trace_long_split);
roi_trace_bgsub_df_split = mat2cell(roi_trace_bgsub_df, ...
    size(roi_trace_long_cat,1),frames_per_file);

roi_traces.roi_trace_bgsub_df_split = roi_trace_bgsub_df_split;
roi_traces.roi_trace_split = roi_trace_long_split;
roi_traces.roi_trace_bg_split = roi_trace_long_bg_split;

disp('Finished creating background-subtracted traces')






