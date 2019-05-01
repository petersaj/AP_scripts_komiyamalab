function AP_make_bgROI_batch(animal)
% AP_make_bgROI_batch(animal)
% From the ground up: create fixed background ROI files, create
% concatenated traces, create analysis files with (corrected) df/f traces
% 
% NOTE: THIS IS OLD, NEW PROCESSING USES AP_make_bgroi

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

% make folder for bg roi files and analysis
bgROI_dir = [data_path filesep animal '_bgROI'];
if ~exist(bgROI_dir,'dir')
    mkdir([data_path filesep animal '_bgROI']);
end

%% Create bgROIs

for curr_day = 1:length(days)
    clearvars -except animal days curr_day bgROI_dir
    
    day = num2str(days{curr_day});
    
    disp(['Creating bgROIs: ' animal ' ' day])
    
    %%%% Load in files
    
    % load summed movie
    summed_movie_filename = ['/usr/local/lab/People/Andy/Data/' ...
        animal filesep day filesep day '_' animal '_summed_50.tif'];
    
    imageinfo=imfinfo(summed_movie_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;

    im = zeros(N*M,numframes,'double');
    for loadframe = 1:numframes
        curr_frame = imread(summed_movie_filename,'tiff',loadframe,'Info',imageinfo);
        im(:,loadframe) = curr_frame(:);
    end
    
    % load ROI file
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roi_filename = [day '_' animal '_summed_50_template.roi'];
    load([analysis_path filesep roi_filename],'-MAT');
    
    %%%% Create background ROIs
    
    % normally from GUI:
    aspect_ratio = 1;
    bg_roi_inner = 0;
    bg_roi_outer = 6;
    bg_roi_outer_full = bg_roi_outer + bg_roi_inner;
    n_polygon = length(polygon.ROI);
    
    
    % loop through ROIs, draw background ROIs with given inner/outer radius
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
        
        % if there are < 10 pixels, use whole bgroi
        if sum(polygon.bgROI_mask{curr_roi}(:)) < 10
            polygon.bgROI_mask{curr_roi}(bg_pixel_outliers) = true;
        end
    end
    
    %%%% Save the background roi
    curr_save_dir = [bgROI_dir filesep day '_' animal '_bgROI.roi'];
    save(curr_save_dir,'polygon');    
    disp(['Created ROIs for ' animal ': ' day]);
end


%% Get concatenated traces

disp('Finished creating ROIs, getting concatenated traces');
%%%% Create concatenated traces from bgROIs 
AP_getConcatTrace_batch(animal,bgROI_dir);

%% Make analysis files

%%%% Go through all analysis files, get bhv, df/f (corrected), 

for curr_day = 1:length(days)
    clearvars -except animal days curr_day bgROI_dir
    
    day = num2str(days{curr_day});
    
    % load raw trace
    analysis_filename = [bgROI_dir filesep day '_' animal '_bgROI_analysis.mat'];
    load([analysis_filename],'-MAT');
    
    % if for some reason im structure already created, use that trace
    if exist('im','var')
        roi_trace_long_bg = im.roi_trace_long_bg;
        roi_trace_long = im.roi_trace_long;
        roi_trace_long_split = im.roi_trace_split;
    end
    
    % get behavior
    bhv = AP_getBehavior_leverBscope(animal,day,roi_trace_long);
    
    % make df/f trace     
    [roi_trace_bg_df roi_trace_bg_baseline] = ...
        AP_baselineEstimation(roi_trace_long_bg,bhv.framerate);
    roi_trace_long_bgsubtract = roi_trace_long - ...
        (roi_trace_long_bg - roi_trace_bg_baseline);
    
    roi_trace_df = ...
        AP_baselineEstimation(roi_trace_long_bgsubtract,bhv.framerate);
    
    im.roi_trace_df = roi_trace_df;
    im.roi_trace_long = roi_trace_long;
    im.roi_trace_split = roi_trace_long_split;
    im.roi_trace_long_bg = roi_trace_long_bg;
    
    save(analysis_filename,'bhv', 'im');
    
    disp(['Finished analyzing: ' animal ' day ' day])
    
end
disp('bgROI Completed')










