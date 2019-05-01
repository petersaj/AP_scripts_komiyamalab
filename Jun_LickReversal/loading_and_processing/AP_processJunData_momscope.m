function AP_processJunData_momscope(animal)
% Adapted from AP_make_bgROI_batch(animal) (AP 2/13/14)
% Uses Jun's file/folder structure
%
% From the ground up: create fixed background ROI files, create
% traces, create analysis files with (corrected) df/f traces

% find days (assume all folders in _roi folder correspond to days as animal_day)
data_path = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_roi'];
data_dir = dir(data_path);
data_names = {data_dir.name};
days = cellfun(@(x) x(end-5:end),data_names(3:end),'UniformOutput',false);

% make folder for bg roi files and analysis
processed_dir = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_processedData'];
if ~exist(processed_dir,'dir')
    mkdir(processed_dir);
end

%% Create summed movies
% ROIs were drawn on raw, so create summed for effective BG roi creation
AP_sumMovie(10,animal,days,'Jun')


%% Create bgROIs
% make folder in processed directory to store bg ROIs
if ~exist([processed_dir filesep animal '_bgROI'],'dir')
    mkdir([processed_dir filesep animal '_bgROI']);
end

for curr_day = 1:length(days)
    clearvars -except animal days curr_day processed_dir data_path
    
    day = num2str(days{curr_day});
    day_path = [data_path filesep animal '_' day];
    
    curr_save_filename = [processed_dir filesep animal '_bgROI' filesep day '_' animal '_bgROI.roi'];
    if exist(curr_save_filename,'file')
        disp(['Already created bgROIs, skipping: ' animal ' ' day])
        continue
    end
    
    disp(['Creating bgROIs: ' animal ' ' day])
    
    %%%% Load in files
    
    % load summed movie
    summed_movie_filename = [data_path filesep animal '_' day ...
        filesep 'summed_movie' filesep animal '_' day '_summed_10.tif'];
    
    imageinfo=imfinfo(summed_movie_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    im = zeros(N*M,numframes,'double');
    for loadframe = 1:numframes
        curr_frame = imread(summed_movie_filename,'tiff',loadframe,'Info',imageinfo);
        im(:,loadframe) = curr_frame(:);
        disp(loadframe);disp(numframes);
    end
    
    % load ROI file
    roi_file = dir([day_path filesep '*.roi']);
    if length(roi_file) > 1
        error([animal ' ' day ' has multiple ROI files']);
    end
    roi_filename = [day_path filesep roi_file.name];
    load([roi_filename],'-MAT');
    
    %%%% Create background ROIs
    
    % normally from GUI:
    aspect_ratio = 4;
    bg_roi_inner = 5;
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
    save(curr_save_filename,'polygon');
    disp(['Created ROIs for ' animal ': ' day]);
end


%% Get concatenated traces, make processed files

disp('Finished creating ROIs, getting concatenated traces');

%%%% Create concatenated traces from bgROIs
AP_getConcatTrace_batch_momscope(animal,[processed_dir filesep animal '_bgROI']);

disp(['Finished processing: ' animal])











