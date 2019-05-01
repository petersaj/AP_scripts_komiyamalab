function AP_make_bgroi(data_path,roi_path)
% AP_make_bgroi
%
% Creates background ROI from ROI files
% Draws an anulus around ROI, drops pixels that appear to have events not
% shared by the main ROI (to avoid subtracting non-contaminating neighbors)
%
% Assumes the following about file structure: 
% - tiff files are organized as data_path/6-digit day
% - a folder called 'summed_movie' exists in each day folder which contains
%   one file
% - the first 6 characters in ROI filenames is the 6-digit day


% Find days from ROI files
roi_dir = dir([roi_path filesep '*.roi']);
roi_files = sort({roi_dir.name});

days = cellfun(@(x) x(1:6),roi_files,'uni',false);

% Make folder for bg roi files and analysis
bgROI_path = [roi_path '_bg'];
if ~exist(bgROI_path,'dir')
    mkdir([roi_path '_bg']);
end

% If there is an ROI label file in the ROI folder, copy to new folder
roi_label_dir = dir([roi_path filesep '*.roilabel']);
for curr_label = 1:length(roi_label_dir)
    copyfile([roi_path filesep roi_label_dir(1).name], ...
        [bgROI_path filesep roi_label_dir(1).name]);
end

for curr_session = 1:length(roi_files)
    
    day = days{curr_session};
    
    disp(['Creating bgROIs: ' roi_files{curr_session}])
    
    %%%% Load in files
    
    % load summed movie (some uncertainty depending on which motion
    % correction used - aki automatically makes one and names it a little
    % differently, also his makes one extra file called "summed_align")
    data_day_path = [data_path filesep day];
    summed_movie_dir = dir([data_day_path filesep 'summed*']);
    summed_movie_path = [data_day_path filesep summed_movie_dir.name];
    summed_movie_file = dir([summed_movie_path filesep '*summed_50.tif']);
    summed_movie_filename = [summed_movie_path filesep summed_movie_file.name];
    
    imageinfo = imfinfo(summed_movie_filename,'tiff');
    numframes = length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;

    im = zeros(N*M,numframes,'double');
    for loadframe = 1:numframes
        curr_frame = imread(summed_movie_filename,'tiff',loadframe,'Info',imageinfo);
        im(:,loadframe) = curr_frame(:);
    end
    
    % load ROI file
    load([roi_path filesep roi_files{curr_session}],'-MAT');
    
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
    curr_save_dir = [bgROI_path filesep roi_files{curr_session}(1:end-4) '_bgROI.roi'];
    save(curr_save_dir,'polygon');    
    disp(['Created background ROIs for ' roi_files{curr_session}]);
end













