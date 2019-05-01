function AP_check_motion_correction(animal,local_copy,manual_threshold)
% AP_check_motion_correction(animal,local_copy,manual_threshold)
% Check correlation of all motion corrected images for outliers
%
% local_copy - if on PC, this specifies if files should be copied to local
% drive first (default is TRUE)
% manual_threshold - manually choose the correlation cutoff for fixing
% frames as each day is evaluated (default is FALSE)


%% Set paths 

% Check if local computer
local_comp = ispc;

if nargin == 1 && local_comp
    local_copy = true;
elseif nargin == 1 && ~local_comp
    local_copy = false;
end

if nargin < 3 || ~islogical(manual_threshold) || isempty(manual_threshold)
   manual_threshold = false; 
end

if ~local_comp
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    save_filename = ['/usr/local/lab/People/Andy/Data/motion_correction_check' filesep animal '_motion_correction_check'];
else
    data_path = ['Z:\People\Andy\Data\' animal];
    save_filename = ['Z:\People\Andy\Data\motion_correction_check' filesep animal '_motion_correction_check'];
end

data_path_dir = dir(data_path);
data_path_dir_days = cellfun(@(x) any(x),regexp({data_path_dir.name},'^\d\d\d\d\d\d$'));

day_paths = cellfun(@(x) [data_path filesep x],{data_path_dir(([data_path_dir.isdir]) & data_path_dir_days).name},'uni',false);

%% Check each frame's motion correction through correlation

% Initialize summary variable
im_corr = struct('prefix',cell(length(day_paths),1),'postfix',cell(length(day_paths),1));

remake_summed_movie = false(1,length(day_paths));

for curr_day = 1:length(day_paths)
    
    % Load TIFF target (from Aki's motion correct)
    target_dir = dir([day_paths{curr_day} filesep 'target' filesep '*.tif']);
    target_file = [day_paths{curr_day} filesep 'target' filesep target_dir.name];
    target_im = imread(target_file,'tiff');
    
    % Get TIFF filenames
    tiff_dir = dir([day_paths{curr_day} filesep '*.tif']);
    tiff_files = cellfun(@(x) [day_paths{curr_day} filesep x],sort({tiff_dir.name}),'uni',false);
    
    % If on local computer and copy desired, create dir for temporary files
    if local_comp && local_copy
        local_dir = 'C:\temp_files';
        % Make this filename unique in case multiple running simultaneously
        local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
        if ~isdir(local_dir)
            mkdir(local_dir)
        end
        % Delete local temp file if it exists
        if exist(local_filename)
            delete(local_filename)
        end
    end
    
    % Loop through TIFF files, get correlation to target image
    for curr_tiff = 1:length(tiff_files)
        
        % Set appropriate image filename
        if ~local_comp || (local_comp && ~local_copy)
            img_filename = [tiff_files{curr_tiff}];
        elseif  local_comp && local_copy
            disp('Copying movie to local drive...')
            tic
            copyfile(tiff_files{curr_tiff},local_filename);
            toc
            disp('Done.')
            img_filename = local_filename;
        end
        
        % Load current file
        imageinfo=imfinfo(img_filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Height;
        N=imageinfo(1).Width;
        
        % Get 2d correlation between current and target frame
        curr_im_corr = nan(numframes,1);
        disp(['Loading & correlating image (' tiff_files{curr_tiff} ')']);
        tic
        for loadframe = 1:numframes
            im_temp = [];
            im_temp = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
            curr_im_corr(loadframe) = corr2(im_temp,target_im);
        end
        toc
        disp('Done.');
        
        % Store prefix image correlations
        im_corr(curr_day).prefix{curr_tiff} = curr_im_corr;
    end
    
    % Create array for postfix correlation (NaN, filled only if fixed)
    im_corr(curr_day).postfix = cellfun(@(x) nan(size(x)),im_corr(curr_day).prefix,'uni',false);    
 
end
       
%% Set bad correlation cutoff for all days

bad_frames = cell(size(im_corr));

for curr_day = 1:length(day_paths)
    
    im_corr_cat = vertcat(im_corr(curr_day).prefix{:});   
    im_corr_cat_diff = abs(diff(im_corr_cat));  
    
    % Set threshold for bad frames
    if ~manual_threshold
        % If automatic, define threshold by standard deviation
        im_corr_cutoff =  median(im_corr_cat)-(mean(im_corr_cat_diff)+3.5*std(im_corr_cat_diff));
    elseif manual_threshold
        % If manual, display correlation and pick threshold
        manual_thresh_fig = figure; hold on;
        title(['Day ' num2str(curr_day)]);
        corr_plot = plot(im_corr_cat,'.k');
        bad_corr_plot = plot(nan(size(im_corr_cat)),'.r');
        thresh_line = line(xlim,[min(ylim) min(ylim)],'color','b');
        good_thresh = false;
        while ~good_thresh
            [~,user_thresh] = ginput(1);
            set(thresh_line,'YData',[user_thresh,user_thresh]);
            bad_corr_y = nan(size(im_corr_cat));
            bad_corr_y(im_corr_cat < user_thresh) = im_corr_cat(im_corr_cat < user_thresh);
            set(bad_corr_plot,'YData',bad_corr_y);
            thresh_respond = input('Good threshold? (y/n)','s');
            if strcmp(thresh_respond,'y')
                good_thresh = true;
            end
            im_corr_cutoff = user_thresh;
        end
        close(manual_thresh_fig);
        
    end
    
    curr_bad_frames = cellfun(@(x) find(x < im_corr_cutoff),im_corr(curr_day).prefix,'uni',false);
    
    % Manually stopping acquisition gives few bad frames at end, ignore
    end_frames = 10;
    curr_bad_frames{end}(curr_bad_frames{end} >= length(im_corr(curr_day).prefix{end})-end_frames) = [];
    
    bad_frames{curr_day} = curr_bad_frames;
    
end

%% Correct mis-aligned frames through 2d cross correlation

for curr_day = 1:length(day_paths)
    
    bad_files = find(cellfun(@(x) ~isempty(x),bad_frames{curr_day}));
    
    % Loop through files with bad frames and re-motion correct (from raw)  
    if ~isempty(bad_files)
        disp('Found bad frames, re-aligning...');
        
        % Flag day to remake summed movie
        remake_summed_movie(curr_day) = true;
        
        for curr_tiff = bad_files
       
            % Get raw movie filename
            [~,curr_name] = fileparts(tiff_files{curr_tiff});
            uncorrected_file = curr_name(1:strfind(curr_name,'_corrected')-1);
            if ~local_comp
                uncorrected_filename = ['/usr/local/lab/Data/ImagingRig3/' ...
                    uncorrected_file(1:6) filesep animal filesep uncorrected_file '.tif'];
            else
                uncorrected_filename = ['Z:\Data\ImagingRig3\' ...
                    uncorrected_file(1:6) filesep animal filesep uncorrected_file '.tif'];
            end
            
            % If on local computer, copy image to temporary file
            if ~local_comp || (local_comp && ~local_copy)
                img_filename = [tiff_files{curr_tiff}];
            elseif local_comp && local_copy
                disp('Copying bad movie to local drive...')
                tic
                copyfile(tiff_files{curr_tiff},local_filename);
                toc
                disp('Done.')
                img_filename = local_filename;
            end
            
            % Get TIFF information
            imageinfo=imfinfo(img_filename,'tiff');
            numframes=length(imageinfo);
            M=imageinfo(1).Height;
            N=imageinfo(1).Width;
            
            % Open TIFF for overwriting and get parsing parameters
            t = Tiff(img_filename,'r+');
           
            % Load in bad frames, re-align, replace image data
            disp('Re-correcting and replacing bad frames...')
            for curr_bad_frame = bad_frames{curr_day}{curr_tiff}'
                
                % Load raw frame
                im_temp = imread(uncorrected_filename,'tiff',curr_bad_frame);
                
                % Align frame via 2d cross correlation
                cc = normxcorr2(im_temp,target_im);
                
                % Get the point of max correlation
                xcorr_midpoint = round(length(cc)/2);
                [max_c,imax] = max(cc(:));
                [ypeak xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [(ypeak-size(im_temp,2)) ...
                    (xpeak-size(im_temp,1))];
                
                % Shift image, set off-screen pixels to black
                im_aligned = circshift(im_temp,corr_offset);
                offscreen_idx = false(512,512);
                offscreen_idx(1:abs(corr_offset(1)),:) = true;
                offscreen_idx(:,1:abs(corr_offset(2))) = true;
                offscreen_idx = circshift(offscreen_idx,corr_offset.*(corr_offset<0));
                im_aligned(offscreen_idx) = 0;
                
                % Get new correlation. If it got worse, then ignore (this
                % happens occasionally and looks like it's because of a
                % stretching of the field - I think this is an imaging
                % artifact i.e. the mirror goes to fast? because it's only
                % in the vertical direction and only for 1 frame at a time)
                new_corr = corr2(im_aligned,target_im);
                if new_corr < im_corr(curr_day).prefix{curr_tiff}(curr_bad_frame)
                    continue
                end
                
                % Replace bad TIFF image with re-aligned image                
                % set current frame
                t.setDirectory(curr_bad_frame);
                % overwrite with re-aligned image
                write(t,im_aligned);
     
                % Store postfix image correlations
                im_corr(curr_day).postfix{curr_tiff}(curr_bad_frame) = new_corr;
            end
       
            % Close the edited TIFF object
            close(t);
            
            % If on local computer, copy re-corrected file back to server
            if local_comp && local_copy
                disp('Copying re-corrected movie to server...')
                tic
                copyfile(img_filename,tiff_files{curr_tiff});
                toc
                disp('Done.')
            end            
        end
        disp('Done.')
    end
    
    % Get rid of temporary file if created
    if local_comp && local_copy && exist(local_filename,'file');
        delete(local_filename)
    end
    
end


    
%% Save summary of pre and post image correlations
disp('Saving image correlation summary')
save(save_filename,'im_corr');

disp(['Finished checking motion correction for ' animal]);


%% For any days that were re-corrected: re-make summed movie
AP_fixed_motion_correction_summed_movies(animal,local_copy)









