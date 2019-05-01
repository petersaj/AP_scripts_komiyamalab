function AP_decimateMovie(smooth_factor,animal,all_days)
% AP_decimateMovie(smooth_factor,animal,day)
% load in multiple image files, decimate (smooth, downsample), save
% all_days must be cell array

% THIS VERSION IS TO REDUCE MEMORY USAGE - not finished/started yet

disp('THIS IS CURRENTLY SPECIFIC TO ANDY - SMOOTHING 50 FRAMES NO STITCHING')

for curr_day = 1:length(all_days)
    
    day= all_days{curr_day};
    
    % make the day a string, if it's just a number
    if ~ischar(day);
        day = num2str(day);
    end
    
    tiff_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep day];
    dir_currfolder = dir(tiff_path);
    dir_filenames = {dir_currfolder.name};
    tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    tiff_filenames = dir_filenames(tiff_files);
    
    % find only tiff files that match filename format for turboreg'd only
    % this is fucked. this should be easier. either matlab can't do it, or the
    % documentation is total shit.
    tiff_corrected_regexp = regexp(tiff_filenames, ...
        '\d\d\d\d\d\d_\w\w\w\w_\d\d\d_\d\d\d_Turboreg.tif');
    tiff_corrected_filenames_indx = cellfun(@(x) ~isempty(x),tiff_corrected_regexp);
    tiff_corrected_filenames = tiff_filenames(tiff_corrected_filenames_indx);
    tiff_corrected_filenames = sort(tiff_corrected_filenames);
    % display them - to know they're in order
    disp(vertcat(tiff_corrected_filenames{:}));
    
    im = [];
    tiff_filename = tiff_corrected_filenames;
    
    % get number of frames
    disp('Calculating total number of frames');
    total_numframes = 0;
    for i = 1:length(tiff_filename);
        img_filename = [tiff_path filesep tiff_filename{i}];
        % get traces from current tiff file
        imageinfo=imfinfo(img_filename,'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        numframes=length(imageinfo);
        total_numframes = total_numframes+numframes;
    end
        
    % initialize full matrix
    disp('Initializing full matrix')
    im = zeros(N*M,total_numframes);
    
    disp('Loading file...')
    fprintf('%2d',0);
    curr_frame = 1;
    for i = 1:length(tiff_filename);
        img_filename = [tiff_path filesep tiff_filename{i}];
        % get traces from current tiff file
        imageinfo=imfinfo(img_filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        clear im_temp im_temp_long
   
        for loadframe = 1:numframes
            im_temp = imread(img_filename,'tiff',loadframe);
            im(:,curr_frame) = im_temp(:);
            curr_frame = curr_frame+1;            
        end
        
        fprintf('%c%c%2d',8,8,round(100*i/length(tiff_filename)));
    end
    
    % go through each column, smooth then downsample
    smooth_kernel = ones(1,smooth_factor)/smooth_factor;
    numframes = size(im,2);
    im = single(im);
    disp('Smoothing movie...')
    im_decimate = conv2(im,smooth_kernel,'same');
    disp('Assuming uint16');
    im_decimate = uint16(im_decimate(:,1:smooth_factor:numframes));
    decimate_frames = size(im_decimate,2);
    
    im_decimate = reshape(im_decimate,N,M,decimate_frames);
    
    savename = [tiff_path filesep tiff_filename{1}(1:11) '_decimated_' ...
        num2str(smooth_factor) '.tif'];
    image_description = imageinfo(1).ImageDescription;
    %     tic
    %     saveastiff(im_decimate,savename,1,1,1,0,image_description);
    %     toc
    disp('Saving images...');
    fprintf('%2d',0); 
    for u = 1:decimate_frames
        if u == 1 && i == 1;
            for windows_pause = 1:10
                try
                    imwrite(im_decimate(:,:,u),savename,'tif','Compression','none','WriteMode','overwrite', ...
                        'Description', imageinfo(1).ImageDescription);
                    fprintf('%c%c%2d',8,8,round(round(100*u/decimate_frames)));
                    break;
                catch me
                    pause(0.2);
                end
                disp('Windows Locked! Didn''t Write')
                
            end
        else
            for windows_pause = 1:10
                try
                    imwrite(im_decimate(:,:,u),savename,'tif','Compression','none','WriteMode','append', ...
                        'Description', imageinfo(1).ImageDescription);
                    fprintf('%c%c%2d',8,8,round(round(100*u/decimate_frames))); 
                    break;
                catch me
                    pause(0.2);
                end
                disp('Windows Locked! Didn''t Write')
            end
        end
    end
end