function AP_maxMovie(sum_frames,animal,all_days,person)
% AP_maxMovie(sum_frames,animal,day,person)
%
% load in multiple image files, MAX groups of sum_frames frames, save
%
% If only the number of sum_frames is given, user prompted for directory
% and all tiff files in directory are summed
%
% all_days must be cell array, or 'full' for all days in folder
%
% person can be: 'Andy', 'Jeff' (but not required, will prompt for folder)

if nargin == 4;
    if strcmp(all_days,'full')
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
        all_days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);
    end
    
    if ~iscell(all_days)
        all_days = num2cell(all_days);
    end
    
else 
    all_days = {1};
end

for curr_day = 1:length(all_days)
    clearvars -except sum_frames animal all_days person curr_day
    
    day = all_days{curr_day};
    
    % make the day a string, if it's just a number
    if ~ischar(day);
        day = num2str(day);
    end
    
    if nargin == 4
        % if person details are given, find corrected files for mouse/day
        switch person
            case 'Andy'
                tiff_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep day];
            case 'Jeff'
                tiff_path = ['/usr/local/lab/People/Jeffrey/ImageData' filesep animal filesep day];
        end
        
        dir_currfolder = dir(tiff_path);
        dir_filenames = {dir_currfolder.name};
        tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
        tiff_filenames = dir_filenames(tiff_files);
        
        % find only tiff files have been turboreg'd
        tiff_corrected_filenames_indx = find(cellfun(@(x) strcmp(x(end-11:end),'Turboreg.tif'),tiff_filenames));
        tiff_corrected_filenames = tiff_filenames(tiff_corrected_filenames_indx);
        tiff_corrected_filenames = sort(tiff_corrected_filenames);
        % display them - to know they're in order
        disp(vertcat(tiff_corrected_filenames{:}));
        
        tiff_filename_all = tiff_corrected_filenames;
        
    else
        % if no person information is given, give prompt for directory and
        % get tiff files there
        tiff_path = uigetdir;
        dir_currfolder = dir(tiff_path);
        dir_filenames = {dir_currfolder.name};
        tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
        tiff_filenames = dir_filenames(tiff_files);
        tiff_filenames = sort(tiff_filenames);
        % display them - to know they're in order
        for disp_filename = 1:length(tiff_filenames)
            disp(vertcat(tiff_filenames{disp_filename}));
        end
        
        tiff_filename_all = tiff_filenames;
    end
    
    % set up save name
    savename = [tiff_path filesep tiff_filename_all{1}(1:end-21) '_max_' ...
        num2str(sum_frames) '.tif'];
    % delete the file if it already exists
    if exist(savename)
        delete(savename)
    end
    imageinfo = imfinfo([tiff_path filesep tiff_filename_all{1}]);
    image_description = imageinfo(1).ImageDescription;
    
    % Loop through each file seperately, sum
    im_sum_all = {};
    disp('Creating maxed image...')
    fprintf('%2d',0);
    for curr_tiff = 1:length(tiff_filename_all)
        
        tiff_filename = tiff_filename_all{curr_tiff};
        
        % get number of frames
        img_filename = [tiff_path filesep tiff_filename];
        % get traces from current tiff file
        imageinfo=imfinfo(img_filename,'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        numframes = length(imageinfo);
        
        % initialize full matrix
        clear im
        im = zeros(N*M,numframes,'uint16');
        
        clear im_temp
        
        for loadframe = 1:numframes
            im_temp = imread(img_filename,'tiff',loadframe);
            im(:,loadframe) = im_temp(:);
        end
        
        % sum every sum_frames frames together
        % create group matrix
        im_sum_frames = floor(size(im,2)/sum_frames);
        frame_groups = [];
        frame_groups = 1:im_sum_frames*sum_frames;
        frame_groups = reshape(frame_groups,sum_frames,im_sum_frames);
        frame_groups = mat2cell(frame_groups,sum_frames,ones(im_sum_frames,1));
        % sum up every sum_frames frames
        im_sum = [];
        im_sum = cellfun(@(x) uint16(max(im(:,x),[],2)/sum_frames),frame_groups,'UniformOutput',false);
        im_sum = cell2mat(im_sum);
        
        % tag on frames that weren't an even multiple
        sum_cut = mod(size(im,2),sum_frames);
        if sum_cut ~= 0
            im_sum_append = sum(im(:,im_sum_frames*sum_frames+1:end),2);
            im_sum_append = im_sum_append/size(im_sum_append,2);
            im_sum = [im_sum im_sum_append];
        end
        
        im_sum = reshape(im_sum,N,M,size(im_sum,2));
        
        
        %         disp('Saving images...');
        %     fprintf('%2d',0);
        for saveframe = 1:im_sum_frames
            if saveframe == 1 && curr_tiff == 1;
                for windows_pause = 1:10
                    try
                        imwrite(im_sum(:,:,saveframe),savename,'tif','Compression','none','WriteMode','overwrite', ...
                            'Description', imageinfo(1).ImageDescription);
                        %                     fprintf('%c%c%2d',8,8,round(round(100*u/im_sum_frames)));
                        break;
                    catch me
                        pause(0.2);
                    end
                    disp('Windows Locked/Error! Didn''t Write')
                    
                end
            else
                for windows_pause = 1:10
                    try
                        imwrite(im_sum(:,:,saveframe),savename,'tif','Compression','none','WriteMode','append', ...
                            'Description', imageinfo(1).ImageDescription);
                        %                     fprintf('%c%c%2d',8,8,round(round(100*u/im_sum_frames)));
                        break;
                    catch me
                        pause(0.2);
                    end
                    disp('Windows Locked/Error! Didn''t Write')
                end
            end
        end   
        fprintf('%c%c%2d',8,8,round(100*curr_tiff/length(tiff_filename_all)));
    end
end