function AP_sumMovie(sum_frames,animal,all_days,person,channels,zdepths)
% AP_sumMovie(sum_frames,animal,day,person,channels,zdepths)
%
% load in multiple image files, sum groups of sum_frames frames, save 
%
% If only the number of sum_frames is given, user prompted for directory
% and all tiff files in directory are summed
%
% all_days must be cell array, or 'full' for all days in folder (Andy only)
%
% person can be: 'Andy', 'Jeff', 'Jun' (but not required, will prompt for folder)
%
% channels: if multiple, assumes channel 1 is green
%
% zdepths: if using fast z, how many depths are used (if multiple, sums
% seperately but alternates final frames)

if nargin > 2 && ~isempty(all_days);
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

% assume 1 channel,zdepth if not specified
if ~exist('channels','var');
    channels = 1;
end

if ~exist('zdepths','var');
    zdepths = 1;
end

for curr_day = 1:length(all_days)
    clearvars -except sum_frames animal all_days person curr_day channels zdepths
    
    day = all_days{curr_day};
    
    % make the day a string, if it's just a number
    if ~ischar(day);
        day = num2str(day);
    end
    
    if nargin >= 4 && ~isempty(person)
        % if person details are given, find corrected files for mouse/day
        switch person
            case 'Andy'
                tiff_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep day];
            case 'Jeff'
                tiff_path = ['/usr/local/lab/People/Jeffrey/ImageData' filesep animal filesep day];
            case 'Jun'
                tiff_path = ['/usr/local/lab/People/Jun/Data' filesep animal ...
                    filesep animal '_roi' filesep animal '_' day];
        end
        
        dir_currfolder = dir(tiff_path);
        dir_filenames = {dir_currfolder.name};
        tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
        tiff_filenames = dir_filenames(tiff_files);
        
        % find only tiff files have been turboreg'd
        tiff_corrected_filenames_indx = find(cellfun(@(x) ~isempty(strfind(x,'Turboreg')),tiff_filenames));
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
    underscore_idx = strfind(tiff_filename_all{1},'_');
    curr_animalday = tiff_filename_all{1}(1:underscore_idx(2)-1);
    savedir = [tiff_path filesep 'summed_movie'];
    if ~exist(savedir,'dir');
        mkdir(savedir);
    end
    savename = [savedir filesep curr_animalday '_summed_' ...
        num2str(sum_frames) '.tif'];
    
    % skip creation if file already exists
    if exist(savename)
        disp(['Sum movie already created: ' animal ' ' day])
        continue
    end
    imageinfo = imfinfo([tiff_path filesep tiff_filename_all{1}]);
    
    % Loop through each file seperately, sum    
    disp('Creating summed image...')
    fprintf('%2d',0);
    for curr_tiff = 1:length(tiff_filename_all)
        im_sum_all = {};
        
        tiff_filename = tiff_filename_all{curr_tiff};
        
        % get number of frames
        img_filename = [tiff_path filesep tiff_filename];
        % get traces from current tiff file
        imageinfo=imfinfo(img_filename,'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        numframes = length(imageinfo);
        
        for curr_zdepth = 1:zdepths;
            % initialize full matrix
            clear im
            
            frame_idx = 1+(channels*(curr_zdepth-1)):channels*zdepths:numframes;
            im = zeros(N*M,length(frame_idx),'uint16');
            
            clear im_temp
            
            
            for loadframe = 1:length(frame_idx)
                curr_frame_idx = frame_idx(loadframe);
                im_temp = imread(img_filename,'tiff',curr_frame_idx);
                im(:,loadframe) = im_temp(:);
            end
            
            % sum every sum_frames frames together
            % create group matrix
            im_sum_frames = floor(size(im,2)/sum_frames);
            frame_groups = [];
            frame_groups = 1:im_sum_frames*sum_frames;
            frame_groups = reshape(frame_groups,sum_frames,im_sum_frames);
            frame_groups = mat2cell(frame_groups,sum_frames,ones(im_sum_frames,1));
            
            % save the frame groups to be put in image header
            imfinfo_frame_groups = cellfun(@num2str,frame_groups,'uni',false);            
            
            % sum up every sum_frames frames
            im_sum = [];
            im_sum = cellfun(@(x) uint16(sum(im(:,x),2)/sum_frames),frame_groups,'UniformOutput',false);
            im_sum = cell2mat(im_sum);
            
            % tag on frames that weren't an even multiple
            sum_cut = mod(size(im,2),sum_frames);
            if sum_cut ~= 0
                im_sum_append = sum(im(:,im_sum_frames*sum_frames+1:end),2);
                im_sum_append = im_sum_append/size(im_sum_append,2);
                im_sum = [im_sum im_sum_append];
                
                imfinfo_frame_groups{end+1} = ...
                    num2str([im_sum_frames*sum_frames+1:size(im,2)]');
            end
            im_sum_all{curr_zdepth,1} = im_sum;                       
        end
        im_sum_interlaced = reshape(cell2mat(im_sum_all),N*M,[]);
        im_sum_reshape = reshape(im_sum_interlaced,N,M,size(im_sum_interlaced,2));
        
        %         disp('Saving images...');
        %     fprintf('%2d',0);
        for saveframe = 1:im_sum_frames
            if saveframe == 1 && curr_tiff == 1;
                for windows_pause = 1:10
                    try
                        imwrite(im_sum_reshape(:,:,saveframe),savename,'tif','Compression','none','WriteMode','overwrite', ...
                            'Description', imfinfo_frame_groups{saveframe});
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
                        imwrite(im_sum_reshape(:,:,saveframe),savename,'tif','Compression','none','WriteMode','append', ...
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