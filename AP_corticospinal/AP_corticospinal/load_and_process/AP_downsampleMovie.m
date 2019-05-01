function AP_downsampleMovie(sum_frames,tiff_path,local_copy,overwrite,channels)
% AP_downsampleMovie(sum_frames,tiff_path,local_copy,overwrite,channels)
% Downsample all movies in tiff_path by averaging sum_frames together
%
% local_copy - if done on local PC, choose whether to make local temp files
% overwrite - specificy whether to replace any existing downsampled movies,
% default false
% channels - if multiple channels, number of channels, default 1

% If no tiff path given, prompt
if ~exist('tiff_path','var') || isempty(tiff_path)
    tiff_path = uigetdir('Select TIFF path');
end

% Get all tiff filenames in selected path
tiff_files = dir([tiff_path filesep '*.tif']);
tiff_filenames = cellfun(@(x) [tiff_path filesep x],sort({tiff_files.name}),'uni',false);

if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% Make save directory and filename
save_dir = [tiff_path filesep 'summed'];
if ~isdir(save_dir)
    mkdir(save_dir)
    % clear out folder if overwrite selected
elseif isdir(save_dir) && overwrite
    rmdir(save_dir,'s');
    mkdir(save_dir);
    % If folder already exists and overwrite not selected, exit
elseif isdir(save_dir) && ~overwrite
    return
end
[~,basename,~] = fileparts(tiff_filenames{1});
save_filename = [save_dir filesep basename '_summed_' num2str(sum_frames) '.tif'];

% If on local computer and copy desired, create dir for temporary files
if ispc && local_copy
    local_dir = 'C:\temp_files';
    % Make this filename unique in case multiple running simultaneously
    local_basename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF')];
    local_filenames = arrayfun(@(x) [local_basename '_' num2str(x) '.tif'],1:length(tiff_filenames),'uni',false);
    if ~isdir(local_dir)
        mkdir(local_dir)
    end
    im_filenames = local_filenames;
else
    im_filenames = tiff_filenames;
end

% If channels not entered, assume 1 channel
if ~exist('channels','var') || isempty(channels)
    channels = 1;
end

% Get imfinfo & total number of frames in path and group for downsampling
all_info = cell(length(tiff_filenames),1);
all_numframes = nan(length(tiff_filenames),1);
disp('Getting header information of all movies...')
fprintf('%2d',0);
for curr_tiff = 1:length(tiff_filenames)
   curr_info = imfinfo(tiff_filenames{curr_tiff},'tiff');
   all_info{curr_tiff} = curr_info;
   all_numframes(curr_tiff) = length(curr_info);
   fprintf('%c%c%2d',8,8,round(100*curr_tiff/length(tiff_filenames)));
end
disp('Done')
M = all_info{1}.Height;
N = all_info{1}.Width;

% Correct number of frames by channels
channel_numframes = all_numframes./channels;

% don't include leftover frames from division
summed_frames = floor(sum(channel_numframes)/sum_frames);
leftover_frames = mod(sum(channel_numframes),sum_frames);
channel_numframes(end) = channel_numframes(end) - leftover_frames;
if channel_numframes(end) == 0
    channel_numframes(end) = [];
end

for curr_channel = 1:channels
    frame_grp{curr_channel} = ...
        mat2cell(cell2mat(arrayfun(@(x) curr_channel:channels:channel_numframes(x), ...
        1:length(channel_numframes),'uni',false)),1, ...
        repmat(sum_frames,1,summed_frames));
end
file_grp = mat2cell(cell2mat(arrayfun(@(x) repmat(x,1,channel_numframes(x)), ...
    1:length(channel_numframes),'uni',false)),1, ...
    repmat(sum_frames,1,summed_frames));

% Load and average frames that compromise each downsampled frame
disp('Creating summed image...')
fprintf('%2d',0);
for curr_summed_frame = 1:summed_frames
       
    % If local copy, check if temp files are available, otherwise copy
    if ispc && local_copy
        used_files = unique(file_grp{curr_summed_frame});        
        % If files are available and not needed, delete them
        unused_files = arrayfun(@(x) exist(local_filenames{x},'file') & ~ismember(x,used_files),1:length(tiff_filenames));
        if any(unused_files)
            for i = find(unused_files)
                delete(local_filenames{i});
            end
        end
        % If files are needed and not available, copy them
        for i = used_files
            % If it's needed and not available, copy it
            if ~exist(local_filenames{i},'file')
                copyfile(tiff_filenames{i},local_filenames{i});
            end
        end
    end
    
    % Loop through channels, get and write average frame
    for curr_channel = 1:channels
        curr_summed_im = zeros(N,M,'uint16');
        for curr_frame = 1:sum_frames
            curr_file = file_grp{curr_summed_frame}(curr_frame);
            curr_file_frame = frame_grp{curr_channel}{curr_summed_frame}(curr_frame);
            curr_im = imread(im_filenames{curr_file},'tiff', ...
                curr_file_frame,'Info',all_info{curr_file});
            curr_summed_im = curr_summed_im + uint16(curr_im./sum_frames);
        end
        
        % write summed image to tiff
        if curr_summed_frame == 1 && channels == 1;
            imwrite(curr_summed_im,save_filename,'tif','Compression','none','WriteMode','overwrite');
        else
            imwrite(curr_summed_im,save_filename,'tif','Compression','none','WriteMode','append');
        end
    end
             
    fprintf('%c%c%2d',8,8,round(100*curr_summed_frame/summed_frames));
end

% If local copy, look for any leftover temp files and delete
if ispc && local_copy
    leftover_files = arrayfun(@(x) exist(local_filenames{x},'file'),1:length(tiff_filenames));
    if any(leftover_files)
        for i = find(leftover_files)
            delete(local_filenames{i});
        end
    end
end