function AP_mcfastz(use_slices)
% AP_mcfastz(use_slices)
%
% Turboreg correct within each slice of fast-z
% Multiple stacks can be selected if part of same file, one ref used
% NOTE: if acquisition split into multiple stacks, assumes complete volumes
% within each file (number of frames within file must be multiple of number
% of slices)
%
% use_slices - specificy slices to be kept if discarding slices
% (which happens if long range and so added slices to allow for travel time)

% Pick stacks to motion correct and average
[tiff_filenames,stack_path] = uigetfile('*.tif','Choose stacks','multiselect','on');

if ~iscell(tiff_filenames)
    tiff_filenames = {tiff_filenames};
end

tiff_filenames = sort(tiff_filenames);

% Get information from the first stack

curr_stack_filename = tiff_filenames{1};

im_info = imfinfo([stack_path filesep curr_stack_filename]);
numframes = length(im_info);

[img_parameter,img_value] = strread(im_info(1).ImageDescription,'%s %s', 'delimiter','=\n');

% Confirm that stack is fast-z: if not, error out
fastz_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.fastZEnable')),img_parameter,'UniformOutput',0));
fastz = str2num(img_value{fastz_indx});
if ~fastz
    error('File is not fast-z');
end

% Find number of volumes
numframes_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.fastZNumVolumes')),img_parameter,'UniformOutput',0));
total_frames = str2num(img_value{numframes_indx});

% Find number of slices
numslices_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.stackNumSlices')),img_parameter,'UniformOutput',0));
numslices = str2num(img_value{numslices_indx});

if ~exist('use_slices','var') || isempty(use_slices);
    use_slices = 1:numslices;
end

% Find pixels to a side
M = im_info(1).Height;
N = im_info(1).Width;

% Save into folder within the uncorrected folder
dir_sep = strfind(stack_path,filesep);
save_dir = [stack_path filesep 'motion_corrected'];

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% Loop through all slices, get average slice for reference, motion correct all
% frames in a slice, then average all slices for final stack

% Add ImageJ java controls for turboreg
spath = javaclasspath('-static');
spath = cell2mat(spath(1));
javafolder = strfind(spath,['java' filesep]);
javafolder = spath(1:javafolder+3);

javaaddpath([javafolder filesep 'jar' filesep 'ij.jar'])
javaaddpath([javafolder filesep 'jar' filesep 'mij.jar'])
javaaddpath([javafolder])

w = waitbar(0,['Creating references for each slice']);

% Make reference image for all slices from the first stack
slice_ref = cell(numslices,1);
for curr_slice = use_slices
    
    numframes_slice = numframes/numslices;
    
    % load in current slice
    curr_slice_raw = repmat(imread( ...
        [stack_path filesep curr_stack_filename],'tiff',1),[1,1,numframes_slice]);
    
    slice_frames = curr_slice:numslices:numframes;
    for frame = 1:numframes_slice
        curr_slice_raw(:,:,frame)=imread( ...
            [stack_path filesep curr_stack_filename],'tiff',slice_frames(frame));
    end
    
    % find the most stable frames to create reference image
    num_stable_frames = 5;
    im_ref = reshape(curr_slice_raw,N*M,numframes_slice);
    im_diff = sum(abs(diff(im_ref,1,2)));
    im_diff_smooth = smooth(im_diff,num_stable_frames);
    [im_diff_min im_diff_min_indx] = min(im_diff_smooth( ...
        floor(num_stable_frames/2):end-floor(num_stable_frames/2)));
    im_diff_min_indx = im_diff_min_indx + floor(num_stable_frames/2);
    stable_frames = im_diff_min_indx - floor(num_stable_frames/2) : ...
        im_diff_min_indx + floor(num_stable_frames/2);
    
    % make reference image from most stable frames
    curr_slice_ref = mean(curr_slice_raw(:,:,stable_frames),3);
    slice_ref{curr_slice} = array2ijStackAK(curr_slice_ref);
    
    waitbar(curr_slice/numslices,w,['Creating references for each slice']);
    
end



for curr_stack = 1:length(tiff_filenames)
    
    waitbar(0,w,['Registering fast z-stack:' ...
        num2str(curr_stack) '/' num2str(length(tiff_filenames))]);
    
    curr_stack_filename = tiff_filenames{curr_stack};
    
    im_info = imfinfo([stack_path filesep curr_stack_filename]);
    
    % Check that total number of frames = slices*volumes
    % If not volumes are not complete, error out
    if mod(length(im_info),numslices) ~= 0
        error(['Not all volumes are complete: stack ' num2str(curr_stack)]);
    end
    
    curr_file_frames = length(im_info);
    
    % Set up save filename
    save_file = [save_dir filesep curr_stack_filename(1:end-4) '_turboreg.tif'];
    
    % Only use frames corresponding to selected slices
    use_frames = sort(cell2mat(cellfun(@(x) x:numslices: ...
        curr_file_frames,num2cell(use_slices),'uni',false)));
    
    % Correct every frame to its appropriate reference
    for curr_frame_idx = 1:length(use_frames)
        
        curr_frame = use_frames(curr_frame_idx);
        
        curr_slice = mod(curr_frame-1,numslices) + 1;
        
        % Set parameters for turboreg
        % Crop 10 frames on all sides: ignore blank flyback problems
        crop_border = 0;
        Cropping = [num2str(crop_border) ' ' num2str(crop_border) ' '...
            num2str(size(curr_slice_ref,2)-crop_border) ' ' ...
            num2str(size(curr_slice_ref,1)-crop_border)];
        Transformation='-translation';
        center_x = num2str(fix(size(curr_slice_ref,2)/2)-1-crop_border);
        center_y = num2str(fix(size(curr_slice_ref,1)/2)-1-crop_border);
        landmarks = [center_x,' ',center_y,' ',center_x,' ',center_y];
        cmdstr=['-align -window s ',Cropping,' -window t ', Cropping,' ', ...
            Transformation, ' ', landmarks,' -hideOutput'];
        
        % Load frame, turboreg
        curr_slice_raw = imread( ...
            [stack_path filesep curr_stack_filename],'tiff',curr_frame);
        
        im_registered = curr_slice_raw;
        source_ij=array2ijStackAK(curr_slice_raw);
        
        al=IJAlign_AK;
        registered_ij = al.doAlign(cmdstr, source_ij, slice_ref{curr_slice});
        registered = ij2arrayAK(registered_ij);
        a = uint16(round(registered));
        im_registered = a;
        
        for windows_lock = 1:100
            if curr_frame_idx == 1;
                try
                    imwrite(im_registered,save_file,'tif','Compression','none','WriteMode','overwrite', ...
                        'Description', im_info(curr_frame).ImageDescription);
                    break;
                catch me
                    pause(0.2);
                end
            else
                try
                    imwrite(im_registered,save_file,'tif','Compression','none','WriteMode','append', ...
                        'Description', im_info(curr_frame).ImageDescription);
                    break;
                catch me
                    pause(0.2);
                end
            end
        end
        waitbar(curr_frame_idx/length(use_frames),w,['Registering fast z-stack:' ...
            num2str(curr_stack) '/' num2str(length(tiff_filenames))]);
        
    end    
end
close(w);
