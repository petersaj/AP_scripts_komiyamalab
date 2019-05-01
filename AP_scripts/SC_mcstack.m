function SC_mcstack()
% Turboreg correct slices and create average stack

% Pick directory with stacks to motion correct
stack_path = uigetdir('*.tif','Pick folder of stacks');

% Find tiff files in directory (assume all tiff files are stacks)
dir_stack = dir(stack_path);
dir_stack_filenames = {dir_stack.name};
tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_stack_filenames);
tiff_filenames = dir_stack_filenames(tiff_files);

for curr_stack = 1:length(tiff_filenames)
    
    curr_stack_filename = tiff_filenames{curr_stack};
    
    warning off
    tiffobj = Tiff([stack_path filesep curr_stack_filename],'r');
    stack_info = tiffobj.getTag('ImageDescription');
    tiffobj.close();
    warning on
    
    [img_parameter img_value] = strread(stack_info,'%s %s', 'delimiter','=\n');
    
    % Find number of frames
    numframes_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.acqNumFrames')),img_parameter,'UniformOutput',0));
    numframes = str2num(img_value{numframes_indx});
    
    % Find number of slices
    numslices_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.stackNumSlices')),img_parameter,'UniformOutput',0));
    numslices = str2num(img_value{numslices_indx});
    
    % Find pixels to a side (assume square)
    numpixels_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanPixelsPerLine')),img_parameter,'UniformOutput',0));
    numpixels = str2num(img_value{numpixels_indx});
    N = numpixels;
    M = numpixels;
    
    % get and setup save directory - just take the directory name from the
    % last two folders from the stack path
    dir_sep = strfind(stack_path,filesep);
    save_dir = ['/usr/local/lab/People/Simon/Data' stack_path(dir_sep(end-2):end)];
    save_file = [save_dir filesep curr_stack_filename(1:end-4) '_turboreg_avg.tif'];
    
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
    
    w = waitbar(0,['Registering and creating average z-stack:' ...
        num2str(curr_stack) '/' num2str(length(tiff_filenames))]);
    
    for curr_slice = 1:numslices
        
        % load in current slice
        curr_slice_raw = zeros(N,M,numframes,'uint16');
        
        tif_frames = (curr_slice-1)*numframes+1:(curr_slice-1)*numframes+numframes;
        for frame = 1:numframes
            curr_slice_raw(:,:,frame)=imread( ...
                [stack_path filesep curr_stack_filename],'tiff',tif_frames(frame));
        end
        
        % find the most stable frames to create reference image
        num_stable_frames = 5;
        im_ref = reshape(curr_slice_raw,N*M,numframes);
        im_diff = sum(abs(diff(im_ref,1,2)));
        im_diff_smooth = smooth(im_diff,num_stable_frames);
        [im_diff_min im_diff_min_indx] = min(im_diff_smooth( ...
            floor(num_stable_frames/2):end-floor(num_stable_frames/2)));
        im_diff_min_indx = im_diff_min_indx + floor(num_stable_frames/2);
        stable_frames = im_diff_min_indx - floor(num_stable_frames/2) : ...
            im_diff_min_indx + floor(num_stable_frames/2);
        
        % make reference image from most stable frames
        curr_slice_ref = mean(curr_slice_raw(:,:,stable_frames),3);
        
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
        
        % Make ImageJ object from reference
        target_ij = array2ijStackAK(curr_slice_ref);
        im_registered = zeros(numpixels,numpixels,numframes,'uint16');
        for curr_frame = 1:numframes;
            curr_slice_frame = curr_slice_raw(:,:,curr_frame);
            source_ij=array2ijStackAK(curr_slice_frame);
            
            al=IJAlign_AK;
            registered_ij = al.doAlign(cmdstr, source_ij, target_ij);
            registered = ij2arrayAK(registered_ij);
            a=uint16(round(registered));
            im_registered(:,:,curr_frame) = a;
        end
        
        im_registered_mean = uint16(mean(im_registered,3));
        
        for windows_lock = 1:100
            if curr_slice == 1;
                try
                    imwrite(im_registered_mean,save_file,'tif','Compression','none','WriteMode','overwrite');
                    break;
                catch me
                    pause(0.2);
                end
            else
                try
                    imwrite(im_registered_mean,save_file,'tif','Compression','none','WriteMode','append');
                    break;
                catch me
                    pause(0.2);
                end
            end
        end
        waitbar(curr_slice/numslices,w,['Registering and creating average z-stack:' ...
        num2str(curr_stack) '/' num2str(length(tiff_filenames))]);
        
    end
    close(w)
    
end
