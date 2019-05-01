%% Motion correct using 2D cross correlation of stable channel
% (eventually should subtract off GCaMP channel in order to get more stable
% other channel)



% NOTE: THIS IS NO FASTER THAN TURBOREG (?? ask Aki for his version, which
% is apparently very fast)


tic;

orig_path = '/usr/local/lab/Data/ImagingRig3/140922/AP111';
save_path = '/usr/local/lab/People/Andy/Data/AP111/test_xcorr_mc';
channels = 2;
border = 50;
smooth_px = 5;

reference_channel = 2;
correct_channel = 1;

% Make save directory if it doesn't exist
if ~exist(save_path,'dir')
    mkdir(save_path);
end

tiff_dir = dir([orig_path filesep '*.tif']);
tiff_files = sort({tiff_dir.name});
source_filenames = cellfun(@(x) [orig_path filesep x],tiff_files,'uni',false);
destination_filenames = cellfun(@(x) [save_path filesep x],tiff_files,'uni',false);

% Create template from average of first file
imageinfo = imfinfo(source_filenames{1},'tiff');
numframes = length(imageinfo);
M = imageinfo(1).Height;
N = imageinfo(1).Width;

reference_im = zeros(N,M,1);
for curr_frame = reference_channel:channels:numframes
    reference_im = reference_im + double(imread(source_filenames{1},'tiff',...
        curr_frame,'Info',imageinfo))./(numframes/channels);
end

im_conv = fspecial('disk',smooth_px);
reference_im_smooth = conv2(reference_im,im_conv,'same');
disp('Motion correcting...')
fprintf('%2d%%',0);
% Get maximum cross correlation between every frame and reference image
for curr_file = 1:length(source_filenames)
    
    imageinfo = imfinfo(source_filenames{curr_file},'tiff');
    numframes = length(imageinfo);
    
    save_filename = [save_path filesep tiff_files{curr_file}(1:end-4) '_xcorrected.tif']; 
    
    for curr_frame = 1:numframes/channels
        
        registered_im = zeros(M,N,'uint16');
        
        % Get corresponding frame numbers
        reference_frame = (curr_frame-1)*channels + reference_channel;
        correct_frame = (curr_frame-1)*channels + correct_channel;
        
        % Load, smooth reference channel frame
        curr_ref_im = imread(source_filenames{curr_file},'tiff',...
            reference_frame,'Info',imageinfo);
        
        curr_ref_im_smooth = conv2(double(curr_ref_im),im_conv,'same');
        
        % Get cross correlation between current frame and reference image
        curr_xc = normxcorr2(reference_im(border:end-border,border:end-border), ...
            curr_ref_im_smooth(border:end-border,border:end-border));
        
        % Get the image offset
        xcorr_midpoint = round(length(curr_xc)/2);
        [max_c,imax] = max(curr_xc(:));
        [ypeak xpeak] = ind2sub(size(curr_xc),imax(1));
        corr_offset = [xcorr_midpoint-ypeak xcorr_midpoint-xpeak];
        
        % Load frame to correct
        curr_correct_im = imread(source_filenames{curr_file},'tiff',...
            correct_frame,'Info',imageinfo);

        registered_im(max(1,1+corr_offset(1)):min(N,N+corr_offset(1)), ...
            max(1,1+corr_offset(2)):min(M,M+corr_offset(2))) = ...
            curr_correct_im(max(1,-corr_offset(1)+1):min(N,N-corr_offset(1)), ...
            max(1,-corr_offset(2)+1):min(M,M-corr_offset(2)));
        
        % Write corrected frame
        if curr_frame == 1
            imwrite(registered_im, ...
                save_filename,'tif','Compression','none','WriteMode','overwrite', ...
                'Description', imageinfo(correct_frame).ImageDescription);
        else
            imwrite(registered_im, ...
                save_filename,'tif','Compression','none','WriteMode','append', ...
                'Description', imageinfo(correct_frame).ImageDescription);
        end
        
    end   
    
    fprintf('%c%c%c%2d%%',8,8,8,round(100*curr_file/length(source_filenames)));
    
end

toc;






