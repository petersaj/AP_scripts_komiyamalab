function im = AP_load_tiff(tiff_file,channels)

% Load tiff file
% im = AP_load_tiff(tiff_file,channels)
%
% 'channels': if multiple channels, 2-element vector of 
% [number of channels, channel to load]

if ~exist('tiff_file','var') || isempty(tiff_file)
    [tiff_filename, tiff_path] = uigetfile('*.tif');
    tiff_file = [tiff_path filesep tiff_filename];
end

if ~exist('channels','var') || isempty(channels)
    channels = [1,1];
end

im_info = imfinfo(tiff_file);
M = im_info(1).Height;
N = im_info(1).Width;
numframes = length(im_info);

numframes_channel = numframes/channels(1);

im = imread(tiff_file,'tiff',1,'Info',im_info);
im = repmat(im,[1,1,numframes_channel]);

load_frames = channels(2):channels(1):numframes;

for curr_frame_idx = 1:length(load_frames);
    curr_frame = load_frames(curr_frame_idx);
    im(:,:,curr_frame_idx) = imread(tiff_file,'tiff',curr_frame,'Info',im_info);
    disp([num2str(curr_frame) '/' num2str(numframes)]);
end