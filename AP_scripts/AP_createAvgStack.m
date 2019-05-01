[stack_file stack_path] = uigetfile('*.tif','Pick stack');
im_info = imfinfo([stack_path stack_file]);

% hardcode number of channels
channels = 2;

% get number of frames to average from header
frame_string_begin = findstr(im_info(1).ImageDescription,'state.acq.numberOfFrames');
frame_string_end = findstr(im_info(1).ImageDescription,'state.acq.linesPerFrame');
n_frames = str2num(im_info(1).ImageDescription(frame_string_begin+25:frame_string_end-2));
n_slices = length(im_info)/(channels*n_frames);

% filename to save is old filename with AVG at the end
[save_path save_file save_ext] = fileparts([stack_path stack_file]);
save_file = [save_path filesep save_file '_AVG' save_ext];

% get image size from first file
M=im_info(1).Width;
N=im_info(1).Height;

% loop through each slice and average frames, write to final stack
w = waitbar(0,'Creating average z-stack')
for slice = 1:n_slices
    
    curr_slice_avg = zeros(N,M);
    curr_slice_raw = zeros(N,M,n_frames);
    
    tif_frames = 1 + (n_frames*channels*(slice-1)):channels: ...
        (n_frames*channels*(slice-1)) + n_frames*channels;
    for frame = 1:n_frames
        curr_slice_raw(:,:,frame)=imread([stack_path stack_file],'tiff',tif_frames(frame));
    end
    
    curr_slice_avg = mean(curr_slice_raw,3);
    curr_slice_avg = uint16(curr_slice_avg);
    
    for windows_lock = 1:100
        if slice == 1;
            try
                imwrite(curr_slice_avg,save_file,'tif','Compression','none','WriteMode','overwrite');
                break;
            catch me
                pause(0.2);
            end
        else
            try
                imwrite(curr_slice_avg,save_file,'tif','Compression','none','WriteMode','append');
                break;
            catch me
                pause(0.2);
            end
        end
    end
    waitbar(slice/n_slices,w,'Creating average z-stack');
end
close(w)