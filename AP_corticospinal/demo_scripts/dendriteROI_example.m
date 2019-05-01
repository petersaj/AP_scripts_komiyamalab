%% This is to make a demo movie to show how dendriteROI works

save_dir = 'Z:\People\Andy\Corticospinal\methods_demo\dendriteROI';

curr_file = 'C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP128\150421_AP128_001_001_summed_50.tif';
im_info = imfinfo(curr_file);

M = im_info(1).Height;
N = im_info(1).Width;
n_frames = length(im_info);

% This is arbitrary, do anything better that comes along
fluor_std_threshold = 5;


% Load in current day, align, make average
im = zeros(M,N,n_frames,'uint16');
for i = 1:n_frames
    curr_frame = imread(curr_file,'tiff',i,'Info',im_info);
    im(:,:,i) = curr_frame;
end
curr_avg = nanmean(double(im),3);

% Make df image, set activity df threshold
im_processed = double(im) - repmat(double(curr_avg),1,1,size(im,3));
df_threshold = std(single(im_processed(:)))*fluor_std_threshold;

% Find centroids of active spots
im_perim = false(size(im));
im_centroid = false(size(im));
im_thresh = false(size(im));

active_centroids = zeros(512,512);

for i = 1:n_frames
    curr_frame = im_processed(:,:,i);
    
    curr_frame_thresh = curr_frame >= df_threshold;   
    
%    % Split thresholded objects by local maxima
%    curr_activity = zeros(M,N);
%    curr_activity(curr_frame_thresh) = curr_frame(curr_frame_thresh);
%    curr_local_max = imextendedmax(curr_activity,df_threshold/10);
%    curr_local_max_borders = watershed(bwdist(curr_local_max)) == 0;
%    curr_frame_thresh(curr_local_max_borders) = 0;
    
    im_thresh(:,:,i) = curr_frame_thresh;
    
    % Eliminiate objects which are too small
    min_px = 3;
    curr_frame_thresh_filt = bwareaopen(curr_frame_thresh,min_px);
        
    stats = regionprops(curr_frame_thresh_filt);
    if length(stats) == 0
        continue
    end
    active_centroids_sub = round(vertcat(stats.Centroid));
    active_centroids_ind = sub2ind([M,N],active_centroids_sub(:,2), ...
        active_centroids_sub(:,1));
    
    active_centroids(active_centroids_ind) = active_centroids(active_centroids_ind)+1;
    
    im_perim(:,:,i) = bwperim(curr_frame_thresh);
    im_centroid(:,:,i) = active_centroids ~= 0;
 
    disp(i);
end

% Make image for showing process
im_compiled = permute(cat(4,cat(4,mat2gray(im,[0 2000]),im_perim),im_centroid),[1 2 4 3]);

active_centroids_max = imextendedmax(active_centroids,2);
%active_centroids_max = imregionalmax(active_centroids);

% Force borders to zero
x_crop_px = 25;
y_crop_px = 15;
border_idx = true(size(active_centroids_max));
border_idx(y_crop_px:end-y_crop_px+1,x_crop_px:end-x_crop_px+1) = false;
active_centroids_max(border_idx) = false;

% Store indicies of each ROI centroid
cc = bwconncomp(active_centroids_max);
active_centroid_rois = cc.PixelIdxList;

% Find shape of ROI
% Find frames where ROI is active (centroid = 1);
frac_active_frames = 0.5;
roi_mask = false(M,N,length(active_centroid_rois));
for curr_roi = 1:length(active_centroid_rois)
    [curr_y curr_x] = ind2sub([M,N],active_centroid_rois{curr_roi});
    active_frames = any(cell2mat(arrayfun(@(x) im_thresh(curr_y(x),curr_x(x),:),1:length(curr_y),'uni',false)),2);
    active_px = nansum(im_thresh(:,:,active_frames),3);
    
    [curr_y curr_x] = ind2sub([M,N],active_centroid_rois{curr_roi});
    
    if ~any(active_px(:))
        continue;
    end
    roi_mask(:,:,curr_roi) = bwselect(active_px >= sum(active_frames)*frac_active_frames,curr_x,curr_y); 
    disp([num2str(curr_roi) '/' num2str(length(active_centroid_rois))]);
end



%% Create demo movies and figures

% write compiled movie to avi
compiled_savename = [save_dir filesep 'dendriteROI_compiled.avi'];

writerObj = VideoWriter(compiled_savename);
writerObj.FrameRate = 100;
open(writerObj);

fig = figure;
set(0,'defaultaxesposition',[0 0 1 1])
ax = imagesc(im_compiled(:,:,:,1));
axis off;
writeVideo(writerObj,getframe(fig));
for i = 2:size(im_compiled,4)
    set(ax,'CData',im_compiled(:,:,:,i))
    writeVideo(writerObj,getframe(fig));
end
close(fig);
close(writerObj);









