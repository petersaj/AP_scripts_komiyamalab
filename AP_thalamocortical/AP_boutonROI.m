%% Load in video 

[im_file im_path] = uigetfile('*.tif');
im_fullfile = [im_path im_file];
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info));
for i = 1:length(im_info)
    im(:,:,i) = imread(im_fullfile,'tiff',i,'Info',im_info);
    disp(i) 
end


%% try auto-rois for boutons/axons

[im_filename im_path]= uigetfile('*.tif');
im_fullfilename = [im_path im_filename];
imageinfo = imfinfo(im_fullfilename);

border = 20;
border_idx = border+1:512-border;

% TRY SMOOTHING THE VIDEO HERE
im = zeros(length(border_idx)^2,length(imageinfo));

im_lumnorm = zeros(length(border_idx),length(border_idx),length(imageinfo));
smooth_disk = fspecial('gaussian',4,1);
h = fspecial('average', [20 20]);

for i = 1:length(imageinfo)
    curr_im = double(imread(im_fullfilename,'tiff',i,'Info',imageinfo));
    curr_im_smooth = imfilter(curr_im,smooth_disk);
    im(:,i) = reshape(curr_im_smooth(border_idx,border_idx),[],1);
    im_lumnorm(:,:,i) = curr_im_smooth(border_idx, border_idx)./ ...
        imfilter(curr_im_smooth(border_idx, border_idx),h);
    disp(i);
end
im_max = max(im,[],2);
im_lumnorm_diff = diff(im_lumnorm.^3,[],3);
im_lumnorm_reshape = reshape(im_lumnorm,[],size(im_lumnorm,3));

max_im_lumnorm = max(im_lumnorm_diff,[],3);

thresh_cutoff = 2*std(max_im_lumnorm(:));
max_lumnorm_thresh = max_im_lumnorm > thresh_cutoff;

% get rid of clusters of less than n pixels
smallest_roi = 10;
cc = bwconncomp(max_lumnorm_thresh);
rois = cc.PixelIdxList;
rois(cellfun(@length,rois) < smallest_roi) = [];

roi_trace = cell2mat(cellfun(@(x) nanmean(im_lumnorm_reshape(x,:),1),rois','uni',false));

% Have some sort of check (correlation?) that all pixels within ROI
% actually do the same thing (if not, split or delete)

% Try to find entire dendrite by looking for correlated pixels?

% ROI-pixel correlation

roi_norm = bsxfun(@minus,roi_trace,mean(roi_trace,2));
pixel_norm = bsxfun(@minus,im_lumnorm_reshape,mean(im_lumnorm_reshape,2));

roi_norm = bsxfun(@times,roi_norm,1./sqrt(sum(roi_norm.^2,2)));
pixel_norm = bsxfun(@times,pixel_norm,1./sqrt(sum(pixel_norm.^2,2)));

roi_pixel_corr = roi_norm*pixel_norm';

% try ICA on the pixel correlation maps
[icasig, A, W] = fastica(roi_pixel_corr,'maxNumIterations',100);
ics = icasig';
is_neg_skew = skewness(ics) < 0;
ics(:, is_neg_skew) = -ics(:, is_neg_skew);

% get correlation of ICs with all pixels
ics_trace = ics'*im_lumnorm_reshape;
ics_trace_norm = bsxfun(@minus,ics_trace,mean(ics_trace,2));
pixel_norm = bsxfun(@minus,im,mean(im,2));

roi_norm = bsxfun(@times,roi_norm,1./sqrt(sum(roi_norm.^2,2)));
pixel_norm = bsxfun(@times,pixel_norm,1./sqrt(sum(pixel_norm.^2,2)));

roi_pixel_corr = roi_norm*pixel_norm';


% segment by thresholding ICs and then picking the largest ROI?
% didn't really work: too noisy
% ic_threshold = nanmean(ics) + 2*nanstd(ics);
% binary_ics = reshape(bsxfun(@gt,ics,ic_threshold),sqrt(size(ics,1)),[],size(ics,2));
% 
% binaryCell = [];
% binaryCell = bwselect(binaryImg,x,y); % select cell that was clicked
% if sum(binaryCell(:)) < 0.5*cell_area
%     disp('No coherent cell detected')
%     return
% end
% first_nonzero = find(binaryCell > 0,1);
% [y_nonzero x_nonzero] = ind2sub([N M],first_nonzero);
% roi_edge = bwtraceboundary(binaryCell,[y_nonzero x_nonzero],'N'); % find boundary, in order
% % create ROI from detected cell
% polygon = get(handles.roiList,'UserData');
% n_polygon = size(get(handles.roiList,'string'),1) + 1;


%%
% try ICA on the pixel correlation maps?
[icasig, A, W] = fastica(roi_pixel_corr);
ics = icasig';
is_neg_skew = skewness(ics) < 0;
ics(:, is_neg_skew) = -ics(:, is_neg_skew);

% try segmenting by std cutoffs?
ic_px = zeros(size(ics));
for curr_ic = 1:size(ics,2);
    curr_thresh_ic = ics(:,curr_ic) > std(ics(:,curr_ic));
    cc = bwconncomp(reshape(curr_thresh_ic,sqrt(length(curr_thresh_ic)),[]));
    ic_conn = cc.PixelIdxList;
    ic_conn(cellfun(@length,ic_conn) < 10) = [];

    ic_conn_px = vertcat(ic_conn{:});
    ic_px(ic_conn_px,curr_ic) = 1;    
end

%%
curr_mov = im_lumnorm_reshape(ic_conn_px,:);
[icasig, A, W] = fastica(curr_mov);




%% PCA-ICA combo
file = '/usr/local/lab/People/Andy/Data/AP107/140807/summed_movie/140807_AP107_summed_50.tif';
im = tifread(file);
imageinfo = imfinfo(file);

%% PREPROCESS

[im_filename im_path]= uigetfile('*.tif');
im_fullfilename = [im_path im_filename];
imageinfo = imfinfo(im_fullfilename);

border = 20;
border_idx = border+1:512-border;

im = zeros(length(border_idx)^2,length(imageinfo));

im_lumnorm = zeros(length(border_idx),length(border_idx),length(imageinfo));
smooth_disk = fspecial('gaussian',4,1);
h = fspecial('average', [20 20]);

for i = 1:length(imageinfo)
    curr_im = double(imread(im_fullfilename,'tiff',i,'Info',imageinfo));
    curr_im_smooth = conv2(curr_im,smooth_disk,'same');
    im(:,i) = reshape(curr_im_smooth(border_idx,border_idx),[],1);
    im_lumnorm(:,:,i) = curr_im_smooth(border_idx, border_idx)./ ...
        imfilter(curr_im_smooth(border_idx, border_idx),h);
    disp(i);
end
im_max = max(im,[],2);
im_lumnorm_diff = diff(im_lumnorm.^3,[],3);

im_data = im_lumnorm;%.^3;

down_sample = [4 4 1];
which_pcs = 1:200;

im_x = 1:size(im_data, 2);
im_y = 1:size(im_data, 1);
im_z = 1:size(im_data, 3);

self.pp.im_y = im_y(1:down_sample(1):length(im_y));
self.pp.im_x = im_x(1:down_sample(2):length(im_x));
self.pp.im_z = im_z(1:down_sample(3):length(im_z));

self.pp.im_data = im_data(self.pp.im_y, self.pp.im_x, self.pp.im_z);
self.pp.data = reshape(self.pp.im_data, [], size(self.pp.im_data, 3));
%% PCA

[pca.scores, pca.pcs, pca.eigs] = princomp(self.pp.data');

figure
for i = 1:100;
subplot(10,10,i)
imagesc(reshape(pca.scores(:,i),sqrt(length(pca.scores)),[]));colormap(gray)
axis off
end

%% ICA (SPATIAL)
which_pcs = which_pcs;
which_pcs(which_pcs > size(pca.scores, 2)) = [];
[icasig, A, W] = fastica(pca.scores(:, which_pcs)');

ica.scores = icasig';

% Then we want the component scores to all skew positively.
is_neg_skew = skewness(ica.scores) < 0;
ica.scores(:, is_neg_skew) = -ica.scores(:, is_neg_skew);
A(:, is_neg_skew) = -A(:, is_neg_skew);

ica.A = A;
ica.W = W;
ics = pca.pcs(:, which_pcs) * A;

figure
for i = 1:100;
subplot(10,10,i)
imagesc(reshape(ica.scores(:,i),sqrt(size(ica.scores,1)),[]));colormap(gray)
axis off
end

%% ICA (TEMPORAL)
which_pcs = which_pcs;
which_pcs(which_pcs > size(pca.scores, 2)) = [];
[icasig, A, W] = fastica(pca.pcs(:,which_pcs)');

ica.scores = icasig';

% Then we want the component scores to all skew positively.
is_neg_skew = skewness(ica.scores) < 0;
ica.scores(:, is_neg_skew) = -ica.scores(:, is_neg_skew);
A(:, is_neg_skew) = -A(:, is_neg_skew);

ica.A = A;
ica.W = W;
ics = pca.scores(:, which_pcs) * A;

figure
for i = 1:100;
subplot(10,10,i)
imagesc(reshape(ics(:,i),sqrt(size(ics,1)),[]));colormap(gray)
axis off
end

%% ICA (SPATIOTEMPORAL)
% it doesn't work here? not under tolerance within limits

pca.im_data = reshape(pca.scores,sqrt(length(pca.scores)), ...
    sqrt(length(pca.scores)),[]);

s.mu = 0;
% spatiotemporal ica
[ics, im_data, A, niter] = CellsortICA(pca.pcs(:, which_pcs)', ...
    pca.im_data(:,:,which_pcs), pca.eigs(which_pcs), [], s.mu);



%% Zebra Auto-ROI (currently using)

[im_file im_path] = uigetfile('*.tif');
im_fullfile = [im_path im_file];
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info));
for i = 1:length(im_info)
    im(:,:,i) = imread(im_fullfile,'tiff',i,'Info',im_info);
    disp(i) 
end

numframes = size(im,3);

crop_px = 20;
crop_idx = crop_px+1:512-crop_px;

im = im(crop_idx,crop_idx,:);

disk_rad = 3;
gaussian_rad_std = [6 1];

% make the luminance normalization filter, filter the avg
h_disk = fspecial('disk', disk_rad);
im_avg = mean(im,3);
im_avg_smooth = double(imfilter(im_avg,h_disk));

% smooth/lumnorm image with gaussian filter
h_gaussian = fspecial('gaussian',gaussian_rad_std(1),gaussian_rad_std(2));
disp('Finished (%):')
fprintf('%2d',0)
im_processed = zeros(size(im));
im_zebraroi = zeros(size(im_avg),'double');
for i = 1:numframes-1
    curr_frame = zeros(size(im_avg),'double');
    % this was the original zebra roi, but highlights blood vessels etc
    %curr_frame = (double(im(:,:,i))-double(im_avg))./(double(im_avg)+double(im_avg_smooth));
    curr_frame = (double(im(:,:,i))-double(im_avg_smooth))./(double(im_avg_smooth));
    curr_frame = curr_frame.^3;
    curr_frame = imfilter(curr_frame,h_gaussian);
    
    im_processed(:,:,i) = curr_frame;
    
    im_zebraroi = im_zebraroi+curr_frame./numframes;
    fprintf('%c%c%2d',8,8,round(100*i/numframes));
end

thresh_cutoff = 2*std(im_zebraroi(:));
im_zebraroi_thresh = im_zebraroi > thresh_cutoff;

% force borders to zero
border_idx = true(size(im_zebraroi_thresh));
border_idx(11:end-10,11:end-10) = false;
im_zebraroi_thresh(border_idx) = false;

% get rid of clusters of less than n pixels
smallest_roi = 10;
cc = bwconncomp(im_zebraroi_thresh);
rois = cc.PixelIdxList;
rois(cellfun(@length,rois) < smallest_roi) = [];

im_processed = reshape(im_processed,size(im,1)^2,size(im,3));

roi_trace = cell2mat(cellfun(@(x) nanmean(im_processed(x,:),1),rois','uni',false));

roi_norm = bsxfun(@minus,roi_trace,mean(roi_trace,2));
pixel_norm = bsxfun(@minus,im_processed,mean(im_processed,2));

roi_norm = bsxfun(@times,roi_norm,1./sqrt(sum(roi_norm.^2,2)));
pixel_norm = bsxfun(@times,pixel_norm,1./sqrt(sum(pixel_norm.^2,2)));

roi_pixel_corr = roi_norm*pixel_norm';

figure;
keep_rois = false(length(rois),1);
for i = 1:size(roi_trace,1)
    
    roi_plot = zeros(size(im,1),size(im,2));
    roi_plot(rois{i}) = 1;
    roi_plot(:,:,2:3) = 0;
    
    roi_alpha = zeros(size(im,1),size(im,2));
    roi_alpha(rois{i}) = 1;
    
    subplot(2,1,1);
    plot(roi_trace(i,:),'k');
    
    subplot(2,2,3);colormap(gray);
    set(gca,'YDir','reverse')
    imagesc(im_avg);caxis([0 500]);
    hold on;
    r = imagesc(roi_plot);
    set(r,'AlphaData',roi_alpha);
    xlim([0 size(im,2)]);
    ylim([0 size(im,1)]);
    hold off;
    
    subplot(2,2,4);colormap(gray);
    set(gca,'YDir','reverse')
    imagesc(reshape(roi_pixel_corr(i,:),sqrt(size(roi_pixel_corr,2)),[]));colormap(gray)
    hold on;
    r = imagesc(roi_plot);
    set(r,'AlphaData',roi_alpha);
    xlim([0 size(im,2)]);
    ylim([0 size(im,1)]);
    hold off;
    
    curr_roi_keep = input('Keep roi?' ,'s');
    if strcmp(curr_roi_keep,'y')
        keep_rois(i) = true;
    end
end

use_rois = rois(keep_rois);

% Group ROIs which are on the same axon

% Make binary masks for averaging
roi_masks = zeros(512,512,size(use_rois,1));
for i = 1:length(use_rois)
   curr_mask = zeros(size(im,1),size(im,2),1);
   curr_mask(use_rois{i}) = true;
   roi_masks(crop_idx,crop_idx,i) = curr_mask;
end
polygon.autosort = roi_masks;

% Make outlined polygons for viewing
polygon.ROI = cell(size(polygon.autosort,3),1);
for i = 1:size(polygon.autosort,3)
    first_nonzero = find(polygon.autosort(:,:,i) > 0,1);
    [y_nonzero x_nonzero] = ind2sub([size(polygon.autosort,1), ...
        size(polygon.autosort,2)],first_nonzero);
    roi_edge = bwtraceboundary(polygon.autosort(:,:,i), ...
        [y_nonzero x_nonzero],'N'); % find boundary, in order   
    polygon.ROI{i} = fliplr(roi_edge); % make x,y coordinates
end



















