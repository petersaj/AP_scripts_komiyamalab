%% Load in video 

[im_file im_path] = uigetfile('*.tif');
im_fullfile = [im_path im_file];
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info),'uint16');
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
file = '/usr/local/lab/People/Andy/Data/AP122/150208/summed/150822_AP122_summed_50.tif';
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



%% Zebra Auto-ROI

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

thresh_cutoff = 1*std(im_zebraroi(:));
im_zebraroi_thresh = im_zebraroi > thresh_cutoff;

% force borders to zero
border_idx = true(size(im_zebraroi_thresh));
border_idx(11:end-10,11:end-10) = false;
im_zebraroi_thresh(border_idx) = false;

% get rid of clusters of less than n pixels
smallest_roi = 1;
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
    imagesc(im_avg);caxis([0 1000]);
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

rois = rois(keep_rois);

roi_masks = zeros(512,512,size(rois,1));
for i = 1:length(rois)
   curr_mask = zeros(size(im,1),size(im,2),1);
   curr_mask(rois{i}) = true;
   roi_masks(crop_idx,crop_idx,i) = curr_mask;
end

polygon.ROI = nan(length(rois),1);
polygon.autosort = roi_masks;

%% Zebra auto detect (no manual stuff)

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
im_processed = zeros(size(im,1)*size(im,2),size(im,3));
im_zebraroi = zeros(size(im_avg),'double');
for i = 1:numframes-1
    curr_frame = zeros(size(im_avg),'double');
    % this was the original zebra roi, but highlights blood vessels etc
    %curr_frame = (double(im(:,:,i))-double(im_avg))./(double(im_avg)+double(im_avg_smooth));
    curr_frame = (double(im(:,:,i))-double(im_avg_smooth))./(double(im_avg_smooth));
    curr_frame = curr_frame.^3;
    %curr_frame = imfilter(curr_frame,h_gaussian);
    
    im_processed(:,i) = curr_frame(:);
    
    im_zebraroi = im_zebraroi+curr_frame./numframes;
    fprintf('%c%c%2d',8,8,round(100*i/numframes));
end

thresh_cutoff = 1*std(im_zebraroi(:));
im_zebraroi_thresh = im_zebraroi > thresh_cutoff;

% force borders to zero
border_idx = true(size(im_zebraroi_thresh));
border_idx(11:end-10,11:end-10) = false;
im_zebraroi_thresh(border_idx) = false;

% get rid of clusters of less than n pixels
smallest_roi = 3;
cc = bwconncomp(im_zebraroi_thresh);
rois = cc.PixelIdxList;
rois(cellfun(@length,rois) < smallest_roi) = [];

% "correlation" map of ROIs (or just z-scored dot product)
im_z = zscore(im,[],3);
im_reshape = reshape(im,[],size(im,3));
roi_trace_z = zscore(cell2mat(cellfun(@(x) nanmean(im_reshape(x,:),1),rois','uni',false)));

correlation_map = roi_trace_z*reshape(im_z,[],size(im_z,3))';

for curr_roi = 1:length(rois)
    
end


%% Detect with velocity (doesn't really look much better)

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

% make the luminance normalization filter, filter the avg
h_disk = fspecial('disk', disk_rad);
im_avg = mean(im,3);
im_avg_smooth = double(imfilter(im_avg,h_disk));

im_processed = cell2mat(arrayfun(@(x) reshape((double(im(:,:,x)) - ...
    double(im_avg_smooth))./(double(im_avg_smooth)),[],1).^3,1:size(im,3),'uni',false));


thresh_cutoff = 1*std(im_zebraroi(:));
im_zebraroi_thresh = im_zebraroi > thresh_cutoff;

% force borders to zero
border_idx = true(size(im_zebraroi_thresh));
border_idx(11:end-10,11:end-10) = false;
im_zebraroi_thresh(border_idx) = false;

% get rid of clusters of less than n pixels
smallest_roi = 1;
cc = bwconncomp(im_zebraroi_thresh);
rois = cc.PixelIdxList;
rois(cellfun(@length,rois) < smallest_roi) = [];

% "correlation" map of ROIs (or just z-scored dot product)
im_z = zscore(im,[],3);
roi_trace_z = zscore(cell2mat(cellfun(@(x) nanmean(im_processed(x,:),1),rois','uni',false)));

correlation_map = roi_trace_z*reshape(im_z,[],size(im_z,3))';

for curr_roi = 1:length(rois)
    
end


%% Try to get rid of point spread

[im_file im_path] = uigetfile('*.tif');
im_fullfile = [im_path im_file];
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info));
for i = 1:length(im_info)
    im(:,:,i) = imread(im_fullfile,'tiff',i,'Info',im_info);
    disp(i) 
end

numframes = size(im,3);

% bin the fluorescence values
n_bins = 10;
fluor_bins = linspace(min(im(:)),max(im(:)),n_bins);

neighbor_pixels = 20;

neighbor_fluor = zeros(n_bins,neighbor_pixels*2+1);
bin_members = zeros(n_bins,1);

for z = 1:size(im,3);
    for y = 1:size(im,2);
        for x = neighbor_pixels+1:size(im,2)-(neighbor_pixels+1)
            
            curr_fluor = im(y,x,z);
            
            if curr_fluor > 0;
                curr_bin = find(histc(curr_fluor,fluor_bins));
                
                bin_members(curr_bin) = bin_members(curr_bin) + 1;
                
                neighbor_fluor(curr_bin,:) = neighbor_fluor(curr_bin,:) + ...
                    im(y,x-neighbor_pixels:x+neighbor_pixels,z);
            end
            
        end
    end
    disp(z);
end

neighbor_fluor = bsxfun(@times,neighbor_fluor,1./bin_members);
neighbor_fluor_norm = cell2mat(arrayfun(@(x) (neighbor_fluor(x,:)- ...
    min(neighbor_fluor(x,:)))./max(neighbor_fluor(x,:)- ...
    min(neighbor_fluor(x,:))),[1:n_bins]','uni',false));

figure; 
subplot(2,1,1); hold on;
col = jet(n_bins);
for i = 1:n_bins
   plot(neighbor_fluor(i,:),'color',col(i,:));    
end
subplot(2,1,2); hold on;
col = jet(n_bins);
for i = 1:n_bins
   plot(neighbor_fluor_norm(i,:),'color',col(i,:));    
end



% try deconvolution

pointspread_std = 2;
filter_size = pointspread_std*2*2+1;
pointspread_filter = fspecial('gaussian',filter_size,pointspread_std);

im_test = im(:,:,100);

[m,n] = size(im_test);
[mb,nb] = size(pointspread_filter);

mm = m + mb - 1;
nn = n + nb - 1;

d = ifft2(fft2(im_test,mm,nn) ./ fft2(pointspread_filter,mm,nn));



%% **CURRENTLY USING**: Draw ROIs on aligned images of all sessions

animal = 'AP135';
save_flag = true;

summed_path = ['C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\' animal];
summed_files = dir([summed_path filesep '*.tif']);
summed_filenames = sort({summed_files.name});
n_sessions = length(summed_filenames);

M = 512;
N = 512;

roi_size = 15;

% To fine-align by cross-correlation (I actually think this makes it worse)
xc_align = false;

disp('Summing activity within sessions')

summed_activity = double(nan(M,N,n_sessions));
summed_max = double(nan(M,N,n_sessions));
summed_frames = nan(n_sessions,1);
% Loop through sessions, get 1) maximum projection, 2) summed activity
for curr_session = 1:n_sessions
    
    % Get max projection of day
    curr_file = [summed_path filesep summed_filenames{curr_session}];
    
    curr_summed_activity = zeros(M,N);
    curr_max = nan(M,N);
    im_info = imfinfo(curr_file);
    n_frames = length(im_info);
    
    % Load in current day    
    im = zeros(M,N,n_frames,'uint16');
    for i = 1:n_frames
        im(:,:,i) = imread(curr_file,'tiff',i,'Info',im_info);
    end
      
    % Make blurred average
    curr_avg = nanmean(im,3);
    h_disk = fspecial('disk', roi_size);
    im_avg_smooth = double(imfilter(curr_avg,h_disk));
    % Sum activity
    for i = 1:n_frames
        curr_frame = zeros(size(curr_avg),'double');
        curr_frame = (double(im(:,:,i))-double(im_avg_smooth))./(double(im_avg_smooth));
        curr_frame = curr_frame.^3;
        
        curr_summed_activity = curr_summed_activity+curr_frame./n_frames;
    end
    
    summed_frames(curr_session) = n_frames;
    summed_activity(:,:,curr_session) = curr_summed_activity;

    summed_max(:,:,curr_session) = max(im,[],3);
    disp(['Session: ' num2str(curr_session)]);
    
    clear im
end

% Register days to each other
disp('Registering maximum projections')
tform_matrix = cell(n_sessions,1);
tform_matrix{1} = eye(3);

summed_max_reg = nan(size(summed_max));
summed_max_reg(:,:,1) = summed_max(:,:,1);
for curr_session = 2:n_sessions

    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.02;%0.0001;
    %optimizer.RelaxationFactor = 0.1;
    %optimizer.GradientMagnitudeTolerance = 1e-6;
    %optimizer.MaximumIterations = 300;
    
    tformEstimate = imregtform(summed_max(:,:,curr_session),summed_max(:,:,1),'affine',optimizer,metric);
    %tformEstimate = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');
    
    im_reg = imwarp(summed_max(:,:,curr_session),tformEstimate,'Outputview',imref2d([M N]));
    
    tform_matrix{curr_session} = tformEstimate.T;
    summed_max_reg(:,:,curr_session) = im_reg;
    
end

% Human-verify registration (if it isn't near-perfect, don't use)
gui_fig = AP_image_scroll(summed_max_reg);
check_alignment = input('Perfect alignment (y/n)?: ','s');

if strcmp(check_alignment,'y')
    
    % Align summed activity
    disp('Aligning activity projections')
    summed_activity_reg = nan(size(summed_activity));
    for curr_session = 1:n_sessions
        curr_tform = affine2d;
        curr_tform.T = tform_matrix{curr_session};
        
        im_reg = imwarp(summed_activity(:,:,curr_session),curr_tform,'Outputview',imref2d([M N]));
        summed_activity_reg(:,:,curr_session) = im_reg;
    end
    
    % Scale each session, get sum total activity
    summed_activity_scaled = arrayfun(@(x) ...
        summed_activity_reg(:,:,x).*summed_frames(x),1:n_sessions,'uni',false);
    summed_activity_total = sum(cat(3,summed_activity_scaled{:}),3)./sum(summed_frames);
    
    % Force borders to zero
    x_crop_px = 25;
    y_crop_px = 15;
    border_idx = true(size(summed_activity_total));
    border_idx(y_crop_px:end-y_crop_px+1,x_crop_px:end-x_crop_px+1) = false;
    
    summed_activity_total(border_idx) = 0;
    
    disp('Thresholding and watershedding')
    % Set threshold for detecting active points
    thresh_cutoff = 1*std(summed_activity_total(:));
    activity_thresh = summed_activity_total > thresh_cutoff;
    
    % Get rid of clusters of less than n pixels, dilate remaining
    dilate_matrix = ones(5);
    smallest_roi = 1;
    activity_thresh_clean = imdilate(bwareaopen(activity_thresh,smallest_roi),dilate_matrix);
    
    % Split clusters by local maxima to separate neighboring dendrites
    local_max = imextendedmax(summed_activity_total,thresh_cutoff/2);
    im_inv = imcomplement(summed_activity_total);
    im_inv_minrestrict = imimposemin(im_inv, ~activity_thresh_clean | local_max);
    im_watershed = watershed(im_inv_minrestrict);
    activity_thresh_split = im_watershed > 1;
    
    save_path = [summed_path filesep animal '_batch_thresh_roi'];
    if ~exist(save_path,'dir')
        mkdir(save_path)
    end
    
    % Define surrounding pixels to pull (needs a lot because sparse)
    disp('Creating ROIs (%):')
    fprintf('%2d',0)
    surround_px = 20;
    for curr_session = 1:n_sessions;
        
        % get connected blobs as single ROIs
        cc = bwconncomp(activity_thresh_split);
        rois = cc.PixelIdxList;
        curr_session_rois = cell(size(rois));
        template_roicenters = nan(length(rois),2);
        
        % make reverse registration transform matrix
        curr_tform = affine2d;
        inverse_tform = inv(tform_matrix{curr_session});
        % inverse is not perfect: force last column
        inverse_tform(1:2,3) = 0;
        inverse_tform(3,3) = 1;
        curr_tform.T = inverse_tform;
        
        % prepare offset matrix from xcorr
        roi_xcorr_offsets = zeros(length(rois),2);
        
        for curr_roi = 1:length(rois);
            
            binaryCell = zeros(size(activity_thresh_split(:,:,1)));
            binaryCell(rois{curr_roi}) = 1;
            
            % reverse register cell to fit session
            binaryCell_reg = +(imwarp(binaryCell,curr_tform,'Outputview',imref2d([M N])) > 0);
            
            first_nonzero = find(binaryCell_reg > 0,1);
            [y_nonzero x_nonzero] = ind2sub([M N],first_nonzero);
            roi_edge = bwtraceboundary(binaryCell_reg,[y_nonzero x_nonzero],'N'); % find boundary, in order
            
            % get difference in distance to center across border
            % function from file exchange
            if size(roi_edge,1) > 10
                num_verticies = 10;
                roi_edge_downsample = reduce_poly(roi_edge',num_verticies)';
            end
            
            % fine-tune position by local correlation compared to day 1 (unless
            % it is day 1, in which case store the ROI centers for local
            % comparison)
            if xc_align
                if curr_session == 1;
                    [~,curr_x,curr_y] = polycenter(roi_edge_downsample(:,2),roi_edge_downsample(:,1));
                    template_roicenters(curr_roi,:) = round([curr_x curr_y]);
                else
                    % get current ROI center
                    [~,curr_x,curr_y] = polycenter(roi_edge_downsample(:,2),roi_edge_downsample(:,1));
                    curr_x = round(curr_x);
                    curr_y = round(curr_y);
                    
                    % grab template ROI image, find location of max cross correlation
                    template_x = template_roicenters(curr_roi,1) - surround_px : ...
                        template_roicenters(curr_roi,1) + surround_px;
                    template_y = template_roicenters(curr_roi,2) - surround_px : ...
                        template_roicenters(curr_roi,2) + surround_px;
                    
                    use_x = template_x > 0 & template_x < N;
                    use_y = template_y > 0 & template_y < M ;
                    
                    % luminance normalize images
                    h = fspecial('disk', roi_size);
                    template_lumnorm = summed_max(:,:,1)./imfilter(summed_max(:,:,1),h);
                    curr_lumnorm = summed_max(:,:,curr_session)./imfilter(summed_max(:,:,curr_session),h);
                    
                    template_roi_max = template_lumnorm(template_y(use_y),template_x(use_x),1);
                    
                    xc = normxcorr2(template_roi_max,curr_lumnorm);
                    
                    % get the point of max correlation
                    xcorr_midpoint = round(length(xc)/2);
                    [max_c,imax] = max(xc(:));
                    [ypeak xpeak] = ind2sub(size(xc),imax(1));
                    corr_offset = [(xpeak-surround_px)-curr_x ...
                        (ypeak-surround_px)-curr_y];
                    
                    roi_xcorr_offsets(curr_roi,:) = corr_offset;
                end
            end
            
            % make x,y coordinates
            curr_session_rois{curr_roi} = fliplr(roi_edge_downsample);
            
        end
        
        % don't use offsets past a reasonable moving distance
        dist_thresh = 10;
        offset_dist = sqrt(sum(roi_xcorr_offsets.^2,2));
        roi_xcorr_offsets(offset_dist > dist_thresh,:) = 0;
        
        % apply offsets from correlation to ROIs
        curr_session_rois_fixed = arrayfun(@(x) curr_session_rois{x} + ...
            repmat(roi_xcorr_offsets(x,:),size(curr_session_rois{x},1),1), ...
            1:length(curr_session_rois),'uni',false);
        
        clear polygon.ROI
        polygon.ROI = curr_session_rois_fixed;
        
        % save
        if save_flag
            curr_save_filename = [save_path filesep summed_filenames{curr_session}(1:end-4) '.roi'];
            save(curr_save_filename,'polygon');
        end
        
        fprintf('%c%c%2d',8,8,round(100*curr_session/n_sessions));
    end
    
    disp('Done')
    
end

%%% Human verification of ROIs

% Draw ROIs on the aligned images for inspection
cc = bwconncomp(activity_thresh_split);
rois = cc.PixelIdxList;
figure(gui_fig)
for curr_roi = 1:length(rois)
    % Get verticies for ROI
    binaryCell = zeros(size(activity_thresh_split(:,:,1)));
    binaryCell(rois{curr_roi}) = 1;
    
    first_nonzero = find(binaryCell > 0,1);
    [y_nonzero x_nonzero] = ind2sub([M N],first_nonzero);
    roi_edge = bwtraceboundary(binaryCell,[y_nonzero x_nonzero],'N'); 
    
    roi_edge(end+1,:) = roi_edge(1,:);
    
    line(roi_edge(:,2),roi_edge(:,1));
    
end

% Draw borders to know what's been excluded
line([x_crop_px,M-x_crop_px+1],[y_crop_px,y_crop_px],'color','r');
line([x_crop_px,x_crop_px],[y_crop_px,N-y_crop_px+1],'color','r');
line([x_crop_px,M-x_crop_px+1],[N-y_crop_px+1,N-y_crop_px+1],'color','r');
line([M-x_crop_px+1,M-x_crop_px+1],[y_crop_px,N-y_crop_px+1],'color','r');





