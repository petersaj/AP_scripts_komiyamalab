%% Test background interpolation subtraction on toy data

a = zeros(50,100,100);
t1 = sin(linspace(0,3*pi,100))*0.5;
t2 = sin(linspace(0,5*pi,100));

a(25,25,:) = t1;
a(25,75,:) = t2;

rois = zeros(size(a,1),size(a,2),2);
rois(25,25,1) = 1;
rois(25,75,2) = 1;
dilate_filt = fspecial('disk',10);
rois = convn(rois,dilate_filt,'same') > 0;

psf = fspecial('gaussian',200,50);
a2 = convn(a,psf,'same');

a2_reshape = reshape(a2,[],size(a2,3));

% Raw ROI values
t1_r = nanmean(a2_reshape(reshape(rois(:,:,1),[],1),:));
t2_r = nanmean(a2_reshape(reshape(rois(:,:,2),[],1),:));

% Interpolate across ROIs
b = a2;
b(repmat(sum(rois,3) > 0,1,1,size(b,3))) = NaN;
for i = 1:size(b,3)
   b(:,:,i) = inpaint_nans(b(:,:,i)); 
end

% Subtract interpolation
is = a2-b;
is_reshape = reshape(is,[],size(is,3));
t1_is = nanmean(is_reshape(reshape(rois(:,:,1),[],1),:));
t2_is = nanmean(is_reshape(reshape(rois(:,:,2),[],1),:));

% Do like old background subtraction
bgrois = zeros(size(a,1),size(a,2),2);
bgrois(25,25,1) = 1;
bgrois(25,75,2) = 1;
dilate_filt = fspecial('disk',16);
bgrois = ((convn(bgrois,dilate_filt,'same') > 0) - rois) > 0;

t1_bg = t1_r - nanmean(a2_reshape(reshape(bgrois(:,:,1),[],1),:));
t2_bg = t2_r - nanmean(a2_reshape(reshape(bgrois(:,:,2),[],1),:));

%% Test background interpolation subtraction on real data

im_fullfile = 'C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP128\150421_AP128_001_001_summed_50.tif';
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info),'uint16');
for i = 1:length(im_info)
    im(:,:,i) = imread(im_fullfile,'tiff',i,'Info',im_info);
    disp(i) 
end

a = load('C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP128\AP128_batch_thresh_roi\150421_AP128_001_001_summed_50.roi','-mat');
curr_polygons = a.polygon.ROI;

curr_masks = cellfun(@(x) poly2mask(x(:,1), ...
        x(:,2),512,512),curr_polygons,'uni',false);
curr_masks_cat = sum(cat(3,curr_masks{:}),3) > 0;

im_bg = zeros(size(im),'uint16');
for i = 1:size(im,3)
    curr_nan_frame = double(im(:,:,i));
    curr_nan_frame(curr_masks_cat) = NaN;
    curr_interp_frame = inpaint_nans(curr_nan_frame,4);
    im_bg(:,:,i) = curr_interp_frame;
    disp(i);
end

im_interpsub = zeros(size(im),'uint16');
for i = 1:size(im,3)
    curr_nan_frame = double(im(:,:,i));
    curr_nan_frame(curr_masks_cat) = NaN;
    curr_interp_frame = inpaint_nans(curr_nan_frame);
    im_interpsub(:,:,i) = im(:,:,i) - uint16(curr_interp_frame);
    disp(i);
end

% Average and interp avg
im_avg = nanmean(double(im),3);

im_avg_roinan = im_avg;
im_avg_roinan(curr_masks_cat) = NaN;
im_avg_interpsub = inpaint_nans(im_avg_roinan);

% Do this like current: subtract background df
im_bgsub = im - (im_bg - repmat(uint16(im_avg_interpsub),1,1,size(im,3)));

im_reshape_df = reshape((double(im) - repmat(im_avg,1,1,size(im,3))) ...
    ./repmat(im_avg,1,1,size(im,3)),[],size(im,3));
im_bgsub_reshape_df = reshape((double(im_bgsub) - repmat(im_avg,1,1,size(im,3))) ...
    ./repmat(im_avg,1,1,size(im,3)),[],size(im,3));

demo_roi = 143;
demo_roi_trace = nanmean(im_reshape_df(reshape(curr_masks{demo_roi},[],1),:));
demo_roi_trace_interpsub = nanmean(im_bgsub_reshape_df(reshape(curr_masks{demo_roi},[],1),:));

figure; hold on;
plot(demo_roi_trace,'k')
plot(demo_roi_trace_interpsub,'r')

% This makes it noisier but also gets rid of artifacts slightly better?
im_bgsub = im - im_bg;

im_reshape_df = reshape((double(im) - repmat(im_avg,1,1,size(im,3))) ...
    ./repmat(im_avg,1,1,size(im,3)),[],size(im,3));
im_bgsub_reshape_df = reshape((double(im_bgsub) - repmat(im_avg_interpsub,1,1,size(im,3))) ...
    ./repmat(im_avg_interpsub,1,1,size(im,3)),[],size(im,3));

demo_roi = 143;
demo_roi_trace = nanmean(im_reshape_df(reshape(curr_masks{demo_roi},[],1),:));
demo_roi_trace_interpsub = nanmean(im_bgsub_reshape_df(reshape(curr_masks{demo_roi},[],1),:));

figure; hold on;
plot(demo_roi_trace,'k')
plot(demo_roi_trace_interpsub,'r')


%% Testing subtraction of interpolated df

im_fullfile = 'C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP128\150421_AP128_001_001_summed_50.tif';
im_info = imfinfo(im_fullfile);

im = zeros(512,512,length(im_info),'uint16');
for i = 1:length(im_info)
    im(:,:,i) = imread(im_fullfile,'tiff',i,'Info',im_info);
    disp(i) 
end

a = load('C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP128\AP128_batch_thresh_roi\150421_AP128_001_001_summed_50.roi','-mat');
curr_polygons = a.polygon.ROI;

curr_masks = cellfun(@(x) poly2mask(x(:,1), ...
        x(:,2),512,512),curr_polygons,'uni',false);
curr_masks_cat = sum(cat(3,curr_masks{:}),3) > 0;
curr_masks_reshape = cell2mat(cellfun(@(x) reshape(x,[],1),curr_masks,'uni',false));

im_bg = zeros(size(im),'uint16');
for i = 1:size(im,3)
    curr_nan_frame = double(im(:,:,i));
    curr_nan_frame(curr_masks_cat) = NaN;
    curr_interp_frame = inpaint_nans(curr_nan_frame,4);
    im_bg(:,:,i) = curr_interp_frame;
    disp(i);
end

% Average and interp avg
im_avg = nanmean(double(im),3);

im_avg_roinan = im_avg;
im_avg_roinan(curr_masks_cat) = NaN;
im_avg_interpsub = inpaint_nans(im_avg_roinan);

% Subtract background df (i.e. only get rid of background signals)
% (sometimes goes negative if neighboring signal is particularly strong)
im_bgsub = im - (im_bg - repmat(uint16(im_avg_interpsub),1,1,size(im,3)));

% NEW TEST: subtract interp, but not MORE THAN ROI df (often doesn't
% subtract enough)
%im_bgsub = im - min(cat(4,(im_bg - repmat(uint16(im_avg_interpsub),1,1,size(im,3))), ...
%    (im - repmat(uint16(im_avg),1,1,size(im,3)))),[],4);

im_reshape_df = reshape((double(im) - repmat(im_avg,1,1,size(im,3))) ...
    ./repmat(im_avg,1,1,size(im,3)),[],size(im,3));
im_bgsub_reshape_df = reshape((double(im_bgsub) - repmat(im_avg,1,1,size(im,3))) ...
    ./repmat(im_avg,1,1,size(im,3)),[],size(im,3));

roi_traces = bsxfun(@times,double(im_reshape_df)'*curr_masks_reshape,1./sum(curr_masks_reshape));
roi_traces_interpsub = bsxfun(@times,double(im_bgsub_reshape_df)'*curr_masks_reshape,1./sum(curr_masks_reshape));

% Find the difference in cells
[~,roi_diff_sort] = sort(nanmean(roi_traces - roi_traces_interpsub,1));

% this comparison actually works better with correlation between the two
% traces










