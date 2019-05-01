function AP_zAlign_daily

% This function is designed to give an estimate of current day's alignment
%
% Input: select 1) an image from today, and 2) your reference images
% *NOTE: reference image names must be in the form:
% '(any filename)_zdistance.tif', i.e. '120612_reference_+10.tif'
%
% Output: a plot showing correlations and offsets for each reference image


% pick today's image
[orig_file orig_path] = uigetfile('*.tif','Today''s image');

% pick reference images
[ref_file ref_path] = uigetfile('*.tif','Reference images','Multiselect','On');

% get the z's from the reference image names
ref_zs = zeros(length(ref_file),1);
for i = 1:length(ref_file)
    underscore_indx = strfind(ref_file{i},'_');
    ref_zs(i) = str2num(ref_file{i}(underscore_indx(end)+1:end-4));
end
    

% check whether one or two channels
orig_info = imfinfo([orig_path orig_file]);
num_channels = length(orig_info);
M=orig_info(1).Width;
N=orig_info(1).Height;
aspect_ratio = N/M;

% set up blur filter
cell_size = 5;
blur = fspecial('average',[round(cell_size*aspect_ratio) cell_size]);

% load in original image, luminance normalize
im_orig = zeros(N,M,num_channels);
im_orig_norm = zeros(N,M,num_channels);
for i = 1:num_channels
    im_orig(:,:,i) = double(imread([orig_path orig_file],'tiff',i));
    im_orig_norm(:,:,i) = im_orig(:,:,i)./imfilter(im_orig(:,:,i),blur);
end

num_refs = length(ref_file);
% load and align each channel
zcorr = zeros(num_refs,num_channels);
xoffset = zeros(num_refs,num_channels);
yoffset = zeros(num_refs,num_channels);
border = 40;
im_orig_border_x = border;
im_orig_border_y = round(border*aspect_ratio);
for channel = 1:num_channels
    % load in references of that channel  
    for curr_ref = 1:num_refs;
        im_orig_norm_curr = im_orig_norm(:,:,channel);
        
        im_ref = zeros(N,M,1);
        im_ref = double(imread([ref_path ref_file{curr_ref}],'tiff',channel));
        im_ref_norm = im_ref./imfilter(im_ref,blur);
        
        im_orig_norm_curr = im_orig_norm_curr( ...
            im_orig_border_y+1:end-im_orig_border_y, ...
            im_orig_border_x+1:end-im_orig_border_x);
        
        temp_xcorr = [];
        temp_xcorr = normxcorr2(im_orig_norm_curr,im_ref_norm);
        
        % check for max value only within max xy displacement range
        xcorr_midpoint_x = round(size(temp_xcorr,2)/2);
        xcorr_midpoint_y = round(size(temp_xcorr,1)/2);
        temp_xcorr_cutoff = true(size(temp_xcorr));
        temp_xcorr_cutoff(xcorr_midpoint_y-im_orig_border_y+2: ...
            xcorr_midpoint_y+im_orig_border_y-1, ...
            xcorr_midpoint_x-im_orig_border_x+2: ...
            xcorr_midpoint_x+im_orig_border_x-1) = false;
        temp_xcorr(temp_xcorr_cutoff) = 0;
        [max_c,imax] = max(temp_xcorr(:));
        [ypeak xpeak] = ind2sub(size(temp_xcorr),imax(1));
        corr_offset = [(xpeak-size(im_orig_norm_curr,2)-im_orig_border_x) ...
            (ypeak-size(im_orig_norm_curr,1)-im_orig_border_x)];
        
        zcorr(curr_ref,channel) = max_c;
        xoffset(curr_ref,channel) = corr_offset(1);
        yoffset(curr_ref,channel) = corr_offset(2);
    end
    disp(['Finished Ch ' num2str(channel)])
end

% Plot z correlations, xy offsets
figure; hold on;

plot_zcorr = [ref_zs zcorr];
plot_zcorr = sortrows(plot_zcorr,1);
plot_xoffset = [ref_zs xoffset];
plot_xoffset = sortrows(plot_xoffset,1);
plot_yoffset = [ref_zs yoffset];
plot_yoffset = sortrows(plot_yoffset,1);

subplot(3,1,1)
plot(plot_zcorr(:,1),plot_zcorr(:,2:end),'Linewidth',3);
title('Correlation');
line([0 0],ylim,'linestyle','--','color','k')

subplot(3,1,2);
plot(plot_xoffset(:,1),plot_xoffset(:,2:end),'Linewidth',3);
title('X Offset');
a = ylim;
ylim(a+[-3 3])
line([0 0],ylim,'linestyle','--','color','k')
line(xlim,[0 0],'linestyle','--','color','k')

subplot(3,1,3);
plot(plot_yoffset(:,1),plot_yoffset(:,2:end),'Linewidth',3);
title('Y Offset');
a = ylim;
ylim(a+[-3 3])
line([0 0],ylim,'linestyle','--','color','k')
line(xlim,[0 0],'linestyle','--','color','k')






