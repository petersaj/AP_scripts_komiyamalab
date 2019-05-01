function AP_zAlign(animal)
% AP_zAlign(animal)
% Align z-stack to all movies for one animal

% find days
data_path = ['/usr/local/lab/People/Andy/Data/' animal];
data_dir = dir(data_path);
data_names = {data_dir.name};
data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
    'UniformOutput',false);
data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
data_folders = data_fullnames(data_isfolders);
data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);

zalign = struct;
skip_points = 15;

for curr_day = 1:length(days)
    clearvars -except animal days zalign curr_day skip_points
    
    % pick original image (single slice)
    orig_path = ['/usr/local/lab/People/Andy/Data/' animal filesep days{curr_day} filesep];
    orig_file = [days{curr_day} '_' animal '_summed_50'];
    %[orig_file orig_path] = uigetfile('*.tif','Pick decimated stack');
    
    % pick stack image (multiple slices)
    stack_path = ['/usr/local/lab/People/Andy/Data/' 'z_alignment_stacks' filesep];
    stack_file = [animal '_zstack'];
    %[stack_file stack_path] = uigetfile('*.tif','Pick stack for alignment');
    
    % load in original image
    stack_info = imfinfo([orig_path orig_file],'tiff');
    M=stack_info(1).Width;
    N=stack_info(1).Height;
    stack_num = length(stack_info);
    im_orig_full = zeros(N,M,stack_num);
    for i = 1:stack_num
        im_orig_full(:,:,i) = double(imread([orig_path orig_file],'tiff',i));
    end
    
    % load in stack
    stack_info = imfinfo([stack_path stack_file],'tiff');
    M=stack_info(1).Width;
    N=stack_info(1).Height;
    stack_num = length(stack_info);
    im_stack = zeros(N,M,stack_num/2);
    
    % luminance normalize stack: divide by averaged area within half-cell
    % length radius (hardcoded as 5 pixels at the moment ~ 5 microns)
    blur = fspecial('average',[5 5]);
    im_stack_norm = zeros(N,M,stack_num/2);
    
    % assume 2 channels for stack
    for i = 1:2:stack_num
        temp_im = [];
        temp_im = double(imread([stack_path stack_file],'tiff',i));
        im_stack(:,:,ceil(i/2)) = temp_im;
        im_stack_norm(:,:,ceil(i/2)) = temp_im./imfilter(temp_im,blur);
    end   
    
    num_z_estimates = floor(size(im_orig_full,3)/skip_points);
      
    max_c_all = cell(num_z_estimates,1);
    corr_offset_all = cell(num_z_estimates,1);
   
    for j = 1:num_z_estimates;
        curr_frame = 1+(j-1)*skip_points;
        im_orig = im_orig_full(:,:,curr_frame);
        
        % luminance normalize image: divide image by averaged area within
        % 2-cell length radius (hardcoded as 20 pixels at the moment)
        blur = fspecial('disk',20);
        im_orig_norm = im_orig./imfilter(im_orig,blur);
        
        % try just doing 2D Corrcoef (xcorr2) on some inside square
        % this assumes that there isn't large rotation
        num_slices = size(im_stack,3);
        im_orig_border = 40;
        im_orig_clipped = im_orig(im_orig_border+1:end-im_orig_border,...
            im_orig_border+1:end-im_orig_border);
        im_orig_norm_clipped = im_orig_norm(im_orig_border+1:end-im_orig_border,...
            im_orig_border+1:end-im_orig_border);
        for i = 1:num_slices
            temp_xcorr = [];
            
            % align by edges (high freq) to prevent aligning of general
            % diffuse neuropil signal
%            blur = fspecial('disk',20);
                        
%             hy = fspecial('sobel');
%             hx = hy';
%             im1y = imfilter(im_orig_clipped,hy,'replicate');
%             im1x = imfilter(im_orig_clipped,hx,'replicate');
%             im1edge = sqrt(im1x.^2 + im1y.^2);
%             im1edgeblur = imfilter(im1edge,blur);
%                        
%             im2y = imfilter(im_stack(:,:,i),hy,'replicate');
%             im2x = imfilter(im_stack(:,:,i),hx,'replicate');
%             im2edge = sqrt(im2x.^2 + im2y.^2);
%             im2edgeblur = imfilter(im2edge,blur);
            
            %temp_xcorr = normxcorr2(im1edgeblur,im2edgeblur);
            
            % try aligning by luminance normalized images only
            temp_xcorr = normxcorr2(im_orig_norm_clipped,im_stack_norm(:,:,i));

%            temp_xcorr = normxcorr2(im_orig_clipped,im_stack(:,:,i));
            % check for max value only within max xy displacement range
            xcorr_midpoint = round(length(temp_xcorr)/2);
            xcorr_midpoint = xcorr_midpoint;
            temp_xcorr_cutoff = true(size(temp_xcorr));
            temp_xcorr_cutoff(xcorr_midpoint-im_orig_border+2: ...
                xcorr_midpoint+im_orig_border-1, ...
                xcorr_midpoint-im_orig_border+2: ...
                xcorr_midpoint+im_orig_border-1) = false;
            temp_xcorr(temp_xcorr_cutoff) = 0;
            [max_c,imax] = max(temp_xcorr(:));
            [ypeak xpeak] = ind2sub(size(temp_xcorr),imax(1));
            corr_offset = [(xpeak-size(im_orig_norm_clipped,2)-im_orig_border) ...
                (ypeak-size(im_orig_norm_clipped,1)-im_orig_border)];
            
            max_c_all{j}(i) = max_c;
            corr_offset_all{j}(i,:) = corr_offset;
        end
        disp([num2str(j) '/' num2str(num_z_estimates) ' of day ' days{curr_day}]);
    end
    
    % find the max of the xcorr, try smoothed, smoothed fit and quadradic fit
    smooth_factor = 10;
    z_slice_smooth = cellfun(@(x) find(smooth(x,smooth_factor,'lowess') == max(smooth(x,smooth_factor,'lowess')),1),max_c_all);
    z_slice_smoothcurve = cellfun(@(x) smooth(x,smooth_factor,'lowess'),max_c_all,'UniformOutput',false);
    z_slice_smoothfit = zeros(length(z_slice_smoothcurve),1);
    for i = 1:size(z_slice_smoothcurve);
        curr_curve = z_slice_smoothcurve{i};
        curr_fit = fit([1:length(curr_curve)]',curr_curve,'cubicinterp');
        curr_fit_interp = curr_fit([1:0.01:length(curr_curve)]);
        % multiply by 100 for interpolation, add 1 for starting frame
        z_slice_smoothfit(i) = find(curr_fit_interp == max(curr_fit_interp))/100+1;   
    end
    b = cellfun(@(x) polyfit([1:length(x)],x,2),max_c_all,'UniformOutput',false);
    z_slice_quadfit = cellfun(@(x) roots([2*x(1) x(2)]),b);
    
    % save in a structure
    zalign.day{curr_day} = days{curr_day};
    zalign.max_c_all{curr_day} = max_c_all;
    zalign.corr_offset_all{curr_day} = corr_offset_all;
    zalign.z_slice_smooth{curr_day} = z_slice_smooth;
    zalign.z_slice_quadfit{curr_day} = z_slice_quadfit;
    zalign.z_slice_smoothcurve{curr_day} = z_slice_smoothfit;
        
    % % input parameters and minimize 1-corrcoef
    % align_fun = zeros(stack_num,1);
    % estimates = zeros(stack_num,2);
    %
    % w = waitbar(0,'Aligning stack to original')
    % for j = 1:stack_num
    %     im_align_z = im_stack(:,:,j);
    %     starting = [0 0]; % dx dy
    %     options = optimset('Algorithm','interior-point');
    %     lb = [-50 -15];
    %     ub = [50 15];
    %     [estimates(j,:) align_fun(j)] = fmincon(@(x) AP_3DaffineAlign(x,im_orig,im_align_z),starting,[],[],[],[],lb,ub,[],options);
    %     waitbar(j/stack_num,w,'Aligning stack to original');
    % end
    % close(w);
    %
    
end

% save the data
save([stack_path animal '_zalign_lumnorm_' num2str(skip_points)],'zalign');

%% Plot the results

plot_zalign = 1;
if plot_zalign == 1;

% Plot the data, find good times
curr_z = zalign.z_slice_smoothcurve;

day_sep = cellfun(@length,curr_z);
day_sep = cumsum(day_sep);

% z drift
z_median = cellfun(@median,curr_z);
z_drift = vertcat(curr_z{:});
z_align_median = nanmedian(z_median);
z_align_sem = 1.98*nanstd(z_median)/sqrt(sum(~isnan(z_median)));
z_good_indx = find(z_median < z_align_median + z_align_sem ... 
    & z_median > z_align_median - z_align_sem);
z_align = mean(z_median(z_good_indx));

z_smooth = cellfun(@(x) smooth(x,20), ...
    curr_z,'UniformOutput',false);
z_smooth_cat = vertcat(z_smooth{:});

z_cutoff_high = z_align+3;
z_cutoff_low = z_align-3;
good_frames = cellfun(@(x) x >= z_cutoff_low & ...
    x <= z_cutoff_high,z_smooth,'UniformOutput',false);

% give a little leway for bad frames: make surrounding frames bad
good_frames_smooth = cellfun(@(x) smooth(+x,5),good_frames,'UniformOutput',false);
good_frames = cellfun(@(x) x >= 0.9, good_frames_smooth,'UniformOutput',false);

% % Once you go bad, they're all bad? this isn't true if fixes...
%badframe_indx = cellfun(@(x) find(x == false,1),good_frames,'UniformOutput',false);

good_frames_cat = vertcat(good_frames{:});
z_smooth_cat_bad = z_smooth_cat;
z_smooth_cat_bad(good_frames_cat) = NaN;

% xy offset
offset_cat = vertcat(zalign.corr_offset_all{:});
offset = zeros(length(z_drift),1);
for i = 1:length(z_drift);
offset(i,1:2) = offset_cat{i}(round(z_drift(i)),:);
end
% max correlation
max_c_cat = vertcat(zalign.max_c_all{:});
max_c = [];
for i = 1:length(z_drift);
max_c = [max_c;max(max_c_cat{i},[],2)];
end

%figure('Name',animal);
figure;

z_plot = subplot(3,1,1);
hold on;
xlim([1 size(z_drift,1)])
fill([xlim fliplr(xlim)], ...
    [z_cutoff_low z_cutoff_low z_cutoff_high z_cutoff_high],'y',...
    'EdgeColor','none');
plot(z_drift,'.k');
plot(z_smooth_cat,'g','linewidth',2);
plot(z_smooth_cat_bad,'r','linewidth',2);
for i = 1:length(day_sep)
    line([day_sep(i) day_sep(i)],ylim,'linestyle','--','color','r');
    text(day_sep(i),33+3*mod(i,2),zalign.day{i},'HorizontalAlignment','Right')
end
title('Z alignment')
ylabel('z (\mum)')

offset_plot = subplot(3,1,2);
plot(offset)
for i = 1:length(day_sep)
    line([day_sep(i) day_sep(i)],ylim,'linestyle','--','color','r');
end
xlim([1 size(offset,1)])
title('X,Y alignment')
ylabel('x,y (~\mum)')

max_c_plot = subplot(3,1,3);
plot(max_c,'.k')
for i = 1:length(day_sep)
    line([day_sep(i) day_sep(i)],ylim,'linestyle','--','color','r');
end
xlim([1 size(max_c,1)])
title('Correlation')
ylabel('corr coef')

linkaxes([z_plot offset_plot max_c_plot],'x');


end