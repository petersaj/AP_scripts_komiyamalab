%% Load in image file tensor - no memory, tiffstack

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','off');
img_filename = [tiff_path tiff_filename];
clear im;
warning off
im = TIFFStack(img_filename);
warning on

%% load in image file

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','off');

img_filename = [tiff_path tiff_filename];

% get traces from current tiff file
imageinfo=imfinfo(img_filename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% try
%     matlabpool
% catch me
% end

im = zeros(N,M,numframes,'uint16');
disp('Loading file....')
for loadframe = 1:numframes
    im(:,:,loadframe) = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
    disp(['Frame ' num2str(loadframe) '/' num2str(numframes)]);
end
disp('Done')

%% load in multiple image files

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');

im_temp = [];
im = [];

% % get total number of frames
% ref_numframes = 0;
% for pfiles = 1:size(filelist,1);
%     ref_iminfo=imfinfo([tiff_path tiff_filename{pfiles}],'tiff');
%     ref_numframes_long=length(ref_iminfo);
%     ref_numframes=ref_numframes + ref_numframes_long/channels;
% end

for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    try
        matlabpool
    catch me
    end
    
    im_temp_long = zeros(N*M,numframes);
    im_temp = zeros(N,M);
    disp('Loading file....')
    
    for loadframe = 1:numframes
        im_temp = double(imread(img_filename,'tiff',loadframe));
        im_temp_long(:,loadframe) = im_temp(:);
    end
    im = [im im_temp_long];
    disp(['Done ' num2str(i) '/' num2str(length(tiff_filename))])
end
clearvars -except im

%% Create total STD image

clear all

try
    matlabpool
catch me
end

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
[avg_filename,avg_path]=uigetfile('*.tif','Choose Average TIFF file','Multiselect','off');

disp('Loading average file....')
im_avg = double(imread([avg_path avg_filename],'tiff'));
disp('Done')
imageinfo=imfinfo([avg_path avg_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
im_std = zeros(N,M);

for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    im = zeros(N,M,numframes);
    im_temp = zeros(N,M);
    disp(['Loading file ' num2str(i) '/' num2str(length(tiff_filename))  '....'])
    parfor loadframe = 1:numframes
        im_temp = double(imread(img_filename,'tiff',loadframe));
        im(:,:,loadframe) = im_temp;
    end
    disp('Done')
    
    im_diff = zeros(N,M);
    for j = 1:numframes;
        im_diff = im_diff + (im(:,:,j) - im_avg).^2;
    end
    
    im_std = im_std + sqrt(im_diff./numframes);
    
end

im_std = uint16(im_std./length(tiff_filename));

save_filename = [tiff_path tiff_filename{end}(1:15) 'STD.tif'];
for windows_lock = 1:100
    try
        imwrite(im_std,save_filename,'tif','Compression','none','WriteMode','overwrite');
        break;
    catch me
        pause(0.2);
    end
end



%% Get trace from only ROI by loading in whole movies and multithreading

[roi_file roi_path] = uigetfile('*','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
w = waitbar(0,'Getting and concatenating traces');
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % Create mask for selected ROI, get trace
    for n_polygon = 1:length(polygon.ROI);
        if ~isfield(polygon,'autosort');
            temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            cellMask(:,n_polygon) = temp_mask(:);
        elseif isfield(polygon,'autosort');
            temp_mask = polygon.autosort(:,:,n_polygon);
            cellMask(:,n_polygon) = temp_mask(:);
        end
    end
    
    im = [];
    im = zeros(N*M,numframes);
    
    try
        matlabpool
    catch me
    end
    
    tic
    disp('Loading file....')
    parfor loadframe = 1:numframes
        im_temp = double(imread(img_filename,'tiff',loadframe));
        % image: columns are frames
        im(:,loadframe) = im_temp(:);
    end
    disp('Done')
    toc
    
    roi_trace = zeros(length(polygon.ROI),numframes);
    w_curr = waitbar(0,'Getting activity from current file');
    traceVal = zeros(length(polygon.ROI),numframes);
    for n_polygon = 1:length(polygon.ROI)
        % sum of pixel brighness / number of pixels, avg fluor in ROI
        traceVal(n_polygon,:) = sum(im(cellMask(:,n_polygon),:)./sum(sum(cellMask(:,n_polygon))));
        waitbar(n_polygon/length(polygon.ROI),w_curr,'Getting activity from current file');
    end
    roi_trace = traceVal;
    close(w_curr);
    
    % create df/f for roi_trace
    roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        try
            norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two gaussians
        catch me
            try
                norm_fit_obj = gmdistribution.fit(cellTrace',1); % there is one gaussian
            catch me
                continue
            end
        end
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f
        cellTrace = cellTrace./baseline_firing - 1;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear roi_trace
    waitbar(i/length(tiff_filename),w,'Getting and concatenating traces');
end
close(w);

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);
% clear to save memory
im = [];


%% create df/f for roi_trace
roi_trace_df = zeros(size(roi_trace));
for i = 1:size(roi_trace,1)
    cellTrace = roi_trace(i,:);
    % don't bother if they're all NaNs
    if ~any(isfinite(cellTrace));
        continue
    end
    try
        norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
    catch me
        norm_fit_obj = gmdistribution.fit(cellTrace',1); % one hidden curve if error
    end
    dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
    x = min(cellTrace):dx:max(cellTrace);
    x = x(1:end-1); % it adds an extra, not sure why but can fix later
    norm_fit = pdf(norm_fit_obj,x');
    baseline_firing = x(norm_fit == max(norm_fit));
    
    % change to df/f make peakin relative to df
    cellTrace = (cellTrace-baseline_firing)./baseline_firing;
    roi_trace_df(i,:) = cellTrace;
end

%% Get distances, correlations by distance between rois

raw_corr_flag = 0;

zoom_factor = 2;

x_microns = (850/512)/zoom_factor;
y_microns = (850/128)/zoom_factor;

% get centers of all ROIs
roi_centers = zeros(size(roi_trace_long),2);
roi_polygon_centers = zeros(size(roi_trace_long),2);
for i = 1:size(roi_trace_long,1);
    % get center x and y from polygon boundaries
    center_x = mean(polygon.ROI{i}(:,1));
    center_y = mean(polygon.ROI{i}(:,2));
    roi_centers(i,1) = center_x*x_microns;
    roi_centers(i,2) = center_y*y_microns;
    roi_polygon_centers(i,1) = center_x;
    roi_polygon_centers(i,2) = -center_y;
end

if raw_corr_flag == 1;
    
    % loop through ROIs, sort by distance and get corrcoefs
    roi_corrcoef = zeros(size(roi_trace_long,1));
    roi_dists_all = [];
    roi_corrcoef_all = [];
    for j = 1:size(roi_trace_long,1);
        % for chosen roi, get distances and sort
        curr_roi = j;
        curr_center = roi_centers(curr_roi,:);
        center_diffs = roi_centers(:,1) - curr_center(1);
        center_diffs(:,2) = roi_centers(:,2) - curr_center(2);
        roi_dists = sqrt(center_diffs(:,1).^2 + center_diffs(:,2).^2);
        roi_dists(:,2) = [1:size(roi_dists,1)];
        roi_dists_sort = sortrows(roi_dists);
        
        roi_dists_all = [roi_dists_all; roi_dists(:,1)];
        
        roi_trace_sort = (roi_trace_long(roi_dists_sort(:,2),:));
        for i = 1:size(roi_trace_sort,1)
            curr_corrcoef = corrcoef(roi_trace_sort(1,:),roi_trace_sort(i,:));
            curr_roi_corrcoef(i) = curr_corrcoef(2);
        end
        
        roi_corrcoef_all = [roi_corrcoef_all; curr_roi_corrcoef'];
        
        roi_corrcoef(curr_roi,:) = curr_roi_corrcoef;
        disp(j/size(roi_trace_long,1));
    end
    
    num_bins = 100;
    [n,bin] = histc(roi_dists_all,linspace(min(roi_dists_all),max(roi_dists_all),num_bins));
    for i = 1:num_bins
        flagBinMembers = (bin == i);
        binMembers = roi_corrcoef_all(flagBinMembers);
        binMean(i) = mean(binMembers);
    end
    figure;plot(binMean,'k','linewidth',2)
    set(gca,'XTick',[1:10:num_bins])
    set(gca,'XTickLabel',round(linspace(min(roi_dists_all),max(roi_dists_all),num_bins/10)));
    xlabel('Distance (\mum)')
    ylabel('Correlation coefficient')
    
end

% show distances vs. correlation but using coincident spike index instead
% look by surround
roi_dists_all = [];
coincidence_indx_all = [];
coincidence_indx_tk_all = [];
cell_pair_all = [];
num_coincidence_all = [];

for curr_roi = 1:size(roi_trace_long,1);
    
    % for chosen roi, get distances
    curr_center = roi_centers(curr_roi,:);
    center_diffs = roi_centers(:,1) - curr_center(1);
    center_diffs(:,2) = roi_centers(:,2) - curr_center(2);
    roi_dists = sqrt(center_diffs(:,1).^2 + center_diffs(:,2).^2);
    roi_dists_all = [roi_dists_all; roi_dists(:,1)];
    
    frame_back = 5;
    frame_forward = 5;
    roi_trace_coincident = zeros(size(roi_trace_long));
    for i = 1:length(spike_frames{curr_roi})
        frame_search = spike_frames{curr_roi}(i)-frame_back:spike_frames{curr_roi}(i)+frame_forward;
        for j = 1:size(roi_trace_long,1)
            event_detect = [];
            event_detect = ismember(spike_frames{j},frame_search);
            if any(event_detect)
                roi_trace_coincident(j,spike_frames{curr_roi}(i)) = 1;
            end
        end
    end
    
    % get coincidence index for each ROI (from Komiyama 2010)
    num_frames = size(roi_trace_long,2);
    num_frames_search = size(roi_trace_long,2)/(frame_back+frame_forward+1);
    coincidence_indx = zeros(size(roi_trace_long,1),1);
    coincidence_indx_tk = zeros(size(roi_trace_long,1),1);
    spike_1 = length(spike_frames{curr_roi});
    spike_1_prob = spike_1/num_frames;
    for i = 1:size(roi_trace_long,1)
        spike_2 = length(spike_frames{i});
        spike_2_prob = spike_2/num_frames;
        % this is what they do in the paper... does this make sense?
        geometric_mean = num_frames*sqrt(spike_1_prob*spike_2_prob);
        coincidence_indx_tk(i) = ((spike_1+spike_2)-(num_frames*spike_1_prob*spike_2_prob))/(geometric_mean);
        % this makes more intuitive sense to me, co-spikes/all spikes
        co_spike = sum(roi_trace_coincident(i,:));
        coincidence_indx(i) = (co_spike*2)/(spike_1+spike_2);
        cell_pair_all = [cell_pair_all; curr_roi i];
        num_coincidence_all = [num_coincidence_all co_spike];
    end
    coincidence_indx_tk_all = [coincidence_indx_tk_all; coincidence_indx_tk];
    coincidence_indx_all = [coincidence_indx_all;coincidence_indx];
    disp(num2str(curr_roi/size(roi_trace_long,1)))
end

% get rid of cells with themselves (distance = 0)
self_compare = find(roi_dists_all == 0);
roi_dists_all(self_compare)= [];
coincidence_indx_all(self_compare) = [];
cell_pair_all(self_compare,:) = [];
num_coincidence_all(self_compare) = [];

num_bins = 40;
[n,bin] = histc(roi_dists_all,linspace(min(roi_dists_all),max(roi_dists_all),num_bins));
binMean = [];
for i = 1:num_bins
    flagBinMembers = (bin == i);
    binMembers = coincidence_indx_all(flagBinMembers);
    binMean(i) = nanmean(binMembers);
end
binMean_notnan = find(~isnan(binMean) == 1);
figure;plot(binMean_notnan,binMean(binMean_notnan),'k','linewidth',2)
set(gca,'XTick',[1:2:num_bins])
set(gca,'XTickLabel',round(linspace(min(roi_dists_all),max(roi_dists_all),num_bins/2)));
xlabel('Distance (\mum)')
ylabel('Correlation coefficient')

% plot the ROIs
figure
hold on;
set(gca,'color','k');
set(gcf,'color','k');
for i = 1:length(polygon.ROI)
    r = patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
    set(r,'linewidth',length(spike_frames{i})/10+0.01);
    set(r,'edgecolor','w');
end

for i = 1:length(mod_cells_early)
    curr_color = mod_cells_early_percent(i)*0.5+0.5;
    r = patch(polygon.ROI{mod_cells_early(i)}(:,1),-polygon.ROI{mod_cells_early(i)}(:,2),[0 curr_color 0]);
    set(r,'linewidth',length(spike_frames{mod_cells_early(i)})/10+0.01);
    set(r,'edgecolor','w');
end

for i = 1:length(mod_cells_middle)
    curr_color = mod_cells_middle_percent(i)*0.5+0.5;
    r = patch(polygon.ROI{mod_cells_middle(i)}(:,1),-polygon.ROI{mod_cells_middle(i)}(:,2),[curr_color curr_color 0.5]);
    set(r,'linewidth',length(spike_frames{mod_cells_middle(i)})/10+0.01);
    set(r,'edgecolor','w');
end

for i = 1:length(mod_cells_late)
    curr_color = mod_cells_late_percent(i)*0.5+0.5
    r = patch(polygon.ROI{mod_cells_late(i)}(:,1),-polygon.ROI{mod_cells_late(i)}(:,2),[curr_color 0.5 0.5]);
    set(r,'linewidth',length(spike_frames{mod_cells_late(i)})/10+0.01);
    set(r,'edgecolor','w');
end

ylim([-128 0]);
xlim([0 512]);

% plot lines on ROIs - thicker = more correlated
for i = 1:length(coincidence_indx_all)
    cell1 = cell_pair_all(i,1);
    cell2 = cell_pair_all(i,2);
    if coincidence_indx_all(i) > 0
        %         line([roi_polygon_centers(cell1,1) roi_polygon_centers(cell2,1)], ...
        %             [roi_polygon_centers(cell1,2) roi_polygon_centers(cell2,2)],'color',[1-coincidence_indx_all(i) 1-coincidence_indx_all(i) 1-coincidence_indx_all(i)]);
        x = [roi_polygon_centers(cell1,1) roi_polygon_centers(cell2,1)];
        y = [roi_polygon_centers(cell1,2) roi_polygon_centers(cell2,2)];
        z = [0 0];
        r = patch(x,y,z);
        set(r,'edgealpha',coincidence_indx_all(i));
        set(r,'linewidth',num_coincidence_all(i)/10);
        set(r,'edgecolor','w');
    end
end

%% Compare peaks of traces

minpeak = 0.1; % min increase (std) to be considered active
roi_peak = {};
for i = 1:size(roi_trace_df,1)
    [max_ca min_ca] = peakdet(roi_trace_df(i,:),minpeak);
    roi_peak{i} = max_ca;
end

%% get calcium transients (the dumb way)

%peakstd = [1.5*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
peakstd = [1.5*ones(1,10)]; % definitions for max trigger, in stds
peak_back = 10;
local_max = 30;
drop_percentage = .8;
spike_frames = {};
movie_frames = 4000;

num_rois = size(roi_trace_long,1);

% this is to catch a dumb bug: if ROI is offscreen, just make noise
for i = 1:num_rois
    if all(roi_trace_long(i,:) == roi_trace_long(i,1))
        roi_trace_long(i,:) = rand(size(roi_trace_long(i,:)));
    end
end

% Estimate noise through smoothing
size_smooth = 90;
smooth_type = 'loess';
size_file = 4000;
roi_trace_df_smooth = zeros(size(roi_trace_df));
for i = 1:size(roi_trace_df,1);
    for j = 1:size(roi_trace_df,2)/size_file;
        roi_trace_df_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_df(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end
smooth_std = zeros(size(roi_trace_df,1),1);
smooth_std = sqrt(sum((abs(roi_trace_df - roi_trace_df_smooth).^2)/size(roi_trace_df_smooth,2),2));

parfor i = 1:num_rois
    % Estimate noise through smoothing
    curr_trace_std = smooth_std(i);
    minpeak = curr_trace_std*peakstd;
    minpeak = 0.5
    spike_frames{i} = ap_schmittDetect(roi_trace_df(i,:),minpeak,peak_back,local_max,drop_percentage);
    disp([num2str(i) '/' num2str(num_rois)])
end

% for i = 1:num_rois
%     % eliminate peaks at file junctions because they're mistimed
%     if ~isempty(spike_frames{i})
%         file_stitch = [];
%         file_stitch = find(mod(spike_frames{i}-1,movie_frames) == 0);
%         peak_frames(file_stitch,:) = [];
%     end
% end

%%%% what's been kind of working
%[peak_frames] = ap_schmittDetect(roi_trace_long(24,:),minpeak,5,20,0.5);
%%%%

%% Concatenatenate traces of .roi already made

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','on');
for i = 1:length(tiff_filename);
    
    load([tiff_path tiff_filename],'-MAT')
    % create df/f for roi_trace
    roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f make peakin relative to df
        cellTrace = cellTrace./baseline_firing;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear polygon roi_mask roi_trace
end

num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_acq = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_acq);

%% Concatenatenate traces from only ROI

auto_align_flag = 0; % do you want to auto align ROI for each file?
if auto_align_flag == 1;
    [align_file align_path] = uigetfile('*.tif','Pick image to align to');
    img_align = double(imread([align_path align_file],'tiff'));
end

[roi_file roi_path] = uigetfile('*','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
w = waitbar(0,'Getting and concatenating traces');
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    if auto_align_flag == 0;
        % Create mask for selected ROI, get trace
        for n_polygon = 1:length(polygon.ROI);
            if ~isfield(polygon,'autosort');
                cellMask(:,:,n_polygon) = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            elseif isfield(polygon,'autosort');
                cellMask(:,:,n_polygon) = polygon.autosort(:,:,n_polygon);
            end
        end
        
        % if auto align selected, make masks from aligned ROIs
    elseif auto_align_flag == 1;
        aligned_ROI = cell(size(polygon.ROI));
        % get current file average
        img_avg = zeros(N,M);
        w_avg = waitbar(0,'Creating average of current file');
        for frame = 1:numframes;
            img_temp(:,:) = double(imread(img_filename,'tiff',frame));
            img_avg(:,:) = img_avg + img_temp./numframes;
            waitbar(frame/numframes,w_avg,'Creating average of current file');
        end
        close(w_avg);
        
        % DON'T HARDCODE THIS IN THE FUTURE
        disp 'xy_ratio hardcoded as 4'
        xy_ratio = 4;
        
        % get estimate for dx dy
        cc = normxcorr2(img_align,img_avg);
        [max_cc, imax] = max(abs(cc(:)));
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        
        % maximize correlation through transformation
        starting = [xpeak ypeak 0]; % dx dy theta
        options = optimset('Algorithm','interior-point');
        lb = [-100 -100/xy_ratio -pi/2];
        ub = [100 100/xy_ratio pi/2];
        % estimate transform of orig -> curr
        display 'Finding best fit...'
        estimates = fmincon(@(x) AP_affineAlign(x,img_align,img_avg,xy_ratio),starting,[],[],[],[],lb,ub,[],options);
        display 'Done'
        % get estimates from best guess
        dx = estimates(1);
        dy = estimates(2);
        theta = estimates(3);
        
        % transform polygons according to image match
        for poly_align = 1:length(polygon.ROI)
            
            udata = [1 N] - (N+1)/2;% input image x
            vdata = [1 M] - (M+1)/2;% input image y
            
            xdata = udata;% output image x
            ydata = vdata;% output image y
            
            % define affine transform matrix
            A = ...
                [cos(theta) sin(theta) 0;
                -sin(theta) cos(theta) 0;
                dx          dy          1];
            
            tform_roi = maketform('affine',A);
            
            % transform polygon, fix for xy_ratio
            aligned_ROI{poly_align} = tformfwd(tform_roi,[polygon.ROI{poly_align}(:,1)...
                polygon.ROI{poly_align}(:,2).*xy_ratio]);%, ...
            % 'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
            % correct final ROIs for xy_ratio
            aligned_ROI{poly_align}(:,2) = aligned_ROI{poly_align}(:,2)./xy_ratio;
            
            cellMask(:,:,poly_align) = poly2mask(aligned_ROI{poly_align}(:,1)',aligned_ROI{poly_align}(:,2)',N,M);
        end
    end
    
    roi_trace = zeros(length(polygon.ROI),numframes);
    w_curr = waitbar(0,'Getting activity from current file');
    for frame=1:numframes
        im_s=double(imread([tiff_path tiff_filename{i}],frame));
        traceVal = [];
        for n_polygon = 1:length(polygon.ROI)
            normalactivityMap = [];
            % sum of pixel brighness / number of pixels, avg fluor in ROI
            normactivityMap = (im_s.*cellMask(:,:,n_polygon))./sum(sum(cellMask(:,:,n_polygon)));
            traceVal(n_polygon,1) = sum(normactivityMap(:));
        end
        roi_trace(:,frame) = traceVal;
        waitbar(frame/numframes,w_curr,'Getting activity from current file');
    end
    close(w_curr);
    
    % create df/f for roi_trace
    roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        try
            norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two gaussians
        catch me
            try
                norm_fit_obj = gmdistribution.fit(cellTrace',1); % there is one gaussian
            catch me
                continue
            end
        end
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f
        cellTrace = cellTrace./baseline_firing - 1;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear roi_trace
    waitbar(i/length(tiff_filename),w,'Getting and concatenating traces');
    
end
close(w);

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);

%% Get and concatenate events
n_channels = 1; % true at the moment for all motion corrected files
n_frames = 1000; %get this later

[bhv_filename bhv_path] = uigetfile('*','Behavior file');
bhv = load([bhv_path bhv_filename],'-MAT');

saved_history = [];
saved_history = bhv.saved_history; % this in for now, somewhere below fix

[acq_filename_all acq_path] = uigetfile('*.xsg','Choose XSG files','Multiselect','on');

w = waitbar(0, 'Getting and concatenating events');

left_down_all = [];
left_up_all = [];
reward_all = [];
unrewarded_all = [];
cue_onset_all = [];
frames_left_down_all = [];
frames_left_up_all = [];
frames_reward_all = [];
frames_unrewarded_all = [];
frames_cue_onset_all = [];
timeDiff = [];
reward_trial = [];
unrewarded_trial = [];

for i = 1:length(acq_filename_all);
    
    % get current acq filename
    acq_filename = [acq_path acq_filename_all{i}];
    
    % this is hard coded until the header problem is fixed
    % % %     % pull out time for each for each frame from header
    % % %     img_info = imfinfo([img_pathname img_filename]);
    % % %     msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
    % % %     msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
    % % %     msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
    % % %     lines = img_info(1).Height;
    % % %     fpms = msPerLine*lines;
    
    msPerLine = 1.24;
    lines = 128;
    fpms = msPerLine*lines;
    
    % Get trials in s since started
    TrialNumList = TK_ReadBitCode([acq_filename]);
    close(gcf)
    % put trial time in ms
    TrialNumList(2,:) = TrialNumList(2,:)*1000;
    
    left_down = [];
    left_up = [];
    reward = [];
    unrewarded = [];
    cue_onset = [];
    
    for trial = TrialNumList(1,:)
        if trial > length(bhv.saved_history.RewardsSection_LastTrialEvents);
            continue
        end
        bhv_current = bhv.saved_history.RewardsSection_LastTrialEvents{trial,1};
        % fixes possible bug of putting in old timepoints at the end
        LastGoodTrial = find(diff(bhv_current(:,3)) < 0,1);
        if ~isempty(LastGoodTrial)
            bhv_current = bhv_current(1:LastGoodTrial,:);
        end
        % get time lag of behavior to acq (t_bhv + timeDiff = t_acq) in ms
        trialStart_bhv = bhv_current(find(bhv_current(:,1) == 101 & bhv_current(:,2) == 0,1,'last'),3)*1000;
        lag_time = TrialNumList(2,find(TrialNumList(1,:) == trial)) - trialStart_bhv;
        timeDiff(trial,1) = lag_time;
        % make a second column for which acquisition you're in to frame shift
        timeDiff(trial,2) = i*ones(length(lag_time),1);
        % time for all lever up and down (left down = 5, left up = 6) in ms
        %left_down = [left_down; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff(trial)];
        %left_up = [left_up; bhv_current(find(bhv_current(:,2) == 6),3)*1000 + timeDiff(trial)];
        reward = [reward; bhv_current(find(bhv_current(:,1) == 44,1),3)*1000 + timeDiff(trial,1)];
        unrewarded = [unrewarded; bhv_current(find(bhv_current(:,1) == 45,1),3)*1000 + timeDiff(trial,1)];
        if ~isempty(find(bhv_current(:,1) == 44,1))
            reward_trial = [reward_trial trial];
        end
        if ~isempty(find(bhv_current(:,1) == 45,1))
            unrewarded_trial = [unrewarded_trial trial];
        end
        cue_onset = [cue_onset; bhv_current(find(bhv_current(:,1) == 41,1),3)*1000 + timeDiff(trial,1)];
    end
    
    % concatenate all events, multiply by img number to get total time
    
    %left_down_all = [left_down_all; left_down];
    %left_up_all = [left_up_all; left_up];
    unrewarded_all = [unrewarded_all; unrewarded];
    reward_all = [reward_all; reward];
    cue_onset_all = [cue_onset_all; cue_onset];
    
    %frames_left_down = ceil(left_down(~isnan(left_down))/fpms)';
    %frames_left_up = ceil(left_up(~isnan(left_up))/fpms)';
    frames_reward = ceil(reward(~isnan(reward))/fpms)';
    frames_unrewarded = ceil(unrewarded(~isnan(unrewarded))/fpms)';
    frames_cue_onset = ceil(cue_onset(~isnan(cue_onset))/fpms)';
    
    % round up for number of channels
    %frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
    %frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
    frames_cue_onset = (1-mod(frames_cue_onset,n_channels))+frames_cue_onset;
    frames_reward = (1-mod(frames_reward,n_channels))+frames_reward;
    frame_unrewarded = (1-mod(frames_unrewarded,n_channels))+frames_unrewarded;
    
    % Get frame for each event
    %frames_left_down_all = [frames_left_down_all frames_left_down+n_frames*(i-1)];
    %frames_left_up_all = [frames_left_up_all frames_left_up+n_frames*(i-1)];
    frames_reward_all = [frames_reward_all frames_reward+n_frames*(i-1)];
    frames_unrewarded_all = [frames_unrewarded_all frames_unrewarded+n_frames*(i-1)];
    frames_cue_onset_all = [frames_cue_onset_all frames_cue_onset+n_frames*(i-1)];
    
    waitbar(i/length(acq_filename_all),w,'Getting and concatenating events');
    
    
end

% get all lever timing at once to get rid of artifacts correctly
AP_getLeftLeverTimesImg
% correct lever times for time delay, keep old ones
left_down_acq = left_down.*1000 + timeDiff(trial_down,1);
trial_down_acq = trial_down;
states_down_acq = states_down;
left_up_acq = left_up.*1000 + timeDiff(trial_up,1);
trial_up_acq = trial_up;
states_up_acq = states_up;
% get rid of lever times when trial not recorded
trial_notRecorded = find(timeDiff(:,1) == 0);
elim_left_down = ismember(trial_down,trial_notRecorded);
elim_left_up = ismember(trial_up,trial_notRecorded);

left_down_acq(elim_left_down) = [];
trial_down_acq(elim_left_down) = [];
states_down_acq(elim_left_down) = [];
left_up_acq(elim_left_up) = [];
trial_up_acq(elim_left_up) = [];
states_up_acq(elim_left_up) = [];
% pick out the left_down that were in the right place (state 41)
left_down_air = left_down_acq(states_down_acq == 41);

% convert lever times to frames
frames_left_down = ceil(left_down_acq(~isnan(left_down_acq))/fpms)';
frames_left_up = ceil(left_up_acq(~isnan(left_up_acq))/fpms)';
% round frames based on channels
frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
% account for concatenation
% get amount of frames to shift based on trial number
frame_shift_left_down = [];
frame_shift_left_up = [];
frame_shift_left_down = (timeDiff(trial_down_acq,2)-1)*n_frames;
frame_shift_left_up = (timeDiff(trial_up_acq,2)-1)*n_frames;
frames_left_down_all = [frames_left_down_all frames_left_down'+frame_shift_left_down];
frames_left_up_all = [frames_left_up_all frames_left_up'+frame_shift_left_up];

% get the frames where left lever down during correct time
frames_left_down_air_all = frames_left_down_all(states_down_acq == 41);

close(w)
cd ..

%% Make traces for plots, etc

%mod_cells = mod_up_cells;
mod_cells = 12;

% hist_flag = 1 if you want to see histograms and raw trace with events
hist_flag = 0;
% reward_trace_flag = 1 if you want to see all traces plotted for each rew
reward_trace_flag = 0;

frame_sample = 60;
reward_trace_concat = zeros(length(frames_reward_all),frame_sample*2+1,length(mod_cells));

%for curr_mod = 1:size(roi_trace_long,1)
for curr_mod = 1:length(mod_cells);
    
    curr_roi = mod_cells(curr_mod);
    %curr_roi = curr_mod;
    
    % get traces time-locked to events
    reward_trace = zeros(length(frames_reward_all),frame_sample*2+1);
    for reward_num = 1:length(frames_reward_all);
        if frames_reward_all(reward_num)-frame_sample > 1 && ...
                frames_reward_all(reward_num)+frame_sample < size(roi_trace_long,2);
            reward_trace(reward_num,:) = roi_trace_long...
                (curr_roi,frames_reward_all(reward_num)-frame_sample:frames_reward_all(reward_num)+frame_sample);
        else
            reward_trace(reward_num,:) = NaN;
        end
    end
    
    %     left_down_trace = [];
    %     for left_down_num = 1:length(frames_left_down_all);
    %         if frames_left_down_all(left_down_num)-frame_sample > 1 && ...
    %                 frames_left_down_all(left_down_num)+frame_sample < size(roi_trace_long,2);
    %             left_down_trace(left_down_num,:) = roi_trace_long...
    %                 (curr_roi,frames_left_down_all(left_down_num)-frame_sample:frames_left_down_all(left_down_num)+frame_sample);
    %         end
    %     end
    %
    %     left_up_trace = [];
    %     for left_up_num = 1:length(frames_left_up_all);
    %         if frames_left_up_all(left_up_num)-frame_sample > 1 && ...
    %                 frames_left_up_all(left_up_num)+frame_sample < size(roi_trace_long,2);
    %             left_up_trace(left_up_num,:) = roi_trace_long...
    %                 (curr_roi,frames_left_up_all(left_up_num)-frame_sample:frames_left_up_all(left_up_num)+frame_sample);
    %         end
    %     end
    
    %     % get the left downs and up that count (down during 41, up during 43)
    %     left_down_cued = find(states_down_acq == 41);
    %     left_up_cued = find(states_up_acq == 43);
    %
    %     left_down_cued_trace = [];
    %     for left_down_cued_num = left_down_cued';
    %         if frames_left_down_all(left_down_cued_num)-frame_sample > 1 && ...
    %                 frames_left_down_all(left_down_cued_num)+frame_sample < size(roi_trace_long,2);
    %             left_down_cued_trace(find(left_down_cued_num == left_down_cued),:) = roi_trace_long...
    %                 (curr_roi,frames_left_down_all(left_down_cued_num)-frame_sample:frames_left_down_all(left_down_cued_num)+frame_sample);
    %         end
    %     end
    %
    %     left_up_cued_trace = [];
    %     for left_up_cued_num = left_up_cued';
    %         if frames_left_up_all(left_up_cued_num)-frame_sample > 1 && ...
    %                 frames_left_up_all(left_up_cued_num)+frame_sample < size(roi_trace_long,2);
    %             left_up_cued_trace(find(left_up_cued_num == left_up_cued),:) = roi_trace_long...
    %                 (curr_roi,frames_left_up_all(left_up_cued_num)-frame_sample:frames_left_up_all(left_up_cued_num)+frame_sample);
    %         end
    %     end
    
    
    %     % get time differences for events around peaks
    %     %peak_frames = roi_peak_max{curr_roi}(:,1);
    %     peak_frames = spike_frames{curr_roi};
    %     frame_search = 5;
    %     left_down_peak_diff = [];
    %     left_down_peak_event = [];
    %     left_up_peak_diff = [];
    %     left_up_peak_event = [];
    %     reward_peak_diff = [];
    %     reward_peak_event = [];
    %
    %     for curr_peak = 1:length(peak_frames)
    %         left_down_peak = [];
    %         left_up_peak = [];
    %         reward_peak = [];
    %
    %         frame_range = peak_frames(curr_peak)-frame_search:peak_frames(curr_peak);
    %
    %         left_down_peak = ismember(frames_left_down_all,frame_range);
    %         left_down_peak_diff = [left_down_peak_diff; frames_left_down_all(left_down_peak)-peak_frames(curr_peak)];
    %         if sum(left_down_peak)~=0
    %             left_down_peak_event(curr_peak) = 1;
    %         end
    %         left_up_peak = ismember(frames_left_up_all,frame_range);
    %         left_up_peak_diff = [left_up_peak_diff; frames_left_up_all(left_up_peak)-peak_frames(curr_peak)];
    %         if sum(left_up_peak)~=0
    %             left_up_peak_event(curr_peak) = 1;
    %         end
    %         reward_peak = ismember(frames_reward_all,frame_range);
    %         reward_peak = find(reward_peak == 1,1);
    %         reward_peak_diff = [reward_peak_diff; frames_reward_all(reward_peak)-peak_frames(curr_peak)];
    %         if sum(reward_peak)~=0
    %             reward_peak_event(curr_peak) = 1;
    %         end
    %
    %     end
    %
    %     % average all peaks
    %     peak_trace = [];
    %     for peak_num = 1:length(peak_frames);
    %         % only take peaks with frame_sample = surrounding frames
    %         if peak_frames(peak_num)-frame_sample > 1 && ...
    %                 peak_frames(peak_num)+frame_sample < size(roi_trace_long,2)
    %             peak_trace(peak_num,:) = roi_trace_long...
    %                 (curr_roi,peak_frames(peak_num)-frame_sample:peak_frames(peak_num)+frame_sample);
    %         end
    %     end
    
    if hist_flag == 1;
        % create stacked histogram for peak events
        
        % get the max peak for bin, round to nearest tenth
        %peak_max_max = ceil(max(roi_trace_long(curr_roi,peak_frames))*10)/10;
        
        %         [hist_all_peak_n hist_all_peak_x] = hist(spike_amplitudes{curr_roi},[1:0.1:peak_max_max]);
        %         [hist_left_down_peak_n hist_left_down_peak_x] = ...
        %             hist(spike_amplitudes{curr_roi}((left_down_peak_event == 1)),[1:0.1:peak_max_max]);
        %         [hist_left_up_peak_n hist_left_up_peak_x] = ...
        %             hist(spike_amplitudes{curr_roi}((left_up_peak_event == 1)),[1:0.1:peak_max_max]);
        %         [hist_reward_peak_n hist_reward_peak_x] = ...
        %             hist(spike_amplitudes{curr_roi}((reward_peak_event == 1)),[1:0.1:peak_max_max]);
        %
        %         figure;
        %         bar(hist_all_peak_x,[hist_all_peak_n; hist_reward_peak_n; hist_left_up_peak_n; hist_left_down_peak_n]');
        %         colormap(summer)
        %         legend('All','Reward','Left up','Left down');
        %         title(['ROI ' num2str(curr_roi)]);
        
        % plot raw trace with ticks: down = black, up = red, reward = blue
        figure;hold on;
        plot(roi_trace_long(curr_roi,:))
        for i = 1:length(frames_left_down_all)
            line([frames_left_down_all(i) frames_left_down_all(i)], [0 0.2],'color','black')
        end
        for i = 1:length(frames_left_up_all)
            line([frames_left_up_all(i) frames_left_up_all(i)], [0.3 0.5],'color','red')
        end
        for i = 1:length(frames_reward_all)
            line([frames_reward_all(i) frames_reward_all(i)], [0.6 0.8],'color','blue')
        end
        plot(spike_frames{curr_roi},roi_trace_long(curr_roi,spike_frames{curr_roi}),'.r');
    end
    
    if reward_trace_flag == 1;
        figure;
        imagesc(reward_trace);
        colormap(gray);
        caxis([0 2]);
        title(num2str(curr_roi));
        x_sec = round(((str2num(get(gca,'XTickLabel')) - frame_sample+1).*(1.24*128/1000)).*10)./10;
        set(gca,'XTickLabel',x_sec);
        xlabel('Time from reward (s)')
        ylabel('Reward number')
    end
    
    reward_trace_concat(:,:,curr_mod) = reward_trace;
    %     % only count trials with events around reward, NaN otherwise
    %     rewarded_trials = surround_reward(curr_roi,:) > 0;
    %     reward_trace_concat(rewarded_trials,:,curr_mod) = NaN;
    
end

% plot reward trace response
sqr_num = ceil(sqrt(size(reward_trace_concat,3)));
figure;
for i = 1:size(reward_trace_concat,3)
    subplot(sqr_num,sqr_num,i);
    imagesc(reward_trace_concat(:,:,i));
    colormap(gray);
    title(num2str(mod_cells(i)));
    line([frame_sample frame_sample],ylim,'linestyle','--','color','k');
    axis off
end

% get ylimit
ymax = max(max(nanmean(reward_trace_concat,1))+0.1);
ymin = min(min(nanmean(reward_trace_concat,1)));

% plot average response
sqr_num = ceil(sqrt(size(reward_trace_concat,3)));
figure;
for i = 1:size(reward_trace_concat,3)
    subplot(sqr_num,sqr_num,i);
    plot(nanmean(reward_trace_concat(:,:,i)),'linewidth',2,'color','k')
    hold on
    nanstd_reward_trace_u = nanmean(reward_trace_concat(:,:,i))+...
        nanstd(reward_trace_concat(:,:,i));
    nanstd_reward_trace_l = nanmean(reward_trace_concat(:,:,i))-...
        nanstd(reward_trace_concat(:,:,i));
    nansem_reward_trace_u = nanmean(reward_trace_concat(:,:,i)) +...
        1.96*((nanstd(reward_trace_concat(:,:,i))/sqrt...
        (size(reward_trace_concat(:,:,i),1))));
    nansem_reward_trace_l = nanmean(reward_trace_concat(:,:,i)) -...
        1.96*((nanstd(reward_trace_concat(:,:,i))/sqrt...
        (size(reward_trace_concat(:,:,i),1))));
    x_reward = [1:1:size(reward_trace_concat(:,:,i),2)];
    
    %     % get bootstrap confidence intervals - this takes longer and I think
    %     it does the same thing
    
    %     num_rep = 1000;
    %     ci = zeros(2,size(reward_trace_concat,2));
    %     for j = 1:size(reward_trace_concat,2)
    %         ci(:,j) = bootci(num_rep,@mean,reward_trace_concat(:,j,i));
    %         disp(['Bootstrapping ' num2str(j/size(reward_trace_concat,2))]);
    %     end
    %     jbfill(x_reward,ci(2,:),ci(1,:),'blue','none',0,0.5);
    
    
    jbfill(x_reward,nansem_reward_trace_u,...
        nansem_reward_trace_l,'blue','none',0,0.5);
    set(gca,'xlim',[0 2*frame_sample+1]);
    x_sec = round(((str2num(get(gca,'XTickLabel')) - frame_sample-1).*(1.24*128/1000)).*10)./10;
    set(gca,'XTickLabel',x_sec);
    ylabel('{\Delta}F/F_0')
    title(num2str(mod_cells(i)));
    ylim([ymin ymax]);
    line([frame_sample frame_sample],ylim,'linestyle','--','color','k');
    axis off
    ylim_curr = ylim;
    %     % plot significantly above 0 points with astericks
    %     sig_increase_x = find(nansem_reward_trace_l > 0);
    %     sig_increase_y = ones(size(sig_increase_x))*(ylim_curr(2)-0.01*diff(ylim_curr));
    %     text(sig_increase_x, sig_increase_y, '*', 'FontSize', 14, 'FontWeight','bold');
    
end

% plot scale bar (assume last box isn't filled)
subplot(sqr_num,sqr_num,sqr_num^2);
line([0 12.6008],[0 0],'linewidth',3,'color','k') % 2 second in frames
line([0 0],[0 0.5],'linewidth',3,'color','k') % 50% DF/F0
ylim([ymin ymax])
xlim([-size(reward_trace_concat,2)/2+1 size(reward_trace_concat,2)/2+1])
axis off

%% Plot events with all peaks

figure;hold on;
for i = 1:size(roi_trace_long,1)
    for j = 1:size(spike_frames{i})
        % color lines by amplitude - darker is larger, -0.2 to make all vis
        if spike_amplitudes{i}(j) > 0; %FIX THIS! spike amplitudes should never be neg, but are sometimes
            color_amp = 1-(spike_amplitudes{i}(j)/(max(spike_amplitudes{i})+0.2));
            line([spike_frames{i}(j) spike_frames{i}(j)], ...
                [i i+0.9],'color',[color_amp color_amp color_amp]);
        end
    end
end

topline = length(spike_frames);
for i = 1:length(frames_left_down_all)
    line([frames_left_down_all(i) frames_left_down_all(i)], [topline+1 topline+1.9],'color','green')
end
for i = 1:length(frames_left_up_all)
    line([frames_left_up_all(i) frames_left_up_all(i)], [topline+2 topline+2.9],'color','red')
end
for i = 1:length(frames_reward_all)
    line([frames_reward_all(i) frames_reward_all(i)], [topline+3 topline+3.9],'color','blue')
end
for i = 1:length(frames_cue_onset_all)
    line([frames_cue_onset_all(i) frames_cue_onset_all(i)], [topline+4 topline+4.9],'color','cyan')
end
ylim([0 topline+5]);

%% Plot stack of traces
figure;
hold on;
mod_up_cells = [1:size(roi_trace_df,1)]';
for i = 1:size(mod_up_cells,1)
    %for i = 1:length(temp)
    plot(roi_trace_df(mod_up_cells(i),:) + i,'-k','linewidth',2);
%     plot(roi_trace_smooth_df(mod_up_cells(i),:) + i,'-r');
    if ~isempty(spike_frames{mod_up_cells(i)});
        plot(spike_frames{mod_up_cells(i)}, i,'.b','MarkerSize',20);
    end
    %      plot(spike_frames{mod_up_cells(i)}, spike_amplitudes{mod_up_cells(i)} + i,'.r');
    %      plot(spike_t0s{mod_up_cells(i)},zeros(1,length(spike_t0s{mod_up_cells(i)}))+i,'.r');
end

% % plot lever down as dotted black lines
% for i = 1:length(frames_left_down_all)
%     line([frames_left_down_all(i) frames_left_down_all(i)],ylim,'linestyle','--','color','k');
% end
%
% % plot lever up as dotted blue lines
% for i = 1:length(frames_left_up_all)
%     line([frames_left_up_all(i) frames_left_up_all(i)],ylim,'linestyle','--','color','b');
% end

% % plot rewards as dashed red lines
% for i = 1:length(frames_reward_all)
%     line([frames_reward_all(i) frames_reward_all(i)],ylim,'linestyle','--','color','r');
% end
% 
% % plot non-rewards as dashed magenta lines
% for i = 1:length(frames_unrewarded_all)
%     line([frames_unrewarded_all(i) frames_unrewarded_all(i)],ylim,'linestyle','--','color','m');
% end

% % plot air onset as dotted blue lines
% for i = 1:length(frames_cue_onset_all)
%     line([frames_cue_onset_all(i) frames_cue_onset_all(i)],ylim,'linestyle',':','color','b');
% end

%% Show which ROIs had events around given timeframe
frame_back = 40;
frame_forward = 40;
% events around reward/left up
surround_reward = zeros(size(roi_trace_df,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(roi_trace_df(curr_roi,spike_frames{curr_roi}(event_detect)));
            surround_reward(curr_roi,curr_reward) = event_detect_amp;
        end
    end
end
figure;imagesc(surround_reward);
colormap(gray);
title('Events around reward/left up')

%% Find reward-related cells (requires last cell be run)
tic
n_sample = 10000;

frame_back = 0;
frame_forward = 40;

surround_random_all = zeros(size(roi_trace_df,1),n_sample);
frames_random = zeros(n_sample, length(frames_reward_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_df,2),length(frames_reward_all));
end

try
    matlabpool
catch me
end

parfor curr_roi = 1:length(spike_frames)
    frames_random_spike = [];
    frames_random_spike = repmat(frames_random,[1 1 length(spike_frames{curr_roi})]) ;
    % get difference between random frames and spike t0s
    for curr_spike = 1:length(spike_frames{curr_roi});
        frames_random_spike(:,:,curr_spike) = frames_random_spike(:,:,curr_spike) - spike_frames{curr_roi}(curr_spike);
    end
    % find any of that matrix that are between frame_back and frame_forward
    frames_random_spike = frames_random_spike > -frame_back & frames_random_spike < frame_forward;
    surround_random_all(curr_roi,:) = sum(any(frames_random_spike,3),2);
    disp(['ROI ' num2str(curr_roi) '/' num2str(length(spike_frames))]);
end

% make rewarded events binary
surround_reward_binary = surround_reward;
surround_reward_binary(surround_reward_binary > 0) = 1;
surround_reward_binary = sum(surround_reward_binary,2);

p = zeros(size(roi_trace_df,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_df,1)
    % tack on test value to end of random distribution
    distr_test = [surround_random_all(curr_roi,:) surround_reward_binary(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end

mod_up_cells = [];
mod_down_cells = [];

mod_up_cells = find(p > 0.975);
mod_down_cells = find(p < 0.025);
toc

%% show up/down/non-modulated cells

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
end

for i = 1:length(mod_up_cells)
    patch(polygon.ROI{mod_up_cells(i)}(:,1),-polygon.ROI{mod_up_cells(i)}(:,2),[0 1 0]);
end

for i = 1:length(mod_down_cells)
    patch(polygon.ROI{mod_down_cells(i)}(:,1),-polygon.ROI{mod_down_cells(i)}(:,2),[1 0 0]);
end
ylim([-128 0]);
xlim([0 512]);
title('Cells correlated with trained task, p > 0.05')

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1-p(i) p(i) 0]);
end
ylim([-128 0]);
xlim([0 512]);
title('Cells correlated with trained task, graded')

%% Plot by kmeans grouping for general trace
% make traces just with delta functions
roi_trace_long_spikes = zeros(size(roi_trace_long));
spike_frames_long = cell(length(spike_frames));
for i = 1:length(spike_frames);
    for j = 1:length(spike_frames{i})
        spike_frames_long{i} =  [spike_frames_long{i} spike_frames{i}(j)-10:spike_frames{i}(j)+10];
    end
end

for i = 1:length(spike_frames_long);
    roi_trace_long_spikes(i,spike_frames_long{i}) = 1;
end

kmeans_var = {};
for cluster_num = 1:12
    kmeans_var{cluster_num} = zeros(length(temp1),1);
    try
        [temp1 centers sumd d] = kmeans(roi_trace_long_spikes,cluster_num);
    catch me
        continue
    end
    for i = 1:length(temp1);
        kmeans_var{cluster_num}(i) = abs(roi_trace_long_spikes(i,:) - centers(temp1(i),:))/(roi_trace_long_spikes(i,:) + centers(temp1(i),:));
    end
    disp([num2str(cluster_num)]);
end

figure;hold on;
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

for i = 1:length(a)
    patch(polygon.ROI{a(i)}(:,1),polygon.ROI{a(i)}(:,2),[1 0 0]);
end
for i = 1:length(b)
    patch(polygon.ROI{b(i)}(:,1),polygon.ROI{b(i)}(:,2),[0 1 0]);
end
for i = 1:length(c)
    patch(polygon.ROI{c(i)}(:,1),polygon.ROI{c(i)}(:,2),[0 0 1]);
end
title('Cells sorted through K Means')

% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    temp_trace(i,:) = roi_trace_long(sort_temp1(i,2),:);
end
figure;imagesc(temp_trace);
figure;imagesc(corrcoef(temp_trace'))
figure;plot(centers(1,:));hold on;plot(centers(2,:),'--r');
title('K-means centers');
legend('Center 1','Center 2');

%% Plot by kmeans grouping for which rewards had events around them
surround_reward_bin = surround_reward;
surround_reward_bin(surround_reward > 0) = 1;

[temp1 centers] = kmeans(surround_reward_bin,2);
figure;hold on;
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

for i = 1:length(a)
    patch(polygon.ROI{a(i)}(:,1),polygon.ROI{a(i)}(:,2),[1 0 0]);
end
for i = 1:length(b)
    patch(polygon.ROI{b(i)}(:,1),polygon.ROI{b(i)}(:,2),[0 1 0]);
end
for i = 1:length(c)
    patch(polygon.ROI{c(i)}(:,1),polygon.ROI{c(i)}(:,2),[0 0 1]);
end
title('Cells sorted through K Means')

% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(surround_reward));
for i = 1:size(surround_reward,1);
    temp_trace(i,:) = surround_reward(sort_temp1(i,2),:);
end
figure;imagesc(temp_trace);
figure;imagesc(corrcoef(temp_trace'))
figure;plot(centers(1,:));hold on;plot(centers(2,:),'--r');
title('K-means centers');
legend('Center 1','Center 2');

% plot the correlations between groups

corrcoef_1  = corrcoef(roi_trace_long(a,:)');
corrcoef_1(corrcoef_1 == 1) = NaN;
mean_corrcoef_1 = nanmean(corrcoef_1(:));
sem_corrcoef_1 = nanstd(corrcoef_1(:))./sqrt(size(roi_trace_long,1));

corrcoef_2  = corrcoef(roi_trace_long(b,:)');
corrcoef_2(corrcoef_2 == 1) = NaN;
mean_corrcoef_2 = nanmean(corrcoef_2(:));
sem_corrcoef_2 = nanstd(corrcoef_2(:))./sqrt(size(roi_trace_long,1));

corrcoef_all  = corrcoef(roi_trace_long');
corrcoef_all(corrcoef_all == 1) = NaN;
mean_corrcoef_all = nanmean(corrcoef_all(:));
sem_corrcoef_all = nanstd(corrcoef_all(:))./sqrt(size(roi_trace_long,1));

means = [mean_corrcoef_1 mean_corrcoef_2 mean_corrcoef_all];
sems = [sem_corrcoef_1 sem_corrcoef_2 sem_corrcoef_all];

figure; hold on;
h = bar(means)
% this can't be copied and pasted, very lame
%ch = get(h,'Children');
%set(ch,'FaceVertexCData',[1 0 0 ; 0 1 0; .5 .5 .5])
errorb(means,sems,'top')
colormap(gray)
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'Group 1','Group 2','All'});
colormap(gray)
ylim([0 0.3])
ylabel('Correlation')

%% Plot by affinity clustering grouping for which rewards had events around them
N = size(surround_reward_bin,1);
x = surround_reward_bin;
M=N*N-N; s=zeros(M,3);
j=1;
for i=1:N
    for k=[1:i-1,i+1:N]
        s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
        j=j+1;
    end;
end
p=median(s(:,3)); % Set preference to median similarity
[idx,netsim,dpsim,expref]=apcluster(s,p,'plot');
fprintf('Number of clusters: %d\n',length(unique(idx)));
fprintf('Fitness (net similarity): %f\n',netsim);
figure; % Make a figures showing the data and the clusters
for i=unique(idx)'
    ii=find(idx==i); h=plot(x(ii,1),x(ii,2),'o'); hold on;
    col=rand(1,3); set(h,'Color',col,'MarkerFaceColor',col);
    xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii));
    line([x(ii,1),xi1]',[x(ii,2),xi2]','Color',col);
end

surround_reward_bin_order = surround_reward_bin;
surround_reward_bin_order(:,size(surround_reward_bin,2)+1) = idx;
suround_reward_bin_order = sort(surround_reward_bin_order,size(surround_reward_bin,2));
figure;
surround_reward_bin_order = surround_reward_bin_order(:,1:end-1);
imagesc(surround_reward_bin_order);
colormap(gray);


%% Find unrewarded-related cells
% find which ROIs had events around given timeframe
frame_back = 0;
frame_forward = 20;
% events around reward/left up
surround_unrewarded = zeros(size(roi_trace_long,1),length(frames_unrewarded_all));
for curr_reward = 1:length(frames_unrewarded_all)
    curr_frame = frames_unrewarded_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_unrewarded(curr_roi,curr_reward) = event_detect_amp;
        end
    end
end
figure;imagesc(surround_unrewarded);
title('Events around reward/left up')

% get behavior-related cells:

% events around left down during air
surround_left_down_air = zeros(size(roi_trace_long,1),length(frames_left_down_air_all));
for curr_left_down_air = 1:length(frames_left_down_air_all)
    curr_frame = frames_left_down_air_all(curr_left_down_air);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_left_down_air(curr_roi,curr_left_down_air) = event_detect_amp;
        end
    end
end


tic
n_sample = 10000;

surround_random_all = zeros(size(roi_trace_long,1),n_sample);
frames_random = zeros(n_sample, length(frames_unrewarded_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_long,2),length(frames_unrewarded_all));
end

w = waitbar(0)
for curr_roi = 1:length(spike_frames)
    frames_random_spike = [];
    frames_random_spike = repmat(frames_random,[1 1 length(spike_frames{curr_roi})]) ;
    % get difference between random frames and spikes
    for curr_spike = 1:length(spike_frames{curr_roi});
        frames_random_spike(:,:,curr_spike) = frames_random_spike(:,:,curr_spike) - spike_frames{curr_roi}(curr_spike);
    end
    % find any of that matrix that are between frame_back and frame_forward
    frames_random_spike = frames_random_spike > -frame_back & frames_random_spike < frame_forward;
    surround_random_all(curr_roi,:) = sum(any(frames_random_spike,3),2);
    disp(['ROI' num2str(curr_roi)]);
    waitbar(curr_roi/length(spike_frames),w);
end
close(w);

% make rewarded events binary
surround_unrewarded_binary = surround_unrewarded;
surround_unrewarded_binary(surround_unrewarded_binary > 0) = 1;
surround_unrewarded_binary = sum(surround_unrewarded_binary,2);

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [surround_random_all(curr_roi,:) surround_unrewarded_binary(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end

mod_up_cells_unrewarded = [];
mod_down_cells_unrewarded = [];

mod_up_cells_unrewarded = find(p > 0.95);
mod_down_cells_unrewarded = find(p < 0.05);
toc


%% CALTRACE EXPORT: Divide caltrace conts y by 1/4 and call them polygon.ROI
polygon.ROI = cell(size(CONTS));
for i = 1:length(CONTS)
    polygon.ROI{i}(:,1) = CONTS{i}(:,1);
    polygon.ROI{i}(:,2) = CONTS{i}(:,2)/4;
end

%% For reward traces, plot stacks
figure;
hold on;
for i = 1:size(reward_trace_concat,1)
    plot(reward_trace_concat(i,:,1) + i,'-k');
end

%% Make raster plot of events for roi_trace_long
figure
hold on;
for i = 1:size(roi_trace_long,1)
    if ~isempty(spike_frames{i});
        max_amp = max([spike_amplitudes{i}]);
        for j = 1:size(spike_frames{i});
            line([spike_frames{i}(j) spike_frames{i}(j)], ...
                [i-0.5 i+0.5], 'color',...
                [1-spike_amplitudes{i}(j)/max_amp 1-spike_amplitudes{i}(j)/max_amp 1-spike_amplitudes{i}(j)/max_amp]);
        end
    end
end


%% cluster which rewarded trials had particular spikes (ONLY FOR MOD UP CELLS!)

frame_back = 10;
frame_forward = 20;
% events around reward/left up
surround_reward = zeros(mod_up_cells,length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for num_roi = 1:length(mod_up_cells)
        curr_roi = mod_up_cells(num_roi);
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            surround_reward(num_roi,curr_reward) = 1;
        end
    end
end

surround_reward_bin = surround_reward;
surround_reward_bin(surround_reward_bin > 0) = 1;
surround_reward_bin(surround_reward_bin < 0) = 0;
surround_reward_bin = surround_reward_bin';

z = linkage(surround_reward_bin,'ward');
c = cluster(z,'MaxClust',3);

temp = [surround_reward_bin c];
temp = sortrows(temp,size(temp,2));
figure;imagesc(temp);colormap(jet);
figure;imagesc([corrcoef(temp(:,1:end-1)') sort(c)]);

%%%%
% reaction time from air onset to lever press

times = saved_history.RewardsSection_LastTrialEvents;
odor_rxn_time_curr = [];
odor_rxn_times = [];
odor_rxn_trials = [];
for curr_trial = reward_trial;
    if ~ismember(curr_trial, trial_notRecorded)
        begin_41 = [];
        down_41 = [];
        
        begin_41 = find(times{curr_trial}(:,1) == 41,1);
        down_41 = find(times{curr_trial}(:,1) == 41 &...
            times{curr_trial}(:,2) == 5,1);
        if ~isempty(down_41);
            odor_rxn_time_curr = times{curr_trial}(down_41,3)...
                - times{curr_trial}(begin_41,3);
            odor_rxn_times = [odor_rxn_times;odor_rxn_time_curr];
            odor_rxn_trials = [odor_rxn_trials; curr_trial];
        end
    end
end

% Hold time for the press that counted

times = saved_history.RewardsSection_LastTrialEvents;
air_triggered_hold = [];
for curr_trial = reward_trial;
    first_down_41 = [];
    first_up_post41 = [];
    curr_air_triggered_hold = [];
    
    first_down_41 = find(times{curr_trial}(:,1) == 41 ...
        & times{curr_trial}(:,2) == 5,1);
    
    first_up_post41 = find((times{curr_trial}(:,1) == 42 ...
        | times{curr_trial}(:,1) == 43) ...
        & times{curr_trial}(:,2) == 6,1);
    
    curr_air_triggered_hold = times{curr_trial}(first_up_post41,3)...
        - times{curr_trial}(first_down_41,3);
    
    if ~isempty(curr_air_triggered_hold);
        air_triggered_hold = [air_triggered_hold;...
            curr_air_triggered_hold];
    end
end

% Get number of presses during the trial before the air came on
AP_getLeftLeverTimes;

times = saved_history.RewardsSection_LastTrialEvents;
pre_air_presses = [];
for curr_trial = reward_trial
    begin_40 = [];
    down_41 = [];
    
    begin_40 = find(times{curr_trial}(:,1) == 40,1);
    down_41 = find(times{curr_trial}(:,1) == 41 &...
        times{curr_trial}(:,2) == 5,1);
    
    begin_40 = times{curr_trial}(begin_40,3);
    down_41 = times{curr_trial}(down_41,3);
    
    pre_air_presses_curr = length(find(left_down > begin_40 & ...
        left_down < down_41));
    
    pre_air_presses = [pre_air_presses;pre_air_presses_curr];
end
%%%%

% get behav params
figure
subplot(3,1,1);
boxplot(odor_rxn_times,c,'notch','on','plotstyle','traditional','whisker',0,'colors',[0 0 0],'symbol','.k');
title('Reaction time')
subplot(3,1,2);
boxplot(air_triggered_hold,c,'notch','on','plotstyle','traditional','whisker',0,'colors',[0 0 0],'symbol','.k');
title('Reaction time')
subplot(3,1,3);
boxplot(pre_air_presses,c,'notch','on','plotstyle','traditional','whisker',0,'colors',[0 0 0],'symbol','.k');

%% try to find modulated cells with xcorr
% this was stupid and a waste of time, it's the same thing as just
% averaging activity around a reward

max_lag = 60; % frames for maximum lag
n_sample = 50;

frames_random = zeros(n_sample, length(frames_reward_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_long,2),length(frames_reward_all));
end

reward_time_real = zeros(1,size(roi_trace_long,2));
reward_time_real(frames_reward_all) = 1;
xcorr_grid = [];
for i = 1:size(roi_trace_long,1)
    xcorr_grid(i,:) = xcorr(roi_trace_long(i,:),reward_time_real,max_lag,'coeff');
end

% % random event times
% xcorr_max_rand = cell(size(roi_trace_long,1));
% for cell_num = 1:size(roi_trace_long,1)
%     xcorr_max_rand{cell_num} = zeros(n_sample,1);
%     for i = 1:n_sample
%     reward_time = zeros(1,size(roi_trace_long,2));
%     reward_time(frames_random(i,:)) = 1;
%     xcorr_max_rand{cell_num}(i) = max(xcorr(roi_trace_long(cell_num,:),reward_time,max_lag));
%     disp(['Finished sample ' num2str(i) ' : ' num2str(cell_num)]);
%     end
% end

% jackknife?
xcorr_max_jack = cell(size(roi_trace_long,1));
for cell_num = 1:size(roi_trace_long,1)
    xcorr_max_jack{cell_num} = zeros(n_sample,1);
    for i = 1:length(frames_reward_all)
        reward_time_one_out = reward_time_real;
        reward_time_one_out(frames_reward_all(i)) = 0;
        xcorr_max_jack{cell_num}(i) = max(xcorr(roi_trace_long(cell_num,:),reward_time_one_out,max_lag,'coeff'));
        disp(['Finished sample ' num2str(i) ' : ' num2str(cell_num)]);
    end
end

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [xcorr_max_jack{curr_roi}; max(xcorr_grid(curr_roi,:))];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end

mod_cells = [];
mod_cells = find(p > 0.95);


%% Get pairwise coincident firing during rewards for each cell

% frames to search for reward
frame_back = 10;
frame_forward = 10;

coincident_pairs = zeros(length(mod_up_cells));

% events around reward/left up
surround_reward = zeros(mod_up_cells,length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for num_roi = 1:length(mod_up_cells)
        curr_roi = mod_up_cells(num_roi);
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            surround_reward(num_roi,curr_reward) = 1;
        end
    end
end

for cell1 = 1:length(mod_up_cells)
    for cell2 = 1:length(mod_up_cells)
        coincident_pairs(cell1,cell2) = surround_reward(cell1,:)*surround_reward(cell2,:)'/...
            min([sum(surround_reward(cell1,:)) sum(surround_reward(cell2,:))]);
    end
end
figure;imagesc(coincident_pairs);title('Coincident pairs')

% get raw correlation, for comparison
correlation_pairs = zeros(length(mod_up_cells));
for cell1 = 1:length(mod_up_cells)
    for cell2 = 1:length(mod_up_cells)
        corrcoef_curr = corrcoef(roi_trace_long(mod_up_cells(cell1),:),roi_trace_long(mod_up_cells(cell2),:));
        correlation_pairs(cell1,cell2) = corrcoef_curr(1,2);
    end
end
figure;imagesc(correlation_pairs);title('Correlation pairs')



%%  find  modulated cells by percent of trials with rewards

% define how many rewards have to have events around them
percent_reward = 0.2;

mod_cells_early = [];
mod_cells_middle = [];
mod_cells_late = [];

mod_cells_early_percent = [];
mod_cells_middle_percent = [];
mod_cells_late_percent = [];

% frames to search for reward (early)
frame_back = 10;
frame_forward = 10;

surround_reward_all = zeros(size(roi_trace_long,1),length(frames_reward_all));

% events around reward/left up
surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:size(roi_trace_long,1)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            surround_reward(curr_roi,curr_reward) = 1;
            surround_reward_all(curr_roi,curr_reward) = 1;
        end
    end
end

sum_rewards_early = sum(surround_reward,2);
mod_cells_early = find(sum_rewards_early > percent_reward*length(frames_reward_all));
mod_cells_early_percent = sum_rewards_early(mod_cells_early)./length(frames_reward_all);

% frames to search for reward (middle)
frame_back = -15;
frame_forward = 30;

% events around reward/left up
surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:size(roi_trace_long,1)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            surround_reward(curr_roi,curr_reward) = 1;
            surround_reward_all(curr_roi,curr_reward) = 2;
        end
    end
end

sum_rewards_middle = sum(surround_reward,2);
mod_cells_middle = find(sum_rewards_middle > percent_reward*length(frames_reward_all));
mod_cells_middle_percent = sum_rewards_middle(mod_cells_middle)./length(frames_reward_all);

% frames to search for reward (late)
frame_back = -30;
frame_forward = 45;

% events around reward/left up
surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:size(roi_trace_long,1)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            surround_reward(curr_roi,curr_reward) = 1;
            surround_reward_all(curr_roi,curr_reward) = 3;
        end
    end
end

sum_rewards_late = sum(surround_reward,2);
mod_cells_late = find(sum_rewards_late > percent_reward*length(frames_reward_all));
mod_cells_late_percent = sum_rewards_late(mod_cells_late)./length(frames_reward_all);

% if there are duplicates, take which category is higher, first for
% early/middle, then middle/late, then early/late
overlap_cells = intersect(mod_cells_early,mod_cells_middle);
for i = overlap_cells';
    cell_indx1 = find(mod_cells_early == i);
    cell_indx2 = find(mod_cells_middle == i);
    check_more = sum_rewards_early(cell_indx1) > sum_rewards_middle(cell_indx2);
    if check_more == 1;
        mod_cells_middle(cell_indx2) = [];
    elseif check_more == 0;
        mod_cells_early(cell_indx1) = [];
    end
end

overlap_cells = intersect(mod_cells_middle,mod_cells_late);
for i = overlap_cells';
    cell_indx1 = find(mod_cells_middle == i);
    cell_indx2 = find(mod_cells_late == i);
    check_more = sum_rewards_middle(cell_indx1) > sum_rewards_late(cell_indx2);
    if check_more == 1;
        mod_cells_late(cell_indx2) = [];
    elseif check_more == 0;
        mod_cells_middle(cell_indx1) = [];
    end
end


overlap_cells = intersect(mod_cells_early,mod_cells_late);
for i = overlap_cells';
    cell_indx1 = find(mod_cells_early == i);
    cell_indx2 = find(mod_cells_late == i);
    check_more = sum_rewards_early(cell_indx1) > sum_rewards_late(cell_indx2);
    if check_more == 1;
        mod_cells_late(cell_indx2) = [];
    elseif check_more == 0;
        mod_cells_early(cell_indx1) = [];
    end
end

%% Check if modulated cell significantly changed relative timing

shift_cells = [];

frame_back = 10;
frame_forward = 10;

for curr_mod_cell = mod_up_cells'
    
    event_time = [];
    surround_reward_trace = [];
    
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        
        % get all of the traces
        if any(curr_frame_search > size(roi_trace_long,2));
            curr_frame_search(curr_frame_search > size(roi_trace_long,2)) = size(roi_trace_long,2);
        end
        
        % get only the traces with events
        event_detect = [];
        event_detect = ismember(spike_frames{curr_mod_cell},curr_frame_search);
        if any(event_detect)
            % the event start time is when it was within 1 std of baseline
            event_time = [event_time (min(spike_frames{curr_mod_cell}(event_detect)) - curr_frame)];
            surround_reward_trace = [surround_reward_trace;roi_trace_long(curr_mod_cell,curr_frame_search)];
        end
    end
    
    % check if there was drift, compare first third to last third
    third_trials = floor(length(event_time)/3);
    h = ttest2(event_time(1:third_trials),event_time(end-third_trials:end));
    
    if h == 1;
        shift_cells = [shift_cells curr_mod_cell];
    end
end


%% Coincident firing by spikes and multiplication

curr_roi = 44;

% find by multiplication
roi_trace_mult = [];
roi_trace_rep = [];
roi_trace_rep = repmat(roi_trace_long(curr_roi,:),size(roi_trace_long,1),1);
roi_trace_mult = sqrt(abs(roi_trace_rep.*roi_trace_long));

% look by surround
frame_back = 5;
frame_forward = 5;
roi_trace_coincident = zeros(size(roi_trace_long));
for i = 1:length(spike_frames{curr_roi})
    frame_search = spike_frames{curr_roi}(i)-frame_back:spike_frames{curr_roi}(i)+frame_forward;
    for j = 1:size(roi_trace_long,1)
        event_detect = [];
        event_detect = ismember(spike_frames{j},frame_search);
        if any(event_detect)
            roi_trace_coincident(j,spike_frames{curr_roi}(i)) = 1;
        end
    end
end
figure;
imagesc(roi_trace_coincident);
colormap(gray);

figure;
hold on;
for i = 1:size(roi_trace_coincident,1)
    plot(find(roi_trace_coincident(i,:) == 1), i*ones(size(find(roi_trace_coincident(i,:) == 1))),'.','MarkerSize',5)
end


% get coincidence index for each ROI (from Komiyama 2010)
num_frames = size(roi_trace_long,2);
num_frames_search = size(roi_trace_long,2)/(frame_back+frame_forward+1);
coincidence_indx = zeros(size(roi_trace_long,1),1);
coincidence_indx_tk = zeros(size(roi_trace_long,1),1);
spike_1 = length(spike_frames{curr_roi});
spike_1_prob = spike_1/num_frames;
for i = 1:size(roi_trace_long,1)
    spike_2 = length(spike_frames{i});
    spike_2_prob = spike_2/num_frames;
    % this is what they do in the paper... does this make sense?
    geometric_mean = num_frames*sqrt(spike_1_prob*spike_2_prob);
    disp('CHECK THIS: SHOULDNt BE 1+2, should be COINC1+ COINC2')
    %``````````````````````CHECK THIS!```````````````````````````%
    coincidence_indx_tk(i) = ((spike_1+spike_2)-(num_frames*spike_1_prob*spike_2_prob))/(geometric_mean);
    % this makes more intuitive sense to me, co-spikes/all spikes
    co_spike = sum(roi_trace_coincident(i,:));
    coincidence_indx(i) = (co_spike*2)/(spike_1+spike_2);
end

figure;
imagesc(roi_trace_mult);
caxis([0 1]);
title(num2str(curr_roi));

z = linkage(roi_trace_mult,'ward');
c = cluster(z,'MaxClust',5);

temp = [roi_trace_mult [1:size(roi_trace_mult,1)]' c];
temp = sortrows(temp,size(temp,2));
figure;imagesc(temp);colormap(jet);
caxis([0 1])

% plot rewards as dashed red lines
for i = 1:length(frames_reward_all)
    line([frames_reward_all(i) frames_reward_all(i)],ylim,'linestyle','--','color','r');
end

% plot non-rewards as dashed magenta lines
for i = 1:length(frames_unrewarded_all)
    line([frames_unrewarded_all(i) frames_unrewarded_all(i)],ylim,'linestyle','--','color','m');
end

%% show up/down/non-modulated cells from dumb-way detection

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
end

for i = 1:length(mod_cells_early)
    patch(polygon.ROI{mod_cells_early(i)}(:,1),-polygon.ROI{mod_cells_early(i)}(:,2),[0 1 0]);
end

for i = 1:length(mod_cells_middle)
    patch(polygon.ROI{mod_cells_middle(i)}(:,1),-polygon.ROI{mod_cells_middle(i)}(:,2),[1 1 0]);
end

for i = 1:length(mod_cells_late)
    patch(polygon.ROI{mod_cells_late(i)}(:,1),-polygon.ROI{mod_cells_late(i)}(:,2),[1 0 0]);
end

ylim([-128 0]);
xlim([0 512]);
title('Cells correlated with lever press')

%% Get reward-related cells the smart way, early/middle/late
frame_back = 10;
frame_forward = 10;
% events around reward/left up
surround_reward_early = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(roi_trace_long(curr_roi,spike_frames{curr_roi}(event_detect)));
            surround_reward_early(curr_roi,curr_reward) = 1;
        end
    end
end

frame_back = -10;
frame_forward = 30;
% events around reward/left up
surround_reward_middle = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(roi_trace_long(curr_roi,spike_frames{curr_roi}(event_detect)));
            surround_reward_middle(curr_roi,curr_reward) = 1;
        end
    end
end

frame_back = -30;
frame_forward = 50;
% events around reward/left up
surround_reward_late = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(roi_trace_long(curr_roi,spike_frames{curr_roi}(event_detect)));
            surround_reward_late(curr_roi,curr_reward) = 1;
        end
    end
end

n_sample = 10000;
disp('Currently hardcoding movie length')

% for generating random spike matrix
frame_back = 10;
frame_forward = 10;

surround_random_all = zeros(size(roi_trace_long,1),n_sample);
frames_random = zeros(n_sample, length(frames_reward_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_long,2),length(frames_reward_all));
end

try
    matlabpool
catch me
end

% find how often random spikes correlate with rewards
parfor curr_roi = 1:length(spike_frames)
    % get difference between random frames and spike frames
    % first, repeat matrix over the course of the frame search
    frames_random_spike = repmat(frames_random,[1 1 frame_back+frame_forward+1]);
    % second, make the random frames span the range of frame searches
    frame_search_linear = [-frame_back:frame_forward];
    for frame_search_fix = 1:size(frames_random_spike,3)
        frames_random_spike(:,:,frame_search_fix) = frames_random_spike...
            (:,:,frame_search_fix)+frame_search_linear(frame_search_fix);
    end
    
    frames_random_spike_coincident = zeros(size(frames_random_spike));
    for i = 1:size(frames_random_spike,3)
        frames_random_spike_coincident(:,:,i) = ismember(frames_random_spike(:,:,i),spike_frames{curr_roi});
    end
    frames_random_spike_coincident = sum(frames_random_spike_coincident,3);
    surround_random_all(curr_roi,:) = sum(frames_random_spike_coincident,2)';
    disp(['ROI ' num2str(curr_roi) '/' num2str(length(spike_frames))]);
end

p_early = zeros(size(roi_trace_long,1),1);
p_middle = zeros(size(roi_trace_long,1),1);
p_late = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test_early = [surround_random_all(curr_roi,:) sum(surround_reward_early(curr_roi),2)];
    distr_test_middle = [surround_random_all(curr_roi,:) sum(surround_reward_middle(curr_roi),2)];
    distr_test_late = [surround_random_all(curr_roi,:) sum(surround_reward_late(curr_roi),2)];
    % rank all values, nonparametric way to get percentile
    distr_rank_early = tiedrank(distr_test_early);
    distr_rank_middle = tiedrank(distr_test_middle);
    distr_rank_late = tiedrank(distr_test_late);
    % percentile is the ratio of the test value to all total ranks
    p_early(curr_roi) = distr_rank_early(end)/length(distr_rank_early);
    p_middle(curr_roi) = distr_rank_middle(end)/length(distr_rank_middle);
    p_late(curr_roi) = distr_rank_late(end)/length(distr_rank_late);
end

mod_cells_early = [];
mod_cells_middle = [];
mod_cells_late = [];

mod_cells_early = find(p_early > 0.95);
mod_cells_middle = find(p_middle > 0.95);
mod_cells_late = find(p_late > 0.95);


%% load and align tdms data (from force lever)

[tdms_filename,tdms_path]=uigetfile('*.tdms','Choose TDMS file','Multiselect','off');
tdms_filename = [tdms_path tdms_filename];

lever_force = TDMS_getStruct(tdms_filename);

figure; hold on;
% plot the raw force data in blue
plot(lever_force.lever.force.data,'.b');
plot(lever_force.lever.force.data);
% plot raw force data over threshold in red
% the cross threshhold shifts by 1 for some reason? fixes this:
threshhold_cross = find(lever_force.lever.threshhold_cross.data == 1)+1;
plot(threshhold_cross, lever_force.lever.force.data(threshhold_cross),'.r');
line(xlim,[lever_force.Props.force_threshold__V_ lever_force.Props.force_threshold__V_], ...
    'linestyle','--','color','k');

% get all raw times that solo said lever activity
AP_getForceLeverTimes;

% get all times for crossing threshold, +2 at end for diff + 1 shift
cross_indx_over = find(diff(lever_force.lever.threshhold_cross.data)== 1)+2;
cross_indx_under = find(diff(lever_force.lever.threshhold_cross.data)== -1)+2;
cross_indx_length = cross_indx_under - cross_indx_over;

% eliminate cross_indx when there was only one point above (solo misses)
cross_indx_single = find(cross_indx_length == 1);
cross_indx_over(cross_indx_single) = [];
cross_indx_under(cross_indx_single) = [];

% if they're not the same size, stop b/c something's up
if length(cross_indx_over) ~= length(left_down)
    disp('Threshhold crossings and Dispatcher lever reports not equal')
    keyboard
end

% translate force linearly: t = b + mx
b = left_down(1) - cross_indx_over(1);
m = (left_down(end) - left_down(1))/(cross_indx_over(end) - cross_indx_over(1));

threshhold_cross_over_ms_ballpark =  (cross_indx_over+b-left_down(1))*m+left_down(1);
threshhold_cross_under_ms_ballpark =  (cross_indx_under+b-left_up(1))*m+left_up(1);

%% Draw ROIs, put numbers

% plot the ROIs
figure
hold on;
set(gca,'color','k');
for i = 1:length(polygon.ROI)
    r = patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
    %set(r,'linewidth',length(spike_frames{i})/10+0.01);
    set(r,'edgecolor','w');
end

% get centers of all ROIs
roi_centers = zeros(size(roi_trace_long),2);
roi_polygon_centers = zeros(size(roi_trace_long),2);
for i = 1:size(roi_trace_long,1);
    % get center x and y from polygon boundaries
    center_x = mean(polygon.ROI{i}(:,1));
    center_y = mean(polygon.ROI{i}(:,2));
    roi_polygon_centers(i,1) = center_x;
    roi_polygon_centers(i,2) = -center_y;
end

% write numbers in the centers for number of spikes
for i = 1:length(polygon.ROI);
    curr_num_spike = length(spike_frames{i});
    roi_numbers_h(i) = text(roi_polygon_centers(i,1),roi_polygon_centers(i,2),num2str(i),'color','white');
end


%% Attempt at autodetecting: find max of smooth
% requires image tensor im
% this breaks laptop: don't use there
Y = size(im,1);
X = size(im,2);
im_max = zeros(Y,X);
w = waitbar(0);
for x_curr = 1:X
    for y_curr = 1:Y
        temp_trace = double(im(X,Y,:));
        temp_trace = smooth(temp_trace,10);
        im_max(y_curr,x_curr) = max(temp_trace);
    end
    waitbar(x_curr/X,w);
end



