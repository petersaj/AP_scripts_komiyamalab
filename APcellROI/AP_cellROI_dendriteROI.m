function AP_cellROI_dendriteROI

%% Load and align summed movies

save_flag = true;

[summed_filenames, summed_path] = uigetfile('*.tif','Choose tiff files','multiselect','on');
if ~iscell (summed_filenames)
    summed_filenames = {summed_filenames};
end
summed_filenames = sort(summed_filenames);

n_sessions = length(summed_filenames);

dir_sep = strfind(summed_path,filesep);
animal = summed_path(dir_sep(end-1)+1:dir_sep(end)-1);

M = 512;
N = 512;

% If more than 1 sessions, align all to the first session
if n_sessions > 1
    disp('Loading movies for alignment')
    
    summed_max = double(nan(M,N,n_sessions));
    summed_mean = double(nan(M,N,n_sessions));
    summed_borders = zeros(M,N,n_sessions);
    % Loop through sessions, get max/avg projection
    for curr_session = 1:n_sessions
        
        % Get max projection of day
        curr_file = [summed_path filesep summed_filenames{curr_session}];
        
        im_info = imfinfo(curr_file);
        n_frames = length(im_info);
        
        % Load in current day
        im = zeros(M,N,n_frames,'uint16');
        for i = 1:n_frames
            im(:,:,i) = imread(curr_file,'tiff',i,'Info',im_info);
        end
        
        % Store mean/max
        summed_max(:,:,curr_session) = max(im,[],3);
        summed_mean(:,:,curr_session) = mean(im,3);
        
        % Get the borders with motion artifacts (all zero rows/columns)
        summed_borders(:,any(all(im == 0,1),3),curr_session) = 1;
        summed_borders(any(all(im == 0,2),3),:,curr_session) = 1;
        
        disp(['Loading for alignment, session: ' num2str(curr_session)]);
        
        clear im
    end
    
    
    % Register days to each other
    disp('Registering maximum projections')
    tform_matrix = cell(n_sessions,1);
    tform_matrix{1} = eye(3);
    
    summed_max_reg = nan(size(summed_max));
    summed_max_reg(:,:,1) = summed_max(:,:,1);
    summed_mean_reg = nan(size(summed_mean));
    summed_mean_reg(:,:,1) = summed_mean(:,:,1);
    for curr_session = 2:n_sessions
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumStepLength = 0.02;
        %optimizer.MaximumStepLength = 0.0001;
        optimizer.RelaxationFactor = 0.1;
        optimizer.GradientMagnitudeTolerance = 1e-5;
        optimizer.MaximumIterations = 300;
%         optimizer = registration.optimizer.OnePlusOneEvolutionary();
%         optimizer.MaximumIterations = 500;
%         optimizer.GrowthFactor = 1.00001;   
%         optimizer.InitialRadius = 1e-5;
        
        % Perform the registration on the maximum image
        tformEstimate = imregtform(summed_max(:,:,curr_session),summed_max(:,:,1),'affine',optimizer,metric);
        %tformEstimate = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');
        
        max_im_reg = imwarp(summed_max(:,:,curr_session),tformEstimate,'Outputview',imref2d([M N]));
        mean_im_reg = imwarp(summed_mean(:,:,curr_session),tformEstimate,'Outputview',imref2d([M N]));
        
        tform_matrix{curr_session} = tformEstimate.T;
        summed_max_reg(:,:,curr_session) = max_im_reg;
        summed_mean_reg(:,:,curr_session) = mean_im_reg;
        
%         %%%%%%%%%%%%%
%         % This is to do correlation, then affine (if above doesn't work)
%         [optimizer, metric] = imregconfig('monomodal');
%         optimizer = registration.optimizer.OnePlusOneEvolutionary();
%         optimizer.MaximumIterations = 500;
%         optimizer.GrowthFactor = 1.00001;   
%         optimizer.InitialRadius = 1e-5;
%         
%         % Register 1) correlation
%         tformEstimate_corr = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');       
%         max_im_reg_corr = imwarp(summed_max(:,:,curr_session),tformEstimate_corr,'Outputview',imref2d([M N]));
%         
%         % Register 2) affine
%         tformEstimate_affine = imregtform(max_im_reg_corr,summed_max(:,:,1),'affine',optimizer,metric);
%         
%         tformEstimate_combined = tformEstimate_corr;
%         tformEstimate_combined.T = tformEstimate_affine.T*tformEstimate_corr.T;
% 
%         max_im_reg = imwarp(summed_max(:,:,curr_session),tformEstimate_combined,'Outputview',imref2d([M N]));
%         mean_im_reg = imwarp(summed_mean(:,:,curr_session),tformEstimate_combined,'Outputview',imref2d([M N]));
%         
%         tform_matrix{curr_session} = tformEstimate_combined.T;
%         summed_max_reg(:,:,curr_session) = max_im_reg;
%         summed_mean_reg(:,:,curr_session) = mean_im_reg;
        
    end
    
    % Human-verify registration (if it isn't near-perfect, don't use)
    gui_fig = AP_image_scroll(summed_max_reg);
    check_alignment = input('Perfect alignment (y/n)?: ','s');
    
    if ~strcmp(check_alignment,'y')
        keyboard
    end
    
    close(gui_fig);
    drawnow;
    
else
    
    tform_matrix = cell(n_sessions,1);
    tform_matrix{1} = eye(3);
    
    curr_file = [summed_path filesep summed_filenames{1}];
    im_info = imfinfo(curr_file);
    n_frames = length(im_info);
    
    % Load in current day
    im = zeros(M,N,n_frames,'uint16');
    for i = 1:n_frames
        im(:,:,i) = imread(curr_file,'tiff',i,'Info',im_info);
    end
    % Get the borders with motion artifacts (all zero rows/columns)
    summed_borders = zeros(M,N);
    summed_borders(:,any(all(im == 0,1),3)) = 1;
    summed_borders(any(all(im == 0,2),3),:) = 1;
    summed_max_reg = max(im,[],3);
    clear im
    
end

%% Get active spots

% Find centroids of active spots, define ROIs around those 

% This is arbitrary, do anything better that comes along
fluor_std_threshold = 5;

active_centroids = zeros(512,512);

for curr_session = 1:length(summed_filenames)
    
    curr_file = [summed_path filesep summed_filenames{curr_session}];
    im_info = imfinfo(curr_file);
    n_frames = length(im_info);
    
    % Load in current day, align, make average
    im = zeros(M,N,n_frames,'uint16');
    for i = 1:n_frames
        curr_frame = imread(curr_file,'tiff',i,'Info',im_info);
        curr_tform = affine2d;
        curr_tform.T = tform_matrix{curr_session};
        im_reg = imwarp(curr_frame,curr_tform,'Outputview',imref2d([M N]));
        im(:,:,i) = im_reg;
    end
    curr_avg = nanmean(double(im),3);
    
    % Make df image, set activity df threshold
    im_processed = double(im) - repmat(double(curr_avg),1,1,size(im,3));
    df_threshold = std(single(im_processed(:)))*fluor_std_threshold;
    
    % Find centroids of active spots
    for i = 1:n_frames
        curr_frame = im_processed(:,:,i);
        
        curr_frame_thresh = curr_frame >= df_threshold;
        
%         % TESTING: this makes it take a lot longer but might be worth it to
%         % separate concurrently active neighboring dendrites
%         % (this didn't seem to help much, so I'm disabling it)
%         % Split thresholded objects by local maxima
%         curr_activity = zeros(M,N);
%         curr_activity(curr_frame_thresh) = curr_frame(curr_frame_thresh);
%         curr_local_max = imextendedmax(curr_activity,df_threshold/10);
%         curr_local_max_borders = watershed(bwdist(curr_local_max)) == 0;
%         curr_frame_thresh(curr_local_max_borders) = 0;
        
        % Eliminiate objects which are too small
        min_px = 3;
        curr_frame_thresh = bwareaopen(curr_frame_thresh,min_px);
        
        stats = regionprops(curr_frame_thresh);
        if length(stats) == 0
            continue
        end
        active_centroids_sub = round(vertcat(stats.Centroid));
        active_centroids_ind = sub2ind([M,N],active_centroids_sub(:,2), ...
            active_centroids_sub(:,1));
        
        active_centroids(active_centroids_ind) = active_centroids(active_centroids_ind)+1;
    end
    
    disp(['Active centroids, session: ' num2str(curr_session)]);
    
end

active_centroids_max = imextendedmax(active_centroids,2);
%active_centroids_max = imregionalmax(active_centroids);

% Store indicies of each ROI centroid
cc = bwconncomp(active_centroids_max);
active_centroid_rois = cc.PixelIdxList;

% % Eliminate active centers that are close to other active centers
% % (this gets rid of good ones and isn't necessary given later checks)
% max_dist = 5;
% close_roi = false(length(active_centroid_rois),1);
% for i = 1:length(active_centroid_rois)
%    temp_map = active_centroids_max;
%    temp_map(active_centroid_rois{i}) = 0;
%    temp_dist = bwdist(temp_map); 
%    close_roi(i) = any(temp_dist(active_centroid_rois{i}) < max_dist);
% end
% active_centroid_rois(close_roi) = [];

%% Get shapes of active objects

% Once centers gotten, go through thresholded images
% and get median shape of objects selected by their centers
roi_active_px = zeros(M,N,length(active_centroid_rois));
roi_active_frames = zeros(length(active_centroid_rois),1);
for curr_session = 1:length(summed_filenames)
    
    curr_file = [summed_path filesep summed_filenames{curr_session}];
    im_info = imfinfo(curr_file);
    n_frames = length(im_info);
    
    % Load in current day, align, make average
    im = zeros(M,N,n_frames,'uint16');
    for i = 1:n_frames
        curr_frame = imread(curr_file,'tiff',i,'Info',im_info);
        curr_tform = affine2d;
        curr_tform.T = tform_matrix{curr_session};
        im_reg = imwarp(curr_frame,curr_tform,'Outputview',imref2d([M N]));
        im(:,:,i) = im_reg;
    end
    curr_avg = nanmean(double(im),3);
    
    % Make df image, set activity df threshold
    im_processed = double(im) - repmat(double(curr_avg),1,1,size(im,3));
    df_threshold = std(single(im_processed(:)))*fluor_std_threshold;
    
    im_thresh = im_processed > df_threshold;
    
    % Find frames where ROI is active (centroid = 1);
    for curr_roi = 1:length(active_centroid_rois)
        [curr_y curr_x] = ind2sub([M,N],active_centroid_rois{curr_roi});
        active_frames = any(cell2mat(arrayfun(@(x) im_thresh(curr_y(x),curr_x(x),:),1:length(curr_y),'uni',false)),2);
        active_px = nansum(im_thresh(:,:,active_frames),3);
        
        roi_active_px(:,:,curr_roi) = roi_active_px(:,:,curr_roi) + active_px;
        roi_active_frames(curr_roi) = roi_active_frames(curr_roi) + sum(active_frames);       
    end
    
    disp(['Active shape, session: ' num2str(curr_session)]);
    
end

% Find the shape for each ROI
frac_active_frames = 0.5;

dilate_size = 1;
dilate_filt = fspecial('disk',dilate_size);

roi_mask = false(M,N,length(active_centroid_rois));
for curr_roi = 1:length(active_centroid_rois)    
    [curr_y curr_x] = ind2sub([M,N],active_centroid_rois{curr_roi});
    curr_mask = bwselect(roi_active_px(:,:,curr_roi) >= roi_active_frames(curr_roi)*frac_active_frames,curr_x,curr_y); 
    
    % If it's no or all pixels, skip it
    if all(curr_mask(:)) || ~any(curr_mask(:));
        continue
    end
    
    % If it's less than 10 pixels, dilate
    if sum(curr_mask(:)) < 10
       curr_mask = conv2(+curr_mask,+dilate_filt,'same') > 0; 
    end
    
    roi_mask(:,:,curr_roi) = curr_mask;     
end

% Sometimes multiple centers are picked up for the same ROI where one is
% inside of the other. Look for large overlaps, delete the larger one.
overlap_frac_thresh = 0.5;

roi_overlap_px = reshape(+roi_mask,[],size(roi_mask,3))'*reshape(+roi_mask,[],size(roi_mask,3));
roi_px = sum(reshape(roi_mask,[],size(roi_mask,3)));
roi_min_px = min(cat(3,repmat(roi_px,size(roi_mask,3),1), ...
    repmat(roi_px,size(roi_mask,3),1)'),[],3);
roi_overlap_frac = tril((roi_overlap_px./roi_min_px),-1) > overlap_frac_thresh;
[roi_1,roi_2] = ind2sub(repmat(size(roi_mask,3),2,1),find(roi_overlap_frac));
remove_overlying_roi = unique([roi_1((roi_px(roi_1) > roi_px(roi_2))); ...
    roi_2((roi_px(roi_2) > roi_px(roi_1)))]);

roi_mask(:,:,remove_overlying_roi) = [];

% Get rid of remaining overlap with buffer zone
dilate_size = 1;
dilate_filt = fspecial('disk',dilate_size);

overlapping_px = conv2(+(sum(roi_mask,3) > 1),+dilate_filt,'same') > 0;
roi_mask(repmat(overlapping_px,1,1,size(roi_mask,3))) = 0;

% Make buffer zone next to immediate neighbors
dilate_size = 1;
dilate_filt = fspecial('disk',dilate_size);

neighboring_px = sum(convn(+roi_mask,dilate_filt,'same') > 0,3) > 1;
roi_mask(repmat(neighboring_px,1,1,size(roi_mask,3))) = 0;

% If ROI was split during this, choose the largest object
roi_bwobj = arrayfun(@(x) bwconncomp(roi_mask(:,:,x)),1:size(roi_mask,3),'uni',false);
split_rois = find(cellfun(@(x) x.NumObjects,roi_bwobj) > 1);
[~,larger_shape] = cellfun(@(x) max(cellfun(@length,x.PixelIdxList)),roi_bwobj(split_rois));
[roi_i,roi_j] = arrayfun(@(x) ind2sub(roi_bwobj{split_rois(x)}.ImageSize, ...
    roi_bwobj{split_rois(x)}.PixelIdxList{larger_shape(x)}(1)),1:length(split_rois));
selected_bw = arrayfun(@(x) bwselect(roi_mask(:,:,split_rois(x)),roi_j(x),roi_i(x)), ...
    1:length(split_rois),'uni',false);

roi_mask(:,:,split_rois) = cat(3,selected_bw{:});

% After this process, get rid of any ROIs which are too small (either in
% total pixel number, or along minor axis i.e. just one row of pixels)
min_px = 5;
roi_px = sum(reshape(roi_mask,[],size(roi_mask,3)),1);
roi_mask(:,:,roi_px < min_px) = [];

% This didn't work because some are a single line but curve
% min_minoraxis = 2;
% roi_minor_axis = cellfun(@(x) x.MinorAxisLength,arrayfun(@(x) ...
%     regionprops(roi_mask(:,:,x),'MinorAxisLength'),1:size(roi_mask,3),'uni',false));
% roi_mask(:,:,roi_minor_axis < min_minoraxis) = [];

% Get rid of any ROIs which do not contain a center of activity
peripheral_roi = ~any(bsxfun(@plus,reshape(roi_mask,[],size(roi_mask,3)),active_centroids_max(:)) > 1,1);
roi_mask(:,:,peripheral_roi) = [];

%% Align ROIs for each session

% Reverse-align ROIs for each session, draw polygon
disp('Creating ROIs (%):')
fprintf('%2d',0)
roi_polygon = cell(n_sessions,1);

roi_dilate = strel('rectangle',[3 3]);

for curr_session = 1:n_sessions;
    
    % make reverse registration transform matrix
    curr_tform = affine2d;
    inverse_tform = inv(tform_matrix{curr_session});
    % inverse is not perfect: force last column
    inverse_tform(1:2,3) = 0;
    inverse_tform(3,3) = 1;
    curr_tform.T = inverse_tform;
    
    roi_polygon{curr_session} = cell(1,size(roi_mask,3));
    for curr_roi = 1:size(roi_mask,3);
        
        binaryCell = roi_mask(:,:,curr_roi);
        
        % reverse register cell to fit session
        binaryCell_reg = +(imwarp(binaryCell,curr_tform,'Outputview',imref2d([M N])) > 0);
        
        % If registering pushes ROI off edge, ignore now and delete later
        if sum(binaryCell_reg(:)) == 0
            continue
        end
        
        % dilate the ROI to trace around the active pixels (make more
        % accurate polygon for poly2mask conversion)
        binaryCell_reg_dilated = imdilate(binaryCell_reg,roi_dilate);
        
        first_nonzero = find(binaryCell_reg_dilated > 0,1);
        [y_nonzero x_nonzero] = ind2sub([M N],first_nonzero);
        roi_edge = bwtraceboundary(binaryCell_reg_dilated,[y_nonzero x_nonzero],'N'); % find boundary, in order
        
%         % reduce the number of points to make border (from file exchange)
%         if size(roi_edge,1) > 10
%             num_verticies = 10;
%             roi_edge = reduce_poly(roi_edge',num_verticies)';
%         end
              
        % make x,y coordinates
        roi_polygon{curr_session}{curr_roi} = fliplr(roi_edge);
        
    end
    fprintf('%c%c%2d',8,8,round(100*curr_session/n_sessions));
end


% Delete any ROIs that are not present / off edge after aligning
offscreen_rois = any(cell2mat(cellfun(@(x) cellfun(@(x) isempty(x), x'), roi_polygon','uni',false)),2);
roi_polygon_use = cellfun(@(x) x(~offscreen_rois),roi_polygon,'uni',false);

%% Check for and eliminate ROIs that overlap with edge artifacts

% Give leeway for summed borders (more for vertical because of edge artifacts)
vertical_leeway = 20;
horizontal_leeway = 10;

vertical_filt = ones(1,vertical_leeway);
horizontal_filt = ones(horizontal_leeway,1);

% Set all direct edges to 1 have general edge leeway even if no artifact
summed_borders(:,1,:) = 1;
summed_borders(:,end,:) = 1;
summed_borders(1,:,:) = 1;
summed_borders(end,:,:) = 1;

summed_borders_leeway = cell2mat(arrayfun(@(x) conv2(summed_borders(:,:,x),vertical_filt,'same') > 0 | ...
    conv2(summed_borders(:,:,x),horizontal_filt,'same') > 0,permute(1:n_sessions,[1 3 2]),'uni',false));

edge_rois = cell(1,n_sessions);
for curr_session = 1:n_sessions;
    edge_rois{curr_session} = false(length(roi_polygon_use{curr_session}),1);
    for curr_roi = 1:length(roi_polygon_use{curr_session})
        curr_mask = poly2mask(roi_polygon_use{curr_session}{curr_roi}(:,1)',roi_polygon_use{curr_session}{curr_roi}(:,2)',M,N);
        edge_rois{curr_session}(curr_roi) = any(reshape(curr_mask & summed_borders_leeway(:,:,curr_session),[],1));
    end  
end
edge_rois = any(horzcat(edge_rois{:}),2);

% Delete any ROIs that are too close to edge artifacts
roi_polygon_final = cellfun(@(x) x(~edge_rois),roi_polygon_use,'uni',false);




%% Display final ROIs over aligned max images
if n_sessions > 1
    gui_fig = AP_image_scroll(summed_max_reg);
else
    gui_fig = figure;imagesc(summed_max_reg);colormap(gray);
end

% Draw ROIs on the aligned images for inspection
figure(gui_fig)
for curr_roi = 1:length(roi_polygon_final{1})
    line(roi_polygon_final{1}{curr_roi}(:,1),roi_polygon_final{1}{curr_roi}(:,2));
end


%% Save the ROI polygons
if save_flag
    save_path = [summed_path animal '_batch_thresh_roi'];
    if ~exist(save_path,'dir')
        mkdir(save_path)
    end
    for curr_session = 1:n_sessions;
        clear polygon
        polygon.ROI = roi_polygon_final{curr_session};
        
        curr_save_filename = [save_path filesep summed_filenames{curr_session}(1:end-4) '.roi'];
        save(curr_save_filename,'polygon');
    end
end

disp('Done.')





