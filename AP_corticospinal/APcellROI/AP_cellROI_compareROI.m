function AP_cellROI_compareROI
% This GUI is to manually check and tweak all ROIs across all sessions

%% Choose number of PCs to save and the downsampling factor

x = inputdlg({'Number of PCs to save','Average around activity (0/1)', ...
    'Image downsample factor', 'Cells to extract','Interpolation subtraction'},'',1,{'0','1','1','all','0'});

% PCs to save
save_pcs = str2num(x{1});
% Average activity
avg_active = logical(str2num(x{2}));
% Downsample ratio
downsample_px = str2num(x{3});
% Cells to extract
use_rois = x{4};
% Interpolation subtraction (point-source ROIs only!)
interp_sub = logical(str2num(x{5}));


%% Select tiff images and corresponding ROI files

[tiff_files, tiff_path] = uigetfile('*.tif','Choose tiff files','multiselect','on');
if ~iscell (tiff_files)
    tiff_files = {tiff_files};
end
tiff_files = sort(tiff_files);
[roi_files, roi_path] = uigetfile('*.roi','Choose ROI files','multiselect','on');
if ~iscell (roi_files)
    roi_files = {roi_files};
end
roi_files = sort(roi_files);

% Error out if the number doesn't match up, or if each tiff file doesn't
% have a corresponding ROI file containing it's filename
if length(tiff_files) ~= length(roi_files)
    error('Number of TIFF files does not match number of ROI files');
elseif ~all(cellfun(@(x) any(cellfun(@(y) any(strfind(x(1:end-4),y(1:end-4))),tiff_files)),roi_files));
    error('TIFF file names do not have corresponding ROI file names')
end

%% Prepare data: loop through each session, pull out ROI snippets

% Get number of ROIs from the first file
curr_polygon_temp = load([roi_path roi_files{1}],'-mat');
num_rois = length(curr_polygon_temp.polygon.ROI);
if strmatch(use_rois,'all')
    use_rois = 1:num_rois;
else
    use_rois = str2num(use_rois);
end
% Calculate size of image to extract (2 * median largest dimension of ROI)
roi_maxlength = cellfun(@(x) max(pdist(x)),curr_polygon_temp.polygon.ROI);
%cell_border = 2*median(round(roi_maxlength));
% AP 150509 - changed to be standard value
cell_border = 50;
clear curr_polygon_temp

% Initialize structure for extracing ROIs
rois = struct('number',num2cell(use_rois)','polygon',cell(length(use_rois),1), ...
    'center',cell(length(use_rois),1),'im',cell(length(use_rois),1));

% Save user settings
user_settings = struct('border',cell_border,'save_pcs',save_pcs,'downsample_px',downsample_px);

load_waitbar = waitbar(0,'Loading summed movie');
for curr_session = 1:length(tiff_files)
    waitbar((curr_session-1)/length(tiff_files),load_waitbar,'Loading summed movie');
    imageinfo = imfinfo([tiff_path tiff_files{curr_session}],'tiff');
    M = imageinfo(1).Width;
    N = imageinfo(1).Height;
    num_frames = length(imageinfo);
    im = repmat({zeros(N,M,'uint16')},num_frames,1);
    for loadframe = 1:num_frames
        curr_frame = imread([tiff_path tiff_files{curr_session}],'tiff',loadframe,'Info',imageinfo);
        im{loadframe} = uint16(curr_frame);
    end
            
    % Load current ROI set
    curr_polygon_temp = load([roi_path roi_files{curr_session}],'-mat');
    curr_polygons = curr_polygon_temp.polygon.ROI;
    
    % Interpolation subtraction
    if interp_sub      
        im_interp_sub = im;
        % Make temporary normal-shape ROI masks
        curr_masks = sum(cell2mat(permute(cellfun(@(x) poly2mask(x(1:downsample_px:end,1), ...
            x(1:downsample_px:end,2),512,512),curr_polygons,'uni',false),[1 3 2])),3) > 0;
        
        % Average and interp avg
        im_avg = nanmean(double(cat(3,im{:})),3);
        
        im_avg_roinan = im_avg;
        im_avg_roinan(curr_masks) = NaN;
        im_avg_interpsub = inpaint_nans(im_avg_roinan,4);
        
        % Make ROIs = NaN, interpolate, subtract
        waitbar((curr_session-1)/length(tiff_files),load_waitbar,'Subtracting interpolated background');
        for i = 1:length(im)
            curr_nan_frame = double(im{i});
            curr_nan_frame(curr_masks) = NaN;
            curr_interp_frame = inpaint_nans(curr_nan_frame);
            im_interp_sub{i} = reshape( ...
                im{i} - min(cat(3,(uint16(curr_interp_frame) - uint16(im_avg_interpsub)), ...
                (uint16(double(im{i})) - uint16(im_avg))),[],3),[],1);
        end 
        im_interp_sub = horzcat(im_interp_sub{:});
    end
    
    % Downsample image (try to save space and not create new variable)
    disk_avg = fspecial('disk',downsample_px);
    for i = 1:length(im);
        curr_conv = conv2(double(im{i}),disk_avg,'same');
        im{i} = reshape(curr_conv(1:downsample_px:end,1:downsample_px:end),[],1);
    end
    im = horzcat(im{:});
    N_ds = ceil(N/downsample_px);
    M_ds = ceil(M/downsample_px);
    cell_border_ds = ceil(cell_border/downsample_px);
    
    % Get cell extraction pixels
    cell_x = round(cell_border_ds/2);
    cell_y = round(cell_border_ds*(N_ds/M_ds)/2);
    
    % Make masks from polygon
    curr_masks = cellfun(@(x) reshape(poly2mask(x(1:downsample_px:end,1), ...
        x(1:downsample_px:end,2),N_ds,M_ds),[],1),curr_polygons,'uni',false);
    
    waitbar((curr_session-1)/length(tiff_files),load_waitbar,'Extracing ROI images');
    for curr_roi_idx = 1:length(use_rois)
        curr_roi = use_rois(curr_roi_idx);
        % Initialize image variable
        rois(curr_roi_idx).im{curr_session} = nan(1+2*cell_y,1+2*cell_x,1);
        
        % Find center of ROIs, save
        [~,cx,cy] = polycenter(curr_polygons{curr_roi}(:,1), ...
            curr_polygons{curr_roi}(:,2));
        cx_ds = round(cx/downsample_px);
        cy_ds = round(cy/downsample_px);
        
        rois(curr_roi_idx).polygon{curr_session} = ...
            [curr_polygons{curr_roi}(:,1)/downsample_px-cx_ds+cell_x+1 ...
            curr_polygons{curr_roi}(:,2)/downsample_px-cy_ds+cell_y+1];
        
        % Store downsampled ROI centers (relative to where they'll be plotted)
        rois(curr_roi_idx).center{curr_session} = [cx_ds-cell_x-1 cy_ds-cell_y-1];
        
        % Define range of image to look at, pull out
        cell_y_range = cy_ds-cell_y:cy_ds+cell_y;
        cell_x_range = cx_ds-cell_x:cx_ds+cell_x;
        
        cell_y_range(cell_y_range < 1 | cell_y_range > N_ds) = [];
        cell_x_range(cell_x_range < 1 | cell_x_range > M_ds) = [];
        
        [cell_x_mesh,cell_y_mesh] = meshgrid(cell_x_range,cell_y_range);
        cell_ind = sub2ind([N_ds M_ds],cell_y_mesh(:),cell_x_mesh(:));
        cell_im = im(cell_ind,:);
        
        % Get image range of ROI
        cell_y_range = cy_ds-cell_y:cy_ds+cell_y;
        cell_x_range = cx_ds-cell_x:cx_ds+cell_x;
        cell_y_use = cell_y_range >= 1 & cell_y_range <= N_ds;
        cell_x_use = cell_x_range >= 1 & cell_x_range <= M_ds;
        
        % Save mean(first image)/max(second image)
        cell_im_mean = nan(1+cell_y*2,1+cell_x*2);
        cell_im_mean(cell_y_use,cell_x_use,:) = ...
            reshape(nanmean(cell_im,2),sum(cell_y_use),sum(cell_x_use));
    
        cell_im_max = nan(1+cell_y*2,1+cell_x*2);
        cell_im_max(cell_y_use,cell_x_use,:) = ...
            reshape(nanmax(cell_im,[],2),sum(cell_y_use),sum(cell_x_use));
        
        curr_mask = poly2mask(rois(curr_roi_idx).polygon{curr_session}(:,1), ...
            rois(curr_roi_idx).polygon{curr_session}(:,2), ...
            sum(cell_y_use),sum(cell_x_use));
        mean_max_px = max(cell_im_mean(curr_mask)) * 1.5;
        max_max_px = max(cell_im_max(curr_mask)) * 1.5;
        
        rois(curr_roi_idx).im{curr_session}(:,:,1) = mat2gray(cell_im_mean,[0 mean_max_px]);
        rois(curr_roi_idx).im{curr_session}(:,:,2) = mat2gray(cell_im_max,[0 max_max_px]);        
        
        % Do PCA on pixels within dilated image mask, apply those time
        % weights to all neighborhood pixels
        if save_pcs > 0
            se = strel('disk',5);
            dilated_mask = imdilate(reshape(curr_masks{curr_roi},N_ds,M_ds),se);
            [coeff,score,latent] = princomp(zscore(double(im(dilated_mask(:),:)),[],2));
            neighborhood_pc_scores = zscore(double(cell_im),[],2)*coeff(:,1:save_pcs);
            % Skew the scores in the same direction
            neighborhood_pc_scores = bsxfun(@times,neighborhood_pc_scores, ...
                sign(skewness(neighborhood_pc_scores)));
            
            % If ICA is wanted (the way I'm doing it, looks worse than PCA):
            %[icasig,A,W] = fastica(zscore(double(im(dilated_mask(:),:)),[],2)','numOfIc',save_pcs);
            %neighborhood_pc_scores = transpose(W*zscore(double(cell_im),[],2)');
            
            % Save the PCs in square form
            cell_pc_score = nan(1+cell_y*2,1+cell_x*2,save_pcs);
            cell_pc_score(cell_y_use,cell_x_use,:) = ...
                reshape(neighborhood_pc_scores,sum(cell_y_use),sum(cell_x_use),save_pcs);
            
            % Save PCs (3:end = PC scores)
            rois(curr_roi_idx).im{curr_session} = cat(3,rois(curr_roi_idx).im{curr_session},cell_pc_score);
        end
        
        % If average activity selected: save average around active portions
        if avg_active
            
            pre_frames = 3;
            updown_frames = 3;
            post_frames = 3;
            
            if interp_sub
                curr_roi_trace = nanmean(im_interp_sub(curr_masks{curr_roi},:));
            else
                curr_roi_trace = nanmean(im(curr_masks{curr_roi},:));               
            end
            
            % Get envelope of first derivative
            vel_env = abs(hilbert(diff(curr_roi_trace)));
            % Median of envelope is approx the noise level
            vel_thresh = 2*median(vel_env);
            vel_overthresh = vel_env > vel_thresh;
            activity_onsets = find(diff([0 vel_overthresh 0]) > 0);
            activity_offsets = find(diff([0 vel_overthresh 0]) < 0);
            
            out_of_bounds = cellfun(@(x,y) (x-pre_frames < 1) | ...
                (y + post_frames > length(curr_roi_trace)),num2cell(activity_onsets), ...
                num2cell(activity_offsets));
            activity_duration = activity_offsets - activity_onsets;
            
            use_act = ~out_of_bounds & activity_duration >= updown_frames;
            
            activity_frames = cellfun(@(x,y) [x-pre_frames:x+updown_frames-1 ...
                y-updown_frames+1:y+post_frames], num2cell(activity_onsets(use_act)), ...
                num2cell(activity_offsets(use_act)),'uni',false);
            
            curr_act_im = reshape(nanmean(cell2mat(permute(cellfun(@(x) cell_im(:,x), ...
                activity_frames,'uni',false),[1 3 2])),3),sum(cell_y_use),sum(cell_x_use),[]);
            
            % Normalize the caxis based on the ROI
            curr_mask = poly2mask(rois(curr_roi_idx).polygon{curr_session}(:,1), ...
                rois(curr_roi_idx).polygon{curr_session}(:,2), ...
                sum(cell_y_use),sum(cell_x_use));
            max_px = max(curr_act_im(repmat(curr_mask,1,1, ...
                size(curr_act_im,3)))) * 1.5;
            
            curr_act_im = mat2gray(curr_act_im,[0 max_px]);
            
            total_frames = pre_frames + 2*updown_frames + post_frames;
            
            if ~isempty(curr_act_im)
                rois(curr_roi_idx).im{curr_session}(cell_y_use,cell_x_use,end+1:end+size(curr_act_im,3)) = curr_act_im;
            else
                rois(curr_roi_idx).im{curr_session}(:,:,end+1:end+size(curr_act_im,3)) = 0;
            end

        end
        
     
    end
end
close(load_waitbar);


%% Create GUI for manually checking ROIs

handles = struct;

% Create figure
roi_compare_gui = figure;

% Set the axes for plots
% Set the X by Y subplots according to 1.6 ratio monitor
subplots_y = round(sqrt(length(rois(1).im)/1.6));
subplots_x = round(length(rois(1).im)/subplots_y);

handles.im_axes = tight_subplot(subplots_y,subplots_x,0,0);
% Plot the first image of first ROI for each session
for curr_session = 1:length(rois(1).im)
    handles.im(curr_session) = ...
        imagesc(rois(1).im{curr_session}(:,:,1),'parent',handles.im_axes(curr_session));
    axis(handles.im_axes,'off');
    colormap(gray);
end

% Plot ROI
for curr_session = 1:length(rois(1).im)
    handles.roi(curr_session) = ...
        impoly(handles.im_axes(curr_session),rois(1).polygon{1,curr_session});
    % This was intended to keep all polygons same shape
    %addNewPositionCallback(handles.roi(curr_session),@(y,z) ...
    %        update_roi_verticies(handles.roi(curr_session),roi_compare_gui));
end

% Save ROI filenames
handles.roi_files = roi_files;

% Set the current ROI/im indicies
handles.curr_im = 1;
handles.curr_roi = 1;

% Package variables for passing to functions
gui_data.handles = handles;
gui_data.rois = rois;
gui_data.user_settings = user_settings;

% Set up scoring (for keeping or rejecting)
gui_data.score = zeros(size(rois));

% Store guidata
guidata(roi_compare_gui, gui_data);

% Set functions for mouse wheel / button press and pass necessary handles
set(roi_compare_gui,'KeyPressFcn', {@change_event, roi_compare_gui}); 
set(roi_compare_gui,'WindowScrollWheelFcn',{@slide_movie, roi_compare_gui});
set(roi_compare_gui,'Name',num2str(gui_data.rois(1).number));

end


%% Slider listener function

function slide_movie(currentObject, eventdata, roi_compare_gui)
% Get guidata
gui_data = guidata(roi_compare_gui);

% Get current event/frame
curr_roi = gui_data.handles.curr_roi;
old_im = gui_data.handles.curr_im;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
new_im = old_im + mouse_wheel_count;
for curr_session = 1:length(gui_data.rois(curr_roi).im);
    num_ims = size(gui_data.rois(curr_roi).im{curr_session},3);
    if new_im > num_ims
        new_im = num_ims;
    elseif new_im < 1
        new_im = 1;
    end
    
    % Update the image and trace marker
    set(gui_data.handles.im(curr_session),'Cdata', ...
        gui_data.rois(curr_roi).im{curr_session}(:,:,new_im));
    
    % Update the caxis based on the ROI
    temp_im = gui_data.rois(curr_roi).im{curr_session}(:,:,new_im);
    curr_mask = poly2mask(gui_data.rois(curr_roi).polygon{curr_session}(:,1), ...
        gui_data.rois(curr_roi).polygon{curr_session}(:,2), ...
        size(gui_data.rois(curr_roi).im{curr_session}(:,:,new_im),1), ...
        size(gui_data.rois(curr_roi).im{curr_session}(:,:,new_im),2));
    %caxis(gui_data.handles.im_axes(curr_session),[0 max(temp_im(curr_mask))]);
end

% Update guidata
gui_data.handles.curr_im = new_im;
guidata(roi_compare_gui, gui_data);
end


function change_event(currentObject, eventdata, roi_compare_gui)
% Get guidata
gui_data = guidata(roi_compare_gui);

curr_roi = gui_data.handles.curr_roi;

% Define the keypress events
change_roi = false;
score_event = false;
show_roi = false;
save_rois = false;
switch eventdata.Key
    case 'leftarrow'
        if curr_roi ~= 1
            new_roi = curr_roi - 1;
            change_roi = true;
            gui_data.handles.curr_roi = new_roi;
        end
    case 'rightarrow'
        if curr_roi ~= length(gui_data.rois)
            new_roi = curr_roi + 1;
            change_roi = true;
            gui_data.handles.curr_roi = new_roi;
        end
    case 'uparrow'
        if gui_data.score(curr_roi) ~= 1
            new_score = gui_data.score(curr_roi) + 1;
            gui_data.score(curr_roi) = new_score;
            score_event = true;
        end
    case 'downarrow'
        if gui_data.score(curr_roi) ~= -1
            new_score = gui_data.score(curr_roi) - 1;
            gui_data.score(curr_roi) = new_score;
            score_event = true;
        end
    case 'r'
        show_roi = true;
    case 's'
        save_rois = true;
end


% Change events if left/right arrow is pressed
if change_roi
    
    % Update background color    
    switch gui_data.score(new_roi)
        case -1
            set(roi_compare_gui,'Color',[0.8 0 0])
        case 0 
            set(roi_compare_gui,'Color',[0.8 0.8 0.8]);
        case 1
            set(roi_compare_gui,'Color',[0 0.8 0])
    end
    
    % Update image / ROI
    gui_data.handles.curr_im = 1;
    for curr_session = 1:length(gui_data.rois(curr_roi).im);
        set(gui_data.handles.im(curr_session),'Cdata', ...
            gui_data.rois(new_roi).im{curr_session}(:,:,gui_data.handles.curr_im));
        if all(size(gui_data.rois(curr_roi).im(:,:,1),1) ~= ...
                size(gui_data.rois(new_roi).im(:,:,1),1));
            % If new image is different size, replace data and fix x/y limits
            new_im_size = size(gui_data.rois(new_roi).im(:,:,1));
            set(gui_data.handles.im_axes,'YLim',[0 new_im_size(1)]+0.5, ...
                'XLim',[0 new_im_size(2)]+0.5);
        end
        % Reposition polygons with new verticies
         setPosition(gui_data.handles.roi(curr_session), ...
             gui_data.rois(new_roi).polygon{curr_session});
         
         % Update the caxis based on the ROI
         temp_im = gui_data.rois(new_roi).im{curr_session}(:,:,gui_data.handles.curr_im);
         curr_mask = poly2mask(gui_data.rois(new_roi).polygon{curr_session}(:,1), ...
             gui_data.rois(new_roi).polygon{curr_session}(:,2), ...
             size(gui_data.rois(new_roi).im{curr_session}(:,:,gui_data.handles.curr_im),1), ...
             size(gui_data.rois(new_roi).im{curr_session}(:,:,gui_data.handles.curr_im),2));
         caxis(gui_data.handles.im_axes(curr_session),[0 max(temp_im(curr_mask))]);
         
    end
    
    
    
    set(roi_compare_gui,'Name',num2str(gui_data.rois(new_roi).number));
end


% Score event if up/down arrow is pressed
if score_event
    switch new_score
        case -1
            set(roi_compare_gui,'Color',[0.8 0 0])
        case 0 
            set(roi_compare_gui,'Color',[0.8 0.8 0.8]);
        case 1
            set(roi_compare_gui,'Color',[0 0.8 0])
            % If the score is set to 1: save the new ROI positions
            for curr_session = 1:length(gui_data.rois(curr_roi).polygon);
                gui_data.rois(curr_roi).polygon{curr_session} = ...
                    getPosition(gui_data.handles.roi(curr_session));
            end
    end    
end

if show_roi
    % Show/hide ROI
    for curr_session = 1:length(gui_data.rois(curr_roi).polygon);
        curr_visibility = get(gui_data.handles.roi(curr_session),'Visible');
        switch curr_visibility
            case 'on'
                new_visibility = 'off';
            case 'off'
                new_visibility = 'on';
        end
        set(gui_data.handles.roi(curr_session),'Visible',new_visibility);
    end
end

if save_rois
    % Save ROIs in the same format as saved by cellROI
    save_path = uigetdir([],'Directory to save fixed ROIs');   
    % Only use rois without a negative score
    delete_rois  = [gui_data.score] == -1;
    for curr_session = 1:length(gui_data.handles.roi_files)
        curr_rois = arrayfun(@(x) gui_data.rois(x).polygon{curr_session},1:length(gui_data.rois),'uni',false);
        curr_centers = arrayfun(@(x) gui_data.rois(x).center{curr_session},1:length(gui_data.rois),'uni',false);
        % Restore ratio and center for each ROI
        clear polygon;
        restored_rois = cellfun(@(poly,center) ...
            gui_data.user_settings.downsample_px*(poly + ...
            repmat(center,size(poly,1),1)),curr_rois, ...
            curr_centers,'uni',false);
        polygon.ROI = restored_rois(~delete_rois);

        curr_savefile = [save_path filesep gui_data.handles.roi_files{curr_session}(1:end-4) '_fixed.roi'];
        save(curr_savefile,'polygon');
    end
    disp([datestr(now) ': Saved fixed ROIs']);
end



% Update guidata
guidata(roi_compare_gui, gui_data);
end



%% ROI update function: if a vertex is changed in one, change in all

% It turns out this is difficult to implement: 
% Could try using line/poly selection system like cellROI, switching ROIs
% might be faster then too if line objects, but selection wouldn't be as
% fast on the fly changes so probably not as good (also, nothing wrong with
% having different polygons if clearly the same cell - the SNR would
% change, but that could change between sessions with window/laser anyway

function update_roi_verticies(curr_roi,roi_compare_gui)
%gui_data = guidata(roi_compare_gui);
% 
% % Get the handle index of the currently moved ROI
% curr_roi_idx = find(gui_data.handles.roi == curr_roi);
% 
% % Get relative vertex positions of current ROI
% verticies_diff = [0 0 ; diff(verticies,[],1)];
% keyboard
% % Set the verticies of all other ROIs the same relative to one another
% for curr_session = setdiff(1:length(gui_data.handles.roi),curr_roi_idx)
%     setPosition(gui_data.handles.roi(curr_session), ...
%         getPosition(gui_data.handles.roi(curr_session)) + verticies_diff);    
% end
% 
% disp('universal update verticies');

end








