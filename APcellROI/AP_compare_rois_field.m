function AP_compare_rois_field
% This is to do large scale, blunt quality control
% Align maximum images across days and overlay ROIs, mark bad ROIs

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

n_sessions = length(tiff_files);
imageinfo = imfinfo([tiff_path tiff_files{1}],'tiff');
M = imageinfo(1).Width;
N = imageinfo(1).Height;

% Load first ROI file as template
template_rois = load([roi_path roi_files{1}],'-mat');

% Load ROI files, create masks
disp('Loading ROIs...');
roi_all_mask = nan(M,N,n_sessions);
for curr_session = 1:length(tiff_files)
    
    curr_polygon_temp = load([roi_path roi_files{curr_session}],'-mat');
    curr_polygons = curr_polygon_temp.polygon.ROI;
    
    cellMask = nan(M,N,length(curr_polygons));
    for curr_roi = 1:length(curr_polygons);
        cellMask(:,:,curr_roi) = ...
            poly2mask(curr_polygons{curr_roi}(:,1)',curr_polygons{curr_roi}(:,2)',N,M);
    end
    
    roi_all_mask(:,:,curr_session) = any(cellMask,3);
end

% Find affine transform matrix based on ROI masks
disp('Aligning ROI masks...');
roi_all_mask_aligned = roi_all_mask;
for curr_session = 2:n_sessions
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength = 0.02;
    optimizer.GradientMagnitudeTolerance = 1e-5;
    optimizer.MaximumIterations = 300;
    
    tformEstimate = imregtform(roi_all_mask(:,:,curr_session),roi_all_mask(:,:,1),'affine',optimizer,metric);
    
    mask_reg = imwarp(roi_all_mask(:,:,curr_session),tformEstimate,'Outputview',imref2d([M N]));
    
    tform_matrix{curr_session} = tformEstimate.T;
    roi_all_mask_aligned(:,:,curr_session) = mask_reg;
    
end

% Loop through sessions, get max/mean projection
disp('Loading images...');
summed_max = double(nan(M,N,n_sessions));
summed_mean = double(nan(M,N,n_sessions));
for curr_session = 1:n_sessions
    
    % Get max projection of day
    curr_file = [tiff_path filesep tiff_files{curr_session}];
    
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
        
    clear im
end

% Align max/mean projections
disp('Aligning images...');
summed_max_aligned = double(nan(M,N,n_sessions));
summed_mean_aligned = double(nan(M,N,n_sessions));

summed_max_aligned(:,:,1) = summed_max(:,:,1);
summed_mean_aligned(:,:,1) = summed_mean(:,:,1);
for curr_session = 2:n_sessions
    
    curr_tform = affine2d;
    curr_tform.T = tform_matrix{curr_session};
    
    summed_max_aligned(:,:,curr_session) = ...
        imwarp(summed_max(:,:,curr_session),curr_tform,'Outputview',imref2d([M N]));
    
    summed_mean_aligned(:,:,curr_session) = ...
        imwarp(summed_mean(:,:,curr_session),curr_tform,'Outputview',imref2d([M N]));
    
end

% Display aligned images
gui_fig = AP_image_scroll(summed_max_aligned);
axis square;

roi_handles = nan(length(template_rois.polygon.ROI),1);
for curr_roi = 1:length(template_rois.polygon.ROI)
    roi_handles(curr_roi) = line(template_rois.polygon.ROI{curr_roi}(:,1), ...
        template_rois.polygon.ROI{curr_roi}(:,2),'color','g');
    set(roi_handles(curr_roi), 'ButtonDownFcn', {@selectROI, gui_fig});
end

% Retain gui data from image scroll
gui_data = guidata(gui_fig);

% Add ROI gui data
gui_data.roi_handles = roi_handles;
gui_data.bad_rois = false(length(template_rois.polygon.ROI),1);

% Set up saving feature ('s' key)
set(gui_fig,'KeyPressFcn', {@gui_fig_keypress, gui_fig}); 

% Update guidata
guidata(gui_fig, gui_data);

disp('Set good/bad ROIs');
end


function selectROI(hObject, eventdata, gui_fig)
% Execute when clicking on ROI polygon

gui_data = guidata(gui_fig);

% Current ROI status
roi_selected = find([gui_data.roi_handles(:)] == hObject);
roi_old_bad_status = gui_data.bad_rois(roi_selected);
roi_new_bad_status = ~roi_old_bad_status;

% Recolor ROI by new bad status
switch roi_new_bad_status
    case true
        set(gui_data.roi_handles(roi_selected),'color','r');
    case false
        set(gui_data.roi_handles(roi_selected),'color','g');
end

% Flip current good/bad status of selected ROI
gui_data.bad_rois(roi_selected) = roi_new_bad_status;

% Update guidata
guidata(gui_fig, gui_data);

end

function gui_fig_keypress(currentObject, eventdata, gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    case 's';
        % Save bad ROIs
        [save_filename,save_path] = uiputfile([pwd filesep 'bad_rois.mat'],'Save bad ROIs');
        
        bad_rois = gui_data.bad_rois;
        save([save_path filesep save_filename],'bad_rois');
        
        disp([datestr(now) ': Saved bad ROIs']);
        
    case 'r'
        % Toggle ROI visibility
        curr_visible = strcmp(get(gui_data.roi_handles(1),'Visible'),'on');
        if curr_visible
            set(gui_data.roi_handles,'Visible','off')
        else
            set(gui_data.roi_handles,'Visible','on')
        end        
        
end

end












