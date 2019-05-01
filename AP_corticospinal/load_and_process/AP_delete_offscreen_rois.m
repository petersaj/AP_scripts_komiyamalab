function AP_delete_offscreen_rois(animal)
% AP_delete_offscreen_rois(animal)
%
% Get maximum motion correction artifacts on borders, delete those ROIs
% from ROI files and from anaylsis files (this is now built into dendrite
% ROI and shouldn't be used for new ROIs)

%% Find ROIs too close to edges

data_path = ['C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\' animal];
summed_files = dir([data_path filesep '*.tif']);
summed_filenames = cellfun(@(x) [data_path filesep x],{summed_files.name},'uni',false);

n_sessions = length(summed_filenames);
M = 512;
N = 512;

% Load all summed movies, look for motion correction artifacts (complete
% line of zeros)
borders = zeros(M,N,n_sessions);
disp('Loading summed movies, checking for borders...')
for curr_session = 1:n_sessions;
    
    im_info = imfinfo(summed_filenames{curr_session});
    n_frames = length(im_info);
        
    for i = 1:n_frames
        curr_im = imread(summed_filenames{curr_session},'tiff',i,'Info',im_info);
        borders(:,all(curr_im == 0,1),curr_session) = 1;
        borders(all(curr_im == 0,2),:,curr_session) = 1;
    end
end

% Give leeway for borders (more for vertical because of edge artifacts)
vertical_leeway = 20;
horizontal_leeway = 10;

vertical_filt = ones(1,vertical_leeway);
horizontal_filt = ones(horizontal_leeway,1);

% Set all direct edges to 1 have general edge leeway even if no artifact
borders(:,1,:) = 1;
borders(:,end,:) = 1;
borders(1,:,:) = 1;
borders(end,:,:) = 1;

borders_leeway = cell2mat(arrayfun(@(x) conv2(borders(:,:,x),vertical_filt,'same') > 0 | ...
    conv2(borders(:,:,x),horizontal_filt,'same') > 0,permute(1:n_sessions,[1 3 2]),'uni',false));
disp('Marking offscreen ROIs...')
% Mark ROIs that ever fall within a motion correction artifact zone
roi_path = [data_path filesep animal '_batch_thresh_roi'];
roi_files = dir([roi_path filesep '*.roi']);
roi_filenames = cellfun(@(x) [roi_path filesep x],{roi_files.name},'uni',false);

if ~all(cellfun(@(x,y) strcmp(x(1:end-3),y(1:end-3)),{summed_files.name},{roi_files.name}))
    error('Summed filenames do not match ROI filenames')
end

bad_rois = cell(1,n_sessions);
for curr_session = 1:n_sessions;
    curr_rois = load(roi_filenames{curr_session},'-MAT');
    bad_rois{curr_session} = false(length(curr_rois.polygon.ROI),1);
    for curr_roi = 1:length(curr_rois.polygon.ROI)
        curr_mask = poly2mask(curr_rois.polygon.ROI{curr_roi}(:,1)',curr_rois.polygon.ROI{curr_roi}(:,2)',M,N);
        bad_rois{curr_session}(curr_roi) = any(reshape(curr_mask & borders_leeway(:,:,curr_session),[],1));
    end
end

edge_rois = any(horzcat(bad_rois{:}),2);


%% Delete ROIs from ROI files and analysis files
% Must be local

% Delete edge rois from ROI files
disp('Deleting edge rois from ROI files...')
for curr_session = 1:n_sessions
    clear polygon
    load(roi_filenames{curr_session},'-MAT');
    polygon.ROI(edge_rois) = [];
    save(roi_filenames{curr_session},'polygon');
    disp(curr_session);
end

% Delete edge rois from analysis files
disp('Deleting edge rois from analysis files...')
analysis_path = ['C:\Users\Andy\Documents\KomiyamaLab\data\analysis_files\' animal];
for curr_session = 1:n_sessions    
    clear im bhv    
    
    curr_analysis_filename = [analysis_path filesep ...
        summed_files(curr_session).name(1:end-4) '_analysis.mat'];
    
    load(curr_analysis_filename);
    
    im.roi_trace_df(edge_rois,:) = [];
    im.roi_trace(edge_rois,:) = [];
    im.roi_trace_bg(edge_rois,:) = [];
    
    save(curr_analysis_filename,'im','bhv');  
    disp(curr_session);
end

disp(['Finished ' animal]);





    