% Testing for event detection
% (requires data loaded in through toolbox)


%% Apply event detection to chosen animal/session
animal = 'JL075';
curr_session = 4;

curr_animal = find(cellfun(@(x) strcmp(animal,x),{mice.name}));
scope = mice(curr_animal).scope;
curr_concat_data = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});
% This is to test different thresholds and manually scoring data?
curr_concat_data_thresh = AP_caEvents_threshtest(curr_concat_data,4);

% This was to test out the event detection used by bscope
%curr_concat_data_thresh = AP_caEvents(curr_concat_data, ...
%    ~data_all(curr_animal).labels.gad,data_all(curr_animal).labels.gad);


% switch scope
%     case 1
%         % MOM event detection
%         curr_concat_data_thresh = AP_caEvents_momscope(curr_concat_data, ...
%             ~data_all(curr_animal).labels.gad,data_all(curr_animal).labels.gad);
%     case 2
%         % Bscope event detection
%         curr_concat_data_thresh = AP_caEvents(curr_concat_data, ...
%             ~data_all(curr_animal).labels.gad,data_all(curr_animal).labels.gad);
% end


%% Get summed movie to look at images of cells

% Find the summed movie with the selected animal/session, load  movie
data_path = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_roi'];
data_dir = dir(data_path);
data_names = {data_dir.name};
days = cellfun(@(x) x(end-5:end),data_names(3:end),'UniformOutput',false);

curr_abs_session = data_all(curr_animal).session(curr_session);

im_dir = [data_path filesep animal '_' days{curr_abs_session}];
im_summed_folder = dir([im_dir filesep 'summed*']);
summed_movie_dir = [im_dir filesep im_summed_folder.name];
summed_movie_files = dir([summed_movie_dir filesep '*summed*.tif']);
summed_movie_filename = [summed_movie_dir filesep summed_movie_files.name];

imageinfo=imfinfo(summed_movie_filename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

im = zeros(N,M,numframes,'double');
disp('Loading summed movie');
for loadframe = 1:numframes
    curr_frame = imread(summed_movie_filename,'tiff',loadframe,'Info',imageinfo);
    im(:,:,loadframe) = curr_frame;
    disp(round(loadframe*100/numframes)/100);
end
disp('Loaded')

% Load ROIs
roi_file = dir([im_dir filesep '*.roi']);
if length(roi_file) ~= 1
    error('Bad ROI file');
end
roi_filename = [im_dir filesep roi_file.name];
load(roi_filename,'-MAT');


%% Plot cell of choice from previous event detection

curr_cell = 125;

switch scope
    case 1
        sum_frames = 10;
    case 2
        sum_frames = 50;
end

figure; hold on
plot((1:size(curr_concat_data,2))/sum_frames,curr_concat_data(curr_cell,:),'k');
plot((1:size(curr_concat_data,2))/sum_frames,curr_concat_data_thresh(curr_cell,:),'r');


%% Pull out active/inactive images of cell
% potentially useful for going over lots of data manually but not sure how
% to use it at the moment
% NOTE: This doesn't correctly align at the moment: I think it's because
% considering borders, the frames don't exactly line up to
% frames/sum_frames

cell_border = 100;
figure;
for curr_cell = 1:size(curr_concat_data,1);
    
    if all(curr_concat_data_thresh(curr_cell,:) == 0)
        continue
    end
    
    cell_x = round(cell_border/2);
    cell_y = round(cell_border*(N/M)/2);
    [~,cx,cy] = polycenter(polygon.ROI{curr_cell}(:,1),polygon.ROI{curr_cell}(:,2));
    cx = round(cx);
    cy = round(cy);
    cell_summed = im(cy-cell_y:cy+cell_y,cx-cell_x:cx+cell_x,:);
    
    % Parse active/inactive portions of the trace
    active_frames = curr_concat_data_thresh(curr_cell,:) ~= 0;
    change_epochs = diff(find(diff([2 active_frames 2]) ~= 0));
    
    parsed_active_frames = mat2cell(active_frames,1,change_epochs);
    parsed_trace = mat2cell(curr_concat_data(curr_cell,:),1,change_epochs);
    % Define each epoch as 1 or 0 (should be alternating, should be one type)
    parsed_state = cellfun(@(x) unique(x),parsed_active_frames);
    
    % Get summed frames corresponding to active and inactive portions
    summed_frames = round((1:size(curr_concat_data,2))/sum_frames);
    summed_frames(summed_frames < 1) = 1;
    summed_frames(summed_frames > size(im,3)) = size(im,3);
    
    parsed_summed_frames = mat2cell(summed_frames,1,change_epochs);
    
    parsed_cell_im = cellfun(@(x) nanmean(cell_summed(:,:,unique(x)),3),parsed_summed_frames,'uni',false);
    
    % Pull out prctiles of active and inactive frames for comparison
    parsed_cell_active_im = prctile(cell_summed(:,:,unique(horzcat(parsed_summed_frames{parsed_state == 1}))),95,3);
    parsed_cell_inactive_im = prctile(cell_summed(:,:,unique(horzcat(parsed_summed_frames{parsed_state == 0}))),20,3);
    
    parsed_cell_df = (parsed_cell_active_im - parsed_cell_inactive_im)./parsed_cell_inactive_im;
    
    % Plot pseudo-df with ROI on top
    imagesc(parsed_cell_df);colormap(gray);
    if polygon.ROI{curr_cell}(1,1) ~= polygon.ROI{curr_cell}(end,1)
        polygon.ROI{curr_cell} = [polygon.ROI{curr_cell};polygon.ROI{curr_cell}(1,:)];
    end
    line(polygon.ROI{curr_cell}(:,1)-cx+cell_x,polygon.ROI{curr_cell}(:,2)-cy+cell_y);
    
    waitforbuttonpress

end

% figure;
% active_epochs = find(parsed_state);
% for curr_activity = 1:length(active_epochs)
%     curr_epoch = active_epochs(curr_activity);
%     curr_caxis = minmax(reshape(cat(3,parsed_cell_im{curr_epoch-1:curr_epoch+1}),1,[]));
%     subplot(1,3,1);colormap(gray);
%     imagesc(parsed_cell_im{curr_epoch-1});caxis(curr_caxis);axis off;
%     subplot(1,3,2);colormap(gray);
%     imagesc(parsed_cell_im{curr_epoch});caxis(curr_caxis);axis off;
%     subplot(1,3,3);colormap(gray);
%     imagesc(parsed_cell_im{curr_epoch+1});caxis(curr_caxis);axis off;
%     waitforbuttonpress;
% end
% 
% 
% % TO DO HERE: try getting just the ROI pixels for every "active" event, get
% % correlation matrix: are these events similar to each other in terms of
% % relative pixel values? or are there groupings that don't fit (i.e.
% % neuropil, other cells, etc)?
% im_reshape = reshape(im,N*M,[]);
% curr_mask = reshape(poly2mask(polygon.ROI{curr_cell}(:,1),polygon.ROI{curr_cell}(:,2),N,M),[],1);
% parsed_cell_px = cellfun(@(x) im_reshape(curr_mask,unique(x)),parsed_summed_frames,'uni',false);

%% Show movie of cell during activity

% different approach: get max image of every active epoch, get the max
% of 5 frames +/- active

summed_frames = round((1:size(curr_concat_data,2))/sum_frames);
summed_frames(summed_frames < 1) = 1;
summed_frames(summed_frames > size(im,3)) = size(im,3);

parsed_summed_frames = mat2cell(summed_frames,1,change_epochs);

parsed_cell_im = cellfun(@(x) nanmean(cell_summed(:,:,unique(x)),3),parsed_summed_frames,'uni',false);

compare_frames = 5;

figure;
active_epochs = find(parsed_state);
for curr_activity = 1:length(active_epochs)
    curr_epoch = active_epochs(curr_activity);
    curr_active_frame = parsed_summed_frames{curr_epoch}(1);
    curr_compare_frames = curr_active_frame-compare_frames:curr_active_frame+compare_frames;
    
    curr_compare_frames(curr_compare_frames < 1 | ...
        curr_compare_frames > size(cell_summed,3)) = [];
    curr_caxis = [0 2000];
    
    i = 0;
    while i == 0;
        for curr_frame = curr_compare_frames
            imagesc(cell_summed(:,:,curr_frame));
            colormap(gray);caxis(curr_caxis);axis off;
            pause(0.1);
        end
    end
    
end



%% Get the pixel values for the surrounding max of each event
thresh_vals = [1 1.5 2:0.2:4];
thresh_px_corr = nan(size(curr_concat_data,1),length(thresh_vals));

im_reshape = reshape(im,N*M,[]);
for curr_thresh = 1:length(thresh_vals)
    
    curr_thresh_val = thresh_vals(curr_thresh);
    curr_concat_data_thresh = AP_caEvents_threshtest(curr_concat_data,curr_thresh_val);
    
    active_px_corr = nan(size(curr_concat_data_thresh,1),1);
    for curr_cell = 1:size(curr_concat_data_thresh,1);
        
        % Parse active/inactive portions of the trace
        active_frames = curr_concat_data_thresh(curr_cell,:) ~= 0;
        change_epochs = diff(find(diff([2 active_frames 2]) ~= 0));
        
        parsed_active_frames = mat2cell(active_frames,1,change_epochs);
        parsed_trace = mat2cell(curr_concat_data(curr_cell,:),1,change_epochs);
        % Define each epoch as 1 or 0 (should be alternating, should be one type)
        parsed_state = cellfun(@(x) unique(x),parsed_active_frames);
        
        summed_frames = round((1:size(curr_concat_data,2))/sum_frames);
        summed_frames(summed_frames < 1) = 1;
        summed_frames(summed_frames > size(im,3)) = size(im,3);
        
        parsed_summed_frames = mat2cell(summed_frames,1,change_epochs);
        
        search_active_frames = cellfun(@(x) x(1):x(1)+6,parsed_summed_frames(parsed_state == 1),'uni',false);
        search_active_frames_use = cellfun(@(x) x(x < size(im,3)),search_active_frames,'uni',false);       
        
        curr_mask = reshape(poly2mask(polygon.ROI{curr_cell}(:,1),polygon.ROI{curr_cell}(:,2),N,M),[],1);
        
        [im_active_f im_active_idx] = cellfun(@(x) nanmax(nanmean(im_reshape(curr_mask,x)),[],2), ...
            search_active_frames_use,'uni',false);
        
        im_active = cellfun(@(x,y) im_reshape(curr_mask,x(1)+y), ...
            parsed_summed_frames(parsed_state == 1),im_active_idx,'uni',false);
        
        if length(parsed_cell_im_active) > 4;
            active_px_corr(curr_cell) = nanmean(AP_itril(corrcoef(horzcat(im_active{:}))));
        end
                
    end
    thresh_px_corr(:,curr_thresh) = active_px_corr;
    disp(curr_thresh);
    
end



%% Manually score activity and compare to thresholds

% different approach: get max image of every active epoch, get the max
% of 5 frames +/- active

summed_frames = round((1:size(curr_concat_data,2))/sum_frames);
summed_frames(summed_frames < 1) = 1;
summed_frames(summed_frames > size(im,3)) = size(im,3);

cell_x = round(cell_border/2);
cell_y = round(cell_border*(N/M)/2);
[~,cx,cy] = polycenter(polygon.ROI{curr_cell}(:,1),polygon.ROI{curr_cell}(:,2));
cx = round(cx);
cy = round(cy);
cell_summed = im(cy-cell_y:cy+cell_y,cx-cell_x:cx+cell_x,:);
compare_frames = 5;

figure;
active_epochs = find(parsed_state);
for curr_activity = 1:length(active_epochs)
    curr_epoch = active_epochs(curr_activity);
    curr_active_frame = parsed_summed_frames{curr_epoch}(1);
    curr_compare_frames = curr_active_frame-compare_frames:curr_active_frame+compare_frames;
    
    curr_compare_frames(curr_compare_frames < 1 | ...
        curr_compare_frames > size(cell_summed,3)) = [];
    curr_caxis = [0 2000];
    
    imagesc(cell_summed(:,:,curr_active_frame));
    colormap(gray);caxis(curr_caxis);axis off;
        
    if polygon.ROI{curr_cell}(1,1) ~= polygon.ROI{curr_cell}(end,1)
        polygon.ROI{curr_cell} = [polygon.ROI{curr_cell};polygon.ROI{curr_cell}(1,:)];
    end
    line(polygon.ROI{curr_cell}(:,1)-cx+cell_x,polygon.ROI{curr_cell}(:,2)-cy+cell_y);    
    
    waitforbuttonpress;
end









