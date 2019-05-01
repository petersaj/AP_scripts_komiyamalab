%% Load and prepare data from corticospinal animals

clear all;
animals = {'AP140' 'AP141' 'AP142' 'AP145' 'AP146' 'AP147' 'AP148' 'AP150'};
data = AP_corticospinal_load_all(animals);
analysis = AP_corticospinal_prepare_loaded(data);
classified_rois = AP_classify_movement_cells_continuous(data,analysis);


%% Load and prepare data from corticospinal no-task animals

clear all;
animals = {'AP152','AP153','AP154','AP156','AP157','AP159','AP160'};

data = AP_corticospinal_load_all(animals);
analysis = AP_corticospinal_prepare_loaded(data);
classified_rois = AP_classify_movement_cells_continuous(data,analysis);


% AP152 has 2 days missing, add nans to classified and shift data
AP152 = find(cellfun(@(x) strcmp('AP152',x),{data(:).animal}));

if ~isempty(AP152) && size(classified_rois(AP152).movement,2) < 14
    
    classified_rois(AP152).movement = +classified_rois(AP152).movement;
    classified_rois(AP152).quiescent = +classified_rois(AP152).quiescent;
    classified_rois(AP152).unclassified_active = +classified_rois(AP152).unclassified_active;
    
    classified_rois(AP152).movement(:,4:14) = classified_rois(AP152).movement(:,2:12);
    classified_rois(AP152).quiescent(:,4:14) = classified_rois(AP152).quiescent(:,2:12);
    classified_rois(AP152).unclassified_active(:,4:14) = classified_rois(AP152).unclassified_active(:,2:12);
    
    classified_rois(AP152).movement(:,2:3) = NaN;
    classified_rois(AP152).quiescent(:,2:3) = NaN;
    classified_rois(AP152).unclassified_active(:,2:3) = NaN;
    
    % set up empty days to switch
    data(AP152).im(14).roi_trace_df = [];
    data(AP152).bhv(14).framerate = [];
    analysis(AP152).lever(14).lever_move_frames = [];
    analysis(AP152).im(14).move_onset_aligned = [];
    
    data(AP152).im([4:14,2:3]) = data(AP152).im([2:12,13:14]);
    data(AP152).bhv([4:14,2:3]) = data(AP152).bhv([2:12,13:14]);
    analysis(AP152).lever([4:14,2:3]) = analysis(AP152).lever([2:12,13:14]);
    analysis(AP152).im([4:14,2:3]) = analysis(AP152).im([2:12,13:14]);  
    
end

%% Load and prepare data from L2/3 (old data in new structures)

clear all;
[data,analysis] = AP_corticospinal_load_layer23;
classified_rois = AP_classify_movement_cells_continuous(data,analysis);


%% 1E: example dendrite/soma images and traces

clearvars -except data analysis classified_rois

% Load example cropped image
crop_img_filename = '/usr/local/lab/People/Andy/Corticospinal/corticospinal_paper/figures/dendrite_soma_examples_figure/151212_AP163_region3_cell1_crop.tif';
im = double(AP_load_tiff(crop_img_filename));

% Plot images and projections
figure;
subplot(3,2,1);
imagesc(im(:,:,1));colormap(gray);
axis off
axis square
title('Dendrite plane')

subplot(3,2,2);
conv_filt = fspecial('average',[5,5]);
imagesc(conv2(im(:,:,end-4),conv_filt,'same'));colormap(gray);
axis off
axis square
title('Soma plane')

subplot(3,2,3);
curr_projection_norm = max(bsxfun(@times,im,1./ ...
    repmat(max(max(im,[],2),[],1),size(im,1),size(im,2))),[],3);
conv_filt = fspecial('average',[5,5]);
curr_projection_norm = conv2(curr_projection_norm,conv_filt,'same');
imagesc(curr_projection_norm);colormap(gray);
axis off
axis square
title('Z-projection');

subplot(3,2,4);
curr_projection = permute(max(im,[],1),[3,2,1]);
conv_filt = fspecial('average',[2,5]);
curr_projection_norm = conv2(bsxfun(@times,curr_projection, ...
    1./max(curr_projection,[],2)),conv_filt,'same');
imagesc(curr_projection_norm);colormap(gray);
axis off
title('X-projection (z max-normalized)')

subplot(3,2,5);
curr_projection = permute(max(im,[],2),[3,1,2]);
conv_filt = fspecial('average',[2,5]);
curr_projection_norm = conv2(bsxfun(@times,curr_projection, ...
    1./max(curr_projection,[],2)),conv_filt,'same');
imagesc(curr_projection_norm);colormap(gray);
axis off
title('Y-projection (z max-normalized)')

% Get and plot trace
subplot(3,2,6); hold on;
curr_example_dayanimal = '151212_AP163';
curr_example_region = 'region3';
curr_example_cell = '1';
curr_example_path = [ ...
    '/usr/local/lab/People/Andy/Data/corticospinal_fastz_imaging/' ...
    curr_example_dayanimal filesep curr_example_region filesep curr_example_cell];

roi_path = [curr_example_path filesep 'motion_corrected/summed'];
roi_dir = dir([roi_path filesep '*.roi']);
roi_filenames = cellfun(@(x) [roi_path filesep x],sort({roi_dir(:).name}),'uni',false);

mc_path = [curr_example_path filesep 'motion_corrected'];
roi_trace = AP_getConcatTrace_continuous_batch_fastz( ...
    roi_filenames,[1,2],mc_path);

dendrite_trace_df = cell2mat(arrayfun(@(x) (roi_trace{1}(x,:) - ...
    prctile(roi_trace{1}(x,:),10))./ ...
    prctile(roi_trace{1}(x,:),10),[1:size(roi_trace{1},1)]','uni',false));

soma_trace_df = cell2mat(arrayfun(@(x) (roi_trace{2}(x,:) - ...
    prctile(roi_trace{2}(x,:),10))./ ...
    prctile(roi_trace{2}(x,:),10),[1:size(roi_trace{2},1)]','uni',false));

dendrite_trace_df_norm = bsxfun(@times,dendrite_trace_df,1./max(dendrite_trace_df,[],2));
soma_trace_df_norm = bsxfun(@times,soma_trace_df,1./max(soma_trace_df,[],2));

plot(soma_trace_df_norm,'k');
plot(bsxfun(@plus,dendrite_trace_df_norm',1:size(dendrite_trace_df_norm,1)),'r');

xlim([1 length(soma_trace_df)]);
title('Soma/Dendrite fluorescence')
axis off



%% 2D: similarity cutoff for sibling branch identification

clearvars -except data analysis classified_rois

pairwise_corr_all = cell(2,2);

% ANIMAL 1

% Set paths
data_path = '/usr/local/lab/People/Andy/Data/AP162';
roi_labels_filename = [data_path filesep 'AP162_roilabels.roilabel'];
data_filename = [data_path filesep 'AP162_rois' filesep ...
    '151209_AP162_001_001_corrected_summed_50_analysis'];

% Load data
data = load(data_filename);

% Mark bad ROIs manually
bad_rois = [17,40,12];

% Load branch labels
curr_labels = load(roi_labels_filename,'-MAT');
branched_rois = cellfun(@(x) ~isempty(x),curr_labels.roi_labels);
num_labels = nan(length(curr_labels.roi_labels),1);
num_labels(branched_rois) = ...
    cellfun(@(x) str2num(x{1}),curr_labels.roi_labels(branched_rois));
num_labels(bad_rois) = [];

pairwise_branch = AP_itril(repmat(num_labels,1,length(num_labels)) == ...
    repmat(num_labels',length(num_labels),1),-1);

% Smooth data
curr_smoothpeak = nan(size(data.im.roi_trace_df));
for curr_roi = 1:size(curr_smoothpeak,1)
    curr_smoothpeak(curr_roi,:) = smooth(data.im. ...
        roi_trace_df(curr_roi,:),100,'loess');
end

% Get normalized dot product
curr_smoothpeak_norm = ...
    bsxfun(@times,curr_smoothpeak,1./sqrt(sum(curr_smoothpeak.^2,2)));

curr_smoothpeak_norm(bad_rois,:) = [];
pairwise_corr = AP_itril(curr_smoothpeak_norm*curr_smoothpeak_norm',-1);

pairwise_corr_all(1,:) = {pairwise_corr(pairwise_branch),pairwise_corr(~pairwise_branch)};

% ANIMAL 2

% Set paths
data_path = '/usr/local/lab/People/Andy/Data/AP163';
roi_labels_filename = [data_path filesep 'AP163_roilabels.roilabel'];
data_filename = [data_path filesep 'AP163_rois' filesep ...
    '151209_AP163_001_001_corrected_summed_50_analysis'];

% Load data
data = load(data_filename);

% Mark bad ROIs manually
bad_rois = [12,31,20,43,17,27,37,38,55,47,48,35];

% Load branch labels
curr_labels = load(roi_labels_filename,'-MAT');
branched_rois = cellfun(@(x) ~isempty(x),curr_labels.roi_labels);
num_labels = nan(length(curr_labels.roi_labels),1);
num_labels(branched_rois) = ...
    cellfun(@(x) str2num(x{1}),curr_labels.roi_labels(branched_rois));
num_labels(bad_rois) = [];

pairwise_branch = AP_itril(repmat(num_labels,1,length(num_labels)) == ...
    repmat(num_labels',length(num_labels),1),-1);

% Smooth data
curr_smoothpeak = nan(size(data.im.roi_trace_df));
for curr_roi = 1:size(curr_smoothpeak,1)
    curr_smoothpeak(curr_roi,:) = smooth(data.im. ...
        roi_trace_df(curr_roi,:),100,'loess');
end

% Get normalized dot product
curr_smoothpeak_norm = ...
    bsxfun(@times,curr_smoothpeak,1./sqrt(sum(curr_smoothpeak.^2,2)));

curr_smoothpeak_norm(bad_rois,:) = [];
pairwise_corr = AP_itril(curr_smoothpeak_norm*curr_smoothpeak_norm',-1);

pairwise_corr_all(2,:) = {pairwise_corr(pairwise_branch),pairwise_corr(~pairwise_branch)};

% Combine pairs across animals
grp_sib = vertcat(pairwise_corr_all{:,1});
grp_nonsib = vertcat(pairwise_corr_all{:,2});

% Plot histograms for each population indivually normalized
bin_edges = linspace(-1,1,100);
bin_centers = [bin_edges(1:end-1) + diff(bin_edges)/2];
bin_edges(end) = Inf;

grp_sib_hist = histc(grp_sib,bin_edges)./length(grp_sib);
grp_nonsib_hist = histc(grp_nonsib,bin_edges)./length(grp_nonsib);

figure; hold on;
bar(bin_centers,grp_sib_hist(1:end-1)./max(grp_sib_hist), ...
    'FaceColor','r','EdgeColor','none','BarWidth',1);
bar(bin_centers,grp_nonsib_hist(1:end-1)./max(grp_nonsib_hist), ...
    'FaceColor','k','EdgeColor','none','BarWidth',1);
xlabel('Normalized dot product')
ylabel('Normalized fraction of pairs')
line([0.8,0.8],ylim,'color','k','linestyle','--','linewidth',2);


%% 2E: example population traces

clearvars -except data analysis classified_rois

example_data_file = '/usr/local/lab/People/Andy/Data/AP141/AP141_batch_thresh_roi/150726_AP141_001_001_summed_50_analysis.mat';
curr_data = load(example_data_file);

use_rois = [4,5,6,7,15,20,41,21,24,28,30,35,37,38,44,46,50,52,13,18];

curr_df = curr_data.im.roi_trace_df(use_rois,:);

% Threshold
event_thresh = 3;
method = 2;
curr_thresh = AP_caEvents_thresh(curr_df,event_thresh,method);

% use 5 minutes of data
use_time = 5*60;
first_frame = 13200;
last_frame = find(curr_data.bhv.frame_times < ...
    curr_data.bhv.frame_times(first_frame) + use_time,1,'last');
frame_time = (curr_data.bhv.frame_times(first_frame:last_frame) - ...
    curr_data.bhv.frame_times(first_frame))/60;

% Plot, color example branch ROIs (x-axis is in minutes)
df_space = 7;
figure; hold on;
plot(frame_time,bsxfun(@plus,curr_df(1:end-2,first_frame:last_frame)', ...
    (1:size(curr_df,1)-2)*df_space),'k');
plot(frame_time,bsxfun(@plus,curr_df(end-1:end,first_frame:last_frame)', ...
    (size(curr_df,1)-1:size(curr_df,1))*df_space),'b');

% Draw scalebar
y_scalebar = 10; % df/f
x_scalebar = 1; % minutes

line([0,0],[0,y_scalebar],'color','k','linewidth',3);
line([0,x_scalebar],[0,0],'color','k','linewidth',3);


%% 3B: movement correlation (CR movements)

num_animals = length(data);

% Plot the trial-by-trial movement correlation
movement_start_time = 1001; % (defined in prepare_processed);
movement_use_time = 3000; % ms 

% (do in loop in case dispatcher file not saved, was in old cohort)
max_sessions = max(arrayfun(@(x) length(analysis(x).lever),1:num_animals));
movement_use_trials = cell(num_animals,1);
for curr_animal = 1:num_animals    
   for curr_session = 1:length(analysis(curr_animal).lever)
      if ~isempty(analysis(curr_animal).lever(curr_session).rewarded_movement)
          
          curr_move = ...
              horzcat(analysis(curr_animal).lever(curr_session).rewarded_movement_fixtime{ ...
              analysis(curr_animal).lever(curr_session).cued_movement_trials});
          
          movement_use_trials{curr_animal}{curr_session} = curr_move;          
      end
   end
end

movement_corr_grid = nan(max_sessions,max_sessions,num_animals);
for curr_animal = 1:num_animals
    for session_x = 1:length(analysis(curr_animal).lever);
        for session_y = 1:length(analysis(curr_animal).lever)            
            if ~isempty(movement_use_trials{curr_animal}{session_x}) && ...
                    ~isempty(movement_use_trials{curr_animal}{session_y})
                
                num_trials_x = size(movement_use_trials{curr_animal}{session_x},2);
                
                curr_move_corr = corrcoef([movement_use_trials{ ...
                    curr_animal}{session_x}(movement_start_time:movement_start_time+movement_use_time,:) ...
                    movement_use_trials{curr_animal}{session_y}( ...
                    movement_start_time:movement_start_time+movement_use_time,:)]);
                
                curr_use_corr = curr_move_corr(1:num_trials_x,num_trials_x+1:end);
                
                movement_corr_grid(session_x,session_y,curr_animal) = nanmedian(curr_use_corr(:));
                
            end
        end
    end
end

% Movement correlation
figure; 
subplot(1,3,1);
imagesc(nanmean(movement_corr_grid,3));colormap(gray);
ylabel('Day');
xlabel('Day');
title('Movement correlation');
axis square;
colorbar;

first_movecorr_diag = cell2mat(arrayfun(@(x) diag(movement_corr_grid(:,:,x)), ...
    1:size(movement_corr_grid,3),'uni',false))';
second_movecorr_diag = cell2mat(arrayfun(@(x) diag(movement_corr_grid(:,:,x),-1), ...
    1:size(movement_corr_grid,3),'uni',false))';

subplot(1,3,2);
errorbar(nanmean(first_movecorr_diag),nanstd(first_movecorr_diag)./ ...
    sqrt(sum(~isnan(first_movecorr_diag))),'k','linewidth',2);
xlim([0 15]);
ylabel('Median movement correlation')
xlabel('Day');
title('Within day');

subplot(1,3,3);
errorbar(nanmean(second_movecorr_diag),nanstd(second_movecorr_diag)./ ...
    sqrt(sum(~isnan(second_movecorr_diag))),'k','linewidth',2);
xlim([0 14]);
ylabel('Median movement correlation')
xlabel('Day v. next day');
title('Adjacent day');

% Statistics
day_grid = repmat(1:14,length(data),1);
[rw,pw] = corrcoef(day_grid(:),first_movecorr_diag(:));

day_grid = repmat(1:13,length(data),1);
[ra,pa] = corrcoef(day_grid(:),second_movecorr_diag(:));


%% 4A: heatmap of activity of all movements (on/off stitch) all cells sorted

% Get activity with all movements > n sec

move_act_onset_all = cell(length(data),14);
move_act_offset_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Get onsets of movement
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 28*2;
        max_move_time = Inf;
        min_pre_move_iti = 28*1;
        min_post_move_iti = 28*1;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Pull out timed activity
        nonmove_frames = 28*1;
        move_frames = 28*1;        
        total_frames = nonmove_frames+move_frames+1;
        
        move_frame_timed_onset_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-nonmove_frames: ...
            use_movement_start_frames(x)+move_frames,1:length(use_movement_start_frames),'uni',false)'); 
        
        move_frame_timed_offset_idx = cell2mat(arrayfun(@(x) use_movement_stop_frames(x)-move_frames: ...
            use_movement_stop_frames(x)+nonmove_frames,1:length(use_movement_stop_frames),'uni',false)');
        
        discard_movements = any(move_frame_timed_onset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_onset_idx < 1,2) | ...
            any(move_frame_timed_offset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_offset_idx < 1,2);
        
        move_frame_timed_onset_idx(discard_movements,:) = [];
        move_frame_timed_onset_idx = reshape(move_frame_timed_onset_idx',[],1);
        
        move_frame_timed_offset_idx(discard_movements,:) = [];
        move_frame_timed_offset_idx = reshape(move_frame_timed_offset_idx',[],1);
        
        move_act_onset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_onset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        move_act_offset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_offset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed;

        
    end
    disp(curr_animal);
end

move_act_onset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_onset_all,'uni',false);
move_act_offset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_offset_all,'uni',false);

use_animal = 3;
use_day = 9;

curr_data = [move_act_onset_all_mean{use_animal,use_day}, ...
    move_act_offset_all_mean{use_animal,use_day}];

curr_data_minsub = bsxfun(@minus,curr_data,min(curr_data,[],2));
curr_data_norm = bsxfun(@times,curr_data_minsub,1./max(curr_data_minsub,[],2));

m = classified_rois(use_animal).movement(:,use_day) == 1;
q = classified_rois(use_animal).quiescent(:,use_day) == 1;
u = ~m & ~q;

curr_data_m = curr_data_norm(m,:);
curr_data_q = curr_data_norm(q,:);
curr_data_u = curr_data_norm(u,:);

% Sort by difference between m/q activity
m_frames = false(size(curr_data_norm,2),1);
m_frames(nonmove_frames+1:nonmove_frames+1+move_frames*2) = true;

curr_data_m_diff = ...
    (nanmean(curr_data_m(:,m_frames),2)-nanmean(curr_data_m(:,~m_frames),2))./ ...
    (nanmean(curr_data_m(:,m_frames),2)+nanmean(curr_data_m(:,~m_frames),2));

curr_data_q_diff = ...
    (nanmean(curr_data_q(:,m_frames),2)-nanmean(curr_data_q(:,~m_frames),2))./ ...
    (nanmean(curr_data_q(:,m_frames),2)+nanmean(curr_data_q(:,~m_frames),2));

curr_data_u_diff = ...
    (nanmean(curr_data_u(:,m_frames),2)-nanmean(curr_data_u(:,~m_frames),2))./ ...
    (nanmean(curr_data_u(:,m_frames),2)+nanmean(curr_data_u(:,~m_frames),2));

[~,sort_idx_m] = sort(curr_data_m_diff);
[~,sort_idx_q] = sort(curr_data_q_diff);
[~,sort_idx_u] = sort(curr_data_u_diff);

curr_data_sort = [curr_data_m(sort_idx_m,:); ...
    curr_data_u(sort_idx_u,:); ...
    curr_data_q(sort_idx_q,:)];

num_rois = cumsum([sum(m),sum(u)]);

figure;
subplot(1,2,1);
imagesc(curr_data_sort(:,1:total_frames));
colormap(gray);
line(xlim,[sum(m),sum(m)],'color','r');
line(xlim,[sum(m)+sum(u),sum(m)+sum(u)],'color','r');
line(repmat(nonmove_frames,2,1),ylim,'color','r');
ylabel('Sorted ROI')
xlabel('Frame');

subplot(1,2,2);
imagesc(curr_data_sort(:,end-(total_frames-1):end));
colormap(gray);
line(xlim,[sum(m),sum(m)],'color','r');
line(xlim,[sum(m)+sum(u),sum(m)+sum(u)],'color','r');
line(repmat(move_frames,2,1),ylim,'color','r');
ylabel('Sorted ROI')
xlabel('Frame');



%% 4B/D: average activity of all neurons aligned to movement onset/offset (all movements)

clearvars -except data analysis classified_rois

move_act_onset_all = cell(length(data),14);
move_act_offset_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Get onsets of movement
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 28*2;
        max_move_time = Inf;
        min_pre_move_iti = 28*1;
        min_post_move_iti = 28*1;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Pull out timed activity
        nonmove_frames = 28*3;
        move_frames = 28*1;        
        total_frames = nonmove_frames+move_frames+1;
        
        move_frame_timed_onset_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-nonmove_frames: ...
            use_movement_start_frames(x)+move_frames,1:length(use_movement_start_frames),'uni',false)'); 
        
        move_frame_timed_offset_idx = cell2mat(arrayfun(@(x) use_movement_stop_frames(x)-move_frames: ...
            use_movement_stop_frames(x)+nonmove_frames,1:length(use_movement_stop_frames),'uni',false)');
        
        discard_movements = any(move_frame_timed_onset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_onset_idx < 1,2) | ...
            any(move_frame_timed_offset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_offset_idx < 1,2);
        
        move_frame_timed_onset_idx(discard_movements,:) = [];
        move_frame_timed_onset_idx = reshape(move_frame_timed_onset_idx',[],1);
        
        move_frame_timed_offset_idx(discard_movements,:) = [];
        move_frame_timed_offset_idx = reshape(move_frame_timed_offset_idx',[],1);
        
        move_act_onset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_onset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        move_act_offset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_offset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed;

        
    end
    disp(curr_animal);
end

move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_all,'uni',false);
move_act_m_onset_mean = cell(size(move_act_onset_mean));
move_act_q_onset_mean = cell(size(move_act_onset_mean));
move_act_ua_onset_mean = cell(size(move_act_onset_mean));

move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_all,'uni',false);
move_act_m_offset_mean = cell(size(move_act_offset_mean));
move_act_q_offset_mean = cell(size(move_act_offset_mean));
move_act_ua_offset_mean = cell(size(move_act_offset_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        move_act_ua_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).unclassified_active(:,curr_day) == 1),1),3);
        
        move_act_m_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);   
        
        move_act_ua_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).unclassified_active(:,curr_day) == 1),1),3);  
        
    end
end

% Get mean activity across days for each animal
move_act_mean_onset_animalavg = nan(length(data),total_frames);
move_act_mean_m_onset_animalavg = nan(length(data),total_frames);
move_act_mean_q_onset_animalavg = nan(length(data),total_frames);
move_act_mean_ua_onset_animalavg = nan(length(data),total_frames);

move_act_mean_offset_animalavg = nan(length(data),total_frames);
move_act_mean_m_offset_animalavg = nan(length(data),total_frames);
move_act_mean_q_offset_animalavg = nan(length(data),total_frames);
move_act_mean_ua_offset_animalavg = nan(length(data),total_frames);

for curr_animal = 1:length(data)
    
    move_act_mean_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_onset_mean{curr_animal,:}),1);
    move_act_mean_m_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_onset_mean{curr_animal,:}),1);
    move_act_mean_q_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_onset_mean{curr_animal,:}),1);
    move_act_mean_ua_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_ua_onset_mean{curr_animal,:}),1);
    
    move_act_mean_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_offset_mean{curr_animal,:}),1);
    move_act_mean_m_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_offset_mean{curr_animal,:}),1);
    move_act_mean_q_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_offset_mean{curr_animal,:}),1);
    move_act_mean_ua_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_ua_offset_mean{curr_animal,:}),1);
    
end

all_avg_onset = nanmean(move_act_mean_onset_animalavg,1);
m_avg_onset = nanmean(move_act_mean_m_onset_animalavg,1);
q_avg_onset = nanmean(move_act_mean_q_onset_animalavg,1);
ua_avg_onset = nanmean(move_act_mean_ua_onset_animalavg,1);

all_sem_onset = nanstd(move_act_mean_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_onset_animalavg)));
m_sem_onset = nanstd(move_act_mean_m_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_m_onset_animalavg)));
q_sem_onset = nanstd(move_act_mean_q_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_q_onset_animalavg)));
ua_sem_onset = nanstd(move_act_mean_ua_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_ua_onset_animalavg)));

all_avg_offset = nanmean(move_act_mean_offset_animalavg,1);
m_avg_offset = nanmean(move_act_mean_m_offset_animalavg,1);
q_avg_offset = nanmean(move_act_mean_q_offset_animalavg,1);
ua_avg_offset = nanmean(move_act_mean_ua_offset_animalavg,1);

all_sem_offset = nanstd(move_act_mean_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_offset_animalavg)));
m_sem_offset = nanstd(move_act_mean_m_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_m_offset_animalavg)));
q_sem_offset = nanstd(move_act_mean_q_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_q_offset_animalavg)));
ua_sem_offset = nanstd(move_act_mean_ua_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_ua_offset_animalavg)));

% Plot all ROIs
df_scalebar = 0.02;
time_scalebar = 14.6; % frames

ymax = max([all_avg_onset+all_sem_onset,all_avg_offset+all_sem_offset]);

figure;
subplot(1,2,1); hold on;
fill([1:length(all_avg_onset),length(all_avg_onset):-1:1], ...
    [all_avg_onset+all_sem_onset,fliplr(all_avg_onset-all_sem_onset)],[0.5,0.5,0.5]);
plot(all_avg_onset,'k','linewidth',2);
ylim([0,ymax]);
line(repmat(nonmove_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Average \DeltaF/F');
title('Movement onset (all ROIs, all animals, all days)');

subplot(1,2,2); hold on;
fill([1:length(all_avg_offset),length(all_avg_offset):-1:1], ...
    [all_avg_offset+all_sem_offset,fliplr(all_avg_offset-all_sem_offset)],[0.5,0.5,0.5]);
plot(all_avg_offset,'k','linewidth',2);
ylim([0,ymax]);
line(repmat(move_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Average \DeltaF/F');
title('Movement onset (all ROIs, all animals, all days)');

line([0,0],[0,df_scalebar],'color','k','linewidth',2);
line([0,time_scalebar],[0,0],'color','k','linewidth',2);

% Plot classified ROIs
df_scalebar = 0.05;
time_scalebar = 14.6; % frames

ymax = max([m_avg_onset+m_sem_onset,m_avg_offset+m_sem_offset, ...
    q_avg_onset+q_sem_onset,q_avg_offset+q_sem_offset]);

figure;
subplot(1,2,1); hold on;
fill([1:length(m_avg_onset),length(m_avg_onset):-1:1], ...
    [m_avg_onset+m_sem_onset,fliplr(m_avg_onset-m_sem_onset)],[0.5,0.5,0.5]);
plot(m_avg_onset,'color',[0,0.8,0],'linewidth',2);
fill([1:length(q_avg_onset),length(q_avg_onset):-1:1], ...
    [q_avg_onset+q_sem_onset,fliplr(q_avg_onset-q_sem_onset)],[0.5,0.5,0.5]);
plot(q_avg_onset,'color',[0.8,0,0],'linewidth',2);
fill([1:length(ua_avg_onset),length(ua_avg_onset):-1:1], ...
    [ua_avg_onset+ua_sem_onset,fliplr(ua_avg_onset-ua_sem_onset)],[0.5,0.5,0.5]);
plot(ua_avg_onset,'color',[0,0,0.8],'linewidth',2);
ylim([0,ymax]);
line(repmat(nonmove_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Average \DeltaF/F');
title('Movement onset (all ROIs, all animals, all days)');

subplot(1,2,2); hold on;
fill([1:length(m_avg_offset),length(m_avg_offset):-1:1], ...
    [m_avg_offset+m_sem_offset,fliplr(m_avg_offset-m_sem_offset)],[0.5,0.5,0.5]);
plot(m_avg_offset,'color',[0,0.8,0],'linewidth',2);
fill([1:length(q_avg_offset),length(q_avg_offset):-1:1], ...
    [q_avg_offset+q_sem_offset,fliplr(q_avg_offset-q_sem_offset)],[0.5,0.5,0.5]);
plot(q_avg_offset,'color',[0.8,0,0],'linewidth',2);
fill([1:length(ua_avg_offset),length(ua_avg_offset):-1:1], ...
    [ua_avg_offset+ua_sem_offset,fliplr(ua_avg_offset-ua_sem_offset)],[0.5,0.5,0.5]);
plot(ua_avg_offset,'color',[0,0,0.8],'linewidth',2);
ylim([0,ymax]);
line(repmat(move_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Average \DeltaF/F');
title('Movement onset (all ROIs, all animals, all days)');

line([0,0],[0,df_scalebar],'color','k','linewidth',2);
line([0,time_scalebar],[0,0],'color','k','linewidth',2);

% Get time to divergence
baseline_frames = 10:20;
framerate = 29;

m_cat = vertcat(move_act_m_onset_mean{:});
q_cat = vertcat(move_act_q_onset_mean{:});

m_cat_minsub = bsxfun(@minus,m_cat,min(m_cat,[],2));
q_cat_minsub = bsxfun(@minus,q_cat,min(q_cat,[],2));

m_cat_norm = bsxfun(@rdivide,m_cat_minsub,max(m_cat_minsub,[],2));
q_cat_norm = bsxfun(@rdivide,q_cat_minsub,max(q_cat_minsub,[],2));

m_cat_norm_mean = nanmean(m_cat_norm,1);
q_cat_norm_mean = nanmean(q_cat_norm,1);

m_ci = bootci(10000,@mean, ...
    reshape(m_cat(:,baseline_frames),[],1));
m_diverge = (nonmove_frames - ...
    (find(m_avg_onset(baseline_frames(end)+1:end) > m_ci(2),1) + ...
    baseline_frames(end)))/framerate;

q_ci = bootci(10000,@mean, ...
    reshape(q_cat(:,baseline_frames),[],1));
q_diverge = (nonmove_frames - ...
    (find(q_avg_onset(baseline_frames(end)+1:end) < q_ci(1),1) + ...
    baseline_frames(end)))/framerate;

% Get time to divergence individually for each animal
baseline_frames = 10:20;
framerate = 29;

m_diverge = nan(length(data),1);
q_diverge = nan(length(data),1);
for curr_animal = 1:length(data)
    
    m_cat = vertcat(move_act_m_onset_mean{curr_animal,:});
    q_cat = vertcat(move_act_q_onset_mean{curr_animal,:});
    
    m_cat_minsub = bsxfun(@minus,m_cat,min(m_cat,[],2));
    q_cat_minsub = bsxfun(@minus,q_cat,min(q_cat,[],2));
    
    m_cat_norm = bsxfun(@rdivide,m_cat_minsub,max(m_cat_minsub,[],2));
    q_cat_norm = bsxfun(@rdivide,q_cat_minsub,max(q_cat_minsub,[],2));
    
    m_cat_norm_mean = nanmean(m_cat_norm,1);
    q_cat_norm_mean = nanmean(q_cat_norm,1);
    
    m_ci = bootci(10000,@mean, ...
        reshape(m_cat(:,baseline_frames),[],1));
    m_diverge(curr_animal) = (nonmove_frames - ...
        (find(m_avg_onset(baseline_frames(end)+1:end) > m_ci(2),1) + ...
        baseline_frames(end)))/framerate;
    
    q_ci = bootci(10000,@mean, ...
        reshape(q_cat(:,baseline_frames),[],1));
    q_diverge(curr_animal) = (nonmove_frames - ...
        (find(q_avg_onset(baseline_frames(end)+1:end) < q_ci(1),1) + ...
        baseline_frames(end)))/framerate;
    
end

%% 4C: histogram of activity rate relative to movement/quiescence

clearvars -except data analysis classified_rois

% Get fraction of activity during session and during movement
activity_mq_rate_ratio = cell(length(data),1);
activity_mq_rate_ratio_class = cell(length(data),1);

for curr_animal = 1:length(data);
    
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    activity_mq_rate_ratio{curr_animal} = cell(n_days,1);
    activity_mq_rate_ratio_class{curr_animal} = cell(n_days,1);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        % Only use "active" ROIs, at least 5 events
        min_events = 5;
        active_rois = sum(diff(curr_act > 0,[],2) == 1,2) > min_events;
        curr_act = curr_act(active_rois,:);
        
        % "Active" classified ROIs
        curr_class = classified_rois(curr_animal).movement(active_rois,curr_day) | ...
            classified_rois(curr_animal).quiescent(active_rois,curr_day);
        
        curr_act_m = (curr_act*lever_active_frames)./sum(lever_active_frames);
        curr_act_q = (curr_act*~lever_active_frames)./sum(~lever_active_frames);
        curr_act_mq = (curr_act_m - curr_act_q)./(curr_act_m + curr_act_q);
        activity_mq_rate_ratio{curr_animal}{curr_day} = curr_act_mq;
        
        activity_mq_rate_ratio_class{curr_animal}{curr_day} = curr_act_mq(curr_class);      
     
    end
    disp(curr_animal);
end

% Histogram of real vs. shuffled (all days/ROIs concatenated)
activity_mq_ratio_cat = cell2mat((cellfun(@(x) vertcat(x{:}),activity_mq_rate_ratio,'uni',false)));
activity_mq_ratio_class_cat = cell2mat((cellfun(@(x) vertcat(x{:}),activity_mq_rate_ratio_class,'uni',false)));

bin_edges = linspace(-1,1,20);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(activity_mq_ratio_cat,bin_centers);
class_dist = hist(activity_mq_ratio_class_cat,bin_centers);

figure; hold on;
bar(bin_centers,real_dist,'w','linewidth',2);
bar(bin_centers,class_dist,'k','linewidth',2)
legend({'Real','Classified'})

ylabel('Number of ROIs');
xlabel('Rate M-Q/M+Q');

%% 4E: Temporal changes across days (onsets)

clearvars -except data analysis classified_rois

move_act_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Get onsets of movement
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 28*2;
        max_move_time = Inf;
        min_pre_move_iti = 28*1;
        min_post_move_iti = 0;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Pull out timed activity
        back_frames = 28*1;
        forward_frames = 28*2;        
        total_frames = back_frames+forward_frames+1;
        
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        move_act_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        % To use only movement onsets
%         move_act_onsets = [zeros(size(move_act_timed,1),1, ...
%             size(move_act_timed,3)),diff(move_act_timed > 0,[],2) == 1];
%         move_act_timed(~move_act_onsets) = 0;
        
        move_act_all{curr_animal,curr_day} = move_act_timed;    
        
    end
    
    disp(curr_animal);

end

% Get distribution of activity across movement ROIs
baseline_frames = [-20:-10] + back_frames + 1;

move_act_total = nan(forward_frames+1,14,length(data));
move_act_distribution = nan(forward_frames+1,14,length(data));
for curr_animal = 1:length(data)
    
    curr_m = classified_rois(curr_animal).movement == 1;
    
    for curr_day = 1:min(length(data(curr_animal).im),14);
        
        curr_avg_baseline = permute(nanmean(nanmean(move_act_all{curr_animal,curr_day}(:, ...
            baseline_frames,curr_m(:,curr_day)),2),1),[3,1,2]);
        
        curr_act = move_act_all{curr_animal,curr_day}(:,back_frames+1:end,curr_m(:,curr_day));
        curr_act_onset = zeros(size(curr_act));
        use_onsets = [false(size(curr_act,1),1,size(curr_act,3)),diff(curr_act,[],2) > 0];
        curr_act_onset(use_onsets) = curr_act(use_onsets) > 0;
        
        total_act_frame = nanmean(nanmean(curr_act_onset,1),3)';
        total_act_frame_baselinesub = total_act_frame;% - nanmean(curr_avg_baseline);
                
        total_act_frame = smooth(total_act_frame,5);
        total_act_frame_baselinesub = smooth(total_act_frame_baselinesub,5);
        
        total_act_frame_dist = total_act_frame./sum(total_act_frame);
                
        move_act_total(:,curr_day,curr_animal) = total_act_frame_baselinesub;
        move_act_distribution(:,curr_day,curr_animal) = ...
            total_act_frame_dist;
        
    end
    
end

% Plot average across animals, group days
df_scale = 0.01;
time_scale = 29*0.5; % Frames

figure; hold on
day_group = {[1,2,3,4],[5,6,7,8],[9,10,11],[12,13,14]};
%day_group = {[1,2,3],[4,5,6],[7,8,9],[10 11 12],[13 14]};
%day_group = {[1,2],[3,4],[5,6],[7,8],[9,10],[11,12],[13,14]};
% Total activity
move_act_total_mean = nanmean(move_act_distribution,3);
move_act_total_grp = nan(size(move_act_total_mean,1),length(day_group));
for i = 1:length(day_group)
   move_act_total_grp(:,i) = nanmean(move_act_total_mean(:,day_group{i}),2); 
end
set(gca,'ColorOrder',copper(length(day_group)));
plot(move_act_total_grp,'linewidth',2);
xlabel('Time');
ylabel('Fraction of activity onsets');
title('Movement ROIs')
legend(cellfun(@(x) num2str(x),day_group,'uni',false));
axis square

line([0,0],[0,df_scale],'color','k','linewidth',2);
line([0,time_scale],[0,0],'color','k','linewidth',2);

% Statistics
day_group = num2cell(1:14);
move_act_total_grp = nan(size(move_act_total,1),length(day_group),length(data));
for i = 1:length(day_group)
   move_act_total_grp(:,i,:) = nanmean(move_act_total(:,day_group{i},:),2); 
end

move_act_total_early = permute(nanmean(move_act_total_grp(1:10,:,:,1)),[2,3,1]);

d = repmat(transpose(1:length(day_group)),1,length(data));
[re,pe] = corrcoef(d(~isnan(move_act_total_early)),move_act_total_early(~isnan(move_act_total_early)));


%% 4F: Temporal changes within ROIs

clearvars -except data analysis classified_rois

pref_act_time = cell(length(data),1);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Get onsets of movement
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 28*2;
        max_move_time = Inf;
        min_pre_move_iti = 28*1;
        min_post_move_iti = 0;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Pull out timed activity
        back_frames = 28*1;
        forward_frames = 28*2;        
        total_frames = back_frames+forward_frames+1;
        
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        curr_m = classified_rois(curr_animal).movement(:,curr_day) == 1;
        
        move_act_timed_m = permute(reshape(curr_act(curr_m, ...
            move_frame_timed_idx)',total_frames,[],sum(curr_m)),[2,1,3]);
        
        baseline_frames = [-20:-10] + back_frames + 1;
        avg_act_m_baseline = permute(nanmean(nanmean( ...
            move_act_timed_m(:,baseline_frames,:),2),1),[3,1,2]);
        
        avg_act_m = permute(nanmean(move_act_timed_m(:,back_frames+1:end,:),1),[3,2,1]);
        [~,avg_act_m_maxtime] = max(avg_act_m,[],2);
        
        curr_pref_act_time = nan(size(curr_m));
        curr_pref_act_time(curr_m) = avg_act_m_maxtime;
        pref_act_time{curr_animal,curr_day} = curr_pref_act_time;
        
    end
    
    disp(curr_animal);

end

% Plot changes in preferred timing
pref_diff = nan(14,14,length(data));
for curr_animal = 1:length(data)
    
    curr_pref = permute(horzcat(pref_act_time{curr_animal,:}),[2,3,1]);
        
    curr_pref_diff = repmat(curr_pref,[1,size(curr_pref,1)]) - ...
        permute(repmat(curr_pref,[1,size(curr_pref,1)]),[2,1,3]);
        
    pref_diff(1:size(curr_pref,1),1:size(curr_pref,1),curr_animal) = ...
        nanmean(curr_pref_diff,3);
    
end

figure;imagesc(nanmean(pref_diff,3));colormap(gray);
xlabel('Day');
ylabel('Day');
title('Timing changes within movement ROIs')
axis square
c = colorbar;
ylabel(c,'Timing y-x (frames)');

% Statistics
pref_timing_forward = nanmean(cell2mat(arrayfun(@(x) AP_itril(pref_diff(:,:,x),-1), ...
    1:size(pref_diff,3),'uni',false)),2);
[p,h,stats] = signrank(pref_timing_forward);



%% 5A: fraction of movement/quiescent ROIs

clearvars -except data analysis classified_rois

m_frac = cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false);
m_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),m_frac,'uni',false));

q_frac = cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false);
q_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),q_frac,'uni',false));

figure; hold on
errorbar(nanmean(m_frac_pad),nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad))),'k','linewidth',2)
errorbar(nanmean(q_frac_pad),nanstd(q_frac_pad)./sqrt(sum(~isnan(q_frac_pad))),'r','linewidth',2)
ylabel('Fraction of ROIs');
xlabel('Day');
legend({'Movement-related' 'Quiescence-related'})

% Statistics

days_1 = 1:6;
days_2 = 6:14;

% Quiescent
p = anova1(q_frac_pad,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(q_frac_pad(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(q_frac_pad(:,days_2),[],1));

% Movement
p = anova1(m_frac_pad,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(m_frac_pad(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(m_frac_pad(:,days_2),[],1));


%% 5B: example ROIs

roi_plot = { ...
    [120], ... % 1
    [65], ... % 2
    [13,222], ... % 3
    [93,26], ... % 4
    [], ... % 5
    [], ... % 6
    [], ... % 7
    [9]};    % 8

% Get activity with all movements > n sec

move_act_onset_all = cell(length(data),14);
move_act_offset_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Get onsets of movement
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 28*2;
        max_move_time = Inf;
        min_pre_move_iti = 28*1;
        min_post_move_iti = 28*1;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Pull out timed activity
        nonmove_frames = 28*3;
        move_frames = 28*1;        
        total_frames = nonmove_frames+move_frames+1;
        
        move_frame_timed_onset_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-nonmove_frames: ...
            use_movement_start_frames(x)+move_frames,1:length(use_movement_start_frames),'uni',false)'); 
        
        move_frame_timed_offset_idx = cell2mat(arrayfun(@(x) use_movement_stop_frames(x)-move_frames: ...
            use_movement_stop_frames(x)+nonmove_frames,1:length(use_movement_stop_frames),'uni',false)');
        
        discard_movements = any(move_frame_timed_onset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_onset_idx < 1,2) | ...
            any(move_frame_timed_offset_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_offset_idx < 1,2);
        
        move_frame_timed_onset_idx(discard_movements,:) = [];
        move_frame_timed_onset_idx = reshape(move_frame_timed_onset_idx',[],1);
        
        move_frame_timed_offset_idx(discard_movements,:) = [];
        move_frame_timed_offset_idx = reshape(move_frame_timed_offset_idx',[],1);
        
        move_act_onset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_onset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        move_act_offset_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_offset_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed;

        
    end
    disp(curr_animal);
end


example_roi_savedir = ['/usr/local/lab/People/Andy/Corticospinal/corticospinal_paper/figures/example_rois/chosen/revision'];

line_fig = figure('Position',[12,537,1456,361]);
set(line_fig,'PaperPositionMode','auto');

for curr_animal = 1:length(data)
    curr_plot_rois = roi_plot{curr_animal};
    for curr_roi = curr_plot_rois
   
        % Plot as line
        figure(line_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_onset_all(curr_animal,:),'uni',false);
        curr_offset = cellfun(@(x) x(:,:,curr_roi),move_act_offset_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            
            % Plot aligned to movement onset
            subplot(2,7,curr_day);
            cla;
            
            curr_onset_mean = nanmean(curr_onset{curr_day},1);
            curr_onset_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            curr_onset_x = 1:length(curr_onset_mean);
            
            curr_offset_mean = nanmean(curr_offset{curr_day},1);
            curr_offset_sem = nanstd(curr_offset{curr_day},[],1)./sqrt(sum(~isnan(curr_offset{curr_day})));
            curr_offset_x = (1:length(curr_offset_mean)) + curr_onset_x(end) + 1;
            
            jbfill(curr_onset_x,curr_onset_mean+curr_onset_sem, ...
                curr_onset_mean-curr_onset_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_onset_x,curr_onset_mean,'k','linewidth',2);
            
            jbfill(curr_offset_x,curr_offset_mean+curr_offset_sem, ...
                curr_offset_mean-curr_offset_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_offset_x,curr_offset_mean,'k','linewidth',2);
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            elseif classified_rois(curr_animal).unclassified_active(curr_roi,curr_day)
                set(gca,'XColor','b');set(gca,'YColor','b');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.2 1.5]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line(repmat(nonmove_frames+1,1,2),ylim,'color','k');
            line(repmat(curr_onset_x(end)+1,1,2),ylim,'color','k');
            line(repmat(curr_onset_x(end)+1+move_frames,1,2),ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);         
            
        end
        
        curr_savefile = [example_roi_savedir filesep 'A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_line'];
        
        saveas(line_fig,curr_savefile,'fig');
        
    end
end

close(line_fig);

% Get images of chosen example ROIs

chosen_roi_numbers = cellfun(@(chosen,idx,grp) cellfun(@(grp) ...
    idx(grp),grp(chosen),'uni',false), ...
    roi_plot,{data(:).roi_idx},{data(:).roi_group},'uni',false);

% Copy over summed movies and roi files for dendrite/soma
data_path = '/usr/local/lab/People/Andy/Data';

for curr_animal = 1:length(roi_plot)
    
    if isempty(roi_plot{curr_animal})
        continue
    end
    
    animal = data(curr_animal).animal;
    
    curr_data_path = [data_path filesep animal];
    curr_data_dir = dir(curr_data_path);
    day_dir = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d\d\d')),{curr_data_dir.name});
    days = {curr_data_dir(day_dir).name};
    day_filenames = sort(cellfun(@(x) [curr_data_path filesep x],days,'uni',false));
    
    roi_path = [curr_data_path filesep animal '_batch_thresh_roi'];
    curr_roi_dir = dir([roi_path filesep '*.roi']);
    roi_filenames = sort(cellfun(@(x) [roi_path filesep x],{curr_roi_dir.name},'uni',false));
    
    plot_rois = vertcat(chosen_roi_numbers{curr_animal}{:});
    plot_roi_grp = cell2mat(arrayfun(@(x) repmat(roi_plot{curr_animal}(x), ...
        length(chosen_roi_numbers{curr_animal}{x}),1), ...
        transpose(1:length(roi_plot{curr_animal})),'uni',false));
    
    roi_im_px = 20;
    roi_im = nan(roi_im_px*2+1,roi_im_px*2+1,length(days),length(plot_rois));
    roi_polygon = cell(length(plot_rois),length(days));
    
    for curr_day = 1:length(day_filenames)
        
        curr_summed_path = [day_filenames{curr_day} filesep 'summed'];
        curr_summed_filename = dir([curr_summed_path filesep '*_summed_50.tif']);
        im = AP_load_tiff([curr_summed_path filesep curr_summed_filename.name]);
        im_mean = nanmean(im,3);
        
        curr_rois = load([roi_filenames{curr_day}],'-MAT');
        
        for curr_roi_idx = 1:length(plot_rois);
            
            curr_roi = plot_rois(curr_roi_idx);
            
            [~,cx,cy] = polycenter(curr_rois.polygon.ROI{curr_roi}(:,1), ...
                curr_rois.polygon.ROI{curr_roi}(:,2));
          
            roi_im(:,:,curr_day,curr_roi_idx) = im_mean(cy-roi_im_px:cy+roi_im_px, ...
                cx-roi_im_px:cx+roi_im_px);     
            
            
            curr_poly = curr_rois.polygon.ROI{curr_roi};
            roi_polygon{curr_roi_idx,curr_day} = ...
                curr_poly - [repmat(cx,size(curr_poly,1),1), ...
                repmat(cy,size(curr_poly,1),1)]+roi_im_px+1;
            
        end
        
    end
    
    % Plot the image with ROI every day
    save_path = ['/usr/local/lab/People/Andy/Corticospinal/corticospinal_paper/figures/example_rois/chosen/im'];
    roi_fig = figure('Position',[3,489,1865,486]);
    
    for curr_roi_idx = 1:length(plot_rois);
        for curr_day = 1:size(roi_im,3)
            
            subplot(2,7,curr_day); cla;
            imagesc(roi_im(:,:,curr_day,curr_roi_idx));
            colormap(gray)
            axis off;
            axis square;
            set(roi_fig,'Name',[animal '_R' num2str(plot_rois(curr_roi_idx)) ...
                '_Rgrp' num2str(plot_roi_grp(curr_roi_idx))]);
            
            line(roi_polygon{curr_roi_idx,curr_day}(:,1), ...
                roi_polygon{curr_roi_idx,curr_day}(:,2));
            
        end
        
        curr_save_filename = [save_path filesep animal ...
            '_r' num2str(plot_rois(curr_roi_idx)) '_rgrp' num2str(plot_roi_grp(curr_roi_idx))];
        saveas(roi_fig,curr_save_filename,'fig');
    end
    close(roi_fig);
    
end


%% 5C: same classification overlap

clearvars -except data analysis classified_rois

% Same classification overlap
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m2m = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + repmat(sum(m,1),size(m,2),1)),m,'uni',false);
m2m_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),m2m,'uni',false);
m2m_cat = cat(3,m2m_pad{:});

q2q = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + repmat(sum(q,1),size(q,2),1)),q,'uni',false);
q2q_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),q2q,'uni',false);
q2q_cat = cat(3,q2q_pad{:});

% Same classification overlap shuffle
n_rep = 1000;

m2m_cat_shuff = nan(14,14,length(data),n_rep);
q2q_cat_shuff = nan(14,14,length(data),n_rep);

m_use = cellfun(@(x) x(any(x,2),:),m,'uni',false);
q_use = cellfun(@(x) x(any(x,2),:),q,'uni',false);

for i = 1:n_rep
    m_shuff = cellfun(@(x) shake(x,1),m_use,'uni',false);    
    m2m_shuff = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + ...
        repmat(sum(m,1),size(m,2),1)),m_shuff,'uni',false);  
    m2m_shuff_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),m2m_shuff,'uni',false);
    m2m_cat_shuff(:,:,:,i) = cat(3,m2m_shuff_pad{:});
    
    q_shuff = cellfun(@(x) shake(x,1),q_use,'uni',false);    
    q2q_shuff = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + ...
        repmat(sum(q,1),size(q,2),1)),q_shuff,'uni',false);   
    q2q_shuff_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),q2q_shuff,'uni',false);
    q2q_cat_shuff(:,:,:,i) = cat(3,q2q_shuff_pad{:});
 
    disp(i);
end

% Real relative to shuffle (in stds from mean)
m2m_cat_shuff_mean = nanmean(m2m_cat_shuff,4);
m2m_cat_shuff_std = nanstd(m2m_cat_shuff,[],4);
m2m_cat_stds = (m2m_cat-m2m_cat_shuff_mean)./m2m_cat_shuff_std;

q2q_cat_shuff_mean = nanmean(q2q_cat_shuff,4);
q2q_cat_shuff_std = nanstd(q2q_cat_shuff,[],4);
q2q_cat_stds = (q2q_cat-q2q_cat_shuff_mean)./q2q_cat_shuff_std;

% Average triangular values in first and second week
m2m_tril = cell2mat(arrayfun(@(x) [nanmean(AP_itril(m2m_cat_stds(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(m2m_cat_stds(8:end,8:end,x),-1))], ...
    transpose(1:size(m2m_cat_stds,3)),'uni',false));

q2q_tril = cell2mat(arrayfun(@(x) [nanmean(AP_itril(q2q_cat_stds(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(q2q_cat_stds(8:end,8:end,x),-1))], ...
    transpose(1:size(q2q_cat_stds,3)),'uni',false));

% Plot the grid
figure; 
subplot(1,3,1);
imagesc(nanmean(m2m_cat_stds,3));
axis square;
caxis([0 5]);
c = colorbar;
ylabel(c,'Z-score');
title('M2M');

subplot(1,3,2);
imagesc(nanmean(q2q_cat_stds,3));
axis square;
caxis([0 5]);
c = colorbar;
ylabel(c,'Z-score');
title('Q2Q');

% Plot the difference between wk1 and wk2 triangular pieces
m2m_plot_mean = nanmean(m2m_tril);
m2m_plot_sem = nanstd(m2m_tril)./sqrt(sum(~isnan(m2m_tril)));

q2q_plot_mean = nanmean(q2q_tril);
q2q_plot_sem = nanstd(q2q_tril)./sqrt(sum(~isnan(q2q_tril)));

subplot(1,3,3); hold on;
bar([m2m_plot_mean;q2q_plot_mean],'linewidth',2); colormap(gray)
errorb([m2m_plot_mean;q2q_plot_mean],[m2m_plot_sem;q2q_plot_sem],'linewidth',2);

set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'M2M','Q2Q'});
legend({'Week 1','Week 2'});
ylabel('Z-score')

% Stats
[p_m2m,~,stats_m2m] = signrank(m2m_tril(:,1),m2m_tril(:,2));
[p_q2q,~,stats_q2q] = signrank(q2q_tril(:,1),q2q_tril(:,2));


%% 5D: classification switching

clearvars -except data analysis classified_rois

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);

% Grid of overlap
m2q = cellfun(@(m,q) (+m'*+q)./(repmat(sum(m,1)',1,size(m,2)) + ...
    repmat(sum(q,1),size(q,2),1)),m_pad,q_pad,'uni',false);
m2q_cat = cat(3,m2q{:});

% Week 1 vs. week 2 classification switch
min_classdays = 1;
m2q = cellfun(@(m,q) ...
    sum((sum(m(:,1:7),2) >= min_classdays) & ...
    (sum(q(:,8:end),2) >= min_classdays))./ ...
    sum((sum(m(:,1:7),2) >= min_classdays)),m_pad,q_pad,'uni',false);
q2m = cellfun(@(m,q) ...
    sum((sum(q(:,1:7),2) >= min_classdays) & ...
    (sum(m(:,8:end),2) >= min_classdays))./ ...
    sum((sum(q(:,1:7),2) >= min_classdays)),m_pad,q_pad,'uni',false);

m2q_mean = cellfun(@nanmean,m2q);
q2m_mean = cellfun(@nanmean,q2m);

figure; 
% Plot overlap grid
subplot(1,2,1);
imagesc(nanmean(m2q_cat,3));colormap(gray);
xlabel('Quiescent');
ylabel('Movement');
c = colorbar;
ylabel(c,'Fraction of overlapping ROIs')
axis square;

% Plot fraction of ROIs classified in week 1 that switched week 2
subplot(1,2,2);hold on;
bar([nanmean(m2q_mean) nanmean(q2m_mean)],'FaceColor','w','linewidth',2);
errorbar([nanmean(m2q_mean) nanmean(q2m_mean)], ...
    [nanstd(m2q_mean)./sqrt(sum(~isnan(m2q_mean))), ...
    nanstd(q2m_mean)./sqrt(sum(~isnan(q2m_mean)))],'.k','linewidth',2);
ylabel('Fraction of ROIs')
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'M2Q','Q2M'});
xlim([0,3]);

% Stats
p = signrank(m2q_mean,q2m_mean);


%% 5E: M2Q ROI average activity rate

clearvars -except data analysis classified_rois

avg_movement_act_m2q = cell(length(data),1);
avg_quiescent_act_m2q = cell(length(data),1);
avg_total_act_m2q = cell(length(data),1);

for curr_animal = 1:length(data)
    
    curr_m2q = ...
        (cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent;
    
    use_cells = any(curr_m2q,2);
        
    curr_avg_movement_act = nan(length(use_cells),14);
    curr_avg_quiescent_act = nan(length(use_cells),14);
    curr_avg_total_act = nan(length(use_cells),14);
    
     for i = 1:14
        curr_avg_movement_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames),2);
        
        curr_avg_quiescent_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,~analysis(curr_animal).lever(i).lever_move_frames),2);
        
        curr_avg_total_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,:),2);
     end
    
    avg_movement_act_m2q{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act_m2q{curr_animal} = curr_avg_quiescent_act;
    avg_total_act_m2q{curr_animal} = curr_avg_total_act;
end

% Plot average activity during classified days
act_all = nan(length(data),6);
for curr_animal = 1:length(data)
    curr_movement_act_m = avg_movement_act_m2q{curr_animal};
    curr_movement_act_m(~classified_rois(curr_animal).movement) = NaN;
    
    curr_movement_act_q = avg_movement_act_m2q{curr_animal};
    curr_movement_act_q(~classified_rois(curr_animal).quiescent) = NaN;
    
    curr_quiescent_act_m = avg_quiescent_act_m2q{curr_animal};
    curr_quiescent_act_m(~classified_rois(curr_animal).movement) = NaN;
    
    curr_quiescent_act_q = avg_quiescent_act_m2q{curr_animal};
    curr_quiescent_act_q(~classified_rois(curr_animal).quiescent) = NaN;
    
    curr_total_act_m = avg_total_act_m2q{curr_animal};
    curr_total_act_m(~classified_rois(curr_animal).movement) = NaN;
    
    curr_total_act_q = avg_total_act_m2q{curr_animal};
    curr_total_act_q(~classified_rois(curr_animal).quiescent) = NaN;
    
    act_all(curr_animal,:) = [ ...
        nanmean(nanmean(curr_movement_act_m,2)), ...
        nanmean(nanmean(curr_movement_act_q,2)), ...
        nanmean(nanmean(curr_quiescent_act_m,2)), ...
        nanmean(nanmean(curr_quiescent_act_q,2)), ...
        nanmean(nanmean(curr_total_act_m,2)), ...
        nanmean(nanmean(curr_total_act_q,2))];
end

act_mean = reshape(nanmean(act_all),2,3)';
act_sem = reshape(nanstd(act_all)./sqrt(sum(~isnan(act_all))),2,3)';
figure;errorb(act_mean,act_sem);colormap(gray)
set(gca,'XTickLabel',{'Movement activity','Quiescent activity','Total activity'});
ylabel('Average activity');
legend({'Movement classified','Quiescent classified'});
title('M2Q cells');

% Stats
pm = signrank(act_all(:,1),act_all(:,2));
pq = signrank(act_all(:,3),act_all(:,4));
pt = signrank(act_all(:,5),act_all(:,6));



%% 6A: Activity/movement correlation of movements (within/across epochs)

clearvars -except data analysis classified_rois

session_epochs = {[1,2,3,4],[11,12,13,14]};

% Get activity and movements
epoch_movement = cell(length(data),length(session_epochs));
epoch_activity = cell(length(data),length(session_epochs));

for curr_animal = 1:length(data)
    disp(['Animal ' num2str(curr_animal)]);
    
    for curr_epoch = 1:length(session_epochs)
        
        use_sessions = session_epochs{curr_epoch};
        
        for curr_session_idx = 1:length(use_sessions)
            
            curr_session = session_epochs{curr_epoch}(curr_session_idx);
            % Skip if animal doesn't have current session
            if length(data(curr_animal).bhv) < curr_session || ...
                    isempty(data(curr_animal).bhv(curr_session).lever_force)
                continue
            end            
            
            % Get movement durations
            [lever_active,lever_force_resample] = AP_parseLeverMovement_updated( ...
                data(curr_animal).bhv(curr_session).lever_force);
            
            movement_starts = find(diff([0;lever_active;0]) == 1);
            movement_stops = find(diff([0;lever_active;0]) == -1);
            movement_durations = movement_stops - movement_starts;
            movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
            pre_movement_iti = [NaN;movement_iti];
            post_movement_iti = [movement_iti;NaN];
            
            % Set criteria for pulling out movements
            use_time = 2000;
            use_frames = (use_time/1000)*28;
            
            min_move_time = 2000;
            max_move_time = Inf;
            min_pre_move_iti = 1000;
            min_post_move_iti = 0;
            use_movements = ...
                movement_durations > min_move_time & ...
                movement_durations < max_move_time & ...
                pre_movement_iti > min_pre_move_iti & ...
                post_movement_iti > min_post_move_iti;
            
            use_movement_start_times = movement_starts(use_movements);
            [~,use_movement_start_frames] = arrayfun(@(x) ...
                min(abs(data(curr_animal).bhv(curr_session).frame_times - ...
                use_movement_start_times(x)/1000)),transpose(1:length(use_movement_start_times)));
            
            % Don't use any movements that don't have full timing in both
            % the lever and the imaging
                                  
            edge_movements = (use_movement_start_times + (use_time-1) > ...
                length(lever_force_resample)) | ...
                (use_movement_start_frames + (use_frames-1) > ...
                size(data(curr_animal).im(curr_session).roi_trace_df,2));
                                
            use_movement_start_times(edge_movements) = [];
            use_movement_start_frames(edge_movements) = [];
            
            movements_idx = repmat(use_movement_start_times,1,use_time) + ...
                repmat(0:use_time-1,length(use_movement_start_times),1);
            timed_movements = lever_force_resample(movements_idx');
            
            timed_activity = permute(cell2mat(arrayfun(@(x) ...
                data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
                use_movement_start_frames(x):use_movement_start_frames(x)+(use_frames-1)), ...
                permute(1:length(use_movement_start_frames),[3,1,2]),'uni',false)),[2,1,3]);
            
            epoch_movement{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_movements;
            epoch_activity{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_activity;
            
            disp(curr_session);
            
        end
    end  
end

% Replace any empty data with empty cells for uniform type
empty_data = cellfun(@(x) isempty(x),epoch_movement);
epoch_movement(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);
epoch_activity(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);

% Get trial-by-trial activity/movement correlation

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

epoch_movement_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_movement_cat{x,:})), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_movement_corr = cellfun(@(allcorr,act1) ...
    {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
    AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1), ...
    reshape(allcorr(1:size(act1,2),size(act1,2)+1:end),[],1)}, ...
    epoch_movement_corr_all,epoch_movement_cat(:,1)','uni',false);
epoch_movement_corr = vertcat(epoch_movement_corr{:});

epoch_activity_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_activity_cat_reshape{x,:}),'rows','pairwise'), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_activity_corr = cellfun(@(allcorr,act1) ...
    {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
    AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1), ...
    reshape(allcorr(1:size(act1,2),size(act1,2)+1:end),[],1)}, ...
    epoch_activity_corr_all,epoch_activity_cat_reshape(:,1)','uni',false);
epoch_activity_corr = vertcat(epoch_activity_corr{:});

% Bin trials by movement correlation
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use)-1;
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),epoch_movement_corr,'uni',false);

% Get mean/sem of activity correlations in each movement correlation bin
actcorr_bin_mean = cellfun(@(act,movebin) grpstats(act,movebin,'nanmean'), ...
    epoch_activity_corr,movecorr_bin_idx,'uni',false);
actcorr_bin_sem = cellfun(@(act,movebin) grpstats(act,movebin,'sem'), ...
    epoch_activity_corr,movecorr_bin_idx,'uni',false);

actcorr_bin_mean_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_mean,'uni',false);
actcorr_bin_sem_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_sem,'uni',false);

for i = 1:numel(actcorr_bin_mean)
    curr_used_bins = unique(movecorr_bin_idx{i});
    actcorr_bin_mean_pad{i}(curr_used_bins) = actcorr_bin_mean{i};
    actcorr_bin_sem_pad{i}(curr_used_bins) = actcorr_bin_sem{i};
end

epoch_combine_mean = nan(length(corr_bin),3);
epoch_combine_sem = nan(length(corr_bin),3);
for i = 1:3
    epoch_combine_mean(:,i) = nanmean(horzcat(actcorr_bin_mean_pad{:,i}),2);
    epoch_combine_sem(:,i) = nanstd(horzcat(actcorr_bin_mean_pad{:,i}),[],2)./ ...
        sqrt(sum(~isnan(horzcat(actcorr_bin_mean_pad{:,i})),2));
end

figure; hold on;
errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);
legend({'Within early','Within late','Across early/late'});
xlabel('Pairwise movement correlation');
ylabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);
xlim([-1,1]);
set(gca,'XTick',-1:0.5:1);

% Statistics

% Time within each group
corr_r = nan(3,1);
corr_p = nan(3,1);
anova_p = nan(3,1);
for i = 1:3
    curr_data = horzcat(actcorr_bin_mean_pad{:,i});
    curr_bins = repmat(transpose(1:length(corr_bin)),1,size(curr_data,2));
    
    curr_data = curr_data(1:end-1,:);
    curr_bins = curr_bins(1:end-1,:);
    
    %[r,p] = corrcoef(curr_bins,zscore(curr_data,[],1));
    [r,p] = corrcoef(curr_bins,curr_data);
    corr_r(i) = r(2);
    corr_p(i) = p(2);
    
    anova_p(i) = anova1(zscore(curr_data(4:end,:),[],1)',[],'off');
end

% Across early/late
early_data = horzcat(actcorr_bin_mean_pad{:,1});
early_data_use = early_data(1:end-1,:);
late_data = horzcat(actcorr_bin_mean_pad{:,2});
late_data_use = late_data(1:end-1,:);
across_data = horzcat(actcorr_bin_mean_pad{:,3});
across_data_use = across_data(1:end-1,:);

p_negcorr = signrank(reshape(early_data_use(1:3,:),[],1),reshape(late_data_use(1:3,:),[],1));
p_poscorr = signrank(reshape(early_data_use(4:6,:),[],1),reshape(late_data_use(4:6,:),[],1));

% Within early late, sign rank for reviewer 3?
p_early = signrank(early_data_use(1,:),early_data_use(end,:));
p_late = signrank(late_data_use(1,:),late_data_use(end,:));
p_across = signrank(across_data_use(1,:),across_data_use(end,:));

% two-way anova? 
[p,tbl,stats] = anova2([early_data_use(:),late_data_use(:),across_data_use(:)],length(data));

% New stats (reviewer 3)
early_data = horzcat(actcorr_bin_mean_pad{:,1});
early_data_use = early_data(1:end-1,:);
late_data = horzcat(actcorr_bin_mean_pad{:,2});
late_data_use = late_data(1:end-1,:);
across_data = horzcat(actcorr_bin_mean_pad{:,3});
across_data_use = across_data(1:end-1,:);
% (slope difference)
curr_fit = arrayfun(@(x) polyfit(1:n_bins,early_data_use(:,x)',1),1:length(data),'uni',false);
early_slope = cellfun(@(x) x(1),curr_fit);

curr_fit = arrayfun(@(x) polyfit(1:n_bins,late_data_use(:,x)',1),1:length(data),'uni',false);
late_slope = cellfun(@(x) x(1),curr_fit);

curr_fit = arrayfun(@(x) polyfit(1:n_bins,across_data_use(:,x)',1),1:length(data),'uni',false);
across_slope = cellfun(@(x) x(1),curr_fit);

p_slope_earlylate = signrank(early_slope,late_slope);
p_slope_earlyacross = signrank(early_slope,across_slope);
p_slope_lateacross = signrank(late_slope,across_slope);
% (plot segment difference)
p_negcorr = signrank(reshape(early_data_use(1:3,:),[],1),reshape(late_data_use(1:3,:),[],1));
p_poscorr = signrank(reshape(early_data_use(4:6,:),[],1),reshape(late_data_use(4:6,:),[],1));


%% 6B: (degeneracy) Pairwise activity within movement template correlation group

clearvars -except data analysis classified_rois

session_epochs = {[1,2,3,4],[11,12,13,14]};

% Get activity and movements
epoch_movement = cell(length(data),length(session_epochs));
epoch_activity = cell(length(data),length(session_epochs));

for curr_animal = 1:length(data)
    disp(['Animal ' num2str(curr_animal)]);
    
    for curr_epoch = 1:length(session_epochs)
        
        use_sessions = session_epochs{curr_epoch};
        
        for curr_session_idx = 1:length(use_sessions)
            
            curr_session = session_epochs{curr_epoch}(curr_session_idx);
            % Skip if animal doesn't have current session
            if length(data(curr_animal).bhv) < curr_session || ...
                    isempty(data(curr_animal).bhv(curr_session).lever_force)
                continue
            end            
            
            % Get movement durations
            [lever_active,lever_force_resample] = AP_parseLeverMovement_updated( ...
                data(curr_animal).bhv(curr_session).lever_force);
            
            movement_starts = find(diff([0;lever_active;0]) == 1);
            movement_stops = find(diff([0;lever_active;0]) == -1);
            movement_durations = movement_stops - movement_starts;
            movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
            pre_movement_iti = [NaN;movement_iti];
            post_movement_iti = [movement_iti;NaN];
            
            % Set criteria for pulling out movements
            use_time = 2000;
            use_frames = (use_time/1000)*28;
            
            min_move_time = 2000;
            max_move_time = Inf;
            min_pre_move_iti = 1000;
            min_post_move_iti = 0;
            use_movements = ...
                movement_durations > min_move_time & ...
                movement_durations < max_move_time & ...
                pre_movement_iti > min_pre_move_iti & ...
                post_movement_iti > min_post_move_iti;
            
            use_movement_start_times = movement_starts(use_movements);
            [~,use_movement_start_frames] = arrayfun(@(x) ...
                min(abs(data(curr_animal).bhv(curr_session).frame_times - ...
                use_movement_start_times(x)/1000)),transpose(1:length(use_movement_start_times)));
            
            % Don't use any movements that don't have full timing in both
            % the lever and the imaging
            edge_movements = (use_movement_start_times + (use_time-1) > ...
                length(lever_force_resample)) | ...
                (use_movement_start_frames + (use_frames-1) > ...
                size(data(curr_animal).im(curr_session).roi_trace_df,2));
                     
            use_movement_start_times(edge_movements) = [];
            use_movement_start_frames(edge_movements) = [];
            
            movements_idx = repmat(use_movement_start_times,1,use_time) + ...
                repmat(0:use_time-1,length(use_movement_start_times),1);
            timed_movements = lever_force_resample(movements_idx');
            
            timed_activity = permute(cell2mat(arrayfun(@(x) ...
                data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
                use_movement_start_frames(x):use_movement_start_frames(x)+(use_frames-1)), ...
                permute(1:length(use_movement_start_frames),[3,1,2]),'uni',false)),[2,1,3]);
            
            epoch_movement{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_movements;
            epoch_activity{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_activity;
            
            disp(curr_session);
            
        end
    end  
end

% Replace any empty data with empty cells for uniform type
empty_data = cellfun(@(x) isempty(x),epoch_movement);
epoch_movement(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);
epoch_activity(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

% Make template get trial-by-trial activity/movement correlation
corr_bin_use = [-Inf -0.5 0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin_valid = corr_bin_use;
corr_bin_valid(1) = -1;
corr_bin_valid(end) = 1;
corr_bin_plot = [corr_bin_valid(1:end-1) + diff(corr_bin_valid)/2 NaN];

all_act_corr = cell(size(epoch_movement_cat,1),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement from median of end sessions movements
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Bin trials by movement correlation
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);

actcorr_bin_mean = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
actcorr_bin_sem = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
movecorr_bin_mean = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
movecorr_bin_sem = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        for curr_bin = 1:n_bins
            
            use_trials =  ...
                find(movecorr_bin_idx{curr_animal,curr_epoch} == curr_bin);
            
%             % Skip if there aren't enough trials
%             if length(use_trials) < 5
%                 continue
%             end
            
            % Get pairwise correlation of activity within bins
            curr_act_corr = AP_itril(corrcoef( ...
                epoch_activity_cat_reshape{curr_animal,curr_epoch}(:,use_trials),'rows','pairwise'),-1);
            actcorr_bin_mean(curr_bin,curr_epoch,curr_animal) = nanmean(curr_act_corr);
            actcorr_bin_sem(curr_bin,curr_epoch,curr_animal) = nanstd(curr_act_corr)./sqrt(sum(~isnan(curr_act_corr)));
            
            % Get pairwise correlation of movement within bins, for use maybe in later sanity check
            curr_move_corr = AP_itril(corrcoef( ...
                epoch_movement_cat{curr_animal,curr_epoch}(:,use_trials),'rows','pairwise'),-1);
            movecorr_bin_mean(curr_bin,curr_epoch,curr_animal) = nanmean(curr_move_corr);
            movecorr_bin_sem(curr_bin,curr_epoch,curr_animal) = nanstd(curr_move_corr)./sqrt(sum(~isnan(curr_move_corr)));
            
            % Store all the activity correlations (total grouping?)
            all_act_corr{curr_animal,curr_epoch,curr_bin} = curr_act_corr;
        end
    end
end

epoch_combine_mean(:,:) = nanmean(actcorr_bin_mean,3);
epoch_combine_sem(:,:) = nanstd(actcorr_bin_mean,[],3)./sqrt(sum(~isnan(actcorr_bin_mean),3));

% few_animals = sum(~isnan(actcorr_bin_mean),3) < 2;
% epoch_combine_mean(few_animals) = NaN;
% epoch_combine_sem(few_animals) = NaN;

figure; 
plot_bins = [1,3];
% Plot pairwise movement correlations
subplot(2,1,1);hold on;
errorbar(repmat(corr_bin_plot(plot_bins)',[1,2]),nanmean(movecorr_bin_mean(plot_bins,:,:),3), ...
    nanmean(movecorr_bin_sem(plot_bins,:,:),3),'linewidth',2);
xlabel('Movement template correlation');
ylabel('Pairwise movement correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));

% Plot pairwise activity correlations
subplot(2,1,2);hold on;
errorbar(repmat(corr_bin_plot(plot_bins)',[1,2]),nanmean(epoch_combine_mean(plot_bins,:,:),3), ...
    nanmean(epoch_combine_sem(plot_bins,:,:),3),'linewidth',2);
xlabel('Movement template correlation');
ylabel('Pairwise activity correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));


% Statistics
p_early = signrank(squeeze(actcorr_bin_mean(plot_bins(1),1,:)),squeeze(actcorr_bin_mean(plot_bins(2),1,:)));
p_late = signrank(squeeze(actcorr_bin_mean(plot_bins(1),2,:)),squeeze(actcorr_bin_mean(plot_bins(2),2,:)));

p_early = signrank(squeeze(movecorr_bin_mean(plot_bins(1),1,:)),squeeze(movecorr_bin_mean(plot_bins(2),1,:)));
p_late = signrank(squeeze(movecorr_bin_mean(plot_bins(1),2,:)),squeeze(movecorr_bin_mean(plot_bins(2),2,:)));


%% 6C: (separability) pairwise movecorr cutoff AND template corr diff cutoff

clearvars -except data analysis classified_rois

session_epochs = {[1,2,3,4],[11,12,13,14]};

% Get activity and movements
epoch_movement = cell(length(data),length(session_epochs));
epoch_activity = cell(length(data),length(session_epochs));

for curr_animal = 1:length(data)
    disp(['Animal ' num2str(curr_animal)]);
    
    for curr_epoch = 1:length(session_epochs)
        
        use_sessions = session_epochs{curr_epoch};
        
        for curr_session_idx = 1:length(use_sessions)
            
            curr_session = session_epochs{curr_epoch}(curr_session_idx);
            % Skip if animal doesn't have current session
            if length(data(curr_animal).bhv) < curr_session || ...
                    isempty(data(curr_animal).bhv(curr_session).lever_force)
                continue
            end            
            
            % Get movement durations
            [lever_active,lever_force_resample] = AP_parseLeverMovement_updated( ...
                data(curr_animal).bhv(curr_session).lever_force);
            
            movement_starts = find(diff([0;lever_active;0]) == 1);
            movement_stops = find(diff([0;lever_active;0]) == -1);
            movement_durations = movement_stops - movement_starts;
            movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
            pre_movement_iti = [NaN;movement_iti];
            post_movement_iti = [movement_iti;NaN];
            
            % Set criteria for pulling out movements
            use_time = 2000;
            use_frames = (use_time/1000)*28;
            
            min_move_time = 2000;
            max_move_time = Inf;
            min_pre_move_iti = 1000;
            min_post_move_iti = 0;
            use_movements = ...
                movement_durations > min_move_time & ...
                movement_durations < max_move_time & ...
                pre_movement_iti > min_pre_move_iti & ...
                post_movement_iti > min_post_move_iti;
            
            use_movement_start_times = movement_starts(use_movements);
            [~,use_movement_start_frames] = arrayfun(@(x) ...
                min(abs(data(curr_animal).bhv(curr_session).frame_times - ...
                use_movement_start_times(x)/1000)),transpose(1:length(use_movement_start_times)));
            
            % Don't use any movements that don't have full timing in both
            % the lever and the imaging
            edge_movements = (use_movement_start_times + (use_time-1) > ...
                length(lever_force_resample)) | ...
                (use_movement_start_frames + (use_frames-1) > ...
                size(data(curr_animal).im(curr_session).roi_trace_df,2));
                     
            use_movement_start_times(edge_movements) = [];
            use_movement_start_frames(edge_movements) = [];
            
            movements_idx = repmat(use_movement_start_times,1,use_time) + ...
                repmat(0:use_time-1,length(use_movement_start_times),1);
            timed_movements = lever_force_resample(movements_idx');
            
            timed_activity = permute(cell2mat(arrayfun(@(x) ...
                data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
                use_movement_start_frames(x):use_movement_start_frames(x)+(use_frames-1)), ...
                permute(1:length(use_movement_start_frames),[3,1,2]),'uni',false)),[2,1,3]);
            
            epoch_movement{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_movements;
            epoch_activity{curr_animal,curr_epoch}{curr_session_idx} = ...
                timed_activity;
            
            disp(curr_session);
            
        end
    end  
end

% Replace any empty data with empty cells for uniform type
empty_data = cellfun(@(x) isempty(x),epoch_movement);
epoch_movement(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);
epoch_activity(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

% Make template movement
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Get pairwise correlations of all trials
movement_pairwise_corr = cellfun(@(x) corrcoef(x),epoch_movement_cat,'uni',false);

% Get differences in pairwise correlation with template
movement_template_corr_diff = cellfun(@(x) abs(bsxfun(@minus,x,x')), ...
    movement_template_corr,'uni',false);

% Select trials by pairwise correlation AND template correlation
movecorr_cutoff = -0.4;
movetemp_cutoff = 1.2;
ortho_movetemp_cutoff = 0.2;

template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movetempdiff > movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

ortho_template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movecorr < movecorr_cutoff & ...
    movetempdiff < ortho_movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        template_dim_pairs = template_axis_use{curr_animal,curr_epoch};
        ortho_template_dim_pairs = ortho_template_axis_use{curr_animal,curr_epoch};
        
        curr_move_corr_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch},'rows','pairwise');
        curr_act_corr_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch},'rows','pairwise');      
        
        % Get pairwise movement correlations along template edges
        move_corr_templatedim_mean = nanmean(curr_move_corr_grid(template_dim_pairs));
        
        % Get pairwise movement correlations in selected middle bin trials
        move_corr_orthotemplatedim_mean = nanmean(curr_move_corr_grid(ortho_template_dim_pairs));
        
        % Get pairwise activity correlations along template edges
        act_corr_templatedim_mean = nanmean(curr_act_corr_grid(template_dim_pairs));
        
        % Get pairwise activity correlations in selected middle bin trials
        act_corr_orthotemplatedim_mean = nanmean(curr_act_corr_grid(ortho_template_dim_pairs));
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_orthotemplatedim_mean];
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_orthotemplatedim_mean];
        
    end
end

movecorr_movedim_mean = nanmean(movecorr_movedim,3);
movecorr_movedim_sem = nanstd(movecorr_movedim,[],3)./sqrt(sum(~isnan(movecorr_movedim),3));

actcorr_movedim_mean = nanmean(actcorr_movedim,3);
actcorr_movedim_sem = nanstd(actcorr_movedim,[],3)./sqrt(sum(~isnan(actcorr_movedim),3));

figure;

subplot(2,1,1);
errorbar(movecorr_movedim_mean,movecorr_movedim_sem,'linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Template dimension','Ortho-template dimension'});
ylabel('Pairwise movement correlations');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));

subplot(2,1,2);
errorbar(actcorr_movedim_mean,actcorr_movedim_sem,'linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Template dimension','Ortho-template dimension'});
ylabel('Pairwise activity correlations');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));


% Statistics
p_early = signrank(squeeze(actcorr_movedim(1,1,:)),squeeze(actcorr_movedim(2,1,:)));
p_late = signrank(squeeze(actcorr_movedim(1,2,:)),squeeze(actcorr_movedim(2,2,:)));

p_early = signrank(squeeze(movecorr_movedim(1,1,:)),squeeze(movecorr_movedim(2,1,:)));
p_late = signrank(squeeze(movecorr_movedim(1,2,:)),squeeze(movecorr_movedim(2,2,:)));




