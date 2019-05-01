%% Smooth df trace

smoothed_roi_df = cell(14,1);
for i = 1:14
    
    roi_trace_df_smooth = nan(size(data.im(i).roi_trace_df));
    for curr_roi = 1:size(roi_trace_df_smooth,1);
        roi_trace_df_smooth(curr_roi,:) = smooth(data.im(i).roi_trace_df(curr_roi,:),30,'loess');
    end
    
    smoothed_roi_df{i} = roi_trace_df_smooth;
    
    disp(i)
    
end


%% Compare raw trace to interpolation subtracted background


curr_data = load('/usr/local/lab/People/Andy/Data/AP140/AP140_batch_thresh_roi/150719_AP140_001_001_summed_50_analysis.mat');

curr_df = AP_baselineEstimation(curr_data.im.roi_trace,30);
curr_interpsub_df = curr_data.im.roi_trace_df;

n_rois = size(curr_df,1);

df_corr = zeros(size(curr_df,1),1);
for i = 1:n_rois
    curr_corr = corrcoef(curr_interpsub_df(i,:),curr_df(i,:));
    df_corr(i) = curr_corr(2);
end

[~,sort_idx] = sort(df_corr);


%% Classify modulated ROIs (the old way)
% NOTE: requires loaded/prepared data
% this is copied into its own script: AP_classify_movement_cells_continuous

classified_rois = struct('movement',cell(length(data),1),'quiescent',cell(length(data),1));
for curr_animal = 1:length(data);
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    n_sessions = length(data(curr_animal).im);
    
    classified_rois(curr_animal).movement = false(n_rois,n_sessions);
    classified_rois(curr_animal).quiescent = false(n_rois,n_sessions);
    
    for curr_session = 1:n_sessions
             
        lever_active_frames = analysis(curr_animal).lever(curr_session).lever_move_frames;
        num_frames = size(data(curr_animal).im(curr_session).roi_trace_df,2);
        
        % Split active frames by movement blocks / quiescent frames
        % (to split up by movement blocks / quiescent blocks)
        %lever_active_frames_nan = +lever_active_frames;
        %lever_active_frames_nan(~lever_active_frames) = NaN;
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
        
        % Get activity during movement/quiescence
        movement_activity = +(data(curr_animal).im(curr_session).roi_trace_thresh(:,1:num_frames) > 0)*lever_active_frames;

        % Get shuffled activity distribution
        num_rep = 10000;
        shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
        lever_active_shuffle = nan(length(lever_active_frames),num_rep);
        for i = 1:num_rep
           lever_active_shuffle(:,i) = ...
               vertcat(lever_active_frames_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        shuffle_movement_activity = +(data(curr_animal).im(curr_session).roi_trace_thresh(:,1:num_frames) > 0)*lever_active_shuffle;
        clear lever_active_shuffle

        movement_rank = tiedrank([movement_activity shuffle_movement_activity]')';
        movement_p = movement_rank(:,1)/(num_rep+1);
        movement_cells = movement_p > 0.975;
        quiescent_cells = movement_p < 0.025;
        
        classified_rois(curr_animal).movement(:,curr_session) = movement_cells;
        classified_rois(curr_animal).quiescent(:,curr_session) = quiescent_cells;
        
        disp(['Classified session ' num2str(curr_session)]);
    end
    
    disp(['Classified ' data(curr_animal).animal]);
end
disp('Finished all')


%% Plot PCA of average activity (movement, quiescent, all)
% requires sigcells

curr_movecells = any(classified_rois(1).movement,2);
curr_quiescentcells = any(classified_rois(1).quiescent,2);

avg_act = nan(size(analysis.im(1).move_onset_aligned,3), ...
    size(analysis.im(1).move_onset_aligned,2),14);
for i = 1:14
    avg_act(:,:,i) = permute(nanmean(analysis.im(i).move_onset_aligned,1),[3 2 1]); 
end

use_cells = curr_movecells;
avg_act_reshape = reshape(permute(avg_act(use_cells,:,:),[2 3 1]),size(avg_act,2)*size(avg_act,3),[]);
[coeff,score,latent] = pca(avg_act_reshape);
avg_act_score = mat2cell(score,repmat(size(avg_act,2),size(avg_act,3),1),size(score,2));

col = jet(14);
figure; hold on
for i = 1:14
plot3(avg_act_score{i}(:,2),avg_act_score{i}(:,3),avg_act_score{i}(:,4), ...
    'color',col(i,:),'linewidth',2)
end


%% Average activity correlation

act_aligned_corr_move = nan(14,14,length(data));
act_total_corr_move = nan(14,14,length(data));

act_aligned_corr_quiescent = nan(14,14,length(data));
act_total_corr_quiescent = nan(14,14,length(data));

for curr_animal = 1:length(data);
    
    % Correlation of aligned movement activity
    % NOTE: this looks the same regardless of alignment...
    curr_movecells = any(classified_rois(curr_animal).movement,2);
    curr_quiescentcells = any(classified_rois(curr_animal).quiescent,2);
    
    use_frames = analysis(curr_animal).surrounding_frames >= 0 &  ...
        analysis(curr_animal).surrounding_frames <= 90;
    
    % MOVEMENT CELLS
    use_cells = find(curr_movecells);
     
    avg_aligned_act = nan(length(use_cells),sum(use_frames),14);
    for i = 1:14
        avg_aligned_act(:,:,i) = permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned(:,use_frames,use_cells),1),[3 2 1]);
    end
    
    act_aligned_reshape = reshape(permute(avg_aligned_act,[2 1 3]),[],size(avg_aligned_act,3));
    act_aligned_corr_move(:,:,curr_animal) = corrcoef(act_aligned_reshape);
    
    % Correlation of average activity
    avg_total_act = nan(length(use_cells),14);
    for i = 1:14
        avg_total_act(:,i) = nanmean(data(curr_animal).im(i).roi_trace_thresh(use_cells,:),2);
    end
    act_total_corr_move(:,:,curr_animal) = corrcoef(avg_total_act);  
    
    % QUIESCENT CELLS
    use_cells = find(curr_quiescentcells);
    
    avg_aligned_act = nan(length(use_cells),sum(use_frames),14);
    for i = 1:14
        avg_aligned_act(:,:,i) = permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned(:,use_frames,use_cells),1),[3 2 1]);
    end
    
    act_aligned_reshape = reshape(permute(avg_aligned_act,[2 1 3]),[],size(avg_aligned_act,3));
    act_aligned_corr_quiescent(:,:,curr_animal) = corrcoef(act_aligned_reshape);
    
    % Correlation of average activity
    avg_total_act = nan(length(use_cells),14);
    for i = 1:14
        avg_total_act(:,i) = nanmean(data(curr_animal).im(i).roi_trace_thresh(use_cells,:),2);
    end
    act_total_corr_quiescent(:,:,curr_animal) = corrcoef(avg_total_act);      
end

figure
subplot(2,2,1);
imagesc(nanmean(act_aligned_corr_move,3));colormap(hot)
xlabel('Day')
ylabel('Day')
title('Average aligned activity correlation: move cells')
subplot(2,2,2);
imagesc(nanmean(act_total_corr_move,3));colormap(hot)
xlabel('Day')
ylabel('Day')
title('Average total activity correlation: move cells')
subplot(2,2,3);
imagesc(nanmean(act_aligned_corr_quiescent,3));colormap(hot)
xlabel('Day')
ylabel('Day')
title('Average aligned activity correlation: quiescent cells')
subplot(2,2,4);
imagesc(nanmean(act_total_corr_quiescent,3));colormap(hot)
xlabel('Day')
ylabel('Day')
title('Average total activity correlation: quiescent cells')


%% Activity during quiescence / movement

avg_total_act_movement_corr_movecells = nan(14,14,length(data));
avg_total_act_quiescent_corr_movecells = nan(14,14,length(data));

avg_total_act_movement_corr_quiescentcells = nan(14,14,length(data));
avg_total_act_quiescent_corr_quiescentcells = nan(14,14,length(data));

for curr_animal = 1:length(data);
    curr_movecells = any(classified_rois(curr_animal).movement,2);
    curr_quiescentcells = any(classified_rois(curr_animal).quiescent,2);
    
    % MOVEMENT CELLS
    use_cells = find(curr_movecells);
    
    % Correlation of average activity
    avg_total_act_movement = nan(length(use_cells),14);
    avg_total_act_quiescent = nan(length(use_cells),14);
    
    for i = 1:14
        avg_total_act_movement(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames),2);
        avg_total_act_quiescent(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,~analysis(curr_animal).lever(i).lever_move_frames),2);
    end
    
    avg_total_act_movement_corr_movecells(:,:,curr_animal) = corrcoef(avg_total_act_movement);
    avg_total_act_quiescent_corr_movecells(:,:,curr_animal) = corrcoef(avg_total_act_quiescent);
    
    % QUIESCENT CELLS
    use_cells = find(curr_quiescentcells);
    
    % Correlation of average activity
    avg_total_act_movement = nan(length(use_cells),14);
    avg_total_act_quiescent = nan(length(use_cells),14);
    
    for i = 1:14
        avg_total_act_movement(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames),2);
        avg_total_act_quiescent(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,~analysis(curr_animal).lever(i).lever_move_frames),2);
    end
    
    avg_total_act_movement_corr_quiescentcells(:,:,curr_animal) = corrcoef(avg_total_act_movement);
    avg_total_act_quiescent_corr_quiescentcells(:,:,curr_animal) = corrcoef(avg_total_act_quiescent);
end

figure; 
subplot(2,2,1);
imagesc(nanmean(avg_total_act_movement_corr_movecells,3));colormap(hot);
xlabel('Day')
ylabel('Day')
title('Average movement activity correlation: move cells')
subplot(2,2,2);
imagesc(nanmean(avg_total_act_quiescent_corr_movecells,3));colormap(hot);
xlabel('Day')
ylabel('Day')
title('Average quiescent activity correlation: move cells')
subplot(2,2,3);
imagesc(nanmean(avg_total_act_movement_corr_quiescentcells,3));colormap(hot);
xlabel('Day')
ylabel('Day')
title('Average movement activity correlation: quiescent cells')
subplot(2,2,4);
imagesc(nanmean(avg_total_act_quiescent_corr_quiescentcells,3));colormap(hot);
xlabel('Day')
ylabel('Day')
title('Average quiescent activity correlation: quiescent cells')



%% Activity timing: std of activity onset from movement onset

act_std_all = cell(1,length(data));
for curr_animal = 1:length(data);
    
    use_frames = analysis(curr_animal).surrounding_frames >= 0 & analysis(curr_animal).surrounding_frames <= 90;

    act_starts = cell(14,size(data(curr_animal).im(1).roi_trace_df,1));
    
    for i = 1:14
        [curr_max,curr_max_idx] = max(analysis(curr_animal).im(i).move_onset_aligned(:,use_frames,:) > 0,[],2);
        curr_act_starts = arrayfun(@(x) curr_max_idx(curr_max(:,:,x),:,x),1:size(curr_max,3),'uni',false);
        act_starts(i,:) = curr_act_starts';
    end
    
    act_std = cellfun(@nanstd,act_starts);
    num_act = cellfun(@length,act_starts);
    min_num_act = 5;
    act_std(num_act < min_num_act) = NaN;
    
    % Get activity of only movement_related cells
    act_std(~classified_rois(curr_animal).movement') = NaN;
    
    act_std_all{curr_animal} = act_std;
end

figure;
framerate = 29; % just estimating for now
act_std_mean_cat = cell2mat(cellfun(@(x) nanmean(x,2)./framerate,act_std_all,'uni',false));
errorbar(nanmean(act_std_mean_cat,2),nanstd(act_std_mean_cat,[],2)./sqrt(sum(~isnan(act_std_mean_cat),2)),'k','linewidth',2);
xlabel('Day');
ylabel('Standard devation of movement-cell activity timing (s)')

%% Get movement aligned to activity onset
n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);

% Raw lever trace
back_time = 3000;
forward_time = 3000;
act_aligned_move = nan(n_rois,back_time+forward_time+1,14);

for curr_session = 1:14;  
    % downsample lever
    downsample_factor = 100;
    lever_force_resample = resample(data(curr_animal).bhv(curr_session).lever_force,1,downsample_factor);
    for curr_cell = 1:n_rois
        curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
        curr_act_start_time = round(data(curr_animal).bhv(curr_session).frame_times(curr_act_starts)*(10000/downsample_factor));
        curr_act_start_time(curr_act_start_time-back_time <= 0 | ...
            curr_act_start_time + forward_time > length(lever_force_resample)) = [];
        
        curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_start_time(x)-back_time: ...
            curr_act_start_time(x)+forward_time,(1:length(curr_act_start_time))','uni',false));
        
        % Lever trace
        curr_act_move = lever_force_resample(curr_act_move_idx);
        
        act_aligned_move(curr_cell,:,curr_session) = nanmean(curr_act_move);
    end
end


% Movement-binarized lever
back_time = 29*30;
forward_time = 29*30;
act_aligned_move_binary = nan(n_rois,back_time+forward_time+1,14);

for curr_session = 1:14;
    for curr_cell = 1:n_rois
        curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
        curr_act_starts(curr_act_starts-back_time <= 0 | ...
            curr_act_starts + forward_time > length(analysis(curr_animal).lever(curr_session).lever_move_frames)) = [];
        
        curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_starts(x)-back_time: ...
            curr_act_starts(x)+forward_time,(1:length(curr_act_starts))','uni',false));
        
        % Movement-binarzed lever
        curr_act_move = analysis(curr_animal).lever(curr_session).lever_move_frames(curr_act_move_idx);
        % Movement-binarized onset/offset
        %lever_diff = smooth(diff([0;analysis(curr_animal).lever(curr_session).lever_move_frames]),10);
        %curr_act_move = lever_diff(curr_act_move_idx);
        
        act_aligned_move_binary(curr_cell,:,curr_session) = nanmean(curr_act_move);
    end
    disp(curr_session);
end

% Get standard deviation of maximum movement probability timing
[~,max_idx] = arrayfun(@(x) max(act_aligned_move_binary(classified_rois(curr_animal).movement(:,x) > 0,:,x),[],2),1:14,'uni',false);
move_timing_std = cellfun(@nanstd,max_idx);


%% Plot fraction of movement/quiescence ROIs

move_frac = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false));
quiescent_frac = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false));
figure; hold on
errorbar(nanmean(move_frac),nanstd(move_frac)./sqrt(sum(~isnan(move_frac))),'k','linewidth',2)
errorbar(nanmean(quiescent_frac),nanstd(quiescent_frac)./sqrt(sum(~isnan(quiescent_frac))),'r','linewidth',2)
ylabel('Fraction of ROIs');
xlabel('Day');
legend({'Movement-related' 'Quiescence-related'})

mq_ratio = move_frac./quiescent_frac;
figure;
errorbar(nanmean(mq_ratio),nanstd(mq_ratio)./sqrt(sum(~isnan(mq_ratio))),'k','linewidth',2);
ylabel('Movement/Quiescent fraction')
xlabel('Day')



%% Compare movement correlation vs. time within day
% Activity looks like it happens in blocks for cells: if trials are sorted
% according to proximity in time vs. proximity in correlation, which one
% makes the difference within a day?

framerate = 28;

% Set up bins for movement correlation and timing difference
n_movecorr_bins = 20;
n_movetime_bins = 20;
max_frames = NaN;
for curr_animal = 1:length(data)
    curr_max_frames = max(cellfun(@(x) size(x,2),{data(curr_animal).im(:).roi_trace_df}));
    max_frames = max(max_frames,curr_max_frames);
end
movecorr_edges = linspace(-1,1,n_movecorr_bins);
movecorr_edges(1) = -inf;movecorr_edges(end) = inf;
movetime_edges = linspace(1,max_frames,n_movetime_bins);
movetime_edges(1) = -inf;movetime_edges(end) = inf;

actcorr_movebin = cell(length(analysis),1);
for curr_animal = 1:length(analysis)
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    movement_corr = cell(n_sessions,1);
    movement_corr(use_sessions) = cellfun(@(movement,trials) AP_itril(corrcoef(cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials),'uni',false)')),-1), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    % Trial-by-trial movement start time difference
    movement_timediff = cell(n_sessions,1);
    movement_timediff(use_sessions) = cellfun(@(movement,trials) AP_itril(pdist2( ...
        movement(trials),movement(trials)),-1), ...
        {analysis(curr_animal).lever(use_sessions).move_onset_frames}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    use_cells = true(size(classified_rois(curr_animal).movement));
    
    activity_corr = cell(n_sessions,1);
    for curr_session = find(use_sessions)
       curr_act = analysis(curr_animal).im(curr_session). ...
           move_onset_aligned(:,use_frames,use_cells(:,curr_session));
       activity_corr{curr_session} = AP_itril(corrcoef(reshape(permute( ...
           curr_act,[2,3,1]),[],size(curr_act,1)),'rows','pairwise'),-1);
    end
    
    % Bin trials by time in trace AND movement correlation   
    actcorr_movebin{curr_animal} = nan(n_movecorr_bins,n_movetime_bins,length(use_sessions));
    for curr_session = find(use_sessions)       
        [~,movecorr_bin] = histc(movement_corr{curr_session},movecorr_edges);
        [~,movetime_bin] = histc(movement_timediff{curr_session},movetime_edges);
        
        [actcorr_bin_mean bin_coords] = grpstats(activity_corr{curr_session}, ...
            {movecorr_bin movetime_bin},{'mean','gname'});
        
        bin_coords = cellfun(@(x) str2num(x),bin_coords);
        
        for i = 1:length(actcorr_bin_mean);
            actcorr_movebin{curr_animal}(bin_coords(i,1),bin_coords(i,2),curr_session) = ...
                actcorr_bin_mean(i);
        end
    end
    
end

% Fit each day to linear combination of movement time and correlation
a = cat(4,actcorr_movebin{:});
b = nanmean(a,4);
coeff_vals = nan(14,3);
coeff_gof = nan(14,1);
for i = 1:14
    c = b(7:end-1,1:end-7,i);
    [x,y] = meshgrid(1:length(c));
    nan_vals = isnan(c(:));
    
    x(nan_vals) = [];
    y(nan_vals) = [];
    c(nan_vals) = [];
            
    [sf,gof] = fit([x(:) y(:)],c(:),'poly11');
    coeff_vals(i,:) = coeffvalues(sf);
    coeff_gof(i) = gof.sse;
end
% Plot coefficient values of each component
figure;
subplot(2,1,1);
plot(coeff_vals(:,2),'k','linewidth',2);
title('Movement time difference coefficient')
xlabel('Day')
ylabel('Coefficient value');
subplot(2,1,2);
plot(coeff_vals(:,3),'k','linewidth',2);
title('Movement correlation coefficient')
xlabel('Day')
ylabel('Coefficient value');


% Compare trials on non-within days
% if time matters here, then it's TIME WITHIN SESSION dependent??
actcorr_movebin_crossday = cell(length(analysis),1);

for curr_animal = 1:length(analysis)
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    day_trials = cellfun(@sum,{analysis(curr_animal).lever(use_sessions).cued_movement_trials});
    move_day = cell2mat(cellfun(@(x,y) repmat(x,y,1),num2cell(1:14)',num2cell(day_trials)','uni',false));
    exclude_comparisons = AP_itril(pdist2(move_day,move_day) == 0,-1);
    
    % Concat all movement parameters
    curr_move = {analysis(curr_animal).lever(use_sessions).rewarded_movement}';
    curr_move = vertcat(curr_move{:});
    
    curr_trials = {analysis(curr_animal).lever(use_sessions).cued_movement_trials}';
    curr_trials = vertcat(curr_trials{:});
    
    movement_corr = {};
    movement_corr{1} = AP_itril(corrcoef(cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),curr_move(curr_trials),'uni',false)')),-1);
    movement_corr{1}(exclude_comparisons) = [];
    
    curr_movetime = {analysis(curr_animal).lever(use_sessions).move_onset_frames}';
    curr_movetime = vertcat(curr_movetime{:});
    
    movement_timediff = {};
    movement_timediff{1} = AP_itril(pdist2( ...
        curr_movetime(curr_trials),curr_movetime(curr_trials)),-1);
    movement_timediff{1}(exclude_comparisons) = [];
    
    % Concat all activity patterns
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    curr_act = {analysis(curr_animal).im(use_sessions).move_onset_aligned};
    curr_act = vertcat(curr_act{:});
    
    activity_corr = {};
    activity_corr{1} = AP_itril(corrcoef(reshape(permute( ...
        curr_act(:,use_frames,use_cells), ...
        [2 3 1]),[],size(curr_act,1)),'rows','pairwise'),-1);
    activity_corr{1}(exclude_comparisons) = [];
    
    % Bin trials by time in trace AND movement correlation
    actcorr_movebin_crossday{curr_animal} = nan(n_movecorr_bins,n_movetime_bins);
    [~,movecorr_bin] = histc(movement_corr{1},movecorr_edges);
    [~,movetime_bin] = histc(movement_timediff{1},movetime_edges);
    
    [actcorr_bin_mean bin_coords] = grpstats(activity_corr{1}, ...
        {movecorr_bin movetime_bin},{'mean','gname'});
    
    bin_coords = cellfun(@(x) str2num(x),bin_coords);
    
    for i = 1:length(actcorr_bin_mean);
        actcorr_movebin_crossday{curr_animal}(bin_coords(i,1),bin_coords(i,2)) = ...
            actcorr_bin_mean(i);
    end
    
end


% Compare trials on n-separated days
actcorr_movebin_ncrossday = cell(length(analysis),1);

for curr_animal = 1:length(analysis)
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    day_trials = cellfun(@sum,{analysis(curr_animal).lever(use_sessions).cued_movement_trials});
    move_day = cell2mat(cellfun(@(x,y) repmat(x,y,1),num2cell(1:14)',num2cell(day_trials)','uni',false));
    
    % Concat all movement parameters
    curr_move = {analysis(curr_animal).lever(use_sessions).rewarded_movement}';
    curr_move = vertcat(curr_move{:});
    
    curr_trials = {analysis(curr_animal).lever(use_sessions).cued_movement_trials}';
    curr_trials = vertcat(curr_trials{:});
    
    movement_corr = {};
    movement_corr{1} = AP_itril(corrcoef(cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),curr_move(curr_trials),'uni',false)')),-1);
     
    curr_movetime = {analysis(curr_animal).lever(use_sessions).move_onset_frames}';
    curr_movetime = vertcat(curr_movetime{:});
    
    movement_timediff = {};
    movement_timediff{1} = AP_itril(pdist2( ...
        curr_movetime(curr_trials),curr_movetime(curr_trials)),-1);    
    
    % Concat all activity patterns
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    curr_act = {analysis(curr_animal).im(use_sessions).move_onset_aligned};
    curr_act = vertcat(curr_act{:});
    
    activity_corr = {};
    activity_corr{1} = AP_itril(corrcoef(reshape(permute( ...
        curr_act(:,use_frames,use_cells), ...
        [2 3 1]),[],size(curr_act,1)),'rows','pairwise'),-1);
    
    
    % Bin trials by time in trace AND movement correlation
    actcorr_movebin_ncrossday{curr_animal} = nan(n_movecorr_bins,n_movetime_bins,13);
    for curr_daycomp = 1:13
        exclude_comparisons = AP_itril(pdist2(move_day,move_day) ~= curr_daycomp,-1);
        
        curr_movement_corr =  movement_corr{1}(~exclude_comparisons);
        curr_movement_timediff = movement_timediff{1}(~exclude_comparisons);
        curr_activity_corr = activity_corr{1}(~exclude_comparisons);
        
        [~,movecorr_bin] = histc(curr_movement_corr,movecorr_edges);
        [~,movetime_bin] = histc(curr_movement_timediff,movetime_edges);
        
        [actcorr_bin_mean bin_coords] = grpstats(curr_activity_corr, ...
            {movecorr_bin movetime_bin},{'mean','gname'});
        
        bin_coords = cellfun(@(x) str2num(x),bin_coords);
        
        for i = 1:length(actcorr_bin_mean);
            actcorr_movebin_ncrossday{curr_animal}(bin_coords(i,1),bin_coords(i,2),curr_daycomp) = ...
                actcorr_bin_mean(i);
        end
    end
end


% Fit each day to linear combination of movement time and correlation
a = cat(4,actcorr_movebin_ncrossday{:});
b = nanmean(a,4);
coeff_vals = nan(14,3);
coeff_gof = nan(14,1);
for i = 1:size(b,3)
    c = b(7:end-1,1:end-7,i);
    [x,y] = meshgrid(1:length(c));
    nan_vals = isnan(c(:));
    
    x(nan_vals) = [];
    y(nan_vals) = [];
    c(nan_vals) = [];
            
    [sf,gof] = fit([x(:) y(:)],c(:),'poly11');
    coeff_vals(i,:) = coeffvalues(sf);
    coeff_gof(i) = gof.sse;
end
% Plot coefficient values of each component
figure;
subplot(2,1,1);
plot(coeff_vals(:,2),'k','linewidth',2);
title('Movement time difference coefficient')
xlabel('Days separated')
ylabel('Coefficient value');
subplot(2,1,2);
plot(coeff_vals(:,3),'k','linewidth',2);
title('Movement correlation coefficient')
xlabel('Days separated')
ylabel('Coefficient value');


%% Correlate activity and movement correlations in each animal
% requires both behavior analysis and activity correlation to be run

% WORKING ON THIS: definite patterns present in both movement_corr_grid and
% activity i.e. avg_total_act_movement_corrcells, figure out what it is


figure;
for curr_animal = 1:length(data)
    subplot(length(data),2,(curr_animal-1)*2+1);
    imagesc(movement_corr_grid(:,:,curr_animal));
    colormap(gray);
    axis off
    
    if curr_animal == 1
        title('Movement corr grid')
    end
    
    subplot(length(data),2,(curr_animal-1)*2+2);
    imagesc(avg_total_act_movement_corr_movecells(:,:,curr_animal));
    colormap(gray);   
    axis off
    
    if curr_animal == 1
        title('Avg movement act movecells corr')
    end
    
end


% Manually enter days for new activity blocks
block_idx = [10,4,8,9,3,6,4,5];
act_pad = padarray(act_aligned_corr_move,[14,14,0],NaN,'pre');
move_pad = padarray(movement_corr_grid,[14,14,0],NaN,'pre');

for i = 1:8
   act_pad(:,:,i) = circshift(act_pad(:,:,i),[-block_idx(i),-block_idx(i)]);  
   move_pad(:,:,i) = circshift(move_pad(:,:,i),[-block_idx(i),-block_idx(i)]);  
end


%% Get spacing of activity (i.e. clumped or evenly distributed)
% s = skewness, ff = fano factor, n = number of ca events

s = cell(8,1);
ff = cell(8,1);
n = cell(8,1);
for curr_animal = 1:8
    s{curr_animal} = cell(14,1);
    ff{curr_animal} = cell(14,1);
    n{curr_animal} = cell(14,1);
    for curr_session = 1:14;
        n_rois = size(data(curr_animal).im(curr_session).roi_trace_thresh,1);
        s{curr_animal}{curr_session} = nan(n_rois,1);
        ff{curr_animal}{curr_session} = nan(n_rois,1);
        n{curr_animal}{curr_session} = nan(n_rois,1);
        for curr_cell = 1:n_rois
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session). ...
                roi_trace_thresh(curr_cell,:)>0]) == 1);
            
            curr_iei = diff(curr_act_starts);
            
            s{curr_animal}{curr_session}(curr_cell) = skewness(curr_iei);
            ff{curr_animal}{curr_session}(curr_cell) = std(curr_iei)./mean(curr_iei);
            n{curr_animal}{curr_session}(curr_cell) = length(curr_iei);
            
        end
    end
end


%% Plot avg movement activity of all ROIs combined

avg_moveroi_act = cell(8,14);
for curr_animal = 1:length(data);
    
    avg_moveroi_act(curr_animal,:) = cellfun(@(x,y) permute(nanmean(x(:,:,y)),[3 2 1]), ...
        {analysis(curr_animal).im(:).move_onset_aligned}, ...
        mat2cell(classified_rois(curr_animal).movement, ...
        size(classified_rois(curr_animal).movement,1),ones(14,1)),'uni',false);
    
end

avg_moveroi_act_cat = cell(1,14);
for i = 1:14
    curr_act = vertcat(avg_moveroi_act{:,i});
    [~,max_idx] = max(curr_act(:,32:end),[],2);
    [~,sort_idx] = sort(max_idx);
    
    avg_moveroi_act_cat{i} = curr_act(sort_idx,:);
end

move_start = find(analysis(1).surrounding_frames > 0,1);

figure;
for i = 1:14
   subplot(2,7,i);
   imagesc(avg_moveroi_act_cat{i});
   colormap(hot)   
   axis off
   line([move_start move_start],ylim,'color','w')
   title(['Day ' num2str(i)]);
end



% Plot average activity of movement cells on one day across days
use_session = 10;
avg_moveroi_act = cell(8,14);
for curr_animal = 1:length(data);    
    avg_moveroi_act(curr_animal,:) = cellfun(@(x) permute(nanmean( ...
        x(:,:,classified_rois(curr_animal).movement(:,use_session))),[3 2 1]), ...
        {analysis(curr_animal).im(:).move_onset_aligned},'uni',false);    
end

avg_moveroi_act_cat = cell(1,14);
curr_act = vertcat(avg_moveroi_act{:,use_session});
[~,max_idx] = max(curr_act(:,32:end),[],2);
[~,sort_idx] = sort(max_idx);
for i = 1:14
    curr_act = vertcat(avg_moveroi_act{:,i});
    avg_moveroi_act_cat{i} = curr_act(sort_idx,:);
end


%% Get classified ROI turnover

% Get fraction of currently classified ROIs that were classified previously
m2m = nan(length(data),14);
for curr_animal = 1:length(data)
m2m(curr_animal,:) = ...
    (sum((cumsum(classified_rois(curr_animal).movement,2) > 1).* ...
    classified_rois(curr_animal).movement))./ ...
    sum(classified_rois(curr_animal).movement);
end
q2q = nan(length(data),14);
for curr_animal = 1:length(data)
q2q(curr_animal,:) = ...
    (sum((cumsum(classified_rois(curr_animal).quiescent,2) > 1).* ...
    classified_rois(curr_animal).quiescent))./ ...
    sum(classified_rois(curr_animal).quiescent);
end
m2q = nan(length(data),14);
for curr_animal = 1:length(data)
m2q(curr_animal,:) = ...
    (sum((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
    classified_rois(curr_animal).quiescent))./ ...
    sum(classified_rois(curr_animal).quiescent);
end
q2m = nan(length(data),14);
for curr_animal = 1:length(data)
q2m(curr_animal,:) = ...
    (sum((cumsum(classified_rois(curr_animal).quiescent,2) >= 1).* ...
    classified_rois(curr_animal).movement))./ ...
    sum(classified_rois(curr_animal).movement);
end

n_rep = 5000;

m2m_shuff_1 = nan(length(data),14,n_rep);
q2q_shuff_1 = nan(length(data),14,n_rep);
m2q_shuff_1 = nan(length(data),14,n_rep);
q2m_shuff_1 = nan(length(data),14,n_rep);

m2m_shuff_2 = nan(length(data),14,n_rep);
q2q_shuff_2 = nan(length(data),14,n_rep);
m2q_shuff_2 = nan(length(data),14,n_rep);
q2m_shuff_2 = nan(length(data),14,n_rep);

% Get shuffled values of fraction previously classified
for curr_rep = 1:n_rep
    shuff_classified_movement = arrayfun(@(x) shake(classified_rois(x).movement,1),1:length(data),'uni',false);
    shuff_classified_quiescent = arrayfun(@(x) shake(classified_rois(x).quiescent,1),1:length(data),'uni',false);
    
    % Get fraction of currently classified ROIs that were classified previously
    for curr_animal = 1:length(data)
        m2m_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_movement{curr_animal},2) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2q_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_quiescent{curr_animal},2) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        m2q_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_movement{curr_animal},2) >= 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2m_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_quiescent{curr_animal},2) >= 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    
    
    shuff_classified_movement = arrayfun(@(x) shake(classified_rois(x).movement,2),1:length(data),'uni',false);
    shuff_classified_quiescent = arrayfun(@(x) shake(classified_rois(x).quiescent,2),1:length(data),'uni',false);
    
    % Get fraction of currently classified ROIs that were classified previously
    for curr_animal = 1:length(data)
        m2m_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_movement{curr_animal},2) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2q_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_quiescent{curr_animal},2) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        m2q_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_movement{curr_animal},2) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2m_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((cumsum(shuff_classified_quiescent{curr_animal},2) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    disp(curr_rep)
end

% Plot fraction of currently classified ROIs previously classified
figure; 
subplot(2,2,1);hold on;
errorbar(nanmean(m2m), ...
    nanstd(m2m)./sqrt(sum(~isnan(m2m))),'k','linewidth',2)
plot(permute(prctile(nanmean(m2m_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(m2m_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')

xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Movement to movement');
subplot(2,2,2);hold on;
errorbar(nanmean(q2q), ...
    nanstd(q2q)./sqrt(sum(~isnan(q2q))),'k','linewidth',2)
plot(permute(prctile(nanmean(q2q_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(q2q_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Quiescent to quiescent');
subplot(2,2,3);hold on;
errorbar(nanmean(m2q), ...
    nanstd(m2q)./sqrt(sum(~isnan(m2q))),'k','linewidth',2)
plot(permute(prctile(nanmean(m2q_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(m2q_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Movement to quiescent');
subplot(2,2,4);hold on;
errorbar(nanmean(q2m), ...
    nanstd(q2m)./sqrt(sum(~isnan(q2m))),'k','linewidth',2)
plot(permute(prctile(nanmean(q2m_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(q2m_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Quiescent to movement');



% Get fraction of currently move ROIs which will not be next day
declassify_move_frac = nan(length(data),13);
for curr_animal = 1:length(data)
    declassify_move_frac(curr_animal,:) = ...
        sum(classified_rois(curr_animal).movement(:,1:end-1) & ...
        ~classified_rois(curr_animal).movement(:,2:end))./ ...
        sum(classified_rois(curr_animal).movement(:,1:end-1));
end
figure;
errorbar(nanmean(declassify_move_frac), ...
    nanstd(declassify_move_frac)./sqrt(sum(~isnan(declassify_move_frac))),'k','linewidth',2)
xlabel('Day');
ylabel('Fraction of movement ROIs not classified next day');



%% Get correlation with activity and movement change over days

% Correlate thresholded activity with binary lever movement
activity_move_corr = cell(length(data),1);
for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    activity_move_corr{curr_animal} = nan(n_rois,n_days);
    for curr_day = 1:length(data(curr_animal).im)
        warning off
        activity_move_corr{curr_animal}(:,curr_day) = ...
            1-pdist2(data(curr_animal).im(curr_day).roi_trace_thresh, ...
            analysis(curr_animal).lever(curr_day).lever_move_frames','correlation');
        warning on
    end
    disp(curr_animal);
end

early_days = 1:3;
late_days = 11:14;

% Assign NaN to non-modulated cells on selected days
actitivty_move_corr_sig = activity_move_corr;
for curr_animal = 1:length(data)
    
    curr_movecells = any(classified_rois(curr_animal). ...
        movement(:,[early_days late_days]),2);
    curr_quiescentcells = any(classified_rois(curr_animal). ...
        quiescent(:,[early_days late_days]),2);
    
    actitivty_move_corr_sig{curr_animal}(~(curr_movecells | ...
        curr_quiescentcells),:) = NaN;
end

% Plot heatmap of changes in correlation
activity_move_corr_cat = vertcat(actitivty_move_corr_sig{:});
activity_move_corr_cat_earlylate = [ ...
    nanmean(activity_move_corr_cat(:,late_days),2), ...
    nanmean(activity_move_corr_cat(:,early_days),2)];

bin_edges = linspace(-1,1,500);
bin_centers = diff(bin_edges)/2+bin_edges(1:end-1);

movecorr_density = hist3(activity_move_corr_cat_earlylate,repmat({bin_centers},2,1));
h = fspecial('gaussian',10,2);
movecorr_density_smooth = imfilter(movecorr_density,h,'same');
figure;imagesc(movecorr_density_smooth);colormap(hot);

set(gca,'YDir','Normal');

zero_corr_bin = find(bin_centers > 0,1);
line(repmat(zero_corr_bin-0.5,2,1),ylim,'color','w')
line(xlim,repmat(zero_corr_bin-0.5,2,1),'color','w')

% Zoom in
filled_box = movecorr_density_smooth > ...
    prctile(movecorr_density_smooth(movecorr_density_smooth ~= 0),60);
min_lim = min([find(any(filled_box,1),1) find(fliplr(any(filled_box,1)),1) ...
    find(any(filled_box,2),1) find(flipud(any(filled_box,1)),1)]);
xlim([min_lim-0.5 length(filled_box)-min_lim+0.5])
ylim([min_lim-0.5 length(filled_box)-min_lim+0.5])

set(gca,'XTickLabel',round(bin_centers(get(gca,'XTick'))*100)/100);
set(gca,'YTickLabel',round(bin_centers(get(gca,'YTick'))*100)/100);

ylabel('Late (11:14) activity-movement correlation');
xlabel('Early (1:3) activity-movement correlation');


% Get fraction of ROIs in each quadrant per animal
m2m = nan(length(data),1);
q2q = nan(length(data),1);
m2q = nan(length(data),1);
q2m = nan(length(data),1);
for curr_animal = 1:length(data)    
    curr_corr = [ ...
        nanmean(actitivty_move_corr_sig{curr_animal}(:,early_days),2), ...
        nanmean(actitivty_move_corr_sig{curr_animal}(:,late_days),2)];
    curr_corr(any(isnan(curr_corr),2),:) = [];
    
    m2m(curr_animal) = sum(curr_corr(:,1) > 0 & curr_corr(:,2) > 0)./size(curr_corr,1);
    q2q(curr_animal) = sum(curr_corr(:,1) < 0 & curr_corr(:,2) < 0)./size(curr_corr,1);
    m2q(curr_animal) = sum(curr_corr(:,1) > 0 & curr_corr(:,2) < 0)./size(curr_corr,1);
    q2m(curr_animal) = sum(curr_corr(:,1) < 0 & curr_corr(:,2) > 0)./size(curr_corr,1);
       
    movecorr_density_animal(:,:,curr_animal) = ...
        hist3(curr_corr,repmat({bin_centers},2,1));   
end
figure; hold on;
corr_plot = [m2m q2q m2q q2m];
bar(nanmean(corr_plot),'FaceColor','None','linewidth',2);
errorbar(nanmean(corr_plot), ...
    nanstd(corr_plot)./sqrt(sum(~isnan(corr_plot))),'.k','linewidth',2);
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'+/+' '-/-' '+/-' '-/+'});
ylabel('Fraction of classified ROIs');


%% Timing of activity relative to movement onset

% Method 1: get mean of activity start timing for CR movements

act_starts_mean_all = cell(length(data),1);
for curr_animal = 1:length(data);
    
    use_frames = analysis(curr_animal).surrounding_frames >= 0 & analysis(curr_animal).surrounding_frames <= 90;

    act_starts = cell(14,size(data(curr_animal).im(1).roi_trace_df,1));
    
    for i = 1:14
        [curr_max,curr_max_idx] = max(analysis(curr_animal).im(i).move_onset_aligned(:,use_frames,:) > 0,[],2);
        curr_act_starts = arrayfun(@(x) curr_max_idx(curr_max(:,:,x),:,x),1:size(curr_max,3),'uni',false);
        act_starts(i,:) = curr_act_starts';
    end
    
    act_starts_mean = cellfun(@nanmedian,act_starts);
    
    % Get activity of only movement_related cells
    act_starts_mean(~classified_rois(curr_animal).movement') = NaN;    
    act_starts_mean_all{curr_animal} = act_starts_mean;
end

early_timing = cellfun(@(x) nanmean(x(1:3,:)),act_starts_mean_all,'uni',false);
late_timing = cellfun(@(x) nanmean(x(11:14,:)),act_starts_mean_all,'uni',false);

timing_bins = linspace(0,90,20);
timing_bins_plot = (timing_bins(1:end-1) + diff(timing_bins)/2)/29;
timing_bins(1) = -inf;
timing_bins(end) = inf;

early_timing_bins = cell2mat(cellfun(@(x) cumsum(histc(x,timing_bins))/ ...
    sum(~isnan(x)),early_timing,'uni',false));
early_timing_bins = early_timing_bins(:,1:end-1);

late_timing_bins = cell2mat(cellfun(@(x) cumsum(histc(x,timing_bins))/ ...
    sum(~isnan(x)),late_timing,'uni',false));
late_timing_bins = late_timing_bins(:,1:end-1);

figure; hold on;
errorbar(timing_bins_plot,nanmean(early_timing_bins),nanstd(early_timing_bins)./ ...
    sqrt(sum(~isnan(early_timing_bins))),'k','linewidth',2);
errorbar(timing_bins_plot,nanmean(late_timing_bins),nanstd(late_timing_bins)./ ...
    sqrt(sum(~isnan(late_timing_bins))),'r','linewidth',2);

legend({'Early (days 1:3)' 'Late (days 11:14)'},'location','se');
ylabel('Cumulative fraction');
xlabel('Time (s)')
title('Distribution of movement-cell mean activity onset')


% Method 2: use maximum of activity-aligned movement onset

back_time = 29*0;
forward_time = 29*3;

max_move = cell(length(data),1);

for curr_animal = 1:length(data)
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    curr_max_move = nan(14,n_rois);
    
    for curr_session = 1:14;
        
        use_rois = find(classified_rois(curr_animal).movement(:,curr_session));
        
        for curr_cell = use_rois';
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
            curr_act_starts(curr_act_starts-back_time <= 0 | ...
                curr_act_starts + forward_time > length(analysis(curr_animal).lever(curr_session).lever_move_frames)) = [];
            
            curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_starts(x)-back_time: ...
                curr_act_starts(x)+forward_time,(1:length(curr_act_starts))','uni',false));
            
            % Movement-binarzed lever
            curr_act_move = analysis(curr_animal).lever(curr_session).lever_move_frames(curr_act_move_idx);
            
            [~,curr_max_move_frame] = max(nanmean(curr_act_move));
            
            curr_max_move(curr_session,curr_cell) = curr_max_move_frame-(back_time+1);
        end
    end
    
    max_move{curr_animal} = curr_max_move;
    disp(curr_animal);
    
end

early_timing = cellfun(@(x) nanmean(x(1:3,:)),max_move,'uni',false);
late_timing = cellfun(@(x) nanmean(x(11:14,:)),max_move,'uni',false);

timing_bins = linspace(-back_time,forward_time,30);
timing_bins_plot = (timing_bins(1:end-1) + diff(timing_bins)/2)/29;
timing_bins(1) = -inf;
timing_bins(end) = inf;

early_timing_bins = cell2mat(cellfun(@(x) cumsum(histc(x,timing_bins))/ ...
    sum(~isnan(x)),early_timing,'uni',false));
early_timing_bins = early_timing_bins(:,1:end-1);

late_timing_bins = cell2mat(cellfun(@(x) cumsum(histc(x,timing_bins))/ ...
    sum(~isnan(x)),late_timing,'uni',false));
late_timing_bins = late_timing_bins(:,1:end-1);

figure; hold on;
errorbar(timing_bins_plot,nanmean(early_timing_bins),nanstd(early_timing_bins)./ ...
    sqrt(sum(~isnan(early_timing_bins))),'k','linewidth',2);
errorbar(timing_bins_plot,nanmean(late_timing_bins),nanstd(late_timing_bins)./ ...
    sqrt(sum(~isnan(late_timing_bins))),'r','linewidth',2);

legend({'Early (days 1:3)' 'Late (days 11:14)'},'location','nw');
ylabel('Cumulative fraction');
xlabel('Time (s)')
title('Distribution of movement time relative to activity onset')


%% PCA/jPCA 
% (this was kind of a bust)

% Get average activity
avg_moveroi_act = cell(8,14);
avg_act = cell(8,14);

for curr_animal = 1:length(data);
    
    avg_moveroi_act(curr_animal,:) = cellfun(@(x,y) permute(nanmean(x(:,:,y)),[3 2 1]), ...
        {analysis(curr_animal).im(:).move_onset_aligned}, ...
        mat2cell(classified_rois(curr_animal).movement, ...
        size(classified_rois(curr_animal).movement,1),ones(14,1)),'uni',false);
    
    avg_act(curr_animal,:) = cellfun(@(x,y) permute(nanmean(x),[3 2 1]), ...
        {analysis(curr_animal).im(:).move_onset_aligned},'uni',false);
    
end

% jPCA on average aligned activity
figure;
for curr_animal = 1:size(avg_act,1)
    for curr_session = 1:size(avg_act,2)
        use_data = avg_act{curr_animal,curr_session};
        curr_plot = subplot(size(avg_act,1),size(avg_act,2), ...
            curr_session + (curr_animal-1)*size(avg_act,2));
        AP_jPCA(use_data',size(use_data,2),'rotation',curr_plot);
        drawnow
    end
end

figure;
for curr_animal = 1:size(avg_act,1)
    curr_plot = subplot(3,3,curr_animal);
    use_data = horzcat(avg_act{curr_animal,:});
    AP_jPCA(use_data',repmat(size(avg_act{curr_animal,1},2),1,14),'rotation',curr_plot);
    drawnow
end



% PCA on average activity over days
figure;
for curr_animal = 1:length(data);
    use_data = horzcat(avg_act{curr_animal,:});
    use_data = zscore(use_data,[],2);
        
    [coeff,score,latent] = princomp(use_data');
    a = mat2cell(score,repmat(size(avg_act{curr_animal,1},2), ...
        size(score,1)/size(avg_act{curr_animal,1},2),1),size(score,2));
    
    subplot(3,3,curr_animal); hold on
    col = jet(14);
    for i = 1:14
        plot3(a{i}(:,1),a{i}(:,2),a{i}(:,3),'color',col(i,:),'linewidth',2)
    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    xlabel('PC1');ylabel('PC2');zlabel('PC3');
    view(3);
    camproj('perspective');
    axis tight;
    
end

% Sanity check:
% PCA on average activity over days (shuffle cells independently by day)
figure;
for curr_animal = 1:length(data);
    
    cell_act = cellfun(@(x) mat2cell(x,ones(size(x,1),1), ...
        size(x,2)),avg_act(curr_animal,:),'uni',false);
    use_data = cell2mat(shake(horzcat(cell_act{:}),2));
    use_data = zscore(use_data,[],2);
    
    [coeff,score,latent] = princomp(use_data');
    a = mat2cell(score,repmat(size(avg_act{curr_animal,1},2), ...
        size(score,1)/size(avg_act{curr_animal,1},2),1),size(score,2));
    
    subplot(3,3,curr_animal); hold on
    col = jet(14);
    for i = 1:14
        plot3(a{i}(:,1),a{i}(:,2),a{i}(:,3),'color',col(i,:),'linewidth',2)
    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    xlabel('PC1');ylabel('PC2');zlabel('PC3');
    view(3);
    camproj('perspective');
    axis tight;
    
end


% Get average of PC1-PC2 activity (this was to show that the major axis of
% change is quiescence cells and not movement cells, but it didn't really
% work)
figure;
for curr_animal = 1:length(data);
    use_data = horzcat(avg_act{curr_animal,:});
    use_data = zscore(use_data,[],2);
        
    [coeff,score,latent] = princomp(use_data');
    a = mat2cell(score,repmat(size(avg_act{curr_animal,1},2), ...
        size(score,1)/size(avg_act{curr_animal,1},2),1),size(score,2));
    
    b = nanmean(cat(3,a{:}),3);
    
    [~,sort_idx] = sort(coeff(:,1) - coeff(:,2));
    
    % Split sorted cells into even number of groups
    n_grps = 10;
    grp_n = diff(ceil(prctile(sort_idx,linspace(0,100,n_grps))));
    grp_n(end) = grp_n(end) + (length(sort_idx)-sum(grp_n));
    sort_idx_grp = mat2cell(sort_idx,grp_n,1);
    
    c = cellfun(@(x) nanmean(reshape(nanmean(use_data(x,:))',[],14),2),sort_idx_grp,'uni',false);
    
    subplot(3,3,curr_animal); hold on
    col = jet(n_grps-1);
    for i = 1:n_grps-1
        plot(c{i},'color',col(i,:),'linewidth',2)
    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);

end


% Sort by PC1-PC2, look at classification
figure;
r = nan(length(data),5);
p = nan(length(data),5);
for curr_animal = 1:length(data);
    
    use_data = horzcat(avg_act{curr_animal,:});
    use_data = zscore(use_data,[],2);
        
    [coeff,score,latent] = princomp(use_data');
    
    [~,sort_idx] = sort(coeff(:,1) - coeff(:,2));
       
    m = nanmean(classified_rois(curr_animal).movement,2);
    q = nanmean(classified_rois(curr_animal).quiescent,2);
    
    [rc,pc] = corrcoef([m-q coeff(:,1:5)]);
    
    r(curr_animal,:) = rc(1,2:end);
    p(curr_animal,:) = pc(1,2:end);
    
    subplot(3,3,curr_animal);
    plot(m-q,coeff(:,2),'.k')
    
    %%%%% Look at which coefficients correlatie with movement and what
    %%%%% direction those coefficients go in (I think it just ends up
    %%%%% saying that move cells go down and quiescent cells go up
    % 
    % also I don't think this makes sense: m/q cells will share temporal
    % activity and so of course will have similar coefficients, that would
    % only make sense if activity was averaged in time <-- do that
 
end

a = mat2cell(score,repmat(size(avg_act{curr_animal,1},2), ...
    size(score,1)/size(avg_act{curr_animal,1},2),1),size(score,2));
b = cell2mat(cellfun(@nanmean,a,'uni',false));


% PCA on average movement/quiescent activity
avg_act_movement = cell(length(data),1);
avg_act_quiescent = cell(length(data),1);
for curr_animal = 1:length(data);
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    avg_act_movement{curr_animal} = nan(n_rois,14);
    avg_act_quiescent{curr_animal} = nan(n_rois,14);
    
    for i = 1:14
        avg_act_movement{curr_animal}(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,analysis(curr_animal).lever(i).lever_move_frames),2);
        avg_act_quiescent{curr_animal}(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,~analysis(curr_animal).lever(i).lever_move_frames),2);
    end

end





%% Get average movement/quiescent activity for m2q cells
% Do they keep their spontaneous/quiescent activity but lower movement act?

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

% Plot average activity of switch cells over days
avg_movement_act_m2q_cat = cell2mat(cellfun(@nanmean,avg_movement_act_m2q,'uni',false));
avg_quiescent_act_m2q_cat = cell2mat(cellfun(@nanmean,avg_quiescent_act_m2q,'uni',false));
avg_total_act_m2q_cat = cell2mat(cellfun(@nanmean,avg_total_act_m2q,'uni',false));

figure; hold on
errorbar(nanmean(avg_movement_act_m2q_cat),nanstd(avg_movement_act_m2q_cat)./ ...
    sqrt(sum(~isnan(avg_movement_act_m2q_cat))),'k','linewidth',2);
errorbar(nanmean(avg_quiescent_act_m2q_cat),nanstd(avg_quiescent_act_m2q_cat)./ ...
    sqrt(sum(~isnan(avg_quiescent_act_m2q_cat))),'r','linewidth',2);

ylabel('Average activity')
xlabel('Day');
legend({'Movement' 'Quiescence'});

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


% Get activty-triggered movement for m2q cells on m/q days
back_time = 29*3;
forward_time = 29*3;
class_movesurround_act = nan(4,back_time+forward_time+1,length(data));
for curr_animal = 1:length(data)
    
    curr_m2q = ...
        (cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent;
    
    use_m2q = any(curr_m2q,2);    
    use_m2m = any(sum(classified_rois(curr_animal).movement,2) > 1) & ...
        ~any(classified_rois(curr_animal).quiescent,2);
    use_q2q = any(sum(classified_rois(curr_animal).quiescent,2) > 1) & ...
        ~any(classified_rois(curr_animal).movement,2);
    
    curr_act_move_all = nan(size(curr_m2q,1),back_time+forward_time+1,14);
    for curr_session = 1:14;
        for curr_cell = find(use_m2m | use_m2q | use_q2q)';
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
            curr_act_starts(curr_act_starts-back_time <= 0 | ...
                curr_act_starts + forward_time > length(analysis(curr_animal).lever(curr_session).lever_move_frames)) = [];
            
            curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_starts(x)-back_time: ...
                curr_act_starts(x)+forward_time,(1:length(curr_act_starts))','uni',false));
            
            % Movement-binarzed lever
            curr_act_move = analysis(curr_animal).lever(curr_session).lever_move_frames(curr_act_move_idx);
            curr_act_move_all(curr_cell,:,curr_session) = ...
                nanmean(curr_act_move);
        end
    end
    
    curr_m2q_act_mclass = cell2mat(arrayfun(@(x) nanmean(curr_act_move_all(x,:, ...
        classified_rois(curr_animal).movement(x,:)),3),find(use_m2q),'uni',false));
    
    curr_m2q_act_qclass = cell2mat(arrayfun(@(x) nanmean(curr_act_move_all(x,:, ...
        classified_rois(curr_animal).quiescent(x,:)),3),find(use_m2q),'uni',false));
    
    curr_m2m_act_mclass = cell2mat(arrayfun(@(x) nanmean(curr_act_move_all(x,:, ...
        classified_rois(curr_animal).movement(x,:)),3),find(use_m2m),'uni',false));
    
    curr_q2q_act_qclass = cell2mat(arrayfun(@(x) nanmean(curr_act_move_all(x,:, ...
        classified_rois(curr_animal).quiescent(x,:)),3),find(use_q2q),'uni',false));
    
    class_movesurround_act(:,:,curr_animal) = [ ...
        nanmean(curr_m2q_act_mclass); ...
        nanmean(curr_m2q_act_qclass); ...
        nanmean(curr_m2m_act_mclass); ...
        nanmean(curr_q2q_act_qclass)];
    
    disp(curr_animal);
    
end

% Plot mean of m2q/m, m2q/q, and m2m/m, and m2q/m-m2m/m
figure;

xvals = ((1:size(class_movesurround_act,2))-(back_time+1))/29;

subplot(2,1,1); hold on;
m2q_q_mean = nanmean(class_movesurround_act(1,:,:),3);
m2q_q_sem =  nanstd(class_movesurround_act(1,:,:),[],3)./sqrt(sum(~isnan( ...
    class_movesurround_act(1,:,:)),3));

m2q_m_mean = nanmean(class_movesurround_act(2,:,:),3);
m2q_m_sem =  nanstd(class_movesurround_act(2,:,:),[],3)./sqrt(sum(~isnan( ...
    class_movesurround_act(1,:,:)),3));

m2m_m_mean = nanmean(class_movesurround_act(3,:,:),3);
m2m_m_sem =  nanstd(class_movesurround_act(3,:,:),[],3)./sqrt(sum(~isnan( ...
    class_movesurround_act(1,:,:)),3));

q2q_q_mean = nanmean(class_movesurround_act(4,:,:),3);
q2q_q_sem =  nanstd(class_movesurround_act(4,:,:),[],3)./sqrt(sum(~isnan( ...
    class_movesurround_act(1,:,:)),3));

plot(xvals,m2q_q_mean,'k','linewidth',2);
plot(xvals,m2q_m_mean,'r','linewidth',2);
plot(xvals,m2m_m_mean,'b','linewidth',2);
plot(xvals,q2q_q_mean,'m','linewidth',2);

jbfill(xvals, ...
    m2q_q_mean+m2q_q_sem,m2q_q_mean-m2q_q_sem,'k'); hold on;
jbfill(xvals, ...
    m2q_m_mean+m2q_m_sem,m2q_m_mean-m2q_m_sem,'r'); hold on;
jbfill(xvals, ...
    m2m_m_mean+m2m_m_sem,m2m_m_mean-m2m_m_sem,'b'); hold on;
jbfill(xvals, ...
    q2q_q_mean+q2q_q_sem,q2q_q_mean-q2q_q_sem,'m'); hold on;
line([0 0],ylim,'color','r','linestyle','--');
xlabel('Time from activity onset');
ylabel('Movement probability')
legend({'M2Q/M' 'M2Q/Q' 'M2M/M' 'Q2Q/Q'});

subplot(2,1,2);
curr_sub = class_movesurround_act(1,:,:) - class_movesurround_act(3,:,:);
curr_mean = nanmean(curr_sub,3);
curr_sem =  nanstd(curr_sub,[],3)./sqrt(sum(~isnan(curr_sub),3));
plot(xvals,curr_mean,'k','linewidth',2);
jbfill(xvals, ...
    curr_mean+curr_sem,curr_mean-curr_sem,'k'); hold on;
line([xvals(1) xvals(end)],[0 0],'color','r','linestyle','--');
line([0 0],ylim,'color','r','linestyle','--');
xlabel('Time from activity onset');
ylabel('\Delta Movement probability');
title('M2Q/M - M2M/M');





%% Correlation of changes across sessions
% this is also in 150327

pop_corr = nan(14,14,length(analysis));
pop_diff_corr = nan(13,13,length(analysis));
pop_diff_corr_allcomp = nan(13,13,length(analysis));

for curr_animal = 1:length(analysis);
    
    curr_movecells = any(classified_rois(curr_animal).movement,2);
    curr_quiescentcells = any(classified_rois(curr_animal).quiescent,2);
    
    use_frames = analysis(curr_animal).surrounding_frames >= 0 &  ...
        analysis(curr_animal).surrounding_frames <= 90;
    
    use_cells = 1:length(curr_movecells);
     
    
    use_sessions = find(cellfun(@(x) ~isempty(x), ...
        {analysis(curr_animal).im(:).move_onset_aligned}));
              
    % Temporal activity
    temp_act = cellfun(@(act) reshape(permute(nanmean(act(:,use_frames,use_cells)),[2 3 1]),[],1), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
%     % Temporal activity, shuffled independently by cell
%     curr_act = cellfun(@(x) permute(nanmean(x,1),[3,2,1]), ...
%         {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);    
%     curr_act_cellseparate = permute(mat2cell(permute(cat(3,curr_act{:}),[2 1 3]), ...
%         size(curr_act{1},2),ones(size(curr_act{1},1),1),ones(length(curr_act),1)),[2 3 1]);
%     curr_act_cellseparate_shuff = shake(curr_act_cellseparate,2);
%     temp_act_cat = cell2mat(curr_act_cellseparate_shuff);
    
    % Mean binary activity (% active trials)
    %temp_act = cellfun(@(act) permute(nanmean(any(act,2),1),[3 2 1]), ...
    %    {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Average activity during movement
%     temp_act = cell(14,1);
%     for i = 1:14
%         temp_act{i} = nanmean(data(curr_animal).im(i). ...
%             roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames),2);
%     end

    temp_act_cat = horzcat(temp_act{:});
%     r = shake(diff(temp_act_cat,[],2),2);
%     m = cumsum(r,2);
%     temp_act_cat = [temp_act_cat(:,1) repmat(temp_act_cat(:,1),1,13) + ...
%         m];
    
    pop_corr(use_sessions,use_sessions,curr_animal) = corrcoef(temp_act_cat);
    
    temp_act_diff_all = cell(size(temp_act_cat,2));
    for pre = 1:size(temp_act_cat,2)
        for post = 1:size(temp_act_cat,2)
            temp_act_diff_all{pre,post} = temp_act_cat(:,post) - ...
                temp_act_cat(:,pre);
        end
    end
        
    diff_inout = cell(size(temp_act_cat,2),1);
    for curr_session = 1:size(temp_act_cat,2)   
        diff_in = horzcat(temp_act_diff_all{:,curr_session});
        diff_out = horzcat(temp_act_diff_all{curr_session,:});
        diff_inout{curr_session} = [diff_in diff_out];
    end
    diff_inout_allcorr = mat2cell(corrcoef(horzcat(diff_inout{:}),'rows','complete'), ...
        repmat(size(temp_act_cat,2)*2,size(temp_act_cat,2),1), ...
        repmat(size(temp_act_cat,2)*2,size(temp_act_cat,2),1));
    %diff_inout_corr_mean = cellfun(@(x) nanmean(AP_itril(x,-size(temp_act_cat,2)+1)), ...
    %    diff_inout_allcorr);
    diff_inout_corr_mean = cellfun(@(x) nanmean(reshape(x(size(temp_act_cat,2)+1:end, ...
        1:size(temp_act_cat,2)),[],1)), ...
        diff_inout_allcorr);
    
    pop_diff_corr_allcomp(use_sessions,use_sessions,curr_animal) = diff_inout_corr_mean;
    
    
    temp_act_diff = diff(temp_act_cat,[],2);
    act_diff_corr = corrcoef(temp_act_diff);
     
    use_sessions_diff = use_sessions(1:end-1);
    pop_diff_corr(use_sessions_diff,use_sessions_diff,curr_animal) = act_diff_corr;
end

figure;

subplot(1,3,1);
imagesc(nanmean(pop_corr,3));
colormap(gray);
caxis([0 1]);
xlabel('Session')
ylabel('Session')
title('Pop corr')

subplot(1,3,2);imagesc(nanmean(pop_diff_corr,3));
colormap(gray);
caxis([-1 1]);
xlabel('Session')
ylabel('Session')
title('Pop diff corr')

subplot(1,3,3);imagesc(nanmean(pop_diff_corr_allcomp,3));
colormap(gray);
caxis([-1 1]);
ylabel('Session')
xlabel('Session')
title('Pop diff corr allcomp');

% Pop_diff_corr diagonals
pop_diff_corr_allcomp_diag = cell(length(analysis),14);
for curr_animal = 1:length(analysis)
   pop_diff_corr_allcomp_diag(curr_animal,:) = ...
       arrayfun(@(x) diag(pop_diff_corr_allcomp(:,:,curr_animal),x),0:13,'uni',false);
end
pop_diff_corr_allcomp_diag_mean = cellfun(@nanmean,pop_diff_corr_allcomp_diag);

figure;
errorbar(0:13,nanmean(pop_diff_corr_allcomp_diag_mean), ...
    nanstd(pop_diff_corr_allcomp_diag_mean)./ ...
    sqrt(sum(~isnan(pop_diff_corr_allcomp_diag_mean))),'k','linewidth',2);
line(xlim,[0 0],'linewidth',2,'color','r','linestyle','--');
xlabel('Session difference');
ylabel('Mean difference correlation');


%% Predict lever from activity (using PCA modes TRIALS)

n_split = 5;

framerate = 28;

lever_predict_corr = cell(length(data),14);
for curr_animal = 1:length(data)
    
    real_active = cell(14,1);
    predict_active = cell(14,1);    
    
    % currently set up to include 1 second before movement   
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > (-framerate) & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    %use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
    %    analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    curr_lever = cell(n_sessions,1);
    curr_lever(use_sessions) = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    curr_act = cell(n_sessions,1);
    curr_act(use_sessions) = cellfun(@(activity) ...
        permute(activity(:,use_frames,use_cells),[2 3 1]), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);  
    
    for curr_session = 1:14;
        
        n_trials = size(curr_act{curr_session},3);
        
        use_cells = 1:size(classified_rois(curr_animal).movement,1);
        use_curr_act = reshape(curr_act{curr_session}(:,use_cells,:),[], ...
            n_trials);
                
        % Decompse lever/activity into principle components
        [lever_coeff,lever_score,lever_latent] = princomp(zscore(curr_lever{curr_session}));
        [act_coeff,act_score,act_latent] = princomp(use_curr_act);
                               
        % Get number of PCs for lever/activity that explain n% variance
        n_pcs_lever = find(cumsum(lever_latent)/sum(lever_latent) > 0.8,1);
        n_pcs_act = find(cumsum(act_latent)/sum(act_latent) > 0.8,1);
        
        % Reconstruct lever given number of PCs
        [~,lever_reconstructed] = pcares(zscore(curr_lever{curr_session}),n_pcs_lever);
       
        n_split_trials = diff(ceil(linspace(0,n_trials,n_split+1)));
        split_trials = mat2cell(randperm(n_trials),1,n_split_trials);
        
        % Predict lever through activity to lever PCs via cross validation
        predict_lever = nan(size(curr_lever{curr_session}));
        for curr_split = 1:n_split
            
            train_trials = horzcat(split_trials{setdiff(1:n_split,curr_split)});
            test_trials = split_trials{curr_split};
            
            train_lever_coeff = lever_coeff(train_trials,1:n_pcs_lever);
            train_act_coeff = act_coeff(train_trials,1:n_pcs_act);
            
            w = train_act_coeff\train_lever_coeff;
            
            predict_lever_coeff = act_coeff(test_trials,1:n_pcs_act)*w;
            predict_lever(:,test_trials) = ...
                lever_score(:,1:n_pcs_lever)*predict_lever_coeff';
            
        end
        
        % Get similarity between predicted and actual lever
        lever_predict_corr{curr_animal,curr_session} = nan(n_trials,1);
        for i = 1:n_trials
            %curr_corr = corrcoef(predict_lever(:,i), ...
            %    curr_lever{curr_session}(:,i));
            curr_corr = corrcoef(predict_lever(:,i), ...
                lever_reconstructed(:,i));
                       
            lever_predict_corr{curr_animal,curr_session}(i) = curr_corr(2);
        end    
    end   
    disp(curr_animal);
end

lever_predict_corr_median = cellfun(@nanmedian,lever_predict_corr);
figure;errorbar(nanmean(lever_predict_corr_median), ...
    nanstd(lever_predict_corr_median)./sqrt(sum( ...
    ~isnan(lever_predict_corr_median))),'k','linewidth',2)
ylabel('Median real/predicted lever correlation');
xlabel('Day')

%% Predict lever from activity (using PCA modes TRIALS - concat all lever)

n_split = 5;

framerate = 28;

lever_predict_corr = cell(length(data),14);
for curr_animal = 1:length(data)
    
    real_active = cell(14,1);
    predict_active = cell(14,1);    
    
    % currently set up to include 1 second before movement   
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > (-framerate) & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    %use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
    %    analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    curr_lever = cell(n_sessions,1);
    curr_lever(use_sessions) = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    curr_act = cell(n_sessions,1);
    curr_act(use_sessions) = cellfun(@(activity) ...
        permute(activity(:,use_frames,use_cells),[2 3 1]), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);  
    
    % Decompose lever
    [lever_coeff,lever_score,lever_latent] = princomp(zscore(horzcat(curr_lever{:})));
    n_pcs_lever = find(cumsum(lever_latent)/sum(lever_latent) > 0.8,1);
    trials_per_session = cellfun(@(x) size(x,2),curr_lever);
    lever_coeff_session = mat2cell(lever_coeff(:,1:n_pcs_lever), ...
        trials_per_session,n_pcs_lever);
    
    % Reconstruct lever given number of PCs
    [~,lever_reconstructed] = pcares(zscore(horzcat(curr_lever{:})),n_pcs_lever);
    lever_reconstructed_session = mat2cell(lever_reconstructed, ...
        size(lever_reconstructed,1),trials_per_session);

    for curr_session = 1:14;
        
        n_trials = size(curr_act{curr_session},3);
        
        use_cells = 1:size(classified_rois(curr_animal).movement,1);
        use_curr_act = reshape(curr_act{curr_session}(:,use_cells,:),[], ...
            n_trials);
                
        % Decompse lever/activity into principle components        
        [act_coeff,act_score,act_latent] = princomp(use_curr_act);
                       
        % Get number of PCs for lever/activity that explain n% variance
        n_pcs_act = find(cumsum(act_latent)/sum(act_latent) > 0.8,1);    
       
        n_split_trials = diff(ceil(linspace(0,n_trials,n_split+1)));
        split_trials = mat2cell(randperm(n_trials),1,n_split_trials);
        
        % Predict lever through activity to lever PCs via cross validation
        predict_lever = nan(size(curr_lever{curr_session}));
        for curr_split = 1:n_split
            
            train_trials = horzcat(split_trials{setdiff(1:n_split,curr_split)});
            test_trials = split_trials{curr_split};
            
            train_lever_coeff = lever_coeff_session{curr_session}(train_trials,1:n_pcs_lever);
            train_act_coeff = act_coeff(train_trials,1:n_pcs_act);
            
            w = train_lever_coeff\train_act_coeff;
            
            predict_lever_coeff = act_coeff(test_trials,1:n_pcs_act)/w;
            predict_lever(:,test_trials) = ...
                lever_score(:,1:n_pcs_lever)*predict_lever_coeff';
            
        end
        
        % Get similarity between predicted and actual lever
        lever_predict_corr{curr_animal,curr_session} = nan(n_trials,1);
        for i = 1:n_trials
            %curr_corr = corrcoef(predict_lever(:,i), ...
            %    curr_lever{curr_session}(:,i));
            curr_corr = corrcoef(predict_lever(:,i), ...
                lever_reconstructed_session{curr_session}(:,i));
                       
            lever_predict_corr{curr_animal,curr_session}(i) = curr_corr(2);
        end    
    end   
    disp(curr_animal);
end

lever_predict_corr_median = cellfun(@nanmedian,lever_predict_corr);
figure;errorbar(nanmean(lever_predict_corr_median), ...
    nanstd(lever_predict_corr_median)./sqrt(sum( ...
    ~isnan(lever_predict_corr_median))),'k','linewidth',2)
ylabel('Median real/predicted lever correlation');
xlabel('Day')

%% Predict lever from activity (using lever PCA, mean act SCORES)

n_split = 5;

framerate = 28;

lever_predict_corr = cell(length(data),14);
for curr_animal = 1:length(data)
    
    real_active = cell(14,1);
    predict_active = cell(14,1);    
    
    % currently set up to include 1 second before movement   
    movement_start_time = 1001;
    movement_use_time = 2000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > (-framerate) & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    %use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
    %    analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    curr_lever = cell(n_sessions,1);
    curr_lever(use_sessions) = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    curr_act = cell(n_sessions,1);
    curr_act(use_sessions) = cellfun(@(activity) ...
        permute(activity(:,use_frames,use_cells),[2 3 1]), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Decompose lever
    [lever_coeff,lever_score,lever_latent] = princomp(zscore(horzcat(curr_lever{:}))');
    n_pcs_lever = find(cumsum(lever_latent)/sum(lever_latent) > 0.8,1);
        
    trials_per_session = cellfun(@(x) size(x,2),curr_lever);
    lever_score_session = mat2cell(lever_score(:,1:n_pcs_lever), ...
        trials_per_session,n_pcs_lever);
    
    % Reconstruct lever given number of PCs
    [~,lever_reconstructed] = pcares(zscore(horzcat(curr_lever{:})),n_pcs_lever);
    lever_reconstructed_session = mat2cell(lever_reconstructed, ...
        size(lever_reconstructed,1),trials_per_session);

    for curr_session = 1:14;
        
        n_trials = size(curr_act{curr_session},3);
        
        use_cells = 1:size(classified_rois(curr_animal).movement,1);
        use_curr_act = reshape(curr_act{curr_session}(1:60,use_cells,:),[], ...
            n_trials);
                
        % Decompse lever/activity into principle components        
        %[act_coeff,act_score,act_latent] = princomp(use_curr_act');
        
        % Get number of PCs for lever/activity that explain n% variance
        %n_pcs_act = find(cumsum(act_latent)/sum(act_latent) > 0.8,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% TEST use all avg activity
        act_score = use_curr_act';
        n_pcs_act = size(act_score,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        n_split_trials = diff(ceil(linspace(0,n_trials,n_split+1)));
        split_trials = mat2cell(randperm(n_trials),1,n_split_trials);
        
        % Predict lever through activity to lever PCs via cross validation
        predict_lever = nan(size(curr_lever{curr_session}));
        for curr_split = 1:n_split
            
            train_trials = horzcat(split_trials{setdiff(1:n_split,curr_split)});
            test_trials = split_trials{curr_split};
                       
            train_lever_score = lever_score_session{curr_session}(train_trials,1:n_pcs_lever);
            train_act_score = act_score(train_trials,1:n_pcs_act);
            
            w = train_lever_score\train_act_score;
            
            predict_lever_score = act_score(test_trials,1:n_pcs_act)/w;
            predict_lever(:,test_trials) = ...
                transpose(predict_lever_score*lever_coeff(:,1:n_pcs_lever)');
            
        end
        
        % Get similarity between predicted and actual lever
        lever_predict_corr{curr_animal,curr_session} = nan(n_trials,1);
        for i = 1:n_trials
            %curr_corr = corrcoef(predict_lever(:,i), ...
            %    curr_lever{curr_session}(:,i));
            curr_corr = corrcoef(predict_lever(:,i), ...
                lever_reconstructed_session{curr_session}(:,i));
                       
            lever_predict_corr{curr_animal,curr_session}(i) = curr_corr(2);
        end    
    end   
    disp(curr_animal);
end

lever_predict_corr_median = cellfun(@nanmedian,lever_predict_corr);
figure;errorbar(nanmean(lever_predict_corr_median), ...
    nanstd(lever_predict_corr_median)./sqrt(sum( ...
    ~isnan(lever_predict_corr_median))),'k','linewidth',2)
ylabel('Median real/predicted lever correlation');
xlabel('Day')

%% Predict lever from activity (using lever CHANNELS, any act)

n_split = 5;

framerate = 28;

% currently set up to include 1 second before movement
movement_baseline_time = 500; % time to estimate baseline, 1:n
movement_start_time = 1001;
movement_use_time = 3000;
    
n_lever_channels = 30;
lever_flip_times = linspace(1,movement_use_time,n_lever_channels);
lever_flip_times = round(lever_flip_times(1:end-1));
lever_channels = zeros(movement_use_time,n_lever_channels);
for curr_channel = 1:n_lever_channels-1
    lever_channels(lever_flip_times(curr_channel):end,curr_channel) = 1; 
end
% Add an offset channel
lever_channels(n_lever_channels,:) = 1;

lever_predict_corr = cell(length(data),14);
for curr_animal = 1:length(data)
    
    real_active = cell(14,1);
    predict_active = cell(14,1);    
  
    % to correspond to movement
    %use_frames = ...
    %    analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    use_frames = analysis(curr_animal).surrounding_frames > framerate*(0.5/1000) & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    curr_lever = cell(n_sessions,1);
    curr_lever(use_sessions) = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time-1),movement(trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    curr_lever_baseline = cell(n_sessions,1);
    curr_lever_baseline(use_sessions) = cellfun(@(movement,trials) cellfun(@(x) ...
        nanmean(x(1:movement_baseline_time)),movement(trials))', ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    curr_act = cell(n_sessions,1);
    curr_act(use_sessions) = cellfun(@(activity) ...
        permute(activity(:,use_frames,use_cells),[2 3 1]), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Decompose lever
    lever_cat = bsxfun(@minus,horzcat(curr_lever{:}),horzcat(curr_lever_baseline{:}));
    lever_channel_score = lever_channels\lever_cat;
    trials_per_session = cellfun(@(x) size(x,2),curr_lever);
    lever_channel_score_session = mat2cell(lever_channel_score', ...
        trials_per_session,n_lever_channels);
    
    % Reconstruct lever given number of PCs
    lever_reconstructed = lever_channels*lever_channel_score;
    lever_reconstructed_session = mat2cell(lever_reconstructed, ...
        size(lever_reconstructed,1),trials_per_session);

    for curr_session = 1:14;
        
        n_trials = size(curr_act{curr_session},3);       
        use_cells = classified_rois(curr_animal).movement(:,curr_session);
        
        act_score = reshape(max(curr_act{curr_session}(:,use_cells,:)),[], ...
             n_trials)';
        % To use just activity onsets
        %act_onsets = reshape( ...
        %    [zeros(1,length(use_cells),n_trials);...
        %    diff(curr_act{curr_session}(:,use_cells,:) > 0)],[],n_trials);
        %act_score(act_onsets' ~= 1) = 0;
               
        n_split_trials = diff(ceil(linspace(0,n_trials,n_split+1)));
        split_trials = mat2cell(randperm(n_trials),1,n_split_trials);
        
        % Predict lever through activity to lever PCs via cross validation
        predict_lever = nan(size(curr_lever{curr_session}));
        for curr_split = 1:n_split
            
            train_trials = horzcat(split_trials{setdiff(1:n_split,curr_split)});
            test_trials = split_trials{curr_split};
                       
            train_lever_score = lever_channel_score_session{curr_session}(train_trials,:);
            train_act_score = act_score(train_trials,:);
            
            w = train_lever_score\train_act_score;
            
            predict_lever_score = act_score(test_trials,:)/w;
            predict_lever(:,test_trials) = ...
                lever_channels*predict_lever_score';
            
        end
        
        % Get similarity between predicted and actual lever
        lever_predict_corr{curr_animal,curr_session} = nan(n_trials,1);
        for i = 1:n_trials
            %curr_corr = corrcoef(predict_lever(:,i), ...
            %    curr_lever{curr_session}(:,i));
            curr_corr = corrcoef(predict_lever(:,i), ...
                lever_reconstructed_session{curr_session}(:,i));
                       
            lever_predict_corr{curr_animal,curr_session}(i) = curr_corr(2);
        end    
    end   
    disp(curr_animal);
end

lever_predict_corr_median = cellfun(@nanmedian,lever_predict_corr);
figure;errorbar(nanmean(lever_predict_corr_median), ...
    nanstd(lever_predict_corr_median)./sqrt(sum( ...
    ~isnan(lever_predict_corr_median))),'k','linewidth',2)
ylabel('Median real/predicted lever correlation');
xlabel('Day')



%% Get / use active cells (in case go away from bone growth, etc)

% Find ROIs that have at least n events per day

active_rois = cell(length(data),1);
for curr_animal = 1:8    
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    n_events = nan(n_rois,14);
    for curr_session = 1:14;       
        n_events(:,curr_session) = ...
            sum(diff([zeros(n_rois,1) data(curr_animal).im(curr_session). ...
            roi_trace_thresh>0],[],2) == 1,2);   
    end
    active_rois{curr_animal} = all(n_events > 5,2);
    disp(curr_animal);
end




%% Predict activity of ROI given activity of all other ROIs?

n_split = 5;
real_predict_corr = nan(length(data),14);

for curr_animal = 1:length(data)
    
    real_active = cell(14,1);
    predict_active = cell(14,1);
    
    for curr_session = 1:14;
        
        curr_movement = classified_rois(curr_animal).movement(:,curr_session);
        curr_quiescent = classified_rois(curr_animal).quiescent(:,curr_session);
        
        curr_active = permute(max(analysis(curr_animal).im(curr_session).move_onset_aligned,[],2),[1,3,2]);
        real_active{curr_session} = curr_active;
      
        n_trials = size(curr_active,1);
        n_split_trials = diff(ceil(linspace(0,n_trials,n_split+1)));
        
        predict_active{curr_session} = nan(size(curr_active));
        for curr_roi = 1:size(curr_active,2)
            
            curr_active_leaveout = curr_active(:,setdiff(1:size(curr_active,2),curr_roi));
            
            split_trials = mat2cell(randperm(n_trials),1,n_split_trials);
            
            for curr_split = 1:n_split
                
                curr_train_trials = horzcat(split_trials{setdiff(1:n_split,curr_split)});
                curr_test_trials = split_trials{curr_split};
                
                curr_weights = +curr_active_leaveout(curr_train_trials,:)\ ...
                    +curr_active(curr_train_trials,curr_roi);
                
                curr_predict = +curr_active_leaveout(curr_test_trials,:)*curr_weights;
                predict_active{curr_session}(curr_test_trials,curr_roi) = curr_predict;
            end
        end
    end
    
    curr_real_predict_corr = cellfun(@(x,y) ...
        corrcoef(+x(:),+y(:)),real_active,predict_active,'uni',false);
    
    real_predict_corr(curr_animal,:) = cellfun(@(x) x(2),curr_real_predict_corr);
    
    disp(curr_animal);
    
end



curr_pre = permute(max(analysis(curr_animal).im(curr_session).move_onset_aligned(:,1:32,:),[],2),[1,3,2]);
curr_move = permute(max(analysis(curr_animal).im(curr_session).move_onset_aligned(:,32:116,:),[],2),[1,3,2]);



%% Correlate pre-move activity with move activity 


framerate = 28;

activity_corr_stages_mean = cell(length(analysis),2);
activity_corr_stages_sem = cell(length(analysis),2);
for curr_animal = 1:length(analysis)
    
    
    % to correspond to movement
    pre_frames = analysis(curr_animal).surrounding_frames < 0;
    move_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*3;

    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    n_rois = size(analysis(curr_animal).im(1).move_onset_aligned,3);
    use_cells_pre = mat2cell(classified_rois(curr_animal).movement, ...
        n_rois,ones(14,1));
    use_cells_move = mat2cell(classified_rois(curr_animal).movement, ...
        n_rois,ones(14,1));
    
    pre_activity_corr = cell(n_sessions,1);
    pre_activity_corr(use_sessions) = cellfun(@(activity,use_cells) ...
        AP_itril(corrcoef(reshape(permute( ...
        max(activity(:,pre_frames,use_cells),[],2), ...
        [2 3 1]),[],size(activity,1)),'rows','pairwise'),-1), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned}, ...
        use_cells_pre,'uni',false);
    
    move_activity_corr = cell(n_sessions,1);
    move_activity_corr(use_sessions) = cellfun(@(activity,use_cells) ...
        AP_itril(corrcoef(reshape(permute( ...
        max(activity(:,move_frames,use_cells),[],2), ...
        [2 3 1]),[],size(activity,1)),'rows','pairwise'),-1), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned}, ...
        use_cells_move,'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = [-1:0.2:1];
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin
    stages = {[1:3] [10:14]};    
    % if no days available from given stage, exclude animal
    if any(~cellfun(@(x) any(ismember(x,1:length(analysis(curr_animal).im))),stages));
        continue
    end       
            
    pre_activity_corr_stages = cellfun(@(x) vertcat(pre_activity_corr{x}),stages,'uni',false);
    move_activity_corr_stages = cellfun(@(x) vertcat(move_activity_corr{x}),stages,'uni',false);
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),pre_activity_corr_stages,'uni',false);
    stages_bins_used = cellfun(@unique,stages_bin_idx,'uni',false);
    
    curr_stage_mean = repmat({nan(length(corr_bin_use),1)},1,2);
    curr_stage_sem = repmat({nan(length(corr_bin_use),1)},1,2);
    for curr_stage = 1:length(stages)
       
        curr_mean = grpstats(move_activity_corr_stages{curr_stage},stages_bin_idx{curr_stage},'mean');
        curr_sem = grpstats(move_activity_corr_stages{curr_stage},stages_bin_idx{curr_stage},'sem');
        
        % NaN correlation gives group 0, ignore
        use_grps = stages_bins_used{curr_stage} > 0;
        
        curr_stage_mean{curr_stage}(stages_bins_used{curr_stage}(use_grps)) = curr_mean(use_grps);
        curr_stage_sem{curr_stage}(stages_bins_used{curr_stage}(use_grps)) = curr_sem(use_grps);
        
    end
    
    activity_corr_stages_mean(curr_animal,:) = curr_stage_mean;
    activity_corr_stages_sem(curr_animal,:) = curr_stage_sem;
    
    disp(curr_animal);
    
end

% Plot stages combine animals
stages_col = {'k','r'};
figure; hold on
for curr_stage = 1:length(stages)
    errorbar(corr_bin_plot,nanmean(horzcat( ...
        activity_corr_stages_mean{:,curr_stage}),2), ...
        nanmean(horzcat( ...
        activity_corr_stages_sem{:,curr_stage}),2), ...
        'color',stages_col{curr_stage});
end


%% (plot stuff for committee meeting)

% avg act over time
avg_act = cell(8,1);
for curr_animal = 1:8;
    for i = 1:14
        avg_act{curr_animal}(:,:,i) = permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned,1),[3 2 1]);
    end
end

a = cellfun(@(x) permute(nanmean(x),[3 2 1]),avg_act,'uni',false);
b = cat(3,a{:});
c = nanmean(b,3);


% get reliability by max time
avg_binary_act = cell(8,1);
for curr_animal = 1:8;
    for i = 1:14
        avg_binary_act{curr_animal}(:,:,i) = ...
            permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned > 0,1),[3 2 1]);
    end
end

[max_reliability,max_reliability_time] = cellfun(@(x) ...
    max(permute(x,[1 3 2]),[],3),avg_binary_act,'uni',false);
reliability_bin = nan(14,121,8);
for curr_animal = 1:8
    for curr_session = 1:14
        a = grpstats(max_reliability{curr_animal}(:,curr_session), ...
            max_reliability_time{curr_animal}(:,curr_session));
        reliability_bin(curr_session, ...
            unique(max_reliability_time{curr_animal}(:,curr_session)),curr_animal) = a;
    end
end
a = nanmedian(reliability_bin,3);
c = ones(1,5)/5;
b = conv2(a,c,'same');
figure; hold on
col = jet(14);
 x = analysis(curr_animal).surrounding_frames/28;
for i = 1:14
    plot(x(2:end),b(i,2:end),'color',col(i,:),'linewidth',2)
end
plot(x(2:end),nanmean(b(:,2:end)),'k','linewidth',2);
xlabel('ROI max time (s)');
ylabel('Average reliability');


% plot by reliability prctile
avg_binary_act = cell(8,1);
for curr_animal = 1:8;
    for i = 1:14
        avg_binary_act{curr_animal}(:,:,i) = ...
            permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned > 0,1),[3 2 1]);
    end
end

[max_reliability,max_reliability_time] = cellfun(@(x) ...
    max(permute(x,[1 3 2]),[],3),avg_binary_act,'uni',false);

reliability_edges = linspace(0,1,10);
reliability_edges(end) = inf;

reliability_bin = cell(8,1);
for curr_animal = 1:8
    for curr_session = 1:14        
      [~,bin] = histc(max_reliability{curr_animal}(:,curr_session),reliability_edges);    
      reliability_bin{curr_animal}(:,curr_session) = bin;
    end
end

reliability_bin_meantrace = nan(14,length(reliability_edges)-1,121,8);
for curr_animal = 1:8
    for curr_session = 1:14        
      for curr_bin = 1:length(reliability_edges)-1;
          curr_rois = reliability_bin{curr_animal}(:,curr_session) == curr_bin;
          if any(curr_rois)            
              reliability_bin_meantrace(curr_session,curr_bin,:,curr_animal) = ...
                  nanmean(avg_act{curr_animal}(curr_rois,:,curr_session),1);
          end
      end
    end
end

a = nanmean(reliability_bin_meantrace,4);
b = permute(nanmean(a,1),[2 3 1]);
figure;
subplot(2,1,1); hold on;
x = analysis(1).surrounding_frames/28;
c = b(any(b,2),:);
col = jet(size(c,1));
for i = 1:size(c,1)
    plot(x,c(i,:),'color',col(i,:),'linewidth',2);    
end
ylabel('Average \DeltaF/F')
xlabel('Time from movement onset (s)')


reliability_bin_meantrace = nan(14,length(reliability_edges)-1,121,8);
for curr_animal = 1:8
    for curr_session = 1:14        
      for curr_bin = 1:length(reliability_edges)-1;
          curr_rois = reliability_bin{curr_animal}(:,curr_session) == curr_bin;
          if any(curr_rois)            
              reliability_bin_meantrace(curr_session,curr_bin,:,curr_animal) = ...
                  nanmean(avg_binary_act{curr_animal}(curr_rois,:,curr_session),1);
          end
      end
    end
end

a = nanmean(reliability_bin_meantrace,4);
b = permute(nanmean(a,1),[2 3 1]);
subplot(2,1,2);hold on;
x = analysis(1).surrounding_frames/28;
c = b(any(b,2),:);
col = jet(size(c,1));
for i = 1:size(c,1)
    plot(x,c(i,:),'color',col(i,:),'linewidth',2);    
end
ylabel('Average reliability')
xlabel('Time from movement onset (s)')












