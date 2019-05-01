%% (ALT, UNUSED) 3C: movement correlation (all movements)

clearvars -except data analysis classified_rois

% Get activity and movements
all_movement = cell(length(data),length(14));

for curr_animal = 1:length(data)
    
    disp(['Animal ' num2str(curr_animal)]);
    
    for curr_session = 1:length(data(curr_animal).im)
        
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
            length(lever_force_resample));
        
        use_movement_start_times(edge_movements) = [];
        use_movement_start_frames(edge_movements) = [];
        
        movements_idx = repmat(use_movement_start_times,1,use_time) + ...
            repmat(0:use_time-1,length(use_movement_start_times),1);
        timed_movements = lever_force_resample(movements_idx');        
        
        all_movement{curr_animal,curr_session} = ...
            timed_movements;
        
        disp(curr_session);
        
    end
end

n_move = cellfun(@(x) size(x,2),all_movement);
all_corr = nan(14,14,length(data));
for curr_animal = 1:length(data)
    
    curr_corr = corrcoef(horzcat(all_movement{curr_animal,:}));
    
    curr_corr_split = ...
        mat2cell(curr_corr,n_move(curr_animal,:),n_move(curr_animal,:));
    
    curr_corr_mean = nan(size(curr_corr_split));
    
    diag_idx = logical(eye(size(curr_corr_split)));
    offdiag_idx = ~eye(size(curr_corr_split));
    
    curr_corr_mean(diag_idx) = ...
        cellfun(@(x) nanmedian(AP_itril(x,-1)),curr_corr_split(diag_idx));
    curr_corr_mean(offdiag_idx) = ...
        cellfun(@(x) nanmedian(x(:)),curr_corr_split(offdiag_idx));
    
    all_corr(:,:,curr_animal) = curr_corr_mean;
    
    disp(curr_animal);
    
end

within_day_corr = cell2mat(arrayfun(@(x) diag(all_corr(:,:,x)),1:length(data),'uni',false))';
adjacent_day_corr = cell2mat(arrayfun(@(x) diag(all_corr(:,:,x),-1),1:length(data),'uni',false))';

figure('Name','Median movement correlation (all movements > 2 sec)');
subplot(1,3,1);
imagesc(nanmean(all_corr,3));
colormap(gray);
xlabel('Day');
ylabel('Day');
axis square
subplot(1,3,2);
errorbar(nanmean(within_day_corr),nanstd(within_day_corr)./sqrt(sum(~isnan(within_day_corr))),'k','linewidth',2)
title('Within-day');
xlabel('Day');
ylabel('Correlation');
subplot(1,3,3);
errorbar(nanmean(adjacent_day_corr),nanstd(adjacent_day_corr)./sqrt(sum(~isnan(adjacent_day_corr))),'k','linewidth',2)
title('Adjacent day');
xlabel('Day');
ylabel('Correlation');

%% (ALT,UNUSED) 4A: example movement/quiescent activity

clearvars -except data analysis classified_rois

use_animal = 2;
use_day = 7;
use_frames = 22500:25020;

% Classified ROIs
curr_m = classified_rois(use_animal).movement(:,use_day);
curr_q = classified_rois(use_animal).quiescent(:,use_day);

% Classified ROI population average
m_avg_act = nanmean(data(use_animal).im(use_day).roi_trace_df(curr_m,use_frames),1);
q_avg_act = nanmean(data(use_animal).im(use_day).roi_trace_df(curr_q,use_frames),1);

% Lever position
move_frames = analysis(use_animal).lever(use_day).lever_move_frames(use_frames);
lever_position = analysis(use_animal).lever(use_day).lever_position_frames(use_frames);

% Cue frames
cue_frames = cell2mat(cellfun(@(x) x.states.cue(1):x.states.cue(end), ...
    data(use_animal).bhv(use_day).bhv_frames,'uni',false)');
cue_frames_use = cue_frames(cue_frames >= use_frames(1) & ...
    cue_frames <= use_frames(end)) - use_frames(1) + 1;
cue_frames_trace = zeros(length(use_frames),1);
cue_frames_trace(cue_frames_use) = 1;


figure; colormap(gray);

subplot(3,1,1);
imagesc(data(use_animal).im(use_day).roi_trace_df(curr_m,use_frames));
caxis([0,5]);
title('Movement-related ROIs');
axis off;

subplot(3,1,2);
imagesc(data(use_animal).im(use_day).roi_trace_df(curr_q,use_frames));
caxis([0,5]);
title('Quiescent-related ROIs');
axis off;

subplot(3,1,3); hold on
plot(cue_frames_trace,'b')
plot((lever_position/max(lever_position))+2,'k')
plot(move_frames+3,'m');
plot(m_avg_act+4,'color',[0,0.7,0]);
plot(q_avg_act+5,'color',[0.7,0,0]);
xlim([0,length(use_frames)])
legend({'Cue','Lever position','Movement', ...
    'Movement-related ROI population average','Quiescent-related ROI population average'});

% Plot scalebar
df_scale = 1;
time_scale = 10; % seconds

line([0,0],[0,df_scale],'color','k','linewidth',2);
line([0,time_scale*data(use_animal).bhv(use_day).framerate],[0,0], ...
    'color','k','linewidth',2);


%% (ALT, UNUSED) 4B: average activity of all neurons aligned to movement onset/offset (CR movements)

clearvars -except data analysis classified_rois

avg_act_onset_all = cell(length(data),1);
avg_act_offset_all = cell(length(data),1);

for curr_animal = 1:length(analysis)

    for curr_day = 1:length(data(curr_animal).im)
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_onset_aligned,1),[3 2 1]);                     
        avg_act_onset_all{curr_animal}(curr_day,:) = nanmean(curr_act(:,:),1);
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_offset_aligned,1),[3 2 1]);                       
        avg_act_offset_all{curr_animal}(curr_day,:) = nanmean(curr_act(:,:),1);
        
    end

end

x = ((1:181) - 91)./28;

figure;
% All days: movement onset aligned
subplot(1,2,1); hold on;
plot(x,nanmean(vertcat(avg_act_onset_all{:})),'k');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (all ROIs, all animals, all days)');

% All days: movement offset aligned
subplot(1,2,2); hold on;
plot(x,nanmean(vertcat(avg_act_offset_all{:})),'k');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (all ROIs, all animals, all days)');

%% (ALT, UNUSED) 4D: average activity of M/Q ROIs aligned to movement onset/offset

clearvars -except data analysis classified_rois

avg_act_onset_m = cell(8,1);
avg_act_onset_q = cell(8,1);
avg_act_offset_m = cell(8,1);
avg_act_offset_q = cell(8,1);
for curr_animal = 1:length(analysis)

    for curr_day = 1:length(data(curr_animal).im)
        
        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_onset_aligned,1),[3 2 1]);                     
        avg_act_onset_m{curr_animal}(curr_day,:) = nanmean(curr_act(m_rois,:),1);
        avg_act_onset_q{curr_animal}(curr_day,:) = nanmean(curr_act(q_rois,:),1);
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_offset_aligned,1),[3 2 1]);                     
        avg_act_offset_m{curr_animal}(curr_day,:) = nanmean(curr_act(m_rois,:),1);
        avg_act_offset_q{curr_animal}(curr_day,:) = nanmean(curr_act(q_rois,:),1);   
        
    end

end

x = ((1:181) - 91)./28;

figure;
% All days: movement onset aligned
subplot(1,2,1); hold on;
plot(x,nanmean(vertcat(avg_act_onset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_onset_q{:})),'color',[0.8,0,0]);
ylabel('Average \DeltaF/F');
legend({'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (all animals, all days)');

% All days: movement offset aligned
subplot(1,2,2); hold on;
plot(x,nanmean(vertcat(avg_act_offset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_offset_q{:})),'color',[0.8,0,0]);
ylabel('Average \DeltaF/F');
legend({'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (all animals, all days)');

%% OLD 6A: Total activity changes

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
 
        move_act_all{curr_animal,curr_day} = move_act_timed;
        
    end
    disp(curr_animal);
end

move_act_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_all,'uni',false);
move_act_m_mean = cell(size(move_act_mean));
move_act_q_mean = cell(size(move_act_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).movement(:,curr_day) == 1)) = 0;
        move_act_m_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
        
        move_act_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).quiescent(:,curr_day) == 1)) = 0;
        move_act_q_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
        
    end
end

% Fill in empty data with nans
empty_data = cellfun(@isempty,move_act_all);
move_act_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);

% Baseline normalize by subtraction
baseline_frames = [-20:-10] + back_frames + 1;

move_act_baseline = cellfun(@(x) nanmean(x(baseline_frames)),move_act_mean);
move_act_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_mean,'uni',false);

move_act_m_baseline_mean = cellfun(@(x) nanmean(x(baseline_frames)),move_act_m_mean);
move_act_m_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean,'uni',false);

move_act_q_baseline_mean = cellfun(@(x) nanmean(x(baseline_frames)),move_act_q_mean);
move_act_q_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean,'uni',false);

% Get (mean) total activity
move_act_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean);
move_act_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean_baselinesub);

move_act_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean);
move_act_m_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean_baselinesub);

move_act_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean);
move_act_q_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean_baselinesub);

% Plot changes in total activity
figure; 

subplot(2,3,1);
curr_data = move_act_m_total_mean_baselinesub;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Movement ROI activity (during movement)');

subplot(2,3,2);
curr_data = move_act_q_total_mean_baselinesub;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Quiescent ROI activity (during movement)');

subplot(2,3,4);
curr_data = move_act_m_baseline_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Movement ROI activity (before movement)');

subplot(2,3,5);
curr_data = move_act_q_baseline_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Quiescent ROI activity (before movement)');

subplot(2,3,3);
curr_data = move_act_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Total activity (during movement)');

subplot(2,3,6);
curr_data = move_act_baseline;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Total activity (before movement)');


% Statistics template
curr_data = move_act_q_baseline_mean;

use_days = 1:4;

curr_days_repmat = repmat(use_days,length(data),1);
curr_data_days = curr_data(:,use_days);
[r,p] = corrcoef(curr_days_repmat(~isnan(curr_data_days)), ...
    curr_data_days(~isnan(curr_data_days)));

p = anova1(move_act_m_baseline_mean)

%% 7B: Activity/movement correlation of movements (across epochs)

clearvars -except data analysis classified_rois

session_epochs = {[1,2,3,4],[5,6,7,8,9,10],[11,12,13,14]};

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
epoch_movement_corr = cellfun(@(allcorr,act1,act2,act3) ...
    mat2cell(allcorr,[size(act1,2),size(act2,2),size(act3,2)], ...
    [size(act1,2),size(act2,2),size(act3,2)]), ...
    epoch_movement_corr_all,epoch_movement_cat(:,1)', ...
    epoch_movement_cat(:,2)',epoch_movement_cat(:,3)','uni',false);
epoch_movement_corr = cellfun(@(x) {x{2,1}(:),x{3,2}(:),x{3,1}(:)}, ...
    epoch_movement_corr,'uni',false);
epoch_movement_corr = vertcat(epoch_movement_corr{:});

epoch_activity_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_activity_cat_reshape{x,:}),'rows','pairwise'), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_activity_corr = cellfun(@(allcorr,act1,act2,act3) ...
    mat2cell(allcorr,[size(act1,2),size(act2,2),size(act3,2)], ...
    [size(act1,2),size(act2,2),size(act3,2)]), ...
    epoch_activity_corr_all,epoch_movement_cat(:,1)', ...
    epoch_movement_cat(:,2)',epoch_movement_cat(:,3)','uni',false);
epoch_activity_corr = cellfun(@(x) {x{2,1}(:),x{3,2}(:),x{3,1}(:)}, ...
    epoch_activity_corr,'uni',false);
epoch_activity_corr = vertcat(epoch_activity_corr{:});

% Bin trials by movement correlation
corr_bin = linspace(-1,1,21);
corr_bin_use = [corr_bin(1:end-1) Inf];
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use);
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
legend({ ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{2})]...
    ['Across sessions ' num2str(session_epochs{2}) ', ' num2str(session_epochs{3})]...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{3})]});




%% 7??: Pairwise activity vs. movement correlation grouped by correlation with template


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


% Get correlation with (median) movement template
template_movement = cellfun(@(move,trials) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Bin trials by correlation with template movement
template_corr_bin = [-Inf -0.5 0.5 Inf];
[~,template_corr_bin_idx] = cellfun(@(x) histc(x,template_corr_bin),movement_template_corr,'uni',false);

% Set up bins for pairwise movement correlation
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

% Make pairwise activity vs. pairwise movement plots within template
% movement correlation bins
epoch_combine_mean = nan(length(corr_bin),3,length(template_corr_bin)-1);
epoch_combine_sem = nan(length(corr_bin),3,length(template_corr_bin)-1);
for curr_template_bin = 1:length(template_corr_bin)-1
    
    use_movement = cellfun(@(move,bins) move(:,bins == curr_template_bin), ...
        epoch_movement_cat,template_corr_bin_idx,'uni',false);
    use_activity = cellfun(@(move,bins) move(:,bins == curr_template_bin), ...
        epoch_activity_cat_reshape,template_corr_bin_idx,'uni',false);
    
    epoch_movement_corr_all = arrayfun(@(x) ...
        corrcoef(horzcat(use_movement{x,:})), ...
        1:size(use_activity,1),'uni',false);
    epoch_movement_corr = cellfun(@(allcorr,act1) ...
        {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
        AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1), ...
        reshape(allcorr(1:size(act1,2),size(act1,2)+1:end),[],1)}, ...
        epoch_movement_corr_all,use_movement(:,1)','uni',false);
    epoch_movement_corr = vertcat(epoch_movement_corr{:});
    
    epoch_activity_corr_all = arrayfun(@(x) ...
        corrcoef(horzcat(use_activity{x,:}),'rows','pairwise'), ...
        1:size(use_activity,1),'uni',false);
    epoch_activity_corr = cellfun(@(allcorr,act1) ...
        {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
        AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1), ...
        reshape(allcorr(1:size(act1,2),size(act1,2)+1:end),[],1)}, ...
        epoch_activity_corr_all,use_activity(:,1)','uni',false);
    epoch_activity_corr = vertcat(epoch_activity_corr{:});
    
    % Bin trials by movement correlation
   
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
        
    for i = 1:3
        epoch_combine_mean(:,i,curr_template_bin) = nanmean(horzcat(actcorr_bin_mean_pad{:,i}),2);
        epoch_combine_sem(:,i,curr_template_bin) = nanstd(horzcat(actcorr_bin_mean_pad{:,i}),[],2)./ ...
            sqrt(sum(~isnan(horzcat(actcorr_bin_mean_pad{:,i})),2));
    end

    
end

figure;
for curr_template_bin = 1:length(template_corr_bin)-1
    subplot(1,length(template_corr_bin)-1,curr_template_bin)
    errorbar(repmat(corr_bin_plot',[1,3,1]), ...
        epoch_combine_mean(:,:,curr_template_bin), ...
        epoch_combine_sem(:,:,curr_template_bin),'linewidth',2);
    legend({'Within early','Within late','Across early/late'});
    xlabel('Pairwise movement correlation');
    ylabel('Pairwise activity correlation');
    legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
        ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);
    title(['Correlation with movement template: ' ...
        num2str(template_corr_bin(curr_template_bin)) ':' ...
        num2str(template_corr_bin(curr_template_bin+1))]);
end

%% 7: test - make plot orthogonal(ish - as good as can get) to 7a

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
    AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1)}, ...
    epoch_movement_corr_all,epoch_movement_cat(:,1)','uni',false);
epoch_movement_corr = vertcat(epoch_movement_corr{:});

epoch_activity_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_activity_cat_reshape{x,:}),'rows','pairwise'), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_activity_corr = cellfun(@(allcorr,act1) ...
    {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
    AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1)}, ...
    epoch_activity_corr_all,epoch_activity_cat_reshape(:,1)','uni',false);
epoch_activity_corr = vertcat(epoch_activity_corr{:});

% Get correlation with (median) movement template
template_movement = cellfun(@(move,trials) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Get mean pairwise template correlation
pairwise_movement_template_corr = ...
    cellfun(@(x) AP_itril(bsxfun(@plus,x,x')./2,-1),movement_template_corr,'uni',false);

% Bin trial pairs by correlation with template movement
template_corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
template_corr_bin = template_corr_bin_use;
template_corr_bin(1) = -1;
template_corr_bin(end) = 1;
template_corr_bin_plot = [template_corr_bin(1:end-1) + diff(template_corr_bin)/2 NaN];
[~,template_corr_bin_idx] = cellfun(@(x) histc(x,template_corr_bin_use), ...
    pairwise_movement_template_corr,'uni',false);

% Bin trials by movement correlation
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),epoch_movement_corr,'uni',false);

% Get mean/sem of activity correlations in each movement correlation bin

[corr_grps,actcorr_bin_mean] = cellfun(@(act,movebin,templatebin) ...
    grpstats(act,{movebin,templatebin},{'gname','nanmean'}), ...
    epoch_activity_corr,movecorr_bin_idx,template_corr_bin_idx,'uni',false);

corr_grps = cellfun(@(x) cellfun(@str2num,x),corr_grps,'uni',false);

actcorr_bin_mean_pad = cellfun(@(x) nan(length(corr_bin), ...
    length(template_corr_bin)),actcorr_bin_mean,'uni',false);

for i = 1:numel(actcorr_bin_mean)    
    for curr_binidx = 1:size(corr_grps{i},1)
        actcorr_bin_mean_pad{i}(corr_grps{i}(curr_binidx,1), ...
            corr_grps{i}(curr_binidx,2)) = actcorr_bin_mean{i}(curr_binidx);       
    end   
end

epoch_combine_mean = nan(length(corr_bin),length(template_corr_bin),2);
epoch_combine_sem = nan(length(corr_bin),length(template_corr_bin),2);
for i = 1:2
    epoch_combine_mean(:,:,i) = nanmean(cat(3,actcorr_bin_mean_pad{:,i}),3);
    epoch_combine_sem(:,:,i) = nanstd(cat(3,actcorr_bin_mean_pad{:,i}),[],3)./ ...
        sqrt(sum(~isnan(cat(3,actcorr_bin_mean_pad{:,i})),3));
end

figure;

subplot(1,2,1); hold on;
set(gca,'ColorOrder',copper(length(template_corr_bin)-1));
errorbar(repmat(template_corr_bin_plot',1,size(epoch_combine_mean,2)), ...
    epoch_combine_mean(:,:,1)',epoch_combine_sem(:,:,1)')
xlabel('Mean pair correlation with template')
ylabel('Pairwise activity correlation');
title(['Sessions ' num2str(session_epochs{1})]);
legend(cellfun(@(x) ['Pairwise movement correlation ' num2str(x)], ...
    num2cell(corr_bin_plot(1:end-1)),'uni',false));

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(length(template_corr_bin)-1));
errorbar(repmat(template_corr_bin_plot',1,size(epoch_combine_mean,2)), ...
    epoch_combine_mean(:,:,2)',epoch_combine_sem(:,:,2)')
xlabel('Mean pair correlation with template')
ylabel('Pairwise activity correlation');
title(['Sessions ' num2str(session_epochs{2})]);
legend(cellfun(@(x) ['Pairwise movement correlation ' num2str(x)], ...
    num2cell(corr_bin_plot(1:end-1)),'uni',false));


%% 7B UNUSED: Activity/movement correlation with template

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
num_rep = 1;

corr_bin = linspace(-1,1,20);
corr_bin_use = [corr_bin(1:end-1) Inf];
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

epoch_combine_mean = nan(length(corr_bin),2,num_rep);
epoch_combine_sem = nan(length(corr_bin),2,num_rep);
for curr_rep = 1:num_rep
    
    % Make template movement/activity from random 50% of trials
    template_trials = [cell(size(epoch_movement_cat,1),1), ...
        cellfun(@(x) sort(randperm(size(x,2),floor(size(x,2)/2))), ...
        epoch_movement_cat(:,2),'uni',false)];
    
    template_movement = cellfun(@(move,trials) nanmean(move(:,trials),2), ...
        epoch_movement_cat(:,2),template_trials(:,2),'uni',false);
    template_activity = cellfun(@(act,trials) nanmean(act(:,trials),2), ...
        epoch_activity_cat_reshape(:,2),template_trials(:,2),'uni',false);
    
    % Get correlation of all trials with template
    movement_template_corr_all = cellfun(@(move,template,trials) ...
        corrcoef([template,move(:,setdiff(1:size(move,2),trials))],'rows','pairwise'), ...
        epoch_movement_cat,repmat(template_movement,1,2),template_trials,'uni',false);
    movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);
    
    activity_template_corr_all = cellfun(@(act,template,trials) ...
        corrcoef([template,act(:,setdiff(1:size(act,2),trials))],'rows','pairwise'), ...
        epoch_activity_cat_reshape,repmat(template_activity,1,2),template_trials,'uni',false);
    activity_template_corr = cellfun(@(x) x(2:end,1),activity_template_corr_all,'uni',false);
    
    % Bin trials by movement correlation   
    [~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);
    
    % Get mean/sem of activity correlations in each movement correlation bin
    actcorr_bin_mean = cellfun(@(act,movebin) grpstats(act,movebin,'nanmean'), ...
        activity_template_corr,movecorr_bin_idx,'uni',false);
    actcorr_bin_sem = cellfun(@(act,movebin) grpstats(act,movebin,'sem'), ...
        activity_template_corr,movecorr_bin_idx,'uni',false);
    
    actcorr_bin_mean_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_mean,'uni',false);
    actcorr_bin_sem_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_sem,'uni',false);
    
    for i = 1:numel(actcorr_bin_mean)
        curr_used_bins = unique(movecorr_bin_idx{i});
        actcorr_bin_mean_pad{i}(curr_used_bins) = actcorr_bin_mean{i};
        actcorr_bin_sem_pad{i}(curr_used_bins) = actcorr_bin_sem{i};
    end
    
    for i = 1:2
        epoch_combine_mean(:,i,curr_rep) = nanmean(horzcat(actcorr_bin_mean_pad{:,i}),2);
        epoch_combine_sem(:,i,curr_rep) = nanstd(horzcat(actcorr_bin_mean_pad{:,i}),[],2)./ ...
            sqrt(sum(~isnan(horzcat(actcorr_bin_mean_pad{:,i})),2));
        
        % Used to use average SEM - doesn't really make sense
        %epoch_combine_sem(:,i) = nanmean(horzcat(actcorr_bin_sem_pad{:,i}),2);
    end
    
    disp(curr_rep);
    
end

figure; hold on;
errorbar(repmat(corr_bin_plot',[1,2]),nanmean(epoch_combine_mean,3), ...
    nanmean(epoch_combine_sem,3),'linewidth',2);
xlabel('Movement template correlation');
ylabel('Activity template correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));





%%%%%% Sanity check that binned movemements have expected pairwise
%%%%%% correlations

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template,trials) ...
    corrcoef([template,move(:,setdiff(1:size(move,2),trials))],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),template_trials,'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);
movement_corr = cellfun(@(x) x(2:end,2:end),movement_template_corr_all,'uni',false);

% Bin trials by movement correlation
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);

movement_corr_within_bin = cell(size(movement_corr));
for i = 1:numel(movement_corr_within_bin)
    for curr_bin = 1:length(corr_bin);
               
        curr_corr = nanmean(AP_itril(movement_corr{i}(movecorr_bin_idx{i} == curr_bin, ...
            movecorr_bin_idx{i} == curr_bin),-1));
        
        movement_corr_within_bin{i}(curr_bin,1) = curr_corr;
        
    end    
end
 
movement_corr_within_bin_mean = nan(length(corr_bin),2);
for i = 1:2
    movement_corr_within_bin_mean(:,i) = ...
        nanmean(horzcat(movement_corr_within_bin{:,i}),2);
end
figure;plot(corr_bin,movement_corr_within_bin_mean);
xlabel('Correlation with template');
ylabel('Pairwise movement correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));


%% TEST: 7c differences in movement/quiescent

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

%%% TEST: use only one class of ROIs

for curr_animal = 1:size(epoch_activity_cat,1)
    for curr_epoch = 1:size(epoch_activity_cat,2)
        
        use_rois = any(classified_rois(curr_animal).unclassified_active(:,session_epochs{curr_epoch}),2);
        curr_act = epoch_activity_cat{curr_animal,curr_epoch}(:,use_rois,:);
        
        epoch_activity_cat_reshape{curr_animal,curr_epoch} = ...
            reshape(curr_act,[],size(curr_act,3));
    end
end

%%%

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
            
            % Skip if there aren't enough trials
            if length(use_trials) < 5
                continue
            end
            
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

few_animals = sum(~isnan(actcorr_bin_mean),3) < 2;
epoch_combine_mean(few_animals) = NaN;
epoch_combine_sem(few_animals) = NaN;

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


%% 7D TEST: Pairwise activity within movement template correlation group

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

all_act_corr = cell(length(data),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Make orthogonal template movement
template_movement_regression = cellfun(@(template,movement) ...
    [template,ones(size(template))]\movement, ...
    template_movement,epoch_movement_cat(:,2),'uni',false);

residual_movement = cellfun(@(template,regression,movement) ...
    movement - [template,ones(size(template))]*regression, ...
    template_movement,template_movement_regression,epoch_movement_cat(:,2),'uni',false);

residual_movement_templates = cellfun(@(x) princomp(x'),residual_movement','uni',false);
ortho_template_movement = cellfun(@(x) x(:,1),residual_movement_templates,'uni',false)';

% %%% TEST: MINIMIZE PC/TEMPLATE SCORE CORR TO GET UNCORR TEMPLATE
% 
% [pc_coeff,pc_score] = cellfun(@(x) princomp(x'),epoch_movement_cat(:,2),'uni',false);
% template_score = cellfun(@(x,m) x'*m,template_movement,epoch_movement_cat(:,2),'uni',false);
% 
% [min_pc_corr,min_pc] = cellfun(@(template_score,pc_score) ...
%     min(abs(1-pdist2(template_score,pc_score(:,1:10)','correlation'))),template_score,pc_score,'uni',false);
% 
% ortho_template_movement = cellfun(@(pc_coeff,min_pc) pc_coeff(:,min_pc), ...
%     pc_coeff,min_pc,'uni',false);
% 
% %%% 

%%% TEST: GET AVERAGE DIFFERENCE BETWEEN PAIRS OF TRIALS?


%%%

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

movement_otemplate_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(ortho_template_movement,1,2),'uni',false);
movement_otemplate_corr = cellfun(@(x) x(2:end,1),movement_otemplate_corr_all,'uni',false);

% Bin trials by movement correlation with template
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);
[~,omovecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_otemplate_corr,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        template_dim_trials = { ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        otemplate_dim_trials = { ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 3)]};        
        
        % Get pairwise movement correlations along these dimensions
        curr_move_corr_templatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        move_corr_templatedim_mean = ...
            nanmean(reshape(curr_move_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_move_corr_otemplatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        move_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_move_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Get pairwise activity correlations along these dimensions
        curr_act_corr_templatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        act_corr_templatedim_mean = ...
            nanmean(reshape(curr_act_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_act_corr_otemplatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        act_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_act_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_otemplatedim_mean];   
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_otemplatedim_mean];      
        
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


%% 7D TEST: compare quadrants (no this will still have the same problem)


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

corr_bin_use = [-Inf -0.0001 0.0001 Inf];
n_bins = length(corr_bin_use);
corr_bin_valid = corr_bin_use;
corr_bin_valid(1) = -1;
corr_bin_valid(end) = 1;
corr_bin_plot = [corr_bin_valid(1:end-1) + diff(corr_bin_valid)/2 NaN];

all_act_corr = cell(length(data),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Make orthogonal template movement
template_movement_regression = cellfun(@(template,movement) ...
    [template,ones(size(template))]\movement, ...
    template_movement,epoch_movement_cat(:,2),'uni',false);

residual_movement = cellfun(@(template,regression,movement) ...
    movement - [template,ones(size(template))]*regression, ...
    template_movement,template_movement_regression,epoch_movement_cat(:,2),'uni',false);

residual_movement_templates = cellfun(@(x) princomp(x'),residual_movement','uni',false);
ortho_template_movement = cellfun(@(x) x(:,1),residual_movement_templates,'uni',false)';

% %%% TEST: MINIMIZE PC/TEMPLATE SCORE CORR TO GET UNCORR TEMPLATE
% 
% [pc_coeff,pc_score] = cellfun(@(x) princomp(x'),epoch_movement_cat(:,2),'uni',false);
% template_score = cellfun(@(x,m) x'*m,template_movement,epoch_movement_cat(:,2),'uni',false);
% 
% [min_pc_corr,min_pc] = cellfun(@(template_score,pc_score) ...
%     min(abs(1-pdist2(template_score,pc_score(:,1:10)','correlation'))),template_score,pc_score,'uni',false);
% 
% ortho_template_movement = cellfun(@(pc_coeff,min_pc) pc_coeff(:,min_pc), ...
%     pc_coeff,min_pc,'uni',false);
% 
% %%% 

%%% TEST: GET AVERAGE DIFFERENCE BETWEEN PAIRS OF TRIALS?


%%%

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

movement_otemplate_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(ortho_template_movement,1,2),'uni',false);
movement_otemplate_corr = cellfun(@(x) x(2:end,1),movement_otemplate_corr_all,'uni',false);

% Bin trials by movement correlation with template
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);
[~,omovecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_otemplate_corr,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        
        
        template_dim_trials = { ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        otemplate_dim_trials = { ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 3)]};      
        
        template_dim_trials = { ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1 & omovecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3 & omovecorr_bin_idx{curr_animal,curr_epoch} == 1)]};
        otemplate_dim_trials = { ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 1 & movecorr_bin_idx{curr_animal,curr_epoch} == 3)], ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 3 & movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        
        
        % Get pairwise movement correlations along these dimensions
        curr_move_corr_templatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        move_corr_templatedim_mean = ...
            nanmean(reshape(curr_move_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_move_corr_otemplatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        move_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_move_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Get pairwise activity correlations along these dimensions
        curr_act_corr_templatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        act_corr_templatedim_mean = ...
            nanmean(reshape(curr_act_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_act_corr_otemplatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        act_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_act_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_otemplatedim_mean];   
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_otemplatedim_mean];      
        
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



%% 7D TEST: make pairwise move corr cutoff in center bin? 
% i.e. no making any specific other axis

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

all_act_corr = cell(length(data),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);


% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Bin trials by movement correlation with template
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);

% Of trials in middle bin, get pairwise correlations
% for now just define cutoff to see if this works
movecorr_cutoff = -0.4;
middlebin_use = cellfun(@(move,bins) ...
    tril(corrcoef(move(:,bins == 2)) < movecorr_cutoff,-1), ...
    epoch_movement_cat,movecorr_bin_idx,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        template_dim_trials = { ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        
        middle_bin_trials = find(movecorr_bin_idx{curr_animal,curr_epoch} == 2);
     
        % Get pairwise movement correlations along template edges
        curr_move_corr_templatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        move_corr_templatedim_mean = ...
            nanmean(reshape(curr_move_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        % Get pairwise movement correlations in selected middle bin trials
        curr_move_corr_middlebin_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            middle_bin_trials),'rows','pairwise');
        move_corr_middlebin_mean = ...
            nanmean(curr_move_corr_middlebin_grid(middlebin_use{curr_animal,curr_epoch}));      
        
        % Get pairwise activity correlations along template edges
        curr_act_corr_templatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        act_corr_templatedim_mean = ...
            nanmean(reshape(curr_act_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        % Get pairwise activity correlations in selected middle bin trials
        curr_act_corr_middlebin_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            middle_bin_trials),'rows','pairwise');
        act_corr_middlebin_mean = ...
            nanmean(curr_act_corr_middlebin_grid(middlebin_use{curr_animal,curr_epoch}));    
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_middlebin_mean];   
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_middlebin_mean];      
        
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


        


%% 7D TEST: differences between movement/quiescent

% Replace any empty data with empty cells for uniform type
empty_data = cellfun(@(x) isempty(x),epoch_movement);
epoch_movement(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);
epoch_activity(empty_data) = cellfun(@(x) {},num2cell(1:sum(empty_data(:))),'uni',false);

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);


for curr_animal = 1:size(epoch_activity_cat,1)
    for curr_epoch = 1:size(epoch_activity_cat,2)
        
        use_rois = any(classified_rois(curr_animal).unclassified_active(:,session_epochs{curr_epoch}),2);
        curr_act = epoch_activity_cat{curr_animal,curr_epoch}(:,use_rois,:);
        
        epoch_activity_cat_reshape{curr_animal,curr_epoch} = ...
            reshape(curr_act,[],size(curr_act,3));
    end
end

% Make template get trial-by-trial activity/movement correlation

corr_bin_use = [-Inf -0.5 0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin_valid = corr_bin_use;
corr_bin_valid(1) = -1;
corr_bin_valid(end) = 1;
corr_bin_plot = [corr_bin_valid(1:end-1) + diff(corr_bin_valid)/2 NaN];

all_act_corr = cell(length(data),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Make orthogonal template movement
template_movement_regression = cellfun(@(template,movement) ...
    [template,ones(size(template))]\movement, ...
    template_movement,epoch_movement_cat(:,2),'uni',false);

residual_movement = cellfun(@(template,regression,movement) ...
    movement - [template,ones(size(template))]*regression, ...
    template_movement,template_movement_regression,epoch_movement_cat(:,2),'uni',false);

residual_movement_templates = cellfun(@(x) princomp(x'),residual_movement','uni',false);
ortho_template_movement = cellfun(@(x) x(:,1),residual_movement_templates,'uni',false)';

% %%% TEST: MINIMIZE PC/TEMPLATE SCORE CORR TO GET UNCORR TEMPLATE
% 
% [pc_coeff,pc_score] = cellfun(@(x) princomp(x'),epoch_movement_cat(:,2),'uni',false);
% template_score = cellfun(@(x,m) x'*m,template_movement,epoch_movement_cat(:,2),'uni',false);
% 
% [min_pc_corr,min_pc] = cellfun(@(template_score,pc_score) ...
%     min(abs(1-pdist2(template_score,pc_score(:,1:10)','correlation'))),template_score,pc_score,'uni',false);
% 
% ortho_template_movement = cellfun(@(pc_coeff,min_pc) pc_coeff(:,min_pc), ...
%     pc_coeff,min_pc,'uni',false);
% 
% %%% 

%%% TEST: GET AVERAGE DIFFERENCE BETWEEN PAIRS OF TRIALS?


%%%

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

movement_otemplate_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(ortho_template_movement,1,2),'uni',false);
movement_otemplate_corr = cellfun(@(x) x(2:end,1),movement_otemplate_corr_all,'uni',false);

% Bin trials by movement correlation with template
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);
[~,omovecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_otemplate_corr,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        template_dim_trials = { ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        otemplate_dim_trials = { ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
            [find(omovecorr_bin_idx{curr_animal,curr_epoch} == 3)]};        
        
        % Get pairwise movement correlations along these dimensions
        curr_move_corr_templatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        move_corr_templatedim_mean = ...
            nanmean(reshape(curr_move_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_move_corr_otemplatedim_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        move_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_move_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Get pairwise activity correlations along these dimensions
        curr_act_corr_templatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [template_dim_trials{1};template_dim_trials{2}]),'rows','pairwise');
        act_corr_templatedim_mean = ...
            nanmean(reshape(curr_act_corr_templatedim_grid(1:length(template_dim_trials{1}), ...
            length(template_dim_trials{1})+1:end),[],1));
        
        curr_act_corr_otemplatedim_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch}(:, ...
            [otemplate_dim_trials{1};otemplate_dim_trials{2}]),'rows','pairwise');
        act_corr_otemplatedim_mean = ...
            nanmean(reshape(curr_act_corr_otemplatedim_grid(1:length(otemplate_dim_trials{1}), ...
            length(otemplate_dim_trials{1})+1:end),[],1));
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_otemplatedim_mean];   
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_otemplatedim_mean];      
        
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


%% 7D TEST: pairwise movecorr cutoff AND template corr diff cutoff

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

all_act_corr = cell(length(data),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

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

% Bin trials by movement correlation with template
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);

% Of trials in middle bin, get pairwise correlations
% for now just define cutoff to see if this works
movecorr_cutoff = -0.4;
movetemp_cutoff = 0.1;
middlebin_use = cellfun(@(movecorr,movetempdiff,bins) ...
    tril(movecorr(bins == 2,bins == 2) < movecorr_cutoff & ...
    movetempdiff(bins == 2,bins == 2) < movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,movecorr_bin_idx,'uni',false);

%%% TEST: don't do this by bin, but rather just by pairwise and template

movecorr_cutoff = -0.4;
movetemp_cutoff = 1;
ortho_movetemp_cutoff = 0.1;

template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movecorr < movecorr_cutoff & ...
    movetempdiff > movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

ortho_template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movecorr < movecorr_cutoff & ...
    movetempdiff < ortho_movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

%%%

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        %         template_dim_trials = { ...
        %             [find(movecorr_bin_idx{curr_animal,curr_epoch} == 1)], ...
        %             [find(movecorr_bin_idx{curr_animal,curr_epoch} == 3)]};
        %
        %         middle_bin_trials = find(movecorr_bin_idx{curr_animal,curr_epoch} == 2);
        
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



%% (ALT, UNUSED) 4E: Temporal changes across days (total activity)

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
                        
        total_act_frame = nanmean(nanmean(curr_act,1),3)';
        total_act_frame_baselinesub = total_act_frame - nanmean(curr_avg_baseline);
        
        total_act_frame_dist = total_act_frame./sum(total_act_frame);
        
        move_act_total(:,curr_day,curr_animal) = total_act_frame_baselinesub;
        move_act_distribution(:,curr_day,curr_animal) = ...
            total_act_frame_dist;
        
    end
    
end

% Plot average across animals, group days
df_scale = 0.05;
time_scale = 29*0.5; % Frames

figure; hold on
day_group = {[1,2,3],[4,5,6],[7,8,9],[10 11 12] [13 14]};
% Total activity
move_act_total_mean = nanmean(move_act_total,3);
move_act_total_grp = nan(size(move_act_total_mean,1),length(day_group));
for i = 1:length(day_group)
   move_act_total_grp(:,i) = nanmean(move_act_total_mean(:,day_group{i}),2); 
end
set(gca,'ColorOrder',copper(length(day_group)));
plot(move_act_total_grp,'linewidth',2);
xlabel('Time');
ylabel('Average \Delta activity per ROI');
title('Movement ROIs')
legend(cellfun(@(x) num2str(x),day_group,'uni',false));
axis square

line([0,0],[0,df_scale],'color','k','linewidth',2);
line([0,time_scale],[0,0],'color','k','linewidth',2);

% Statistics
move_act_total_grp = nan(size(move_act_distribution,1),length(day_group),length(data));
for i = 1:length(day_group)
   move_act_total_grp(:,i,:) = nanmean(move_act_distribution(:,day_group{i},:),2); 
end

move_act_total_sec1 = permute(nanmean(move_act_total_grp(1:10,:,:,1)),[2,3,1]);
move_act_total_sec2 = permute(nanmean(move_act_total_grp(29:end,:,:,1)),[2,3,1]);

d = repmat(transpose(1:length(day_group)),1,length(data));
[re,pe] = corrcoef(d(~isnan(move_act_total_sec1)),move_act_total_sec1(~isnan(move_act_total_sec1)));
[rl,pl] = corrcoef(d(~isnan(move_act_total_sec2)),move_act_total_sec2(~isnan(move_act_total_sec2)));

sec1_act_early = squeeze(nanmean(nanmean(move_act_distribution(1:15,1:4,:),1),2));
sec1_act_late = squeeze(nanmean(nanmean(move_act_distribution(1:15,11:14,:),1),2));

actdistr_early = squeeze(nanmean(nanmean(move_act_distribution(:,1:4,:),2),3));
actdistr_late = squeeze(nanmean(nanmean(move_act_distribution(:,11:14,:),2),3));


% Stats - shuffle test?


actdistr_early = squeeze(nanmean(nanmean(move_act_distribution(:,1:4,:),2),3));
actdistr_late = squeeze(nanmean(nanmean(move_act_distribution(:,11:14,:),2),3));
real_diff = actdistr_early - actdistr_late;

n_shuff = 10000;
shuff_diff = nan(length(real_diff),n_shuff);
for curr_shuff = 1:n_shuff
    move_act_distribution_shuff = shake(move_act_distribution,1);
    actdistr_early = squeeze(nanmean(nanmean(move_act_distribution_shuff(:,1:4,:),2),3));
    actdistr_late = squeeze(nanmean(nanmean(move_act_distribution_shuff(:,11:14,:),2),3));
    shuff_diff(:,curr_shuff) = actdistr_early - actdistr_late;
end

figure; hold on
plot(real_diff,'k');
plot(prctile(shuff_diff,[2.5,97.5],2),'r');
% no this is crap too fails with 5% alpha









