%% Load and prepare data from control animals

animals = {'AP152','AP153','AP154','AP156','AP157','AP159','AP160'};

data = AP_corticospinal_load_all(animals);
analysis = AP_corticospinal_prepare_loaded(data);
classified_rois = AP_classify_movement_cells_continuous(data,analysis);


% AP152 has 2 days missing, add nans to classified and shift data
AP152 = find(cellfun(@(x) strcmp('AP152',x),{data(:).animal}));

if ~isempty(AP152) && size(classified_rois(AP152).movement,2) < 14
    
    classified_rois(AP152).movement = +classified_rois(AP152).movement;
    classified_rois(AP152).quiescent = +classified_rois(AP152).quiescent;
    
    classified_rois(AP152).movement(:,4:14) = classified_rois(AP152).movement(:,2:12);
    classified_rois(AP152).quiescent(:,4:14) = classified_rois(AP152).quiescent(:,2:12);
    
    classified_rois(AP152).movement(:,2:3) = NaN;
    classified_rois(AP152).quiescent(:,2:3) = NaN;
    
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

%% Get statistics of movement (all movements)

clearvars -except data analysis classified_rois

move_stats = struct( ...
    'frac_moving',cell(size(data)), ...
    'move_push_movefrac',cell(size(data)), ...
    'move_push_globalfrac',cell(size(data)), ...
    'move_duration',cell(size(data)), ...
    'quiescent_duration',cell(size(data)), ...
    'move_amplitude',cell(size(data)));

for curr_animal = 1:length(data)
    
    move_stats(curr_animal).frac_moving = nan(1,14);
    move_stats(curr_animal).move_push_movefrac = nan(1,14);
    move_stats(curr_animal).move_push_globalfrac = nan(1,14);
    move_stats(curr_animal).move_duration = nan(1,14);
    move_stats(curr_animal).quiescent_duration = nan(1,14);
    move_stats(curr_animal).move_amplitude = nan(1,14);
    
    for curr_session = 1:length(data(curr_animal).im)
        
        % Parse movement/nonmovement
        lever_move = data(curr_animal). ...
            bhv(curr_session).imaged_downsampled_lever_force;
        % (if no data, skip - this happened in one control animal)
        if isempty(lever_move)
            continue
        end
        [lever_active,~,~,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated(lever_move,true);
        
        % Split movements/nonmovements
        boundary_times = find(diff([Inf;lever_active;Inf]) ~= 0);
        lever_move_split = mat2cell(lever_move,diff(boundary_times));
        lever_envelope_split = mat2cell(lever_velocity_envelope_smooth,diff(boundary_times));

        lever_move_split_move = cellfun(@any,mat2cell(lever_active,diff(boundary_times)));
        
               
        % fraction time spent moving
        frac_moving = nanmean(lever_active);
        
        % when moving, fraction time spent above or below prior rest
        move_push_movefrac = nanmedian(cellfun(@(x) nanmean(x > x(1)),lever_move_split(lever_move_split_move)));
        move_push_globalfrac = sum(cellfun(@(x) sum(x > x(1)), ...
            lever_move_split(lever_move_split_move)))./sum(cellfun(@length, ...
            lever_move_split(lever_move_split_move)));

        % duration of move/quiescent epochs
        move_duration = nanmedian(cellfun(@length, ...
            lever_move_split(lever_move_split_move)));
        
        quiescent_duration = nanmedian(cellfun(@length, ...
            lever_move_split(~lever_move_split_move)));
        
        % average amplitude of movements
        move_amplitude = nanmedian(cellfun(@nanmedian,lever_envelope_split(lever_move_split_move)));
        
        move_stats(curr_animal).frac_moving(curr_session) = frac_moving;
        move_stats(curr_animal).move_push_movefrac(curr_session) = move_push_movefrac;
        move_stats(curr_animal).move_push_globalfrac(curr_session) = move_push_globalfrac;
        move_stats(curr_animal).move_duration(curr_session) = move_duration;
        move_stats(curr_animal).quiescent_duration(curr_session) = quiescent_duration;
        move_stats(curr_animal).move_amplitude(curr_session) = move_amplitude;
        
    end
    disp(['Move stats: animal ' num2str(curr_animal)]);
end

% If AP152 is included (which has days 2-3 missing), fix
AP152 = find(cellfun(@(x) strcmp('AP152',x),{data(:).animal}));
if ~isempty(AP152)
   move_stats(AP152).frac_moving(:,4:14) = move_stats(AP152).frac_moving(:,2:12);
   move_stats(AP152).frac_moving(:,2:3) = NaN;
   
   move_stats(AP152).move_push_movefrac(:,4:14) = move_stats(AP152).move_push_movefrac(:,2:12);
   move_stats(AP152).move_push_movefrac(:,2:3) = NaN;
   
   move_stats(AP152).move_push_globalfrac(:,4:14) = move_stats(AP152).move_push_globalfrac(:,2:12);
   move_stats(AP152).move_push_globalfrac(:,2:3) = NaN;
   
   move_stats(AP152).move_duration(:,4:14) = move_stats(AP152).move_duration(:,2:12);
   move_stats(AP152).move_duration(:,2:3) = NaN;
   
   move_stats(AP152).quiescent_duration(:,4:14) = move_stats(AP152).quiescent_duration(:,2:12);
   move_stats(AP152).quiescent_duration(:,2:3) = NaN;
   
   move_stats(AP152).move_amplitude(:,4:14) = move_stats(AP152).move_amplitude(:,2:12);
   move_stats(AP152).move_amplitude(:,2:3) = NaN;
end

% Plot each feature across days

figure; 

subplot(2,3,1);
curr_stat = vertcat(move_stats.frac_moving);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Fraction of time moving')

subplot(2,3,2);
curr_stat = vertcat(move_stats.move_push_movefrac);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Fraction of movement pushing (vs. pulling) within movements')

subplot(2,3,3);
curr_stat = vertcat(move_stats.move_push_globalfrac);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Fraction of movement pushing (vs. pulling) whole day')

subplot(2,3,4);
curr_stat = vertcat(move_stats.move_duration);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Movement duration')

subplot(2,3,5);
curr_stat = vertcat(move_stats.quiescent_duration);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Quiescence duration')

subplot(2,3,6);
curr_stat = vertcat(move_stats.move_amplitude);
errorbar(nanmean(curr_stat),nanstd(curr_stat)./ ...
    sqrt(sum(~isnan(curr_stat))),'k','linewidth',2);
title('Move amplitude')



%% Align activity to onset of movement

for curr_animal = 1:length(data)
    for curr_session = 1:length(data(curr_animal).im)
        
        lever_active = analysis(curr_animal).lever(curr_session).lever_move_frames;
        
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 15;
        max_move_time = Inf;
        min_pre_move_iti = 10;
        min_post_move_iti = 10;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        
        frames_back = 90;
        frames_forward = 90;
        valid_frames = (use_movement_start_frames - frames_back > 0) & ...
            (use_movement_start_frames + frames_forward < size(data(curr_animal).im(curr_session).roi_trace_df,2));
        
        
        curr_aligned_act = cell2mat(cellfun(@(x) permute(data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
            x-frames_back:x+frames_forward),[3,2,1]),num2cell( ...
            use_movement_start_frames(valid_frames)),'uni',false));
        
        analysis(curr_animal).im(curr_session).allmove_onset_aligned = ...
            curr_aligned_act;
        
    end
end



%% Plot fraction of movement/quiescence ROIs

% Fraction of movement/quiescent cells
move_frac_raw = cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false);
quiescent_frac_raw = cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false);

% Pad any smaller arrays to 14 days with NaNs
move_frac = cell2mat(cellfun(@(x) padarray(x,[0,14-length(x)], ...
    NaN,'post'),move_frac_raw,'uni',false));
quiescent_frac = cell2mat(cellfun(@(x) padarray(x,[0,14-length(x)], ...
    NaN,'post'),quiescent_frac_raw,'uni',false));

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


%% Plot the fraction of classified ROIs on one day overlapping with another

m2m_overlap_grid = nan(14,14,length(data));
q2q_overlap_grid = nan(14,14,length(data));
m2q_overlap_grid = nan(14,14,length(data));
q2m_overlap_grid = nan(14,14,length(data));
for curr_animal = 1:length(data)
   for curr_day1 = 1:size(classified_rois(curr_animal).movement,2);
       for curr_day2 = 1:size(classified_rois(curr_animal).movement,2);
           % m2m
           curr_overlap = +classified_rois(curr_animal).movement(:,curr_day1)'* ...
               +classified_rois(curr_animal).movement(:,curr_day2);
           curr_total = sum(classified_rois(curr_animal).movement(:,curr_day1));
           
           m2m_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % q2q
           curr_overlap = +classified_rois(curr_animal).quiescent(:,curr_day1)'* ...
               +classified_rois(curr_animal).quiescent(:,curr_day2);
           curr_total = sum(classified_rois(curr_animal).quiescent(:,curr_day1));
           
           q2q_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % m2q
           curr_overlap = +classified_rois(curr_animal).movement(:,curr_day1)'* ...
               +classified_rois(curr_animal).quiescent(:,curr_day2);
           curr_total = sum(classified_rois(curr_animal).movement(:,curr_day1));
           
           m2q_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % q2m
           curr_overlap = +classified_rois(curr_animal).quiescent(:,curr_day1)'* ...
               +classified_rois(curr_animal).movement(:,curr_day2);
           curr_total = sum(classified_rois(curr_animal).quiescent(:,curr_day1));
           
           q2m_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
       end
   end    
end

figure;
subplot(2,2,1);
imagesc(nanmean(m2m_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('M2M')
subplot(2,2,2);
imagesc(nanmean(q2q_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('Q2Q')
subplot(2,2,3);
imagesc(nanmean(m2q_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('M2Q')
subplot(2,2,4);
imagesc(nanmean(q2m_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('Q2M')


%% Activity-movement analysis: pull out activity/movements by duration
% Only use movements not overlapping with reward, "ITI movements"

%session_epochs = {[1:4],[11:14]};
session_epochs = num2cell(1:14);

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
            use_time = 1000;
            use_frames = (use_time/1000)*28;
            
            min_move_time = 1000;
            max_move_time = Inf;
            min_pre_move_iti = 500;
            min_post_move_iti = 500;
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
            
            % Don't use any movements that overlap with reward/punish
            reward_trials = cellfun(@(x) ~isempty(x.states.reward), ...
                data(curr_animal).bhv(curr_session).bhv_frames);
            
            reward_frames = cell2mat(cellfun(@(x) x.states.reward(1):x.states.reward(2), ...
                data(curr_animal).bhv(curr_session).bhv_frames(reward_trials),'uni',false)');
            punish_frames = cell2mat(cellfun(@(x) x.states.punish(1):x.states.punish(2), ...
                data(curr_animal).bhv(curr_session).bhv_frames(~reward_trials),'uni',false)');
            event_frames = [reward_frames,punish_frames];
            
            frames_idx = repmat(use_movement_start_frames,1,use_frames) + ...
                repmat(0:use_frames-1,length(use_movement_start_frames),1);
            reward_movements = any(ismember(frames_idx,event_frames),2);
            
            use_movement_start_times(edge_movements | reward_movements) = [];
            use_movement_start_frames(edge_movements | reward_movements) = [];
            
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
corr_bin = linspace(-1,1,20);
corr_bin_use = [corr_bin(1:end-1) Inf];
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
    epoch_combine_sem(:,i) = nanmean(horzcat(actcorr_bin_sem_pad{:,i}),2);
end

figure; hold on;
errorbar(epoch_combine_mean,epoch_combine_sem);
legend({'Within early','Within late','Across early/late'});
set(gca,'XTick',1:2:length(corr_bin_plot));
set(gca,'XTickLabel',round(corr_bin_plot(1:2:end)*10)/10);
xlabel('Movement correlation');
ylabel('Activity correlation');















