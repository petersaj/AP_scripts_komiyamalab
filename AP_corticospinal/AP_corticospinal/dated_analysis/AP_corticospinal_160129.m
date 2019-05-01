%% Average activity of all neurons aligned to movement onset/offset

avg_act_onset_m = cell(8,1);
avg_act_onset_q = cell(8,1);
avg_act_onset_all = cell(8,1);
avg_act_offset_m = cell(8,1);
avg_act_offset_q = cell(8,1);
avg_act_offset_all = cell(8,1);
for curr_animal = 1:length(analysis)

    for curr_day = 1:length(data(curr_animal).im)
        
        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_onset_aligned,1),[3 2 1]);                     
        avg_act_onset_m{curr_animal}(curr_day,:) = nanmean(curr_act(m_rois,:),1);
        avg_act_onset_q{curr_animal}(curr_day,:) = nanmean(curr_act(q_rois,:),1);
        avg_act_onset_all{curr_animal}(curr_day,:) = nanmean(curr_act(:,:),1);
        
        curr_act = permute(nanmean(analysis(curr_animal).im(curr_day).move_offset_aligned,1),[3 2 1]);                     
        avg_act_offset_m{curr_animal}(curr_day,:) = nanmean(curr_act(m_rois,:),1);
        avg_act_offset_q{curr_animal}(curr_day,:) = nanmean(curr_act(q_rois,:),1);   
        avg_act_offset_all{curr_animal}(curr_day,:) = nanmean(curr_act(:,:),1);
        
    end

end

x = ((1:181) - 91)./28;

figure;
% All days: movement onset aligned
subplot(4,2,1); hold on;
plot(x,nanmean(vertcat(avg_act_onset_all{:})),'k');
plot(x,nanmean(vertcat(avg_act_onset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_onset_q{:})),'r');
ylabel('Average \DeltaF/F');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset');

% All days: movement offset aligned
subplot(4,2,2); hold on;
plot(x,nanmean(vertcat(avg_act_offset_all{:})),'k');
plot(x,nanmean(vertcat(avg_act_offset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_offset_q{:})),'r');
ylabel('Average \DeltaF/F');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset');

% Within days: movement onset aligned (movement ROIs)
subplot(4,2,3); hold on;
linecolor= [linspace(0,0,14);linspace(0,0.7,14);linspace(0,0,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_onset_all_dayavg = nanmean(cat(3,avg_act_onset_m{:}),3);
plot(x,avg_act_onset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (movement ROIs)');

% Within days: movement offset aligned (movement ROIs)
subplot(4,2,4); hold on;
linecolor= [linspace(0,0,14);linspace(0,0.7,14);linspace(0,0,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_offset_all_dayavg = nanmean(cat(3,avg_act_offset_m{:}),3);
plot(x,avg_act_offset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (movement ROIs)');

% Within days: movement onset aligned (quiescent ROIs)
subplot(4,2,5); hold on;
linecolor= [linspace(0,0.7,14);linspace(0,0,14);linspace(0,0,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_onset_all_dayavg = nanmean(cat(3,avg_act_onset_q{:}),3);
plot(x,avg_act_onset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (quiesent ROIs)');

% Within days: movement offset aligned (quiescent ROIs)
subplot(4,2,6); hold on;
linecolor= [linspace(0,0.7,14);linspace(0,0,14);linspace(0,0,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_offset_all_dayavg = nanmean(cat(3,avg_act_offset_q{:}),3);
plot(x,avg_act_offset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (quiescent ROIs)');

% Within days: movement onset aligned (all ROIs)
subplot(4,2,7); hold on;
linecolor= [linspace(0,0.7,14);linspace(0,0.7,14);linspace(0,0.7,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_onset_all_dayavg = nanmean(cat(3,avg_act_onset_all{:}),3);
plot(x,avg_act_onset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (movement ROIs)');

% Within days: movement offset aligned (all ROIs)
subplot(4,2,8); hold on;
linecolor= [linspace(0,0.7,14);linspace(0,0.7,14);linspace(0,0.7,14)]';
set(gca,'ColorOrder',linecolor);
avg_act_offset_all_dayavg = nanmean(cat(3,avg_act_offset_all{:}),3);
plot(x,avg_act_offset_all_dayavg');
ylabel('Average \DeltaF/F');
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (movement ROIs)');



%% Same classification overlap (corrected for # classified)

% Same classification overlap
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m2m = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + repmat(sum(m,1),size(m,2),1)),m,'uni',false);
m2m_cat = cat(3,m2m{:});
m2m_tril = cell2mat(arrayfun(@(x) [nanmean(AP_itril(m2m_cat(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(m2m_cat(8:end,8:end,x),-1))], ...
    transpose(1:size(m2m_cat,3)),'uni',false));

q2q = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + repmat(sum(q,1),size(q,2),1)),q,'uni',false);
q2q_cat = cat(3,q2q{:});
q2q_tril = cell2mat(arrayfun(@(x) [nanmean(AP_itril(q2q_cat(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(q2q_cat(8:end,8:end,x),-1))], ...
    transpose(1:size(q2q_cat,3)),'uni',false));

% Same classification overlap shuffle
n_rep = 1000;

m2m_cat_shuff = nan(14,14,length(data),n_rep);
q2q_cat_shuff = nan(14,14,length(data),n_rep);

m_use = cellfun(@(x) x(any(x,2),:),m,'uni',false);
q_use = cellfun(@(x) x(any(x,2),:),q,'uni',false);

for i = 1:n_rep
    m_shuff = cellfun(@(x) shake(x,1),m_use,'uni',false);    
    m2m_shuff = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + repmat(sum(m,1),size(m,2),1)),m_shuff,'uni',false);   
    m2m_cat_shuff(:,:,:,i) = cat(3,m2m_shuff{:});
    
    q_shuff = cellfun(@(x) shake(x,1),q_use,'uni',false);    
    q2q_shuff = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + repmat(sum(q,1),size(q,2),1)),q_shuff,'uni',false);   
    q2q_shuff_cat = cat(3,q2q_shuff{:});
    q2q_cat_shuff(:,:,:,i) = cat(3,q2q_shuff{:});
 
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
p_m2m = signrank(m2m_tril(:,1),m2m_tril(:,2));
p_q2q = signrank(q2q_tril(:,1),q2q_tril(:,2));


%% Changes in distribution of activity across movement (different ways)

figure; 

% Get average activity within movement ROIs
back_time = 0;
forward_frames = 29*3;

m_avg_act = cell(length(data),14);
for curr_animal = 1:length(data)        
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < forward_frames;    
    for curr_session = 1:length(data(curr_animal).im);               
        m = find(classified_rois(curr_animal).movement(:,curr_session));    
        m_avg_act{curr_animal,curr_session} = ...
            permute(nanmean(analysis(curr_animal).im(curr_session). ...
            move_onset_aligned(:,use_frames,m) > 0,1),[3,2,1]);
    end
end
session_grp = {[1:4],[6:9],[11:14]};
for i = 1:length(session_grp);
    subplot(3,4,i);
    curr_cat = vertcat(m_avg_act{:,session_grp{i}});
    [~,max_idx] = max(curr_cat,[],2);
    [~,sort_idx] = sort(max_idx);
    imagesc(curr_cat(sort_idx,:));colormap(hot);  
    xlabel('Frames from movement onset');
    ylabel('Sorted movement ROI')
    title(num2str(session_grp{i}));
end

% Get average activity across ROIs
back_time = 0;
forward_frames = 29*3;

m_avg_act_avg = nan(14,86,length(data));
for curr_animal = 1:length(data)        
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < forward_frames;    
    for curr_session = 1:length(data(curr_animal).im);               
        m = find(classified_rois(curr_animal).movement(:,curr_session));    
        m_avg_act_avg(curr_session,:,curr_animal) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session). ...
            move_onset_aligned(:,use_frames,m) > 0,1),[3,2,1]));
    end
end
subplot(3,4,4);
imagesc(nanmean(m_avg_act_avg,3));colormap(hot);
xlabel('Frames from movement onset');
ylabel('Day');
title('Average activity across movement ROIs');

% Mean of activity start timing for CR movements
act_starts_mean_all = cell(length(data),1);
for curr_animal = 1:length(data);
    
    use_frames = analysis(curr_animal).surrounding_frames >= 0 & analysis(curr_animal).surrounding_frames <= 90;

    act_starts = cell(14,size(data(curr_animal).im(1).roi_trace_df,1));
    
    for i = 1:length(data(curr_animal).im);
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

subplot(3,2,3); hold on;
errorbar(timing_bins_plot,nanmean(early_timing_bins),nanstd(early_timing_bins)./ ...
    sqrt(sum(~isnan(early_timing_bins))),'k','linewidth',2);
errorbar(timing_bins_plot,nanmean(late_timing_bins),nanstd(late_timing_bins)./ ...
    sqrt(sum(~isnan(late_timing_bins))),'r','linewidth',2);

legend({'Early (days 1:3)' 'Late (days 11:14)'},'location','se');
ylabel('Cumulative fraction');
xlabel('Time (s)')
title('Distribution of activity onset mean')


% Maximum of activity-aligned movement onset
forward_frames = 29*3;

m_avg_act_max = cell(length(data),14);

for curr_animal = 1:length(data)        
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < forward_frames;    
    for curr_session = 1:length(data(curr_animal).im);               
        m = find(classified_rois(curr_animal).movement(:,curr_session));    
        m_avg_act = permute(nanmean(analysis(curr_animal).im(curr_session). ...
            move_onset_aligned(:,use_frames,m),1),[3,2,1]);
        [~,m_avg_act_max{curr_animal,curr_session}] = max(m_avg_act,[],2);
    end
end

subplot(3,2,4); hold on;
session_grp = {[1:4],[6:9],[11:14]};
col = [linspace(0,1,length(session_grp))',zeros(length(session_grp),2)];
for i = 1:length(session_grp);
    curr_cat = vertcat(m_avg_act_max{:,session_grp{i}});
    plot(sort(curr_cat)./29,linspace(0,1,length(curr_cat)),'color',col(i,:),'linewidth',2);    
end
legend(cellfun(@num2str,session_grp,'uni',false));
ylabel('Cumulative fraction');
xlabel('Time from movement onset (s)');
title('Distribution of frame with maximum average activity');


% Distribution of activity across movement (non-cumulative)
forward_frames = 29*3;

m_total_act = cell(length(data),14);

for curr_animal = 1:length(data)        
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < forward_frames;    
    for curr_session = 1:length(data(curr_animal).im);               
        m = find(classified_rois(curr_animal).movement(:,curr_session));    
        m_roi_act = permute(nansum(analysis(curr_animal).im(curr_session). ...
            move_onset_aligned(:,use_frames,m) > 0,1),[3,2,1]);
        m_total_act{curr_animal,curr_session} = sum(m_roi_act,1)./sum(sum(m_roi_act,1));        
    end
end

subplot(3,2,5); hold on;
session_grp = {[1:4],[6:9],[11:14]};
col = [linspace(0,1,length(session_grp))',zeros(length(session_grp),2)];
for i = 1:length(session_grp);
    curr_cat = vertcat(m_total_act{:,session_grp{i}});
    plot((1:size(curr_cat,2))/29,nanmean(curr_cat,1),'color',col(i,:),'linewidth',2);    
end
legend(cellfun(@num2str,session_grp,'uni',false));
ylabel('Fraction of total activity');
xlabel('Time from movement onset (s)');
title('Distribution activity during movement');

% Distribution of activity from normalized onset to offset all movements
move_act_normdist_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
        % Get  activity
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:) > 0);
        
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
        min_move_time = 15;
        max_move_time = Inf;
        min_pre_move_iti = 10;
        min_post_move_iti = 10;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);

        % Get activity during movement, and +/- time equal to movement time
        surround_move_frames = cellfun(@(start,stop) ...
            (start-round((stop-start)/2)):(stop+round((stop-start)/2)), ...
            num2cell(use_movement_start_frames), ...
            num2cell(use_movement_stop_frames),'uni',false);
        
        bad_frames = cellfun(@(x) any(x < 0) | any(x > size(curr_act,2)), ...
            surround_move_frames);
        
        surround_move_frames(bad_frames) = [];
        
        move_act_m = cellfun(@(x) nanmean(curr_act(curr_m,x),1), ...
            surround_move_frames,'uni',false);
        
        % Get distribution of activity across movement, resampled
        n_t_bins = 60;
        move_act_m_resample = cellfun(@(x) resample(x,n_t_bins,length(x)), ...
            move_act_m,'uni',false);
        move_act_m_resample_cat = vertcat(move_act_m_resample{:});
        move_act_normdist = bsxfun(@times,move_act_m_resample_cat,1./sum(move_act_m_resample_cat,2));
        
        move_act_normdist_all{curr_animal,curr_day} = nanmean(move_act_normdist,1);
        
    end
    disp(curr_animal);
end

subplot(3,2,6); hold on;
session_grp = {[1:4],[6:9],[11:14]};
col = [linspace(0,1,length(session_grp))',zeros(length(session_grp),2)];
for i = 1:length(session_grp);
    curr_cat = vertcat(move_act_normdist_all{:,session_grp{i}});
    plot(linspace(-1,2,size(curr_cat,2)),nanmean(curr_cat,1),'color',col(i,:),'linewidth',2);    
end
line([0,0],ylim,'linestyle','--','color','k');
line([1,1],ylim,'linestyle','--','color','k');
legend(cellfun(@num2str,session_grp,'uni',false));
xlabel('Normalized time from movement onset');
ylabel('Normalized activity');
title('Distribution of activity across all time-normalized movements');


%% Average activity vs. average lever
    
avg_move = nan(5001,14,length(data));
m_avg_act = nan(181,14,length(data));
for curr_animal = 1:length(data)

    for curr_session = 1:14;
        cat_move = horzcat(analysis(curr_animal).lever(curr_session).rewarded_movement_fixtime{:});
        avg_move(:,curr_session,curr_animal) = nanmean(bsxfun(@minus,cat_move,nanmean(cat_move(500:900,:),1)),2);
        
        m = find(classified_rois(curr_animal).movement(:,curr_session));
        m_avg_act(:,curr_session,curr_animal) = ...
            nanmean(nanmean(analysis(curr_animal).im(curr_session). ...
            move_onset_aligned(:,:,m) > 0,1),3)';
    end
    
    act_t = (-90:90)./29;
    move_t = (-1000:4000)/1000;
    
end

figure;
for curr_animal = 1:length(data)
    subplot(2,4,curr_animal);hold on;
    set(gca,'ColorOrder',jet(14));
    curr_data = bsxfun(@times,m_avg_act(:,:,curr_animal),1./nanmean(m_avg_act(1:50,:,curr_animal),1));
    plot(curr_data);
end

curr_animal = 2;
for i = 1:14
    figure('Position',[3,392,571,583]);
    subplot(2,1,1);
    plot(move_t,avg_move(:,i,curr_animal),'k')
    ylim([-0.3,0.3]);
    xlim([-0.5,3]);
    subplot(2,1,2);
    plot(act_t,m_avg_act(:,i,curr_animal),'k')
    ylim([0,0.2]);
    xlim([-0.5,3]);
end

%% Average activity of movement-aligned ROIs across days (for defense)

avg_act_m = cell(length(data),14);
back_frames = 28*1.5;
forward_frames = 28*3;

for curr_animal = 1:length(data)
    
    use_rois = any(classified_rois(curr_animal).movement,2);
    curr_move_rois = classified_rois(curr_animal).movement(use_rois,:);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        use_onset_frames = analysis(curr_animal).surrounding_frames > -back_frames & ...
            analysis(curr_animal).surrounding_frames < forward_frames;
        
        use_offset_frames = analysis(curr_animal).surrounding_frames < forward_frames & ...
            analysis(curr_animal).surrounding_frames > -back_frames;
        
        curr_onset_act = permute(nanmean(analysis(curr_animal). ...
            im(curr_day).move_onset_aligned(:,use_onset_frames,use_rois),1),[3,2,1]);
        
        curr_offset_act = permute(nanmean(analysis(curr_animal). ...
            im(curr_day).move_onset_aligned(:,use_offset_frames,use_rois),1),[3,2,1]);     
        
        curr_onset_act(~curr_move_rois(:,curr_day),:) = 0;  
        
        avg_act_m{curr_animal,curr_day} = curr_onset_act;
        
    end 
end

avg_act_m_animalcat = cell(1,14);
for i = 1:14
   avg_act_m_animalcat{i} = vertcat(avg_act_m{:,i}); 
end

days_combine = {[1:4],[6:9],[11:14]};
avg_act_m_dayavg = cell(1,length(days_combine));
for i = 1:length(days_combine)
    avg_act_m_dayavg{i} = nanmean(cat(3,avg_act_m_animalcat{days_combine{i}}),3);
end

max_act = max(cat(3,avg_act_m_dayavg{:}),[],3);
[~,max_idx] = max(max_act,[],2);
[~,sort_idx] = sort(max_idx);

figure; 
for i = 1:length(days_combine)
    subplot(1,length(days_combine),i);
    imagesc(avg_act_m_dayavg{i}(sort_idx,:)); colormap(hot);
    title(num2str(days_combine{i}));
    xlabel('Frames');
    ylabel('Sorted ROI');
    line([back_frames+1,back_frames+1],ylim,'color','w','linewidth',2);
end

%% Positional preference for quiescent ROIs

% Get average position and velocity when each ROI is active
qroi_position = cell(length(data),1);
qroi_position_shuff = cell(length(data),1);

for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    qroi_position{curr_animal} = cell(n_days,1);
    qroi_position_shuff{curr_animal} = cell(n_days,1);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = smooth(analysis(curr_animal). ...
            lever(curr_day).lever_position_frames,100,'loess');
        
        % Normalize lever parameters
        lever_position_norm = lever_position_frames./std(lever_position_frames);
        
        % Real data
        curr_q = classified_rois(curr_animal).quiescent(:,curr_day);
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        
        % Only use "active" ROIs, at least 5 events
        min_events = 5;
        active_rois = sum(diff(curr_act > 0,[],2) == 1,2) > 5;
        curr_act = curr_act(active_rois,:);
        
        curr_roi_position = (curr_act*lever_position_norm)./sum(curr_act,2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_position_split = mat2cell(lever_position_norm,diff(boundary_frames));
        
        num_rep = 1000;
        shuffle_perms = shake(repmat([1:length(boundary_frames)-1]',1,num_rep),1);
        lever_position_shuffle = nan(length(lever_position_norm),num_rep);
        for i = 1:num_rep
            lever_position_shuffle(:,i) = ...
                vertcat(lever_position_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        curr_roi_position_shuff = bsxfun(@times,curr_act*lever_position_shuffle,1./sum(curr_act,2));

        % Store
        qroi_position{curr_animal}{curr_day} = curr_roi_position;  
        qroi_position_shuff{curr_animal}{curr_day} = curr_roi_position_shuff(:);
        
    end
    disp(curr_animal);
end

figure; hold on;

% Histogram of real vs. shuffled position (all concatenated)
qroi_position_cat = cell2mat((cellfun(@(x) vertcat(x{:}),qroi_position,'uni',false)));
qroi_position_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),qroi_position_shuff,'uni',false)));

bin_edges = linspace(-2,2,50);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(qroi_position_cat,bin_centers)./length(qroi_position_cat);
shuff_dist = hist(qroi_position_shuff_cat,bin_centers)./length(qroi_position_shuff_cat);

plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
legend({'Shuffled','Real'})

ylabel('Number of ROIs');
xlabel('Lever displacement from median (stds)');

%% Positional preference for quiescent ROIs (shuffle from quiescence only)

% Get average position and velocity when each ROI is active
qroi_position = cell(length(data),1);
qroi_position_shuff = cell(length(data),1);

for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    qroi_position{curr_animal} = cell(n_days,1);
    qroi_position_shuff{curr_animal} = cell(n_days,1);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = smooth(analysis(curr_animal). ...
            lever(curr_day).lever_position_frames,100,'loess');
        
        % Normalize lever parameters
        lever_position_norm = lever_position_frames./std(lever_position_frames);
        
        % Real data
        curr_q = classified_rois(curr_animal).quiescent(:,curr_day);
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        
        % Split data, remove movement epochs
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_position_split = mat2cell(lever_position_norm,diff(boundary_frames));
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_frames));
        act_split = mat2cell(curr_act,size(curr_act,1),diff(boundary_frames));
        
        lever_move_epoch = cellfun(@any,lever_active_split);
        
        act_split_q = act_split(~lever_move_epoch);
        lever_position_q = lever_position_split(~lever_move_epoch);
        
        act_q_cat = horzcat(act_split_q{:});
        lever_position_q_cat = vertcat(lever_position_q{:});
        
        % Only use "active" ROIs, at least 5 events
        min_events = 5;
        active_rois = sum(diff(act_q_cat > 0,[],2) == 1,2) > 5;
        act_q_cat = act_q_cat(active_rois,:);
        
        curr_roi_position = (act_q_cat*lever_position_q_cat)./sum(act_q_cat,2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)        
        num_rep = 1000;
        shuffle_perms = shake(repmat([1:length(lever_position_q)]',1,num_rep),1);
        lever_position_shuffle = nan(length(lever_position_q_cat),num_rep);
        for i = 1:num_rep
            lever_position_shuffle(:,i) = ...
                vertcat(lever_position_q{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        curr_roi_position_shuff = bsxfun(@times,act_q_cat*lever_position_shuffle,1./sum(act_q_cat,2));

        % Store
        qroi_position{curr_animal}{curr_day} = curr_roi_position;  
        qroi_position_shuff{curr_animal}{curr_day} = curr_roi_position_shuff(:);
        
    end
    disp(curr_animal);
end

figure; hold on;

% Histogram of real vs. shuffled position (all concatenated)
qroi_position_cat = cell2mat((cellfun(@(x) vertcat(x{:}),qroi_position,'uni',false)));
qroi_position_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),qroi_position_shuff,'uni',false)));

bin_edges = linspace(-2,2,50);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(qroi_position_cat,bin_centers)./length(qroi_position_cat);
shuff_dist = hist(qroi_position_shuff_cat,bin_centers)./length(qroi_position_shuff_cat);

plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
legend({'Shuffled','Real'})

ylabel('Number of ROIs');
xlabel('Lever displacement from median (stds)');



%% Classify qROIs by position (shuffle only quiescent epochs)


for curr_animal = 1:length(data);
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    n_days = length(data(curr_animal).im);
    
    q_position_class(curr_animal).position_up = false(n_rois,n_days);
    q_position_class(curr_animal).position_down = false(n_rois,n_days);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = analysis(curr_animal).lever(curr_day).lever_position_frames;
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);    
        
        % Split data, remove movement epochs
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_position_split = mat2cell(lever_position_frames,diff(boundary_frames));
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_frames));
        act_split = mat2cell(curr_act,size(curr_act,1),diff(boundary_frames));
        
        lever_move_epoch = cellfun(@any,lever_active_split);
        
        act_split_q = act_split(~lever_move_epoch);
        lever_position_q = lever_position_split(~lever_move_epoch);
        
        roi_position = (horzcat(act_split_q{:})*vertcat(lever_position_q{:}))./ ...
            sum(horzcat(act_split_q{:}),2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)
   
        num_rep = 100;
        shuffle_perms = shake(repmat([1:length(lever_position_q)]',1,num_rep),1);
        lever_position_shuffle = nan(length(vertcat(lever_position_q{:})),num_rep);
        for i = 1:num_rep
            lever_position_shuffle(:,i) = ...
                vertcat(lever_position_q{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        roi_position_shuff = bsxfun(@times,horzcat(act_split_q{:})* ...
            lever_position_shuffle,1./sum(horzcat(act_split_q{:}),2));
        
        roi_position_rank = tiedrank([roi_position,roi_position_shuff]')';
        roi_position_p = roi_position_rank(:,1)./(num_rep+1);        
        
        % Store
        q_position_class(curr_animal).position_down(:,curr_day) = ...
            roi_position_p > 0.975;
        
        q_position_class(curr_animal).position_up(:,curr_day) = ...
            roi_position_p < 0.025;
        
    end
    disp(curr_animal);
end


any_q_positionclass = cellfun(@(x,y) x|y, ...
    {q_position_class(:).position_down}, ...
    {q_position_class(:).position_up},'uni',false);

q_position_frac = cellfun(@(p,q) sum(p&q,1)./sum(q,1), ...
    any_q_positionclass,{classified_rois(:).quiescent},'uni',false);

q_position_up_frac = cellfun(@(p,q) sum(p&q,1)./sum(q,1), ...
    {q_position_class(:).position_up},{classified_rois(:).quiescent},'uni',false);

q_position_down_frac = cellfun(@(p,q) sum(p&q,1)./sum(q,1), ...
    {q_position_class(:).position_down},{classified_rois(:).quiescent},'uni',false);

figure;
subplot(1,3,1); hold on
plot(vertcat(q_position_frac{:})','k');
plot(nanmean(vertcat(q_position_frac{:}))','r','linewidth',2);
xlabel('Day');
ylabel('Fraction of Q rois that are position-selective');
title('Significant relative to quiescent epochs only');

subplot(1,3,2); hold on;
plot(vertcat(q_position_up_frac{:})','k');
plot(nanmean(vertcat(q_position_up_frac{:}))','r','linewidth',2);
xlabel('Day');
ylabel('Fraction of Q rois that are position up-selective');
title('Significant relative to quiescent epochs only');

subplot(1,3,3); hold on;
plot(vertcat(q_position_down_frac{:})','k');
plot(nanmean(vertcat(q_position_down_frac{:}))','r','linewidth',2);
xlabel('Day');
ylabel('Fraction of Q rois that are position down-selective');
title('Significant relative to quiescent epochs only');

%% Temporal activity during all movements, speed binned

% Distribution of activity from normalized onset to offset all movements
move_act_m_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
        % Get  activity
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:) > 0);
        
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
        min_move_time = 15;
        max_move_time = Inf;
        min_pre_move_iti = 10;
        min_post_move_iti = 10;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        use_movement_start_frames = movement_starts(use_movements);
        use_movement_stop_frames = movement_stops(use_movements);
        
        % Alternate: timed activity
        back_frames = 28*1;
        forward_frames = 28*3;
        total_frames = back_frames+forward_frames+1;
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        move_act_m_timed = permute(reshape(curr_act(curr_m, ...
            move_frame_timed_idx)',total_frames,[],sum(curr_m)),[2,1,3]);
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_m_all{curr_animal,curr_day} = move_act_m_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

move_act_m_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_m_all,'uni',false);
move_position_mean = cellfun(@(x) nanmean(x,1),move_position_all,'uni',false);
move_speed_mean = cellfun(@(x) nanmean(x,1),move_speed_all,'uni',false);

% Plot average activity/position/speed over days
figure; 

move_act_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_m_day(:,curr_day) = nanmean(vertcat(move_act_m_mean{:,curr_day}),1);
    
end
move_act_m_day_minusbaseline = bsxfun(@minus,move_act_m_day,mean(move_act_m_day(1:20,:),1));
subplot(1,4,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_m_day_minusbaseline);
ylabel('Baseline-subtracted binarized activity');
xlabel('Frames');

move_position_day = nan(total_frames,14);
for curr_day = 1:14
    move_position_day(:,curr_day) = nanmean(vertcat(move_position_mean{:,curr_day}),1);
    
end
move_position_day_minusbaseline = bsxfun(@minus,move_position_day,mean(move_position_day(1:20,:),1));
subplot(1,4,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_position_day_minusbaseline);
ylabel('Baseline-subtracted movement');
xlabel('Frames');

move_speed_day = nan(total_frames,14);
for curr_day = 1:14
    move_speed_day(:,curr_day) = nanmean(vertcat(move_speed_mean{:,curr_day}),1);
    
end
move_speed_day_minusbaseline = bsxfun(@minus,move_speed_day,mean(move_speed_day(1:20,:),1));
subplot(1,4,3); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_speed_day_minusbaseline);
ylabel('Baseline-subtracted speed');
xlabel('Frames');

% Plot max speed vs. max activity
max_move_act_m_minusbaseline = nan(size(move_act_m_mean));
max_move_speed_minusbaseline = nan(size(move_act_m_mean));
use_animaldays = cellfun(@(x) ~isempty(x),move_act_m_mean);

max_move_act_m_minusbaseline(use_animaldays) = ...
    cellfun(@(x) max(x-nanmean(x(1:20))),move_act_m_mean(use_animaldays));
max_move_speed_minusbaseline(use_animaldays) = ...
    cellfun(@(x) max(x-nanmean(x(1:20))),move_speed_mean(use_animaldays));

max_move_act_m_day_minusbaseline = max(move_act_m_day_minusbaseline,[],1);
max_move_speed_day_minusbaseline = max(move_speed_day_minusbaseline,[],1);

subplot(1,4,4); hold on;
plot(max_move_speed_minusbaseline,max_move_act_m_minusbaseline,'.')
plot(max_move_speed_day_minusbaseline,max_move_act_m_day_minusbaseline,'^k','MarkerSize',10)
xlabel('Movement-ROI activity');
ylabel('Speed');


% Bin trials by speed within the first second
use_frames = back_frames+1:back_frames+1+28;
bin_param = cellfun(@(x) max(x(:,:),[],2),move_speed_all,'uni',false);
speed_bins = linspace(min(vertcat(bin_param{:})), ...
    prctile(vertcat(bin_param{:}),95),20);
speed_bins(end) = Inf;
[~,initial_speed_binned] = cellfun(@(x) histc(x,speed_bins),bin_param,'uni',false);

act_speed_binned = cell(size(move_act_m_all,1),size(move_act_m_all,2),length(speed_bins)-1);
speed_bin_mean = nan(size(move_act_m_all,1),size(move_act_m_all,2),length(speed_bins)-1);
for curr_bin = 1:length(speed_bins)-1
    
    act_speed_binned(:,:,curr_bin) = ...
        cellfun(@(act,bins) nanmean(nanmean(act(bins == curr_bin,:,:),3),1), ...
        move_act_m_all,initial_speed_binned,'uni',false);
    
    speed_bin_mean(:,:,curr_bin) = ...
        cellfun(@(speed,bins) nanmean(speed(bins == curr_bin,:)), ...
        bin_param,initial_speed_binned);
    
end

act_speed_binned_day = nan(total_frames,14,length(speed_bins)-1);
for curr_day = 1:14
    for curr_bin = 1:length(speed_bins)-1
        act_speed_binned_day(:,curr_day,curr_bin) = ...
            nanmean(vertcat(act_speed_binned{:,curr_day,curr_bin}),1)';    
    end
end

% Plot by day, colored by speed bin
figure; 
sq = ceil(sqrt(length(speed_bins)-1));
for curr_bin = 1:length(speed_bins)-1
    subplot(sq,sq,curr_bin); hold on;
    set(gca,'ColorOrder',jet(14));
    plot(act_speed_binned_day(:,:,curr_bin));
    title(['Speed bin: ' num2str(curr_bin)]);
end

% Plot by speed bin, colored by day
figure; 
for curr_day = 1:14
    subplot(2,7,curr_day); hold on;
    set(gca,'ColorOrder',jet(length(speed_bins)-1));
    plot(permute(act_speed_binned_day(:,curr_day,:),[1,3,2]));
    title(['Day: ' num2str(curr_day)]);
end



%% Activity during all (long enough) movements, and corresponding movements (full)

% Get activity and movements
used_movement = cell(length(data),14);
used_activity = cell(length(data),14);

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
            back_time = 1000;
            forward_time = 3000;
            
            back_frames = (back_time/1000)*28;
            forward_frames = (forward_time/1000)*28;
            
            min_move_time = 1000;
            max_move_time = 5000;
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
            edge_movements = (use_movement_start_times + (forward_time) > ...
                length(lever_force_resample)) | ...
                (use_movement_start_frames + (forward_frames) > ...
                size(data(curr_animal).im(curr_session).roi_trace_df,2)) & ...
                (use_movement_start_times - (back_time) < 1) | ...
                (use_movement_start_frames - (back_time) < 1);
            
            % Don't use any movements that overlap with reward/punish
%             reward_trials = cellfun(@(x) ~isempty(x.states.reward), ...
%                 data(curr_animal).bhv(curr_session).bhv_frames);
%             
%             reward_frames = cell2mat(cellfun(@(x) x.states.reward(1):x.states.reward(2), ...
%                 data(curr_animal).bhv(curr_session).bhv_frames(reward_trials),'uni',false)');
%             punish_frames = cell2mat(cellfun(@(x) x.states.punish(1):x.states.punish(2), ...
%                 data(curr_animal).bhv(curr_session).bhv_frames(~reward_trials),'uni',false)');
%             event_frames = [reward_frames,punish_frames];
%             
%             frames_idx = repmat(use_movement_start_frames,1,forward_frames) + ...
%                 repmat(0:forward_frames-1,length(use_movement_start_frames),1);
%             reward_movements = any(ismember(frames_idx,event_frames),2);
%             
%             use_movement_start_times(edge_movements | reward_movements) = [];
%             use_movement_start_frames(edge_movements | reward_movements) = [];
            
            use_movement_start_times(edge_movements) = [];
            use_movement_start_frames(edge_movements) = [];
            
            movements_idx = repmat(use_movement_start_times,1,1+back_time+forward_time) + ...
                repmat(-back_time:forward_time,length(use_movement_start_times),1);
            timed_movements = lever_force_resample(movements_idx');
            
            timed_activity = permute(cell2mat(arrayfun(@(x) ...
                data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
                use_movement_start_frames(x)-back_frames: ...
                use_movement_start_frames(x)+forward_frames), ...
                permute(1:length(use_movement_start_frames),[3,1,2]),'uni',false)),[2,1,3]);
            
            used_movement{curr_animal,curr_session} = ...
                timed_movements;
            used_activity{curr_animal,curr_session} = ...
                timed_activity;
            
            disp(curr_session);
            
        end
end







%% Prepare/save events for quality control (must be done on non-grouped)

save_path = ['/usr/local/lab/People/Andy/Data/corticospinal_qc'];

for curr_animal = 1:8  
    ca_events = cell(length(data(curr_animal).im),1);
    for curr_session = 1:length(data(curr_animal).im)
        curr_data = data(curr_animal).im(curr_session).roi_trace_thresh;
        
        curr_data_eventmarkers = diff((curr_data > 0),[],2);
        
        event_onsets = arrayfun(@(x) find(curr_data_eventmarkers(x,:) == 1), ...
            1:size(curr_data_eventmarkers,1),'uni',false);
        event_offsets = arrayfun(@(x) find(curr_data_eventmarkers(x,:) == -1), ...
            1:size(curr_data_eventmarkers,1),'uni',false);
        
        ca_events{curr_session} = cellfun(@(on,off) [on;off],event_onsets,event_offsets,'uni',false);
        
        curr_savename = [save_path filesep data(curr_animal).animal '_events.mat'];
        
        save(curr_savename,'ca_events');       
    end
    disp(['Finished animal ' data(curr_animal).animal]);
end

disp('done')


%% Variability across movement

% Activity
act_pop_var = cell(length(data),14);
act_pop_avg = cell(length(data),14);
for curr_animal = 1:length(analysis)
    
    for curr_day = 1:length(data(curr_animal).im)
        
        curr_data = analysis(curr_animal).im(curr_day).move_onset_aligned > 0;

        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);   
        
        curr_data_m = curr_data(:,:,m_rois);
        
        % don't use trials with any NaNs (where did these come from?)
        nan_trials = any(any(isnan(curr_data_m),2),3);
        curr_data_m = curr_data_m(~nan_trials,:,:);
        
        curr_pop_var = nan(size(curr_data,1),1);
        for curr_frame = 1:size(curr_data_m,2);            
            curr_frame_data = zscore(permute(curr_data_m(:,curr_frame,:),[3,1,2]),[],2);
            curr_pop_var(curr_frame) = trace(cov(curr_frame_data));
        end
        
        act_pop_var{curr_animal,curr_day} = curr_pop_var;
        
        curr_data_m_norm = nan(size(curr_data_m));
        for i = 1:size(curr_data_m,3)
            curr_data_m_norm(:,:,i) = reshape(zscore(reshape( ...
                curr_data_m(:,:,i),[],1)),size(curr_data_m,1),size(curr_data_m,2));
        end
        act_pop_avg{curr_animal,curr_day} = nanmean(nanmean(curr_data_m,3),1);

    end
    disp(curr_animal);
end

% Plot animal-mean variance over days
act_pop_var_day = nan(181,14);
for i = 1:14
    act_pop_var_day(:,i) = nanmean(horzcat(act_pop_var{:,i}),2);
end
figure; hold on
set(gca,'ColorOrder',jet(14));
plot(act_pop_var_day);
ylabel('Population variance');
xlabel('Frame');
title('Movement ROIs');

% Plot animal-mean variance over days
act_pop_avg_day = nan(181,14);
for i = 1:14
    act_pop_avg_day(:,i) = nanmean(vertcat(act_pop_avg{:,i}),1)';
end
figure; hold on
set(gca,'ColorOrder',jet(14));
plot(act_pop_avg_day);
ylabel('Population mean');
xlabel('Frame');
title('Movement ROIs');

% Movement


%% Correlation of ALL movements of a certain duration

% Get activity and movements
all_movement = cell(length(data),length(14));
all_activity = cell(length(data),length(14));

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
        use_frames = (use_time/1000)*28;
        
        min_move_time = 2000;
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
        
        use_movement_start_times(edge_movements) = [];
        use_movement_start_frames(edge_movements) = [];
        
        movements_idx = repmat(use_movement_start_times,1,use_time) + ...
            repmat(0:use_time-1,length(use_movement_start_times),1);
        timed_movements = lever_force_resample(movements_idx');
        
        timed_activity = permute(cell2mat(arrayfun(@(x) ...
            data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
            use_movement_start_frames(x):use_movement_start_frames(x)+(use_frames-1)), ...
            permute(1:length(use_movement_start_frames),[3,1,2]),'uni',false)),[2,1,3]);
        
        all_movement{curr_animal,curr_session} = ...
            timed_movements;
        all_activity{curr_animal,curr_session} = ...
            timed_activity;
        
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
        cellfun(@(x) nanmean(AP_itril(x,-1)),curr_corr_split(diag_idx));
    curr_corr_mean(offdiag_idx) = ...
        cellfun(@(x) nanmean(x(:)),curr_corr_split(offdiag_idx));
    
    all_corr(:,:,curr_animal) = curr_corr_mean;
    
    disp(curr_animal);
    
end

within_day_corr = cell2mat(arrayfun(@(x) diag(all_corr(:,:,x),-1),1:length(data),'uni',false))';
adjacent_day_corr = cell2mat(arrayfun(@(x) diag(all_corr(:,:,x),-1),1:length(data),'uni',false))';

figure;
subplot(1,3,1);
imagesc(nanmean(all_corr,3));
colormap(gray);
title('All movements');
xlabel('Day');
ylabel('Day');
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

%% Get example lever presses (all movements)

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
        back_time = 1000;
        forward_time = 2000;
        
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
        
        % Don't use any movements that don't have full timing
        edge_movements = (use_movement_start_times+forward_time > ...
            length(lever_force_resample)) & ...
            (use_movement_start_times-back_time < 0);
        
        use_movement_start_times(edge_movements) = [];
        use_movement_start_frames(edge_movements) = [];
        
        movements_idx = repmat(use_movement_start_times,1,back_time+forward_time+1) + ...
            repmat(-back_time:forward_time,length(use_movement_start_times),1);
        timed_movements = lever_force_resample(movements_idx');
        
        all_movement{curr_animal,curr_session} = ...
            timed_movements;
        
        disp(curr_session);
        
    end
end

% Subtract baseline for each movement
all_movement_qsub = cellfun(@(x) bsxfun(@minus,x,nanmean(x(1:back_time-10,:),1)),all_movement,'uni',false);
all_movement_mean = cellfun(@(x) nanmean(x,2),all_movement_qsub,'uni',false);

% Plot mean movements of all animals
figure; 
for curr_animal = 1:length(data)
    
    subplot(3,3,curr_animal); hold on;
    set(gca,'ColorOrder',jet(14));
    
    t = -back_time:forward_time;
    plot(t,horzcat(all_movement_mean{curr_animal,:}),'linewidth',2);    
    
    axis square
    xlim([-back_time,forward_time]);
    title(['Animal ' num2str(curr_animal)])
    
end

% Plot all movements for one animal across days
plot_animal = 1;
figure('Name',['Animal ' num2str(plot_animal)]);
for curr_day = 1:14
    
    subplot(2,7,curr_day); hold on;
    
    t = -back_time:forward_time;
    
    plot(t,all_movement_qsub{plot_animal,curr_day},'color',[0.5 0.5 0.5]);
    plot(t,all_movement_mean{plot_animal,curr_day},'k','linewidth',2);
    
    axis square
    xlim([-back_time,forward_time]);
    
end

% Plot random subset of movements for one animal across days
plot_animal = 1;
figure('Name',['Animal ' num2str(plot_animal)]);
for curr_day = 1:14
    
    subplot(2,7,curr_day); hold on;
    
    n_plot = 5;
    curr_rand_plot = randperm(size(all_movement_qsub{plot_animal,curr_day},2),n_plot);
    t = -back_time:forward_time;
    
    plot(t,all_movement_qsub{plot_animal,curr_day}(:,curr_rand_plot),'color',[0.5 0.5 0.5]);
    plot(t,all_movement_mean{plot_animal,curr_day},'k','linewidth',2);
    
    axis square
    xlim([-back_time,forward_time]);
    
end

% Plot mean movement with std for one animal across days
plot_animal = 1;
figure('Name',['Animal ' num2str(plot_animal)]);
for curr_day = 1:14
    
    subplot(2,7,curr_day); hold on;
    
    curr_mean = all_movement_mean{plot_animal,curr_day};
    curr_std = nanstd(all_movement_qsub{plot_animal,curr_day},[],2);
    t = -back_time:forward_time;
    
    plot(t,curr_mean-curr_std,'color',[0.5 0.5 0.5]);
    plot(t,curr_mean+curr_std,'color',[0.5 0.5 0.5]);
    plot(t,curr_mean,'k','linewidth',2);
    
    axis square
    xlim([-back_time,forward_time]);
    
end


%% Temporal & total activity during all (long enough) movements

% Distribution of activity from normalized onset to offset all movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        
        % To zero out all non-classified activity
%         curr_class = (classified_rois(curr_animal).movement(:,curr_day) == 1) | ...
%             (classified_rois(curr_animal).quiescent(:,curr_day) == 1);
%         move_act_timed(:,:,~curr_class) = 0;
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end


% Set all zeros to NaN (to isolate amplitude component)
move_act_amplitude_all = move_act_all;
for i = 1:numel(move_act_amplitude_all)
    % Use this for amplitude of all discrete events
    [m,n,p] = size(move_act_amplitude_all{i});
    event_onsets = [zeros(m,1,p),diff(move_act_amplitude_all{i} > 0,[],2) == 1];
    move_act_amplitude_all{i}(~event_onsets) = NaN;
    
    % Use this for amplitude of all active frames
    %move_act_amplitude_all{i}(move_act_amplitude_all{i} == 0) = NaN;
end

move_act_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_all,'uni',false);
% To use only classified cells on each day to build total 
% move_act_mean = cell(size(move_act_mean));
% for curr_animal = 1:length(data)
%     for curr_day = 1:length(data(curr_animal).im)
%         use_cells = (classified_rois(curr_animal).movement(:,curr_day) == 1) | ...
%             (classified_rois(curr_animal).quiescent(:,curr_day) == 1);
%         curr_data = move_act_all{curr_animal,curr_day};
%         curr_data(:,:,~use_cells) = 0;
%         move_act_mean{curr_animal,curr_day} = ...
%             nanmean(nanmean(curr_data,1),3);
%     end
% end
move_act_binary_mean = cellfun(@(x) nanmean(nanmean(x > 0,1),3),move_act_all,'uni',false);
move_act_amplitude_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_amplitude_all,'uni',false);

move_act_m_mean = cell(size(move_act_mean));
move_act_m_mean_allfrac = cell(size(move_act_mean));
move_act_binary_m_mean = cell(size(move_act_mean));
move_act_amplitude_m_mean = cell(size(move_act_mean));

move_act_q_mean = cell(size(move_act_mean));
move_act_q_mean_allfrac = cell(size(move_act_mean));
move_act_binary_q_mean = cell(size(move_act_mean));
move_act_amplitude_q_mean = cell(size(move_act_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im)
        
        move_act_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).movement(:,curr_day) == 1)) = 0;
        move_act_m_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
                
        move_act_binary_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        move_act_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).quiescent(:,curr_day) == 1)) = 0;
        move_act_q_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
        
        move_act_binary_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
    end
end

% Fill in empty data with nans
empty_data = cellfun(@isempty,move_act_all);
move_act_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean_allfrac(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean_allfrac(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_all(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);

% Baseline normalize by subtraction
baseline_frames = [-20:-10] + back_frames + 1;

move_act_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_mean,'uni',false);
move_act_binary_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_mean,'uni',false);
move_act_amplitude_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_mean,'uni',false);

move_act_m_baseline_mean = cellfun(@(x) nanmean(x(baseline_frames)),move_act_m_mean);
move_act_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean,'uni',false);
move_act_m_mean_allfrac = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean_allfrac,'uni',false);
move_act_binary_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_m_mean,'uni',false);
move_act_amplitude_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_m_mean,'uni',false);

move_act_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean,'uni',false);
move_act_q_mean_allfrac = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean_allfrac,'uni',false);
move_act_binary_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_q_mean,'uni',false);
move_act_amplitude_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_q_mean,'uni',false);

% Get (mean) total activity
move_act_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean);
move_act_binary_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_mean);
move_act_fracrois = cellfun(@(x) nanmean(nanmean(any(x(:,back_frames+1:end,:),2),3)),move_act_all);

move_act_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean);
move_act_m_allfrac_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean_allfrac);
move_act_binary_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_m_mean);

move_act_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean);
move_act_q_allfrac_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean_allfrac);
move_act_binary_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_q_mean);

% Get mean amplitude (average across trial, not frame)
move_act_amplitude_total_mean = cellfun(@(x) ...
    nanmean(nanmean(nanmean(x(:,back_frames+1:end,:),2),1),3),move_act_amplitude_all);

move_act_amplitude_total_mean_baseline = cellfun(@(x) ...
    nanmean(nanmean(nanmean(x(:,baseline_frames,:),2),1),3),move_act_amplitude_all);

move_act_amplitude_m_total_mean = nan(size(move_act_amplitude_all));
move_act_amplitude_m_total_mean_baseline = nan(size(move_act_amplitude_all));

move_act_amplitude_q_total_mean = nan(size(move_act_amplitude_all));
move_act_amplitude_q_total_mean_baseline = nan(size(move_act_amplitude_all));

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im)
        
        move_act_amplitude_m_total_mean(curr_animal,curr_day) = ...
            nanmean(nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,back_frames+1:end, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),2),1),3);
        
        move_act_amplitude_m_total_mean_baseline(curr_animal,curr_day) = ...
            nanmean(nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,baseline_frames, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),2),1),3);
        
        move_act_amplitude_q_total_mean(curr_animal,curr_day) = ...
            nanmean(nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,back_frames+1:end, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),2),1),3);
        
        move_act_amplitude_q_total_mean_baseline(curr_animal,curr_day) = ...
            nanmean(nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,baseline_frames, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),2),1),3);
        
    end
end

move_act_amplitude_total_mean = move_act_amplitude_total_mean - ...
    move_act_amplitude_m_total_mean_baseline;
move_act_amplitude_m_total_mean = move_act_amplitude_m_total_mean - ...
    move_act_amplitude_m_total_mean_baseline;
move_act_amplitude_q_total_mean = move_act_amplitude_q_total_mean - ...
    move_act_amplitude_q_total_mean_baseline;

% Plot average temporal activity/position/speed over days

move_position_mean = cellfun(@(x) nanmean(x,1),move_position_all,'uni',false);
move_speed_mean = cellfun(@(x) nanmean(x,1),move_speed_all,'uni',false);

figure; 

move_act_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_day(:,curr_day) = nanmean(vertcat(move_act_mean{:,curr_day}),1);
end
subplot(3,3,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_m_day(:,curr_day) = nanmean(vertcat(move_act_m_mean{:,curr_day}),1);
end
subplot(3,3,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_m_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_q_day(:,curr_day) = nanmean(vertcat(move_act_q_mean{:,curr_day}),1);
end
subplot(3,3,3); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_q_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_day(:,curr_day) = nanmean(vertcat(move_act_binary_mean{:,curr_day}),1);
end
subplot(3,3,4); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_m_day(:,curr_day) = nanmean(vertcat(move_act_binary_m_mean{:,curr_day}),1);
end
subplot(3,3,5); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_m_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_q_day(:,curr_day) = nanmean(vertcat(move_act_binary_q_mean{:,curr_day}),1);
end
subplot(3,3,6); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_q_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_mean{:,curr_day}),1);
end
subplot(3,3,7); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_m_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_m_mean{:,curr_day}),1);
end
subplot(3,3,8); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_m_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_q_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_q_mean{:,curr_day}),1);
end
subplot(3,3,9); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_q_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');


% Plot average total activity

figure; 

subplot(3,3,1);
curr_data = move_act_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Total activity');

subplot(3,3,2);
curr_data = move_act_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Movement ROI activity');

subplot(3,3,3);
curr_data = move_act_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Quiescent ROI activity');

subplot(3,3,4);
curr_data = move_act_binary_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity');
title('Total activity');

subplot(3,3,5);
curr_data = move_act_binary_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity');
title('Movement ROI activity');

subplot(3,3,6);
curr_data = move_act_binary_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity');
title('Quiescent ROI activity');

subplot(3,3,7);
curr_data = move_act_amplitude_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity');
title('Total activity');

subplot(3,3,8);
curr_data = move_act_amplitude_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity');
title('Movement ROI activity');

subplot(3,3,9);
curr_data = move_act_amplitude_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity');
title('Quiescent ROI activity');


% Plot fraction of active ROIs per movement
figure;

subplot(1,2,1);
errorbar(nanmean(move_act_fracrois,1), ...
    nanstd(move_act_fracrois,[],1)./sqrt(sum(~isnan(move_act_fracrois),1)),'k','linewidth',2)
xlabel('Day');
ylabel('Fraction of ROIs active during movement');

move_act_binary_m_notime = nan(size(move_act_all));
for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im)

        move_act_binary_m_notime(curr_animal,curr_day) = ...
            nanmean(nanmean(any(move_act_all{curr_animal,curr_day}(:,back_frames+1:end, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1) > 0,2),3),1);
      
    end
end
subplot(1,2,2);
errorbar(nanmean(move_act_binary_m_notime,1), ...
    nanstd(move_act_binary_m_notime,[],1)./sqrt(sum(~isnan(move_act_binary_m_notime),1)),'k','linewidth',2)
xlabel('Day');
ylabel('Fraction of movement ROIs active during movement');


% Plot total movement and quiescent activity relative to whole population
figure;
subplot(1,2,1);
curr_data = move_act_m_allfrac_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity (relative to all ROIs)');
title('Movement ROI activity');

subplot(1,2,2);
curr_data = move_act_q_allfrac_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity (relative to all ROIs)');
title('Quiescent ROI activity');

% Plot average number of usable movements
n_move = cellfun(@(x) size(x,1),move_act_all);
figure;
errorbar(nanmean(n_move), ...
    nanstd(n_move)./sqrt(sum(~isnan(n_move))),'k','linewidth',2);
ylabel('Number of used movements');
xlabel('Day');

% Plot temporal activity of m/q relative to all activity
figure;
move_act_m_mean_allfrac_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_m_mean_allfrac_day(:,curr_day) = nanmean(vertcat(move_act_m_mean_allfrac{:,curr_day}),1);
end
subplot(1,2,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_m_mean_allfrac_day);
ylabel('\Delta Activity (relative to all ROIs)');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_q_mean_allfrac_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_q_mean_allfrac_day(:,curr_day) = nanmean(vertcat(move_act_q_mean_allfrac{:,curr_day}),1);
end
subplot(1,2,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_q_mean_allfrac_day);
ylabel('\Delta Activity (relative to all ROIs)');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

% Statistics 

% Total activity
days_1 = 1:6;
days_2 = 6:14;

p = anova1(move_act_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_total_mean(:,days_2),[],1));

% Total movement activity
days_1 = 1:6;
days_2 = 6:14;

p = anova1(move_act_m_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_m_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_m_total_mean(:,days_2),[],1));

% Total quiescent activity
days_1 = 1:3;
days_2 = 3:14;

p = anova1(move_act_q_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_q_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_q_total_mean(:,days_2),[],1));

% Binarized movement activity
days_1 = 1:3;
days_2 = 3:7;

p = anova1(move_act_binary_m_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_binary_m_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_binary_m_total_mean(:,days_2),[],1));

% Binarized quiescent activity
days_1 = 1:14;
days_2 = 3:7;

p = anova1(move_act_binary_q_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_binary_q_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_binary_q_total_mean(:,days_2),[],1));


% Amplitude movement activity
days_1 = 3:7;

p = anova1(move_act_amplitude_m_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_amplitude_m_total_mean(:,days_1),[],1));


% Amplitude quiescent activity
days_1 = 1:14;

p = anova1(move_act_amplitude_q_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_amplitude_q_total_mean(:,days_1),[],1));


% Fraction active ROIs during movement
days_1 = 1:3;
days_2 = 3:7;

p = anova1(move_act_fracrois,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_fracrois(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_fracrois(:,days_2),[],1));

p = anova1(move_act_binary_m_notime,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_binary_m_notime(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_binary_m_notime(:,days_2),[],1));

% Amplitude of active
days_1 = 1:6;
days_2 = 6:14;

p = anova1(move_act_amplitude_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_amplitude_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_amplitude_total_mean(:,days_2),[],1));


% Total M/Q activity relative to all ROIs
days_1 = 1:14;
days_2 = 1:14;

p = anova1(move_act_m_allfrac_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_m_allfrac_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_m_allfrac_total_mean(:,days_2),[],1));

p = anova1(move_act_q_allfrac_total_mean,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(move_act_q_allfrac_total_mean(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(move_act_q_allfrac_total_mean(:,days_2),[],1));


%%  (same as above, stripped down)


% Distribution of activity from normalized onset to offset all movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

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
        min_move_time = 28*1;
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
        forward_frames = 28*1;        
        total_frames = back_frames+forward_frames+1;
        
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        move_act_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        % To zero out all non-classified activity
%         curr_class = (classified_rois(curr_animal).movement(:,curr_day) == 1) | ...
%             (classified_rois(curr_animal).quiescent(:,curr_day) == 1);
%         move_act_timed(:,:,~curr_class) = 0;
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end


% Set all zeros to NaN (to isolate amplitude component)
move_act_amplitude_all = move_act_all;
for i = 1:numel(move_act_amplitude_all)
    % Use this for amplitude of all discrete events
    [m,n,p] = size(move_act_amplitude_all{i});
    event_onsets = [zeros(m,1,p),diff(move_act_amplitude_all{i} > 0,[],2) == 1];
    move_act_amplitude_all{i}(~event_onsets) = NaN;
    
    % Use this for amplitude of all active frames
    %move_act_amplitude_all{i}(move_act_amplitude_all{i} == 0) = NaN;
end

move_act_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_all,'uni',false);
% To use only classified cells on each day to build total 
% move_act_mean = cell(size(move_act_mean));
% for curr_animal = 1:length(data)
%     for curr_day = 1:length(data(curr_animal).im)
%         use_cells = (classified_rois(curr_animal).movement(:,curr_day) == 1) | ...
%             (classified_rois(curr_animal).quiescent(:,curr_day) == 1);
%         curr_data = move_act_all{curr_animal,curr_day};
%         curr_data(:,:,~use_cells) = 0;
%         move_act_mean{curr_animal,curr_day} = ...
%             nanmean(nanmean(curr_data,1),3);
%     end
% end
move_act_binary_mean = cellfun(@(x) nanmean(nanmean(x > 0,1),3),move_act_all,'uni',false);
move_act_amplitude_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_amplitude_all,'uni',false);

move_act_m_mean = cell(size(move_act_mean));
move_act_m_mean_allfrac = cell(size(move_act_mean));
move_act_binary_m_mean = cell(size(move_act_mean));
move_act_amplitude_m_mean = cell(size(move_act_mean));

move_act_q_mean = cell(size(move_act_mean));
move_act_q_mean_allfrac = cell(size(move_act_mean));
move_act_binary_q_mean = cell(size(move_act_mean));
move_act_amplitude_q_mean = cell(size(move_act_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).movement(:,curr_day) == 1)) = 0;
        move_act_m_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
                
        move_act_binary_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        move_act_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        curr_data = move_act_all{curr_animal,curr_day};
        curr_data(:,:,~(classified_rois(curr_animal).quiescent(:,curr_day) == 1)) = 0;
        move_act_q_mean_allfrac{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_data,1),3);
        
        move_act_binary_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
    end
end

% Fill in empty data with nans
empty_data = cellfun(@isempty,move_act_all);
move_act_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean_allfrac(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean_allfrac(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);

% Baseline normalize by subtraction
baseline_frames = [-20:-10] + back_frames + 1;

move_act_baseline = cellfun(@(x) nanmean(x(baseline_frames)),move_act_mean);
move_act_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_mean,'uni',false);

move_act_m_baseline_mean = cellfun(@(x) nanmean(x(baseline_frames)),move_act_m_mean);
move_act_m_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean,'uni',false);
move_act_m_mean_allfrac_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean_allfrac,'uni',false);

move_act_q_baseline_mean = cellfun(@(x) nanmean(x(baseline_frames)),move_act_q_mean);
move_act_q_mean_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean,'uni',false);
move_act_q_mean_allfrac_baselinesub = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean_allfrac,'uni',false);

% Get (mean) total activity
move_act_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean);
move_act_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean_baselinesub);

move_act_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean);
move_act_m_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean_baselinesub);
move_act_m_allfrac_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean_allfrac_baselinesub);

move_act_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean);
move_act_q_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean_baselinesub);
move_act_q_allfrac_total_mean_baselinesub = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean_allfrac_baselinesub);

% Plot average total  activity
figure; 

subplot(2,2,1);
curr_data = move_act_m_total_mean_baselinesub;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Movement ROI activity (during movement)');

subplot(2,2,2);
curr_data = move_act_q_total_mean_baselinesub;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Quiescent ROI activity (during movement)');

subplot(2,2,3);
curr_data = move_act_m_baseline_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Movement ROI activity (before movement)');

subplot(2,2,4);
curr_data = move_act_q_baseline_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Quiescent ROI activity (before movement)');

% Plot average temporal activity
figure; 

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_m_mean_baselinesub{:,curr_day}),1);
end
subplot(2,2,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_q_mean_baselinesub{:,curr_day}),1);
end
subplot(2,2,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_m_mean{:,curr_day}),1);
end
subplot(2,2,3); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('Activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_q_mean{:,curr_day}),1);
end
subplot(2,2,4); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('Activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');


% Plot changes in total activity
figure; 

subplot(3,2,1);
curr_data = move_act_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Total activity (during movement)');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_mean{:,curr_day}),1);
end
subplot(3,2,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('Activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

subplot(3,2,3);
curr_data = move_act_total_mean_baselinesub;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity');
title('Total activity (during movement)');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_mean_baselinesub{:,curr_day}),1);
end
subplot(3,2,4); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

subplot(3,2,5);
curr_data = move_act_baseline;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('Activity');
title('Total activity (before movement)');

curr_data_day = nan(total_frames,14);
for curr_day = 1:14
    curr_data_day(:,curr_day) = nanmean(vertcat(move_act_mean_baselinesub{:,curr_day}),1);
end
subplot(3,2,6); hold on;
set(gca,'ColorOrder',jet(14));
plot(curr_data_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');



% Statistics template
curr_data = move_act_total_mean;

use_days = 2:10;

curr_days_repmat = repmat(use_days,length(data),1);
curr_data_days = curr_data(:,use_days);
[r,p] = corrcoef(curr_days_repmat(~isnan(curr_data_days)), ...
    curr_data_days(~isnan(curr_data_days)))

%% Temporal & total activity during cued-rewarded movements

% Distribution of activity from normalized onset to offset all movements
move_act_all = cell(length(data),14);
for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:));
        
        % Get lever movement
        lever_active = analysis(curr_animal).lever(curr_day).lever_move_frames;
       
        % Pull out timed activity
        surrounding_frames = analysis(curr_animal).surrounding_frames;
        
        back_frames = 28*1;
        forward_frames = 28*1;        
        total_frames = back_frames+forward_frames+1;
        
        use_frames = find(surrounding_frames == 0) + [-back_frames:forward_frames];
       
        move_act_timed = analysis(curr_animal).im(curr_day).move_onset_aligned(:, ...
            use_frames,:);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        
    end
    disp(curr_animal);
end


% Set all zeros to NaN (to isolate amplitude component)
move_act_amplitude_all = move_act_all;
for i = 1:numel(move_act_amplitude_all)
   move_act_amplitude_all{i}(move_act_amplitude_all{i} == 0) = NaN; 
end

move_act_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_all,'uni',false);
move_act_binary_mean = cellfun(@(x) nanmean(nanmean(x > 0,1),3),move_act_all,'uni',false);
move_act_amplitude_mean = cellfun(@(x) nanmean(nanmean(x > 0,1),3),move_act_amplitude_all,'uni',false);

move_act_m_mean = cell(size(move_act_mean));
move_act_binary_m_mean = cell(size(move_act_mean));
move_act_amplitude_m_mean = cell(size(move_act_mean));

move_act_q_mean = cell(size(move_act_mean));
move_act_binary_q_mean = cell(size(move_act_mean));
move_act_amplitude_q_mean = cell(size(move_act_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im)
        
        move_act_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
                
        move_act_binary_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_m_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
        
        move_act_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        move_act_binary_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1) > 0,1),3);
        
        move_act_amplitude_q_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_amplitude_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
    end
end

% Fill in empty data with nans
empty_data = cellfun(@isempty,move_act_all);
move_act_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_m_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_binary_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);
move_act_amplitude_q_mean(empty_data) = repmat({nan(1,total_frames)},sum(empty_data(:)),1);

% Baseline normalize by subtraction
baseline_frames = [-20:-10] + back_frames + 1;

move_act_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_mean,'uni',false);
move_act_binary_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_mean,'uni',false);
move_act_amplitude_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_mean,'uni',false);

move_act_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_m_mean,'uni',false);
move_act_binary_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_m_mean,'uni',false);
move_act_amplitude_m_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_m_mean,'uni',false);

move_act_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_q_mean,'uni',false);
move_act_binary_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_binary_q_mean,'uni',false);
move_act_amplitude_q_mean = cellfun(@(x) x-nanmean(x(baseline_frames)),move_act_amplitude_q_mean,'uni',false);

% Get (mean) total activity
move_act_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_mean);
move_act_binary_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_mean);
move_act_amplitude_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_amplitude_mean);
move_act_fracrois = cellfun(@(x) nanmean(nanmean(any(x(:,back_frames+1:end,:),2),3)),move_act_all);

move_act_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_m_mean);
move_act_binary_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_m_mean);
move_act_amplitude_m_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_amplitude_m_mean);

move_act_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_q_mean);
move_act_binary_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_binary_q_mean);
move_act_amplitude_q_total_mean = cellfun(@(x) nanmean(x(back_frames+1:end)),move_act_amplitude_q_mean);


% Plot average temporal activity over days

figure; 

move_act_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_day(:,curr_day) = nanmean(vertcat(move_act_mean{:,curr_day}),1);
end
subplot(3,3,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_m_day(:,curr_day) = nanmean(vertcat(move_act_m_mean{:,curr_day}),1);
end
subplot(3,3,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_m_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_q_day(:,curr_day) = nanmean(vertcat(move_act_q_mean{:,curr_day}),1);
end
subplot(3,3,3); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_q_day);
ylabel('\Delta Activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_day(:,curr_day) = nanmean(vertcat(move_act_binary_mean{:,curr_day}),1);
end
subplot(3,3,4); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_m_day(:,curr_day) = nanmean(vertcat(move_act_binary_m_mean{:,curr_day}),1);
end
subplot(3,3,5); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_m_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_binary_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_binary_q_day(:,curr_day) = nanmean(vertcat(move_act_binary_q_mean{:,curr_day}),1);
end
subplot(3,3,6); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_binary_q_day);
ylabel('\Delta Binarized activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_mean{:,curr_day}),1);
end
subplot(3,3,7); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Total activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_m_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_m_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_m_mean{:,curr_day}),1);
end
subplot(3,3,8); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_m_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Movement ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');

move_act_amplitude_q_day = nan(total_frames,14);
for curr_day = 1:14
    move_act_amplitude_q_day(:,curr_day) = nanmean(vertcat(move_act_amplitude_q_mean{:,curr_day}),1);
end
subplot(3,3,9); hold on;
set(gca,'ColorOrder',jet(14));
plot(move_act_amplitude_q_day);
ylabel('\Delta Amplitude activity');
xlabel('Frames');
title('Quiescent ROI activity');
line([back_frames+1,back_frames+1],ylim,'color','k');


% Plot average total activity

figure; 

subplot(3,3,1);
curr_data = move_act_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity')
title('Total activity');

subplot(3,3,2);
curr_data = move_act_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity')
title('Movement ROI activity');

subplot(3,3,3);
curr_data = move_act_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Activity')
title('Quiescent ROI activity');

subplot(3,3,4);
curr_data = move_act_binary_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity')
title('Total activity');

subplot(3,3,5);
curr_data = move_act_binary_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity')
title('Movement ROI activity');

subplot(3,3,6);
curr_data = move_act_binary_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Binarized activity')
title('Quiescent ROI activity');

subplot(3,3,7);
curr_data = move_act_amplitude_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity')
title('Total activity');

subplot(3,3,8);
curr_data = move_act_amplitude_m_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity')
title('Movement ROI activity');

subplot(3,3,9);
curr_data = move_act_amplitude_q_total_mean;
errorbar(nanmean(curr_data), ...
    nanstd(curr_data)./sqrt(sum(~isnan(curr_data))),'k','linewidth',2)
xlabel('Day');
ylabel('\Delta Amplitude activity')
title('Quiescent ROI activity');

% Plot fraction of active ROIs per movement
figure; 
errorbar(nanmean(move_act_fracrois,1), ...
    nanstd(move_act_fracrois,[],1)./sqrt(sum(~isnan(move_act_fracrois),1)),'k','linewidth',2)
xlabel('Day');
ylabel('Fraction of ROIs active during movement');


%% Fraction of non-switching classified ROIs

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);

m_frac = cell2mat(cellfun(@(x) nanmean(x,1),m_pad','uni',false));
q_frac = cell2mat(cellfun(@(x) nanmean(x,1),q_pad','uni',false));

m2q_pad = cell(length(data),1);
for curr_animal = 1:length(data)
    
    m2q_pad{curr_animal} = ...
        (cumsum(m_pad{curr_animal},2) >= 1).* ...
        q_pad{curr_animal};
    
end

q2m_pad = cell(length(data),1);
for curr_animal = 1:length(data)
    
    q2m_pad{curr_animal} = ...
        (cumsum(q_pad{curr_animal},2) >= 1).* ...
        m_pad{curr_animal};
    
end


m2q_frac = cell2mat(cellfun(@(x) nanmean(x,1),m2q_pad,'uni',false));
q2m_frac = cell2mat(cellfun(@(x) nanmean(x,1),q2m_pad,'uni',false));

figure; hold on
errorbar(nanmean(m_frac),nanstd(m_frac)./sqrt(sum(~isnan(m_frac))),'k','linewidth',2)
errorbar(nanmean(q_frac),nanstd(q_frac)./sqrt(sum(~isnan(q_frac))),'r','linewidth',2)
errorbar(nanmean(m2q_frac),nanstd(m2q_frac)./sqrt(sum(~isnan(m2q_frac))),'b','linewidth',2)
errorbar(nanmean(q2m_frac),nanstd(q2m_frac)./sqrt(sum(~isnan(q2m_frac))),'m','linewidth',2)

ylabel('Fraction of ROIs');
xlabel('Day');
legend({'M','Q','M to Q','Q to M'})

% Plot all ROIs
figure;

[~,sort_idx] = sort(nansum(vertcat(m_pad{:})-vertcat(q_pad{:}),2));
class_full = vertcat(m_pad{:})-vertcat(q_pad{:});

subplot(1,3,1);
imagesc(class_full(sort_idx,:));colormap(gray)
ylabel('Sorted ROIs');
xlabel('Day');
title('All ROIs')

subplot(1,3,2);
imagesc(class_full(any(vertcat(m2q_pad{:}),2),:));colormap(gray)
ylabel('ROIs');
xlabel('Day');
title('All M2Q')

subplot(1,3,3);
imagesc(class_full(any(vertcat(q2m_pad{:}),2),:));colormap(gray)
ylabel('ROIs');
xlabel('Day');
title('All Q2M')


%% Classification by significant activity pre/post movement

% Distribution of activity from normalized onset to offset all movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        min_move_time = 28*1;
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
        forward_frames = 28*1;        
        total_frames = back_frames+forward_frames+1;
        
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        move_act_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Get significantly different average m/q activity
prepost_move_class = cell(length(data),1);

for curr_animal = 1:length(data);
    
    prepost_move_class{curr_animal} = ...
        zeros(length(data(curr_animal).roi_group),length(data(curr_animal).im));
    
    for curr_day = 1:length(data(curr_animal).im);
        
        curr_data = move_act_all{curr_animal,curr_day};
        
        curr_q_act = permute(nanmean(curr_data(:,1:back_frames,:),2),[1,3,2]);
        curr_m_act = permute(nanmean(curr_data(:,back_frames+1:end,:),2),[1,3,2]);
        
        avg_q_act = nanmean(curr_q_act,1);
        avg_m_act = nanmean(curr_m_act,1);
        
        avg_act_diff = avg_m_act - avg_q_act;
        
        % Build distribution by shuffling m/q periods within trial
        
        cat_act = cat(3,curr_q_act,curr_m_act);
        
        n_shuff = 1000;
        avg_act_diff_shuff = nan(n_shuff,size(curr_data,3));
        for curr_shuff = 1:n_shuff            
            % too slow
            %curr_data_shuff = permute(shake(cat_act,3),[3,1,2]);
            
            curr_rand = rand(size(cat_act,1),size(cat_act,2)) > 0.5;
            curr_data_1_shuff = nan(size(cat_act,1),size(cat_act,2));
            curr_data_2_shuff = nan(size(cat_act,1),size(cat_act,2));
            
            curr_data_1_shuff(curr_rand) = curr_q_act(curr_rand);
            curr_data_1_shuff(~curr_rand) = curr_m_act(~curr_rand);
            
            curr_data_2_shuff(~curr_rand) = curr_q_act(~curr_rand);
            curr_data_2_shuff(curr_rand) = curr_m_act(curr_rand);
            
            avg_act_diff_shuff(curr_shuff,:) = nanmean(curr_data_1_shuff,1) - ...
                nanmean(curr_data_2_shuff,1);   
        end
        
        avg_act_rank = tiedrank([avg_act_diff;avg_act_diff_shuff]);
        avg_act_p = avg_act_rank(1,:)./(n_shuff+1);
        
        curr_mq = (avg_act_p > 0.975) - (avg_act_p < 0.025);
        
        prepost_move_class{curr_animal}(:,curr_day) = curr_mq;

     end
    disp(curr_animal);
end


% Plot fraction of classified ROIs
m_frac = cell2mat(cellfun(@(x) nanmean(x == 1,1),prepost_move_class,'uni',false));
q_frac = cell2mat(cellfun(@(x) nanmean(x == -1,1),prepost_move_class,'uni',false));

m_frac_old = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false));
q_frac_old = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false));

figure; hold on;
errorbar(nanmean(m_frac_old,1),nanstd(m_frac_old,[],1)./sqrt(sum(~isnan(m_frac_old))),'k','linewidth',2)
errorbar(nanmean(q_frac_old,1),nanstd(q_frac_old,[],1)./sqrt(sum(~isnan(q_frac_old))),'r','linewidth',2)

errorbar(nanmean(m_frac,1),nanstd(m_frac,[],1)./sqrt(sum(~isnan(m_frac))),'k','linewidth',2)
errorbar(nanmean(q_frac,1),nanstd(q_frac,[],1)./sqrt(sum(~isnan(q_frac))),'r','linewidth',2)


%% Plot individual ROI activity to find examples

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m2u = cellfun(@(m,q) sum(m,2) >= 3 & sum(q,2) == 0,m,q,'uni',false);
m2q = cellfun(@(m,q) sum(m,2) >= 2 & sum(q,2) >= 2,m,q,'uni',false);
q2q = cellfun(@(m,q) sum(m,2) == 0 & sum(q,2) >= 2,m,q,'uni',false);

% Get activity with all movements > 1 sec

% Distribution of activity from normalized onset to offset all movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

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
        min_move_time = 28*1;
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
        back_frames = 28*3;
        forward_frames = 28*3;        
        total_frames = back_frames+forward_frames+1;
        
        move_frame_timed_idx = cell2mat(arrayfun(@(x) use_movement_start_frames(x)-back_frames: ...
            use_movement_start_frames(x)+forward_frames,1:length(use_movement_start_frames),'uni',false)');  
        
        discard_movements = any(move_frame_timed_idx > size(curr_act,2),2) | ...
            any(move_frame_timed_idx < 1,2);
        
        move_frame_timed_idx(discard_movements,:) = [];
        move_frame_timed_idx = reshape(move_frame_timed_idx',[],1);
        
        move_act_timed = permute(reshape(curr_act(:, ...
            move_frame_timed_idx)',total_frames,[],size(curr_act,1)),[2,1,3]);

        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Plot ROIs to pick good ones
h = figure('Position',[12,537,1456,361]);
for curr_animal = 1:length(data)
    curr_plot_rois = find(m2u{curr_animal})';
    for curr_roi = curr_plot_rois

        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        % To plot heatmap
%         curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
%         
%         for curr_day = 1:length(curr_onset)
%             subplot(2,7,curr_day);
%             imagesc(curr_onset{curr_day});colormap(hot);
%             line([back_frames+1,back_frames+1],ylim,'color','w');
%             
%             if classified_rois(curr_animal).movement(curr_roi,curr_day)
%                 set(gca,'XColor','g');set(gca,'YColor','g');
%             elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
%                 set(gca,'XColor','r');set(gca,'YColor','r');
%             else
%                 set(gca,'XColor','k');set(gca,'YColor','k');
%             end
%             
%             caxis([0 3]);
%             set(gca,'XTick',[]);
%             set(gca,'YTick',[]);
%             set(gca,'LineWidth',3);
%         end
        
        % To plot line plot
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            cla;
            curr_mean = nanmean(curr_onset{curr_day},1);
            curr_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_mean,'k','linewidth',2);
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.2 0.5]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line([back_frames+1,back_frames+1],ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        waitforbuttonpress;
        
    end
end

% Plot and save selected ROIs
m2u_plot = { ...
    [], ... % 1
    [135], ... % 2
    [7,134,227,267], ... % 3
    [], ... % 4
    [100], ... % 5
    [], ... % 6
    [24,149,189,239,286,292], ... % 7
    [9,14,24]}; % 8

u2m_plot = { ...
    [120], ... % 1
    [9,142,156], ... % 2
    [34,114,167,179,181,217,232,250,263,265,271], ... % 3
    [40,93], ... % 4
    [5], ... % 5
    [13], ... % 6
    [245], ... % 7
    [45,141,147,167]}; % 8

m2q_plot = { ...
    [], ... % 1
    [189,205,212,215], ... % 2
    [49,61,74,210,221,222], ... % 3
    [3,18,29,47,65,79], ... % 4
    [], ... % 5
    [], ... % 6
    [], ... % 7
    [47]}; % 8

q2q_plot = { ...
    [7,8,11,12,27,28,29,32,34], ... % 1
    [2,8,18,22,24,30,32,36,50,64,65,77,196], ... % 2
    [8,13,26,30,46,50,67,72], ... % 3
    [7,12,15,17,26,84], ... % 4
    [], ... % 5
    [], ... % 6
    [5,20,33,154,196,285], ... % 7
    [16,27,38,54,56,64]}; % 8


example_roi_savedir = ['/usr/local/lab/People/Andy/Corticospinal/corticospinal_paper/figures/example_rois'];

heatmap_fig = figure('Position',[12,201,993,697]);
set(heatmap_fig,'PaperPositionMode','auto');
line_fig = figure('Position',[12,537,1456,361]);
set(line_fig,'PaperPositionMode','auto');

for curr_animal = 1:length(data)
    curr_plot_rois = m2u_plot{curr_animal};
    for curr_roi = curr_plot_rois

        % Plot as heatmap
        figure(heatmap_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            imagesc(curr_onset{curr_day});colormap(hot);
            line([back_frames+1,back_frames+1],ylim,'color','w');
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            caxis([0 3]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'M2U' filesep 'M2U_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_heatmap'];
        
        saveas(heatmap_fig,curr_savefile,'png');
        
        % Plot as line
        figure(line_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            cla;
            curr_mean = nanmean(curr_onset{curr_day},1);
            curr_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_mean,'k','linewidth',2);
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.2 1.5]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line([back_frames+1,back_frames+1],ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'M2U' filesep 'M2U_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_line'];
        
        saveas(line_fig,curr_savefile,'png');
        
    end
end

for curr_animal = 1:length(data)
    curr_plot_rois = u2m_plot{curr_animal};
    for curr_roi = curr_plot_rois
        
        % Plot as heatmap
        figure(heatmap_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            imagesc(curr_onset{curr_day});colormap(hot);
            line([back_frames+1,back_frames+1],ylim,'color','w');
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            caxis([0 3]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'U2M' filesep 'U2M_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_heatmap'];
        
        saveas(heatmap_fig,curr_savefile,'png');
        
        % Plot as line
        figure(line_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            cla;
            curr_mean = nanmean(curr_onset{curr_day},1);
            curr_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_mean,'k','linewidth',2);
                        
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.2 1.5]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line([back_frames+1,back_frames+1],ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'U2M' filesep 'U2M_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_line'];
        
        saveas(line_fig,curr_savefile,'png');
        
    end
end

for curr_animal = 1:length(data)
    curr_plot_rois = m2q_plot{curr_animal};
    for curr_roi = curr_plot_rois

        % Plot as heatmap
        figure(heatmap_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            imagesc(curr_onset{curr_day});colormap(hot);
            line([back_frames+1,back_frames+1],ylim,'color','w');
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            caxis([0 3]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'M2Q' filesep 'M2Q_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_heatmap'];
        
        saveas(heatmap_fig,curr_savefile,'png');
        
        % Plot as line
        figure(line_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            cla;
            curr_mean = nanmean(curr_onset{curr_day},1);
            curr_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_mean,'k','linewidth',2);
                        
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.2 1.5]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line([back_frames+1,back_frames+1],ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'M2Q' filesep 'M2Q_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_line'];
        
        saveas(line_fig,curr_savefile,'png');
        
    end
end


for curr_animal = 1:length(data)
    curr_plot_rois = q2q_plot{curr_animal};
    for curr_roi = curr_plot_rois

        % Plot as heatmap
        figure(heatmap_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            imagesc(curr_onset{curr_day});colormap(hot);
            line([back_frames+1,back_frames+1],ylim,'color','w');
            
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            caxis([0 3]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'Q2Q' filesep 'Q2Q_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_heatmap'];
        
        saveas(heatmap_fig,curr_savefile,'png');
        
        % Plot as line
        figure(line_fig);
        set(gcf,'Name',['Animal ' num2str(curr_animal) ', ROI grp ' num2str(curr_roi)]);
        
        curr_onset = cellfun(@(x) x(:,:,curr_roi),move_act_all(curr_animal,:),'uni',false);
        
        for curr_day = 1:length(curr_onset)
            subplot(2,7,curr_day);
            cla;
            curr_mean = nanmean(curr_onset{curr_day},1);
            curr_sem = nanstd(curr_onset{curr_day},[],1)./sqrt(sum(~isnan(curr_onset{curr_day})));
            jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5],'none');
            hold on;
            plot(curr_mean,'k','linewidth',2);
                        
            if classified_rois(curr_animal).movement(curr_roi,curr_day)
                set(gca,'XColor','g');set(gca,'YColor','g');
            elseif classified_rois(curr_animal).quiescent(curr_roi,curr_day)
                set(gca,'XColor','r');set(gca,'YColor','r');
            else
                set(gca,'XColor','k');set(gca,'YColor','k');
            end
            
            ylim([-0.1 1]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            set(gca,'LineWidth',3);
            line([back_frames+1,back_frames+1],ylim,'color','k');

            title(['A ' num2str(curr_animal) ', R ' num2str(curr_roi)]);
        end
        
        curr_savefile = [example_roi_savedir filesep 'Q2Q' filesep 'Q2Q_A' ...
            num2str(curr_animal) '_R' num2str(curr_roi) '_activity_line'];
        
        saveas(line_fig,curr_savefile,'png');
        
    end
end

close(heatmap_fig);
close(line_fig);

%% Plot chosen example ROIs, only line plot, stitched onset/offset

roi_plot = { ...
    [120], ... % 1
    [65,135], ... % 2
    [13,222], ... % 3
    [93,26], ... % 4
    [], ... % 5
    [], ... % 6
    [], ... % 7
    []};    % 8

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


example_roi_savedir = ['/usr/local/lab/People/Andy/Corticospinal/corticospinal_paper/figures/example_rois/chosen'];

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


%% Standard deviation of activity onset for all (selected) movements

% Distribution of activity from normalized onset to offset all movements
act_onset_std = cell(length(data),14);

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
        back_frames = 28*0;
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

        % Get activity onsets
        [onset_act,onset_frame] = max(move_act_timed > 0,[],2) ;
        onset_frame(onset_act == 0) = NaN;
        
        % Get std of activity onset time in each ROI
        curr_act_onset_std = permute(nanstd(onset_frame,[],1),[3,1,2]);   
        
        % Ignore ROIs with less than minimum active trials
        min_trials = 10;
        curr_act_onset_std(permute(sum(onset_act) < min_trials,[3,1,2])) = NaN;
        
        act_onset_std{curr_animal,curr_day} = curr_act_onset_std;
        
    end
    disp(curr_animal);
end

% Plot changes in activity onset std over time
act_onset_std_all = cellfun(@nanmean,act_onset_std);

act_onset_std_m = nan(size(act_onset_std));
act_onset_std_q = nan(size(act_onset_std));
for curr_animal = 1:length(data);
    curr_std = horzcat(act_onset_std{curr_animal,1:14});
    
    curr_std_m = curr_std;
    curr_std_m(~(classified_rois(curr_animal).movement(:,1:size(curr_std,2)) == 1)) = NaN;
    act_onset_std_m(curr_animal,1:size(curr_std,2)) = nanmean(curr_std_m,1);
    
    curr_std_q = curr_std;
    curr_std_q(~(classified_rois(curr_animal).quiescent(:,1:size(curr_std,2)) == 1)) = NaN;
    act_onset_std_q(curr_animal,1:size(curr_std,2)) = nanmean(curr_std_q,1);
    
end

figure; hold on;
errorbar(nanmean(act_onset_std_all,1),nanstd(act_onset_std_all,[],1)./ ...
    sqrt(sum(~isnan(act_onset_std_all))),'k','linewidth',2)

errorbar(nanmean(act_onset_std_m,1),nanstd(act_onset_std_m,[],1)./ ...
    sqrt(sum(~isnan(act_onset_std_m))),'g','linewidth',2)

errorbar(nanmean(act_onset_std_q,1),nanstd(act_onset_std_q,[],1)./ ...
    sqrt(sum(~isnan(act_onset_std_q))),'r','linewidth',2)

legend({'All ROIs','Movement ROIs','Quiescent ROIs'});
ylabel('Standard deviation (frames)');
xlabel('Day');
title('Activity onset variability');


%% Activity correlation for all (selected) movements

avg_act = cell(length(data),14);

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
        back_frames = 28*0;
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
        
        
        avg_act{curr_animal,curr_day} = permute(nanmean(move_act_timed,1),[3,2,1]);    
        
    end
    
    disp(curr_animal);

end

% Get temporal correlation of movement-related ROIs
long_corr = nan(14,14,length(data));
roi_corr = nan(14,14,length(data));
for curr_animal = 1:length(data)
    
    curr_m = classified_rois(curr_animal).movement == 1;  
    
    curr_data = nan(size(curr_m,1),total_frames,14);
    for curr_day = 1:min(length(data(curr_animal).im),14);                
        curr_data(:,:,curr_day) = avg_act{curr_animal,curr_day};
        curr_data(~curr_m(:,curr_day),:,curr_day) = NaN;   
    end
    
    % To correlated across all movement-related ROIs together
    curr_data_long = reshape(permute(curr_data,[2,1,3]),[],size(curr_data,3));
    curr_data_long_corr = corrcoef(curr_data_long,'rows','pairwise');
    long_corr(:,:,curr_animal) = curr_data_long_corr;
    
    % To average together correlation within movement-related ROIs
    curr_roi_corr = arrayfun(@(x) corrcoef(permute(curr_data(x,:,:),[2,3,1]), ...
        'rows','pairwise'),1:size(curr_data,1),'uni',false);
    curr_roi_corr_mean = nanmean(cat(3,curr_roi_corr{:}),3);
    roi_corr(:,:,curr_animal) = curr_roi_corr_mean;
    
end

figure;
subplot(1,2,1);
long_corr_plot = nanmean(long_corr,3);
long_corr_plot(logical(eye(14))) = NaN;
imagesc(long_corr_plot);colormap(hot);
xlabel('Day');
ylabel('Day');
title('Correlation across movement ROIs');
axis square;

subplot(1,2,2);
roi_corr_plot = nanmean(roi_corr,3);
roi_corr_plot(logical(eye(14))) = NaN;
imagesc(nanmean(roi_corr_plot,3));colormap(hot);
xlabel('Day');
ylabel('Day');
title('Correlation within movement ROIs');
axis square;

%% Distribution of activity for all (selected) movements

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

% Plot average across animals
figure;
% Total activity
subplot(1,2,1); hold on
set(gca,'ColorOrder',jet(14));
plot(nanmean(move_act_total,3),'linewidth',2);
xlabel('Day');
ylabel('Average \Delta activity per ROI');
title('Movement ROIs')
axis square
% Distribution
subplot(1,2,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(nanmean(move_act_distribution,3),'linewidth',2);
xlabel('Day');
ylabel('Fraction of total activity');
title('Movement ROIs')
axis square

% Plot average across animals, group days
figure;
day_group = {[1,2,3,4],[11,12,13,14]};
% Total activity
move_act_total_mean = nanmean(move_act_total,3);
move_act_total_grp = nan(size(move_act_total_mean,1),length(day_group));
for i = 1:length(day_group)
   move_act_total_grp(:,i) = nanmean(move_act_total_mean(:,day_group{i}),2); 
end
subplot(1,2,1); hold on
set(gca,'ColorOrder',copper(length(day_group)));
plot(move_act_total_grp,'linewidth',2);
xlabel('Day');
ylabel('Average \Delta activity per ROI');
title('Movement ROIs')
legend(cellfun(@(x) num2str(x),day_group,'uni',false));
axis square
% Distribution
move_act_distribution_mean = nanmean(move_act_distribution,3);
move_act_distribution_grp = nan(size(move_act_distribution_mean,1),length(day_group));
for i = 1:length(day_group)
   move_act_distribution_grp(:,i) = nanmean(move_act_distribution_mean(:,day_group{i}),2); 
end
subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(length(day_group)));
plot(move_act_distribution_grp,'linewidth',2);
xlabel('Day');
ylabel('Fraction of total activity');
title('Movement ROIs')
legend(cellfun(@(x) num2str(x),day_group,'uni',false));
axis square


% Statistics
move_act_total_grp = nan(size(move_act_distribution,1),length(day_group),length(data));
for i = 1:length(day_group)
   move_act_total_grp(:,i,:) = nanmean(move_act_distribution(:,day_group{i},:),2); 
end

move_act_total_sec1 = permute(nanmean(move_act_total_grp(1:28,:,:,1)),[2,3,1]);
move_act_total_sec2 = permute(nanmean(move_act_total_grp(29:end,:,:,1)),[2,3,1]);

d = repmat(transpose(1:length(day_group)),1,length(data));
[r,p] = corrcoef(d(~isnan(move_act_total_sec1)),move_act_total_sec1(~isnan(move_act_total_sec1)));
[r,p] = corrcoef(d(~isnan(move_act_total_sec2)),move_act_total_sec2(~isnan(move_act_total_sec2)));


%% Activity and fraction of ROIs with activity at different time bins

binned_frac_m = cell(length(data),14);
binned_act_m = cell(length(data),14);
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
        min_move_time = 28*1;
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
        forward_frames = 28*1;        
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
        
        n_time_bins = 3;
        time_bin_edges = linspace(1,1+forward_frames,n_time_bins+1);
        time_bin_edges_inf = time_bin_edges;
        time_bin_edges_inf(end) = Inf;
        [time_bin_n,time_bin] = histc(avg_act_m_maxtime,time_bin_edges_inf);
        
        binned_act_mean = nan(n_time_bins,1);
        for curr_bin = 1:n_time_bins
            
            curr_bin_act = ...
                nanmean(nanmean(avg_act_m(time_bin == curr_bin, ...
                time_bin_edges(curr_bin):time_bin_edges(curr_bin+1)),2));
            
            curr_bin_act_baselinesub = curr_bin_act - ...
                nanmean(avg_act_m_baseline(time_bin == curr_bin));
            
            binned_act_mean(curr_bin) = curr_bin_act_baselinesub;
                
        end
        
        binned_frac_m{curr_animal,curr_day} = time_bin_n(1:end-1)./sum(curr_m);
        binned_act_m{curr_animal,curr_day} = binned_act_mean;        
        
        curr_pref_act_time = nan(size(curr_m));
        curr_pref_act_time(curr_m) = avg_act_m_maxtime;
        pref_act_time{curr_animal,curr_day} = curr_pref_act_time;
        
    end
    
    disp(curr_animal);

end

% Get mean values across animals

binned_frac_m_day = nan(n_time_bins,14);
binned_act_m_day = nan(n_time_bins,14);
for curr_day = 1:14
    binned_frac_m_day(:,curr_day) = nanmean(horzcat(binned_frac_m{:,curr_day}),2);
    binned_act_m_day(:,curr_day) = nanmean(horzcat(binned_act_m{:,curr_day}),2);  
end

figure;

subplot(1,2,1); hold on;
set(gca,'ColorOrder',copper(n_time_bins));
plot(binned_frac_m_day','linewidth',2);
xlabel('Day');
ylabel('Fraction of movement ROIs in bin');
legend(cellfun(@(x) ['Bin ' num2str(x)],num2cell(1:n_time_bins),'uni',false));

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(n_time_bins));
plot(binned_act_m_day','linewidth',2);
xlabel('Day');
ylabel('\Delta Activity of movement ROIs in bin');
legend(cellfun(@(x) ['Bin ' num2str(x)],num2cell(1:n_time_bins),'uni',false));

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


%% Activity/movement correlation of ALL movements (quick - downsampled)

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Get correlation between movement and activity
stages = {[1:4],[11:14]};

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(data);
    
    curr_act = cellfun(@(x) reshape(permute(x,[2,3,1]), ...
        [],size(x,1)),move_act_all(curr_animal,:),'uni',false);
    curr_move = cellfun(@(x) x',move_position_all(curr_animal,:),'uni',false);
    
    n_trials = cellfun(@(x) size(x,2),curr_move);
    
    movement_corr = mat2cell(corrcoef(horzcat(curr_move{:})),n_trials,n_trials);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    activity_corr = mat2cell(corrcoef(horzcat(curr_act{:})),n_trials,n_trials);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
        
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
       
    % Use sessions within stage that exist 
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages);
        for curr_stage_2 = 1:length(stages);
            
            sessions_1 = stages{curr_stage_1};
            sessions_2 = stages{curr_stage_2};
                  
            curr_num = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'numel');
            curr_mean = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'mean');
            curr_sem = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'sem');
            
            bins_used = unique(vertcat(stages_bin_idx{sessions_1,sessions_2}));
            
            % Cutoff of too small number of trial combinations
            min_trial_pairs = 50;
            below_min_bins = curr_num < min_trial_pairs;
            curr_mean(below_min_bins) = [];
            curr_sem(below_min_bins) = [];
            bins_used(below_min_bins) = [];
            
            bins_unused = setdiff(1:length(corr_bin),bins_used);
            
            curr_stage_mean{curr_stage_1,curr_stage_2}(bins_used) = curr_mean;
            curr_stage_sem{curr_stage_1,curr_stage_2}(bins_used) = curr_sem;
            
            % fill in unused values with nans
            curr_stage_mean{curr_stage_1,curr_stage_2}(bins_unused) = NaN;
            curr_stage_sem{curr_stage_1,curr_stage_2}(bins_unused) = NaN;
            
        end
    end
    
    activity_corr_stages_mean(:,:,curr_animal) = curr_stage_mean;
    activity_corr_stages_sem(:,:,curr_animal) = curr_stage_sem;
    
    disp(curr_animal);
end

% Combine animals
activity_corr_stages_mean_animals = cell(length(stages));
activity_corr_stages_sem_animals = cell(length(stages));
for stage_1 = 1:length(stages)
    for stage_2 = 1:length(stages)
   
        curr_mean = vertcat(activity_corr_stages_mean{stage_1,stage_2,:});
        activity_corr_stages_mean_animals{stage_1,stage_2} = nanmean(curr_mean,1);
        
%         % Average SEM across animals (this doesn't make sense, but it's
%         what I did last time...)
%         curr_sem = vertcat(activity_corr_stages_sem{stage_1,stage_2,:});
%         activity_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem,1);
        
        % SEM of mean across animals
        activity_corr_stages_sem_animals{stage_1,stage_2} = ...
            nanstd(curr_mean)./sqrt(sum(~isnan(curr_mean),1));
        
    end
end

figure; hold on;
col = repmat(linspace(0,0.7,length(stages))',1,3);
for i = 1:length(stages)
    errorbar(corr_bin_plot,activity_corr_stages_mean_animals{i,i}, ...
        activity_corr_stages_sem_animals{1,1},'color',col(i,:),'linewidth',2);
end
errorbar(corr_bin_plot,activity_corr_stages_mean_animals{1,end}, ...
        activity_corr_stages_sem_animals{1,end},'color','b','linewidth',2);
ylabel('Pairwise activity correlation');
xlabel('Pairwise movement correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],stages,'uni',false), ...
    ['Across sessions ' num2str(stages{1}) ', ' num2str(stages{end})]]);


%% Activity/movement correlation of ALL movements (full sample movements)

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
    epoch_combine_sem(:,i) = nanstd(horzcat(actcorr_bin_mean_pad{:,i}),[],2)./ ...
        sqrt(sum(~isnan(horzcat(actcorr_bin_mean_pad{:,i})),2));
    
    % Used to use average SEM - doesn't really make sense
    %epoch_combine_sem(:,i) = nanmean(horzcat(actcorr_bin_sem_pad{:,i}),2);
end

figure; hold on;
errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);
legend({'Within early','Within late','Across early/late'});
xlabel('Pairwise movement correlation');
ylabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);



%% PCA representation of distance/variance of activity by movement corr

% STARTED THIS BUT WORK IN PROGRESS

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Get correlation between movement and activity
stages = {[1:4],[11:14]};


for curr_animal = 1:length(data);
    
    curr_act = cellfun(@(x) reshape(permute(x,[2,3,1]), ...
        [],size(x,1)),move_act_all(curr_animal,:),'uni',false);
    curr_move = cellfun(@(x) x',move_position_all(curr_animal,:),'uni',false);
    
    n_trials = cellfun(@(x) size(x,2),curr_move);
    
    movement_corr = mat2cell(corrcoef(horzcat(curr_move{:})),n_trials,n_trials);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    [coeff,scores,latents] = princomp(horzcat(curr_act{:})');
    
    
    activity_corr = mat2cell(corrcoef(horzcat(curr_act{:})),n_trials,n_trials);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
        
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
       
    % Use sessions within stage that exist 
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages);
        for curr_stage_2 = 1:length(stages);
            
            sessions_1 = stages{curr_stage_1};
            sessions_2 = stages{curr_stage_2};
                  
            curr_num = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'numel');
            curr_mean = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'mean');
            curr_sem = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'sem');
            
            bins_used = unique(vertcat(stages_bin_idx{sessions_1,sessions_2}));
            
            % Cutoff of too small number of trial combinations
            min_trial_pairs = 50;
            below_min_bins = curr_num < min_trial_pairs;
            curr_mean(below_min_bins) = [];
            curr_sem(below_min_bins) = [];
            bins_used(below_min_bins) = [];
            
            bins_unused = setdiff(1:length(corr_bin),bins_used);
            
            curr_stage_mean{curr_stage_1,curr_stage_2}(bins_used) = curr_mean;
            curr_stage_sem{curr_stage_1,curr_stage_2}(bins_used) = curr_sem;
            
            % fill in unused values with nans
            curr_stage_mean{curr_stage_1,curr_stage_2}(bins_unused) = NaN;
            curr_stage_sem{curr_stage_1,curr_stage_2}(bins_unused) = NaN;
            
        end
    end
    

    
    disp(curr_animal);
end

 
%% Timing of m/q to diverge from baseline around movement

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

move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_all,'uni',false);
move_act_m_onset_mean = cell(size(move_act_onset_mean));
move_act_q_onset_mean = cell(size(move_act_onset_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        baseline_frames = 1:15;
        baseline_m = nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,baseline_frames, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),2),3);           
        
    end
end

move_act_m_onset_mean_day = nan(min(length(data(curr_animal).im),14),total_frames);
move_act_q_onset_mean_day = nan(min(length(data(curr_animal).im),14),total_frames);
for curr_day = 1:min(length(data(curr_animal).im),14)
    
    curr_m_act = vertcat(move_act_m_onset_mean{:,curr_day});
    move_act_m_onset_mean_day(curr_day,:) = nanmean(curr_m_act,1);
       
    curr_q_act = vertcat(move_act_q_onset_mean{:,curr_day});
    move_act_q_onset_mean_day(curr_day,:) = nanmean(curr_q_act,1);
    
end

baseline_frames = 10:15;

m_baseline_mean = nanmean(move_act_m_onset_mean_day(:,baseline_frames),2);
m_baseline_std = nanstd(move_act_m_onset_mean_day(:,baseline_frames),[],2);
[~,m_over_idx] = max(move_act_m_onset_mean_day(:,baseline_frames(end)+1:end) > ...
    repmat(m_baseline_mean+2*m_baseline_std,1,total_frames-baseline_frames(end)),[],2);

framerate = 29;
m_over_idx_time = (nonmove_frames - (m_over_idx + baseline_frames(end)))./framerate;


q_baseline_mean = nanmean(move_act_q_onset_mean_day(:,baseline_frames),2);
q_baseline_std = nanstd(move_act_q_onset_mean_day(:,baseline_frames),[],2);
[~,q_under_idx] = max(move_act_q_onset_mean_day(:,baseline_frames(end)+1:end) < ...
    repmat(q_baseline_mean-2*q_baseline_std,1,total_frames-baseline_frames(end)),[],2);

framerate = 29;
q_under_idx_time = (nonmove_frames - (q_under_idx + baseline_frames(end)))./framerate;

figure;

subplot(2,2,1); hold on; set(gca,'ColorOrder',copper(14));
move_act_m_onset_mean_day_norm = ...
    bsxfun(@rdivide,move_act_m_onset_mean_day,max(move_act_m_onset_mean_day,[],2)./2)-1;
plot(move_act_m_onset_mean_day_norm');
axis off
line([nonmove_frames,nonmove_frames],ylim);
title('Average mROI activity across days');

subplot(2,2,2); hold on; set(gca,'ColorOrder',copper(14));
move_act_q_onset_mean_day_norm = ...
    bsxfun(@rdivide,move_act_q_onset_mean_day,max(move_act_q_onset_mean_day,[],2)./2)-1;
plot(move_act_q_onset_mean_day_norm');
axis off
line([nonmove_frames,nonmove_frames],ylim);
title('Average qROI activity across days');

subplot(2,2,3);
plot(m_over_idx_time,'k','linewidth',2);
xlabel('Day');
ylabel('Time m above baseline');

subplot(2,2,4);
plot(q_under_idx_time,'k','linewidth',2);
xlabel('Day');
ylabel('Time q below baseline');


%% Template analysis for specific stabilization (TK requested) (quick - downsampled movements)

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Make "templates" based on average of end days
template_days = [13,14];
template_days_position = arrayfun(@(x) vertcat(move_position_all{x, ...
    template_days}),1:length(data),'uni',false);
move_template = cellfun(@(x) nanmean(x,1),template_days_position,'uni',false);

template_corr = cell(size(move_position_all));
for curr_animal = 1:length(data)
    for curr_session = 1:14
        
        if isempty(move_template{curr_animal}) || ...
                isempty(move_position_all{curr_animal,curr_session})
            continue
        end
        
        template_corr{curr_animal,curr_session} = ...
            1-pdist2(move_template{curr_animal}, ...
            move_position_all{curr_animal,curr_session},'correlation');     
    end
end

% Bin trials according to correlation with template
n_bins = 5;
template_corr_edges = linspace(-1,1,n_bins+1);
bin_centers = template_corr_edges(1:end-1) + diff(template_corr_edges)/2;
template_corr_edges(end) = inf;
[~,template_corr_bin] = cellfun(@(x) ...
    histc(x,template_corr_edges),template_corr,'uni',false);

binned_move_corr = repmat({nan(n_bins,1)},size(move_position_all,1),size(move_position_all,2));
binned_act_corr = repmat({nan(n_bins,1)},size(move_position_all,1),size(move_position_all,2));
for curr_animal = 1:length(data)
    for curr_session = 1:length(data(curr_animal).im);
        for curr_bin = 1:n_bins
            
            use_movements = template_corr_bin{curr_animal,curr_session} == curr_bin;
            
            % Get correlation of movements within bins
            binned_move_corr{curr_animal,curr_session}(curr_bin) = nanmean(AP_itril(corrcoef( ...
                move_position_all{curr_animal,curr_session}(use_movements,:)'),-1));
            
            curr_act = reshape(permute(move_act_all{curr_animal,curr_session}( ...
                use_movements,:,:),[2,3,1]),[],sum(use_movements));
            binned_act_corr{curr_animal,curr_session}(curr_bin) = nanmean(AP_itril(corrcoef( ...
                curr_act,'rows','pairwise'),-1));
            
            
        end
    end
end

% Make sure that movement correlations within bins don't change
binned_move_corr_mean = nan(n_bins,14);
binned_act_corr_mean = nan(n_bins,14);

for curr_session = 1:14
    binned_move_corr_mean(:,curr_session) = ...
        nanmean(horzcat(binned_move_corr{:,curr_session}),2);   
    binned_act_corr_mean(:,curr_session) = ...
        nanmean(horzcat(binned_act_corr{:,curr_session}),2);   
end

figure;

subplot(2,1,1); hold on;
set(gca,'ColorOrder',copper(n_bins))
plot(binned_move_corr_mean')
ylabel('Movement correlation');
xlabel('Day');

subplot(2,1,2); hold on;
set(gca,'ColorOrder',copper(n_bins))
plot(binned_act_corr_mean')
ylabel('Activity correlation');
xlabel('Day');

legend(cellfun(@(x) ['Corr w/ template: ' num2str(x)],num2cell(bin_centers),'uni',false));

% THIS IS JUST TEMPORARY TO CHECK ALL ACTIVITY CORRELATIONS
binned_act_corr = cell(size(move_position_all,1),size(move_position_all,2));
for curr_animal = 1:length(data)
    for curr_session = 1:length(data(curr_animal).im);
            
            curr_act = reshape(permute(move_act_all{curr_animal,curr_session}( ...
                :,:,:),[2,3,1]),[],size(move_act_all{curr_animal,curr_session},1));
            binned_act_corr{curr_animal,curr_session} = nanmean(AP_itril(corrcoef( ...
                curr_act,'rows','pairwise'),-1));
            
            
    end
end
a = nan(1,14);
for i = 1:14
a(i) = nanmean(horzcat(binned_act_corr{:,i}));
end
figure;plot(a,'k')

%% Activity/movement correlation of ALL movements (quick - downsampled) across sessions

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Get correlation between movement and activity by days separated
stages = {[0:4],[5:9],[10:13]};
%stages = num2cell(0:13);

activity_corr_stages_mean = cell(length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(analysis));
for curr_animal = 1:length(data);
    
    curr_act = cellfun(@(x) reshape(permute(x,[2,3,1]), ...
        [],size(x,1)),move_act_all(curr_animal,:),'uni',false);
    curr_move = cellfun(@(x) x',move_position_all(curr_animal,:),'uni',false);
    
    n_trials = cellfun(@(x) size(x,2),curr_move);
    
    movement_corr = mat2cell(corrcoef(horzcat(curr_move{:})),n_trials,n_trials);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    activity_corr = mat2cell(corrcoef(horzcat(curr_act{:})),n_trials,n_trials);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
        
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,22);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
       
    % Use sessions within stage that exist 
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages),1);
    curr_stage_sem = cell(length(stages),1);
    for curr_stage = 1:length(stages);
        
        curr_sessions = false(size(activity_corr));
        curr_session_idx = reshape(1:numel(curr_sessions),length(curr_sessions),[]);
        for i = 1:length(stages{curr_stage})
            curr_sessions(diag(curr_session_idx,stages{curr_stage}(i))) = true;
        end
        
        curr_num = grpstats(vertcat(activity_corr{curr_sessions}), ...
            vertcat(stages_bin_idx{curr_sessions}),'numel');
        curr_mean = grpstats(vertcat(activity_corr{curr_sessions}), ...
            vertcat(stages_bin_idx{curr_sessions}),'mean');
        curr_sem = grpstats(vertcat(activity_corr{curr_sessions}), ...
            vertcat(stages_bin_idx{curr_sessions}),'sem');
        
        bins_used = unique(vertcat(stages_bin_idx{curr_sessions}));
        
        % Cutoff of too small number of trial combinations
        min_trial_pairs = 50;
        below_min_bins = curr_num < min_trial_pairs;
        curr_mean(below_min_bins) = [];
        curr_sem(below_min_bins) = [];
        bins_used(below_min_bins) = [];
        
        bins_unused = setdiff(1:length(corr_bin),bins_used);
        
        curr_stage_mean{curr_stage}(bins_used) = curr_mean;
        curr_stage_sem{curr_stage}(bins_used) = curr_sem;
        
        % fill in unused values with nans
        curr_stage_mean{curr_stage}(bins_unused) = NaN;
        curr_stage_sem{curr_stage}(bins_unused) = NaN;
        
    end
    
    activity_corr_stages_mean(:,curr_animal) = curr_stage_mean;
    activity_corr_stages_sem(:,curr_animal) = curr_stage_sem;
    
    disp(curr_animal);
end

% Combine animals
activity_corr_stages_mean_animals = cell(length(stages),1);
activity_corr_stages_sem_animals = cell(length(stages),1);
for curr_stage = 1:length(stages)
   
        curr_mean = vertcat(activity_corr_stages_mean{curr_stage,:});
        activity_corr_stages_mean_animals{curr_stage} = nanmean(curr_mean,1);
        
%         % Average SEM across animals (this doesn't make sense, but it's
%         what I did last time...)
%         curr_sem = vertcat(activity_corr_stages_sem{stage_1,stage_2,:});
%         activity_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem,1);
        
        % SEM of mean across animals
        activity_corr_stages_sem_animals{curr_stage} = ...
            nanstd(curr_mean)./sqrt(sum(~isnan(curr_mean),1));
  
end

% Plot 
figure; hold on;
col = repmat(linspace(0,0.7,length(stages))',1,3);
for i = 1:length(stages)
    errorbar(corr_bin_plot,activity_corr_stages_mean_animals{i}, ...
        activity_corr_stages_sem_animals{1,1},'color',col(i,:),'linewidth',2);
end
ylabel('Pairwise activity correlation');
xlabel('Pairwise movement correlation');
legend(cellfun(@(x) [num2str(x) ' days separated'],stages,'uni',false));

% Plot minus the 0-bin value
bin_0 = find(corr_bin_plot == 0);
activity_corr_stages_mean_animals_0sub = cellfun(@(x) x - ...
    x(bin_0),activity_corr_stages_mean_animals,'uni',false);
figure; hold on;
col = repmat(linspace(0,0.7,length(stages))',1,3);
for i = 1:length(stages)
    plot(corr_bin_plot,activity_corr_stages_mean_animals_0sub{i}, ...
        'color',col(i,:),'linewidth',2);
end
ylabel('Pairwise activity correlation');
xlabel('Pairwise movement correlation');
legend(cellfun(@(x) [num2str(x) ' days separated'],stages,'uni',false));


%% Get activity axis of dissimilar movements: is it consistent across days?
% NOTE: is this any different from just getting axes of variance?

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Get correlation between movement and activity
stages = {[1:4],[11:14]};

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(data);
    
    curr_act = cellfun(@(x) reshape(permute(x,[2,3,1]), ...
        [],size(x,1)),move_act_all(curr_animal,:),'uni',false);
    curr_move = cellfun(@(x) x',move_position_all(curr_animal,:),'uni',false);
    
    n_trials = cellfun(@(x) size(x,2),curr_move);
    
    movement_corr = mat2cell(corrcoef(horzcat(curr_move{:})),n_trials,n_trials);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    activity_corr = mat2cell(corrcoef(horzcat(curr_act{:})),n_trials,n_trials);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
        
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
       
    % Use sessions within stage that exist 
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    
end


%% Pairwise activity correlation given movement correlation within epochs (quick - downsampled movements)

% Get movements and activity for all selected movements
move_act_all = cell(length(data),14);
move_position_all = cell(length(data),14);
move_speed_all = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:length(data(curr_animal).im);
            
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
        back_frames = 28*0;
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
        
        position_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_position_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
        
        speed_timed = permute(reshape( ...
            analysis(curr_animal).lever(curr_day).lever_speed_frames( ...
            move_frame_timed_idx),total_frames,[]),[2,1,3]);
  
        move_act_all{curr_animal,curr_day} = move_act_timed;
        move_position_all{curr_animal,curr_day} = position_timed;
        move_speed_all{curr_animal,curr_day} = speed_timed;
        
    end
    disp(curr_animal);
end

% Make templates with 50% of trials from end days
template_days = [11:14];
template_days_position = arrayfun(@(x) vertcat(move_position_all{x, ...
    template_days}),1:length(data),'uni',false);
move_template = cellfun(@(x) nanmean(x,1),template_days_position,'uni',false);

template_corr = cell(size(move_position_all,1),size(move_position_all,2));
for curr_animal = 1:length(data)
    for curr_session = 1:14
        
        if isempty(move_template{curr_animal}) || ...
                isempty(move_position_all{curr_animal,curr_session})
            continue
        end
        
        template_corr{curr_animal,curr_session} = ...
            1-pdist2(move_template{curr_animal}, ...
            move_position_all{curr_animal,curr_session},'correlation');     
    end
end

% Bin trials according to correlation with template
n_bins = 21;
template_corr_edges = linspace(-1,1,n_bins+1);
bin_centers = template_corr_edges(1:end-1) + diff(template_corr_edges)/2;
template_corr_edges(end) = inf;
[~,template_corr_bin] = cellfun(@(x) ...
    histc(x,template_corr_edges),template_corr,'uni',false);

epochs = {[1:4],[11:14]};

binned_move_corr = repmat({nan(n_bins,1)},size(move_position_all,1),length(epochs));
binned_act_corr = repmat({nan(n_bins,1)},size(move_position_all,1),length(epochs));
for curr_animal = 1:length(data)
    for curr_epoch = 1:length(epochs);
        for curr_bin = 1:n_bins
            
            use_movements = horzcat(template_corr_bin{curr_animal,epochs{curr_epoch}}) == curr_bin;
            
            % Get correlation of movements within bins
            cat_movements = vertcat(move_position_all{curr_animal,epochs{curr_epoch}});
            binned_move_corr{curr_animal,curr_epoch}(curr_bin) = nanmean(AP_itril(corrcoef( ...
                cat_movements(use_movements,:)'),-1));
            
            cat_act = vertcat(move_act_all{curr_animal,epochs{curr_epoch}});
            curr_act = reshape(permute(cat_act( ...
                use_movements,:,:),[2,3,1]),[],sum(use_movements));
            binned_act_corr{curr_animal,curr_epoch}(curr_bin) = nanmean(AP_itril(corrcoef( ...
                curr_act,'rows','pairwise'),-1));
            
            
        end
    end
end

% Make sure that movement correlations within bins don't change
binned_move_corr_mean = nan(n_bins,length(epochs));
binned_act_corr_mean = nan(n_bins,length(epochs));

for curr_epoch = 1:length(epochs)
    binned_move_corr_mean(:,curr_epoch) = ...
        nanmean(horzcat(binned_move_corr{:,curr_epoch}),2);   
    binned_act_corr_mean(:,curr_epoch) = ...
        nanmean(horzcat(binned_act_corr{:,curr_epoch}),2);   
end

figure;

subplot(2,1,1); hold on;
set(gca,'ColorOrder',copper(length(epochs)))
plot(bin_centers,binned_move_corr_mean)
ylabel('Pairwise movement correlation');
xlabel('Movement correlation with template');

subplot(2,1,2); hold on;
set(gca,'ColorOrder',copper(length(epochs)))
plot(bin_centers,binned_act_corr_mean)
ylabel('Pairwise activity correlation');
xlabel('Movement correlation with template');

legend(cellfun(@(x) ['Sessions ' num2str(x)],epochs,'uni',false));



