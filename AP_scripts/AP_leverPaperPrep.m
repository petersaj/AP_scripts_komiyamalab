%% Plot average movement v quiescent activity

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

still_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z+1),1:length(x)-1,'uni',false), ...
    move_stops,move_starts,'uni',false);

movement_activity = cell(8,15);
still_activity = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % get movement activity
        move_time_frames = cellfun(@(x) x, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        still_time_frames = cellfun(@(x) x, ...
            still_epoch_frames{curr_animal,curr_session},'uni',false);

        curr_move_activity = cellfun(@(x) ...
            pyr_activity_all{curr_animal,curr_session}(:,x), ...
            move_time_frames,'uni',false);
        
        curr_still_activity = cellfun(@(x) ...
            pyr_activity_all{curr_animal,curr_session}(:,x), ...
            still_time_frames,'uni',false);

        cat_move_activity = horzcat(curr_move_activity{:});
        cat_still_activity = horzcat(curr_still_activity{:});
        % sometimes have nans, make those zero
                    
        movement_activity{curr_animal,curr_session} = nanmean(cat_move_activity,2);
        still_activity{curr_animal,curr_session} = nanmean(cat_still_activity,2);
        
    end
end

% plot movement and quiescent activity over days
use_animals = [1 2 4 5 6 7 8];
movement_activity_cat = cell(14,1);
still_activity_cat = cell(14,1);
for i = 1:14
   movement_activity_cat{i} = vertcat(movement_activity{use_animals,i}); 
   still_activity_cat{i} = vertcat(still_activity{use_animals,i});     
end
mm = cellfun(@nanmean,movement_activity_cat);
ms = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),movement_activity_cat);
sm = cellfun(@nanmean,still_activity_cat);
ss = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),still_activity_cat);
figure; hold on
errorbar(mm,ms,'k','linewidth',2)
errorbar(sm,ss,'r','linewidth',2)
xlabel('Sessions')
ylabel('Average activity')
legend({'Movement','Quiescence'})

% cat all, bar movement and quiescent activity
use_animals = [1 2 4 5 6 7 8];
movement_activity_cat = vertcat(movement_activity{use_animals,1:14}); 
still_activity_cat = vertcat(still_activity{use_animals,1:14});     
mm = nanmean(movement_activity_cat);
ms = nanstd(movement_activity_cat)/sqrt(sum(~isnan(movement_activity_cat)));
sm = nanmean(still_activity_cat);
ss = nanstd(still_activity_cat)/sqrt(sum(~isnan(still_activity_cat)));
figure; hold on
bar([mm sm],'FaceColor','w','linewidth',2,'BarWidth',0.5);
errorbar([mm sm],[ms ss],'.k','linewidth',2)
xlim([0 3])
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Moving' 'Quiescent'})
ylabel('Average activity')



%% Gad/pyr class (in)variance: plot example classification grid and cells

% plot class grids
use_animals = [1 2 4 5 6 7 8];

figure;
persistence = cell(8,2);
for curr_animal = use_animals   
    
    curr_gad_class = vertcat(gad_class_all{curr_animal,:})';
    day_center = nanmean(curr_gad_class.*meshgrid(1:size(curr_gad_class,2), ...
        1:size(curr_gad_class,1)),2);
    gad_class_cells = find(any(curr_gad_class,2));
    [asdf sort_idx] = sort(day_center(gad_class_cells));
    
    subplot(8,2,curr_animal*2-1);
    imagesc(curr_gad_class(gad_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Gad')
    end
    
    curr_pyr_class = vertcat(pyr_class_all{curr_animal,:})';
    day_center = nanmean(curr_pyr_class.*meshgrid(1:size(curr_pyr_class,2), ...
        1:size(curr_pyr_class,1)),2);
    pyr_class_cells = find(any(curr_pyr_class,2));
    [asdf sort_idx] = sort(day_center(pyr_class_cells));
    
    subplot(8,2,curr_animal*2);
    imagesc(curr_pyr_class(pyr_class_cells(sort_idx),:));colormap(gray);
    
    if curr_animal == use_animals(1)
        title('Pyr')
    end
    
    persistence{curr_animal,1} = nanmean(curr_gad_class(gad_class_cells,:),2);
    persistence{curr_animal,2} = nanmean(curr_pyr_class(pyr_class_cells,:),2);
    
end
avg_persistence = cellfun(@nanmean,persistence);
figure; hold on;
bar(nanmean(avg_persistence),'FaceColor','w','linewidth',2,'BarWidth',0.5);
errorbar(nanmean(avg_persistence), ...
    nanstd(avg_persistence)./sqrt(sum(~isnan(avg_persistence))),'.k','linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Gad' 'Pyr'})
ylabel('Persistence')
xlim([0 3])

figure;hold on
plot(sort(a{2}),linspace(0,1,length(a{2})),'color',[0 0.7 0],'linewidth',2)
plot(sort(a{1}),linspace(0,1,length(a{1})),'color',[0.8 0 0],'linewidth',2)
ylabel('Cumulative fraction of cells')
xlabel('Fraction of days active')
title('Gad cells are less transient')
legend({'Pyr' 'Gad'},'location','se')


% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 60;fs_f = 200;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end

% get rewarded move starts for plotting activity (using another matlab)
curr_animal = 4;
sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
rewarded_move_start = cell(max(sessions),1);
for session = sessions
   
    curr_rewardmoves = move_epoch_frames{curr_animal,session}( ...
        rewarded_movements{curr_animal,session});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves);
    rewarded_move_start{session} = curr_moveonsets;
end



%% Excitatory/inhibitory balance consistency

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity with lever 
pyr_move_active = cell(8,15);
gad_move_active = cell(8,15);
movement_lengths = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds or whole movement
        move_time_framespan = 28*3;
        %move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
        %    move_epoch_frames{curr_animal,curr_session},'uni',false);
        move_time_frames = cellfun(@(x) x, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_rewardmoves_leeway);
        
    end

    curr_animal
    
end


% get pyr:gad ratio for movements

use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);

figure; hold on
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
col = jet(length(day_combine));
for i = 1:length(day_combine)

    pc = horzcat(p{use_animals,day_combine{i}});
    gc = horzcat(g{use_animals,day_combine{i}});
    
    bin_edges = [0:0.02:0.25];
    %bin_edges = [0:0.1:1-0.1 Inf];
    [n bins] = histc(pc,bin_edges);
    gc_mean = grpstats(gc,bins);
    gc_sem = grpstats(gc,bins,'sem');
    use_bins = unique(bins) ~= 0;
    use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
    errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins), ...
        'color',col(i,:),'linewidth',2);
        
end

pc = horzcat(p{use_animals,:});
gc = horzcat(g{use_animals,:});
bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'--k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Rewarded movement based activity')

legend([cellfun(@num2str,day_combine,'uni',false) 'All'])



%% Pyr population old/new bar graph makeup

% plot class grids
pyr_oldfrac = nan(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    curr_pyr_class = vertcat(pyr_class_all{curr_animal,:})';
    curr_pyr_old = ~diff(curr_pyr_class,[],2) & curr_pyr_class(:,2:end);

    pyr_oldfrac(curr_animal,sessions(2:end)) = ...
        sum(curr_pyr_old)./sum(curr_pyr_class(:,2:end));
end

figure; hold on;
mean_oldfrac = nanmean(pyr_oldfrac(use_animals,2:14));
bar(2:14,[mean_oldfrac' 1-mean_oldfrac'],'stack')

use_animals = [1 2 4 5 6 7 8];
errorbar(2:14,nanmean(pyr_oldfrac(use_animals,2:14)), ...
    nanstd(pyr_oldfrac(use_animals,2:14))./ ...
    sqrt(sum(~isnan(pyr_oldfrac(use_animals,2:14)))),'.k','linewidth',2)

ylabel('Fraction of movement cells')
xlabel('Session')
legend({'Old' 'New'})


%% Temporal alignment examples

% plot aligned examples
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
moveonset_aligned_pyr = cell(size(movement_trace_all));
moveonset_aligned_gad = cell(size(movement_trace_all));
for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 60;fs_f = 200;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    curr_pyr_activity = pyr_activity_all{curr_animalday};
    curr_gad_activity = gad_activity_thresh_all{curr_animalday};
    
    % store activity 
    curr_moveonset_aligned_pyr = cellfun(@(x) curr_pyr_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   curr_moveonset_aligned_gad = cellfun(@(x) curr_gad_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
   
    moveonset_aligned_pyr{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_pyr{:}), [3 2 1]);
    moveonset_aligned_gad{curr_animalday} = permute( ...
        cat(3,curr_moveonset_aligned_gad{:}), [3 2 1]);
    
    curr_animalday/max(animaldays)

end


% plot class grid
curr_animal = 5;
figure;imagesc(vertcat(pyr_class_all{curr_animal,:})');colormap(gray)

% plot example mean avg day plot
curr_animal = 2;
curr_session = 11;
a = permute(nanmean(moveonset_aligned_pyr{curr_animal,curr_session}),[3 2 1]);
figure;imagesc(a)

% plot the trials from animal/cell
curr_animal = 2;
curr_pyrcell = 134;
figure; 
sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
for session = sessions
    subplot(4,4,session);colormap(gray);
    imagesc(moveonset_aligned_pyr{curr_animal,session}(:,:,curr_pyrcell));
    axis off;
end
set(gcf,'Name',['Animal ' num2str(curr_animal) ', cell ' num2str(curr_pyrcell)]);


%% Lever template from average


% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);




% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step

% For this purpose - template is just half of the trials in last 3 days
lever_avg_template_trials = cell(8,15);
lever_avg_template = cell(8,1);
lever_avg_template_corr = cell(8,15);
lever_avg_avg_template_corr = nan(8,15);
activity_template_all = cell(8,1);
activity_template_corr = cell(8,15);
for curr_animal = 1:8
   
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    % pick random subset of trials from last 3 days to make lever template
    template_sessions = sessions(end-2:end);
    curr_lever_template_trials = cell(1,max(sessions));
    curr_lever_template_trials(template_sessions) = ...
        cellfun(@(x) randsample(size(x,1),round(size(x,1)/2)), ...
        lever_used_crop(curr_animal,template_sessions),'uni',false);
    lever_avg_template_trials(curr_animal,1:max(sessions)) = ...
        curr_lever_template_trials;
    
    % make lever template
    lever_avg_use = cellfun(@(x,y) x(y,:), ...
        lever_used_crop(curr_animal,1:max(sessions)), ...
        curr_lever_template_trials,'uni',false);
    curr_template = nanmean(vertcat(lever_avg_use{:})); 
    lever_avg_template{curr_animal} = curr_template;
    
    % get the lever correlation for trials not used to make template
    lever_avg_template_corr_grid = cellfun(@(x) ...
        corrcoef([x' curr_template']), ...
        lever_used_crop(curr_animal,sessions),'uni',false);    
    lever_avg_template_corr(curr_animal,sessions) = cellfun(@(x,y) ...
        x(end,setdiff(1:size(x,2)-1,y)),lever_avg_template_corr_grid, ...
        curr_lever_template_trials(sessions),'uni',false);
    
    % get lever correlation avg of day with template
    curr_lever_avg = cellfun(@(x,y) nanmean(x(setdiff(1:size(x,1),y),:)), ...
        lever_used_crop(curr_animal,sessions), ...
        curr_lever_template_trials(sessions),'uni',false);
    curr_lever_avg_cat = vertcat(curr_lever_avg{:});
    lever_avg_avg_template_corr_grid = corrcoef( ...
        [curr_lever_avg_cat' curr_template']);
    lever_avg_avg_template_corr(curr_animal,sessions) = ...
        lever_avg_avg_template_corr_grid(end,1:end-1);
    
    % make activity template from same trials used for lever template
    movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
        movement_activity_crop(curr_animal,sessions), ...
        curr_lever_template_trials(sessions),'uni',false);
    
    activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
    activity_template_all{curr_animal} = activity_template;
    
    % get the activity correlation for trials not used to make template
    act_avg_template_corr_grid = cellfun(@(x) ...
        corrcoef([x activity_template]), ...
        movement_activity_crop(curr_animal,sessions),'uni',false);    
    activity_template_corr(curr_animal,sessions) = cellfun(@(x,y) ...
        x(end,setdiff(1:size(x,2)-1,y)),act_avg_template_corr_grid, ...
        curr_lever_template_trials(sessions),'uni',false);    
    
    curr_animal
end

% plot each animal avg/trial template correlation
avg_trial_template_corr = cellfun(@nanmean, lever_avg_template_corr);
figure;
for i = 1:8
    subplot(3,3,i); hold on
   plot(lever_avg_avg_template_corr(i,:),'k');
   plot(avg_trial_template_corr(i,:),'r');
   title(i)
end
subplot(3,3,9); hold on
plot(1,1,'.k')
plot(1,1,'.r');
legend({'Avg-template' 'Trial-template'})



% plot the activity v lever template correlations for each animal
figure; hold on;
for earlylate = 1:2
    
    if earlylate == 1
        curr_sessions = 1:3;
        plot_col = 'k';
    elseif earlylate == 2
        curr_sessions = 8:14;
        plot_col = 'r';
    end
    
    use_animals = [1 2 3 4 5 6 7 8];
    for curr_animal = use_animals;
        total_lever_cc = cell(8,1);
        total_act_cc = cell(8,1);
        total_act = cell(8,1);
        
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = intersect(curr_sessions,sessions);
        
        cat_lever_cc = horzcat(lever_avg_template_corr{curr_animal,use_sessions});
        cat_act_cc = horzcat(activity_template_corr{curr_animal,use_sessions});
        
        curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
            movement_activity_crop(curr_animal,use_sessions), ...
            lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
        curr_act_cat = horzcat(curr_act{:});
               
        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;
        total_act{curr_animal} = curr_act_cat;
        
        total_lever_cc_cat = horzcat(total_lever_cc{:});
        total_act_cc_cat = horzcat(total_act_cc{:});
        
        corr_bins = [-1:0.1:1-0.1 Inf];
        [n bins] = histc(total_lever_cc_cat,corr_bins);
        m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
        s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
        
        subplot(3,3,curr_animal); hold on;
        errorbar(corr_bins(1:end-1),nanmean(m,1), ...
            nanmean(s,1),plot_col,'linewidth',2);
        xlabel('Lever template correlation')
        ylabel('Activity template correlation')
        title(curr_animal);        
    end
end

% combine all animals to get total activity and correlations
%use_animals = [1 6 7 8];
use_animals = [4 2];
%use_animals = [1 2 4 6 7 8];
corr_fig = figure; hold on;
act_fig = figure;
for earlylate = 1:2
    
    % they need to be sorted for down the road
    use_animals = sort(use_animals);
    
    if earlylate == 1
        curr_sessions = 1:3;
        plot_col = 'k';
    elseif earlylate == 2
        curr_sessions = 8:14;
        plot_col = 'r';
    end
    
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_act = cell(8,1);    
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = intersect(curr_sessions,sessions);

        cat_lever_cc = horzcat(lever_avg_template_corr{curr_animal,use_sessions});
        cat_act_cc = horzcat(activity_template_corr{curr_animal,use_sessions});
        
        curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
            movement_activity_crop(curr_animal,use_sessions), ...
            lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
        curr_act_cat = horzcat(curr_act{:});
               
        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;
        total_act{curr_animal} = curr_act_cat;        
    end
      
    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    
    corr_bins = [-1:0.1:1-0.1 Inf];
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    figure(corr_fig);
    errorbar(corr_bins(1:end-1),nanmean(m,1),nanmean(s,1),plot_col,'linewidth',2);
    xlabel('Lever template correlation')
    ylabel('Activity template correlation')
    title(['Animals ' num2str(use_animals)])
    %title(['Sessions ' num2str(curr_sessions(1)) ':' ...
    %    num2str(curr_sessions(end))]);
    
    % plot the average activity of cells in bins

    act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
    
    [n bin_split] = cellfun(@(x) histc(x,act_corr_bins),total_lever_cc,'uni',false);
    act_bins = horzcat(bin_split{:});
    act_grp_all = cell(8,1);
    for curr_animal = use_animals
        act_grp = nan(size(total_act{curr_animal},1),length(act_corr_bins));
        act_grp(:,unique(bin_split{curr_animal})) = ...
            grpstats([total_act{curr_animal}]',bin_split{curr_animal})';
        act_grp_reshape = reshape(act_grp,move_time_framespan+2*leeway+1, ...
            [],length(act_corr_bins));
        act_grp_all{curr_animal} = act_grp_reshape;
    end
    
    act_grp_cat = cat(2,act_grp_all{:});
    
    activity_template_cat = vertcat(activity_template_all{use_animals});
    activity_template_reshape = reshape(activity_template_cat, ...
        move_time_framespan+2*leeway+1,[])';
    [c ci] = max(activity_template_reshape,[],2);
    use_cells = find(c > 0.01);
    [asdf sort_idx] = sort(ci(use_cells));
    
    activity_template_reshape_norm = bsxfun(@times,activity_template_reshape,1./c);
    
    figure(act_fig);
    for i = 1:length(act_corr_bins)-1
        subplot(2,length(act_corr_bins),i+((earlylate-1)*length(act_corr_bins)));
        
        imagesc(act_grp_cat(:,use_cells(sort_idx),i)');
        
        colormap(hot);
        caxis([0 0.2]);
        axis off
        title([num2str(act_corr_bins(i)) ...
            ':' num2str(act_corr_bins(i+1))]);
    end
    subplot(2,length(act_corr_bins),length(act_corr_bins)+ ...
        ((earlylate-1)*length(act_corr_bins)));
    imagesc(activity_template_reshape(use_cells(sort_idx),:));
    caxis([0 0.2]);
    colormap(hot);
    axis off
    title('Template')
    colormap(hot);
    
    set(gcf,'Name',['Animals ' num2str(use_animals)]);
    %set(gcf,'Name',['Sessions ' num2str(curr_sessions(1)) ':' ...
    %    num2str(curr_sessions(end))]);
end
figure(corr_fig);
legend({'Early' 'Late'})

% plot the average lever trace for each used animal in each bin
use_animals = 1:8;
for earlylate = 1:2
    
    if earlylate == 1
        curr_sessions = 1:3;
        figure;
        col = [0 0 1];
    elseif earlylate == 2
        curr_sessions = 8:14;
        col = [1 0 0];
    end
        
    total_lever_cc = cell(8,1);
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = intersect(curr_sessions,sessions);

        cat_lever_cc = horzcat(lever_avg_template_corr{curr_animal,use_sessions});

        total_lever_cc{curr_animal} = cat_lever_cc;      
    end
      
    total_lever_cc_cat = horzcat(total_lever_cc{:});
       
    act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
    
    [n act_bins] = histc(total_lever_cc_cat,act_corr_bins);
    
    bin_split = mat2cell(act_bins,1,cellfun(@length,total_lever_cc));
    
    lever_used_bin = nan(length(act_corr_bins)-1,302,8);
    lever_used_bin_sem = nan(length(act_corr_bins)-1,302,8);
    for curr_animal = use_animals
        
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = intersect(curr_sessions,sessions);
        
        lever_used_crop_nontemplate = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
            lever_used_crop(curr_animal,use_sessions), ...
            lever_avg_template_trials(curr_animal,use_sessions),'uni',false);
        lever_cat = vertcat(lever_used_crop_nontemplate{:});
        
        lever_used_bin(unique(bin_split{curr_animal}),:,curr_animal) = ...
            grpstats(lever_cat,bin_split{curr_animal});
        lever_used_bin_sem(unique(bin_split{curr_animal}),:,curr_animal) = ...
            grpstats(lever_cat,bin_split{curr_animal},'sem');
        
    end
    
    for curr_animal = use_animals
        for i = 1:length(act_corr_bins)-1
            subplot(length(act_corr_bins),8,8*(i-1)+curr_animal);hold on;
            plot(lever_used_bin(i,:,curr_animal),'color',col,'linewidth',2);
            
            if earlylate == 1
                plot(lever_avg_template{curr_animal}, ...
                    '--k','linewidth',2)
            end           
            
            if i == 1
                title(curr_animal)
            end
            
            if curr_animal == 1
                ylabel(num2str(act_corr_bins(i)));
            end

            axis off
        end
    end
end


%% 10k/rewarded move classification check
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

pyr_class_new = cell(8,15);
gad_class_new = cell(8,15);

pyr_class_rewardmove_new = cell(8,15);
gad_class_rewardmove_new = cell(8,15);

for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
  
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    clear classified_cells
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        'final_class/' animal '_classified_cells_caEvents_10k'])
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_new(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_new(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    
    clear classified_cells
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        'final_class/' animal '_classified_cells_caEvents_rewardedmove_10k'])
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_rewardmove_new(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_rewardmove_new(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    
    curr_animal
    
end

old_class_cat = cell(14,1);
new_class_cat = cell(14,1);
new_rewardmove_class_cat = cell(14,1);

use_animals = [1 2 4 5 6 7 8];
for i = 1:14
   old_class_cat{i} = horzcat(pyr_class_all{use_animals,i});
   new_class_cat{i} = vertcat(pyr_class_new{use_animals,i});
   new_rewardmove_class_cat{i} = vertcat(pyr_class_rewardmove_new{use_animals,i});
end
figure; hold on;
plot(cellfun(@nanmean,old_class_cat),'k','linewidth',2);
plot(cellfun(@nanmean,new_class_cat),'r','linewidth',2);
plot(cellfun(@nanmean,new_rewardmove_class_cat),'b','linewidth',2);
xlabel('Session')
ylabel('Fraction of cells')
title('Movement-related cells')
legend({'Old' 'New' 'New rewarded move'})

%% Spatial clustering

% get all roi distances and gad/pyr unfilled indicies

% set pixels to microns ratio
um_px_x = 590/512;
um_px_y = 500/512;

roi_distances = cell(8,1);
pyr_unfilled_idx = cell(8,1);
gad_unfilled_idx = cell(8,1);

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    pyr_unfilled_idx{curr_animal} = find(pyr_cells & ~filled_cells);
    gad_unfilled_idx{curr_animal} = find(gad_cells & ~filled_cells);
    
    % load ROIs and get ROI centroids and distances
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    roi_dir = dir(roi_path);
    roi_filenames = {roi_dir(:).name};
    roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
    roi_files_sort = sort(roi_filenames(roi_files));
    load([roi_path filesep roi_files_sort{end}],'-MAT');
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for curr_roi  = 1:num_rois
        [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
            polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
    end
    roi_center_x = roi_center_x*um_px_x;
    roi_center_y = roi_center_y*um_px_y;
    
    roi_dist = nan(length(roi_center_x));
    for k = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
                (roi_center_y(k) - roi_center_y(j))^2);
            roi_dist(k,j) = curr_dist;
        end
    end
    
    roi_distances{curr_animal} = roi_dist;
    
    disp(['Distances: animal ' animal]) 
    
end


% % get the timing preference of cells for timing vs. clustering
% moveonset_aligned = cell(8,15);
% animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';
% pyr_time_diff_grid = cell(8,15);
% aligned_max = cell(8,15);
% for curr_animalday = animaldays
%     movetime_trace = nan(size(movement_trace_all{curr_animalday}));
%     curr_moves = move_epoch_frames{curr_animalday};
%     
%     curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
%         rewarded_movements{curr_animalday});
%     
%     curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
%     curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
%     
%     curr_rewardframes = cellfun(@(x) intersect(x, ...
%         reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
%     curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
%     
%     % Movement onset/offset reward +/- frames
%     fs_b = 5;fs_f = 85;
%     curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
%     curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
%     curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
%        
%     elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
%         length(movetime_trace)),curr_moveonset_surround, ...
%         curr_reward_surround,curr_moveoffset_surround);
%     curr_moveonset_surround(elim_moves) = [];
%     curr_reward_surround(elim_moves) = [];
%     curr_moveoffset_surround(elim_moves) = [];
%      
%     % to get just activity onsets (ignore sustained, to peak)
%     event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
%         diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
%     pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
%     pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
%     
%     curr_pyrclass = pyr_class_all{curr_animalday};
%       
%     use_activity = pyr_activity_all{curr_animalday};%(curr_pyrclass,:);
%     
%     % store activity 
%     curr_moveonset_aligned = cellfun(@(x) use_activity(:,x), ...
%        curr_moveonset_surround,'uni',false);
%     moveonset_aligned{curr_animalday} = nanmean(cat(3,curr_moveonset_aligned{:}),3);
%     [c ci] = max(moveonset_aligned{curr_animalday},[],2);
%     
%     pyr_time_diff_grid{curr_animalday} = ...
%         abs(repmat(ci,1,length(ci)) - repmat(ci',length(ci),1));
%     aligned_max{curr_animalday} = c;
% 
%     curr_animalday/max(animaldays)
% 
% end


% check clustering of movement cells
max_cutoff = 0.1;

pyr_move_dist = cell(8,15);
pyr_move_time_diff = cell(8,15);
pyr_move_dist_p = nan(8,15);
aligned_max_grid = cell(8,15);
cell_dist_p = cell(8,15);
shuff_dist = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for session = sessions
        
        curr_move_cells = pyr_unfilled_idx{curr_animal}( ...
            pyr_class_all{curr_animal,session});
        curr_dist = AP_itril(roi_distances{curr_animal}( ...
            curr_move_cells,curr_move_cells),-1);
        pyr_move_dist{curr_animal,session} = curr_dist;
        
        curr_aligned_max = aligned_max{curr_animal,session};
        aligned_max_grid{curr_animal,session} = ...
            AP_itril(+(curr_aligned_max(pyr_class_all{curr_animal,session})>max_cutoff)* ...
            +(curr_aligned_max(pyr_class_all{curr_animal,session})>max_cutoff)',-1);
        
        curr_move_cell_dist_mean = nanmean(curr_dist);
        
        pyr_move_time_diff{curr_animal,session} = ... 
            AP_itril(pyr_time_diff_grid{curr_animal,session}( ...
            pyr_class_all{curr_animal,session}, ...
            pyr_class_all{curr_animal,session}),-1);
    
        % get probability distribution for random cells
        num_shuff = 1000;
        rand_pyr_unfilled = pyr_unfilled_idx{curr_animal}( ...
            randi(length(pyr_unfilled_idx{curr_animal}), ...
            length(curr_move_cells),num_shuff));
        
        rand_pyr_reshape = permute(rand_pyr_unfilled,[1 3 2]);
        rand_pyr_grid_1 = repmat(rand_pyr_reshape,1, ...
            size(rand_pyr_reshape,1));
        rand_pyr_grid_2 = permute(rand_pyr_grid_1,[2 1 3]);
        
        rand_pyr_idx = sub2ind(size(roi_distances{curr_animal}), ...
            rand_pyr_grid_1,rand_pyr_grid_2);
        rand_pyr_idx_tril = arrayfun(@(x) AP_itril( ...
            rand_pyr_idx(:,:,x),-1),1:num_shuff,'uni',false);
        
        rand_pyr_dist = cellfun(@(x) nanmean(roi_distances{curr_animal}(x)) ,...
            rand_pyr_idx_tril);
        
        shuff_dist{curr_animal,session} = rand_pyr_dist;
                
        curr_dist_rank = tiedrank([rand_pyr_dist curr_move_cell_dist_mean]);
        
        pyr_move_dist_p(curr_animal,session) = curr_dist_rank(end)/(num_shuff+1);
        
        cell_dist_p{curr_animal,session} = nan(size(curr_dist));
        for j = 1:length(curr_dist)
            curr_cell_rank = tiedrank([rand_pyr_dist curr_dist(j)]);
            cell_dist_p{curr_animal,session}(j) = ...
                curr_cell_rank(end)/(num_shuff+1);
        end           
        
    end
    
    curr_animal
    
end

% plot mean/std and real values of distances
figure;
for curr_animal = [1 2 4 5 6 7 8];
    curr_sessions = find(cellfun(@(x) ~isempty(x), shuff_dist(curr_animal,:)));
    curr_shuff_mean = cellfun(@nanmean,shuff_dist(curr_animal,curr_sessions));
    curr_shuff_ci_low = cellfun(@(x) abs(nanmean(x) - prctile(x,[2.5])), ...
        shuff_dist(curr_animal,curr_sessions));
    curr_shuff_ci_high = cellfun(@(x) abs(nanmean(x) - prctile(x,[97.5])), ...
        shuff_dist(curr_animal,curr_sessions));
    subplot(3,3,curr_animal); hold on;
    errorbar(curr_sessions,curr_shuff_mean,curr_shuff_ci_low, ...
        curr_shuff_ci_high,'ok','linewidth',2,'linestyle','none', ...
        'MarkerSize',5,'MarkerFaceColor','k')
    plot(curr_sessions,cellfun(@nanmean,pyr_move_dist(curr_animal,curr_sessions)), ...
        'or','MarkerSize',5,'MarkerFaceColor','r');
end


% % plot grids of significantly close/far days
% figure;
% subplot(2,1,1);
% imagesc(pyr_move_dist_p < 0.025)
% colormap(gray)
% xlabel('Session')
% ylabel('Animal')
% title('Significantly close')
% subplot(2,1,2);
% imagesc(pyr_move_dist_p > 0.975)
% xlabel('Session')
% ylabel('Animal')
% title('Significantly far')
% 
% % plot timing vs. distance
% use_animals = [1 2 4 5 6 7 8];
% use_sessions = 8:14;
% 
% use_pairs = vertcat(aligned_max_grid{use_animals,use_sessions})>0;
% %dist_cat = vertcat(pyr_move_dist{use_animals,use_sessions});
% dist_cat = vertcat(cell_dist_p{use_animals,use_sessions});
% time_cat = vertcat(pyr_move_time_diff{use_animals,use_sessions});
% %use_pairs = true(size(dist_cat));
% time_bin = [0:5:90];
% [n bins] = histc(time_cat(use_pairs),time_bin);
% dist_binned_m = grpstats(dist_cat(use_pairs),bins);
% dist_binned_s = grpstats(dist_cat(use_pairs),bins,'sem');
% figure;errorbar(dist_binned_m,dist_binned_s,'k','linewidth',2)
% set(gca,'XTick',1:2:length(time_bin)-1)
% set(gca,'XTickLabel',round(time_bin(1:2:end)/2.8)/10);

































