%% Template-based correlations

stages = {[1:3] [10:14]};

stages_activity_template_corr_mean = cell(length(analysis),2);
stages_activity_template_corr_sem = cell(length(analysis),2);

for curr_animal = 1:length(analysis)
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    curr_stages = cellfun(@(x) intersect(find(use_sessions),x),stages,'uni',false);
    
    % Concatenate movement/activity into learning stages
    % if no days available from given stage, exclude animal
    if any(~cellfun(@(x) any(ismember(x,1:length(analysis(curr_animal).im))),curr_stages));
        continue
    end  
    stages_col = {'k','r'};
    
    stages_movement = cellfun(@(sessions) cell2mat(cellfun(@(lever,trials) horzcat(cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),lever(trials  & ...
        cellfun(@(y) ~isempty(y),lever)),'uni',false)')), ...
        {analysis(curr_animal).lever(sessions).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(sessions).cued_movement_trials},'uni',false)),curr_stages,'uni',false);
    
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
  
    framerate = 28;
    stages_activity = cellfun(@(sessions) cell2mat(cellfun(@(activity) reshape(permute( ...
        activity(:,analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000),use_cells), ...
        [2 3 1]),[],size(activity,1)),{analysis(curr_animal).im(sessions).move_onset_aligned}, ...
        'uni',false)),curr_stages,'uni',false);
    
    % Subtract the mean from activity (to eliminate non-specific activity)
    stages_activity_meansub = cellfun(@(x) bsxfun(@minus, ...
        x,nanmean(x,2)),stages_activity,'uni',false);
    
    % Make template from the last (expert) stage
    template_trials = randperm(size(stages_movement{end},2),round(size(stages_movement{end},2)/2));
    
    movement_template = nanmean(stages_movement{end}(:,template_trials),2);
    activity_template = nanmean(stages_activity{end}(:,template_trials),2);
    
    % Get correlations between trials and templates (excluding trials used to
    % create the templates)
    valid_trials = cellfun(@(x) 1:size(x,2),stages_movement,'uni',false);
    valid_trials{end}(template_trials) = [];
    
    stages_movement_template_corr_grid = cellfun(@(movement,trials) ...
        corrcoef([movement_template,movement(:,trials)]),stages_movement,valid_trials,'uni',false);
    
    stages_activity_template_corr_grid = cellfun(@(activity,trials) ...
        corrcoef([activity_template,activity(:,trials)],'rows','pairwise'), ...
        stages_activity,valid_trials,'uni',false);
    
    stages_movement_template_corr = cellfun(@(x) x(1,2:end), ...
        stages_movement_template_corr_grid,'uni',false);
    
    stages_activity_template_corr = cellfun(@(x) x(1,2:end), ...
        stages_activity_template_corr_grid,'uni',false);
    
    % Bin movement correlations
    corr_bin = [-1:0.2:1];
    
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    [~,bin_idx] = cellfun(@(x) histc(x,corr_bin_use),stages_movement_template_corr,'uni',false);
    bins_used = cellfun(@unique,bin_idx,'uni',false);
    
%     stages_activity_template_corr_mean = cellfun(@(corr,bins) grpstats(corr,bins,'mean'), ...
%         stages_activity_template_corr,bin_idx,'uni',false);
%     stages_activity_template_corr_sem = cellfun(@(corr,bins) grpstats(corr,bins,'sem'), ...
%         stages_activity_template_corr,bin_idx,'uni',false);
%     
    
    curr_stage_mean = repmat({nan(length(corr_bin_use),1)},1,2);
    curr_stage_sem = repmat({nan(length(corr_bin_use),1)},1,2);
    for curr_stage = 1:length(stages)
       
        curr_mean = grpstats(stages_activity_template_corr{curr_stage},bin_idx{curr_stage},'mean');
        curr_sem = grpstats(stages_activity_template_corr{curr_stage},bin_idx{curr_stage},'sem');
        
        curr_stage_mean{curr_stage}(bins_used{curr_stage}) = curr_mean;
        curr_stage_sem{curr_stage}(bins_used{curr_stage}) = curr_sem;
        
    end
    
    stages_activity_template_corr_mean(curr_animal,:) = curr_stage_mean;
    stages_activity_template_corr_sem(curr_animal,:) = curr_stage_sem;
    
    disp(curr_animal);
end


% 
% % Plot stages
% figure; hold on
% for curr_stage = 1:length(stages)
%     errorbar(corr_bin_plot(bins_used{curr_stage}), ...
%         stages_activity_template_corr_mean{curr_stage}, ...
%         stages_activity_template_corr_sem{curr_stage},'color',stages_col{curr_stage});
% end


% Plot stages combine animals
stages_col = {'k','r'};
figure; hold on
for curr_stage = 1:length(stages)
    errorbar(corr_bin_plot,nanmean(horzcat( ...
        stages_activity_template_corr_mean{:,curr_stage}),2), ...
        nanmean(horzcat( ...
        stages_activity_template_corr_sem{:,curr_stage}),2), ...
        'color',stages_col{curr_stage});
end


%% Trial-by-trial within/across sessions

stages = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
%stages = {[1:3],[10:14]};
framerate = 28;

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % custom
    %use_frames = 15:105;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).im);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials & ...
        cellfun(@(y) ~isempty(y),movement)),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    n_session_movements = cellfun(@(x) size(x,2),use_movement);
    movement_corr = cell(n_sessions);
    movement_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_movement{:})),n_session_movements,n_session_movements);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    use_activity = cellfun(@(activity) reshape(permute( ...
        activity(:,use_frames,use_cells), ...
        [2 3 1]),[],size(activity,1)), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Don't use cells that are ever NaN (unless the whole trial is a NaN -
    % then it'll just be ignored later)
    cat_activity = horzcat(use_activity{:});
    bad_cells = any(isnan(cat_activity(:,~all(isnan(cat_activity),1))),2);
    use_activity = cellfun(@(x) x(~bad_cells,:),use_activity,'uni',false);
    
    activity_corr = cell(n_sessions);
    activity_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_activity{:})),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    % Use sessions within stage that exist
    stages_sessions = cellfun(@(x) intersect(find(use_sessions),x), ...
        stages(use_stages),'uni',false);
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages_sessions);
        for curr_stage_2 = 1:length(stages_sessions);
            
            sessions_1 = stages_sessions{curr_stage_1};
            sessions_2 = stages_sessions{curr_stage_2};
                  
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
        activity_corr_stages_mean_animals{stage_1,stage_2} = nanmean(curr_mean);
        
        curr_sem = vertcat(activity_corr_stages_sem{stage_1,stage_2,:});
        activity_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem);
        
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


%% Trial-by-trial within/across sessions (FLIP: bin by activity, get move)

%stages = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
%stages = {[1:3],[10:14]};
stages = num2cell(1:14);
framerate = 28;

move_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
move_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % custom
    %use_frames = 15:105;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).im);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials & ...
        cellfun(@(y) ~isempty(y),movement)),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    n_session_movements = cellfun(@(x) size(x,2),use_movement);
    movement_corr = cell(n_sessions);
    movement_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_movement{:})),n_session_movements,n_session_movements);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    
    use_activity = cellfun(@(activity) reshape(permute( ...
        activity(:,use_frames,use_cells), ...
        [2 3 1]),[],size(activity,1)), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Don't use cells that are ever NaN (unless the whole trial is a NaN -
    % then it'll just be ignored later)
    cat_activity = horzcat(use_activity{:});
    bad_cells = any(isnan(cat_activity(:,~all(isnan(cat_activity),1))),2);
    use_activity = cellfun(@(x) x(~bad_cells,:),use_activity,'uni',false);
    
    activity_corr = cell(n_sessions);
    activity_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_activity{:})),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,51);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    % Use sessions within stage that exist
    stages_sessions = cellfun(@(x) intersect(find(use_sessions),x), ...
        stages(use_stages),'uni',false);
    
    % Get rid of trials with nans in activity correlation
    movement_corr = cellfun(@(x,y) x(~isnan(y)),movement_corr,activity_corr,'uni',false);
    activity_corr = cellfun(@(x) x(~isnan(x)),activity_corr,'uni',false);
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),activity_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages_sessions);
        for curr_stage_2 = 1:length(stages_sessions);
            
            sessions_1 = stages_sessions{curr_stage_1};
            sessions_2 = stages_sessions{curr_stage_2};
                  
            curr_num = grpstats(vertcat(movement_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'numel');
            curr_mean = grpstats(vertcat(movement_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'mean');
            curr_sem = grpstats(vertcat(movement_corr{sessions_1,sessions_2}), ...
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
    
    move_corr_stages_mean(:,:,curr_animal) = curr_stage_mean;
    move_corr_stages_sem(:,:,curr_animal) = curr_stage_sem;
    
    disp(curr_animal);
    
end

% Combine animals
move_corr_stages_mean_animals = cell(length(stages));
move_corr_stages_sem_animals = cell(length(stages));
for stage_1 = 1:length(stages)
    for stage_2 = 1:length(stages)
   
        curr_mean = vertcat(move_corr_stages_mean{stage_1,stage_2,:});
        move_corr_stages_mean_animals{stage_1,stage_2} = nanmean(curr_mean);
        
        curr_sem = vertcat(move_corr_stages_sem{stage_1,stage_2,:});
        move_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem);
        
    end
end

figure; hold on;
col = repmat(linspace(0,0.7,length(stages))',1,3);
for i = 1:length(stages)
    errorbar(corr_bin_plot,move_corr_stages_mean_animals{i,i}, ...
        move_corr_stages_sem_animals{1,1},'color',col(i,:),'linewidth',2);
end
errorbar(corr_bin_plot,move_corr_stages_mean_animals{1,end}, ...
        move_corr_stages_mean_animals{1,end},'color','b','linewidth',2);
ylabel('Pairwise movement correlation');
xlabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],stages,'uni',false), ...
    ['Across sessions ' num2str(stages{1}) ', ' num2str(stages{end})]]);



%% Trial-by-trial within/across sessions (movement on/off concat)

stages = {[1:4],[11:14]};

framerate = 28;

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);
    
    cat_time = 1500;
    
    % to correspond to movement
    use_frames_onset = analysis(curr_animal).surrounding_frames >= 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(cat_time/1000);
    use_frames_offset = analysis(curr_animal).surrounding_frames > ...
        -framerate*(cat_time/1000) & ...
        analysis(curr_animal).surrounding_frames < 0;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).im);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    rewarded_movement_duration = cellfun(@(x) cellfun(@length,x), ....
        {analysis(curr_animal).lever(use_sessions).rewarded_movement},'uni',false);
    
    time_cutoff_movements = cellfun(@(x) x > cat_time+1, ...
        rewarded_movement_duration,'uni',false);
    
    good_trials = cellfun(@(x) cellfun(@(x) ~isempty(x),x), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement},'uni',false);
    
    % Trial-by-trial movement correlation
    use_movement = cellfun(@(movement,cue_trials,time_trials,good_trials) cell2mat(cellfun(@(x) ...
        x([1:cat_time,end-cat_time:end]),movement(cue_trials & time_trials & ...
        good_trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials}, ...
        time_cutoff_movements,good_trials,'uni',false);
    
    n_session_movements = cellfun(@(x) size(x,2),use_movement);
    movement_corr = cell(n_sessions);
    movement_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_movement{:})),n_session_movements,n_session_movements);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    use_act_trials = cellfun(@(cued,timed,good) timed(cued & good), ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials}, ...
        time_cutoff_movements,good_trials,'uni',false);
    
    use_activity = cellfun(@(on_act,off_act,trials) reshape(permute( ...
        [on_act(trials,use_frames_onset,use_cells), ...
        off_act(trials,use_frames_offset,use_cells)] > 0, ...
        [2 3 1]),[],sum(trials)), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned}, ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned}, ...
        use_act_trials,'uni',false);
    
    activity_corr = cell(n_sessions);
    activity_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_activity{:})),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,20);
    corr_bin_use = corr_bin;
    corr_bin(end) = Inf;
    corr_bin(1) = -Inf;
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    % Use sessions within stage that exist
    stages_sessions = cellfun(@(x) intersect(find(use_sessions),x), ...
        stages(use_stages),'uni',false);
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages_sessions);
        for curr_stage_2 = 1:length(stages_sessions);
            
            sessions_1 = stages_sessions{curr_stage_1};
            sessions_2 = stages_sessions{curr_stage_2};
                  
            curr_mean = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'mean');
            curr_sem = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'sem');
            
            bins_used = unique(vertcat(stages_bin_idx{sessions_1,sessions_2}));
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
        activity_corr_stages_mean_animals{stage_1,stage_2} = nanmean(curr_mean);
        
        curr_sem = vertcat(activity_corr_stages_sem{stage_1,stage_2,:});
        activity_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem);
        
    end
end

figure; hold on;
errorbar(corr_bin_plot,activity_corr_stages_mean_animals{1,1}, ...
    activity_corr_stages_sem_animals{1,1},'k');
errorbar(corr_bin_plot,activity_corr_stages_mean_animals{2,2}, ...
    activity_corr_stages_sem_animals{2,2},'b');
errorbar(corr_bin_plot,activity_corr_stages_mean_animals{1,2}, ...
    activity_corr_stages_sem_animals{1,2},'r');

legend({'Early','Late','Across'})


%% Template-based correlations (movement on/off concat)

stages = {[1:4] [11:14]};

framerate = 28;

stages_activity_template_corr_mean = cell(length(analysis),2);
stages_activity_template_corr_sem = cell(length(analysis),2);

for curr_animal = 1:length(analysis)
    
    cat_time = 1500;
    
    % to correspond to movement
    use_frames_onset = analysis(curr_animal).surrounding_frames >= 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(cat_time/1000);
    use_frames_offset = analysis(curr_animal).surrounding_frames > ...
        -framerate*(cat_time/1000) & ...
        analysis(curr_animal).surrounding_frames < 0;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).lever);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    curr_stages = cellfun(@(x) intersect(find(use_sessions),x),stages,'uni',false);
    
    rewarded_movement_duration = cellfun(@(x) cellfun(@length,x), ....
        {analysis(curr_animal).lever(use_sessions).rewarded_movement},'uni',false);
    
    time_cutoff_movements = cellfun(@(x) x > cat_time+1, ...
        rewarded_movement_duration,'uni',false);
    
    good_trials = cellfun(@(x) cellfun(@(x) ~isempty(x),x), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement},'uni',false);
    
    % Concatenate movement/activity into learning stages
    % if no days available from given stage, exclude animal
    if any(~cellfun(@(x) any(ismember(x,1:length(analysis(curr_animal).im))),curr_stages));
        continue
    end  
    stages_col = {'k','r'};
    
    stages_movement = cellfun(@(sessions) ...
        cell2mat(cellfun(@(movement,cue_trials,time_trials,good_trials) ...
        horzcat(cell2mat(cellfun(@(x) ...
        x([1:cat_time,end-cat_time:end]), ...
        movement(cue_trials & time_trials & good_trials),'uni',false)')), ...
        {analysis(curr_animal).lever(sessions).rewarded_movement}, ...
        {analysis(curr_animal).lever(sessions).cued_movement_trials}, ...
        time_cutoff_movements(sessions),good_trials(sessions),'uni',false)),curr_stages,'uni',false);
    
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    use_act_trials = cellfun(@(cued,timed,good) timed(cued & good), ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials}, ...
        time_cutoff_movements,good_trials,'uni',false);
    
    stages_activity = cellfun(@(sessions) cell2mat( ...
        cellfun(@(on_act,off_act,trials) reshape(permute( ...
        [on_act(trials,use_frames_onset,use_cells), ...
        off_act(trials,use_frames_offset,use_cells)] > 0, ...
        [2 3 1]),[],sum(trials)), ...
        {analysis(curr_animal).im(sessions).move_onset_aligned}, ...
        {analysis(curr_animal).im(sessions).move_onset_aligned}, ...
        use_act_trials(sessions),'uni',false)),curr_stages,'uni',false);
    
    % Subtract the mean from activity (to eliminate non-specific activity)
    stages_activity_meansub = cellfun(@(x) bsxfun(@minus, ...
        x,nanmean(x,2)),stages_activity,'uni',false);
    
    % Make template from the last (expert) stage
    template_trials = randperm(size(stages_movement{end},2),round(size(stages_movement{end},2)/2));
    
    movement_template = nanmean(stages_movement{end}(:,template_trials),2);
    activity_template = nanmean(stages_activity_meansub{end}(:,template_trials),2);
    
    % Get correlations between trials and templates (excluding trials used to
    % create the templates)
    valid_trials = cellfun(@(x) 1:size(x,2),stages_movement,'uni',false);
    valid_trials{end}(template_trials) = [];
    
    stages_movement_template_corr_grid = cellfun(@(movement,trials) ...
        corrcoef([movement_template,movement(:,trials)]),stages_movement,valid_trials,'uni',false);
    
    stages_activity_template_corr_grid = cellfun(@(activity,trials) ...
        corrcoef([activity_template,activity(:,trials)],'rows','pairwise'), ...
        stages_activity_meansub,valid_trials,'uni',false);
    
    stages_movement_template_corr = cellfun(@(x) x(1,2:end), ...
        stages_movement_template_corr_grid,'uni',false);
    
    stages_activity_template_corr = cellfun(@(x) x(1,2:end), ...
        stages_activity_template_corr_grid,'uni',false);
    
    % Bin movement correlations
    corr_bin = [-1:0.2:1];
    
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    [~,bin_idx] = cellfun(@(x) histc(x,corr_bin_use),stages_movement_template_corr,'uni',false);
    bins_used = cellfun(@unique,bin_idx,'uni',false);
    
%     stages_activity_template_corr_mean = cellfun(@(corr,bins) grpstats(corr,bins,'mean'), ...
%         stages_activity_template_corr,bin_idx,'uni',false);
%     stages_activity_template_corr_sem = cellfun(@(corr,bins) grpstats(corr,bins,'sem'), ...
%         stages_activity_template_corr,bin_idx,'uni',false);
%     
    
    curr_stage_mean = repmat({nan(length(corr_bin_use),1)},1,2);
    curr_stage_sem = repmat({nan(length(corr_bin_use),1)},1,2);
    for curr_stage = 1:length(stages)
       
        curr_mean = grpstats(stages_activity_template_corr{curr_stage},bin_idx{curr_stage},'mean');
        curr_sem = grpstats(stages_activity_template_corr{curr_stage},bin_idx{curr_stage},'sem');
        
        curr_stage_mean{curr_stage}(bins_used{curr_stage}) = curr_mean;
        curr_stage_sem{curr_stage}(bins_used{curr_stage}) = curr_sem;
        
    end
    
    stages_activity_template_corr_mean(curr_animal,:) = curr_stage_mean;
    stages_activity_template_corr_sem(curr_animal,:) = curr_stage_sem;
    
    disp(curr_animal);
end


% 
% % Plot stages
% figure; hold on
% for curr_stage = 1:length(stages)
%     errorbar(corr_bin_plot(bins_used{curr_stage}), ...
%         stages_activity_template_corr_mean{curr_stage}, ...
%         stages_activity_template_corr_sem{curr_stage},'color',stages_col{curr_stage});
% end


% Plot stages combine animals
stages_col = {'k','r'};
figure; hold on
for curr_stage = 1:length(stages)
    errorbar(corr_bin_plot,nanmean(horzcat( ...
        stages_activity_template_corr_mean{:,curr_stage}),2), ...
        nanmean(horzcat( ...
        stages_activity_template_corr_sem{:,curr_stage}),2), ...
        'color',stages_col{curr_stage});
end


%% Trial-trial correlation mean/slope across days

stages = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
stages = num2cell(1:14);
framerate = 28;

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);
    
    movement_start_time = 1001;
    movement_use_time = 3000;
    
    % to correspond to movement
    use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
        analysis(curr_animal).surrounding_frames < framerate*(movement_use_time/1000);
    
    % custom
    %use_frames = 15:105;
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).im);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials & ...
        cellfun(@(y) ~isempty(y),movement)),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(use_sessions).cued_movement_trials},'uni',false);
    
    n_session_movements = cellfun(@(x) size(x,2),use_movement);
    movement_corr = cell(n_sessions);
    movement_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_movement{:})),n_session_movements,n_session_movements);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    %use_cells = end_up_rois{curr_animal};
    %use_cells = any(classified_rois(curr_animal).movement,2) & ...
    %    ~any(classified_rois(curr_animal).quiescent,2);
    
    use_activity = cellfun(@(activity) reshape(permute( ...
        activity(:,use_frames,use_cells), ...
        [2 3 1]),[],size(activity,1)), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Don't use cells that are ever NaN (unless the whole trial is a NaN -
    % then it'll just be ignored later)
    cat_activity = horzcat(use_activity{:});
    bad_cells = any(isnan(cat_activity(:,~all(isnan(cat_activity),1))),2);
    use_activity = cellfun(@(x) x(~bad_cells,:),use_activity,'uni',false);
    
    activity_corr = cell(n_sessions);
    activity_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_activity{:})),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    % Use sessions within stage that exist
    stages_sessions = cellfun(@(x) intersect(find(use_sessions),x), ...
        stages(use_stages),'uni',false);
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages_sessions);
        for curr_stage_2 = 1:length(stages_sessions);
            
            sessions_1 = stages_sessions{curr_stage_1};
            sessions_2 = stages_sessions{curr_stage_2};
                  
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

% Get slope of movement/activity correlation in positive domain
positive_bins = find(corr_bin > 0);
corr_fit = cellfun(@(x) nan(1,2),activity_corr_stages_mean,'uni',false);
use_animaldays = cellfun(@(x) ~isempty(x),activity_corr_stages_mean);
corr_fit(use_animaldays) = cellfun(@(x) polyfit(positive_bins(~isnan(x(positive_bins))), ...
    x(positive_bins(~isnan(x(positive_bins)))),1), ...
    activity_corr_stages_mean(use_animaldays),'uni',false);
corr_slope = cellfun(@(x) x(1),corr_fit);

figure;
subplot(1,2,1);
imagesc(nanmean(corr_slope,3));
colormap(hot);
axis square;
xlabel('Day');
ylabel('Day');
title('Slope of positive movement-correlation activity');

% Get mean of positive domain correlations
positive_bins = find(corr_bin > 0 & all(~isnan(vertcat(activity_corr_stages_mean{:})),1));
positive_bin_mean = nan(size(activity_corr_stages_mean));
use_animaldays = cellfun(@(x) ~isempty(x),activity_corr_stages_mean);
positive_bin_mean(use_animaldays) = cellfun(@(x) nanmean(x(positive_bins)), ...
    activity_corr_stages_mean(use_animaldays));

subplot(1,2,2);
imagesc(nanmean(positive_bin_mean,3));
colormap(hot)
axis square;
xlabel('Day');
ylabel('Day');
title('Mean of positive movement-correlation activity');


%% Trial-trial correlation mean/slope across days (ITI movements)

% Get ITI movement/activity


% Get activity and movements
epoch_movement = cell(length(data),14);
epoch_activity = cell(length(data),14);

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
        
        % Set criteria for pulling out movements (in ms)
        use_time = 3000;
        use_frames = (use_time/1000)*28;
        
        min_move_time = 1000;
        max_move_time = 3500;
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
        
        epoch_movement{curr_animal,curr_session} = ...
            timed_movements;
        epoch_activity{curr_animal,curr_session} = ...
            timed_activity;
        
        disp(curr_session);
        
    end
end

% Get movement/activity correlation

stages = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
stages = num2cell(1:14);
framerate = 28;

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);   
    
    good_sessions = cellfun(@(x) ~isempty(x),epoch_activity(curr_animal,:));
    
    % Define trials to use
    n_sessions = find(good_sessions, 1, 'last' );
    use_sessions = find(good_sessions(1:n_sessions));
    
    % Trial-by-trial movement correlation
    use_movement = epoch_movement(curr_animal,1:n_sessions);
    
    n_session_movements = cellfun(@(x) size(x,2),use_movement);
    movement_corr = ...
        mat2cell(corrcoef(horzcat(use_movement{:})),n_session_movements,n_session_movements);
    movement_corr(find(eye(size(movement_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),movement_corr(find(eye(size(movement_corr)))),'uni',false);
    movement_corr(find(~eye(size(movement_corr)))) = ...
        cellfun(@(x) x(:),movement_corr(find(~eye(size(movement_corr)))),'uni',false);
    
    % Trial-by-trial activity correlation
    use_cells = true(size(analysis(curr_animal).im(1).move_onset_aligned,3),1);
    %use_cells = end_up_rois{curr_animal};
    %use_cells = any(classified_rois(curr_animal).movement,2) & ...
    %    ~any(classified_rois(curr_animal).quiescent,2);
    
    use_activity = cellfun(@(x) reshape(x,[],size(x,3)), ...
        epoch_activity(curr_animal,1:n_sessions),'uni',false);
    
    % Don't use cells that are ever NaN (unless the whole trial is a NaN -
    % then it'll just be ignored later)
    cat_activity = horzcat(use_activity{:});
    bad_cells = any(isnan(cat_activity(:,~all(isnan(cat_activity),1))),2);
    use_activity(use_sessions) = cellfun(@(x) x(~bad_cells,:),use_activity(use_sessions),'uni',false);
    
    activity_corr = ...
        mat2cell(corrcoef(horzcat(use_activity{:})),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,21);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    % Use sessions within stage that exist
    stages_sessions = cellfun(@(x) intersect(find(use_sessions),x), ...
        stages(use_stages),'uni',false);
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = 1:length(stages_sessions);
        for curr_stage_2 = 1:length(stages_sessions);
            
            sessions_1 = stages_sessions{curr_stage_1};
            sessions_2 = stages_sessions{curr_stage_2};
                  
            curr_mean = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'mean');
            curr_sem = grpstats(vertcat(activity_corr{sessions_1,sessions_2}), ...
                vertcat(stages_bin_idx{sessions_1,sessions_2}),'sem');
            
            bins_used = unique(vertcat(stages_bin_idx{sessions_1,sessions_2}));
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

% Get slope of movement/activity correlation in positive domain
positive_bins = find(corr_bin > 0);
corr_fit = cellfun(@(x) nan(1,2),activity_corr_stages_mean,'uni',false);
use_animaldays = cellfun(@(x) ~isempty(x),activity_corr_stages_mean);
corr_fit(use_animaldays) = cellfun(@(x) polyfit(positive_bins(~isnan(x(positive_bins))), ...
    x(positive_bins(~isnan(x(positive_bins)))),1), ...
    activity_corr_stages_mean(use_animaldays),'uni',false);
corr_slope = cellfun(@(x) x(1),corr_fit);

figure;
subplot(1,2,1);
imagesc(nanmean(corr_slope,3));
colormap(hot);
axis square;
xlabel('Day');
ylabel('Day');
title('Slope of positive movement-correlation activity');

% Get mean of positive domain correlations
%positive_bins = find(corr_bin > 0 & all(~isnan(vertcat(activity_corr_stages_mean{:})),1));
positive_bins = find(corr_bin > 0);
positive_bin_mean = nan(size(activity_corr_stages_mean));
use_animaldays = cellfun(@(x) ~isempty(x),activity_corr_stages_mean);
positive_bin_mean(use_animaldays) = cellfun(@(x) nanmean(x(positive_bins)), ...
    activity_corr_stages_mean(use_animaldays));

subplot(1,2,2);
imagesc(nanmean(positive_bin_mean,3));
colormap(hot)
axis square;
xlabel('Day');
ylabel('Day');
title('Mean of positive movement-correlation activity');

% Plot line plots
activity_corr_stages_mean_animals = cell(length(stages));
activity_corr_stages_sem_animals = cell(length(stages));
for stage_1 = 1:length(stages)
    for stage_2 = 1:length(stages)
   
        curr_mean = vertcat(activity_corr_stages_mean{stage_1,stage_2,:});
        activity_corr_stages_mean_animals{stage_1,stage_2} = nanmean(curr_mean);
        
        curr_sem = vertcat(activity_corr_stages_sem{stage_1,stage_2,:});
        activity_corr_stages_sem_animals{stage_1,stage_2} = nanmean(curr_sem);
        
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


%% Mutual information?





























