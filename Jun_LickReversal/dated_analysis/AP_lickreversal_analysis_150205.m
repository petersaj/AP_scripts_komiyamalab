%% GENERAL USE: Get activity by epoch (following cells use this)

similarity = 'correlation'; % 'correlation', 'pca'
trial_time = [0 2]; % min/max trial time to look at, seconds from odor onset
sig_cells_only = false; % force nonsignificant cells to have value of 0

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});


% Identify significantly active cells from baseline if selected
if sig_cells_only
    % Go through each epoch, compare to baseline (significant & above)
    trial_epochs = 4;
    baseline_multiepoch_sigdiff = cell(length(analysis),1);
    for curr_animal = 1:length(animals)
        for curr_session = 1:length(analysis.epoch_activity_CL{curr_animal});
            for curr_reversal = 1:length(analysis.epoch_activity_CL{curr_animal}{curr_session})
                for curr_epoch = 1:trial_epochs;
                    
                    frames = analysis.epoch_frames{curr_animal}{curr_session};
                    curr_epoch_frames = frames(1):frames(2) > (curr_epoch-1)*frames(2)/trial_epochs & ...
                        frames(1):frames(2) < curr_epoch*frames(2)/trial_epochs;
                    
                    baseline_frames = frames(1):frames(2) < 0;
                    
                    curr_activity = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_reversal};
                    
                    baseline_multiepoch_sigdiff{curr_animal}{curr_session}{curr_reversal}{curr_epoch} = ...
                        AP_signrank_matrix( ...
                        permute(nanmean(curr_activity(:,baseline_frames,:),2),[1 3 2]), ...
                        permute(nanmean(curr_activity(:,curr_epoch_frames,:),2),[1 3 2])) ...
                        < 0.01 & ...
                        nanmean(permute(nanmean(curr_activity(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
                        nanmean(permute(nanmean(curr_activity(:,curr_epoch_frames,:),2),[1 3 2]),1);
                    
                end
            end
        end
        disp(curr_animal);
    end
end
 

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);
all_sessions = cell(length(animals),1);

% Get activity
for curr_animal = 1:length(mice)
        
    curr_slot = 1;
    for curr_session = 1:length(data_all(curr_animal).im);
        
        % define reversal (use only the first so that it's settled in)
        for curr_rev = 1:length(analysis.epoch_contingencies{curr_animal}{curr_session});
            
            revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
                size(analysis.odor_trials{curr_animal}{curr_session},1)];
            % if there are < 10 trials in reversal, assume early/not real
            % reversal and use the next one
            rev_trials = diff(revs);
            if rev_trials(curr_rev) < 10
                curr_rev = curr_rev + 1;
            end
            
            % use time in the trial defined by user
            use_frames = data_all(curr_animal).im(curr_session).framerate*trial_time;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > use_frames(1) & curr_frames <= use_frames(2);
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (:,curr_odor_frames,:);
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(nanmean(any(curr_act,2),1),[3 1 2]);
            % if selected, force cells not significantly active during
            % selected trial time to be 0
            if sig_cells_only 
                trial_epochs_ends = linspace(0,4,trial_epochs+1); % the ends of the epochs (sec)
                trial_epochs_use = trial_epochs_ends(2:end) > floor(trial_time(1)) & ...
                    trial_epochs_ends(2:end) <= ceil(trial_time(2));
                
                curr_sigcells = ...
                    any(vertcat(baseline_multiepoch_sigdiff{curr_animal} ...
                    {curr_session}{curr_rev}{trial_epochs_use}),1)';                              
                
                curr_active_trials(~curr_sigcells) = 0;
            end
            
            curr_contingency = analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);           
            % if control animal, don't use anything past session 10 or if
            % contingency 2 has appeared (because reversal in these
            % animals, but it's late - I'm not using this data for now)
            if ismember(curr_animal,find(ctrl_animals));
                if curr_session > 10 || curr_contingency == 2
                    continue
                end
            end
            
            all_active_trials{curr_animal,curr_slot} = curr_active_trials;
            all_contingencies{curr_animal}(curr_slot) = curr_contingency;
            all_sessions{curr_animal}(curr_slot) = curr_session;
            
            % 'Late' for reversal animals = anything including and after
            % the second reversal
            if ismember(curr_animal,find(rev_animals));
                all_latereversal{curr_animal}(curr_slot) = ...
                    sum(diff(all_contingencies{curr_animal}) ~= 0) > 1;
            else
                % 'Late for control animals = anything after the 6th session
                all_latereversal{curr_animal}(curr_slot) = ...
                    curr_session >= 6;
            end
            curr_slot = curr_slot + 1;
            
        end
    end
end




%% Slide 1-2) similarity between contingency epochs
% generalized for these features:
% - measurement of similarity
% - time used to gauge active on trial
% - activity on trial threshold or compared to baseline

% If similarity measure is PC distance, get number of PCs to use
if strcmp(similarity,'pca')
    pc_var = cell(length(mice),1);
    for curr_animal = 1:length(mice);
        
        curr_act = horzcat(all_active_trials{curr_animal,:});
        [coeff score latent] = princomp(curr_act');
        pc_var{curr_animal} = cumsum(latent)/sum(latent);
    end
    % use the median number of PCs that give > 80% variance
    pcs = ceil(median(cellfun(@(x) find(x > 0.8,1),pc_var)));
    disp(['PCs used: ' num2str(pcs)]);
end

% Distance of the first n PCs WITHIN EARLY/LATE
activity_distance = cell(length(mice),1);
for curr_animal = 1:length(mice);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    switch similarity       
        case 'pca'
            % euclidean distance between selected PCs
            [coeff score latent] = princomp(curr_act');
            pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
            pc_dist_norm = (pc_dist-nanmean(AP_itril(pc_dist,-1)))./nanstd(AP_itril(pc_dist,-1));
            use_similarity = -pc_dist_norm;
            
        case 'correlation'
            % Pearson's correlation
            activity_corr = corrcoef(curr_act);
            use_similarity = activity_corr;
    end
    
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    
    % split/save distances, dim1 = animal, dim2 = early/late, 
    % dim3 = 1v1 / 2v2 / 1v2
    
    activity_distance{curr_animal,1,1} = ...
        AP_itril(use_similarity(~all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),-1);
    
    activity_distance{curr_animal,2,1} = ...
        AP_itril(use_similarity(all_latereversal{curr_animal,:} & cont1_session, ...
        all_latereversal{curr_animal,:} & cont1_session),-1);
    
    activity_distance{curr_animal,1,2} = ...
        AP_itril(use_similarity(~all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),-1);
    
    activity_distance{curr_animal,2,2} = ...
        AP_itril(use_similarity(all_latereversal{curr_animal,:} & cont2_session, ...
        all_latereversal{curr_animal,:} & cont2_session),-1);
    
    activity_distance{curr_animal,1,3} = ...
        reshape(use_similarity(~all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,2,3} = ...
        reshape(use_similarity(all_latereversal{curr_animal,:} & cont1_session, ...
        all_latereversal{curr_animal,:} & cont2_session),[],1);
    
end

activity_distance_mean_withinrev = cellfun(@nanmean,activity_distance);


% Distance of the first n PCs ACROSS EARLY/LATE
activity_distance = cell(length(mice),1);
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    switch similarity       
        case 'pca'
            % euclidean distance between selected PCs
            [coeff score latent] = princomp(curr_act');
            pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
            pc_dist_norm = (pc_dist-nanmean(AP_itril(pc_dist,-1)))./nanstd(AP_itril(pc_dist,-1));
            use_similarity = -pc_dist_norm;
            
        case 'correlation'
            % Pearson's correlation
            activity_corr = corrcoef(curr_act);
            use_similarity = activity_corr;
    end
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    
    % split/save distances, dim1 = animal, 
    % dim2 = 1v1 / 2v2 / 1v2 / 2v1 (former = early, latter = late)
      
    activity_distance{curr_animal,1} = ...
        reshape(use_similarity(all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),[],1);
    
    activity_distance{curr_animal,2} = ...
        reshape(use_similarity(all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,3} = ...
        reshape(use_similarity(all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,4} = ...
        reshape(use_similarity(all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),[],1);
    
end

activity_distance_mean_acrossrev = cellfun(@nanmean,activity_distance);


% Plot by contingency and area 
alm_distance_withinrev = activity_distance_mean_withinrev(alm_animals & rev_animals,:,:);
alm_distance_acrossrev = activity_distance_mean_acrossrev(alm_animals & rev_animals,:,:);

alm_distance_withinrev_reshape = permute(alm_distance_withinrev,[2 1 3]);
alm_distance_acrossrev_reshape = cat(1,permute(alm_distance_acrossrev(:,1:3),[3 1 2]), ...
    cat(3,zeros(1,size(alm_distance_acrossrev,1),2),permute(alm_distance_acrossrev(:,4),[3 1 2])));

alm_distance_boxplot = cat(1,alm_distance_withinrev_reshape,alm_distance_acrossrev_reshape);


pmm_distance_withinrev = activity_distance_mean_withinrev(pmm_animals & rev_animals,:,:);
pmm_distance_acrossrev = activity_distance_mean_acrossrev(pmm_animals & rev_animals,:,:);

pmm_distance_withinrev_reshape = permute(pmm_distance_withinrev,[2 1 3]);
pmm_distance_acrossrev_reshape = cat(1,permute(pmm_distance_acrossrev(:,1:3),[3 1 2]), ...
    cat(3,zeros(1,size(pmm_distance_acrossrev,1),2),permute(pmm_distance_acrossrev(:,4),[3 1 2])));

pmm_distance_boxplot = cat(1,pmm_distance_withinrev_reshape,pmm_distance_acrossrev_reshape);

figure;
p1 = subplot(2,5,1:4);
aboxplot(alm_distance_boxplot,'colorgrad','green_down')
ylabel('Similarity')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early 1/late 2' 'Across early 2/late 1'});
title('ALM')

p2 = subplot(2,5,6:9);
aboxplot(pmm_distance_boxplot,'colorgrad','green_down')
ylabel('Similarity')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early 1/late 2' 'Across early 2/late 1'});
title('PMM')

% Plot control as early v late contingency 1
alm_distance_withinrev_ctrl = activity_distance_mean_withinrev(alm_animals & ctrl_animals,:,1);
alm_distance_acrossrev_ctrl = activity_distance_mean_acrossrev(alm_animals & ctrl_animals,1);
alm_distance_ctrl_boxplot = [alm_distance_withinrev_ctrl alm_distance_acrossrev_ctrl];

pmm_distance_withinrev_ctrl = activity_distance_mean_withinrev(pmm_animals & ctrl_animals,:,1);
pmm_distance_acrossrev_ctrl = activity_distance_mean_acrossrev(pmm_animals & ctrl_animals,1);
pmm_distance_ctrl_boxplot = [pmm_distance_withinrev_ctrl pmm_distance_acrossrev_ctrl];

subplot(2,5,5); hold on;
aboxplot(alm_distance_ctrl_boxplot,'colorgrad','green_down');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early' 'Late' 'Across'});
ylim(get(p1,'ylim'));
title('ALM control');

subplot(2,5,10); hold on;
aboxplot(pmm_distance_ctrl_boxplot,'colorgrad','green_down');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early' 'Late' 'Across'});
ylim(get(p1,'ylim'));
title('PMM control');

set(gcf,'Name',['User settings: similarity=' similarity ...
    ', trial time=[' num2str(trial_time) '], sig cells = ' num2str(+sig_cells_only)]);



%% Get the mean activity in different conditions

% at the moment requires the last cell be run to avoid copy/paste mess

mean_activity = {};
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    
    % split/save distances, dim1 = animal, 
    % dim2 = contingency (1/2) % dim3 = early/late
      
    mean_activity{curr_animal,1,1} = ...
        nanmean(curr_act(:,~all_latereversal{curr_animal,:} & cont1_session));
    mean_activity{curr_animal,2,1} = ...
        nanmean(curr_act(:,~all_latereversal{curr_animal,:} & cont2_session));
    mean_activity{curr_animal,1,2} = ...
        nanmean(curr_act(:,all_latereversal{curr_animal,:} & cont1_session));
    mean_activity{curr_animal,2,2} = ...
        nanmean(curr_act(:,all_latereversal{curr_animal,:} & cont2_session));
    
end
mean_activity_grouped = cellfun(@nanmedian,mean_activity);


alm_rev_meanact = mean_activity_grouped(alm_animals & rev_animals,:,:);
pmm_rev_meanact = mean_activity_grouped(pmm_animals & rev_animals,:,:);

figure;aboxplot(permute(alm_rev_meanact,[2 1 3]))
title('ALM');
set(gca,'XTickLabel',{'Early' 'Late'});
ylabel('Mean reliability')
legend({'Cont 1' 'Cont 2'});

figure;aboxplot(permute(pmm_rev_meanact,[2 1 3]))
title('PMM')
set(gca,'XTickLabel',{'Early' 'Late'});
ylabel('Mean reliability')
legend({'Cont 1' 'Cont 2'});



%% Slide 3) Similarity of activity by time (sessions)

% If similarity measure is PC distance, get number of PCs to use
if strcmp(similarity,'pca')
    pc_var = cell(length(mice),1);
    for curr_animal = 1:length(mice);
        
        curr_act = horzcat(all_active_trials{curr_animal,:});
        [coeff score latent] = princomp(curr_act');
        pc_var{curr_animal} = cumsum(latent)/sum(latent);
    end
    % use the median number of PCs that give > 80% variance
    pcs = ceil(median(cellfun(@(x) find(x > 0.8,1),pc_var)));
    disp(['PCs used: ' num2str(pcs)]);
end

% Distance by session / contingency / time
max_sessions = max(cellfun(@length,analysis.epoch_activity_CL))+1;
max_rev = max(cellfun(@(x) sum(diff(x) ~= 0),all_contingencies))+1;
rev_grp_all = nan(length(mice),max_sessions,max_rev);
for curr_animal = 1:length(mice);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    switch similarity       
        case 'pca'
            % euclidean distance between selected PCs
            [coeff score latent] = princomp(curr_act');
            pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
            pc_dist_norm = (pc_dist-nanmean(AP_itril(pc_dist,-1)))./nanstd(AP_itril(pc_dist,-1));
            use_similarity = -pc_dist_norm;
            
        case 'correlation'
            % Pearson's correlation
            activity_corr = corrcoef(curr_act);
            use_similarity = activity_corr;
    end
    
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    % Get dfference between variables
    time_difference = pdist2(all_sessions{curr_animal}',all_sessions{curr_animal}');
    n_reversals = cumsum(abs([0 diff(all_contingencies{curr_animal})]));
    reversal_difference = pdist2(n_reversals',n_reversals');
    contingency_difference = pdist2(all_contingencies{curr_animal}',all_contingencies{curr_animal}');
    earlylate_difference = pdist2(all_latereversal{curr_animal}',all_latereversal{curr_animal}');

    
    time_aoc = AP_itril(time_difference,-1);
    similarity_aoc = AP_itril(use_similarity,-1);
    contingency_aoc = AP_itril(contingency_difference,-1);
    reversal_aoc = AP_itril(reversal_difference,-1);
    earlylate_aoc = AP_itril(earlylate_difference,-1);
    
    % for ancova
    % aoctool(time_aoc,similarity_aoc,reversal_aoc)
    
    % group by reversal distance
    for curr_rev_dist = unique(reversal_aoc)';
        curr_rev_grp = reversal_aoc == curr_rev_dist;
        rev_grp_mean_similarity = grpstats(similarity_aoc(curr_rev_grp),time_aoc(curr_rev_grp));
        rev_grp_all(curr_animal,1+unique(time_aoc(curr_rev_grp)),curr_rev_dist+1) = ...
            rev_grp_mean_similarity;    
    end

end


figure; 

rev_grp_alm = permute(nanmean(rev_grp_all(rev_animals & alm_animals,:,:)),[2 3 1]);
subplot(2,1,1); hold on;
n_plot = find(any(~isnan(rev_grp_alm),2),1,'last');
col = jet(n_plot);
for i = 1:n_plot
   plot(rev_grp_alm(i,:),'color',col(i,:)); 
end
title('ALM (color = \DeltaReversal)');
xlabel('Session difference')
ylabel('Similarity')

rev_grp_pmm = permute(nanmean(rev_grp_all(rev_animals & pmm_animals,:,:)),[2 3 1]);
subplot(2,1,2); hold on;
n_plot = find(any(~isnan(rev_grp_pmm),2),1,'last');
col = jet(n_plot);
for i = 1:n_plot
   plot(rev_grp_pmm(i,:),'color',col(i,:)); 
end
title('PMM (color = \DeltaReversal)');
xlabel('Session difference')
ylabel('Similarity')


%% GENERAL USE: Get d' (based on all sessions)

%%%%% NEW (150619) d' was previously being calculated incorrectly, and I
%%%%% don't think specifically d' necessary anyway. Changed it to just hit
%%%%% rate - false alarm rate (has nice range of -1:1)

% define sliding window for d'
d_prime_trials = 10;
slide_window = ones(d_prime_trials,1);

d_prime = cell(size(mice));
for curr_animal = 1:length(mice);
    
    n_session_trials = cellfun(@(x) size(x,1),analysis.condition_trials{curr_animal});
    n_sessions = length(analysis.condition_trials{curr_animal});
    
    curr_conditions = vertcat(analysis.condition_trials{curr_animal}{:});
    curr_conditions_slide = conv2(+curr_conditions,slide_window,'same');
    hit_rate = curr_conditions_slide(:,1)./ ...
        (curr_conditions_slide(:,1) + curr_conditions_slide(:,4));
    fa_rate = curr_conditions_slide(:,3)./ ...
        (curr_conditions_slide(:,3) + curr_conditions_slide(:,2));
    curr_d_prime = hit_rate - fa_rate;
    d_prime{curr_animal} = mat2cell(curr_d_prime,n_session_trials,1);
    
end



%% Slide 4) d' vs. correlation

% take moving average of activity in trials and maybe manual inspection:
% the point here is to see if correlation drops around reversal but then
% comes back up with d'

% or maybe more directly, get correlation v. d' across days


% d': CL cont 1 vs CL cont 1
n_bins = 10;

d_prime_act_corr = nan(n_bins,n_bins,length(mice));
for curr_animal = 1:length(mice)
    
    curr_CL_d_prime = ...
        cell2mat(cellfun(@(x,y) x(y(:,1)), d_prime{curr_animal}, ...
        analysis.condition_trials{curr_animal}','uni',false));
    
    curr_act = cell2mat(cellfun(@(x,y) reshape(permute(vertcat( ...
        x{:}),[2 3 1]),size(x{end},2)*size(x{end},3),[]), ...
        analysis.epoch_activity_CL{curr_animal},'uni',false));
    
    curr_contingency = cell2mat(cellfun(@(a,b) cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
        num2cell(a),b,'uni',false)),analysis.epoch_contingencies{curr_animal}', ...
        analysis.epoch_activity_CL{curr_animal}','uni',false));
    
    %d_prime_edges = linspace(prctile(curr_CL_d_prime,10),max(curr_CL_d_prime),n_bins);
    d_prime_edges = prctile(curr_CL_d_prime,linspace(1,100,10));
    d_prime_edges(1) = -Inf;
    d_prime_edges(end) = Inf;
    [n bin] = histc(curr_CL_d_prime,d_prime_edges);
    
    % correlation by bin
    for curr_bin_1 = 1:n_bins
        for curr_bin_2 = 1:n_bins
            curr_total_corr = corrcoef([curr_act(:,bin == curr_bin_1) ...
                curr_act(:,bin == curr_bin_2)],'rows','complete');
            
            if curr_bin_1 ~= curr_bin_2
                curr_corr = curr_total_corr( ....
                    1:sum(bin == curr_bin_1),sum(bin == curr_bin_1)+1:end);
            else
                curr_corr = AP_itril(curr_total_corr,-1);
            end
            
            curr_median_corr = nanmedian(curr_corr(:));
            
            d_prime_act_corr(curr_bin_1,curr_bin_2,curr_animal) = curr_median_corr;           
        end
    end   
    disp(curr_animal);
end

%%% d': Licking CL cont 1 vs CL cont 1
d_prime_lick_corr = nan(n_bins,n_bins,length(mice));
for curr_animal = 1:length(mice)
    
    curr_CL_d_prime = ...
        cell2mat(cellfun(@(x,y) x(y(:,1)), d_prime{curr_animal}, ...
        analysis.condition_trials{curr_animal}','uni',false));
    
    curr_lick = cell2mat(cellfun(@(x,y) x(y(:,1),:), analysis.lick_rate{curr_animal}', ...
        analysis.condition_trials{curr_animal}','uni',false))';    
    
    curr_contingency = cell2mat(cellfun(@(a,b) cell2mat(cellfun(@(x,y) repmat(x,size(y,1),1), ...
        num2cell(a),b,'uni',false)),analysis.epoch_contingencies{curr_animal}', ...
        analysis.epoch_activity_CL{curr_animal}','uni',false));
    
    d_prime_edges = linspace(prctile(curr_CL_d_prime,10),max(curr_CL_d_prime),n_bins);
    d_prime_edges(1) = -Inf;
    [n bin] = histc(curr_CL_d_prime,d_prime_edges);
    
    % correlation by bin
    for curr_bin_1 = 1:n_bins
        for curr_bin_2 = 1:n_bins
            curr_total_corr = corrcoef([curr_lick(:,bin == curr_bin_1) ...
                curr_lick(:,bin == curr_bin_2)]);
            
            if curr_bin_1 ~= curr_bin_2
                curr_corr = curr_total_corr( ....
                    1:sum(bin == curr_bin_1),sum(bin == curr_bin_1)+1:end);
            else
                curr_corr = AP_itril(curr_total_corr,-1);
            end
            
            curr_median_corr = nanmedian(curr_corr(:));
            
            d_prime_lick_corr(curr_bin_1,curr_bin_2,curr_animal) = curr_median_corr;           
        end
    end   
    disp(curr_animal);
end

%%% d': IL vs CL 
n_bins = 10;

d_prime_act_corr = nan(n_bins,n_bins,length(mice));
for curr_animal = 1:length(mice)
    
    curr_d_prime = ...
        vertcat(d_prime{curr_animal}{:});
    
    curr_condition = vertcat(analysis.condition_trials{curr_animal}{:});
    
    curr_act = cell2mat(cellfun(@(act) reshape(permute(vertcat( ...
        act{:}),[2 3 1]),size(act{end},2)*size(act{end},3),[]), ...
        analysis.epoch_activity_odor{curr_animal},'uni',false));
    
    curr_contingency = vertcat(analysis.odor_trials{curr_animal}{:});
    
    d_prime_edges = linspace(prctile(curr_d_prime,10),max(curr_d_prime),n_bins);
    d_prime_edges(1) = -Inf;
    [n bin] = histc(curr_d_prime,d_prime_edges);
    
    % correlation by bin
    for curr_bin_1 = 1:n_bins
        for curr_bin_2 = 1:n_bins
            
            curr_use_1 = bin == curr_bin_1 & ...
                curr_condition(:,1) & curr_contingency(:,2) == 1;
            
            curr_use_2 = bin == curr_bin_2 & ...
                curr_condition(:,3) & curr_contingency(:,2) == 2;
            
            curr_total_corr = corrcoef([curr_act(:,curr_use_1) ...
                curr_act(:,curr_use_2)],'rows','complete');
            
                curr_corr = curr_total_corr( ....
                    1:sum(curr_use_1),sum(curr_use_1)+1:end);
           
            curr_median_corr = nanmedian(curr_corr(:));
            
            d_prime_act_corr(curr_bin_1,curr_bin_2,curr_animal) = curr_median_corr;           
        end
    end   
    disp(curr_animal);
end

%%% d': CL cont 1 vs CL cont 2
n_bins = 10;

d_prime_act_corr = nan(n_bins,n_bins,length(mice));
for curr_animal = 1:length(mice)
    
    curr_d_prime = ...
        vertcat(d_prime{curr_animal}{:});
    
    curr_condition = vertcat(analysis.condition_trials{curr_animal}{:});
    
    curr_act = cell2mat(cellfun(@(act) reshape(permute(vertcat( ...
        act{:}),[2 3 1]),size(act{end},2)*size(act{end},3),[]), ...
        analysis.epoch_activity_odor{curr_animal},'uni',false));
    
    curr_contingency = vertcat(analysis.odor_trials{curr_animal}{:});
    
    d_prime_edges = linspace(prctile(curr_d_prime,10),max(curr_d_prime),n_bins);
    d_prime_edges(1) = -Inf;
    [n bin] = histc(curr_d_prime,d_prime_edges);
    
    % correlation by bin
    for curr_bin_1 = 1:n_bins
        for curr_bin_2 = 1:n_bins
            
            curr_use_1 = bin == curr_bin_1 & ...
                curr_condition(:,1) & curr_contingency(:,2) == 1;
            
            curr_use_2 = bin == curr_bin_2 & ...
                curr_condition(:,1) & curr_contingency(:,2) == 2;
            
            curr_total_corr = corrcoef([curr_act(:,curr_use_1) ...
                curr_act(:,curr_use_2)],'rows','complete');
            
                curr_corr = curr_total_corr( ....
                    1:sum(curr_use_1),sum(curr_use_1)+1:end);
           
            curr_median_corr = nanmedian(curr_corr(:));
            
            d_prime_act_corr(curr_bin_1,curr_bin_2,curr_animal) = curr_median_corr;           
        end
    end   
    disp(curr_animal);
end


%% GENERAL USE: get activity by epoch, with d' cutoff

% define sliding window for d'
d_prime_trials = 10;
slide_window = ones(d_prime_trials,1);

d_prime = cell(size(mice));
for curr_animal = 1:length(mice);
    
    n_session_trials = cellfun(@(x) size(x,1),analysis.condition_trials{curr_animal});
    n_sessions = length(analysis.condition_trials{curr_animal});
    
    curr_conditions = vertcat(analysis.condition_trials{curr_animal}{:});
    curr_conditions_slide = conv2(+curr_conditions,slide_window,'same');
    hit_rate = curr_conditions_slide(:,1)./ ...
        (curr_conditions_slide(:,1) + curr_conditions_slide(:,4));
    fa_rate = curr_conditions_slide(:,3)./ ...
        (curr_conditions_slide(:,3) + curr_conditions_slide(:,2));
    curr_d_prime = hit_rate - fa_rate;
    d_prime{curr_animal} = mat2cell(curr_d_prime,n_session_trials,1);
    
end

n_bins = 10;
d_prime_use_trials_CL_epoch = cell(length(mice),1);
d_prime_use_trials_CR_epoch = cell(length(mice),1);
for curr_animal = 1:length(mice)
    
    % I think all this going in circles isn't actually necessary, can split
    % d' up earlier instead of unpacking and repacking
    curr_d_prime = ...
        vertcat(d_prime{curr_animal}{:});
    
    curr_condition = vertcat(analysis.condition_trials{curr_animal}{:});
       
    curr_odor = vertcat(analysis.odor_trials{curr_animal}{:});
    curr_contingency = curr_odor(:,2);
    
    d_prime_edges = linspace(-1,1,n_bins);
    d_prime_edges(1) = -Inf;
    d_prime_edges(end) = Inf;
    [n bin] = histc(curr_d_prime,d_prime_edges);
    
    d_prime_cutoff = d_prime_edges(end-3);
    d_prime_use_trials = curr_d_prime >= d_prime_cutoff;
    
    % CL
    d_prime_use_trials_CL = d_prime_use_trials(curr_condition(:,1));
    
    d_prime_use_trials_CL_session = mat2cell(d_prime_use_trials_CL, ...
        cellfun(@(x) size(vertcat(x{:}),1),analysis.epoch_activity_CL{curr_animal}));
    
    d_prime_use_trials_CL_epoch{curr_animal} = cellfun(@(x,y) mat2cell(x,cellfun(@(z) size(z,1),y)), ...
        d_prime_use_trials_CL_session', ...
        analysis.epoch_activity_CL{curr_animal},'uni',false);
    
    % CR
    d_prime_use_trials_CR = d_prime_use_trials(curr_condition(:,2));
    
    d_prime_use_trials_CR_session = mat2cell(d_prime_use_trials_CR, ...
        cellfun(@(x) size(vertcat(x{:}),1),analysis.epoch_activity_CR{curr_animal}));
    
    d_prime_use_trials_CR_epoch{curr_animal} = cellfun(@(x,y) mat2cell(x,cellfun(@(z) size(z,1),y)), ...
        d_prime_use_trials_CR_session', ...
        analysis.epoch_activity_CR{curr_animal},'uni',false);
    
    

    disp(curr_animal);
end



% I don't know why this was included here

similarity = 'pca'; % 'correlation', 'pca'
trial_time = [0 2]; % min/max trial time to look at, seconds from odor onset
sig_cells_only = false; % force nonsignificant cells to have value of 0

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});


% Identify significantly active cells from baseline if selected
if sig_cells_only
    % Go through each epoch, compare to baseline (significant & above)
    trial_epochs = 4;
    baseline_multiepoch_sigdiff = cell(length(animals),1);
    for curr_animal = 1:length(animals)
        for curr_session = 1:length(analysis.epoch_activity_CL{curr_animal});
            for curr_reversal = 1:length(analysis.epoch_activity_CL{curr_animal}{curr_session})
                for curr_epoch = 1:trial_epochs;
                    
                    frames = analysis.epoch_frames{curr_animal}{curr_session};
                    curr_epoch_frames = frames(1):frames(2) > (curr_epoch-1)*frames(2)/trial_epochs & ...
                        frames(1):frames(2) < curr_epoch*frames(2)/trial_epochs;
                    
                    baseline_frames = frames(1):frames(2) < 0;
                    
                    curr_activity = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_reversal}( ...
                        d_prime_use_trials_CL_epoch{curr_animal}{curr_session}{curr_reversal},:,:);
                    
                    baseline_multiepoch_sigdiff{curr_animal}{curr_session}{curr_reversal}{curr_epoch} = ...
                        AP_signrank_matrix( ...
                        permute(nanmean(curr_activity(:,baseline_frames,:),2),[1 3 2]), ...
                        permute(nanmean(curr_activity(:,curr_epoch_frames,:),2),[1 3 2])) ...
                        < 0.01 & ...
                        nanmean(permute(nanmean(curr_activity(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
                        nanmean(permute(nanmean(curr_activity(:,curr_epoch_frames,:),2),[1 3 2]),1);
                    
                end
            end
        end
        disp(curr_animal);
    end
end
 

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);
all_sessions = cell(length(animals),1);

% Get activity
for curr_animal = 1:length(mice)
        
    curr_slot = 1;
    for curr_session = 1:length(data_all(curr_animal).im);
        
        % define reversal (use only the first so that it's settled in)
        for curr_rev = 1:length(analysis.epoch_contingencies{curr_animal}{curr_session});
            
            revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
                size(analysis.odor_trials{curr_animal}{curr_session},1)];
            % if there are < 10 trials in reversal, assume early/not real
            % reversal and use the next one
            rev_trials = diff(revs);
            if rev_trials(curr_rev) < 10
                curr_rev = curr_rev + 1;
            end
            
            % use time in the trial defined by user
            use_frames = data_all(curr_animal).im(curr_session).framerate*trial_time;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > use_frames(1) & curr_frames <= use_frames(2);
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (d_prime_use_trials_CL_epoch{curr_animal}{curr_session}{curr_rev},curr_odor_frames,:);
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(nanmean(any(curr_act,2),1),[3 1 2]);
            % if selected, force cells not significantly active during
            % selected trial time to be 0
            if sig_cells_only 
                trial_epochs_ends = linspace(0,4,trial_epochs+1); % the ends of the epochs (sec)
                trial_epochs_use = trial_epochs_ends(2:end) > floor(trial_time(1)) & ...
                    trial_epochs_ends(2:end) <= ceil(trial_time(2));
                
                curr_sigcells = ...
                    any(vertcat(baseline_multiepoch_sigdiff{curr_animal} ...
                    {curr_session}{curr_rev}{trial_epochs_use}),1)';                              
                
                curr_active_trials(~curr_sigcells) = 0;
            end
            
            curr_contingency = analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);           
            % if control animal, don't use anything past session 10 or if
            % contingency 2 has appeared (because reversal in these
            % animals, but it's late - I'm not using this data for now)
            if ismember(curr_animal,find(ctrl_animals));
                if curr_session > 10 || curr_contingency == 2
                    continue
                end
            end
            
            all_active_trials{curr_animal,curr_slot} = curr_active_trials;
            all_contingencies{curr_animal}(curr_slot) = curr_contingency;
            all_sessions{curr_animal}(curr_slot) = curr_session;
            
            % 'Late' for reversal animals = anything including and after
            % the second reversal
            if ismember(curr_animal,find(rev_animals));
                all_latereversal{curr_animal}(curr_slot) = ...
                    sum(diff(all_contingencies{curr_animal}) ~= 0) > 1;
            else
                % 'Late for control animals = anything after the 6th session
                all_latereversal{curr_animal}(curr_slot) = ...
                    curr_session >= 6;
            end
            curr_slot = curr_slot + 1;
            
        end
    end
end



%% Slide 5+6) Jun requested analysis: frac active, etc
% 1) % hit cells across days in ctrl/rev
% -- for rev, between rev 1-2, beween rev 2-3, after 3rd reversal
% 2) % CR cells across days in ctrl/rev
%
% (also try using d' cutoff here)

% Convert epoch frames into seconds
framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);
epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% Get significant cells (fraction active during CL compared to baseline)
baseline_CL_sigcells = ...
    cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
    AP_chisquare_matrix( ...
    permute(any(odor(:,frames < 0,:),2),[1 3 2]), ...
    permute(any(odor(:,frames > 0 & ...
    frames < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,frames,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

baseline_CR_sigcells = ...
    cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
    AP_chisquare_matrix( ...
    permute(any(odor(:,frames < 0,:),2),[1 3 2]), ...
    permute(any(odor(:,frames > 0 & ...
    frames < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,frames,'uni',false), ...
    analysis.epoch_activity_CR,epoch_seconds,...
    'uni',false);

% Get significant cells (fraction active during CL compared to baseline)
% + only use d' cutoff trials
baseline_dprime_CL_sigcells = ...
        cellfun(@(odor,frames,trials) cellfun(@(odor,frames,trials) cellfun(@(odor,trials) ...
        AP_chisquare_matrix( ...
        permute(any(odor(trials,frames < 0,:),2),[1 3 2]), ...
        permute(any(odor(trials,frames > 0 & ...
        frames < 2,:),2),[1 3 2])) ...
        < 0.05,odor,trials,'uni',false),odor,frames,trials,'uni',false), ...
        analysis.epoch_activity_CL,epoch_seconds,...
        d_prime_use_trials_CL_epoch', ...
        'uni',false);
    
baseline_dprime_CR_sigcells = ...
        cellfun(@(odor,frames,trials) cellfun(@(odor,frames,trials) cellfun(@(odor,trials) ...
        AP_chisquare_matrix( ...
        permute(any(odor(trials,frames < 0,:),2),[1 3 2]), ...
        permute(any(odor(trials,frames > 0 & ...
        frames < 2,:),2),[1 3 2])) ...
        < 0.05,odor,trials,'uni',false),odor,frames,trials,'uni',false), ...
        analysis.epoch_activity_CR,epoch_seconds,...
        d_prime_use_trials_CR_epoch', ...
        'uni',false);

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});


use_animals = find(rev_animals & pmm_animals);
use_cont = 1;
trial_cutoff = 20; % if less than this trials, don't use day
    
curr_baseline_dprime_sigcells = cell(length(analysis),1);
curr_baseline_sigcells = cell(length(analysis),1);
for curr_animal_idx = 1:length(use_animals);
    
    curr_animal = use_animals(curr_animal_idx);    
    
    % d' cutoff
    curr_sig = baseline_dprime_CR_sigcells{curr_animal};
    curr_cont = analysis.epoch_contingencies{curr_animal};
    
    curr_use_epochs = cellfun(@(x) x == use_cont,curr_cont,'uni',false);
    curr_sig_use = cellfun(@(x,y) any(vertcat(x{y}),1),curr_sig,curr_use_epochs,'uni',false);
    
    curr_baseline_dprime_sigcells{curr_animal_idx} = ...
        NaN(length(curr_sig_use),size(vertcat(curr_sig_use{:}),2));
    use_sessions = find(cellfun(@(x) ~isempty(x),curr_sig_use));
    curr_baseline_dprime_sigcells{curr_animal_idx}(use_sessions,:) = vertcat(curr_sig_use{:});
    
    n_trials = cellfun(@(x,y) sum(vertcat(x{y}),1), ...
        d_prime_use_trials_CR_epoch{curr_animal},curr_use_epochs,'uni',false);
    under_trial_cutoff = cellfun(@(x) x < trial_cutoff,n_trials(use_sessions));
    curr_baseline_dprime_sigcells{curr_animal_idx}(use_sessions(under_trial_cutoff),:) = NaN;
    
    
    % all trials
    curr_sig = baseline_CR_sigcells{curr_animal};
    curr_cont = analysis.epoch_contingencies{curr_animal};
    
    curr_use_epochs = cellfun(@(x) x == use_cont,curr_cont,'uni',false);
    curr_sig_use = cellfun(@(x,y) any(vertcat(x{y}),1),curr_sig,curr_use_epochs,'uni',false);
    
    curr_baseline_sigcells{curr_animal_idx} = ...
        NaN(length(curr_sig_use),size(vertcat(curr_sig_use{:}),2));
    use_sessions = cellfun(@(x) ~isempty(x),curr_sig_use);
    curr_baseline_sigcells{curr_animal_idx}(use_sessions,:) = vertcat(curr_sig_use{:});  
    
    n_trials = cellfun(@(x,y) sum(vertcat(x{y}),1), ...
        d_prime_use_trials_CR_epoch{curr_animal},curr_use_epochs,'uni',false);
    under_trial_cutoff = cellfun(@(x) x < trial_cutoff,n_trials(use_sessions));
    curr_baseline_sigcells{curr_animal_idx}(use_sessions(under_trial_cutoff),:) = NaN;
    
end

curr_baseline_dprime_sigcells_mean = cellfun(@(x) ...
    nanmean(x,2),curr_baseline_dprime_sigcells,'uni',false);

curr_baseline_sigcells_mean = cellfun(@(x) ...
    nanmean(x,2),curr_baseline_sigcells,'uni',false);

figure; hold on;

max_sessions = max(cellfun(@length,curr_baseline_dprime_sigcells_mean));
data_nanpad = cell2mat(cellfun(@(x) [x;nan(max_sessions-length(x),1)],curr_baseline_dprime_sigcells_mean,'uni',false))';
errorbar(nanmean(data_nanpad),nanstd(data_nanpad)./sqrt(sum(~isnan(data_nanpad))),'r','linewidth',2)

max_sessions = max(cellfun(@length,curr_baseline_sigcells_mean));
data_nanpad = cell2mat(cellfun(@(x) [x;nan(max_sessions-length(x),1)],curr_baseline_sigcells_mean,'uni',false))';
errorbar(nanmean(data_nanpad),nanstd(data_nanpad)./sqrt(sum(~isnan(data_nanpad))),'k','linewidth',2)
    
ylabel('Fraction of active cells');
xlabel('Session');
legend({'D'' cutoff' 'All'})

% History of active cells
max_sessions = max(cellfun(@(x) size(x,1),curr_baseline_dprime_sigcells));

curr_baseline_dprime_sigcells_history = cell2mat(cellfun(@(x) ...
    padarray(sum(x(2:end,:) == 1 & cumsum(x(1:end-1,:)>0,1) > 0,2)./sum(x(2:end,:) == 1,2), ...
    max_sessions-size(x,1),NaN,'post'),curr_baseline_dprime_sigcells,'uni',false))';

curr_baseline_sigcells_history = cell2mat(cellfun(@(x) ...
    padarray(sum(x(2:end,:) == 1 & cumsum(x(1:end-1,:)>0,1) > 0,2)./sum(x(2:end,:) == 1,2), ...
    max_sessions-size(x,1),NaN,'post'),curr_baseline_sigcells,'uni',false))';

figure; hold on;
errorbar(nanmean(curr_baseline_dprime_sigcells_history), ...
    nanstd(curr_baseline_dprime_sigcells_history)./ ...
    sqrt(sum(~isnan(curr_baseline_dprime_sigcells_history))),'r','linewidth',2)

errorbar(nanmean(curr_baseline_sigcells_history), ...
    nanstd(curr_baseline_sigcells_history)./ ...
    sqrt(sum(~isnan(curr_baseline_sigcells_history))),'k','linewidth',2)

ylabel('Fraction of active cells active previously');
xlabel('Session');
legend({'D'' cutoff' 'All'})

% Sig cell correlation
max_sessions = max(cellfun(@(x) size(x,1),curr_baseline_dprime_sigcells));

curr_baseline_dprime_sigcells_corr = cellfun(@(x) padarray(corrcoef(x'), ...
    repmat(max_sessions-size(x,1),1,2),NaN,'post'),curr_baseline_dprime_sigcells,'uni',false);

curr_baseline_sigcells_corr = cellfun(@(x) padarray(corrcoef(x'), ...
    repmat(max_sessions-size(x,1),1,2),NaN,'post'),curr_baseline_sigcells,'uni',false);

figure;
subplot(2,2,1);
imagesc(nanmean(cat(3,curr_baseline_dprime_sigcells_corr{:}),3));
colormap(gray); caxis([0 1]);
title('D'' cutoff');
xlabel('Session');
ylabel('Session');
subplot(2,2,2);
imagesc(nanmean(cat(3,curr_baseline_sigcells_corr{:}),3));
colormap(gray); caxis([0 1]);
title('All trials')
xlabel('Session');
ylabel('Session');
subplot(2,2,3);
curr_diag = cell2mat(cellfun(@(x) diag(x,-1),curr_baseline_dprime_sigcells_corr,'uni',false));
errorbar(nanmean(curr_diag,2),nanstd(curr_diag,[],2)./sqrt(sum(~isnan(curr_diag),2)),'k','linewidth',2)
xlim([0 max_sessions+1])
title('1st diagonal')
ylabel('Fraction cells');
xlabel('Session');
subplot(2,2,4);
curr_diag = cell2mat(cellfun(@(x) diag(x,-1),curr_baseline_sigcells_corr,'uni',false));
errorbar(nanmean(curr_diag,2),nanstd(curr_diag,[],2)./sqrt(sum(~isnan(curr_diag),2)),'k','linewidth',2)
xlim([0 max_sessions+1])
title('1st diagonal')
ylabel('Fraction cells');
xlabel('Session');

%% For Jun: get timing of cell activity relative to lick onset
% Find activity/lick time in CL trials, and relative timing
%
% output: 
% activity_time (time of max activity relative to odor onset), 
% lick_time (time of lick relative to odor onset), 
% activity_lick_time (time of max activity relative to lick within odor)
%
% nested cell array: animal/session/reversal 
% Time in matricies: trials (rows) by sessions (columns)

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

activity_time = cell(length(analysis.epoch_activity_CL),1);
lick_time = cell(length(analysis.epoch_activity_CL),1);
for curr_animal = 1:length(analysis.epoch_activity_CL)
    for curr_session = 1:length(analysis.epoch_activity_CL{curr_animal})
        
        curr_lick = mat2cell(analysis.lick_rate{curr_animal}{curr_session}( ...
                analysis.condition_trials{curr_animal}{curr_session}(:,1),:), ...
                cellfun(@(x) size(x,1),analysis.epoch_activity_CL{curr_animal}{curr_session}), ...
                size(analysis.lick_rate{curr_animal}{curr_session},2));
            
        for curr_rev = 1:length(analysis.epoch_activity_CL{curr_animal}{curr_session});
            
            use_epoch_seconds = epoch_seconds{curr_animal}{curr_session} > 0 & ...
                epoch_seconds{curr_animal}{curr_session} < 2;
            curr_epoch_seconds = epoch_seconds{curr_animal}{curr_session}(use_epoch_seconds);
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev}( ...
                :,use_epoch_seconds,:);
            curr_act(isnan(curr_act)) = 0;
            [max_act,max_idx] = max(curr_act,[],2);
            act_idx = permute(max_idx,[1 3 2]);
            act_time = curr_epoch_seconds(act_idx);
            act_time(permute(max_act,[1 3 2]) == 0) = NaN;
            
            lick_thresh = 3; % Hz threshold for licking
            [max_lick,curr_lick_thresh] = max(curr_lick{curr_rev} > lick_thresh,[],2);
            curr_lick_thresh(max_lick == 0) = NaN;
            curr_lick_thresh_sec = curr_lick_thresh/100;
            
            activity_time{curr_animal}{curr_session}{curr_rev} = act_time;
            lick_time{curr_animal}{curr_session}{curr_rev} = curr_lick_thresh_sec;
                        
        end
    end
end
    

activity_lick_time = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x,y) ...
    bsxfun(@minus,x,y),x,y,'uni',false),x,y,'uni',false), ...
    activity_time,lick_time,'uni',false);
                

%% Slide 9) TK requested: day-back correlation of population activity across epochs

% Get correlation of mean odor activity
    
framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% % all trials
% odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
%     reshape(permute(nanmean(x(:,y > 0 & y < 2,:),1),[2 3 1]),[],1),x,'uni',false), ...
%     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);

% d' cutoff
odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
    reshape(permute(nanmean(x(z,y > 0 & y < 2,:),1),[2 3 1]),[],1),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);


odor_activity_CL_cat = cellfun(@(x) cell2mat(vertcat(x{:})'),odor_activity_CL,'uni',false);

epoch_sessions_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_sessions,'uni',false);
epoch_contingency_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_contingencies,'uni',false);

% There are random occurrences of the wrong contingency where there was not
% meant to be a reversal. This means that some epochs have very few trials
% but are only identifiable by this fact. Identify where these epochs occur
% and pull them out indivudally.
% --- but this happens a lot? and there isn't a clear cutoff ??? this makes
% this impossible.
% THE ONLY THING I CAN DO: if there under an arbitrary cutoff, then ignore
% it entirely. I'm making the cutoff now based what looks like a line
% separating number of CL trials.

num_CL_trials = cellfun(@(x) cellfun(@(x) size(x,1),vertcat(x{:})), ...
    analysis.epoch_activity_CL,'uni',false);
num_trial_cutoff = 20;
use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_CL_trials,'uni',false);

odor_activity_CL_cat_use = cellfun(@(x,y) x(:,y),odor_activity_CL_cat, ...
    use_epochs,'uni',false);
epoch_sessions_cat_use = cellfun(@(x,y) x(y),epoch_sessions_cat, ...
    use_epochs,'uni',false);
epoch_contingency_cat_use = cellfun(@(x,y) x(y),epoch_contingency_cat, ...
    use_epochs,'uni',false);

% Get reversals (after selecting out epochs)
reversals = cellfun(@(x) find(diff(x) ~= 0),epoch_contingency_cat_use,'uni',false);

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Get day-before correlation
act_corr = cellfun(@(x) corrcoef(x,'rows','complete'), ...
    odor_activity_CL_cat_use,'uni',false);
act_corr_oneback = cellfun(@(x) diag(x,-1),act_corr,'uni',false);

% Plot correlation aligned to n reversal across rev animals
max_rev = max(cellfun(@length,reversals(rev_animals)));
sq_plot = ceil(sqrt(max_rev));
figure;
for curr_rev = 1:max_rev
    subplot(sq_plot,sq_plot,curr_rev);hold on;
    for curr_animal = find(rev_animals);
        if pmm_animals(curr_animal)
            col = 'k';
        else
            col = 'r';
        end
        if length(reversals{curr_animal}) >= curr_rev
            plot([1:size(act_corr_oneback{curr_animal})]- ...
                reversals{curr_animal}(curr_rev),act_corr_oneback{curr_animal},col);
        end
    end
    title(['Reversal ' num2str(curr_rev)]);
end

% Plot average correlation aligned to n reversal across rev animals
max_rev = max(cellfun(@length,reversals(rev_animals)));

sq_plot = ceil(sqrt(max_rev));
figure;
for curr_rev = 1:max_rev
    
    subplot(sq_plot,sq_plot,curr_rev); hold on;
    
    curr_epochs_pmm = [];
    curr_corr_pmm = [];
    
    curr_epochs_alm = [];
    curr_corr_alm = [];
    
    for curr_animal = find(rev_animals);
        if length(reversals{curr_animal}) >= curr_rev
            if pmm_animals(curr_animal)
                curr_epochs_pmm = [curr_epochs_pmm; ([1:size(act_corr_oneback{curr_animal})]'- ...
                    reversals{curr_animal}(curr_rev))];
                curr_corr_pmm = [curr_corr_pmm;act_corr_oneback{curr_animal}];
            else              
                curr_epochs_alm = [curr_epochs_alm; ([1:size(act_corr_oneback{curr_animal})]'- ...
                    reversals{curr_animal}(curr_rev))];
                curr_corr_alm = [curr_corr_alm;act_corr_oneback{curr_animal}];
            end
        end
    end   
    
    mean_corr_pmm = grpstats(curr_corr_pmm,curr_epochs_pmm,'mean');
    plot_epochs_pmm = unique(curr_epochs_pmm);
    plot(plot_epochs_pmm,mean_corr_pmm,'color','k','linewidth',2);
    
    mean_corr_alm = grpstats(curr_corr_alm,curr_epochs_alm,'mean');
    plot_epochs_alm = unique(curr_epochs_alm);
    plot(plot_epochs_alm,mean_corr_alm,'color','r','linewidth',2);
    
    line([0 0],ylim,'color','k','linestyle','--');
    
    title(['Reversal ' num2str(curr_rev)]);
end


% Plot grid aligned to reversal
max_rev = max(cellfun(@length,reversals(rev_animals)));

max_epochs = max(cellfun(@length,act_corr));
act_corr_pad_forward = cellfun(@(x) padarray(x,[max_epochs-length(x) ...
    max_epochs-length(x)],NaN,'post'),act_corr,'uni',false);
act_corr_pad = cellfun(@(x) padarray(x,[max_epochs ...
    max_epochs],NaN,'pre'),act_corr_pad_forward,'uni',false);

pmm_grid_fig = figure;
set(pmm_grid_fig,'Name','PMM');
alm_grid_fig = figure;
set(alm_grid_fig,'Name','ALM');
sq_plot = ceil(sqrt(max_rev));

for curr_rev = 1:5
    % shift all correlation grids to align by reversal
    use_pmm = rev_animals & pmm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_pmm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_pmm),reversals(use_pmm),'uni',false),[1 3 2])),3);
    
    figure(pmm_grid_fig);
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(1,5,curr_rev);
    imagesc(curr_pmm_grid_avg);
    colormap(gray);
    line(xlim,[max_epochs-1.5 max_epochs-1.5],'color','r');
    line([max_epochs-1.5 max_epochs-1.5],ylim,'color','r');
    title(['Reversal ' num2str(curr_rev)]);
    
    use_alm = rev_animals & alm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_alm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_alm),reversals(use_alm),'uni',false),[1 3 2])),3);
    
    figure(alm_grid_fig);
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(1,5,curr_rev);
    imagesc(curr_alm_grid_avg);
    colormap(gray);
    line(xlim,[max_epochs-1.5 max_epochs-1.5],'color','r');
    line([max_epochs-1.5 max_epochs-1.5],ylim,'color','r');
    title(['Reversal ' num2str(curr_rev)]);
end

%% TK requested: like last cell, but trial-trial correlations

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% % all trials
% odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
%     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
%     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);

% % d' cutoff
% odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
%     reshape(permute(x(z,y > 0 & y < 2,:),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
%     x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);

% d' cutoff binary activity
odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
    permute(any(x(z,y > 0 & y < 2,:),2),[3 1 2]),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);

odor_activity_CL_cat = cellfun(@(x) vertcat(x{:})',odor_activity_CL,'uni',false);

epoch_sessions_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_sessions,'uni',false);
epoch_contingency_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_contingencies,'uni',false);

% There are random occurrences of the wrong contingency where there was not
% meant to be a reversal. This means that some epochs have very few trials
% but are only identifiable by this fact. Identify where these epochs occur
% and pull them out indivudally.
% --- but this happens a lot? and there isn't a clear cutoff ??? this makes
% this impossible.
% THE ONLY THING I CAN DO: if there under an arbitrary cutoff, then ignore
% it entirely. I'm making the cutoff now based what looks like a line
% separating number of CL trials.

num_cutoff_CL_trials = cellfun(@(x) cellfun(@(x) sum(x),vertcat(x{:})), ...
    d_prime_use_trials_CL_epoch','uni',false);
num_trial_cutoff = 20;
use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_cutoff_CL_trials,'uni',false);

odor_activity_CL_cat_use = cellfun(@(x,y) x(:,y),odor_activity_CL_cat, ...
    use_epochs,'uni',false);
epoch_sessions_cat_use = cellfun(@(x,y) x(y),epoch_sessions_cat, ...
    use_epochs,'uni',false);
epoch_contingency_cat_use = cellfun(@(x,y) x(y),epoch_contingency_cat, ...
    use_epochs,'uni',false);

% Get reversals (after selecting out epochs)
reversals = cellfun(@(x) find(diff(x) ~= 0),epoch_contingency_cat_use,'uni',false);

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Get day-before TRIAL-TRIAL correlation
act_corr_trial = cellfun(@(x) mat2cell(corrcoef(horzcat(x{:}),'rows','complete'), ...
    cellfun(@(x) size(x,2),x),cellfun(@(x) size(x,2),x)), ...
    odor_activity_CL_cat_use,'uni',false);
act_corr = cellfun(@(x) cellfun(@(x) nanmedian(x(:)),x),act_corr_trial,'uni',false);
act_corr_oneback = cellfun(@(x) diag(x,-1),act_corr,'uni',false);

% Plot correlation aligned to n reversal across rev animals
max_rev = max(cellfun(@length,reversals(rev_animals)));
sq_plot = ceil(sqrt(max_rev));
figure;
for curr_rev = 1:max_rev
    subplot(sq_plot,sq_plot,curr_rev);hold on;
    for curr_animal = find(rev_animals);
        if pmm_animals(curr_animal)
            col = 'k';
        else
            col = 'r';
        end
        if length(reversals{curr_animal}) >= curr_rev
            plot([1:size(act_corr_oneback{curr_animal})]- ...
                reversals{curr_animal}(curr_rev),act_corr_oneback{curr_animal},col);
        end
    end
    title(['Reversal ' num2str(curr_rev)]);
end

% Plot average correlation aligned to n reversal across rev animals
max_rev = max(cellfun(@length,reversals(rev_animals)));

sq_plot = ceil(sqrt(max_rev));
figure;
for curr_rev = 1:max_rev
    
    subplot(sq_plot,sq_plot,curr_rev); hold on;
    
    curr_epochs_pmm = [];
    curr_corr_pmm = [];
    
    curr_epochs_alm = [];
    curr_corr_alm = [];
    
    for curr_animal = find(rev_animals);
        if length(reversals{curr_animal}) >= curr_rev
            if pmm_animals(curr_animal)
                % Only plot within current and last reversal
                curr_corr_oneback = act_corr_oneback{curr_animal};
                if length(reversals{curr_animal}) > curr_rev
                    curr_corr_oneback(reversals{curr_animal}(curr_rev+1):end) = NaN;
                end
                if curr_rev >= 2
                    curr_corr_oneback(1:reversals{curr_animal}(curr_rev-1)) = NaN;
                end
                
                curr_epochs_pmm = [curr_epochs_pmm; ([1:size(act_corr_oneback{curr_animal})]'- ...
                    reversals{curr_animal}(curr_rev))];
                curr_corr_pmm = [curr_corr_pmm;curr_corr_oneback];
            else       
                % Only plot within current and last reversal
                curr_corr_oneback = act_corr_oneback{curr_animal};
                if length(reversals{curr_animal}) > curr_rev
                    curr_corr_oneback(reversals{curr_animal}(curr_rev+1):end) = NaN;
                end
                if curr_rev >= 2
                    curr_corr_oneback(1:reversals{curr_animal}(curr_rev-1)) = NaN;
                end
                
                curr_epochs_alm = [curr_epochs_alm; ([1:size(act_corr_oneback{curr_animal})]'- ...
                    reversals{curr_animal}(curr_rev))];
                curr_corr_alm = [curr_corr_alm;curr_corr_oneback];
            end
        end
    end   
    
    mean_corr_pmm = grpstats(curr_corr_pmm,curr_epochs_pmm,'mean');
    plot_epochs_pmm = unique(curr_epochs_pmm);
    plot(plot_epochs_pmm,mean_corr_pmm,'color','k','linewidth',2);
    
    mean_corr_alm = grpstats(curr_corr_alm,curr_epochs_alm,'mean');
    plot_epochs_alm = unique(curr_epochs_alm);
    plot(plot_epochs_alm,mean_corr_alm,'color','r','linewidth',2);
    
    line([0 0],ylim,'color','k','linestyle','--');
    
    title(['Reversal ' num2str(curr_rev)]);
end

% Plot average control animals
curr_epochs_pmm = [];
curr_corr_pmm = [];

curr_epochs_alm = [];
curr_corr_alm = [];
for curr_animal = find(ctrl_animals);
    if pmm_animals(curr_animal)
        % Don't include 'control' reversals
        curr_corr_oneback = act_corr_oneback{curr_animal};
        if length(reversals{curr_animal}) > 0
            curr_corr_oneback(reversals{curr_animal}(1):end) = NaN;
        end
        curr_epochs_pmm = [curr_epochs_pmm; ([1:size(act_corr_oneback{curr_animal})]')];
        curr_corr_pmm = [curr_corr_pmm;curr_corr_oneback];
    else
        % Don't include 'control' reversals
        curr_corr_oneback = act_corr_oneback{curr_animal};
        if length(reversals{curr_animal}) > 0
            curr_corr_oneback(reversals{curr_animal}(1):end) = NaN;
        end
        curr_epochs_alm = [curr_epochs_alm; ([1:size(act_corr_oneback{curr_animal})]')];
        curr_corr_alm = [curr_corr_alm;curr_corr_oneback];
    end
end
figure; hold on;
mean_corr_pmm = grpstats(curr_corr_pmm,curr_epochs_pmm,'mean');
plot_epochs_pmm = unique(curr_epochs_pmm);
plot(plot_epochs_pmm,mean_corr_pmm,'color','k','linewidth',2,'linestyle','--');

mean_corr_alm = grpstats(curr_corr_alm,curr_epochs_alm,'mean');
plot_epochs_alm = unique(curr_epochs_alm);
plot(plot_epochs_alm,mean_corr_alm,'color','r','linewidth',2,'linestyle','--');

line([0 0],ylim,'color','k','linestyle','--');

title(['Control']);


% Plot correlation relative to last epoch before reversal
figure; hold on
% Rev
curr_epochs_pmm = [];
curr_corr_pmm = [];

curr_epochs_alm = [];
curr_corr_alm = [];
for curr_animal = find(rev_animals);
    if pmm_animals(curr_animal)
        curr_corr_prerev = act_corr{curr_animal}(:,reversals{curr_animal}(1));
        if length(reversals{curr_animal}) > 0
            curr_corr_prerev(reversals{curr_animal}(2):end) = NaN;
        end
        curr_epochs_pmm = [curr_epochs_pmm; ([1:length(curr_corr_prerev)]'-(reversals{curr_animal}(1)))];
        curr_corr_pmm = [curr_corr_pmm;curr_corr_prerev];
    else
        curr_corr_prerev = act_corr{curr_animal}(:,reversals{curr_animal}(1));
        if length(reversals{curr_animal}) > 0
            curr_corr_prerev(reversals{curr_animal}(2):end) = NaN;
        end
        curr_epochs_alm = [curr_epochs_alm; ([1:length(curr_corr_prerev)]'-(reversals{curr_animal}(1)))];
        curr_corr_alm = [curr_corr_alm;curr_corr_prerev];
    end
end
mean_corr_pmm = grpstats(curr_corr_pmm,curr_epochs_pmm,'mean');
plot_epochs_pmm = unique(curr_epochs_pmm);
plot(plot_epochs_pmm,mean_corr_pmm,'color','k','linewidth',2);

mean_corr_alm = grpstats(curr_corr_alm,curr_epochs_alm,'mean');
plot_epochs_alm = unique(curr_epochs_alm);
plot(plot_epochs_alm,mean_corr_alm,'color','r','linewidth',2);

% Ctrl
curr_epochs_pmm = [];
curr_corr_pmm = [];

curr_epochs_alm = [];
curr_corr_alm = [];
use_session = 6;
for curr_animal = find(ctrl_animals);
    if pmm_animals(curr_animal)
        curr_corr_prerev = act_corr{curr_animal}(:,use_session);
        if length(reversals{curr_animal}) > 0
            curr_corr_prerev(reversals{curr_animal}(1):end) = NaN;
        end
        curr_epochs_pmm = [curr_epochs_pmm; ([1:length(curr_corr_prerev)]'-use_session)];
        curr_corr_pmm = [curr_corr_pmm;curr_corr_prerev];
    else
        curr_corr_prerev = act_corr{curr_animal}(:,use_session);
        if length(reversals{curr_animal}) > 0
            curr_corr_prerev(reversals{curr_animal}(1):end) = NaN;
        end
        curr_epochs_alm = [curr_epochs_alm; ([1:length(curr_corr_prerev)]'-use_session)];
        curr_corr_alm = [curr_corr_alm;curr_corr_prerev];
    end
end
mean_corr_pmm = grpstats(curr_corr_pmm,curr_epochs_pmm,'mean');
plot_epochs_pmm = unique(curr_epochs_pmm);
plot(plot_epochs_pmm,mean_corr_pmm,'color','k','linewidth',2,'linestyle','--');

mean_corr_alm = grpstats(curr_corr_alm,curr_epochs_alm,'mean');
plot_epochs_alm = unique(curr_epochs_alm);
plot(plot_epochs_alm,mean_corr_alm,'color','r','linewidth',2,'linestyle','--');


line([0 0],ylim,'color','k','linestyle','--');

title(['Trial-trial correlation with pre-reversal epoch']);


% Correlation prerev/postrev/across
pre_median = nan(size(animals));
post_median = nan(size(animals));
cross_median = nan(size(animals));
for curr_animal = 1:length(animals)
    if ismember(find(rev_animals),curr_animal)
        curr_pre_rev = 1:reversals{curr_animal}(1);
        curr_post_rev = reversals{curr_animal}(1):reversals{curr_animal}(2);
    else
        curr_pre_rev = 1:6;
        curr_post_rev = 6:10;
    end
    
    curr_pre_corr = act_corr_trial{curr_animal}(curr_pre_rev,curr_pre_rev);
    pre_use_epochs = logical(tril(ones(size(curr_pre_corr)),-1));
    pre_corr_trials = cell2mat(cellfun(@(x) x(:),curr_pre_corr(pre_use_epochs),'uni',false));
    pre_median(curr_animal) = nanmedian(pre_corr_trials);
    
    curr_post_corr = act_corr_trial{curr_animal}(curr_post_rev,curr_post_rev);
    post_use_epochs = logical(tril(ones(size(curr_post_corr)),-1));
    post_corr_trials = cell2mat(cellfun(@(x) x(:),curr_post_corr(post_use_epochs),'uni',false));
    post_median(curr_animal) = nanmedian(post_corr_trials);
    
    curr_cross_corr = act_corr_trial{curr_animal}(curr_pre_rev,curr_post_rev);
    cross_corr_trials = cell2mat(cellfun(@(x) x(:),curr_cross_corr(:),'uni',false));
    cross_median(curr_animal) = nanmedian(cross_corr_trials);
end
figure;
bar([...
    nanmean(pre_median(pmm_animals & rev_animals)) ...
    nanmean(pre_median(pmm_animals & ctrl_animals)) ...
    nanmean(pre_median(alm_animals & rev_animals)) ...
    nanmean(pre_median(alm_animals & ctrl_animals)); ...
    nanmean(post_median(pmm_animals & rev_animals)) ...
    nanmean(post_median(pmm_animals & ctrl_animals)) ...
    nanmean(post_median(alm_animals & rev_animals)) ...
    nanmean(post_median(alm_animals & ctrl_animals)); ...
    nanmean(cross_median(pmm_animals & rev_animals)) ...
    nanmean(cross_median(pmm_animals & ctrl_animals)) ...
    nanmean(cross_median(alm_animals & rev_animals)) ...
    nanmean(cross_median(alm_animals & ctrl_animals))]);
legend({'PMM rev' 'PMM ctrl' 'ALM rev' 'ALM ctrl'})
set(gca,'XTickLabel',{'Pre-rev' 'Post-rev' 'Cross'})
ylabel('Median trial-trial correlation')


% Plot grid aligned to reversal
max_rev = max(cellfun(@length,reversals(rev_animals)));

max_epochs = max(cellfun(@length,act_corr));
act_corr_pad_forward = cellfun(@(x) padarray(x,[max_epochs-length(x) ...
    max_epochs-length(x)],NaN,'post'),act_corr,'uni',false);
act_corr_pad = cellfun(@(x) padarray(x,[max_epochs ...
    max_epochs],NaN,'pre'),act_corr_pad_forward,'uni',false);

pmm_grid_fig = figure;
set(pmm_grid_fig,'Name','PMM');
alm_grid_fig = figure;
set(alm_grid_fig,'Name','ALM');
sq_plot = ceil(sqrt(max_rev));

for curr_rev = 1:5
    % shift all correlation grids to align by reversal
    use_pmm = rev_animals & pmm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_pmm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_pmm),reversals(use_pmm),'uni',false),[1 3 2])),3);
    
    figure(pmm_grid_fig);
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(1,5,curr_rev);
    imagesc(curr_pmm_grid_avg);
    colormap(gray);
    line(xlim,[max_epochs-1.5 max_epochs-1.5],'color','r');
    line([max_epochs-1.5 max_epochs-1.5],ylim,'color','r');
    title(['Reversal ' num2str(curr_rev)]);
    
    use_alm = rev_animals & alm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_alm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_alm),reversals(use_alm),'uni',false),[1 3 2])),3);
    
    figure(alm_grid_fig);
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(1,5,curr_rev);
    imagesc(curr_alm_grid_avg);
    colormap(gray);
    line(xlim,[max_epochs-1.5 max_epochs-1.5],'color','r');
    line([max_epochs-1.5 max_epochs-1.5],ylim,'color','r');
    title(['Reversal ' num2str(curr_rev)]);
end

% Plot control grid
figure;
subplot(2,1,1)
curr_animals = find(ctrl_animals & pmm_animals);
curr_grid = cat(3,act_corr_pad_forward{curr_animals});
for i = 1:length(curr_animals)
    curr_animal = curr_animals(i);
    if length(reversals{curr_animal}) > 0
        curr_grid(reversals{curr_animal}(1):end,:,i) = NaN;
        curr_grid(:,reversals{curr_animal}(1):end,i) = NaN;
    end
end
imagesc(nanmean(curr_grid,3));colormap(gray);
title('PMM control')
subplot(2,1,2)
curr_animals = find(ctrl_animals & alm_animals);
curr_grid = cat(3,act_corr_pad_forward{curr_animals});
for i = 1:length(curr_animals)
    curr_animal = curr_animals(i);
    if length(reversals{curr_animal}) > 0
        curr_grid(reversals{curr_animal}(1):end,:,i) = NaN;
        curr_grid(:,reversals{curr_animal}(1):end,i) = NaN;
    end
end
imagesc(nanmean(curr_grid,3));colormap(gray);
title('ALM control')


%% Contingency selectivity by cell

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% % all trials
% odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
%     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
%     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);

% d' cutoff
odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
    permute(any(x(z,y > 0 & y < 2,:),2),[1 3 2]),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);


odor_activity_CL_cat = cellfun(@(x) vertcat(x{:})',odor_activity_CL,'uni',false);

epoch_sessions_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_sessions,'uni',false);
epoch_contingency_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_contingencies,'uni',false);

% There are random occurrences of the wrong contingency where there was not
% meant to be a reversal. This means that some epochs have very few trials
% but are only identifiable by this fact. Identify where these epochs occur
% and pull them out indivudally.
% --- but this happens a lot? and there isn't a clear cutoff ??? this makes
% this impossible.
% THE ONLY THING I CAN DO: if there under an arbitrary cutoff, then ignore
% it entirely. I'm making the cutoff now based what looks like a line
% separating number of CL trials.

num_cutoff_CL_trials = cellfun(@(x) cellfun(@(x) cellfun(@(x) sum(x),x),x,'uni',false), ...
    d_prime_use_trials_CL_epoch','uni',false);
trial_cutoff = 20;
cont_selectivity = cell(length(animals),1);

for curr_animal = 1:length(animals)
    curr_contingencies = cellfun(@(x,y) x(y >= trial_cutoff),analysis.epoch_contingencies{curr_animal}, ...
        num_cutoff_CL_trials{curr_animal},'uni',false);
    
    curr_epochs = cellfun(@(x,y) x(y >= trial_cutoff),odor_activity_CL{curr_animal}, ...
        num_cutoff_CL_trials{curr_animal},'uni',false);
    rev_sessions = find(cellfun(@(x) length(x) > 1,curr_epochs));
    
    curr_selectivity = [];
    for curr_session_idx = 1:length(rev_sessions)
       curr_session = rev_sessions(curr_session_idx);
       cont1_reliability = vertcat(curr_epochs{curr_session}{curr_contingencies{curr_session} == 1});
       cont2_reliability = vertcat(curr_epochs{curr_session}{curr_contingencies{curr_session} == 2});
       
       curr_selectivity = nanmean(cont1_reliability,1) - nanmean(cont2_reliability,1);
       cont_selectivity{curr_animal}(:,curr_session_idx) = curr_selectivity';    
    end    
end

figure('Name','ALM');
use_animals = find(alm_animals & rev_animals);
for i = 1:length(use_animals);
   subplot(1,length(use_animals),i);
   [~,sort_cells] = sort(nanmean(cont_selectivity{use_animals(i)},2));
   imagesc(cont_selectivity{use_animals(i)}(sort_cells,:));
   colormap(gray);
end
figure('Name','PMM');
use_animals = find(pmm_animals & rev_animals);
for i = 1:length(use_animals);
    if isempty(cont_selectivity{use_animals(i)})
        continue
    end
   subplot(1,length(use_animals),i);
   [~,sort_cells] = sort(nanmean(cont_selectivity{use_animals(i)},2));
   imagesc(cont_selectivity{use_animals(i)}(sort_cells,:));
   colormap(gray);
end




%% Get active cells for quality control

test = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
    permute(sum(any(x(:,y > 0 & y < 2,:),2),1),[3 2 1]),x,'uni',false), ...
    x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);


test2 = cellfun(@(x) cellfun(@(x) sum(horzcat(x{:}),2),x,'uni',false), ...
    test,'uni',false);


test3 = cellfun(@(x) horzcat(x{:}),test2,'uni',false);

min_trials = 5;
active_cells = cellfun(@(x) find(any(x >= 5,2)),test3,'uni',false);
active_cell_sessions = cellfun(@(x) arrayfun(@(y) find(x(y,:) >= 5), ...
    1:size(x,1),'uni',false),test3,'uni',false);























