%% Grab reversal epoch activity (within reversals, within day if none)

% split aligned activity into pyr/gad
aligned_df_gad = cellfun(@(x,y) cellfun(@(z) z(:,:,mice(y).interneurons),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);
aligned_df_pyr = cellfun(@(x,y) cellfun(@(z) z(:,:,setdiff(1:size(z,3),mice(y).interneurons)),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);

% Get reversal trials (reversal = first trial of new contingency)
all_reversals = cellfun(@(x) cellfun(@(y) find(diff(y(:,2)) ~= 0)+1,x,'uni',false), ...
   odor_trials_all,'uni',false); 

% Get the start/ends of epoch (either within reversal or entire day if not)
epoch_boundaries = cellfun(@(rev,odor) cellfun(@(rev,odor) ...
    mat2cell([1;reshape([rev-1 rev]',[],1);size(odor,1)], ...
    repmat(2,length(rev)+1,1)),rev,odor,'uni',false), ...
    all_reversals,odor_trials_all,'uni',false);

% Get the contingencies and days corresponding to reversals
epoch_contingencies = cellfun(@(epoch,odor) cellfun(@(epoch,odor) cellfun(@(epoch)...
    odor(epoch(1),2),epoch),epoch,odor,'uni',false), ...
    epoch_boundaries,odor_trials_all,'uni',false);

epoch_sessions = cellfun(@(epoch) cellfun(@(epoch,sessions) ...
    repmat(sessions,length(epoch),1),epoch,num2cell(1:length(epoch)), ...
    'uni',false),epoch_boundaries,'uni',false);

% Pull out the activity associated with epochs
epoch_activity_CL = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,1)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),aligned_df_pyr, ...
    condition_trials_all,epoch_boundaries,'uni',false);

epoch_activity_CR = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,2)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),aligned_df_pyr, ...
    condition_trials_all,epoch_boundaries,'uni',false);

% Concatenate all trials, associate with contingentices / sessions
epoch_activity_CL_cat = cellfun(@(df) cell2mat(cellfun(@(df) ...
    cat(1,df{:}),df,'uni',false)'),epoch_activity_CL,'uni',false);

contingency_CL_cat = cellfun(@(bhv,odor) cell2mat(cellfun(@(bhv,odor) odor(bhv(:,1),2), ...
    bhv,odor,'uni',false)'), condition_trials_all,odor_trials_all,'uni',false);

sessions_CL_cat = cellfun(@(bhv) cell2mat(cellfun(@(bhv,sessions) ...
    repmat(sessions,sum(bhv(:,1)),1),bhv,num2cell(1:length(bhv)), ...
    'uni',false)'),condition_trials_all,'uni',false);


epoch_activity_CR_cat = cellfun(@(df) cell2mat(cellfun(@(df) ...
    cat(1,df{:}),df,'uni',false)'),epoch_activity_CR,'uni',false);

contingency_CR_cat = cellfun(@(bhv,odor) cell2mat(cellfun(@(bhv,odor) odor(bhv(:,2),2), ...
    bhv,odor,'uni',false)'), condition_trials_all,odor_trials_all,'uni',false);

sessions_CR_cat = cellfun(@(bhv) cell2mat(cellfun(@(bhv,sessions) ...
    repmat(sessions,sum(bhv(:,2)),1),bhv,num2cell(1:length(bhv)), ...
    'uni',false)'),condition_trials_all,'uni',false);


%% CATEGORY: Population correlation across reversals

%% 1) Same-contingency population correlation: fraction binary active trials

epoch_activity_CL_binary = cellfun(@(df) ...
    permute(any(df,2),[1 3 2]), ...
    epoch_activity_CL_cat,'uni',false);

rev = cellfun(@(contingency) find(diff(contingency) ~= 0)+1, ...
    contingency_CL_cat,'uni',false);
reversal_animals = cellfun(@(x) ~isempty(x),rev);

rev_diff = cellfun(@(rev,contingency) [rev(1);diff(rev); ...
    length(contingency)-rev(end)],rev(reversal_animals), ...
    contingency_CL_cat(reversal_animals),'uni',false);

split_activity_fracbinary_CL_corr = cellfun(@(activity,rev_diff) ...
    corrcoef(cell2mat(cellfun(@(x) nanmean(x,1),mat2cell(activity,rev_diff,size(activity,2)), ...
    'uni',false))'),epoch_activity_CL_binary(reversal_animals), ...
    rev_diff,'uni',false);

contingency_1_meancorr = cellfun(@(act_corr) nanmean(AP_itril( ...
    act_corr(1:2:end,1:2:end),-1)),split_activity_fracbinary_CL_corr);

contingency_2_meancorr = cellfun(@(act_corr) nanmean(AP_itril( ...
    act_corr(2:2:end,2:2:end),-1)),split_activity_fracbinary_CL_corr);

contingency_1_2_meancorr = cellfun(@(act_corr) nanmean(reshape( ...
    act_corr(1:2:end,2:2:end),[],1)),split_activity_fracbinary_CL_corr);

figure; hold on
bar([nanmean(contingency_1_meancorr) nanmean(contingency_2_meancorr) ...
    nanmean(contingency_1_2_meancorr)],'FaceColor','none','linewidth',2,'BarWidth',0.5)
errorbar([nanmean(contingency_1_meancorr) nanmean(contingency_2_meancorr) ...
    nanmean(contingency_1_2_meancorr)], ...
    [nanstd(contingency_1_meancorr)./sqrt(sum(~isnan(contingency_1_meancorr))) ...
    nanstd(contingency_2_meancorr)./sqrt(sum(~isnan(contingency_2_meancorr))) ...
    nanstd(contingency_1_2_meancorr)./sqrt(sum(~isnan(contingency_1_2_meancorr)))],'.k','linewidth',2);
plot(1,contingency_1_meancorr,'.r');
plot(2,contingency_2_meancorr,'.r');
plot(3,contingency_1_2_meancorr,'.r');
xlim([0 4])
xlabel('Contingency')
ylabel('Mean correlation')
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'1-1' '2-2' '1-2'})
title('Population = mean binary activity within contingency (across days). ANOVA p = 0.83')


%% 2) Correlations of active trials across sessions, reversals, random /control

rev_animals = [1 2 3 6];
ctrl_animals = [4 5];

min_trials = 10; % pre-post trial count

epoch_activity_CL_binary = cellfun(@(df) ...
    permute(any(df,2),[1 3 2]), ...
    epoch_activity_CL_cat,'uni',false);

rev = cellfun(@(contingency) find(diff(contingency) ~= 0)+1, ...
    contingency_CL_cat,'uni',false);

ses = cellfun(@(sessions) find(diff(sessions) ~= 0)+1, ...
    sessions_CL_cat,'uni',false);

% Correlation across reversals
event = rev;
rev_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    rev_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1) || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event)+1:event{curr_animal}(curr_event)+min_trials], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        rev_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
rev_event_corr_mean = cellfun(@nanmean,rev_event_corr);

% Correlation across sessions
event = ses;
ses_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    ses_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1) || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event)+1:event{curr_animal}(curr_event)+min_trials], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        ses_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
ses_event_corr_mean = cellfun(@nanmean,ses_event_corr);

% Correlation across random points (# as max(rev,ses))
event = cellfun(@(x,rev,ses) randi(size(x,1), ...
    10000*max(length(ses),length(rev)),1), ...
    epoch_activity_CL_binary,ses,rev,'uni',false);
rand_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    rand_event_corr_total = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1)          
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        rand_event_corr_total(curr_event) = curr_corr(2);
    end
    rand_event_corr{curr_animal} = prctile(cellfun(@nanmean,mat2cell(rand_event_corr_total, ...
        repmat(length(rand_event_corr_total)/10000,10000,1))),[2.5 97.5]);
end

rand_event_corr_conf_rev = nanmean(cell2mat(rand_event_corr(rev_animals)'));
rand_event_corr_conf_ctrl = nanmean(cell2mat(rand_event_corr(ctrl_animals)'));

figure; hold on;
errorbar([nanmean(rev_event_corr_mean(rev_animals)) ...
    nanmean(ses_event_corr_mean(rev_animals))], ...
    [nanstd(rev_event_corr_mean(rev_animals))./sqrt(sum(~isnan(rev_event_corr_mean(rev_animals)))) ...
    nanstd(ses_event_corr_mean(rev_animals))./sqrt(sum(~isnan(ses_event_corr_mean(rev_animals))))],'.k','linewidth',2);
line([1 2],repmat(rand_event_corr_conf_rev(1),1,2),'color','r','linestyle','--');
line([1 2],repmat(rand_event_corr_conf_rev(2),1,2),'color','r','linestyle','--');

errorbar(3,nanmean(ses_event_corr_mean(ctrl_animals)), ...
       nanstd(ses_event_corr_mean(ctrl_animals))./sqrt(sum(~isnan(ses_event_corr_mean(ctrl_animals)))),'.k','linewidth',2);
line([2 3],repmat(rand_event_corr_conf_ctrl(1),1,2),'color','r','linestyle','--');
line([2 3],repmat(rand_event_corr_conf_ctrl(2),1,2),'color','r','linestyle','--');
   
xlim([0 4])
ylabel('Mean correlation')
xlabel('Event type')
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Reversal','Session','Session (Ctrl)'});
title('Population correlation (fraction binary active / 30 trials) across events')

%% (Same as above but CR)
rev_animals = [1 2 3 6];
ctrl_animals = [4 5];

min_trials = 30; % pre-post trial count

epoch_activity_CR_binary = cellfun(@(df) ...
    permute(any(df,2),[1 3 2]), ...
    epoch_activity_CR_cat,'uni',false);

rev = cellfun(@(contingency) find(diff(contingency) ~= 0)+1, ...
    contingency_CR_cat,'uni',false);

ses = cellfun(@(sessions) find(diff(sessions) ~= 0)+1, ...
    sessions_CR_cat,'uni',false);

% Correlation across reversals
event = rev;
rev_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CR_binary)
    activity = epoch_activity_CR_binary{curr_animal};
    rev_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1) || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event)+1:event{curr_animal}(curr_event)+min_trials], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        rev_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
rev_event_corr_mean = cellfun(@nanmean,rev_event_corr);

% Correlation across sessions
event = ses;
ses_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CR_binary)
    activity = epoch_activity_CR_binary{curr_animal};
    ses_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1) || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event)+1:event{curr_animal}(curr_event)+min_trials], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        ses_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
ses_event_corr_mean = cellfun(@nanmean,ses_event_corr);

% Correlation across random points (# as max(rev,ses))
event = cellfun(@(x,rev,ses) randi(size(x,1), ...
    10000*max(length(ses),length(rev)),1), ...
    epoch_activity_CR_binary,ses,rev,'uni',false);
rand_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CR_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    rand_event_corr_total = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1)          
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:)));
        
        rand_event_corr_total(curr_event) = curr_corr(2);
    end
    rand_event_corr{curr_animal} = prctile(cellfun(@nanmean,mat2cell(rand_event_corr_total, ...
        repmat(length(rand_event_corr_total)/10000,10000,1))),[2.5 97.5]);
end

rand_event_corr_conf_rev = nanmean(cell2mat(rand_event_corr(rev_animals)'));
rand_event_corr_conf_ctrl = nanmean(cell2mat(rand_event_corr(ctrl_animals)'));

figure; hold on;
errorbar([nanmean(rev_event_corr_mean(rev_animals)) ...
    nanmean(ses_event_corr_mean(rev_animals))], ...
    [nanstd(rev_event_corr_mean(rev_animals))./sqrt(sum(~isnan(rev_event_corr_mean(rev_animals)))) ...
    nanstd(ses_event_corr_mean(rev_animals))./sqrt(sum(~isnan(ses_event_corr_mean(rev_animals))))],'.k','linewidth',2);
line([1 2],repmat(rand_event_corr_conf_rev(1),1,2),'color','r','linestyle','--');
line([1 2],repmat(rand_event_corr_conf_rev(2),1,2),'color','r','linestyle','--');

errorbar(3,nanmean(ses_event_corr_mean(ctrl_animals)), ...
       nanstd(ses_event_corr_mean(ctrl_animals))./sqrt(sum(~isnan(ses_event_corr_mean(ctrl_animals)))),'.k','linewidth',2);
line([2 3],repmat(rand_event_corr_conf_ctrl(1),1,2),'color','r','linestyle','--');
line([2 3],repmat(rand_event_corr_conf_ctrl(2),1,2),'color','r','linestyle','--');
   
xlim([0 4])
ylabel('Mean correlation')
xlabel('Event type')
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Reversal','Session','Session (Ctrl)'});
title('Population correlation (fraction binary active / 30 trials) across events')


%% (Like above but CL, and pre vs. next pre (if similar once things settle down))

rev_animals = [1 2 3 6];
ctrl_animals = [4 5];

min_trials = 30; % pre trial count

epoch_activity_CL_binary = cellfun(@(df) ...
    permute(any(df,2),[1 3 2]), ...
    epoch_activity_CL_cat,'uni',false);

rev = cellfun(@(contingency) find(diff(contingency) ~= 0)+1, ...
    contingency_CL_cat,'uni',false);

ses = cellfun(@(sessions) find(diff(sessions) ~= 0)+1, ...
    sessions_CL_cat,'uni',false);

% Correlation across reversals
event = rev;
rev_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    rev_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})-1
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event+1)-min_trials < 1 || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event+1)-min_trials:event{curr_animal}(curr_event+1)-1], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event+1)-min_trials: ...
            event{curr_animal}(curr_event+1)-1,:)));
        
        rev_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
rev_event_corr_mean = cellfun(@nanmean,rev_event_corr);

% Correlation across sessions
event = ses;
ses_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    ses_event_corr{curr_animal} = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})-1
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event+1)-min_trials < 1 || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event+1)-min_trials:event{curr_animal}(curr_event+1)-1], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event+1)-min_trials: ...
            event{curr_animal}(curr_event+1)-1,:)));
        
        ses_event_corr{curr_animal}(curr_event) = curr_corr(2);
    end
end
ses_event_corr_mean = cellfun(@nanmean,ses_event_corr);

% Correlation across random points (# as max(rev,ses))
num_rand = 1000;
event = cellfun(@(x,rev,ses) randi(size(x,1), ...
    num_rand*max(length(ses),length(rev)),1), ...
    epoch_activity_CL_binary,ses,rev,'uni',false);
rand_event_corr = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    activity = epoch_activity_CL_binary{curr_animal};
    rand_event_corr_total = nan(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})-1
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event+1)-min_trials < 1  
            continue
        end
        curr_corr = corrcoef(nanmean(activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:)), ...
            nanmean(activity(event{curr_animal}(curr_event+1)-min_trials: ...
            event{curr_animal}(curr_event+1)-1,:)));
        
        rand_event_corr_total(curr_event) = curr_corr(2);
    end
    rand_event_corr{curr_animal} = prctile(cellfun(@nanmean,mat2cell(rand_event_corr_total, ...
        repmat(length(rand_event_corr_total)/num_rand,num_rand,1))),[2.5 97.5]);
end

rand_event_corr_conf_rev = nanmean(cell2mat(rand_event_corr(rev_animals)'));
rand_event_corr_conf_ctrl = nanmean(cell2mat(rand_event_corr(ctrl_animals)'));

figure; hold on;
errorbar([nanmean(rev_event_corr_mean(rev_animals)) ...
    nanmean(ses_event_corr_mean(rev_animals))], ...
    [nanstd(rev_event_corr_mean(rev_animals))./sqrt(sum(~isnan(rev_event_corr_mean(rev_animals)))) ...
    nanstd(ses_event_corr_mean(rev_animals))./sqrt(sum(~isnan(ses_event_corr_mean(rev_animals))))],'.k','linewidth',2);
line([1 2],repmat(rand_event_corr_conf_rev(1),1,2),'color','r','linestyle','--');
line([1 2],repmat(rand_event_corr_conf_rev(2),1,2),'color','r','linestyle','--');

errorbar(3,nanmean(ses_event_corr_mean(ctrl_animals)), ...
       nanstd(ses_event_corr_mean(ctrl_animals))./sqrt(sum(~isnan(ses_event_corr_mean(ctrl_animals)))),'.k','linewidth',2);
line([2 3],repmat(rand_event_corr_conf_ctrl(1),1,2),'color','r','linestyle','--');
line([2 3],repmat(rand_event_corr_conf_ctrl(2),1,2),'color','r','linestyle','--');
   
xlim([0 4])
ylabel('Mean correlation')
xlabel('Event type')
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Reversal','Session','Session (Ctrl)'});
title('Population correlation (fraction binary active / 30 trials) BEFORE events')


%% 3) Average activity over sessions

rev_animals = [1 2 3 6];
ctrl_animals = [4 5];

session_avg_activity_CL = cellfun(@(x) cellfun(@(x) ...
    nanmean(nanmean(any(vertcat(x{:}),2),1)),x), ...
    epoch_activity_CL,'uni',false);

session_avg_activity_CR = cellfun(@(x) cellfun(@(x) ...
    nanmean(nanmean(any(vertcat(x{:}),2),1)),x), ...
    epoch_activity_CR,'uni',false);

maxdays = max(cellfun(@length,session_avg_activity_CL));
session_avg_activity_CL_grid = nan(length(session_avg_activity_CL),maxdays);
for i = 1:length(session_avg_activity_CL)
   session_avg_activity_CL_grid(i,1:length(session_avg_activity_CL{i})) = ...
       session_avg_activity_CL{i};
end

session_avg_activity_CR_grid = nan(length(session_avg_activity_CR),maxdays);
for i = 1:length(session_avg_activity_CR)
   session_avg_activity_CR_grid(i,1:length(session_avg_activity_CR{i})) = ...
       session_avg_activity_CR{i};
end

figure; hold on;
plot(nanmean(session_avg_activity_CL_grid(rev_animals,:)),'k','linewidth',2)
plot(nanmean(session_avg_activity_CR_grid(rev_animals,:)),'r','linewidth',2)
plot(nanmean(session_avg_activity_CL_grid(ctrl_animals,:)),'--k','linewidth',2)
plot(nanmean(session_avg_activity_CR_grid(ctrl_animals,:)),'--r','linewidth',2)

plot(session_avg_activity_CL_grid(rev_animals,:)','k')
plot(session_avg_activity_CR_grid(rev_animals,:)','r')
plot(session_avg_activity_CL_grid(ctrl_animals,:)','--k')
plot(session_avg_activity_CR_grid(ctrl_animals,:)','--r')

ylabel('Average fraction of active cells');
xlabel('Session');
legend({'CL rev','CR rev','CL ctrl','CR ctrl'});

day_grid = repmat(1:maxdays,size(session_avg_activity_CL_grid,1),1);

temp_y = session_avg_activity_CL_grid(rev_animals,:);
temp_x = day_grid(rev_animals,:);
[r p] = corrcoef(temp_x(~isnan(temp_y)),temp_y(~isnan(temp_y)));

temp_y = session_avg_activity_CR_grid(rev_animals,:);
temp_x = day_grid(rev_animals,:);
[r p] = corrcoef(temp_x(~isnan(temp_y)),temp_y(~isnan(temp_y)));

temp_y = session_avg_activity_CL_grid(ctrl_animals,:);
temp_x = day_grid(rev_animals,:);
[r p] = corrcoef(temp_x(~isnan(temp_y)),temp_y(~isnan(temp_y)));

temp_y = session_avg_activity_CR_grid(ctrl_animals,:);
temp_x = day_grid(rev_animals,:);
[r p] = corrcoef(temp_x(~isnan(temp_y)),temp_y(~isnan(temp_y)));



session_avgcell_activity_CL = cellfun(@(x) cell2mat(cellfun(@(x) ...
    permute(nanmean(any(vertcat(x{:}),2),1),[3 2 1]),x,'uni',false)), ...
    epoch_activity_CL,'uni',false);

session_avgcell_activity_CR = cellfun(@(x) cell2mat(cellfun(@(x) ...
    permute(nanmean(any(vertcat(x{:}),2),1),[3 2 1]),x,'uni',false)), ...
    epoch_activity_CR,'uni',false);



a = cellfun(@(x) diff(x,[],2),session_avgcell_activity_CR(rev_animals),'uni',false);
b = cellfun(@(x) diff(x,[],2),session_avgcell_activity_CL(rev_animals),'uni',false);
a2 = cellfun(@(x) x(:),a,'uni',false);
a3 = vertcat(a2{:});
b2 = cellfun(@(x) x(:),b,'uni',false);
b3 = vertcat(b2{:});
figure;plot(a3,b3,'.k')
ylabel('CL diff')
xlabel('CR diff')
title('Difference in fraction active trials for every cell/animal (rev)')
line([-1 1],[-1 1]);
line([0 0],[-1 1]);
line([-1 1],[0 0]);

%% 4) CL v CR fraction activity

CL_binary = cellfun(@(x) cellfun(@(x) permute(any(vertcat(x{:}),2), ...
    [1 3 2]),x,'uni',false),epoch_activity_CL,'uni',false);
CR_binary = cellfun(@(x) cellfun(@(x) permute(any(vertcat(x{:}),2), ...
    [1 3 2]),x,'uni',false),epoch_activity_CR,'uni',false);

CL_CR_avg_diff = cellfun(@(x,y) cellfun(@(x,y) nanmean(x,1)-nanmean(y,1), ...
    x,y,'uni',false),CL_binary,CR_binary,'uni',false);



%% 5) CR population avg activity across (first) reversal

use_sessions = 1:9;
rev_animals = [1 2 3 6];
ctrl_animals = [4 5];

avg_cat_activity = cellfun(@(x) cellfun(@(x) ...
    permute(nanmean(nanmean(vertcat(x{:}),2),1),[3 1 2]),x,'uni',false), ...
    epoch_activity_CR,'uni',false);

use_sessions_activity_corr = cellfun(@(x) ...
    corrcoef(horzcat(x{use_sessions})),avg_cat_activity,'uni',false);

figure;
for i = 1:length(rev_animals)
    subplot(2,max(length(rev_animals),length(ctrl_animals)),i);
    imagesc(use_sessions_activity_corr{rev_animals(i)});
    colormap(gray);
end

for i = 1:length(ctrl_animals)
    subplot(2,max(length(rev_animals),length(ctrl_animals)),i+length(rev_animals));
    imagesc(use_sessions_activity_corr{ctrl_animals(i)});
    colormap(gray);
end




%% CATEGORY: Single cells

%% 1) CR single cell differences (avoid motion/licking artifacts?)

rev_animals = [1 2 3 6];

min_trials = 30; % pre-post trial count

epoch_activity_CL_binary = cellfun(@(df) ...
    permute(any(df,2),[1 3 2]), ...
    epoch_activity_CL_cat,'uni',false);

rev = cellfun(@(contingency) find(diff(contingency) ~= 0)+1, ...
    contingency_CL_cat,'uni',false);

ses = cellfun(@(sessions) find(diff(sessions) ~= 0)+1, ...
    sessions_CL_cat,'uni',false);

% Get fraction of +/- trials across reversals
event = rev;
pre_event_binary = cell(size(event));
post_event_binary = cell(size(event));
rev_newcontingency = cell(size(event));
for curr_animal = 1:length(epoch_activity_CL_binary)
    curr_event_counter = 1;
    activity = epoch_activity_CL_binary{curr_animal};
    %pre_event_binary{curr_animal} = cell(size(event{curr_animal}));
    %post_event_binary{curr_animal} = cell(size(event{curr_animal}));
    for curr_event = 1:length(event{curr_animal})
        % exclude if minimum forward or back includes another event or
        % exceeds dimensions
        if event{curr_animal}(curr_event)-min_trials < 1 || ...
                event{curr_animal}(curr_event) + min_trials > size(activity,1) || ...
                ~isempty(intersect([event{curr_animal}(curr_event)-min_trials:event{curr_animal}(curr_event)-1 ...
                event{curr_animal}(curr_event)+1:event{curr_animal}(curr_event)+min_trials], ...
                [rev{curr_animal};ses{curr_animal}]))            
            continue
        end
        pre_event_binary{curr_animal}{curr_event_counter} = ...
            activity(event{curr_animal}(curr_event)-min_trials: ...
            event{curr_animal}(curr_event)-1,:);
        post_event_binary{curr_animal}{curr_event_counter} = ...
            activity(event{curr_animal}(curr_event): ...
            event{curr_animal}(curr_event)+min_trials,:);   
        
        rev_newcontingency{curr_animal}(curr_event_counter) = ...
            contingency_CL_cat{curr_animal}(event{curr_animal}(curr_event));
        
        curr_event_counter = curr_event_counter + 1;
    end
end

% Find reversal significance via bootstrap
pre_event_binary_bootci = cellfun(@(activity) ...
    cellfun(@(activity) bootci(10000,@(x) nanmean(x,1),activity), ...
    activity,'uni',false),pre_event_binary(rev_animals),'uni',false);
post_event_binary_bootci = cellfun(@(activity) ...
    cellfun(@(activity) bootci(10000,@(x) nanmean(x,1),activity), ...
    activity,'uni',false),post_event_binary(rev_animals),'uni',false);

rev_newcontingency_revanimals = rev_newcontingency(rev_animals);

% Significance = outside of conf interval AND over 20% in either condition
reversal_sigdiff_cells = cellfun(@(pre,post) cellfun(@(pre,post) ...
    -1*(pre(1,:) > post(2,:)) + 1*(post(1,:) > pre(2,:)) .* ...
    (median(pre) > 0.2 | median(post) > 0.2),pre,post,'uni',false), ...
    pre_event_binary_bootci,post_event_binary_bootci,'uni',false);


%% 2) Get CL-specific cells

session_binary_activity_CL = cellfun(@(x) cellfun(@(x) ...
    permute(any(vertcat(x{:}),2),[3 1 2]),x,'uni',false), ...
    epoch_activity_CL,'uni',false);

session_binary_activity_CR = cellfun(@(x) cellfun(@(x) ...
    permute(any(vertcat(x{:}),2),[3 1 2]),x,'uni',false), ...
    epoch_activity_CR,'uni',false);


p_CL_CR = cellfun(@(x,y) cellfun(@(x,y) AP_chisquare_matrix(x',y'), ...
    x,y,'uni',false),session_binary_activity_CL,session_binary_activity_CR,'uni',false);

% Sig more CL than CR, at least 20% in CL
CL_sigcells = cellfun(@(x,y) cellfun(@(x,y) x < 0.05 & nanmean(y,2)' > ...
    0.2,x,y,'uni',false),p_CL_CR,session_binary_activity_CL,'uni',false);


%% 4) Pre/post fraction scatterplot of control animals

ctrl_animals = [4 5];

% make fake reversals for control animals (halfway through CL trials)
fake_reversals = cellfun(@(x) cellfun(@(x) find(cumsum(x(:,1)) == ...
    round(sum(x(:,1))/2),1),x,'uni',false), ...
    condition_trials_all(ctrl_animals),'uni',false);

fake_boundaries = cellfun(@(rev,odor) cellfun(@(rev,odor) ...
    mat2cell([1;reshape([rev-1 rev]',[],1);size(odor,1)], ...
    repmat(2,length(rev)+1,1)),rev,odor,'uni',false), ...
    fake_reversals,odor_trials_all(ctrl_animals),'uni',false);

% pull out the activity associated with fake reversal epochs
fake_reversal_activity_CL = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,1)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),aligned_df_pyr(ctrl_animals), ...
    condition_trials_all(ctrl_animals),fake_boundaries,'uni',false);

fake_reversal_activity_CL_binary = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    permute(any(x,2),[1 3 2]),x,'uni',false),x,'uni',false),fake_reversal_activity_CL,'uni',false);

% NOTE: below is based on just one reversal

% test for significance
p_fake_reversal = cellfun(@(x) cellfun(@(x) AP_chisquare_matrix(x{1},x{2}), ...
    x,'uni',false),fake_reversal_activity_CL_binary,'uni',false);

% Significance = outside of conf interval AND over 20% in either condition
fake_reversal_sigdiff_cells = cellfun(@(act,p) cellfun(@(act,p) ...
    (p < 0.05) .*...
    (-1*(nanmean(act{1},1) > nanmean(act{2},1)) + 1*(nanmean(act{2},1) > nanmean(act{1},1))) .* ...
    (max([nanmean(act{1},1);nanmean(act{2},1)]) > 0.2),act,p,'uni',false), ...
    fake_reversal_activity_CL_binary,p_fake_reversal,'uni',false);


figure; hold on;
dcm_obj = datacursormode;
set(dcm_obj,'UpdateFcn',@AP_datatip_callback);

animals = [1 2];
sessions = [6 10];
for i = 1:length(animals)
    for j = 1:length(sessions)
        
        animal = animals(i);
        session = sessions(j);
        
        subplot(length(sessions),length(animals),i+(j-1)*length(animals));
        
        sig_color = zeros(length(fake_reversal_sigdiff_cells{animal}{session}),3);
        sig_color(fake_reversal_sigdiff_cells{animal}{session}~=0,1) = 1;
        
        scatter(nanmean(fake_reversal_activity_CL_binary{animal}{session}{1},1), ...
            nanmean(fake_reversal_activity_CL_binary{animal}{session}{2},1), ...
            30*ones(size(sig_color,1),1),sig_color,'fill');
        
        line([0 1],[0 1]);
        xlabel('Pre active fraction')
        ylabel('Post active fraction')
        title(['Animal ' num2str(ctrl_animals(i)) ', Session ' num2str(session)]);
        
    end
end

%% 5) Fraction of significantly different pre/post cells


%%%% Reversal: real reversals
reversal_animals = [1 2 3 6];

%%% CL

rev_activity_CL_binary = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    permute(any(x,2),[1 3 2]),x,'uni',false),x,'uni',false), ...
    epoch_activity_CL(reversal_animals),'uni',false);

p_reversal_CL = cell(size(rev_activity_CL_binary));
reversal_sigdiff_cells_CL = cell(size(rev_activity_CL_binary));
for curr_animal = 1:length(rev_activity_CL_binary)
    for curr_session = 1:length(rev_activity_CL_binary{curr_animal})
        revs = length(rev_activity_CL_binary{curr_animal}{curr_session})-1;
        if revs > 0
            for curr_rev = 1:length(revs)
                curr_pre = rev_activity_CL_binary{curr_animal}{curr_session}{curr_rev};
                curr_post = rev_activity_CL_binary{curr_animal}{curr_session}{curr_rev+1};
                
                p_reversal_CL{curr_animal}{curr_session}{curr_rev} = ...
                    AP_chisquare_matrix(curr_pre,curr_post);
                
                reversal_sigdiff_cells_CL{curr_animal}{curr_session}{curr_rev} = ...
                    (p_reversal_CL{curr_animal}{curr_session}{curr_rev} < 0.05) .*...
                    (-1*(nanmean(curr_pre,1) > nanmean(curr_post,1)) + ...
                    1*(nanmean(curr_post,1) > nanmean(curr_pre,1))) .* ...
                    (max([nanmean(curr_pre,1);nanmean(curr_post,1)]) > 0.2);
            end
        end
    end
end

%%% CR

rev_activity_CR_binary = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    permute(any(x,2),[1 3 2]),x,'uni',false),x,'uni',false), ...
    epoch_activity_CR(reversal_animals),'uni',false);

p_reversal_CR = cell(size(rev_activity_CR_binary));
reversal_sigdiff_cells_CR = cell(size(rev_activity_CR_binary));
for curr_animal = 1:length(rev_activity_CR_binary)
    for curr_session = 1:length(rev_activity_CR_binary{curr_animal})
        revs = length(rev_activity_CR_binary{curr_animal}{curr_session})-1;
        if revs > 0
            for curr_rev = 1:length(revs)
                curr_pre = rev_activity_CR_binary{curr_animal}{curr_session}{curr_rev};
                curr_post = rev_activity_CR_binary{curr_animal}{curr_session}{curr_rev+1};
                
                p_reversal_CR{curr_animal}{curr_session}{curr_rev} = ...
                    AP_chisquare_matrix(curr_pre,curr_post);
                
                reversal_sigdiff_cells_CR{curr_animal}{curr_session}{curr_rev} = ...
                    (p_reversal_CR{curr_animal}{curr_session}{curr_rev} < 0.05) .*...
                    (-1*(nanmean(curr_pre,1) > nanmean(curr_post,1)) + ...
                    1*(nanmean(curr_post,1) > nanmean(curr_pre,1))) .* ...
                    (max([nanmean(curr_pre,1);nanmean(curr_post,1)]) > 0.2);
            end
        end
    end
end

%%%% Control: fake reversals
ctrl_animals = [4 5];

% make fake reversals for control animals (halfway through CL trials)
fake_reversals = cellfun(@(x) cellfun(@(x) find(cumsum(x(:,1)) == ...
    round(sum(x(:,1))/2),1),x,'uni',false), ...
    condition_trials_all(ctrl_animals),'uni',false);

fake_boundaries = cellfun(@(rev,odor) cellfun(@(rev,odor) ...
    mat2cell([1;reshape([rev-1 rev]',[],1);size(odor,1)], ...
    repmat(2,length(rev)+1,1)),rev,odor,'uni',false), ...
    fake_reversals,odor_trials_all(ctrl_animals),'uni',false);


%%% CL

% pull out the activity associated with fake reversal epochs
fake_reversal_activity_CL = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,1)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),aligned_df_pyr(ctrl_animals), ...
    condition_trials_all(ctrl_animals),fake_boundaries,'uni',false);

fake_reversal_activity_CL_binary = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    permute(any(x,2),[1 3 2]),x,'uni',false),x,'uni',false),fake_reversal_activity_CL,'uni',false);

% NOTE: below is based on just one reversal

% test for significance
p_fake_reversal_CL = cellfun(@(x) cellfun(@(x) AP_chisquare_matrix(x{1},x{2}), ...
    x,'uni',false),fake_reversal_activity_CL_binary,'uni',false);

% Significance = outside of conf interval AND over 20% in either condition
fake_reversal_sigdiff_cells_CL = cellfun(@(act,p) cellfun(@(act,p) ...
    (p < 0.05) .*...
    (-1*(nanmean(act{1},1) > nanmean(act{2},1)) + 1*(nanmean(act{2},1) > nanmean(act{1},1))) .* ...
    (max([nanmean(act{1},1);nanmean(act{2},1)]) > 0.2),act,p,'uni',false), ...
    fake_reversal_activity_CL_binary,p_fake_reversal_CL,'uni',false);


%%% CR

% pull out the activity associated with fake reversal epochs
fake_reversal_activity_CR = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,2)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),aligned_df_pyr(ctrl_animals), ...
    condition_trials_all(ctrl_animals),fake_boundaries,'uni',false);

fake_reversal_activity_CR_binary = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    permute(any(x,2),[1 3 2]),x,'uni',false),x,'uni',false),fake_reversal_activity_CR,'uni',false);

% NOTE: below is based on just one reversal

% test for significance
p_fake_reversal_CR = cellfun(@(x) cellfun(@(x) AP_chisquare_matrix(x{1},x{2}), ...
    x,'uni',false),fake_reversal_activity_CR_binary,'uni',false);

% Significance = outside of conf interval AND over 20% in either condition
fake_reversal_sigdiff_cells_CR = cellfun(@(act,p) cellfun(@(act,p) ...
    (p < 0.05) .*...
    (-1*(nanmean(act{1},1) > nanmean(act{2},1)) + 1*(nanmean(act{2},1) > nanmean(act{1},1))) .* ...
    (max([nanmean(act{1},1);nanmean(act{2},1)]) > 0.2),act,p,'uni',false), ...
    fake_reversal_activity_CR_binary,p_fake_reversal_CR,'uni',false);





% Plot the fraction of significant cells in each case

rev_first = num2cell([6 6 6 7]);
frac_up_rev_CL = cellfun(@(x,y) nanmean(x{y}{1} == 1),reversal_sigdiff_cells_CL,rev_first);
frac_down_rev_CL = cellfun(@(x,y) nanmean(x{y}{1} == -1),reversal_sigdiff_cells_CL,rev_first);

frac_up_rev_CR = cellfun(@(x,y) nanmean(x{y}{1} == 1),reversal_sigdiff_cells_CR,rev_first);
frac_down_rev_CR = cellfun(@(x,y) nanmean(x{y}{1} == -1),reversal_sigdiff_cells_CR,rev_first);

ctrl_day = 6;
frac_up_ctrl_CL = cellfun(@(x) nanmean(x{ctrl_day} == 1),fake_reversal_sigdiff_cells_CL);
frac_down_ctrl_CL = cellfun(@(x) nanmean(x{ctrl_day} == -1),fake_reversal_sigdiff_cells_CL);

frac_up_ctrl_CR = cellfun(@(x) nanmean(x{ctrl_day} == 1),fake_reversal_sigdiff_cells_CR);
frac_down_ctrl_CR = cellfun(@(x) nanmean(x{ctrl_day} == -1),fake_reversal_sigdiff_cells_CR);

plot_data = {frac_up_rev_CL frac_up_ctrl_CL frac_up_rev_CR frac_up_ctrl_CR;...
             frac_down_rev_CL frac_down_ctrl_CL frac_down_rev_CR frac_down_ctrl_CR};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    std(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev CL' 'Ctrl CL' 'Rev CR' 'Ctrl CR'});
ylabel('Fraction of significantly modulated cells');



%% (from above: plot the average timecourse from the sig cells above)

figure;
for i = 1:length(reversal_animals)
   subplot(2,length(reversal_animals),i);
   curr_sig = find(reversal_sigdiff_cells_CR{i}{6}{1} == 1);
   plot(permute(nanmean(epoch_activity_CR{reversal_animals(i)}{6}{2}(:,:,curr_sig)),[2 3 1]));
   ylim([0 1]);
   xlim([0 15])
end
title('CR')

for i = 1:length(reversal_animals)
   subplot(2,length(reversal_animals),i+length(reversal_animals));
   curr_sig = find(reversal_sigdiff_cells_CL{i}{6}{1} == 1);
   plot(permute(nanmean(epoch_activity_CL{reversal_animals(i)}{6}{2}(:,:,curr_sig)),[2 3 1]));
   ylim([0 1]);
   xlim([0 15]);
end
title('CL')



%% CATEGORY: Activity timing


%% Timing of cells active in at least 20% of CL trials

high_active_cells = cellfun(@(df) ...
    nanmean(permute(any(df,2),[1 3 2])) > 0.2, ...
    epoch_activity_CL_cat,'uni',false);

epoch_activity_CL_cat_mean = cellfun(@(x) permute(nanmean(x,1),[3 2 1]), ...
    epoch_activity_CL_cat,'uni',false);
















