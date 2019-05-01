% This script was started after quality control: removing contaminated
% ROIs, lowering threshold

%% Define performance cutoff


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




%% Trial-trial correlation barplot (normalized)

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% % all trials
% odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
%     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
%     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);
% 
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



% Correlation prerev/postrev/across
pre_median = nan(size(animals));
post_median = nan(size(animals));
cross_median = nan(size(animals));
for curr_animal = 1:length(animals)
    if ismember(curr_animal,find(rev_animals))
        curr_pre_rev = 1:reversals{curr_animal}(1);
        curr_post_rev = reversals{curr_animal}(1)+1:reversals{curr_animal}(2);
    else
        curr_pre_rev = 1:6;
        curr_post_rev = 7:10;
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

% Normalize each animal to pre-median
pre_median_n = pre_median./pre_median;
post_median_n = post_median./pre_median;
cross_median_n = cross_median./pre_median;

figure;
bar([...
    nanmean(pre_median_n(pmm_animals & rev_animals)) ...
    nanmean(pre_median_n(pmm_animals & ctrl_animals)) ...
    nanmean(pre_median_n(alm_animals & rev_animals)) ...
    nanmean(pre_median_n(alm_animals & ctrl_animals)); ...
    nanmean(post_median_n(pmm_animals & rev_animals)) ...
    nanmean(post_median_n(pmm_animals & ctrl_animals)) ...
    nanmean(post_median_n(alm_animals & rev_animals)) ...
    nanmean(post_median_n(alm_animals & ctrl_animals)); ...
    nanmean(cross_median_n(pmm_animals & rev_animals)) ...
    nanmean(cross_median_n(pmm_animals & ctrl_animals)) ...
    nanmean(cross_median_n(alm_animals & rev_animals)) ...
    nanmean(cross_median_n(alm_animals & ctrl_animals))]);
legend({'PMM rev' 'PMM ctrl' 'ALM rev' 'ALM ctrl'})
set(gca,'XTickLabel',{'Pre-rev' 'Post-rev' 'Cross'})
ylabel('Median trial-trial correlation')
title('Values normalized for each animal to pre-median')


%% First reversal = new patterns of activity?

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% Activity type:
% 1 = trials, temporal
% 2 = trials, binary
% 3 = trials, average
% 4 = average, temporal
% 5 = average, binary

for activity_type = 1:5
    
    % % all trials
    % odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
    %     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
    %     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);
    
    switch activity_type
        
        case 1
            % d' cutoff temporal trials
            odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(x(z,y > 0 & y < 2,:),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);
            act_string = 'tr/time';
            
        case 2
            % d' cutoff binary trial activity
            odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(any(x(z,y > 0 & y < 2,:),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);
            act_string = 'tr/bin';
            
        case 3
            % d' cutoff average trial activity
            odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(x(z,y > 0 & y < 2,:),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);
            act_string = 'tr/avg';
            
        case 4
            % d' cutoff average temporal activity
            odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(nanmean(x(z,y > 0 & y < 2,:),1),[2 3 1]),[],1),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);
            act_string = 'avg/time';
            
        case 5
            % d' cutoff average reliability
            odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(any(x(z,y > 0 & y < 2,:),2),1),[3 2 1]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);
            act_string = 'avg/bin';
            
    end
    
    % Get rid of epochs with too few trials (e.g. accidental reversals)
    num_trial_cutoff = 20;
    use_epochs = cellfun(@(x) cellfun(@(x) cellfun(@(x) sum(x) > num_trial_cutoff,x), ...
        x,'uni',false),d_prime_use_trials_CL_epoch','uni',false);
    
    odor_activity_CL_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
        odor_activity_CL,use_epochs,'uni',false) ;
    
    epoch_contingencies_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
        analysis.epoch_contingencies,use_epochs,'uni',false) ;
    
    
    % Find the first reversal
    first_reversal = cellfun(@(x) find(cellfun(@(x) any(x == 2),x),1), ...
        epoch_contingencies_use,'uni',false);
    
    % Find the second reversal
    second_reversal = cellfun(@(x) find((cellfun(@(x) any(x == 1),x(2:end)) ...
        + cellfun(@(x) any(x == 2),x(1:end-1))) == 2,1)+1, ...
        epoch_contingencies_use,'uni',false);
    
    
    % Get animal group indicies
    rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
    ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
    pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
    alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});
    
    % Extract first pre/post reversal
    for area = 1:2
        
        switch area
            case 1
                use_animals = alm_animals;
                area_string = 'ALM';
            case 2
                use_animals = pmm_animals;
                area_string = 'PMM';
        end
        
        
        % Reversal
        rev_act = cellfun(@(x,y,z) {x(1:y-1) x(y+1:z-1)}, ...
            odor_activity_CL_use(use_animals & rev_animals), ...
            first_reversal(use_animals & rev_animals), ...
            second_reversal(use_animals & rev_animals), 'uni',false);
        
        % get rid of animals where there isn't data in both (this happens
        % in one pmm animal where second reversal was very early???)
        discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),rev_act);
        rev_act(discard_animals) = [];
        
        rev_act_trialcat = cellfun(@(x) cellfun(@(x) ...
            cell2mat(horzcat(x{:})),x,'uni',false),rev_act,'uni',false);
        
        rev_act_trialcorr = cellfun(@(x) corrcoef(cell2mat(x),'rows','complete'), ...
            rev_act_trialcat,'uni',false);
        
        rev_act_trialcorr_median = cellfun(@(x,y) ...
            [nanmedian(AP_itril(x(1:y(1),1:y(1)),-1)),...
            nanmedian(AP_itril(x(y(1)+1:end,y(1)+1:end),-1)), ...
            nanmedian(reshape(x(1:y(1),y(1)+1:end),[],1))], ...
            rev_act_trialcorr, ...
            cellfun(@(x) cellfun(@(x) size(x,2),x),rev_act_trialcat,'uni',false), ...
            'uni',false);
        
        % Control
        % use days 6 and 10 as the reversal days
        ctrl_act = cellfun(@(x) {x(1:5) x(7:9)}, ...
            odor_activity_CL_use(use_animals & ctrl_animals), ...
            'uni',false);
        
        % get rid of animals where there isn't data in both (this happens
        % in one pmm animal where second reversal was very early???)
        discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_act);
        ctrl_act(discard_animals) = [];
        
        ctrl_act_trialcat = cellfun(@(x) cellfun(@(x) ...
            cell2mat(horzcat(x{:})),x,'uni',false),ctrl_act,'uni',false);
        
        ctrl_act_trialcorr = cellfun(@(x) corrcoef(cell2mat(x),'rows','complete'), ...
            ctrl_act_trialcat,'uni',false);
        
        ctrl_act_trialcorr_median = cellfun(@(x,y) ...
            [nanmedian(AP_itril(x(1:y(1),1:y(1)),-1)),...
            nanmedian(AP_itril(x(y(1)+1:end,y(1)+1:end),-1)), ...
            nanmedian(reshape(x(1:y(1),y(1)+1:end),[],1))], ...
            ctrl_act_trialcorr, ...
            cellfun(@(x) cellfun(@(x) size(x,2),x),ctrl_act_trialcat,'uni',false), ...
            'uni',false);
        
        % Control reversals
        ctrl_rev_act = cellfun(@(x,y) {x(y-6:y-1) x(y+1:end)}, ...
            odor_activity_CL_use(use_animals & ctrl_animals), ...
            first_reversal(use_animals & ctrl_animals), 'uni',false);
        
        % get rid of animals where there isn't data in both (this happens
        % in one pmm animal where second reversal was very early???)
        discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_rev_act);
        ctrl_rev_act(discard_animals) = [];
        
        ctrl_rev_act_trialcat = cellfun(@(x) cellfun(@(x) ...
            cell2mat(horzcat(x{:})),x,'uni',false),ctrl_rev_act,'uni',false);
        
        ctrl_rev_act_trialcorr = cellfun(@(x) corrcoef(cell2mat(x),'rows','complete'), ...
            ctrl_rev_act_trialcat,'uni',false);
        
        ctrl_rev_act_trialcorr_median = cellfun(@(x,y) ...
            [nanmedian(AP_itril(x(1:y(1),1:y(1)),-1)),...
            nanmedian(AP_itril(x(y(1)+1:end,y(1)+1:end),-1)), ...
            nanmedian(reshape(x(1:y(1),y(1)+1:end),[],1))], ...
            ctrl_rev_act_trialcorr, ...
            cellfun(@(x) cellfun(@(x) size(x,2),x),ctrl_rev_act_trialcat,'uni',false), ...
            'uni',false);
        
        % Plot
        rev_act_corr = vertcat(rev_act_trialcorr_median{:});
        rev_act_corr_norm = bsxfun(@times,rev_act_corr,1./rev_act_corr(:,1));
        
        ctrl_act_corr = vertcat(ctrl_act_trialcorr_median{:});
        ctrl_act_corr_norm = bsxfun(@times,ctrl_act_corr,1./ctrl_act_corr(:,1));
        
        ctrl_rev_act_corr = vertcat(ctrl_rev_act_trialcorr_median{:});
        ctrl_rev_act_corr_norm = bsxfun(@times,ctrl_rev_act_corr,1./ctrl_rev_act_corr(:,1));
        
        plot_max = max(reshape([rev_act_corr_norm;ctrl_act_corr_norm;ctrl_rev_act_corr_norm],[],1));
        plot_min = min(reshape([rev_act_corr_norm;ctrl_act_corr_norm;ctrl_rev_act_corr_norm],[],1));
        
        figure;
        subplot(1,3,1); hold on;        
        plot(rev_act_corr_norm')
        plot(nanmean(rev_act_corr_norm),'k','linewidth',2);
        set(gca,'XTick',1:3)
        set(gca,'XTickLabel',{'Pre' 'Post' 'Cross'})
        ylabel('Normalized median correlation');
        ylim([plot_min plot_max])
        rev_title = title([area_string ' Rev ' act_string]);
        
        subplot(1,3,2); hold on;       
        plot(ctrl_act_corr_norm')
        plot(nanmean(ctrl_act_corr_norm),'k','linewidth',2);
        set(gca,'XTick',1:3)
        set(gca,'XTickLabel',{'Pre' 'Post' 'Cross'})
        ylabel('Normalized median correlation');
        ylim([plot_min plot_max])
        title([area_string ' Ctrl ' act_string]);
        
        subplot(1,3,3); hold on;        
        plot(ctrl_rev_act_corr_norm')
        plot(nanmean(ctrl_rev_act_corr_norm),'k','linewidth',2);
        set(gca,'XTick',1:3)
        set(gca,'XTickLabel',{'Pre' 'Post' 'Cross'})
        ylabel('Normalized median correlation');
        ylim([plot_min plot_max])
        title([area_string ' Ctrl rev ' act_string]);
        
        % Stats for rev v control cross
        rev_ctrl_cross = [rev_act_corr_norm(:,3);ctrl_act_corr_norm(:,3)];
        rev_ctrl_cross_grp = [ones(size(rev_act_corr_norm,1),1); ...
            2*ones(size(ctrl_act_corr_norm,1),1)];
        n_perm = 1000;
        grp_diff_shuff = nan(n_perm,1);
        for i = 1:n_perm
            grp_diff_shuff(i) = diff(grpstats(rev_ctrl_cross,shake(rev_ctrl_cross_grp)));
        end
        grp_diff = diff(grpstats(rev_ctrl_cross,rev_ctrl_cross_grp));
        grp_diff_rank = tiedrank([grp_diff;grp_diff_shuff]);
        p = grp_diff_rank(1)./(n_perm+1);
        set(rev_title,'String',[get(rev_title,'String') ' p = ' num2str(round(p*100)/100)]);
        
    end
end



%% Correlation by session separation

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% d' cutoff temporal trials
odor_activity_CL = cellfun(@(x,y,z) cellfun(@(x,y,z) cellfun(@(x,z) ...
    reshape(permute(x(z,y > 0 & y < 2,:),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch','uni',false);


% Get rid of epochs with too few trials (e.g. accidental reversals)
num_trial_cutoff = 20;
use_epochs = cellfun(@(x) cellfun(@(x) cellfun(@(x) sum(x) > num_trial_cutoff,x), ...
    x,'uni',false),d_prime_use_trials_CL_epoch','uni',false);

odor_activity_CL_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    odor_activity_CL,use_epochs,'uni',false);

epoch_contingencies_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    analysis.epoch_contingencies,use_epochs,'uni',false);

epoch_sessions_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    analysis.epoch_sessions,use_epochs,'uni',false);


% Find the first reversal
first_reversal = cellfun(@(x) find(cellfun(@(x) any(x == 2),x),1), ...
    epoch_contingencies_use,'uni',false);

% Find the second reversal
second_reversal = cellfun(@(x) find((cellfun(@(x) any(x == 1),x(2:end)) ...
    + cellfun(@(x) any(x == 2),x(1:end-1))) == 2,1)+1, ...
    epoch_contingencies_use,'uni',false);


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});


% Reversal
use_animals = alm_animals;

rev_act = cellfun(@(x,y,z) {x(1:y-1) x(y+1:z-1)}, ...
    odor_activity_CL_use(use_animals & rev_animals), ...
    first_reversal(use_animals & rev_animals), ...
    second_reversal(use_animals & rev_animals), 'uni',false);

rev_sessions = cellfun(@(x,y,z) {x(1:y-1) x(y+1:z-1)}, ...
    epoch_sessions_use(use_animals & rev_animals), ...
    first_reversal(use_animals & rev_animals), ...
    second_reversal(use_animals & rev_animals), 'uni',false);

% get rid of animals where there isn't data in both (this happens
% in one pmm animal where second reversal was very early???)
discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),rev_act);
rev_act(discard_animals) = [];
rev_sessions(discard_animals) = [];

rev_act_trialcat = cellfun(@(x) cellfun(@(x) ...
    cell2mat(horzcat(x{:})),x,'uni',false),rev_act,'uni',false);

rev_act_trialcorr = cellfun(@(x) corrcoef(cell2mat(x),'rows','complete'), ...
    rev_act_trialcat,'uni',false);



% Group trial correlations by sessions separated
rev_sessions_trials = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x,y) cellfun(@(y)...
    repmat(x,1,size(y,2)),y,'uni',false),x,y,'uni',false),x,y,'uni',false),rev_sessions,rev_act,'uni',false);

rev_sessions_trialcat = cellfun(@(x) cellfun(@(x) ...
    cell2mat(horzcat(x{:})),x,'uni',false),rev_sessions_trials,'uni',false);

rev_sessions_trialdiff = cellfun(@(x) abs(repmat(cell2mat(x),length(cell2mat(x)),1) - ...
    repmat(cell2mat(x),length(cell2mat(x)),1)'), ...
    rev_sessions_trialcat,'uni',false);

rev_act_trialcorr_median = cellfun(@(x,y,z) ...
    {grpstats(AP_itril(x(1:y(1),1:y(1)),-1),AP_itril(z(1:y(1),1:y(1)),-1),'median'),...
    grpstats(AP_itril(x(y(1)+1:end,y(1)+1:end),-1),AP_itril(z(y(1)+1:end,y(1)+1:end),-1)), ...
    grpstats(reshape(x(1:y(1),y(1)+1:end),[],1),reshape(z(1:y(1),y(1)+1:end),[],1))}, ...
    rev_act_trialcorr, ...
    cellfun(@(x) cellfun(@(x) size(x,2),x),rev_act_trialcat,'uni',false), ...
    rev_sessions_trialdiff,'uni',false);

rev_act_trialcorr_sessiondiff = cellfun(@(x,y,z) ...
    {unique(AP_itril(z(1:y(1),1:y(1)),-1)),...
    unique(AP_itril(z(y(1)+1:end,y(1)+1:end),-1)), ...
    unique(reshape(z(1:y(1),y(1)+1:end),[],1))}, ...
    rev_act_trialcorr, ...
    cellfun(@(x) cellfun(@(x) size(x,2),x),rev_act_trialcat,'uni',false), ...
    rev_sessions_trialdiff,'uni',false);

figure; hold on;
ls = {'-','--','-.'};
col = lines(length(rev_act_trialcorr_median));
for curr_animal = 1:length(rev_act_trialcorr_median)
    for curr_cond = 1:3
        plot(rev_act_trialcorr_sessiondiff{curr_animal}{curr_cond}, ...
            rev_act_trialcorr_median{curr_animal}{curr_cond}, ...
            'color',col(curr_animal,:),'linestyle',ls{curr_cond});
    end
end

xlabel('Sessions separated')
ylabel('Median correlation')




%% Individual cells difference across reversals - reliability


framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% To select cells significantly active from baseline

% % Chi square
% baseline_sigcells = ...
%     cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
%     AP_chisquare_matrix( ...
%     permute(any(odor(:,seconds < 0,:),2),[1 3 2]), ...
%     permute(any(odor(:,seconds >= 0 & ...
%     seconds < 2,:),2),[1 3 2])) ...
%     < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
%     analysis.epoch_activity_CL,epoch_seconds,...
%     'uni',false);


% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);


% d' cutoff binary trial activity
odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
    permute(any(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds, ...
    d_prime_use_trials_CL_epoch',baseline_sigcells_use,'uni',false);


% Get rid of epochs with too few trials (e.g. accidental reversals)
num_trial_cutoff = 20;
use_epochs = cellfun(@(x) cellfun(@(x) cellfun(@(x) sum(x) > num_trial_cutoff,x), ...
    x,'uni',false),d_prime_use_trials_CL_epoch','uni',false);

odor_activity_CL_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    odor_activity_CL,use_epochs,'uni',false) ;

epoch_contingencies_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    analysis.epoch_contingencies,use_epochs,'uni',false) ;


% Find the first reversal
first_reversal = cellfun(@(x) find(cellfun(@(x) any(x == 2),x),1), ...
    epoch_contingencies_use,'uni',false);

% Find the second reversal
second_reversal = cellfun(@(x) find((cellfun(@(x) any(x == 1),x(2:end)) ...
    + cellfun(@(x) any(x == 2),x(1:end-1))) == 2,1)+1, ...
    epoch_contingencies_use,'uni',false);


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Extract first pre/post reversal
for area = 1:2
    
    switch area
        case 1
            use_animals = alm_animals;
            area_string = 'ALM';
        case 2
            use_animals = pmm_animals;
            area_string = 'PMM';
    end
    
    
    % Reversal
    rev_act = cellfun(@(x,y,z) {x(1:y-1) x(y+1:z-1)}, ...
        odor_activity_CL_use(use_animals & rev_animals), ...
        first_reversal(use_animals & rev_animals), ...
        second_reversal(use_animals & rev_animals), 'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),rev_act);
    rev_act(discard_animals) = [];
    
    rev_act_trialcat = cellfun(@(x) cellfun(@(x) ...
        cell2mat(horzcat(x{:})),x,'uni',false),rev_act,'uni',false);
    
    rev_act_contdiff = cellfun(@(x) diff(cell2mat(cellfun(@(x) ...
        nanmean(x,2),x,'uni',false)),[],2),rev_act_trialcat,'uni',false);
    
    rev_act_contdiff_cat = vertcat(rev_act_contdiff{:});
    
    
    % Control
    % use days 6 and 10 as the reversal days
    ctrl_act = cellfun(@(x) {x(1:5) x(7:9)}, ...
        odor_activity_CL_use(use_animals & ctrl_animals), ...
        'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_act);
    ctrl_act(discard_animals) = [];
    
    ctrl_act_trialcat = cellfun(@(x) cellfun(@(x) ...
        cell2mat(horzcat(x{:})),x,'uni',false),ctrl_act,'uni',false);
    
    ctrl_act_contdiff = cellfun(@(x) diff(cell2mat(cellfun(@(x) ...
        nanmean(x,2),x,'uni',false)),[],2),ctrl_act_trialcat,'uni',false);
    
    ctrl_act_contdiff_cat = vertcat(ctrl_act_contdiff{:});
    
    
    % Control reversals
    ctrl_rev_act = cellfun(@(x,y) {x(y-6:y-1) x(y+1:end)}, ...
        odor_activity_CL_use(use_animals & ctrl_animals), ...
        first_reversal(use_animals & ctrl_animals), 'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_rev_act);
    ctrl_rev_act(discard_animals) = [];
    
    ctrl_rev_act_trialcat = cellfun(@(x) cellfun(@(x) ...
        cell2mat(horzcat(x{:})),x,'uni',false),ctrl_rev_act,'uni',false);
    
    ctrl_rev_act_contdiff = cellfun(@(x) diff(cell2mat(cellfun(@(x) ...
        nanmean(x,2),x,'uni',false)),[],2),ctrl_rev_act_trialcat,'uni',false);
    
    ctrl_rev_act_contdiff_cat = vertcat(ctrl_rev_act_contdiff{:});
    
    
    % Plot
    cell_data = [rev_act_contdiff_cat; ctrl_act_contdiff_cat; ...
        ctrl_rev_act_contdiff_cat];
    
    cell_grps = [ones(length(rev_act_contdiff_cat),1); ...
        2*ones(length(ctrl_act_contdiff_cat),1); ...
        3*ones(length(ctrl_rev_act_contdiff_cat),1)];
    
    figure;
    subplot(2,1,1);
    boxplot(cell_data,cell_grps,'notch','on');
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Rev' 'Ctrl' 'Ctrl rev'})
    ylabel('Post - pre rev')
    title([area_string])
    subplot(2,1,2);
    boxplot(abs(cell_data),cell_grps,'notch','on');
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Rev' 'Ctrl' 'Ctrl rev'})
    ylabel('Post - pre rev')
    title([area_string ' abs'])
    
end



%% Individual cells difference across reversals - temporal


framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);


% To select cells significantly active from baseline

% % Chi square
% baseline_sigcells = ...
%     cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
%     AP_chisquare_matrix( ...
%     permute(any(odor(:,seconds < 0,:),2),[1 3 2]), ...
%     permute(any(odor(:,seconds >= 0 & ...
%     seconds < 2,:),2),[1 3 2])) ...
%     < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
%     analysis.epoch_activity_CL,epoch_seconds,...
%     'uni',false);


% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);


% d' cutoff temporal trials
odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
    x(z,y > 0 & y < 2,sigcells),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds, ...
    d_prime_use_trials_CL_epoch',baseline_sigcells_use,'uni',false);


% Get rid of epochs with too few trials (e.g. accidental reversals)
num_trial_cutoff = 20;
use_epochs = cellfun(@(x) cellfun(@(x) cellfun(@(x) sum(x) > num_trial_cutoff,x), ...
    x,'uni',false),d_prime_use_trials_CL_epoch','uni',false);

odor_activity_CL_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    odor_activity_CL,use_epochs,'uni',false) ;

epoch_contingencies_use = cellfun(@(x,y) cellfun(@(x,y) x(y),x,y,'uni',false), ...
    analysis.epoch_contingencies,use_epochs,'uni',false) ;


% Find the first reversal
first_reversal = cellfun(@(x) find(cellfun(@(x) any(x == 2),x),1), ...
    epoch_contingencies_use,'uni',false);

% Find the second reversal
second_reversal = cellfun(@(x) find((cellfun(@(x) any(x == 1),x(2:end)) ...
    + cellfun(@(x) any(x == 2),x(1:end-1))) == 2,1)+1, ...
    epoch_contingencies_use,'uni',false);


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Extract first pre/post reversal
for area = 1:2
    
    switch area
        case 1
            use_animals = alm_animals;
            area_string = 'ALM';
        case 2
            use_animals = pmm_animals;
            area_string = 'PMM';
    end
    
    
    % Reversal
    rev_act = cellfun(@(x,y,z) {x(1:y-1) x(y+1:z-1)}, ...
        odor_activity_CL_use(use_animals & rev_animals), ...
        first_reversal(use_animals & rev_animals), ...
        second_reversal(use_animals & rev_animals), 'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),rev_act);
    rev_act(discard_animals) = [];
    
    rev_act_trialcat_mean = cellfun(@(x) cellfun(@(x) ...
        permute(nanmean(cell2mat(vertcat(x{:})),1),[2 3 1]),x,'uni',false),rev_act,'uni',false);
    
    rev_act_corr = cellfun(@(x) 1-diag(pdist2(x{1}',x{2}','correlation')), ...
        rev_act_trialcat_mean,'uni',false);
    
    rev_act_corr_cat = vertcat(rev_act_corr{:});
    rev_act_corr_cat(isnan(rev_act_corr_cat)) = [];
    
    % Control
    % use days 6 and 10 as the reversal days
    ctrl_act = cellfun(@(x) {x(1:5) x(7:9)}, ...
        odor_activity_CL_use(use_animals & ctrl_animals), ...
        'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_act);
    ctrl_act(discard_animals) = [];
    
    ctrl_act_trialcat_mean = cellfun(@(x) cellfun(@(x) ...
        permute(nanmean(cell2mat(vertcat(x{:})),1),[2 3 1]),x,'uni',false),ctrl_act,'uni',false);
    
    ctrl_act_corr = cellfun(@(x) 1-diag(pdist2(x{1}',x{2}','correlation')), ...
        ctrl_act_trialcat_mean,'uni',false);
    
    ctrl_act_corr_cat = vertcat(ctrl_act_corr{:});
    ctrl_act_corr_cat(isnan(ctrl_act_corr_cat)) = [];
    
    % Control reversals
    ctrl_rev_act = cellfun(@(x,y) {x(y-6:y-1) x(y+1:end)}, ...
        odor_activity_CL_use(use_animals & ctrl_animals), ...
        first_reversal(use_animals & ctrl_animals), 'uni',false);
    
    % get rid of animals where there isn't data in both (this happens
    % in one pmm animal where second reversal was very early???)
    discard_animals =  cellfun(@(x) any(cellfun(@(x) isempty(x),x)),ctrl_rev_act);
    ctrl_rev_act(discard_animals) = [];
    
    ctrl_rev_act_trialcat_mean = cellfun(@(x) cellfun(@(x) ...
        permute(nanmean(cell2mat(vertcat(x{:})),1),[2 3 1]),x,'uni',false),ctrl_rev_act,'uni',false);
    
    ctrl_rev_act_corr = cellfun(@(x) 1-diag(pdist2(x{1}',x{2}','correlation')), ...
        ctrl_rev_act_trialcat_mean,'uni',false);
    
    ctrl_rev_act_corr_cat = vertcat(ctrl_rev_act_corr{:});
    ctrl_rev_act_corr_cat(isnan(ctrl_rev_act_corr_cat)) = [];
    
    
    % Plot
    cell_data = [rev_act_corr_cat; ctrl_act_corr_cat; ...
        ctrl_rev_act_corr_cat];
    
    cell_grps = [ones(length(rev_act_corr_cat),1); ...
        2*ones(length(ctrl_act_corr_cat),1); ...
        3*ones(length(ctrl_rev_act_corr_cat),1)];
    
    figure;
    subplot(2,1,1);
    boxplot(cell_data,cell_grps,'notch','on');
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Rev' 'Ctrl' 'Ctrl rev'})
    ylabel('Post - pre rev')
    title([area_string])
    subplot(2,1,2);
    boxplot(abs(cell_data),cell_grps,'notch','on');
    set(gca,'XTick',1:3)
    set(gca,'XTickLabel',{'Rev' 'Ctrl' 'Ctrl rev'})
    ylabel('Post - pre rev')
    title([area_string ' abs'])
    
end



%% Slide 11) Correlation grid

% Get correlation of mean odor activity
    
framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% Activity type:
% 1 = trials, temporal
% 2 = trials, binary
% 3 = trials, average
% 4 = average, temporal
% 5 = average, binary

% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);

for activity_type = 1:5;
    
    % % all trials
    % odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
    %     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
    %     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);
    
    switch activity_type
        
        case 1
            % d' cutoff temporal trials
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(x(z,y > 0 & y < 2,sigcells),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            
            act_string = 'tr/time';
            
        case 2
            % d' cutoff binary trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(any(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/bin';
            
        case 3
            % d' cutoff average trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/avg';
            
        case 4
            % d' cutoff average temporal activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(nanmean(x(z,y > 0 & y < 2,sigcells),1),[2 3 1]),[],1),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/time';
            
        case 5
            % d' cutoff average reliability
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(any(x(z,y > 0 & y < 2,sigcells),2),1),[3 2 1]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/bin';
            
    end
    
    if ismember(activity_type,[1 2 3])
        odor_activity_CL_cat = cellfun(@(x) vertcat(x{:})',odor_activity_CL,'uni',false);
    elseif ismember(activity_type,[4 5]);
        odor_activity_CL_cat = cellfun(@(x) cell2mat(vertcat(x{:})'),odor_activity_CL,'uni',false);
    end
    
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
    
%     num_CL_trials = cellfun(@(x) cellfun(@(x) size(x,1),vertcat(x{:})), ...
%         analysis.epoch_activity_CL,'uni',false);
%     num_trial_cutoff = 20;
%     use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_CL_trials,'uni',false);
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
    
    % Get activity correlation
    if ismember(activity_type,[1 2 3]);
        act_corr_trial = cellfun(@(x) mat2cell(corrcoef(horzcat(x{:}),'rows','complete'), ...
            cellfun(@(x) size(x,2),x),cellfun(@(x) size(x,2),x)), ...
            odor_activity_CL_cat_use,'uni',false);
        act_corr = cellfun(@(x) cellfun(@(x) nanmedian(x(:)),x),act_corr_trial,'uni',false);
        for i = 1:length(act_corr)
             act_corr{i}(logical(eye(size(act_corr{i})))) = ...
                 cellfun(@(x) nanmedian(AP_itril(x,-1)), ...
                 act_corr_trial{i}(logical(eye(size(act_corr{i})))));
        end
    elseif ismember(activity_type,[4 5]);
        act_corr = cellfun(@(x) corrcoef(x,'rows','complete'), ...
            odor_activity_CL_cat_use,'uni',false);
    end
    
    % Get pre/post/cross rev 1 comparison
    rev_1 = nan(1,length(reversals));
    rev_1(rev_animals) = cellfun(@(x) x(1),reversals(rev_animals));
    rev_1(ctrl_animals) = 5;
    rev_2 = nan(1,length(reversals));
    for curr_animal = find(rev_animals);
        if length(reversals{curr_animal}) >= 2;
            rev_2(curr_animal) = reversals{curr_animal}(2);
        else
            rev_2(curr_animal) = length(act_corr{curr_animal});
        end
    end
    rev_2(ctrl_animals) = 9;
    if ismember(activity_type,[1 2 3]);
        pre_mean_corr = cellfun(@(x,rev) nanmedian(AP_itril(x(1:rev,1:rev))), ...
            act_corr,num2cell(rev_1));
        post_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(AP_itril(x(rev1+1:rev2,rev1+1:rev2))), ...
            act_corr,num2cell(rev_1),num2cell(rev_2));
    elseif ismember(activity_type,[4 5]);
        pre_mean_corr = cellfun(@(x,rev) nanmedian(AP_itril(x(1:rev,1:rev),-1)), ...
            act_corr,num2cell(rev_1));
        post_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(AP_itril(x(rev1+1:rev2,rev1+1:rev2),-1)), ...
            act_corr,num2cell(rev_1),num2cell(rev_2));
    end
    cross_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(reshape(x(1:rev1,rev1+1:rev2),[],1)), ...
        act_corr,num2cell(rev_1),num2cell(rev_2));
    contingency_corr = [pre_mean_corr;post_mean_corr;cross_mean_corr];
    contingency_corr = bsxfun(@times,contingency_corr,1./contingency_corr(1,:));
    
    % can make this a loop
    curr_rev = 1;
    
    figure; 
    subplot(2,2,1);
    plot(contingency_corr(:,pmm_animals & rev_animals));
    title(['PMM Reversal ' num2str(curr_rev)]);
    ylabel('Normalized median correlation');
    subplot(2,2,2);
    plot(contingency_corr(:,pmm_animals & ctrl_animals));
    title(['PMM Control ' num2str(curr_rev)]);
    ylabel('Normalized median correlation');
    subplot(2,2,3);
    plot(contingency_corr(:,alm_animals & rev_animals));
    title(['ALM Reversal ' num2str(curr_rev)]);
    ylabel('Normalized median correlation');
    subplot(2,2,4);
    plot(contingency_corr(:,alm_animals & ctrl_animals));
    title(['ALM Control ' num2str(curr_rev)]);
    ylabel('Normalized median correlation');
    
    set(gcf,'Name',act_string);
    
    % Plot grid aligned to reversal
    max_rev = max(cellfun(@length,reversals(rev_animals)));
    
    max_epochs = max(cellfun(@length,act_corr));
    act_corr_pad_forward = cellfun(@(x) padarray(x,[max_epochs-length(x) ...
        max_epochs-length(x)],NaN,'post'),act_corr,'uni',false);
    act_corr_pad = cellfun(@(x) padarray(x,[max_epochs ...
        max_epochs],NaN,'pre'),act_corr_pad_forward,'uni',false);
    
    figure;
    
    % shift all correlation grids to align by reversal
    use_pmm = rev_animals & pmm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_pmm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_pmm),reversals(use_pmm),'uni',false),[1 3 2])),3);
    
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(2,2,1);
    imagesc(curr_pmm_grid_avg);
    caxis([0 1])
    colormap(gray);
    line(xlim,[max_epochs-0.5 max_epochs-0.5],'color','r');
    line([max_epochs-0.5 max_epochs-0.5],ylim,'color','r');
    curr_1 = find(~all(isnan(curr_pmm_grid_avg),1),1)-0.5;
    curr_x = [find(~all(isnan(curr_pmm_grid_avg),1),1)-0.5 ...
        size(curr_pmm_grid_avg,2)-find(~all(isnan(fliplr(curr_pmm_grid_avg)),1),1)+1.5];
    curr_y = [find(~all(isnan(curr_pmm_grid_avg),2),1)-0.5 ...
        size(curr_pmm_grid_avg,2)-find(~all(isnan(flipud(curr_pmm_grid_avg)),2),1)+1.5];
    xlim([curr_1 curr_1+10])
    ylim([curr_1 curr_1+10])
    title(['PMM Reversal ' num2str(curr_rev)]);
    
    % Plot control grid
    subplot(2,2,2)
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
    caxis([0 1])
    curr_x = [find(~all(isnan(nanmean(curr_grid,3)),1),1)-0.5 ...
        size(nanmean(curr_grid,3),2)-find(~all(isnan(fliplr(nanmean(curr_grid,3))),1),1)+1.5];
    curr_y = [find(~all(isnan(nanmean(curr_grid,3)),2),1)-0.5 ...
        size(nanmean(curr_grid,3),2)-find(~all(isnan(flipud(nanmean(curr_grid,3))),2),1)+1.5];
    xlim(curr_x)
    ylim(curr_y)
    title('PMM control')
    
    use_alm = rev_animals & alm_animals & ...
        cellfun(@length,reversals) >= curr_rev;
    curr_alm_grid_avg = nanmean(cell2mat(permute(cellfun(@(grid,rev) circshift(grid, ...
        [-rev(curr_rev),-rev(curr_rev)]), ...
        act_corr_pad(use_alm),reversals(use_alm),'uni',false),[1 3 2])),3);
    
    %subplot(sq_plot,sq_plot,curr_rev);
    subplot(2,2,3);
    imagesc(curr_alm_grid_avg);
    caxis([0 1])
    colormap(gray);
    line(xlim,[max_epochs-0.5 max_epochs-0.5],'color','r');
    line([max_epochs-0.5 max_epochs-0.5],ylim,'color','r');
    curr_1 = find(~all(isnan(curr_alm_grid_avg),1),1)-0.5;
    curr_x = [find(~all(isnan(curr_alm_grid_avg),1),1)-0.5 ...
        size(curr_alm_grid_avg,2)-find(~all(isnan(fliplr(curr_alm_grid_avg)),1),1)+1.5];
    curr_y = [find(~all(isnan(curr_alm_grid_avg),2),1)-0.5 ...
        size(curr_alm_grid_avg,2)-find(~all(isnan(flipud(curr_alm_grid_avg)),2),1)+1.5];
    xlim([curr_1 curr_1+10])
    ylim([curr_1 curr_1+10])
    title(['ALM Reversal ' num2str(curr_rev)]);
    
    subplot(2,2,4)
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
    caxis([0 1])
    curr_x = [find(~all(isnan(nanmean(curr_grid,3)),1),1)-0.5 ...
        size(nanmean(curr_grid,3),2)-find(~all(isnan(fliplr(nanmean(curr_grid,3))),1),1)+1.5];
    curr_y = [find(~all(isnan(nanmean(curr_grid,3)),2),1)-0.5 ...
        size(nanmean(curr_grid,3),2)-find(~all(isnan(flipud(nanmean(curr_grid,3))),2),1)+1.5];
    xlim(curr_x)
    ylim(curr_y)
    title('ALM control')
    
    set(gcf,'Name',act_string);
    
    
end

for curr_animal = find(pmm_animals & rev_animals);
figure;
imagesc(act_corr{curr_animal});colormap(gray);caxis([0 1]);
ylim([0.5 reversals{curr_animal}(2)+0.5]);xlim([0.5 reversals{curr_animal}(2)+0.5])
line(xlim,[reversals{curr_animal}(1)+0.5 reversals{curr_animal}(1)+0.5],'color','r');
line([reversals{curr_animal}(1)+0.5 reversals{curr_animal}(1)+0.5],ylim,'color','r');
title(mice(curr_animal).name)
end

for curr_animal = find(alm_animals & ctrl_animals);
figure;
imagesc(act_corr{curr_animal});colormap(gray);caxis([0 1]);
ylim([0.5 10.5]);xlim([0.5 10.5])
title(mice(curr_animal).name)
end



%% Pre/post/cross for cell = n

% Get correlation of mean odor activity

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);


% % all trials
% odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
%     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
%     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);

% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);


% d' cutoff temporal trials
% tr x time x cell
odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
    x(z,y > 0 & y < 2,sigcells),x,z,'uni',false), ...
    x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds, ...
    d_prime_use_trials_CL_epoch',baseline_sigcells_use,'uni',false);
act_string = 'tr/time';

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

%     num_CL_trials = cellfun(@(x) cellfun(@(x) size(x,1),vertcat(x{:})), ...
%         analysis.epoch_activity_CL,'uni',false);
%     num_trial_cutoff = 20;
%     use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_CL_trials,'uni',false);

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


% Get activity correlation for each cell & session
cell_corr = cell(length(odor_activity_CL_cat_use),1);
for curr_animal = 1:length(odor_activity_CL_cat_use)
    
    n_cells = size(odor_activity_CL_cat_use{curr_animal}{1},3);
    curr_cell_corr = nan(length(odor_activity_CL_cat_use{curr_animal}), ...
        length(odor_activity_CL_cat_use{curr_animal}),n_cells);
    
    for curr_cell = 1:n_cells
        curr_act = cellfun(@(x) x(:,:,curr_cell)', ...
            odor_activity_CL_cat_use{curr_animal},'uni',false);
        
        curr_corr_trial = mat2cell(corrcoef(horzcat(curr_act{:}),'rows','complete'), ...
            cellfun(@(x) size(x,2),curr_act),cellfun(@(x) size(x,2),curr_act));
        
        curr_corr = cellfun(@(x) nanmedian(x(:)),curr_corr_trial);
        
        curr_corr(logical(eye(size(curr_corr)))) = ...
            cellfun(@(x) nanmedian(AP_itril(x,-1)), ...
            curr_corr_trial(logical(eye(size(curr_corr)))));
        
        curr_cell_corr(:,:,curr_cell) = curr_corr;
    end   
    cell_corr{curr_animal} = curr_cell_corr;
end

% Get pre/post/cross rev 1 comparison
rev_1 = nan(1,length(reversals));
rev_1(rev_animals) = cellfun(@(x) x(1),reversals(rev_animals));
rev_1(ctrl_animals) = 5;
rev_2 = nan(1,length(reversals));
for curr_animal = find(rev_animals);
    if length(reversals{curr_animal}) >= 2;
        rev_2(curr_animal) = reversals{curr_animal}(2);
    else
        rev_2(curr_animal) = length(act_corr{curr_animal});
    end
end
rev_2(ctrl_animals) = 9;

% For whole block
% pre_mean_corr = cellfun(@(x,rev) arrayfun(@(y) ...
%     nanmedian(AP_itril(x(1:rev,1:rev,y))),1:size(x,3)), ...
%     cell_corr',num2cell(rev_1),'uni',false);
% post_mean_corr = cellfun(@(x,rev1,rev2) arrayfun(@(y) ...
%     nanmedian(AP_itril(x(rev1+1:rev2,rev1+1:rev2,y))),1:size(x,3)), ...
%     cell_corr',num2cell(rev_1),num2cell(rev_2),'uni',false);
% cross_mean_corr = cellfun(@(x,rev1,rev2) arrayfun(@(y) ...
%     nanmedian(reshape(x(1:rev1,rev1+1:rev2,y),[],1)),1:size(x,3)), ...
%     cell_corr',num2cell(rev_1),num2cell(rev_2),'uni',false);

% For n-back from reversal
n_back = 2;

pre_mean_corr = cellfun(@(x,rev) arrayfun(@(y) ...
    nanmedian(AP_itril(x(rev-n_back:rev,rev-n_back:rev,y))),1:size(x,3)), ...
    cell_corr',num2cell(rev_1),'uni',false);
post_mean_corr = cellfun(@(x,rev1,rev2) arrayfun(@(y) ...
    nanmedian(AP_itril(x(rev2-n_back:rev2,rev2-n_back:rev2,y))),1:size(x,3)), ...
    cell_corr',num2cell(rev_1),num2cell(rev_2),'uni',false);
cross_mean_corr = cellfun(@(x,rev1,rev2) arrayfun(@(y) ...
    nanmedian(reshape(x(rev1-n_back:rev1,rev2-n_back:rev2,y),[],1)),1:size(x,3)), ...
    cell_corr',num2cell(rev_1),num2cell(rev_2),'uni',false);

contingency_corr = [pre_mean_corr;post_mean_corr;cross_mean_corr];

figure;
subplot(2,2,1);
curr_plot = cell2mat(contingency_corr(:,pmm_animals & rev_animals));
errorbar(nanmedian(curr_plot'),nanstd(curr_plot')./sqrt(sum(~isnan(curr_plot'))),'k','linewidth',2)
%boxplot(curr_plot','notch','on');
title(['PMM Reversal ' num2str(curr_rev)]);
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Pre','Post','Cross'});
subplot(2,2,2);
curr_plot = cell2mat(contingency_corr(:,pmm_animals & ctrl_animals));
errorbar(nanmedian(curr_plot'),nanstd(curr_plot')./sqrt(sum(~isnan(curr_plot'))),'k','linewidth',2)
%boxplot(curr_plot','notch','on');
title(['PMM Control ' num2str(curr_rev)]);
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Pre','Post','Cross'});
subplot(2,2,3);
curr_plot = cell2mat(contingency_corr(:,alm_animals & rev_animals));
errorbar(nanmedian(curr_plot'),nanstd(curr_plot')./sqrt(sum(~isnan(curr_plot'))),'k','linewidth',2)
%boxplot(curr_plot','notch','on');
title(['ALM Reversal ' num2str(curr_rev)]);
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Pre','Post','Cross'});
subplot(2,2,4);
curr_plot = cell2mat(contingency_corr(:,alm_animals & ctrl_animals));
errorbar(nanmedian(curr_plot'),nanstd(curr_plot')./sqrt(sum(~isnan(curr_plot'))),'k','linewidth',2)
%boxplot(curr_plot','notch','on');
title(['ALM Control ' num2str(curr_rev)]);
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Pre','Post','Cross'});

set(gcf,'Name',[act_string ' cell = n']);


%% Plot 1st and 2nd diagonal of individual animal correlations

% Get correlation of mean odor activity
    
framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% Activity type:
% 1 = trials, temporal
% 2 = trials, binary
% 3 = trials, average
% 4 = average, temporal
% 5 = average, binary

% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);

for activity_type = 1:3;
    
    % % all trials
    % odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
    %     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
    %     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);
    
    switch activity_type
        
        case 1
            % d' cutoff temporal trials
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(x(z,y > 0 & y < 2,sigcells),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            
            act_string = 'tr/time';
            
        case 2
            % d' cutoff binary trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(any(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/bin';
            
        case 3
            % d' cutoff average trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/avg';
            
        case 4
            % d' cutoff average temporal activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(nanmean(x(z,y > 0 & y < 2,sigcells),1),[2 3 1]),[],1),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/time';
            
        case 5
            % d' cutoff average reliability
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(any(x(z,y > 0 & y < 2,sigcells),2),1),[3 2 1]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/bin';
            
    end
    
    if ismember(activity_type,[1 2 3])
        odor_activity_CL_cat = cellfun(@(x) vertcat(x{:})',odor_activity_CL,'uni',false);
    elseif ismember(activity_type,[4 5]);
        odor_activity_CL_cat = cellfun(@(x) cell2mat(vertcat(x{:})'),odor_activity_CL,'uni',false);
    end
    
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
    
%     num_CL_trials = cellfun(@(x) cellfun(@(x) size(x,1),vertcat(x{:})), ...
%         analysis.epoch_activity_CL,'uni',false);
%     num_trial_cutoff = 20;
%     use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_CL_trials,'uni',false);
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
    
    % Get activity correlation
    if ismember(activity_type,[1 2 3]);
        act_corr_trial = cellfun(@(x) mat2cell(corrcoef(horzcat(x{:}),'rows','complete'), ...
            cellfun(@(x) size(x,2),x),cellfun(@(x) size(x,2),x)), ...
            odor_activity_CL_cat_use,'uni',false);
        act_corr = cellfun(@(x) cellfun(@(x) nanmedian(x(:)),x),act_corr_trial,'uni',false);
        for i = 1:length(act_corr)
             act_corr{i}(logical(eye(size(act_corr{i})))) = ...
                 cellfun(@(x) nanmedian(AP_itril(x,-1)), ...
                 act_corr_trial{i}(logical(eye(size(act_corr{i})))));
        end
    elseif ismember(activity_type,[4 5]);
        act_corr = cellfun(@(x) corrcoef(x,'rows','complete'), ...
            odor_activity_CL_cat_use,'uni',false);
    end
    
    % Plot first (within-day) and second (next day) diagonals
    first_diag = cellfun(@(x) diag(x),act_corr,'uni',false);
    second_diag = cellfun(@(x) diag(x,-1),act_corr,'uni',false);
    
    rev_1 = nan(1,length(reversals));
    rev_1(rev_animals) = cellfun(@(x) x(1),reversals(rev_animals));
    rev_1(ctrl_animals) = 5;
    rev_2 = nan(1,length(reversals));
    rev_2(rev_animals) = cellfun(@(x) x(2),reversals(rev_animals));
    rev_2(ctrl_animals) = 9;
    
    
    figure;
    subplot(2,2,1); hold on;
    curr_animals = find(pmm_animals & rev_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal):rev_2(curr_animal)-(rev_1(curr_animal)+1);
        plot(x,first_diag{curr_animal}(1:rev_2(curr_animal)), ...
            'color',col(curr_animal_idx,:));
    end
    title('PMM rev 1')
    subplot(2,2,2); hold on;
    curr_animals = find(pmm_animals & ctrl_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal):rev_2(curr_animal)-(rev_1(curr_animal)+1);
        plot(x,first_diag{curr_animal}(1:rev_2(curr_animal)), ...
            'color',col(curr_animal_idx,:));
    end
    title('PMM ctrl 1')
    subplot(2,2,3); hold on;
    curr_animals = find(alm_animals & rev_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal):rev_2(curr_animal)-(rev_1(curr_animal)+1);
        plot(x,first_diag{curr_animal}(1:rev_2(curr_animal)), ...
            'color',col(curr_animal_idx,:));
    end
    title('ALM rev 1')
    subplot(2,2,4); hold on;
    curr_animals = find(alm_animals & ctrl_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal):rev_2(curr_animal)-(rev_1(curr_animal)+1);
        plot(x,first_diag{curr_animal}(1:rev_2(curr_animal)), ...
            'color',col(curr_animal_idx,:));
    end
    title('ALM ctrl 1')
      
    set(gcf,'Name',act_string);
    
    figure;
    subplot(2,2,1); hold on;
    curr_animals = find(pmm_animals & rev_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal)+1:rev_2(curr_animal)-(rev_1(curr_animal))-1;
        plot(x,second_diag{curr_animal}(1:rev_2(curr_animal)-1), ...
            'color',col(curr_animal_idx,:));
    end
    title('PMM rev 2')
    subplot(2,2,2); hold on;
    curr_animals = find(pmm_animals & ctrl_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal)+1:rev_2(curr_animal)-(rev_1(curr_animal))-1;
        plot(x,second_diag{curr_animal}(1:rev_2(curr_animal)-1), ...
            'color',col(curr_animal_idx,:));
    end
    title('PMM ctrl 2')
    subplot(2,2,3); hold on;
    curr_animals = find(alm_animals & rev_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal)+1:rev_2(curr_animal)-(rev_1(curr_animal))-1;
        plot(x,second_diag{curr_animal}(1:rev_2(curr_animal)-1), ...
            'color',col(curr_animal_idx,:));
    end
    title('ALM rev 2')
    subplot(2,2,4); hold on;
    curr_animals = find(alm_animals & ctrl_animals);
    col = lines(length(curr_animals));
    for curr_animal_idx = 1:length(curr_animals);
        curr_animal = curr_animals(curr_animal_idx);
        x = -rev_1(curr_animal)+1:rev_2(curr_animal)-(rev_1(curr_animal))-1;
        plot(x,second_diag{curr_animal}(1:rev_2(curr_animal)-1), ...
            'color',col(curr_animal_idx,:));
    end
    title('ALM ctrl 2')
    
    set(gcf,'Name',act_string);
    
end


%% Pre/post/cross n-back from reversal

% Get correlation of mean odor activity
    
framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

% Activity type:
% 1 = trials, temporal
% 2 = trials, binary
% 3 = trials, average
% 4 = average, temporal
% 5 = average, binary

% Rank sum
baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% For now, just use cells that are ever significant
baseline_sigcells_use = cellfun(@(x) any(cell2mat(vertcat(x{:})),1), ...
    baseline_sigcells,'uni',false);

for activity_type = 1;
    
    % % all trials
    % odor_activity_CL = cellfun(@(x,y) cellfun(@(x,y) cellfun(@(x) ...
    %     reshape(permute(x(:,y > 0 & y < 2,:),[2 3 1]),[],size(x,1)),x,'uni',false), ...
    %     x,y,'uni',false),analysis.epoch_activity_CL,epoch_seconds,'uni',false);
    
    switch activity_type
        
        case 1
            % d' cutoff temporal trials
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(x(z,y > 0 & y < 2,sigcells),[2 3 1]),[],sum(z)),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            
            act_string = 'tr/time';
            
        case 2
            % d' cutoff binary trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(any(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/bin';
            
        case 3
            % d' cutoff average trial activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(x(z,y > 0 & y < 2,sigcells),2),[3 1 2]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'tr/avg';
            
        case 4
            % d' cutoff average temporal activity
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                reshape(permute(nanmean(x(z,y > 0 & y < 2,sigcells),1),[2 3 1]),[],1),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/time';
            
        case 5
            % d' cutoff average reliability
            odor_activity_CL = cellfun(@(x,y,z,sigcells) cellfun(@(x,y,z) cellfun(@(x,z) ...
                permute(nanmean(any(x(z,y > 0 & y < 2,sigcells),2),1),[3 2 1]),x,z,'uni',false), ...
                x,y,z,'uni',false),analysis.epoch_activity_CL,epoch_seconds,d_prime_use_trials_CL_epoch', ...
                baseline_sigcells_use,'uni',false);
            act_string = 'avg/bin';
            
    end
    
    if ismember(activity_type,[1 2 3])
        odor_activity_CL_cat = cellfun(@(x) vertcat(x{:})',odor_activity_CL,'uni',false);
    elseif ismember(activity_type,[4 5]);
        odor_activity_CL_cat = cellfun(@(x) cell2mat(vertcat(x{:})'),odor_activity_CL,'uni',false);
    end
    
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
    
%     num_CL_trials = cellfun(@(x) cellfun(@(x) size(x,1),vertcat(x{:})), ...
%         analysis.epoch_activity_CL,'uni',false);
%     num_trial_cutoff = 20;
%     use_epochs = cellfun(@(x) x >= num_trial_cutoff,num_CL_trials,'uni',false);
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
    
    % Get activity correlation
    if ismember(activity_type,[1 2 3]);
        act_corr_trial = cellfun(@(x) mat2cell(corrcoef(horzcat(x{:}),'rows','complete'), ...
            cellfun(@(x) size(x,2),x),cellfun(@(x) size(x,2),x)), ...
            odor_activity_CL_cat_use,'uni',false);
        act_corr = cellfun(@(x) cellfun(@(x) nanmedian(x(:)),x),act_corr_trial,'uni',false);
        for i = 1:length(act_corr)
             act_corr{i}(logical(eye(size(act_corr{i})))) = ...
                 cellfun(@(x) nanmedian(AP_itril(x,-1)), ...
                 act_corr_trial{i}(logical(eye(size(act_corr{i})))));
        end
    elseif ismember(activity_type,[4 5]);
        act_corr = cellfun(@(x) corrcoef(x,'rows','complete'), ...
            odor_activity_CL_cat_use,'uni',false);
    end
    
    % Get pre/post/cross rev 1 comparison
    rev_1 = nan(1,length(reversals));
    rev_1(rev_animals) = cellfun(@(x) x(1),reversals(rev_animals));
    rev_1(ctrl_animals) = 5;
    rev_2 = nan(1,length(reversals));
    for curr_animal = find(rev_animals);
        if length(reversals{curr_animal}) >= 2;
            rev_2(curr_animal) = reversals{curr_animal}(2);
        else
            rev_2(curr_animal) = length(act_corr{curr_animal});
        end
    end
    rev_2(ctrl_animals) = 9;
    
    n_back = 2;
    
    % To use diagonal
%     if ismember(activity_type,[1 2 3]);
%         pre_mean_corr = cellfun(@(x,rev) nanmedian(AP_itril(x((rev-n_back):rev,(rev-n_back):rev))), ...
%             act_corr,num2cell(rev_1));
%         post_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(AP_itril(x((rev2-n_back):rev2,(rev2-n_back):rev2))), ...
%             act_corr,num2cell(rev_1),num2cell(rev_2));
%     elseif ismember(activity_type,[4 5]);
        pre_mean_corr = cellfun(@(x,rev) nanmedian(AP_itril(x((rev-n_back):rev,(rev-n_back):rev),-1)), ...
            act_corr,num2cell(rev_1));
        post_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(AP_itril(x((rev2-n_back):rev2,(rev2-n_back):rev2),-1)), ...
            act_corr,num2cell(rev_1),num2cell(rev_2));
%     end
    cross_mean_corr = cellfun(@(x,rev1,rev2) nanmedian(reshape(x((rev1-n_back):rev1,rev2-(n_back):rev2),[],1)), ...
        act_corr,num2cell(rev_1),num2cell(rev_2));
    contingency_corr = [pre_mean_corr;post_mean_corr;cross_mean_corr];
    %contingency_corr = bsxfun(@times,contingency_corr,1./contingency_corr(1,:));
    
    figure; 
    subplot(2,2,1);
    plot(contingency_corr(:,pmm_animals & rev_animals));
    title(['PMM Reversal ' num2str(curr_rev)]);
    ylabel('median correlation');
    subplot(2,2,2);
    plot(contingency_corr(:,pmm_animals & ctrl_animals));
    title(['PMM Control ' num2str(curr_rev)]);
    ylabel('median correlation');
    subplot(2,2,3);
    plot(contingency_corr(:,alm_animals & rev_animals));
    title(['ALM Reversal ' num2str(curr_rev)]);
    ylabel('median correlation');
    subplot(2,2,4);
    plot(contingency_corr(:,alm_animals & ctrl_animals));
    title(['ALM Control ' num2str(curr_rev)]);
    ylabel('median correlation');
    
    set(gcf,'Name',[act_string ' ' num2str(n_back) ' back,no diag']);

end



























