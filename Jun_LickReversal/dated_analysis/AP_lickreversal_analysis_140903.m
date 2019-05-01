%% Script info

% This script is after toolbox included labels, analysis, 
% upgraded mom/bscope used, and ALM included (old scripts no longer work)


%% Fraction of significantly different pre/post cells

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% get cells which are significantly different CR/CL within epoch: chi
% square p < 0.05 and active at least 20% or trials in one condition

CL_CR_sigdiff = cellfun(@(CL,CR) cellfun(@(CL,CR) cellfun(@(CL,CR) ...
    AP_chisquare_matrix( ...
    permute(any(CL,2),[1 3 2]),permute(any(CR,2),[1 3 2])) < 0.05 & ...
    max([nanmean(permute(any(CL,2),[1 3 2]),1);nanmean(permute(any(CR,2),[1 3 2]),1)],[],1) > 0.2, ...
    CL,CR,'uni',false),CL,CR,'uni',false),analysis.epoch_activity_CL,...
    analysis.epoch_activity_CR,'uni',false);

% get cells which are significantly active compared to baseline: chi
% square p < 0.05 and active on at least 10% of trials

odor_sigdiff = cellfun(@(odor,baseline) cellfun(@(odor,baseline)...
    AP_chisquare_matrix( ...
    permute(any(odor,2),[1 3 2]),permute(any(baseline,2),[1 3 2])) < 0.05 & ...
    nanmean(permute(any(odor,2),[1 3 2]),1) > 0.2, ...
    odor,baseline,'uni',false),analysis.aligned_df_all,...
    analysis.baseline_df_all,'uni',false);





CL_frac_active = cellfun(@(CL) cellfun(@(CL) cellfun(@(CL) ...
    permute(any(CL,2),[1 3 2]),CL,'uni',false),CL,'uni',false), ...
    analysis.epoch_activity_CL,'uni',false);


%% Temporal activity of CL/CR cells

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% get cells which are significantly different CR/CL within epoch: chi
% square p < 0.05 and active at least x% or trials in one condition

CL_CR_sigdiff = cellfun(@(CL,CR) cellfun(@(CL,CR) cellfun(@(CL,CR) ...
    (1*(nanmean(permute(any(CL,2),[1 3 2]),1) > nanmean(permute(any(CR,2),[1 3 2]),1)) + ...
    -1*(nanmean(permute(any(CL,2),[1 3 2]),1) < nanmean(permute(any(CR,2),[1 3 2]),1))).* ...
    AP_chisquare_matrix( ...
    permute(any(CL,2),[1 3 2]),permute(any(CR,2),[1 3 2])) < 0.05 & ...
    max([nanmean(permute(any(CL,2),[1 3 2]),1);nanmean(permute(any(CR,2),[1 3 2]),1)],[],1) > 0.1, ...
    CL,CR,'uni',false),CL,CR,'uni',false),analysis.epoch_activity_CL,...
    analysis.epoch_activity_CR,'uni',false);


% get cells which are significantly active compared to baseline (paired
% trials for baseline / trials, so signrank of average df/f for given time)
% and active on at least x% of trials

odor_sigdiff = cellfun(@(odor,baseline) cellfun(@(odor,baseline) cellfun(@(odor,baseline) ...
    AP_signrank_matrix( ...
    permute(nanmean(odor,2),[1 3 2]),permute(nanmean(baseline,2),[1 3 2])) < 0.05 & ...
    nanmean(permute(any(odor,2),[1 3 2]),1) > 0.1, ...
    odor,baseline,'uni',false),odor,baseline,'uni',false), ...
    analysis.epoch_activity_odor,analysis.epoch_activity_baseline,...
    'uni',false);


CL_avg_activity = cellfun(@(CL) cellfun(@(CL) cellfun(@(CL) ...
    permute(nanmean(CL,1),[3 2 1]),CL,'uni',false),CL,'uni',false), ...
    analysis.epoch_activity_CL,'uni',false);

CR_avg_activity = cellfun(@(CR) cellfun(@(CR) cellfun(@(CR) ...
    permute(nanmean(CR,1),[3 2 1]),CR,'uni',false),CR,'uni',false), ...
    analysis.epoch_activity_CR,'uni',false);

CL_CR_corr = cell(size(CL_avg_activity));
for curr_animal = 1:length(mice)
    for curr_session = 1:length(CL_avg_activity{curr_animal})
        for curr_rev = 1:length(CL_avg_activity{curr_animal}{curr_session})
            CL_zscore = reshape(zscore( ...
                CL_avg_activity{curr_animal}{curr_session}{curr_rev}'),[],1);
            CR_zscore = reshape(zscore( ...
                CR_avg_activity{curr_animal}{curr_session}{curr_rev}'),[],1);
            
            nan_points = any(isnan([CL_zscore CR_zscore]),2);
            
            CL_zscore_nonan = CL_zscore(~nan_points);
            CR_zscore_nonan = CR_zscore(~nan_points);
            
            curr_corr = corrcoef(CL_zscore_nonan,CR_zscore_nonan);
            CL_CR_corr{curr_animal}{curr_session}(curr_rev) = curr_corr(2);
        end
    end
end



%% Testing some visualization / PCA stuff

a = permute(analysis.epoch_activity_CL{8}{5}{1},[2 1 3]);
b = permute(analysis.epoch_activity_CR{8}{5}{1},[2 1 3]);

a2 = reshape(a,size(a,1)*size(a,2),size(a,3));
b2 = reshape(b,size(b,1)*size(b,2),size(b,3));

c = zscore([a2;b2]);

[coeff score latent] = princomp(c);

c1 = reshape(score(1:size(a2,1),1:3),size(a,1),size(a,2),size(a,3));
c2 = reshape(score(size(a2,1)+1:end,1:3),size(a,1),size(a,2),size(a,3));

%%%%%%%%%
a = permute(nanmean(analysis.epoch_activity_CL{8}{10}{1},1),[2 3 1]);
b = permute(nanmean(analysis.epoch_activity_CR{8}{10}{1},1),[2 3 1]);
a2 = permute(nanmean(analysis.epoch_activity_CL{8}{10}{2},1),[2 3 1]);
b2 = permute(nanmean(analysis.epoch_activity_CR{8}{10}{2},1),[2 3 1]);

c = zscore([a;a2;b;b2]);

[coeff score latent] = princomp(c);

d = permute(reshape(score(:,1:3)',3,[],4),[2 1 3]);

figure; hold on;
plot3(d(:,1,1),d(:,2,1),d(:,3,1),'k','linewidth',2)
plot3(d(:,1,2),d(:,2,2),d(:,3,2),'r','linewidth',2)
plot3(d(:,1,3),d(:,2,3),d(:,3,3),'k')
plot3(d(:,1,4),d(:,2,4),d(:,3,4),'r')


%% Slide 2-4) Scatter/barplots of reversal up/down modulation

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

CL_avg_activity = cellfun(@(act) cellfun(@(act) cellfun(@(act) ...
    permute(nanmean(act,2),[1 3 2]),act,'uni',false),act,'uni',false), ...
    analysis.epoch_activity_CL,'uni',false);

CR_avg_activity = cellfun(@(act) cellfun(@(act) cellfun(@(act) ...
    permute(nanmean(act,2),[1 3 2]),act,'uni',false),act,'uni',false), ...
    analysis.epoch_activity_CR,'uni',false);

first_reversals = nan(size(mice));
last_reversals = nan(size(mice));

first_reversals(rev_animals) = cellfun(@(x) find(cellfun(@(x) length(x) > 1,x),1),CL_avg_activity(rev_animals));
last_reversals(rev_animals) = cellfun(@(x) find(cellfun(@(x) length(x) > 1,x),1,'last'),CL_avg_activity(rev_animals));

first_reversals(ctrl_animals) = 6;
last_reversals(ctrl_animals) = cellfun(@(x) find(cellfun(@(x) length(x) == 1,x),1,'last'),CL_avg_activity(ctrl_animals));

odor_sigdiff = cellfun(@(odor,baseline) cellfun(@(odor,baseline) cellfun(@(odor,baseline) ...
    AP_signrank_matrix( ...
    permute(nanmean(odor,2),[1 3 2]),permute(nanmean(baseline,2),[1 3 2])) < 0.05 & ...
    nanmean(permute(any(odor,2),[1 3 2]),1) > 0.1, ...
    odor,baseline,'uni',false),odor,baseline,'uni',false), ...
    analysis.epoch_activity_odor,analysis.epoch_activity_baseline,...
    'uni',false);

% First reversal

% CL
CL_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CL, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = first_reversals(curr_animal);
    if any(intersect(curr_animal,find(rev_animals)))
        pre = CL_avg_activity{curr_animal}{curr_session}{1};
        post = CL_avg_activity{curr_animal}{curr_session}{2};
        col = 'b';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1} | ...
        odor_sigdiff{curr_animal}{curr_session}{2};
    else
        num_trials = size(CL_avg_activity{curr_animal}{curr_session}{1},1);
        half_trials = ceil(num_trials/2);
        pre = CL_avg_activity{curr_animal}{curr_session}{1}(1:half_trials,:);
        post = CL_avg_activity{curr_animal}{curr_session}{1}(half_trials+1:end,:);
        col = 'm';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1};
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CL_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
    
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% CR
CR_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CR, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = first_reversals(curr_animal);
    if any(intersect(curr_animal,find(rev_animals)))
        pre = CR_avg_activity{curr_animal}{curr_session}{1};
        post = CR_avg_activity{curr_animal}{curr_session}{2};
        col = 'b';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1} | ...
        odor_sigdiff{curr_animal}{curr_session}{2};
    else
        num_trials = size(CR_avg_activity{curr_animal}{curr_session}{1},1);
        half_trials = ceil(num_trials/2);
        pre = CR_avg_activity{curr_animal}{curr_session}{1}(1:half_trials,:);
        post = CR_avg_activity{curr_animal}{curr_session}{1}(half_trials+1:end,:);
        col = 'm';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1};
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CR_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
        
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% Bar plot
% CL
frac_up_rev_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CL frac_up_rev_alm_CL frac_up_ctrl_pmm_CL frac_up_ctrl_alm_CL;...
             frac_down_rev_pmm_CL frac_down_rev_alm_CL frac_down_ctrl_pmm_CL frac_down_ctrl_alm_CL};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    std(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CL' 'Rev ALM CL' 'Ctrl PMM CL' 'Ctrl ALM CL'});
ylabel('Fraction of significantly modulated cells');
title('First reversal')

% CR
frac_up_rev_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CR frac_up_rev_alm_CR frac_up_ctrl_pmm_CR frac_up_ctrl_alm_CR;...
             frac_down_rev_pmm_CR frac_down_rev_alm_CR frac_down_ctrl_pmm_CR frac_down_ctrl_alm_CR};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    std(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CR' 'Rev ALM CR' 'Ctrl PMM CR' 'Ctrl ALM CR'});
ylabel('Fraction of significantly modulated cells');
title('First reversal')

% Last reversal

% Make exceptions for particular animals with messed up data :(
last_reversals(6) = 12;
last_reversals(7) = NaN; % no other within-session reversals

% CL
CL_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CL, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    if curr_animal == 7;
        continue
    end
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = last_reversals(curr_animal);
    if any(intersect(curr_animal,find(rev_animals)))
        pre = CL_avg_activity{curr_animal}{curr_session}{1};
        post = CL_avg_activity{curr_animal}{curr_session}{2};
        col = 'b';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1} | ...
        odor_sigdiff{curr_animal}{curr_session}{2};
    else
        num_trials = size(CL_avg_activity{curr_animal}{curr_session}{1},1);
        half_trials = ceil(num_trials/2);
        pre = CL_avg_activity{curr_animal}{curr_session}{1}(1:half_trials,:);
        post = CL_avg_activity{curr_animal}{curr_session}{1}(half_trials+1:end,:);
        col = 'm';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1};
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CL_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
    
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% CR
CR_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CR, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    if curr_animal == 7
        continue
    end
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = last_reversals(curr_animal);
    if any(intersect(curr_animal,find(rev_animals)))
        pre = CR_avg_activity{curr_animal}{curr_session}{1};
        post = CR_avg_activity{curr_animal}{curr_session}{2};
        col = 'b';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1} | ...
        odor_sigdiff{curr_animal}{curr_session}{2};
    else
        num_trials = size(CR_avg_activity{curr_animal}{curr_session}{1},1);
        half_trials = ceil(num_trials/2);
        pre = CR_avg_activity{curr_animal}{curr_session}{1}(1:half_trials,:);
        post = CR_avg_activity{curr_animal}{curr_session}{1}(half_trials+1:end,:);
        col = 'm';
        curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session}{1};
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CR_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
        
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% Bar plot
% CL
frac_up_rev_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CL frac_up_rev_alm_CL frac_up_ctrl_pmm_CL frac_up_ctrl_alm_CL;...
             frac_down_rev_pmm_CL frac_down_rev_alm_CL frac_down_ctrl_pmm_CL frac_down_ctrl_alm_CL};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    nanstd(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CL' 'Rev ALM CL' 'Ctrl PMM CL' 'Ctrl ALM CL'});
ylabel('Fraction of significantly modulated cells');
title('Last reversal')

% CR
frac_up_rev_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CR frac_up_rev_alm_CR frac_up_ctrl_pmm_CR frac_up_ctrl_alm_CR;...
             frac_down_rev_pmm_CR frac_down_rev_alm_CR frac_down_ctrl_pmm_CR frac_down_ctrl_alm_CR};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    nanstd(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CR' 'Rev ALM CR' 'Ctrl PMM CR' 'Ctrl ALM CR'});
ylabel('Fraction of significantly modulated cells');
title('Last reversal')


% number of CR/CL trials for each animal pre and post reversal
figure;set(gcf,'Name','Blue = CL,Red = CR');
CL_CR_numtrials = cell(size(mice));
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = first_reversals(curr_animal);
    if any(intersect(curr_animal,find(rev_animals)))
        CL_pre = size(CL_avg_activity{curr_animal}{curr_session}{1},1);
        CL_post = size(CL_avg_activity{curr_animal}{curr_session}{2},1);
        
        CR_pre = size(CR_avg_activity{curr_animal}{curr_session}{1},1);
        CR_post = size(CR_avg_activity{curr_animal}{curr_session}{2},1);
    else
        num_trials = size(CL_avg_activity{curr_animal}{curr_session}{1},1);
        CL_pre = ceil(num_trials/2);
        CL_post = floor(num_trials/2);
        
        num_trials = size(CR_avg_activity{curr_animal}{curr_session}{1},1);
        CR_pre = ceil(num_trials/2);
        CR_post = floor(num_trials/2);
    end

    bar([CL_pre CR_pre;CL_pre CL_post]);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'Pre' 'Post'});
end

%% Slide 5) Scatter/barplots of reversal up/down modulation, but using pre/post sessions for reversal

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

CL_avg_activity = cellfun(@(act) cellfun(@(act) cellfun(@(act) ...
    permute(nanmean(act,2),[1 3 2]),act,'uni',false),act,'uni',false), ...
    analysis.epoch_activity_CL,'uni',false);

CR_avg_activity = cellfun(@(act) cellfun(@(act) cellfun(@(act) ...
    permute(nanmean(act,2),[1 3 2]),act,'uni',false),act,'uni',false), ...
    analysis.epoch_activity_CR,'uni',false);

first_reversals = nan(size(mice));
last_reversals = nan(size(mice));

first_reversals(rev_animals) = cellfun(@(x) find(cellfun(@(x) length(x) > 1,x),1),CL_avg_activity(rev_animals));
last_reversals(rev_animals) = cellfun(@(x) find(cellfun(@(x) length(x) > 1,x),1,'last'),CL_avg_activity(rev_animals));

first_reversals(ctrl_animals) = 6;
last_reversals(ctrl_animals) = cellfun(@(x) find(cellfun(@(x) length(x) == 1,x),1,'last'),CL_avg_activity(ctrl_animals));

odor_sigdiff = cellfun(@(odor,baseline) cellfun(@(odor,baseline) cellfun(@(odor,baseline) ...
    AP_signrank_matrix( ...
    permute(nanmean(odor,2),[1 3 2]),permute(nanmean(baseline,2),[1 3 2])) < 0.05 & ...
    nanmean(permute(any(odor,2),[1 3 2]),1) > 0.1, ...
    odor,baseline,'uni',false),odor,baseline,'uni',false), ...
    analysis.epoch_activity_odor,analysis.epoch_activity_baseline,...
    'uni',false);

% First reversal

% CL
CL_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CL, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = first_reversals(curr_animal);
    % sometimes the next day starts out wrong contingency, make sure it's 2
    if any(intersect(curr_animal,find(rev_animals)))
        post_rev = find(analysis.epoch_contingencies{curr_animal}{curr_session+1} == 2);
    else
        post_rev = 1;
    end
    
    pre = CL_avg_activity{curr_animal}{curr_session-1}{end};
    post = CL_avg_activity{curr_animal}{curr_session+1}{post_rev};
    
    curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session+1}{end} | ...
        odor_sigdiff{curr_animal}{curr_session+1}{post_rev};
    
    if any(intersect(curr_animal,find(rev_animals)))
        col = 'b';
    else
        col = 'm';
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CL_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
    
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% CR
CR_first_rev_mod = cell(size(mice));
figure;set(gcf,'Name','CR, b = rev, -- = pmm')
sq_animals = ceil(sqrt(length(mice)));
for curr_animal = 1:length(mice)
    subplot(sq_animals,sq_animals,curr_animal);hold on;
    
    curr_session = first_reversals(curr_animal);
    % sometimes the next day starts out wrong contingency, make sure it's 2
    if any(intersect(curr_animal,find(rev_animals)))
        post_rev = find(analysis.epoch_contingencies{curr_animal}{curr_session+1} == 2);
    else
        post_rev = 1;
    end
    
    pre = CR_avg_activity{curr_animal}{curr_session-1}{end};
    post = CR_avg_activity{curr_animal}{curr_session+1}{post_rev};
    
    curr_odor_sigdiff = odor_sigdiff{curr_animal}{curr_session+1}{end} | ...
        odor_sigdiff{curr_animal}{curr_session+1}{post_rev};
    
    if any(intersect(curr_animal,find(rev_animals)))
        col = 'b';
    else
        col = 'm';
    end
    if any(intersect(curr_animal,find(pmm_animals)))
        ls = '--';
    else
        ls = '-';
    end
    
    sig_cells = AP_ranksum_matrix(pre,post) < 0.05 & curr_odor_sigdiff;
    CR_first_rev_mod{curr_animal} = (1*(nanmean(pre(:,curr_odor_sigdiff)) ...
        < nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff))) + ...
        (-1*(nanmean(pre(:,curr_odor_sigdiff)) > ...
        nanmean(post(:,curr_odor_sigdiff)) & sig_cells(curr_odor_sigdiff)));
        
    plot(nanmean(pre(:,~sig_cells)),nanmean(post(:,~sig_cells)),'.k','MarkerSize',5);
    plot(nanmean(pre(:,sig_cells)),nanmean(post(:,sig_cells)),'.r','MarkerSize',5);
    lim = [min([nanmean(pre) nanmean(post)]) max([nanmean(pre) nanmean(post)])];
    xlim(lim);
    ylim(lim);
    line(lim,lim,'color',col,'linestyle',ls);
end

% Bar plot
% CL
frac_up_rev_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CL = cellfun(@(x) nanmean(x == 1),CL_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CL = cellfun(@(x) nanmean(x == -1),CL_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CL frac_up_rev_alm_CL frac_up_ctrl_pmm_CL frac_up_ctrl_alm_CL;...
             frac_down_rev_pmm_CL frac_down_rev_alm_CL frac_down_ctrl_pmm_CL frac_down_ctrl_alm_CL};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    nanstd(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CL' 'Rev ALM CL' 'Ctrl PMM CL' 'Ctrl ALM CL'});
ylabel('Fraction of significantly modulated cells');
title('First reversal')

% CR
frac_up_rev_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_down_rev_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & pmm_animals));
frac_up_rev_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(rev_animals & alm_animals));
frac_down_rev_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(rev_animals & alm_animals));

frac_up_ctrl_pmm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_down_ctrl_pmm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & pmm_animals));
frac_up_ctrl_alm_CR = cellfun(@(x) nanmean(x == 1),CR_first_rev_mod(ctrl_animals & alm_animals));
frac_down_ctrl_alm_CR = cellfun(@(x) nanmean(x == -1),CR_first_rev_mod(ctrl_animals & alm_animals));

plot_data = {frac_up_rev_pmm_CR frac_up_rev_alm_CR frac_up_ctrl_pmm_CR frac_up_ctrl_alm_CR;...
             frac_down_rev_pmm_CR frac_down_rev_alm_CR frac_down_ctrl_pmm_CR frac_down_ctrl_alm_CR};

figure;
errorb(cellfun(@nanmean,plot_data)',cellfun(@(x) ...
    nanstd(x)./sqrt(sum(~isnan(x))),plot_data)');
legend({'Up' 'Down'},'location','NW');
set(gca,'XTickLabel',{'Rev PMM CR' 'Rev ALM CR' 'Ctrl PMM CR' 'Ctrl ALM CR'});
ylabel('Fraction of significantly modulated cells');
title('First reversal')


%% NOTE ON ABOVE!!! "BASELINE" NO LONGER SEPERATE VARIABLE


%% Slide 6-39) Get significantly modulated cells by epoch/trial type

% by request of TK/Jun: to find epoch with the most significant cells? or
% something related to finding a threshold? 
% Either way: compare baseline to each of 0-4 seconds in 1 second periods
% (signrank) or compare CL to CR trials (ranksum)

% NOTE: This requires a baseline of 1 sec and an epoch of 0-4 sec

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

epochs = 4;
baseline_multiepoch_sigdiff = cell(epochs,1);
condition_multiepoch_sigdiff = cell(epochs,1);
for curr_epoch = 1:epochs

    baseline_multiepoch_sigdiff{curr_epoch} = ...
        cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
        AP_signrank_matrix( ...
        permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
        permute(nanmean(odor(:,frames(1):frames(2) > (curr_epoch-1)*frames(2)/epochs & ...
        frames(1):frames(2) < curr_epoch*frames(2)/epochs,:),2),[1 3 2])) ...
        < 0.01/epochs, ...
        odor,'uni',false),odor,frames,'uni',false), ...
        analysis.epoch_activity_odor,analysis.epoch_frames,...
        'uni',false);
    
    condition_multiepoch_sigdiff{curr_epoch} = ...
        cellfun(@(CL,CR,frames) cellfun(@(CL,CR,frames) cellfun(@(CL,CR) ...
        AP_ranksum_matrix( ...
        permute(nanmean(CL(:,frames(1):frames(2) > (curr_epoch-1)*frames(2)/epochs & ...
        frames(1):frames(2) < curr_epoch*frames(2)/epochs,:),2),[1 3 2]), ...
        permute(nanmean(CR(:,frames(1):frames(2) > (curr_epoch-1)*frames(2)/epochs & ...
        frames(1):frames(2) < curr_epoch*frames(2)/epochs,:),2),[1 3 2])) ...
        < 0.01/epochs, ...
        CL,CR,'uni',false),CL,CR,frames,'uni',false), ...
        analysis.epoch_activity_CL,analysis.epoch_activity_CR,analysis.epoch_frames,...
        'uni',false);
    
    disp(['Epoch: ' num2str(curr_epoch)]);
end

% Plot total fraction of significant cells for each epoch for each measure
baseline_multiepoch_sigdiff_animals = cellfun(@(x) horzcat(x{:}),baseline_multiepoch_sigdiff,'uni',false);
baseline_multiepoch_sigdiff_sessions = cellfun(@(x) vertcat(x{:}),baseline_multiepoch_sigdiff_animals,'uni',false);
baseline_multiepoch_sigdiff_revs = cellfun(@(x) horzcat(x{:}),baseline_multiepoch_sigdiff_sessions,'uni',false);
baseline_multiepoch_sigdiff_totalmean = cellfun(@nanmean,baseline_multiepoch_sigdiff_revs);

condition_multiepoch_sigdiff_animals = cellfun(@(x) horzcat(x{:}),condition_multiepoch_sigdiff,'uni',false);
condition_multiepoch_sigdiff_sessions = cellfun(@(x) vertcat(x{:}),condition_multiepoch_sigdiff_animals,'uni',false);
condition_multiepoch_sigdiff_revs = cellfun(@(x) horzcat(x{:}),condition_multiepoch_sigdiff_sessions,'uni',false);
condition_multiepoch_sigdiff_totalmean = cellfun(@nanmean,condition_multiepoch_sigdiff_revs);

overlap_multiepoch_sigdiff_totalmean = cellfun(@(x,y) nanmean(x&y), ...
    baseline_multiepoch_sigdiff_revs,condition_multiepoch_sigdiff_revs);

figure;
bar([baseline_multiepoch_sigdiff_totalmean ...
    condition_multiepoch_sigdiff_totalmean overlap_multiepoch_sigdiff_totalmean]);
legend({'Baseline' 'Condition' 'Both'},'location','nw')
set(gca,'XTickLabel',{'0-1s' '1-2s' '2-3s' '3-4s'});
ylabel('Fraction of total significant cells / all cells');

% plot first 36 significant cells in session 3 in both measures all animals
use_session = 3;
for curr_animal = 1:length(mice)
    
   curr_fig = figure('Name',mice(curr_animal).name);  
   
   curr_baseline_sigdiff = cell2mat(cellfun(@(x) ...
       x{curr_animal}{use_session}{1},baseline_multiepoch_sigdiff,'uni',false));
   
   curr_condition_sigdiff = cell2mat(cellfun(@(x) ...
       x{curr_animal}{use_session}{1},condition_multiepoch_sigdiff,'uni',false));
   
   curr_sig_cells = find(any([curr_baseline_sigdiff;curr_condition_sigdiff],1),36);
   for curr_sig_cell = 1:length(curr_sig_cells);
       
       curr_cell = curr_sig_cells(curr_sig_cell);
       subplot(6,6,curr_sig_cell); hold on;
       
       curr_CL = analysis.epoch_activity_CL{curr_animal}{use_session}{1}(:,:,curr_cell);
       curr_CR = analysis.epoch_activity_CR{curr_animal}{use_session}{1}(:,:,curr_cell);
       
       plot(nanmean(curr_CL),'k','linewidth',2);
       plot(nanmean(curr_CL)+nanstd(curr_CL)./sqrt(sum(~isnan(curr_CL))),'--k');
       plot(nanmean(curr_CL)-nanstd(curr_CL)./sqrt(sum(~isnan(curr_CL))),'--k');
       
       plot(nanmean(curr_CR),'r','linewidth',2);
       plot(nanmean(curr_CR)+nanstd(curr_CR)./sqrt(sum(~isnan(curr_CR))),'--r');
       plot(nanmean(curr_CR)-nanstd(curr_CR)./sqrt(sum(~isnan(curr_CR))),'--r');
       
       sig_linestyle = cell(epochs,1);
       sig_linestyle(curr_baseline_sigdiff(:,curr_cell)) = {'--'};
       sig_linestyle(~curr_baseline_sigdiff(:,curr_cell)) = {'-'};
       
       sig_color = cell(epochs,1);
       sig_color(curr_condition_sigdiff(:,curr_cell)) = {'g'};
       sig_color(~curr_condition_sigdiff(:,curr_cell)) = {'b'};       
       
       curr_framelimit = analysis.epoch_frames{curr_animal}{use_session};
       curr_frames = curr_framelimit(1):curr_framelimit(2);
       for curr_epoch = 1:epochs
           line(repmat(find(curr_frames == ceil((curr_epoch-1)* ...
               curr_framelimit(2)/epochs)),2,1),ylim, ...
               'linestyle',sig_linestyle{curr_epoch}, ...
               'color',sig_color{curr_epoch});
       end
       
       set(gca,'XTick',[]);
       title(curr_cell);
       
   end
      
   saveas(curr_fig,['/usr/local/lab/People/Andy/Jun_lickreversal/140903_analysis/example_activity/example_sigcells_animal_' num2str(curr_animal) 'NEWTHRESH.fig'])
   close(curr_fig);
   
end

























