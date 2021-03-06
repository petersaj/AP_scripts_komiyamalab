%% Trial-trial activity correlation of task-related cells

rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Get "task-related" cells by those which are active in significantly
% (just use chi square instead of McNemar statistic which is the right one
% for paired data)

baseline_sigcells = ...
        cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
        AP_chisquare_matrix( ...
        permute(any(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
        permute(any(odor(:,frames(1):frames(2) > 0 & ...
        frames(1):frames(2) < 2,:),2),[1 3 2])) ...
        < 0.05, odor,'uni',false),odor,frames,'uni',false), ...
        analysis.epoch_activity_odor,analysis.epoch_frames,...
        'uni',false);
    
prereversal_sig_CL = ...
    cellfun(@(act,sig) cellfun(@(act,sig) act{1}(:,:,sig{1}),act(1:10),sig(1:10),'uni',false), ...
    analysis.epoch_activity_CL,baseline_sigcells,'uni',false);

prereversal_corr = cellfun(@(x) cellfun(@(y) ...
    AP_itril(corrcoef(reshape(permute(y,[2 3 1]),[],size(y,1))),-1),x, ...
    'uni',false),prereversal_sig_CL,'uni',false);



%% Pull out early-active cells

% Get significantly modulated cells

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Go through each epoch, compare 1) to baseline, 2) across conditions
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
        < 0.01 & ...
        nanmean(permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
        nanmean(permute(nanmean(odor(:,frames(1):frames(2) > (curr_epoch-1)*frames(2)/epochs & ...
        frames(1):frames(2) < curr_epoch*frames(2)/epochs,:),2),[1 3 2]),1), ...
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


early_active_cells = cellfun(@(e1,e2,e3,e4) cellfun(@(e1,e2,e3,e4) ...
    cellfun(@(e1,e2,e3,e4) ...
    (e1|e2) & ~e4, ... % The active epoch contstraints
    e1,e2,e3,e4,'uni',false), ...
    e1,e2,e3,e4,'uni',false),baseline_multiepoch_sigdiff{1}, ... 
    baseline_multiepoch_sigdiff{2}, baseline_multiepoch_sigdiff{3}, ...
    baseline_multiepoch_sigdiff{4},'uni',false);


 
%% Plot the early active cells 

use_session = 5;

for curr_animal = 8;
    
   curr_fig = figure('Name',mice(curr_animal).name);  
   
   curr_baseline_sigdiff = cell2mat(cellfun(@(x) ...
       x{curr_animal}{use_session}{1},baseline_multiepoch_sigdiff,'uni',false));
   
   curr_condition_sigdiff = cell2mat(cellfun(@(x) ...
       x{curr_animal}{use_session}{1},condition_multiepoch_sigdiff,'uni',false));
   
   curr_sig_cells = find(early_active_cells{curr_animal}{use_session}{1});
   
   for curr_sig_cell = 1:length(curr_sig_cells);
       
       sqr_plot = ceil(sqrt(length(curr_sig_cells)));
       
       curr_cell = curr_sig_cells(curr_sig_cell);
       subplot(sqr_plot,sqr_plot,curr_sig_cell); hold on;
       
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
      
   %saveas(curr_fig,['/usr/local/lab/People/Andy/Jun_lickreversal/141007_analysis/early_sigcells/early_sigcells_animal_' num2str(curr_animal) '_thresh4.fig'])
   %close(curr_fig);
   
end




%% Slide 1) Get significance during trial times, get overlap between trial times


% Get significantly modulated cells

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% Go through each epoch, compare 1) to baseline, 2) across conditions
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
        < 0.01 & ...
        nanmean(permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
        nanmean(permute(nanmean(odor(:,frames(1):frames(2) > (curr_epoch-1)*frames(2)/epochs & ...
        frames(1):frames(2) < curr_epoch*frames(2)/epochs,:),2),[1 3 2]),1), ...
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

baseline_sig_epoch_overlap_1 = [nanmean(vertcat(baseline_multiepoch_sigdiff_revs{1:3}) & ...
    vertcat(baseline_multiepoch_sigdiff_revs{2:4}),2);0];
baseline_sig_epoch_overlap_2 = [nanmean(vertcat(baseline_multiepoch_sigdiff_revs{1:2}) & ...
    vertcat(baseline_multiepoch_sigdiff_revs{3:4}),2);0;0];
baseline_sig_epoch_overlap_3 = [nanmean(vertcat(baseline_multiepoch_sigdiff_revs{1}) & ...
    vertcat(baseline_multiepoch_sigdiff_revs{4}),2);0;0;0];

condition_multiepoch_sigdiff_animals = cellfun(@(x) horzcat(x{:}),condition_multiepoch_sigdiff,'uni',false);
condition_multiepoch_sigdiff_sessions = cellfun(@(x) vertcat(x{:}),condition_multiepoch_sigdiff_animals,'uni',false);
condition_multiepoch_sigdiff_revs = cellfun(@(x) horzcat(x{:}),condition_multiepoch_sigdiff_sessions,'uni',false);
condition_multiepoch_sigdiff_totalmean = cellfun(@nanmean,condition_multiepoch_sigdiff_revs);

condition_sig_epoch_overlap_1 = [nanmean(vertcat(condition_multiepoch_sigdiff_revs{1:3}) & ...
    vertcat(condition_multiepoch_sigdiff_revs{2:4}),2);0];
condition_sig_epoch_overlap_2 = [nanmean(vertcat(condition_multiepoch_sigdiff_revs{1:2}) & ...
    vertcat(condition_multiepoch_sigdiff_revs{3:4}),2);0;0];
condition_sig_epoch_overlap_3 = [nanmean(vertcat(condition_multiepoch_sigdiff_revs{1}) & ...
    vertcat(condition_multiepoch_sigdiff_revs{4}),2);0;0;0];

overlap_multiepoch_sigdiff_totalmean = cellfun(@(x,y) nanmean(x&y), ...
    baseline_multiepoch_sigdiff_revs,condition_multiepoch_sigdiff_revs);

figure;
subplot(1,2,1);
bar([baseline_multiepoch_sigdiff_totalmean baseline_sig_epoch_overlap_1 ...
    baseline_sig_epoch_overlap_2 baseline_sig_epoch_overlap_3],'grouped');
legend({'Baseline' 'Epoch+1' 'Epoch+2' 'Epoch+3'},'location','se')
set(gca,'XTickLabel',{'0-1s' '1-2s' '2-3s' '3-4s'});
ylabel('Fraction of total significant cells / all cells');

subplot(1,2,2);
bar([condition_multiepoch_sigdiff_totalmean condition_sig_epoch_overlap_1 ...
    condition_sig_epoch_overlap_2 condition_sig_epoch_overlap_3],'grouped');
legend({'Condition' 'Epoch+1' 'Epoch+2' 'Epoch+3'},'location','se')
set(gca,'XTickLabel',{'0-1s' '1-2s' '2-3s' '3-4s'});
ylabel('Fraction of total significant cells / all cells');




%% Slide 2-3) Average correlation of population activity across days

avg_pop_corr_cr = cell(size(mice));
avg_pop_corr_cl = cell(size(mice));

for curr_animal = 1:length(mice)
    
    curr_mean_cl_act = [];
    curr_mean_cr_act = [];
    
    for curr_session = 1:length(analysis.epoch_activity_CL{curr_animal})        
        for curr_contingency = 1:length(analysis.epoch_activity_CL{curr_animal}{curr_session});
            
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1): ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_frames_use = curr_frames > 0 & curr_frames < 2*(max(curr_frames)/4);
            
            curr_mean_cl_act = [curr_mean_cl_act ...
                reshape(permute(nanmean( ...
                analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_contingency}(:, ...
                curr_frames_use,:),1),[2 3 1]),[],1)];
            
            curr_mean_cr_act = [curr_mean_cr_act ...
                reshape(permute(nanmean( ...
                analysis.epoch_activity_CR{curr_animal}{curr_session}{curr_contingency}(:, ...
                curr_frames_use,:),1),[2 3 1]),[],1)];
                   
        end
    end 
    
    curr_contingencies = vertcat(analysis.epoch_contingencies{curr_animal}{:});
    curr_sessions = vertcat(analysis.epoch_sessions{curr_animal}{:});
    
    % Group trial type within day across contingencies
    curr_mean_cl_act_grp = nan(size(curr_mean_cl_act,1),max(curr_sessions));
    curr_mean_cr_act_grp = nan(size(curr_mean_cl_act,1),max(curr_sessions));
    
    curr_mean_cl_act_grp(:,unique(curr_sessions)) = grpstats(curr_mean_cl_act',curr_sessions)';
    curr_mean_cr_act_grp(:,unique(curr_sessions)) = grpstats(curr_mean_cr_act',curr_sessions)';
    
    % Correlate activity
    avg_pop_corr_cl{curr_animal} = corrcoef(curr_mean_cl_act_grp,'rows','pairwise');
    avg_pop_corr_cr{curr_animal} = corrcoef(curr_mean_cr_act_grp,'rows','pairwise');
    
    disp(curr_animal)
end

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% pad activity correlations with nans

max_sessions = max(cellfun(@length,avg_pop_corr_cl));
avg_pop_corr_cl_pad = cellfun(@(x) padarray(x,repmat(max_sessions-length(x),2,1),nan,'post'), ...
    avg_pop_corr_cl,'uni',false);
avg_pop_corr_cr_pad = cellfun(@(x) padarray(x,repmat(max_sessions-length(x),2,1),nan,'post'), ...
    avg_pop_corr_cr,'uni',false);

cl_cat = cat(3,avg_pop_corr_cl_pad{:});
cr_cat = cat(3,avg_pop_corr_cr_pad{:});

rev_pmm_cl = nanmean(cl_cat(:,:,rev_animals & pmm_animals),3);
rev_alm_cl = nanmean(cl_cat(:,:,rev_animals & alm_animals),3);
ctrl_pmm_cl = nanmean(cl_cat(:,:,ctrl_animals & pmm_animals),3);
ctrl_alm_cl = nanmean(cl_cat(:,:,ctrl_animals & alm_animals),3);

rev_pmm_cr = nanmean(cr_cat(:,:,rev_animals & pmm_animals),3);
rev_alm_cr = nanmean(cr_cat(:,:,rev_animals & alm_animals),3);
ctrl_pmm_cr = nanmean(cr_cat(:,:,ctrl_animals & pmm_animals),3);
ctrl_alm_cr = nanmean(cr_cat(:,:,ctrl_animals & alm_animals),3);

% Plot CL
figure;

subplot(2,2,1);
imagesc(rev_pmm_cl);colormap(gray);
title('rev pmm cl')

subplot(2,2,2);
imagesc(rev_alm_cl);colormap(gray);
title('rev alm cl')

subplot(2,2,3);
imagesc(ctrl_pmm_cl);colormap(gray);
title('ctrl pmm cl')

subplot(2,2,4);
imagesc(ctrl_alm_cl);colormap(gray);
title('ctrl alm cl')

% Plot the first diagonal of the CL correlation
cl_cat_diag = cell2mat(arrayfun(@(x) diag(cl_cat(:,:,x),-1),1:size(cl_cat,3),'uni',false))';
figure; hold on;
% errorbar(nanmean(cl_cat_diag(rev_animals & pmm_animals,:)), ...
%     nanstd(cl_cat_diag(rev_animals & pmm_animals,:))./ ...
%     sqrt(sum(~isnan(cl_cat_diag(rev_animals & pmm_animals,:)))), ...
%     'k','linewidth',2,'linestyle','-')
% errorbar(nanmean(cl_cat_diag(ctrl_animals & pmm_animals,:)), ...
%     nanstd(cl_cat_diag(rev_animals & pmm_animals,:))./ ...
%     sqrt(sum(~isnan(cl_cat_diag(rev_animals & pmm_animals,:)))), ...
%     'k','linewidth',2,'linestyle','--')
% errorbar(nanmean(cl_cat_diag(rev_animals & alm_animals,:)), ...
%     nanstd(cl_cat_diag(rev_animals & pmm_animals,:))./ ...
%     sqrt(sum(~isnan(cl_cat_diag(rev_animals & pmm_animals,:)))), ...
%     'r','linewidth',2,'linestyle','-')
% errorbar(nanmean(cl_cat_diag(ctrl_animals & alm_animals,:)), ...
%     nanstd(cl_cat_diag(rev_animals & pmm_animals,:))./ ...
%     sqrt(sum(~isnan(cl_cat_diag(rev_animals & pmm_animals,:)))), ...
%     'r','linewidth',2,'linestyle','--')
plot(nanmean(cl_cat_diag(rev_animals & pmm_animals,:)), ...
    'k','linewidth',2,'linestyle','-')
plot(nanmean(cl_cat_diag(ctrl_animals & pmm_animals,:)), ...
    'k','linewidth',2,'linestyle','--')
plot(nanmean(cl_cat_diag(rev_animals & alm_animals,:)), ...
    'r','linewidth',2,'linestyle','-')
plot(nanmean(cl_cat_diag(ctrl_animals & alm_animals,:)), ...
    'r','linewidth',2,'linestyle','--')

ylabel('Correlation');
xlabel('Session N with Session N+1')
legend({'rev pmm' 'ctrl pmm' 'rev alm' 'ctrl alm'});

% Plot CR
figure;

subplot(2,2,1);
imagesc(rev_pmm_cr);colormap(gray);
title('rev pmm cr')

subplot(2,2,2);
imagesc(rev_alm_cr);colormap(gray);
title('rev alm cr')

subplot(2,2,3);
imagesc(ctrl_pmm_cr);colormap(gray);
title('ctrl pmm cr')

subplot(2,2,4);
imagesc(ctrl_alm_cr);colormap(gray);
title('ctrl alm cr')



%% Visualize average activity of all cells / sessions in given animal



curr_animal = 1;
for curr_session = 1:length(analysis.epoch_activity_CL{curr_animal});
    % plot only cells sigdiff from baseline in epochs 1 or 2
    %plot_cells = baseline_multiepoch_sigdiff{1}{curr_animal}{curr_session}{1} | ...
    %   baseline_multiepoch_sigdiff{1}{curr_animal}{curr_session}{1};
    
    curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1): ...
        analysis.epoch_frames{curr_animal}{curr_session}(2);
    curr_frames_use = true(size(curr_frames));
    
    plot_cells = setdiff(1:size( ...
        analysis.epoch_activity_CL{curr_animal}{curr_session}{1},3), ...
        mice(curr_animal).interneurons);
    
    figure;imagesc(permute(nanmean( ...
        analysis.epoch_activity_CL{curr_animal}{curr_session}{1} ...
        (:,curr_frames_use,plot_cells),1), ...
        [3 2 1]));colormap(hot);caxis([0 0.5])
    
    line(repmat(find(curr_frames == 0),2,1),ylim,'color','b','linewidth',2);
    line(repmat(find(curr_frames == round(max(curr_frames)/2)),2,1), ...
        ylim,'color','b','linewidth',2);
    
    figure;imagesc(permute(nanmean( ...
        analysis.epoch_activity_CR{curr_animal}{curr_session}{1} ...
        (:,curr_frames_use,plot_cells),1), ...
        [3 2 1]));colormap(hot);caxis([0 0.5])
    
    line(repmat(find(curr_frames == 0),2,1),ylim,'color','b','linewidth',2);
    line(repmat(find(curr_frames == round(max(curr_frames)/2)),2,1), ...
        ylim,'color','b','linewidth',2);
end



%% Get lag of activity over sessions / timing changes of population?

% Can't concatenate activity to get lag, can't overlap between cells

xcorr_lags = cell(length(mice),1);
for curr_animal = 1:length(mice)
    
    plot_cells = setdiff(1:size( ...
        analysis.epoch_activity_CL{curr_animal}{1}{1},3), ...
        mice(curr_animal).interneurons);
    
    curr_frames = analysis.epoch_frames{curr_animal}{2}(1): ...
        analysis.epoch_frames{curr_animal}{2}(2);
    curr_frames_use = curr_frames > 0;
    
    curr_avg_cl = cellfun(@(x) reshape(zscore(permute(nanmean(x{1}(:,curr_frames_use,plot_cells),1),[2 3 1])),[],1), ...
        analysis.epoch_activity_CL{curr_animal},'uni',false);
    
    % get average of last 3 days
    curr_compare_act = nanmean(horzcat(curr_avg_cl{end-2:end}),2);
    
    % get nan values in any day (to ignore)
    nan_vals = any(isnan(horzcat(curr_avg_cl{:})),2);
    
    % set nan values to zero (ok because zscored)
    curr_compare_act(nan_vals) = 0;
    act_nonan = horzcat(curr_avg_cl{:});
    act_nonan(nan_vals,:) = 0;
    
    curr_framerates = arrayfun(@(x) data_all(curr_animal).im(x).framerate, ...
        1:size(act_nonan,2));
    
    [~,curr_xcorr_lagvals] = arrayfun(@(x) xcorr(act_nonan(:,x), ...
        curr_compare_act,round(curr_framerates(x)*2)),1:size(act_nonan,2),'uni',false);
    [~,curr_xcorr_lags] = arrayfun(@(x) max(xcorr(act_nonan(:,x), ...
        curr_compare_act,round(curr_framerates(x)*2))),1:size(act_nonan,2));
    
    xcorr_lags{curr_animal} = arrayfun(@(x) curr_xcorr_lagvals{x}( ...
        curr_xcorr_lags(x))/curr_framerates(x),1:length(curr_xcorr_lags));
    
end


figure; hold on;
for curr_animal = find(pmm_animals)
    plot(xcorr_lags{curr_animal});    
end



%% Slide 4) Get fraction of trials active per cell pre/post first reversal

% get pre-reversal fraction active
cl_frac_active = cellfun(@(x,mice) cellfun(@(x) cellfun(@(x) ...
    permute(nanmean(any(x,2),1),[3 1 2]),x,'uni',false),x,'uni',false), ...
    analysis.epoch_activity_CL,'uni',false);

cr_frac_active = cellfun(@(x,mice) cellfun(@(x) cellfun(@(x) ...
    permute(nanmean(any(x,2),1),[3 1 2]),x,'uni',false),x,'uni',false), ...
    analysis.epoch_activity_CR,'uni',false);

% concatenate the fraction active
cl_frac_active_cat = cellfun(@(x) cell2mat(cellfun(@(x) ...
    horzcat(x{:}),x,'uni',false)),cl_frac_active,'uni',false);

cr_frac_active_cat = cellfun(@(x) cell2mat(cellfun(@(x) ...
    horzcat(x{:}),x,'uni',false)),cr_frac_active,'uni',false);

contingencies_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_contingencies,'uni',false);

sessions_cat = cellfun(@(x) vertcat(x{:}),analysis.epoch_sessions,'uni',false);

% get significantly active cells as more active in 0-2 than -2-0
% baseline_multiepoch_sigdiff = ...
%     cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
%     AP_signrank_matrix( ...
%     permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
%     permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
%     frames(1):frames(2) < frames(2)/2,:),2),[1 3 2])) ...
%     < 0.01 & ...
%     nanmean(permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
%     nanmean(permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
%     frames(1):frames(2) < frames(2)/2,:),2),[1 3 2]),1), ...
%     odor,'uni',false),odor,frames,'uni',false), ...
%     analysis.epoch_activity_odor,analysis.epoch_frames,...
%     'uni',false);

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

% plot the average fraction cells active
cl_frac_active_sessionavg = cell(length(mice),1);
cr_frac_active_sessionavg = cell(length(mice),1);

for curr_animal = 1:length(mice)
    use_cells = setdiff(1:size(cl_frac_active_cat{curr_animal}), ...
        mice(curr_animal).interneurons);
    %plot(sessions_cat{curr_animal},nanmean(cl_frac_active_cat{curr_animal}(use_cells,:)),'k')
    cl_frac_active_sessionavg{curr_animal} = grpstats( ...
        nanmean(cl_frac_active_cat{curr_animal}(use_cells,:)),sessions_cat{curr_animal});
    cr_frac_active_sessionavg{curr_animal} = grpstats( ...
        nanmean(cr_frac_active_cat{curr_animal}(use_cells,:)),sessions_cat{curr_animal});
end

figure;
p1 = subplot(2,1,1); hold on;
xlabel('Session');
ylabel('Fraction active cells per trial');
title('CL')
p2 = subplot(2,1,2); hold on;
xlabel('Session');
ylabel('Fraction active cells per trial');
title('CR')
for curr_cond = 1:2
    for curr_area = 1:2
        switch curr_cond
            case 1
                ls = '-';
                use_cond = rev_animals;
            case 2
                ls = '--';
                use_cond = ctrl_animals;
        end
        
        switch curr_area
            case 1
                col = 'k';
                use_area = pmm_animals;
            case 2
                col = 'r';
                use_area = alm_animals;
        end
               
        use_animals = use_cond & use_area;
        
        session_mean = grpstats(vertcat(cl_frac_active_sessionavg{use_animals}), ...
            horzcat(analysis.data_days{use_animals}));
        session_sem = grpstats(vertcat(cl_frac_active_sessionavg{use_animals}), ...
            horzcat(analysis.data_days{use_animals}),'sem');
        
        %errorbar(session_mean,session_sem,'linestyle',ls,'color',col,'linewidth',2);
        subplot(p1)
        plot(session_mean,'linestyle',ls,'color',col,'linewidth',2)
        
        session_mean = grpstats(vertcat(cr_frac_active_sessionavg{use_animals}), ...
            horzcat(analysis.data_days{use_animals}));
        session_sem = grpstats(vertcat(cr_frac_active_sessionavg{use_animals}), ...
            horzcat(analysis.data_days{use_animals}),'sem');
        
        %errorbar(session_mean,session_sem,'linestyle',ls,'color',col,'linewidth',2);
        subplot(p2)
        plot(session_mean,'linestyle',ls,'color',col,'linewidth',2)
          
    end
end


%% Slide 5) Compare CL and CR
% because visual inspection looked like CL and CR activity became similar
% around time of reversal (which makes sense if they start to lick, or if
% there's some kind of active inhibition or similar signal)

cl_cr_coeff = cell(size(mice));
for curr_animal = 1:length(mice)
    
    plot_cells = setdiff(1:size( ...
        analysis.epoch_activity_CL{curr_animal}{1}{1},3), ...
        mice(curr_animal).interneurons);
    
    curr_frames = analysis.epoch_frames{curr_animal}{2}(1): ...
        analysis.epoch_frames{curr_animal}{2}(2);
    curr_frames_use = curr_frames > 0 & curr_frames < max(curr_frames)/2;
    
    curr_avg_cl = cellfun(@(x) reshape(zscore(permute(nanmean(x{1}(:, ...
        curr_frames_use,plot_cells),1),[2 3 1])),[],1), ...
        analysis.epoch_activity_CL{curr_animal},'uni',false);
    
    curr_avg_cr = cellfun(@(x) reshape(zscore(permute(nanmean(x{1}(:, ...
        curr_frames_use,plot_cells),1),[2 3 1])),[],1), ...
        analysis.epoch_activity_CR{curr_animal},'uni',false);
    
    for curr_session = 1:length(curr_avg_cl)
        % zero nans (ok because zscored)
        curr_cl = curr_avg_cl{curr_session};
        curr_cl(isnan(curr_cl)) = 0;
        
        curr_cr = curr_avg_cr{curr_session};
        curr_cr(isnan(curr_cr)) = 0;
        
        curr_coeff = corrcoef(curr_cl,curr_cr);
        
        cl_cr_coeff{curr_animal}(curr_session) = curr_coeff(2);
    end
    
end

max_sessions = max(cellfun(@length,cl_cr_coeff));
cl_cr_coeff_pad = cellfun(@(x) padarray(x,[0 max_sessions-length(x)],nan,'post'), ...
    cl_cr_coeff,'uni',false);

cl_cr_coeff_cat = vertcat(cl_cr_coeff_pad{:});

% get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

cl_cr_coeff_cat_minsub = bsxfun(@minus,cl_cr_coeff_cat,min(cl_cr_coeff_cat,[],2));
cl_cr_coeff_cat_norm = bsxfun(@times,cl_cr_coeff_cat_minsub,1./max(cl_cr_coeff_cat_minsub,[],2));

figure; hold on;
use_animals = pmm_animals & rev_animals;
errorbar(nanmean(cl_cr_coeff_cat_norm(use_animals,:)), ...
    nanstd(cl_cr_coeff_cat_norm(use_animals,:))./ ...
    sqrt(sum(~isnan(cl_cr_coeff_cat_norm(use_animals,:)))),'k','linewidth',2,'linestyle','-');

use_animals = pmm_animals & ctrl_animals;
errorbar(nanmean(cl_cr_coeff_cat_norm(use_animals,:)), ...
    nanstd(cl_cr_coeff_cat_norm(use_animals,:))./ ...
    sqrt(sum(~isnan(cl_cr_coeff_cat_norm(use_animals,:)))),'k','linewidth',2,'linestyle','--');

use_animals = alm_animals & rev_animals;
errorbar(nanmean(cl_cr_coeff_cat_norm(use_animals,:)), ...
    nanstd(cl_cr_coeff_cat_norm(use_animals,:))./ ...
    sqrt(sum(~isnan(cl_cr_coeff_cat_norm(use_animals,:)))),'r','linewidth',2,'linestyle','-');

use_animals = alm_animals & ctrl_animals;
errorbar(nanmean(cl_cr_coeff_cat_norm(use_animals,:)), ...
    nanstd(cl_cr_coeff_cat_norm(use_animals,:))./ ...
    sqrt(sum(~isnan(cl_cr_coeff_cat_norm(use_animals,:)))),'r','linewidth',2,'linestyle','--');

ylabel('CL-CR correlation')
xlabel('Session')
legend({'rev pmm' 'ctrl pmm' 'rev alm' 'ctrl alm'})

% group pre/post, plot as bar
cl_cr_coeff_pre = nanmean(cl_cr_coeff_cat_norm(:,3:5),2);
cl_cr_coeff_post = nanmean(cl_cr_coeff_cat_norm(:,7:9),2);

animal_grp = nan(size(mice));
animal_grp(pmm_animals & rev_animals) = 1;
animal_grp(pmm_animals & ctrl_animals) = 2;
animal_grp(alm_animals & rev_animals) = 3;
animal_grp(alm_animals & ctrl_animals) = 4;

plot_pre_mean = grpstats(cl_cr_coeff_pre,animal_grp);
plot_pre_sem = grpstats(cl_cr_coeff_pre,animal_grp,'sem');

plot_post_mean = grpstats(cl_cr_coeff_post,animal_grp);
plot_post_sem = grpstats(cl_cr_coeff_post,animal_grp,'sem');

figure;
bar([plot_pre_mean plot_post_mean],'grouped')
ylabel('CL-CR correlation')
set(gca,'XTickLabel',{'rev pmm' 'ctrl pmm' 'rev alm' 'ctrl alm'});
legend({'Sessions 3-5','Sessions 7-9'});

% get significantly active cells 
% (in baseline seconds -2-0 vs odor period 0-2)
cl_baseline_sigdiff = ...
    cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
    AP_signrank_matrix( ...
    permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2])) ...
    < 0.05 & ...
    nanmean(permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
    nanmean(permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2]),1), ...
    odor,'uni',false),odor,frames,'uni',false), ...
    analysis.epoch_activity_CL,analysis.epoch_frames,...
    'uni',false);

cr_baseline_sigdiff = ...
    cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
    AP_signrank_matrix( ...
    permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2])) ...
    < 0.05 & ...
    nanmean(permute(nanmean(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
    nanmean(permute(nanmean(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2]),1), ...
    odor,'uni',false),odor,frames,'uni',false), ...
    analysis.epoch_activity_CR,analysis.epoch_frames,...
    'uni',false);

cl_baseline_sigdiff_prepost = cellfun(@(x) ...
    [nanmean(cell2mat(cellfun(@(x) x{1},x(3:5),'uni',false)')); ...
    nanmean(cell2mat(cellfun(@(x) x{1},x(7:9),'uni',false)'))], ...
    cl_baseline_sigdiff,'uni',false);

cr_baseline_sigdiff_prepost = cellfun(@(x) ...
    [nanmean(cell2mat(cellfun(@(x) x{1},x(3:5),'uni',false)')); ...
    nanmean(cell2mat(cellfun(@(x) x{1},x(7:9),'uni',false)'))], ...
    cr_baseline_sigdiff,'uni',false);


% Find significant cells by fraction active, get cells which are active by
% baseline but NOT by condition
baseline_sigdiff = ...
    cellfun(@(odor,frames) cellfun(@(odor,frames) cellfun(@(odor) ...
    AP_chisquare_matrix( ...
    permute(any(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]), ...
    permute(any(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2])) ...
    < 0.05 & ...
    nanmean(permute(any(odor(:,frames(1):frames(2) < 0,:),2),[1 3 2]),1) < ...
    nanmean(permute(any(odor(:,frames(1):frames(2) >= 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2]),1), ...
    odor,'uni',false),odor,frames,'uni',false), ...
    analysis.epoch_activity_odor,analysis.epoch_frames,...
    'uni',false);

condition_sigdiff = ...
    cellfun(@(CL,CR,frames) cellfun(@(CL,CR,frames) cellfun(@(CL,CR) ...
    AP_chisquare_matrix( ...
    permute(any(CL(:,frames(1):frames(2) > 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2]), ...
    permute(any(CR(:,frames(1):frames(2) > 0 & ...
    frames(1):frames(2) < frames(2)/2,:),2),[1 3 2])) ...
    < 0.05, ...
    CL,CR,'uni',false),CL,CR,frames,'uni',false), ...
    analysis.epoch_activity_CL,analysis.epoch_activity_CR,analysis.epoch_frames,...
    'uni',false);







%% Slide 6) Activity aligned/sorted by lick onset


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

plot_corr = true;

max_sessions = max(cellfun(@length,analysis.odor_trials));

odor_lick_align_slope_ratio = nan(length(mice),max_sessions);

for curr_animal = [3 6];%1:length(mice)
    
    disp(curr_animal);
    
    for curr_session = [5 9];%1:length(data_all(curr_animal).im)
        
        % define reversal
        curr_rev = 1;
        revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
            size(analysis.odor_trials{curr_animal}{curr_session},1)];
        % if there are < 10 trials in reversal, assume early/not real
        % reversal and use the next one
        rev_trials = diff(revs);
        if rev_trials(curr_rev) < 10
            curr_rev = curr_rev + 1;
        end
        
        % only use first 2 seconds (odor delivery) of activity
        odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
        
        curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (:,curr_odor_frames,:);
        % skip this if less than 20 trials
        if size(curr_act,1) < 20
            continue
        end
        curr_act_zscore = reshape(permute(zscore(curr_act,[],2),[2 3 1]),[],size(curr_act,1));
        
        % lick isn't segregated by reversal, do that here (after running the load
        % again because odor trials weren't previously stored)      
        curr_trials = false(size(analysis.odor_trials{curr_animal}{curr_session},1),1);
        curr_trials(revs(curr_rev):revs(curr_rev+1)-1) = true;
        
        curr_lick = analysis.lick_rate{curr_animal}{curr_session}( ...
            curr_trials & analysis.condition_trials{curr_animal}{curr_session}(:,1),:);
        
        % Get the first time where licking is over 3 hz, sort
        [~,lick_start] = max(curr_lick > 3,[],2);
        [~,lick_start_sort] = sort(lick_start);
        
        % Get rid of weird outlier lick trials (for now just cutoff)
        early_cutoff = 50; % depends on lick rate resolution, now 1/100 sec (*10 = ms)
        early_lick_trials = lick_start < early_cutoff;
        
        % Plot activity correlation v lick time start difference, excluding early
        
        lick_start_diff = pdist2(lick_start(~early_lick_trials),lick_start(~early_lick_trials));
        lick_corr = corrcoef(curr_lick(~early_lick_trials,:)'>3);
        
        act_corr = corrcoef(curr_act_zscore(:,~early_lick_trials),'rows','pairwise');
        
        lick_start_diff_tril = AP_itril(lick_start_diff,-1);
        act_corr_tril = AP_itril(act_corr,-1);
        
        lick_bins = [0:10:100];
        lick_bins_plot = lick_bins(1:end-1) + diff(lick_bins)/2;
        lick_bins(end) = inf;
        
        [lick_bins_histc lick_use_bins] = histc(lick_start_diff_tril,lick_bins);
        act_corr_mean = grpstats(act_corr_tril,lick_use_bins);
        act_corr_sem = grpstats(act_corr_tril,lick_use_bins,'sem');
        
        % aligning activity to lick onset (v odor onset)
        
        % how to align pre/post (currently 1.5 pre / 0.5 post, which is similar to
        % odor onset aligned)
        post_frames = round(data_all(curr_animal).im(curr_session).framerate*0.5);
        pre_frames = round(data_all(curr_animal).im(curr_session).framerate*1.5);
        
        % convert lick times into frames
        lick_start_frames = round(lick_start/100* ...
            data_all(curr_animal).im(curr_session).framerate);
        
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        curr_start_frame = find(curr_frames > 0,1);
        
        act_use_frames = arrayfun(@(x) (curr_start_frame + lick_start_frames(x) - pre_frames): ...
            (curr_start_frame + lick_start_frames(x) + post_frames),1:length(lick_start_frames),'uni',false);
        
        curr_act_lickalign = cell2mat( ...
            arrayfun(@(x) analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (x,act_use_frames{x},:),find(~early_lick_trials),'uni',false));
        
        curr_act_lickalign_zscore = reshape(permute(zscore(curr_act_lickalign,[],2), ...
            [2 3 1]),[],size(curr_act_lickalign,1));
        
        act_lickalign_corr = corrcoef(curr_act_lickalign_zscore,'rows','pairwise');
        act_lickalign_corr_tril = AP_itril(act_lickalign_corr,-1);
        
        act_lickalign_corr_mean = grpstats(act_lickalign_corr_tril,lick_use_bins);
        act_lickalign_corr_sem = grpstats(act_lickalign_corr_tril,lick_use_bins,'sem');
        
        
        % align randomly for control
        num_shuff = 1;
        act_shuffle_lickalign_corr_mean = nan(length(unique(lick_use_bins)),num_shuff);
        act_shuffle_lickalign_corr_sem = nan(length(unique(lick_use_bins)),num_shuff);
        for curr_shuff = 1:num_shuff;
            shuffle_act_use_frames(~early_lick_trials) = shake(act_use_frames(~early_lick_trials));
            
            curr_shuffle_act_lickalign = cell2mat( ...
                arrayfun(@(x) analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (x,shuffle_act_use_frames{x},:),find(~early_lick_trials),'uni',false));
            
            curr_act_shuffle_lickalign_zscore = reshape(permute(zscore(curr_shuffle_act_lickalign,[],2), ...
                [2 3 1]),[],size(curr_shuffle_act_lickalign,1));
            
            act_shuffle_lickalign_corr = corrcoef(curr_act_shuffle_lickalign_zscore,'rows','pairwise');
            act_shuffle_lickalign_corr_tril = AP_itril(act_shuffle_lickalign_corr,-1);
            
            act_shuffle_lickalign_corr_mean(:,curr_shuff) = grpstats(act_shuffle_lickalign_corr_tril,lick_use_bins);
            act_shuffle_lickalign_corr_sem(:,curr_shuff) = grpstats(act_shuffle_lickalign_corr_tril,lick_use_bins,'sem');
            
            %disp(curr_shuff);
        end
        
        
        if plot_corr
            figure; hold on;
            errorbar(lick_bins_plot(unique(lick_use_bins)), ...
                act_corr_mean,act_corr_sem,'k','linewidth',2);
            errorbar(lick_bins_plot(unique(lick_use_bins)), ...
                act_lickalign_corr_mean,act_lickalign_corr_sem,'r','linewidth',2);
            errorbar(lick_bins_plot(unique(lick_use_bins)), ...
                nanmean(act_shuffle_lickalign_corr_mean,2), ...
                nanmean(act_shuffle_lickalign_corr_sem,2),'b','linewidth',2);
            
            legend({'Odor-aligned','Lick-aligned','Random lick-aligned'});
            xlabel('Difference in lick onset time (10 ms)')
            ylabel('Population correlation');
            
        end
        
        % fit the lines and get the slopes
        lick_bins_filled = lick_bins_plot(unique(lick_use_bins))';
        use_bins = 1:min(length(lick_bins_filled),4);
        
        fit_act = polyfit(lick_bins_filled(use_bins),act_corr_mean(use_bins),1);
        
        fit_act_lickalign = polyfit(lick_bins_filled(use_bins), ...
            act_lickalign_corr_mean(use_bins),1);
        
        fit_act_shuffle_lickalign = polyfit(lick_bins_filled(use_bins), ...
            nanmean(act_shuffle_lickalign_corr_mean(use_bins,:),2),1);
        
        % get the ratio of slopes between odor and lick aligned
        odor_lick_align_slope_ratio(curr_animal,curr_session) = fit_act(1)/fit_act_lickalign(1);
        
    end
end

figure; hold on;
plot(nanmedian(odor_lick_align_slope_ratio(pmm_animals & rev_animals,:)),'k')
plot(nanmedian(odor_lick_align_slope_ratio(alm_animals & rev_animals,:)),'r')
plot(nanmedian(odor_lick_align_slope_ratio(pmm_animals & ctrl_animals,:)),'--k')
plot(nanmedian(odor_lick_align_slope_ratio(alm_animals & ctrl_animals,:)),'--r')
line(xlim,[1 1])
legend({'pmm rev','alm rev','pmm ctrl','alm ctrl'})
xlabel('Session')
ylabel('Odor-aligned corr slope / lick-aligned corr slope')

% % try backwards - align a cell,and then align lick? or see how other cells
% % align?
% a = curr_act(~early_lick_trials,:,179);
% [c lags] = arrayfun(@(x) xcorr(a(x,:),nanmean(a),10),1:size(a,1),'uni',false);
% [~,max_c_idx] = cellfun(@(x) max(x),c,'uni',false);
% max_lags = cellfun(@(x,y) x(y),lags,max_c_idx);
% a_align = cell2mat(arrayfun(@(x) circshift(a(x,:),[0 -max_lags(x)]),[1:size(a,1)]','uni',false));
% 
% b = curr_act(~early_lick_trials,:,89);
% b_align = cell2mat(arrayfun(@(x) circshift(b(x,:),[0 -max_lags(x)]),[1:size(b,1)]','uni',false));
% 
% figure;imagesc(b);colormap(gray);
% figure;imagesc(b_align);colormap(gray);
% 
% lag_align_licktime = round((max_lags/data_all(curr_animal).im(curr_session).framerate)*100);
% use_lick = curr_lick(~early_lick_trials,:);
% lick_align = cell2mat(arrayfun(@(x) circshift(use_lick(x,:), ...
%     [0 -lag_align_licktime(x)]),[1:size(use_lick,1)]','uni',false));


% TO DO: somehow see if there are odor components aligned to odor and lick
% components aligned to lick?

% TO DO: look at the difference in alignment, seems clear time component,
% but also clear secondary non-time component



%% Find lick/odor-aligned cells (somehow...)

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

max_sessions = max(cellfun(@length,analysis.odor_trials));

for curr_animal = [10 16];%1:length(mice)
    
    disp(curr_animal);
    
    for curr_session = 1:length(data_all(curr_animal).im)
        
        % define reversal
        curr_rev = 1;
        revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
            size(analysis.odor_trials{curr_animal}{curr_session},1)];
        % if there are < 10 trials in reversal, assume early/not real
        % reversal and use the next one
        rev_trials = diff(revs);
        if rev_trials(curr_rev) < 10
            curr_rev = curr_rev + 1;
        end
        
        % only use first 2 seconds (odor delivery) of activity
        odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        %%% TEMP CHANGE
        curr_odor_frames = curr_frames > 0;%& curr_frames <= odor_frames;
        
        curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (:,curr_odor_frames,:);
        % skip this if less than 20 trials
        if size(curr_act,1) < 20
            continue
        end
        curr_act_zscore = reshape(permute(zscore(curr_act,[],2),[2 3 1]),[],size(curr_act,1));
        
        % lick isn't segregated by reversal, do that here (after running the load
        % again because odor trials weren't previously stored)      
        curr_trials = false(size(analysis.odor_trials{curr_animal}{curr_session},1),1);
        curr_trials(revs(curr_rev):revs(curr_rev+1)-1) = true;
        
        curr_lick = analysis.lick_rate{curr_animal}{curr_session}( ...
            curr_trials & analysis.condition_trials{curr_animal}{curr_session}(:,1),:);
        
        % Get the first time where licking is over 3 hz, sort
        [~,lick_start] = max(curr_lick > 3,[],2);
        [~,lick_start_sort] = sort(lick_start);
        
        % Get rid of weird outlier lick trials (for now just cutoff)
        early_cutoff = 50; % depends on lick rate resolution, now 1/100 sec (*10 = ms)
        early_lick_trials = lick_start < early_cutoff;
        
        
        % aligning activity to lick onset (v odor onset)
        
        % how to align pre/post (currently 1.5 pre / 0.5 post, which is similar to
        % odor onset aligned)
        post_frames = round(data_all(curr_animal).im(curr_session).framerate*5);
        pre_frames = round(data_all(curr_animal).im(curr_session).framerate*1);
        
        % convert lick times into frames
        lick_start_frames = round(lick_start/100* ...
            data_all(curr_animal).im(curr_session).framerate);
        
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        curr_start_frame = find(curr_frames > 0,1);
        
        act_use_frames = arrayfun(@(x) (curr_start_frame + lick_start_frames(x) - pre_frames): ...
            (curr_start_frame + lick_start_frames(x) + post_frames),1:length(lick_start_frames),'uni',false);
        
        curr_act_lickalign = cell2mat( ...
            arrayfun(@(x) analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (x,act_use_frames{x},:),find(~early_lick_trials),'uni',false));
        
        curr_act_lickalign_zscore = reshape(permute(zscore(curr_act_lickalign,[],2), ...
            [2 3 1]),[],size(curr_act_lickalign,1));
        
    end
end




%% TEST (paired with above cell)

% Plot aligned activity odor / lick over days for cell
curr_cell = 54;
curr_animal = 10;
for curr_session = 1:length(data_all(curr_animal).im)
    
    % define reversal
    curr_rev = 1;
    revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
        size(analysis.odor_trials{curr_animal}{curr_session},1)];
    % if there are < 10 trials in reversal, assume early/not real
    % reversal and use the next one
    rev_trials = diff(revs);
    if rev_trials(curr_rev) < 10
        curr_rev = curr_rev + 1;
    end
    
    % only use first 2 seconds (odor delivery) of activity
    odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
    curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
        analysis.epoch_frames{curr_animal}{curr_session}(2);
    curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
    
    curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
        (:,curr_odor_frames,:);
    % skip this if less than 20 trials
    if size(curr_act,1) < 20
        continue
    end
    curr_act_zscore = reshape(permute(zscore(curr_act,[],2),[2 3 1]),[],size(curr_act,1));
    
    % lick isn't segregated by reversal, do that here (after running the load
    % again because odor trials weren't previously stored)
    curr_trials = false(size(analysis.odor_trials{curr_animal}{curr_session},1),1);
    curr_trials(revs(curr_rev):revs(curr_rev+1)-1) = true;
    
    curr_lick = analysis.lick_rate{curr_animal}{curr_session}( ...
        curr_trials & analysis.condition_trials{curr_animal}{curr_session}(:,1),:);
    
    % Get the first time where licking is over 3 hz, sort
    [~,lick_start] = max(curr_lick > 3,[],2);
    [~,lick_start_sort] = sort(lick_start);
    
    % Get rid of weird outlier lick trials (for now just cutoff)
    early_cutoff = 50; % depends on lick rate resolution, now 1/100 sec (*10 = ms)
    early_lick_trials = lick_start < early_cutoff;
    
    
    % aligning activity to lick onset (v odor onset)
    
    % time for surrounding align pre/post 
    post_frames = round(data_all(curr_animal).im(curr_session).framerate*1.5);
    pre_frames = round(data_all(curr_animal).im(curr_session).framerate*0.5);
    
    % convert lick times into frames
    lick_start_frames = round(lick_start/100* ...
        data_all(curr_animal).im(curr_session).framerate);
    
    curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
        analysis.epoch_frames{curr_animal}{curr_session}(2);
    curr_start_frame = find(curr_frames > 0,1);
    
    act_use_frames = arrayfun(@(x) (curr_start_frame + lick_start_frames(x) - pre_frames): ...
        (curr_start_frame + lick_start_frames(x) + post_frames),1:length(lick_start_frames),'uni',false);
    
    curr_act_lickalign = cell2mat( ...
        arrayfun(@(x) analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
        (x,act_use_frames{x},:),find(~early_lick_trials),'uni',false));
    
    curr_act_lickalign_zscore = reshape(permute(zscore(curr_act_lickalign,[],2), ...
        [2 3 1]),[],size(curr_act_lickalign,1));
    
    figure; hold on;
    plot(nanmean(curr_act_zscore(:,:,curr_cell)),'k');
    plot(nanmean(curr_act_lickalign_zscore(:,:,curr_cell)),'k');
end



%% Slide 7) Plot the pca space thing for looking at spread/distance of rev activity


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

max_sessions = max(cellfun(@length,analysis.odor_trials));

all_active_trials = cell(length(animals),max_sessions);
all_contingencies = nan(length(animals),max_sessions);
for curr_animal = 1:length(mice)
    
    disp(curr_animal);
    
    for curr_session = 1:length(data_all(curr_animal).im);
        
        % define reversal (use only the first so that it's settled in)
        curr_rev = 1;
        revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
            size(analysis.odor_trials{curr_animal}{curr_session},1)];
        % if there are < 10 trials in reversal, assume early/not real
        % reversal and use the next one
        rev_trials = diff(revs);
        if rev_trials(curr_rev) < 10
            curr_rev = curr_rev + 1;
        end
        
        % only use first 2 seconds (odor delivery) of activity
        odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
        
        curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (:,curr_odor_frames,:);
        % skip this if less than 20 trials
        if size(curr_act,1) < 20
            continue
        end
        
        curr_active_trials = permute(nanmean(any(curr_act,2),1),[3 1 2]);
        
        all_active_trials{curr_animal,curr_session} = curr_active_trials;
        all_contingencies(curr_animal,curr_session) =  ...
            analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);
    end
end


% Distance of the first n PCs
for curr_animal = find(alm_animals & rev_animals);%1:length(mice)
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    
    % euclidean distance between PC1+2 
    pcs = 2;
    pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
    
    % Get the sessions that are contingency 1 and 2
    used_sessions = find(~isnan(all_contingencies(curr_animal,:)));
    cont1_sessions = find(all_contingencies(curr_animal,:) == 1);
    cont2_sessions = find(all_contingencies(curr_animal,:) == 2);
    
    cont1_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 1;
    cont2_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 2;
    
    pc_dist_xcont = pc_dist(cont1_session_idx,cont2_session_idx);    
    
    figure;imagesc(pc_dist_xcont);colormap(gray);

end


% Distance of the first n PCs
alm_rev_cont_dist = cell(3,1);
alm_rev_animals = find(alm_animals & rev_animals);
for curr_animal_idx = 1:length(alm_rev_animals);
    
    curr_animal = alm_rev_animals(curr_animal_idx);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    
    % euclidean distance between PC1+2 
    pcs = 3;
    pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
    
    % Get the sessions that are contingency 1 and 2
    used_sessions = find(~isnan(all_contingencies(curr_animal,:)));
    cont1_sessions = find(all_contingencies(curr_animal,:) == 1);
    cont2_sessions = find(all_contingencies(curr_animal,:) == 2);
    
    cont1_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 1;
    cont2_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 2;
    
    pc_dist_xcont = pc_dist(cont1_session_idx,cont2_session_idx);    
    
    alm_rev_cont_dist{curr_animal_idx} = pc_dist_xcont;
end
max_1 = max(cellfun(@(x) size(x,1),alm_rev_cont_dist));
max_2 = max(cellfun(@(x) size(x,2),alm_rev_cont_dist));
alm_rev_cont_dist_norm = cellfun(@(x) (x-min(x(:)))/max(x(:)),alm_rev_cont_dist,'uni',false);
alm_rev_cont_dist_pad = cellfun(@(x) padarray(x, ...
    [max_1-size(x,1) max_2-size(x,2)],nan,'post'),alm_rev_cont_dist_norm,'uni',false);
alm_rev_cont_dist_cat = cat(3,alm_rev_cont_dist_pad{:});
figure;imagesc(nanmean(alm_rev_cont_dist_cat,3));colormap(gray);

% Distance of the first n PCs
pmm_rev_cont_dist = cell(3,1);
pmm_rev_animals = find(pmm_animals & rev_animals);
for curr_animal_idx = 1:length(pmm_rev_animals);
    
    curr_animal = pmm_rev_animals(curr_animal_idx);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    
    % euclidean distance between PC1+2 
    pcs = 3;
    pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
    
    % Get the sessions that are contingency 1 and 2
    used_sessions = find(~isnan(all_contingencies(curr_animal,:)));
    cont1_sessions = find(all_contingencies(curr_animal,:) == 1);
    cont2_sessions = find(all_contingencies(curr_animal,:) == 2);
    
    cont1_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 1;
    cont2_session_idx = all_contingencies(curr_animal, ...
        ~isnan(all_contingencies(curr_animal,:))) == 2;
    
    pc_dist_xcont = pc_dist(cont1_session_idx,cont2_session_idx);    
    
    pmm_rev_cont_dist{curr_animal_idx} = pc_dist_xcont;
end
max_1 = max(cellfun(@(x) size(x,1),pmm_rev_cont_dist));
max_2 = max(cellfun(@(x) size(x,2),pmm_rev_cont_dist));
pmm_rev_cont_dist_norm = cellfun(@(x) (x-min(x(:)))/max(x(:)),pmm_rev_cont_dist,'uni',false);
pmm_rev_cont_dist_pad = cellfun(@(x) padarray(x, ...
    [max_1-size(x,1) max_2-size(x,2)],nan,'post'),pmm_rev_cont_dist_norm,'uni',false);
pmm_rev_cont_dist_cat = cat(3,pmm_rev_cont_dist_pad{:});
figure;imagesc(nanmean(pmm_rev_cont_dist_cat,3));colormap(gray);


%% Slide 7) PCA contingency epoch activity distance


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);

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
            
            % only use first 2 seconds (odor delivery) of activity
            odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (:,curr_odor_frames,:);
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(nanmean(any(curr_act,2),1),[3 1 2]);
            
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


% Get reasonable number of PCs
pc_var = cell(length(mice),1);
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    pc_var{curr_animal} = cumsum(latent)/sum(latent);
end

pcs = 5;


% Distance of the first n PCs WITHIN EARLY/LATE
activity_distance = cell(length(mice),1);
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    
    % euclidean distance between selected PCs 
    pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));  
    %pc_dist_norm = pc_dist./max(pc_dist(:));
    pc_dist_norm = (pc_dist-nanmean(AP_itril(pc_dist,-1)))./nanstd(AP_itril(pc_dist,-1));
    
    use_dist = pc_dist_norm;
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    
    % split/save distances, dim1 = animal, dim2 = early/late, 
    % dim3 = 1v1 / 2v2 / 1v2
    
    activity_distance{curr_animal,1,1} = ...
        AP_itril(use_dist(~all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),-1);
    
    activity_distance{curr_animal,2,1} = ...
        AP_itril(use_dist(all_latereversal{curr_animal,:} & cont1_session, ...
        all_latereversal{curr_animal,:} & cont1_session),-1);
    
    activity_distance{curr_animal,1,2} = ...
        AP_itril(use_dist(~all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),-1);
    
    activity_distance{curr_animal,2,2} = ...
        AP_itril(use_dist(all_latereversal{curr_animal,:} & cont2_session, ...
        all_latereversal{curr_animal,:} & cont2_session),-1);
    
    activity_distance{curr_animal,1,3} = ...
        reshape(use_dist(~all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,2,3} = ...
        reshape(use_dist(all_latereversal{curr_animal,:} & cont1_session, ...
        all_latereversal{curr_animal,:} & cont2_session),[],1);
    
end

activity_distance_mean_withinrev = cellfun(@nanmean,activity_distance);


% Distance of the first n PCs ACROSS EARLY/LATE
activity_distance = cell(length(mice),1);
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');
    
    % euclidean distance between selected PCs 
    pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));   
    %pc_dist_norm = pc_dist./max(pc_dist(:));
    pc_dist_norm = (pc_dist-nanmean(AP_itril(pc_dist,-1)))./nanstd(AP_itril(pc_dist,-1));
    
    use_dist = pc_dist_norm;
    
    cont1_session = all_contingencies{curr_animal} == 1;
    cont2_session = all_contingencies{curr_animal} == 2;
    
    
    % split/save distances, dim1 = animal, 
    % dim2 = 1v1 / 2v2 / 1v2 / 2v1 (former = early, latter = late)
      
    activity_distance{curr_animal,1} = ...
        reshape(use_dist(all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),[],1);
    
    activity_distance{curr_animal,2} = ...
        reshape(use_dist(all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,3} = ...
        reshape(use_dist(all_latereversal{curr_animal,:} & cont1_session, ...
        ~all_latereversal{curr_animal,:} & cont2_session),[],1);
    
    activity_distance{curr_animal,4} = ...
        reshape(use_dist(all_latereversal{curr_animal,:} & cont2_session, ...
        ~all_latereversal{curr_animal,:} & cont1_session),[],1);
    
end

activity_distance_mean_acrossrev = cellfun(@nanmean,activity_distance);

% Plot by contingency and area 
alm_distance_acrossrev = activity_distance_mean_acrossrev(alm_animals & rev_animals,:,:);
alm_distance_acrossrev_mean = nanmean(alm_distance_acrossrev);
alm_distance_acrossrev_sem = nanstd(alm_distance_acrossrev)./sqrt(sum(~isnan(alm_distance_acrossrev)));

alm_distance_withinrev = permute(activity_distance_mean_withinrev(alm_animals & rev_animals,:,:),[3 2 1]);
alm_distance_withinrev_mean = nanmean(alm_distance_withinrev,3);
alm_distance_withinrev_sem = nanstd(alm_distance_withinrev,[],3)./sqrt(sum(~isnan(alm_distance_withinrev),3));

alm_distance_plot_mean = [[alm_distance_withinrev_mean alm_distance_acrossrev_mean(1:3)'] ...
    [nan(2,1);alm_distance_acrossrev_mean(4)]];
alm_distance_plot_sem = [[alm_distance_withinrev_sem alm_distance_acrossrev_sem(1:3)'] ...
    [nan(2,1);alm_distance_acrossrev_sem(4)]];


pmm_distance_withinrev = permute(activity_distance_mean_withinrev(pmm_animals & rev_animals,:,:),[3 2 1]);
pmm_distance_withinrev_mean = nanmean(pmm_distance_withinrev,3);
pmm_distance_withinrev_sem = nanstd(pmm_distance_withinrev,[],3)./sqrt(sum(~isnan(pmm_distance_withinrev),3));

pmm_distance_acrossrev = permute(activity_distance_mean_acrossrev(pmm_animals & rev_animals,:,:),[3 2 1]);
pmm_distance_acrossrev_mean = nanmean(pmm_distance_acrossrev,3);
pmm_distance_acrossrev_sem = nanstd(pmm_distance_acrossrev,[],3)./sqrt(sum(~isnan(pmm_distance_acrossrev),3));

pmm_distance_plot_mean = [[pmm_distance_withinrev_mean pmm_distance_acrossrev_mean(1:3)'] ...
    [nan(2,1);pmm_distance_acrossrev_mean(4)]];
pmm_distance_plot_sem = [[pmm_distance_withinrev_sem pmm_distance_acrossrev_sem(1:3)'] ...
    [nan(2,1);pmm_distance_acrossrev_sem(4)]];

figure; 
p1 = subplot(2,5,1:4);
errorb(alm_distance_plot_mean,alm_distance_plot_sem);
colormap(gray); 
ylabel('Distance')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early/late' 'Across early/late alt'},'location','nw');
title('ALM')

p2 = subplot(2,5,6:9);
errorb(pmm_distance_plot_mean,pmm_distance_plot_sem);
colormap(gray); 
ylabel('Distance')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early/late' 'Across early/late alt'},'location','nw');
title('PMM')

% Plot control early v late contingency 1
alm_distance_withinrev_ctrl = activity_distance_mean_withinrev(alm_animals & ctrl_animals,:,1);
alm_distance_withinrev_ctrl_mean = nanmean(alm_distance_withinrev_ctrl);
alm_distance_withinrev_ctrl_sem = nanstd(alm_distance_withinrev_ctrl)./ ...
    sqrt(sum(~isnan(alm_distance_withinrev_ctrl)));

alm_distance_acrossrev_ctrl = activity_distance_mean_acrossrev(alm_animals & ctrl_animals,1);
alm_distance_acrossrev_ctrl_mean = nanmean(alm_distance_acrossrev_ctrl);
alm_distance_acrossrev_ctrl_sem = nanstd(alm_distance_acrossrev_ctrl)./ ...
    sqrt(sum(~isnan(alm_distance_acrossrev_ctrl)));

alm_distance_ctrl_plot_mean = [alm_distance_withinrev_ctrl_mean alm_distance_acrossrev_ctrl_mean];
alm_distance_ctrl_plot_sem = [alm_distance_withinrev_ctrl_sem alm_distance_acrossrev_ctrl_sem];


pmm_distance_withinrev_ctrl = activity_distance_mean_withinrev(pmm_animals & ctrl_animals,:,1);
pmm_distance_withinrev_ctrl_mean = nanmean(pmm_distance_withinrev_ctrl);
pmm_distance_withinrev_ctrl_sem = nanstd(pmm_distance_withinrev_ctrl)./ ...
    sqrt(sum(~isnan(pmm_distance_withinrev_ctrl)));

pmm_distance_acrossrev_ctrl = activity_distance_mean_acrossrev(pmm_animals & ctrl_animals,1);
pmm_distance_acrossrev_ctrl_mean = nanmean(pmm_distance_acrossrev_ctrl);
pmm_distance_acrossrev_ctrl_sem = nanstd(pmm_distance_acrossrev_ctrl)./ ...
    sqrt(sum(~isnan(pmm_distance_acrossrev_ctrl)));

pmm_distance_ctrl_plot_mean = [pmm_distance_withinrev_ctrl_mean pmm_distance_acrossrev_ctrl_mean];
pmm_distance_ctrl_plot_sem = [pmm_distance_withinrev_ctrl_sem pmm_distance_acrossrev_ctrl_sem];

subplot(2,5,5); hold on;
bar(alm_distance_ctrl_plot_mean);
errorb(alm_distance_ctrl_plot_mean,alm_distance_ctrl_plot_sem);
colormap(gray); 
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early' 'Late' 'Across'});
ylim(get(p1,'ylim'));
title('ALM control');

subplot(2,5,10); hold on;
bar(pmm_distance_ctrl_plot_mean);
errorb(pmm_distance_ctrl_plot_mean,pmm_distance_ctrl_plot_sem);
colormap(gray); 
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early' 'Late' 'Across'});
ylim(get(p2,'ylim'));
title('PMM control');


% Try plotting this as a box plot: might make more sense because can
% normalize and not get confused by box height which is arbitrary

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
ylabel('Normalized distance')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early/late' 'Across early/late alt'},'location','nw');
title('ALM')

p2 = subplot(2,5,6:9);
aboxplot(pmm_distance_boxplot,'colorgrad','green_down')
ylabel('Normalized distance')
set(gca,'XTickLabel',{'1 v 1' '2 v 2' '1 v 2'})
legend({'Within early' 'Within late' 'Across early/late' 'Across early/late alt'},'location','nw');
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



%% Slide 8) PCA distance of trial activity

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);

for curr_animal = 1:length(mice)
        
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
            
            % only use first 2 seconds (odor delivery) of activity
            odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (:,curr_odor_frames,:);
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(any(curr_act,2),[1 3 2]);
            
            curr_contingency = analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);           
            % if control animal, don't use anything past session 10 or if
            % contingency 2 has appeared (because reversal in these
            % animals, but it's late - I'm not using this data for now)
            if ismember(curr_animal,find(ctrl_animals));
                if curr_session > 10 || curr_contingency == 2
                    continue
                end
            end
            
            all_active_trials{curr_animal} = [all_active_trials{curr_animal};curr_active_trials];
            all_contingencies{curr_animal} = [all_contingencies{curr_animal}; ...
                repmat(curr_contingency, size(curr_active_trials,1),1)];
            
            % 'Late' for reversal animals = anything including and after
            % the second reversal
            if ismember(curr_animal,find(rev_animals));
                all_latereversal{curr_animal} = [all_latereversal{curr_animal}; ...
                    repmat(sum(diff(all_contingencies{curr_animal}) ~= 0) > 1, ...
                    size(curr_active_trials,1),1)];
            else
                % 'Late for control animals = anything after the 6th session
                all_latereversal{curr_animal} = [all_latereversal{curr_animal}; ...
                    repmat(curr_session >= 6, ...
                    size(curr_active_trials,1),1)];
            end
            
            
                       
        end
    end
end


% Get reasonable number of PCs
pc_var = cell(length(mice),1);
for curr_animal = 1:length(mice);
        
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff,score,latent] = princomp(curr_act);
    pc_var{curr_animal} = cumsum(latent)/sum(latent);
end

pcs = 10;


% Plot example animal
curr_animal = 12;
curr_act = horzcat(all_active_trials{curr_animal,:});
[coeff,score,latent] = princomp(curr_act);
figure; 

p1 = subplot(2,1,1); hold on;

plot(score(all_contingencies{curr_animal} == 1 & ...
    ~all_latereversal{curr_animal},1), ...
    score(all_contingencies{curr_animal} == 1 ...
    & ~all_latereversal{curr_animal},2),'.k')

plot(score(all_contingencies{curr_animal} == 2 & ...
    ~all_latereversal{curr_animal},1), ...
    score(all_contingencies{curr_animal} == 2 ...
    & ~all_latereversal{curr_animal},2),'.r')

xlabel('PC1')
ylabel('PC2')
legend({'Contingency 1' 'Contingency 2'})

p2 = subplot(2,1,2); hold on;

plot(score(all_contingencies{curr_animal} == 1 & ...
    all_latereversal{curr_animal},1), ...
    score(all_contingencies{curr_animal} == 1 ...
    & all_latereversal{curr_animal},2),'ok')

plot(score(all_contingencies{curr_animal} == 2 & ...
    all_latereversal{curr_animal},1), ...
    score(all_contingencies{curr_animal} == 2 ...
    & all_latereversal{curr_animal},2),'or')

xlabel('PC1')
ylabel('PC2')
legend({'Contingency 1' 'Contingency 2'})

linkaxes([p1 p2],'xy');


% Get the distance between clusters and variance within cluster
cluster_dist = cell(length(mice),1);
for curr_animal = 1:length(mice)
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act);
    
    trial_pc_dist = pdist2(score(:,1:pcs),score(:,1:pcs));
    
    % Normalize distance
    trial_pc_dist_norm = (trial_pc_dist-mean(AP_itril(trial_pc_dist,-1))./std(AP_itril(trial_pc_dist,-1)));
    use_dist = trial_pc_dist_norm;
    
    % dim1: mice, dim2: within 1, within 2, across 1/2, dim3: early/late
    cluster_dist{curr_animal,1,1} = AP_itril(use_dist( ...
        all_contingencies{curr_animal} == 1 & ...
        ~all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 1 & ...
        ~all_latereversal{curr_animal}),-1);
    
    cluster_dist{curr_animal,2,1} = AP_itril(use_dist( ...
        all_contingencies{curr_animal} == 2 & ...
        ~all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 2 & ...
        ~all_latereversal{curr_animal}),-1);
    
    cluster_dist{curr_animal,3,1} = reshape(use_dist( ...
        all_contingencies{curr_animal} == 1 & ...
        ~all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 2 & ...
        ~all_latereversal{curr_animal}),[],1);
    
    cluster_dist{curr_animal,1,2} = AP_itril(use_dist( ...
        all_contingencies{curr_animal} == 1 & ...
        all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 1 & ...
        all_latereversal{curr_animal}),-1);
    
    cluster_dist{curr_animal,2,2} = AP_itril(use_dist( ...
        all_contingencies{curr_animal} == 2 & ...
        all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 2 & ...
        all_latereversal{curr_animal}),-1);
    
    cluster_dist{curr_animal,3,2} = reshape(use_dist( ...
        all_contingencies{curr_animal} == 1 & ...
        all_latereversal{curr_animal}, ...
        all_contingencies{curr_animal} == 2 & ...
        all_latereversal{curr_animal}),[],1);
    
end

cluster_dist_median = cellfun(@nanmedian,cluster_dist);

cluster_dist_alm_rev = cluster_dist_median(alm_animals & rev_animals,:,:);
cluster_dist_alm_rev_mean = permute(nanmean(cluster_dist_alm_rev,1),[3 2 1]);
cluster_dist_alm_rev_sem = permute(nanstd(cluster_dist_alm_rev)./ ...
    sqrt(sum(~isnan(cluster_dist_alm_rev))),[3 2 1]);

cluster_dist_pmm_rev = cluster_dist_median(pmm_animals & rev_animals,:,:);
cluster_dist_pmm_rev_mean = permute(nanmean(cluster_dist_pmm_rev,1),[3 2 1]);
cluster_dist_pmm_rev_sem = permute(nanstd(cluster_dist_pmm_rev)./ ...
    sqrt(sum(~isnan(cluster_dist_pmm_rev))),[3 2 1]);

figure;

subplot(2,1,1);
errorb(cluster_dist_alm_rev_mean,cluster_dist_alm_rev_sem);
colormap(gray);
set(gca,'XTickLabel',{'Early','Late'})
title('ALM');
legend({'1-1' '2-2' '1-2'});
ylabel('Median distance (trials)')

subplot(2,1,2);
errorb(cluster_dist_pmm_rev_mean,cluster_dist_pmm_rev_sem);
colormap(gray);
set(gca,'XTickLabel',{'Early','Late'})
title('PMM');
legend({'1-1' '2-2' '1-2'});
ylabel('Median distance (trials)')

% Alternatively: plot as boxplot
cluster_dist_alm_rev = cluster_dist_median(alm_animals & rev_animals,:,:);
cluster_dist_alm_rev_boxplot = permute(cluster_dist_alm_rev,[2 1 3]);

cluster_dist_pmm_rev = cluster_dist_median(pmm_animals & rev_animals,:,:);
cluster_dist_pmm_rev_boxplot = permute(cluster_dist_pmm_rev,[2 1 3]);

figure; 

subplot(2,1,1);
aboxplot(cluster_dist_alm_rev_boxplot,'colorgrad','green_down');
colormap(gray);
set(gca,'XTickLabel',{'Early','Late'})
title('ALM');
legend({'1-1' '2-2' '1-2'});
ylabel('Median normalized distance (trials)')

subplot(2,1,2);
aboxplot(cluster_dist_pmm_rev_boxplot,'colorgrad','green_down');
colormap(gray);
set(gca,'XTickLabel',{'Early','Late'})
title('PMM');
legend({'1-1' '2-2' '1-2'});
ylabel('Median normalized distance (trials)')


% Classify trial by contingency, early/late
% order: early contingency, late contingency, early/late
classifier_success = nan(length(mice),3);

% restrict the animals to run, one animal is "reversal" but only once!
for curr_animal = 1:length(animals);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act);
    
    % Classify by SVMsplit trials 80% train / 20% test groups 
    % (leaves out trials for symmetry of group number
    
    % Classify contingencies from each other (only if reversal animal)
    if ismember(curr_animal,find(rev_animals))
        
        % Early contingency
        curr_trials = find(~all_latereversal{curr_animal});
        curr_group = all_contingencies{curr_animal};
        
        num_test_trials = floor(length(curr_trials)*0.2);
        num_sets = floor(length(curr_trials)/num_test_trials);
        test_sets = mat2cell(curr_trials(randperm(length(curr_trials),num_test_trials*num_sets)), ...
            repmat(num_test_trials,1,num_sets));
        
        group_success = cell(length(test_sets),1);
        for curr_test_set = 1:length(test_sets)
            svmstruct = svmtrain(score(curr_trials,1:pcs),...
                curr_group(curr_trials));
            
            svm_group = svmclassify(svmstruct,score(vertcat(test_sets{setdiff(1: ...
                length(test_sets),curr_test_set)}),1:pcs));
            
            real_group = curr_group(vertcat(test_sets{setdiff(1: ...
                length(test_sets),curr_test_set)}));
            
            group_success{curr_test_set} = svm_group == real_group;
        end
        classifier_success(curr_animal,1) = nanmean(vertcat(group_success{:}));
        
        % Late contingencies (only if late reversals)
        if find(sum(diff(all_contingencies{curr_animal}) ~= 0) > 2)
            curr_trials = find(all_latereversal{curr_animal});
            curr_group = all_contingencies{curr_animal};
            
            num_test_trials = floor(length(curr_trials)*0.2);
            num_sets = floor(length(curr_trials)/num_test_trials);
            test_sets = mat2cell(curr_trials(randperm(length(curr_trials),num_test_trials*num_sets)), ...
                repmat(num_test_trials,1,num_sets));
            
            group_success = cell(length(test_sets),1);
            for curr_test_set = 1:length(test_sets)
                svmstruct = svmtrain(score(curr_trials,1:pcs),...
                    curr_group(curr_trials));
                
                svm_group = svmclassify(svmstruct,score(vertcat(test_sets{setdiff(1: ...
                    length(test_sets),curr_test_set)}),1:pcs));
                
                real_group = curr_group(vertcat(test_sets{setdiff(1: ...
                    length(test_sets),curr_test_set)}));
                
                group_success{curr_test_set} = svm_group == real_group;
            end
            classifier_success(curr_animal,2) = nanmean(vertcat(group_success{:}));
        end
    end
    
    % Classify early/late independent of contingency
    curr_trials = [1:size(score,1)]';
    curr_group = all_latereversal{curr_animal};
    
    num_test_trials = floor(length(curr_trials)*0.2);
    num_sets = floor(length(curr_trials)/num_test_trials);
    test_sets = mat2cell(curr_trials(randperm(length(curr_trials),num_test_trials*num_sets)), ...
        repmat(num_test_trials,1,num_sets));
    
    group_success = cell(length(test_sets),1);
    for curr_test_set = 1:length(test_sets)
        svmstruct = svmtrain(score(curr_trials,1:pcs),...
            curr_group(curr_trials));
        
        svm_group = svmclassify(svmstruct,score(vertcat(test_sets{setdiff(1: ...
            length(test_sets),curr_test_set)}),1:pcs));
        
        real_group = curr_group(vertcat(test_sets{setdiff(1: ...
            length(test_sets),curr_test_set)}));
        
        group_success{curr_test_set} = svm_group == real_group;
    end
    classifier_success(curr_animal,3) = nanmean(vertcat(group_success{:}));

end

% Plot classifier results
classifier_success_alm = classifier_success(alm_animals & rev_animals,:);
classifier_success_alm_mean = nanmean(classifier_success_alm);
classifier_success_alm_sem = nanstd(classifier_success_alm)./ ...
    sqrt(sum(~isnan(classifier_success_alm)));

classifier_success_pmm = classifier_success(pmm_animals & rev_animals,:);
classifier_success_pmm_mean = nanmean(classifier_success_pmm);
classifier_success_pmm_sem = nanstd(classifier_success_pmm)./ ...
    sqrt(sum(~isnan(classifier_success_pmm)));

classifier_plot_mean = [classifier_success_alm_mean;classifier_success_pmm_mean]';
classifier_plot_sem = [classifier_success_alm_sem;classifier_success_pmm_sem]';

figure;
subplot(1,3,1:2)
errorb(classifier_plot_mean,classifier_plot_sem);
colormap(gray);
ylim([0.5 1])
set(gca,'XTickLabel',{'Early contingency','Late contingency','Early/late'})
ylabel('Classifier success')
legend({'ALM' 'PMM'});
title('Reversal')

% Plot control early/late classification
classifier_success_alm_ctrl = classifier_success(alm_animals & ctrl_animals,3);
classifier_success_alm_ctrl_mean = nanmean(classifier_success_alm_ctrl);
classifier_success_alm_ctrl_sem = nanstd(classifier_success_alm_ctrl)./ ...
    sqrt(sum(~isnan(classifier_success_alm_ctrl)));

classifier_success_pmm_ctrl = classifier_success(pmm_animals & ctrl_animals,3);
classifier_success_pmm_ctrl_mean = nanmean(classifier_success_pmm_ctrl);
classifier_success_pmm_ctrl_sem = nanstd(classifier_success_pmm_ctrl)./ ...
    sqrt(sum(~isnan(classifier_success_pmm_ctrl)));

classifier_ctrl_plot_mean = [classifier_success_alm_ctrl_mean;classifier_success_pmm_ctrl_mean];
classifier_ctrl_plot_sem = [classifier_success_alm_ctrl_sem;classifier_success_pmm_ctrl_sem];

subplot(1,3,3); hold on;
bar(classifier_ctrl_plot_mean);
errorb(classifier_ctrl_plot_mean,classifier_ctrl_plot_sem);
colormap(gray);
ylim([0.5 1])
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'ALM','PMM'})
ylabel('Classifier success')
title('Control')


%% Slide 9) Temporal profiles of differences

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);

% Get all activity within animal
all_animal_act = cell(length(animals),1);

for curr_animal = 1:length(mice)
        
    curr_slot = 1;
    curr_animal_act = [];
    
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
            
            % only use first 2 seconds (odor delivery) of activity
            odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (:,curr_odor_frames,:);
                        
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(nanmean(any(curr_act,2),1),[3 1 2]);
            
            all_active_trials{curr_animal,curr_slot} = curr_active_trials;
            all_contingencies{curr_animal}(curr_slot) = ...
                analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);
            
            all_latereversal{curr_animal}(curr_slot) = ...
                sum(diff(all_contingencies{curr_animal}) ~= 0) > 1;
            
            curr_slot = curr_slot + 1;
            
            curr_animal_act = [curr_animal_act;nanmean(curr_act > 0,1)];
        end
    end
    all_animal_act{curr_animal} = curr_animal_act;
end

% Example animal plot
% Plot mean activity for contingencies, early/late
curr_animal = 9;
curr_act = horzcat(all_active_trials{curr_animal,:});
[coeff score latent] = princomp(curr_act');

[~,sort_idx] = sort(abs(coeff(:,1)));
active_cells = find(permute(any(nanmean(all_animal_act{curr_animal},1) > 0.01,2),[3 1 2]));
plot_cells = sort_idx(ismember(sort_idx,active_cells));


early_cont1 = permute(nanmean(all_animal_act{curr_animal}(~all_latereversal{curr_animal} & ...
    all_contingencies{curr_animal} == 1,:,plot_cells),1),[3 2 1]);

early_cont2 = permute(nanmean(all_animal_act{curr_animal}(~all_latereversal{curr_animal} & ...
    all_contingencies{curr_animal} == 2,:,plot_cells),1),[3 2 1]);

late_cont1 = permute(nanmean(all_animal_act{curr_animal}(all_latereversal{curr_animal} & ...
    all_contingencies{curr_animal} == 1,:,plot_cells),1),[3 2 1]);

late_cont2 = permute(nanmean(all_animal_act{curr_animal}(all_latereversal{curr_animal} & ...
    all_contingencies{curr_animal} == 2,:,plot_cells),1),[3 2 1]);

figure; colormap(gray);

subplot(3,3,1)
imagesc(early_cont1);
title('Contingency 1, Early')

subplot(3,3,2)
imagesc(early_cont2);
title('Contingency 2, Early')

subplot(3,3,4)
imagesc(late_cont1);
title('Contingency 1, Late')

subplot(3,3,5)
imagesc(late_cont2);
title('Contingency 2, Late')

% Plot subtractions of activity
subplot(3,3,3);
imagesc(early_cont2 - early_cont1);
title('Early 2 - Early 1');

subplot(3,3,6);
imagesc(late_cont2 - late_cont1);
title('Late 2 - Late 1');

subplot(3,3,7);
imagesc(late_cont1 - early_cont1);
title('Late 1 - Early 1');

subplot(3,3,8);
imagesc(late_cont2 - early_cont2);
title('Late 2 - Early 2');

% Plot the average temporal components of subtraction
figure;

subplot(2,2,1); hold on;
curr_sub = zscore(early_cont2 - early_cont1,[],2);
neg_sub = curr_sub;
neg_sub(neg_sub >= 0) = 0;
pos_sub = curr_sub;
pos_sub(pos_sub < 0) = 0;
plot(nanmean(neg_sub),'r');
plot(nanmean(pos_sub),'r');
plot(nanmean(curr_sub),'k');
line(xlim,[0 0],'color','k','linestyle','--');
title('Early 2 - Early 1')
ylabel('Difference of zscored activity');
xlabel('Time (frames)');

subplot(2,2,2); hold on;
curr_sub = zscore(late_cont2 - late_cont1,[],2);
neg_sub = curr_sub;
neg_sub(neg_sub >= 0) = 0;
pos_sub = curr_sub;
pos_sub(pos_sub < 0) = 0;
plot(nanmean(neg_sub),'r');
plot(nanmean(pos_sub),'r');
plot(nanmean(curr_sub),'k');
line(xlim,[0 0],'color','k','linestyle','--');
title('Late 2 - Late 1')
ylabel('Difference of zscored activity');
xlabel('Time (frames)');

subplot(2,2,3); hold on;
curr_sub = zscore(late_cont1 - early_cont1,[],2);
neg_sub = curr_sub;
neg_sub(neg_sub >= 0) = 0;
pos_sub = curr_sub;
pos_sub(pos_sub < 0) = 0;
plot(nanmean(neg_sub),'r');
plot(nanmean(pos_sub),'r');
plot(nanmean(curr_sub),'k');
line(xlim,[0 0],'color','k','linestyle','--');
title('Late 1 - Early 1')
ylabel('Difference of zscored activity');
xlabel('Time (frames)');

subplot(2,2,4); hold on;
curr_sub = zscore(late_cont2 - early_cont2,[],2);
neg_sub = curr_sub;
neg_sub(neg_sub >= 0) = 0;
pos_sub = curr_sub;
pos_sub(pos_sub < 0) = 0;
plot(nanmean(neg_sub),'r');
plot(nanmean(pos_sub),'r');
plot(nanmean(curr_sub),'k');
line(xlim,[0 0],'color','k','linestyle','--');
title('Late 2 - Early 2')
ylabel('Difference of zscored activity');
xlabel('Time (frames)');


% Average temporal difference profiles

% Example animal plot
% Plot mean activity for contingencies, early/late
temporal_difference = cell(length(mice),1);
for curr_animal = 1:length(mice);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    [coeff score latent] = princomp(curr_act');    
    
    early_cont1 = permute(nanmean(all_animal_act{curr_animal}(~all_latereversal{curr_animal} & ...
        all_contingencies{curr_animal} == 1,:,:),1),[3 2 1]);
    
    early_cont2 = permute(nanmean(all_animal_act{curr_animal}(~all_latereversal{curr_animal} & ...
        all_contingencies{curr_animal} == 2,:,:),1),[3 2 1]);
    
    late_cont1 = permute(nanmean(all_animal_act{curr_animal}(all_latereversal{curr_animal} & ...
        all_contingencies{curr_animal} == 1,:,:),1),[3 2 1]);
    
    late_cont2 = permute(nanmean(all_animal_act{curr_animal}(all_latereversal{curr_animal} & ...
        all_contingencies{curr_animal} == 2,:,:),1),[3 2 1]);
    
    e2_e1 = nanmean(zscore(early_cont2 - early_cont1,[],2));
    l2_l1 = nanmean(zscore(late_cont2 - late_cont1,[],2));
    l1_e1 = nanmean(zscore(late_cont1 - early_cont1,[],2));
    l2_e2 = nanmean(zscore(late_cont2 - early_cont2,[],2));
    
    temporal_difference{curr_animal} = [e2_e1;l2_l1;l1_e1;l2_e2];
    
end

% ALM
temporal_difference_alm_rev = cat(3,temporal_difference{alm_animals & rev_animals});
temporal_difference_alm_rev_mean = nanmean(temporal_difference_alm_rev,3);
temporal_difference_alm_rev_sem = nanstd(temporal_difference_alm_rev,[],3)./ ...
    sqrt(sum(~isnan(temporal_difference_alm_rev),3));

figure; set(gcf,'Name','ALM');
for i = 1:4
   subplot(2,2,i); hold on
   plot(temporal_difference_alm_rev_mean(i,:),'k');
   plot(temporal_difference_alm_rev_mean(i,:) - ...
       temporal_difference_alm_rev_sem(i,:),'--k');
   plot(temporal_difference_alm_rev_mean(i,:) + ...
       temporal_difference_alm_rev_sem(i,:),'--k');
   
   line(xlim,[0 0],'color','r');
end

% PMM
% for now exclude the bscope animal, can downsample later or something
scope = [mice.scope];

temporal_difference_pmm_rev = cat(3,temporal_difference{pmm_animals & rev_animals & scope == 1});
temporal_difference_pmm_rev_mean = nanmean(temporal_difference_pmm_rev,3);
temporal_difference_pmm_rev_sem = nanstd(temporal_difference_pmm_rev,[],3)./ ...
    sqrt(sum(~isnan(temporal_difference_pmm_rev),3));

figure; set(gcf,'Name','PMM');
for i = 1:4
   subplot(2,2,i); hold on
   plot(temporal_difference_pmm_rev_mean(i,:),'k');
   plot(temporal_difference_pmm_rev_mean(i,:) - ...
       temporal_difference_pmm_rev_sem(i,:),'--k');
   plot(temporal_difference_pmm_rev_mean(i,:) + ...
       temporal_difference_pmm_rev_sem(i,:),'--k');
   
   line(xlim,[0 0],'color','r');
end













































