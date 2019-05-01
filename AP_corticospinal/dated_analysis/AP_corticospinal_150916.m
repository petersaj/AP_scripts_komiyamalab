%% Classify ROIs

classified_rois = AP_classify_movement_cells_continuous(data,analysis);


%% Compare movement to activity correlation / move-act corr slope

% requires: behavior analysis, 
% across day move-act corr (prelim_move)
% avg act corr (150527)

figure; 

subplot(3,3,1);
imagesc(nanmean(movement_corr_grid,3));
colormap(hot);
title('Movement')

subplot(3,3,2);
imagesc(nanmean(act_aligned_corr_move,3));
colormap(hot);
title('Aligned activity')

subplot(3,3,3);
imagesc(corr_slope);
colormap(hot);
title('Act-move slope');

move_corr = nanmean(movement_corr_grid,3);

actmove_corrslope_tril = AP_itril(corr_slope,-1);
move_corr_tril = AP_itril(nanmean(movement_corr_grid,3),-1);
act_corr_tril = AP_itril(nanmean(act_aligned_corr_move,3),-1);

subplot(3,3,4);
plot(move_corr_tril,actmove_corrslope_tril,'.k');
xlabel('Move corr');
ylabel('Act-move slope')

subplot(3,3,5);
plot(move_corr_tril,act_corr_tril,'.k');
xlabel('Move corr');
ylabel('Act corr')

subplot(3,3,6);
plot(actmove_corrslope_tril,act_corr_tril,'.k');
xlabel('Act-move slope');
ylabel('Act corr')

% Bin movement corr for plotting
move_corr_bin = linspace(min(move_corr_tril),max(move_corr_tril),10);
move_corr_bincenters = move_corr_bin(1:end-1) + diff(move_corr_bin)/2;
move_corr_bin(end) = Inf;
[~,bins] = histc(move_corr_tril,move_corr_bin);
act_meanci = grpstats(act_corr_tril,bins,'std');

% Correlation between move-act corr slope and movement correlation
curr_r = corrcoef(move_corr_tril,actmove_corrslope_tril);
r = curr_r(2);

n_shuff = 1000;
r_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    a = AP_itril(shake(nanmean(movement_corr_grid,3),'diag'),-1);
    b = AP_itril(shake(corr_slope,'diag'),-1);
    curr_r = corrcoef(a,b);
    r_shuff(curr_shuff) = curr_r(2);
    
    disp(curr_shuff);
end

r_rank = tiedrank([r;r_shuff]);
r_p = r_rank(1)/length([r;r_shuff]);

subplot(3,3,7); hold on;
hist(r_shuff,20);
line([r,r],ylim,'color','r')
title('Move corr, act-move slope')

% Correlation between avg act correlation and movement correlation
curr_r = corrcoef(move_corr_tril,act_corr_tril);
r = curr_r(2);
n_shuff = 1000;
r_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    a = AP_itril(shake(nanmean(movement_corr_grid,3),'diag'),-1);
    b = AP_itril(shake(nanmean(act_aligned_corr_move,3),'diag'),-1);    
    curr_r = corrcoef(a,b);
    r_shuff(curr_shuff) = curr_r(2);   
    
    disp(curr_shuff);
end

r_rank = tiedrank([r;r_shuff]);
r_p = r_rank(1)/length([r;r_shuff]);

subplot(3,3,8); hold on;
hist(r_shuff,20);
line([r,r],ylim,'color','r')
title('Move corr, act corr')

% Correlation between avg act correlation and move-act corr slope
curr_r = corrcoef(actmove_corrslope_tril,act_corr_tril);
r = curr_r(2);

n_shuff = 1000;
r_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    a = shake(corr_slope,'diag');
    b = shake(nanmean(act_aligned_corr_move,3),'diag');
    curr_r = corrcoef(AP_itril(a,-1),AP_itril(b,-1));
    r_shuff(curr_shuff) = curr_r(2);
    disp(curr_shuff);
end

r_rank = tiedrank([r;r_shuff]);
r_p = r_rank(1)/length([r;r_shuff]);

subplot(3,3,9); hold on;
hist(r_shuff,20);
line([r,r],ylim,'color','r')
title('Act-move slope, act corr')

set(gcf,'Name','Diagonal shuffle r')

%% Quiescent activity over days
% Do all the changes keep a constant baseline activity?

avg_movement_act = cell(length(data),1);
avg_quiescent_act = cell(length(data),1);

for curr_animal = 1:length(data)
        
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    
    curr_avg_movement_act = nan(n_rois,14);
    curr_avg_quiescent_act = nan(n_rois,14);
    
     for i = 1:14
        curr_avg_movement_act(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,analysis(curr_animal).lever(i).lever_move_frames),2);
        curr_avg_quiescent_act(:,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,~analysis(curr_animal).lever(i).lever_move_frames),2);
     end
    
    avg_movement_act{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act{curr_animal} = curr_avg_quiescent_act;
        
    disp(curr_animal);
end


figure;

subplot(2,2,1); hold on;
curr_act = cell2mat(cellfun(@nanmean,avg_movement_act,'uni',false));
plot(curr_act');
errorbar(nanmean(curr_act),nanstd(curr_act)./sqrt(sum(~isnan(curr_act))),'k','linewidth',2)
title('Average movement activity');
xlabel('Day')
ylabel('Average \DeltaF/F')

subplot(2,2,2); hold on;
curr_act = cell2mat(cellfun(@nanmean,avg_quiescent_act,'uni',false));
plot(curr_act');
errorbar(nanmean(curr_act),nanstd(curr_act)./sqrt(sum(~isnan(curr_act))),'k','linewidth',2)
title('Average quiescent activity')
xlabel('Day')
ylabel('Average \DeltaF/F')

subplot(2,2,3); hold on;
curr_act = cell2mat(cellfun(@nanmean,avg_movement_act,'uni',false));
curr_act = bsxfun(@times,curr_act,1./curr_act(:,1));
plot(curr_act');
errorbar(nanmean(curr_act),nanstd(curr_act)./sqrt(sum(~isnan(curr_act))),'k','linewidth',2)
title('Normalized movement activity');
xlabel('Day')
ylabel('Average \DeltaF/F')

subplot(2,2,4); hold on;
curr_act = cell2mat(cellfun(@nanmean,avg_quiescent_act,'uni',false));
curr_act = bsxfun(@times,curr_act,1./curr_act(:,1));
plot(curr_act');
errorbar(nanmean(curr_act),nanstd(curr_act)./sqrt(sum(~isnan(curr_act))),'k','linewidth',2)
title('Normalized quiescent activity')
xlabel('Day')
ylabel('Average \DeltaF/F')


%% Quiescent activity over time/trials

avg_movement_act = cell(length(data),1);
avg_quiescent_act = cell(length(data),1);

for curr_animal = 1:length(data)
        
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    
    curr_avg_movement_act = nan(n_rois,14);
    curr_avg_quiescent_act = nan(n_rois,14);
    
    for i = 1:14
        curr_m = classified_rois(curr_animal).movement(:,i);
        curr_q = classified_rois(curr_animal).quiescent(:,i);
        
        
        curr_avg_movement_act(curr_m,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(curr_m,analysis(curr_animal).lever(i).lever_move_frames),2);
        curr_avg_quiescent_act(curr_q,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(curr_q,~analysis(curr_animal).lever(i).lever_move_frames),2);
    end
    
    avg_movement_act{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act{curr_animal} = curr_avg_quiescent_act;
        
    disp(curr_animal);
end

m_q_r = nan(curr_animal,curr_day);
m_q_p = nan(curr_animal,curr_day);
for curr_animal = 1:8   
    
    %figure;
    
    for curr_day = 1:14;
        
        % Get quiescent activity that's not within ~1 sec of movement
        move_leeway = 1;
        move_leeway_filt = ones(1,move_leeway);
        
        curr_move = analysis(curr_animal).lever(curr_day).lever_move_frames;
        curr_move_leeway = conv(+curr_move,move_leeway_filt,'same') > 0;
        
        % curr_quiescent_act = ...
        %     nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(:,~curr_move_leeway,1));
        %
        % curr_move_act = ...
        %     nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(:,curr_move_leeway,1));
        
        curr_quiescent_act = ...
            nanmean(data(curr_animal).im(curr_day).roi_trace_thresh( ...
            classified_rois(curr_animal).quiescent(:,curr_day),~curr_move_leeway,1));
        
        curr_move_act = ...
            nanmean(data(curr_animal).im(curr_day).roi_trace_thresh( ...
            classified_rois(curr_animal).movement(:,curr_day),curr_move_leeway,1));
        
        % smooth activity
        curr_q_smooth = smooth(curr_quiescent_act,5000);
        curr_m_smooth = smooth(curr_move_act,5000);
        
        % Interpolate across gaps for move/quiescent
        curr_q_interp = nan(size(curr_move_leeway));
        curr_q_interp(~curr_move_leeway) = curr_q_smooth;
        curr_q_interp(1) = curr_q_interp(find(~curr_move_leeway,1));
        curr_q_interp(end) = curr_q_interp(find(~curr_move_leeway,1,'last'));
        curr_q_interp(isnan(curr_q_interp)) = interp1(find(~isnan(curr_q_interp)), ...
            curr_q_interp(~isnan(curr_q_interp)),find(isnan(curr_q_interp)));

        curr_m_interp = nan(size(curr_move_leeway));
        curr_m_interp(curr_move_leeway) = curr_m_smooth;
        curr_m_interp(1) = curr_m_interp(find(curr_move_leeway,1));
        curr_m_interp(end) = curr_m_interp(find(curr_move_leeway,1,'last'));
        curr_m_interp(isnan(curr_m_interp)) = interp1(find(~isnan(curr_m_interp)), ...
            curr_m_interp(~isnan(curr_m_interp)),find(isnan(curr_m_interp)));
        
        [curr_r,curr_p] = corrcoef(curr_q_interp,curr_m_interp);
        
        m_q_r(curr_animal,curr_day) = curr_r(2);
        m_q_p(curr_animal,curr_day) = curr_p(2);  
        
        %subplot(4,4,curr_day); hold on;
        %plot(curr_q_interp,'k');
        %plot(curr_m_interp,'r');
        
    end
end
disp('done')

% no figure from this yet: doesn't look like it goes down over the day but
% does peak after movements, not clear if coming back to "baseline" or
% actual peak in activity


%% Fano factor / variability around movement

%%% Within ROIs

% Movement onset fano factor: all ROIs, movement rois, quiescent rois
ff_onset_all = nan(14,181,8);
ff_onset_m = nan(14,181,8);
ff_onset_q = nan(14,181,8);
for curr_animal = 1:length(analysis)
    
    curr_ff = nan(size(analysis(curr_animal).im(1).move_onset_aligned,3), ...
        size(analysis(curr_animal).im(1).move_onset_aligned,2),14);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        curr_data = analysis(curr_animal).im(curr_day).move_onset_aligned;

        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);        
        
        ff_onset_all(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,:),1)./ ...
            mean(curr_data(:,:,:),1),[3 2 1]),1);
        
        ff_onset_m(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,m_rois),1)./ ...
            mean(curr_data(:,:,m_rois),1),[3 2 1]),1);
        
        ff_onset_q(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,q_rois),1)./ ...
            mean(curr_data(:,:,q_rois),1),[3 2 1]),1);
        
    end
    
end
ff_onset_all_avg = nanmean(bsxfun(@times,ff_onset_all,1./nanmedian(ff_onset_all(1,1:10,:),2)),3);
ff_onset_m_avg = nanmean(bsxfun(@times,ff_onset_m,1./nanmedian(ff_onset_m(1,1:10,:),2)),3);
ff_onset_q_avg = nanmean(bsxfun(@times,ff_onset_q,1./nanmedian(ff_onset_q(1,1:10,:),2)),3);

% Movement offset fano factor: all ROIs, movement rois, quiescent rois
ff_offset_all = nan(14,181,8);
ff_offset_m = nan(14,181,8);
ff_offset_q = nan(14,181,8);
for curr_animal = 1:length(analysis)
    
    curr_ff = nan(size(analysis(curr_animal).im(1).move_offset_aligned,3), ...
        size(analysis(curr_animal).im(1).move_offset_aligned,2),14);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        curr_data = analysis(curr_animal).im(curr_day).move_offset_aligned;

        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);        
        
        ff_offset_all(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,:),1)./ ...
            mean(curr_data(:,:,:),1),[3 2 1]),1);
        
        ff_offset_m(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,m_rois),1)./ ...
            mean(curr_data(:,:,m_rois),1),[3 2 1]),1);
        
        ff_offset_q(curr_day,:,curr_animal) = ...
            nanmean(permute(var(curr_data(:,:,q_rois),1)./ ...
            mean(curr_data(:,:,q_rois),1),[3 2 1]),1);
        
    end
    
end
ff_offset_all_avg = nanmean(bsxfun(@times,ff_offset_all,1./nanmedian(ff_offset_all(1,1:10,:),2)),3);
ff_offset_m_avg = nanmean(bsxfun(@times,ff_offset_m,1./nanmedian(ff_offset_m(1,1:10,:),2)),3);
ff_offset_q_avg = nanmean(bsxfun(@times,ff_offset_q,1./nanmedian(ff_offset_q(1,1:10,:),2)),3);


%%% Across ROIs

% Movement onset fano factor: all ROIs, movement rois, quiescent rois
ff_onset_all_trialavg = nan(14,181,8);
ff_onset_m_trialavg = nan(14,181,8);
ff_onset_q_trialavg = nan(14,181,8);
for curr_animal = 1:length(analysis)
    
    for curr_day = 1:length(data(curr_animal).im)
        
        curr_data = analysis(curr_animal).im(curr_day).move_onset_aligned;

        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);    
        
        avg_all_trialact = nanmean(curr_data,3);
        avg_m_trialact = nanmean(curr_data(:,:,m_rois),3);
        avg_q_trialact = nanmean(curr_data(:,:,q_rois),3);
                
        ff_onset_all_trialavg(curr_day,:,curr_animal) = ...
            var(avg_all_trialact,1)./ ...
            mean(avg_all_trialact,1);
        
        ff_onset_m_trialavg(curr_day,:,curr_animal) = ...
            var(avg_m_trialact,1)./ ...
            mean(avg_m_trialact,1);
        
        ff_onset_q_trialavg(curr_day,:,curr_animal) = ...
            var(avg_q_trialact,1)./ ...
            mean(avg_q_trialact,1);
        
    end
    
end
ff_onset_all_trialavg_avg = nanmean(bsxfun(@times,ff_onset_all_trialavg,1./nanmedian(ff_onset_all_trialavg(1,1:10,:),2)),3);
ff_onset_m_trialavg_avg = nanmean(bsxfun(@times,ff_onset_m_trialavg,1./nanmedian(ff_onset_m_trialavg(1,1:10,:),2)),3);
ff_onset_q_trialavg_avg = nanmean(bsxfun(@times,ff_onset_q_trialavg,1./nanmedian(ff_onset_q_trialavg(1,1:10,:),2)),3);

% Movement offset fano factor: all ROIs, movement rois, quiescent rois
ff_offset_all_trialavg = nan(14,181,8);
ff_offset_m_trialavg = nan(14,181,8);
ff_offset_q_trialavg = nan(14,181,8);
for curr_animal = 1:length(analysis)
    
    for curr_day = 1:length(data(curr_animal).im)
        
        curr_data = analysis(curr_animal).im(curr_day).move_offset_aligned;

        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);    
        
        avg_all_trialact = nanmean(curr_data,3);
        avg_m_trialact = nanmean(curr_data(:,:,m_rois),3);
        avg_q_trialact = nanmean(curr_data(:,:,q_rois),3);       
        
        ff_offset_all_trialavg(curr_day,:,curr_animal) = ...
            var(avg_all_trialact,1)./ ...
            mean(avg_all_trialact,1);
        
        ff_offset_m_trialavg(curr_day,:,curr_animal) = ...
            var(avg_m_trialact,1)./ ...
            mean(avg_m_trialact,1);
        
        ff_offset_q_trialavg(curr_day,:,curr_animal) = ...
            var(avg_q_trialact,1)./ ...
            mean(avg_q_trialact,1);
        
    end
    
end
ff_offset_all_trialavg_avg = nanmean(bsxfun(@times,ff_offset_all_trialavg,1./nanmedian(ff_offset_all_trialavg(1,1:10,:),2)),3);
ff_offset_m_trialavg_avg = nanmean(bsxfun(@times,ff_offset_m_trialavg,1./nanmedian(ff_offset_m_trialavg(1,1:10,:),2)),3);
ff_offset_q_trialavg_avg = nanmean(bsxfun(@times,ff_offset_q_trialavg,1./nanmedian(ff_offset_q_trialavg(1,1:10,:),2)),3);

% Plot fano factor of roi classes over days 

x = ((1:181) - 91)./28;

%%% Within ROIs

figure;
col = jet(14);

subplot(2,3,1); hold on;
for i = 1:14
plot(x,ff_onset_all_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged all ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,2); hold on;
for i = 1:14
plot(x,ff_onset_m_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged m ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,3); hold on;
for i = 1:14
plot(x,ff_onset_q_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged q ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,4); hold on;
for i = 1:14
plot(x,ff_offset_all_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (all ROIs)');
ylabel('Fano factor (averaged normalized)');

subplot(2,3,5); hold on;
for i = 1:14
plot(x,ff_offset_m_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (averaged m ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,6); hold on;
for i = 1:14
plot(x,ff_offset_q_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (averaged q ROIs)');
ylabel('Fano factor (normalized)');

%%% Across ROIs

figure;
col = jet(14);

subplot(2,3,1); hold on;
for i = 1:14
plot(x,ff_onset_all_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged all ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,2); hold on;
for i = 1:14
plot(x,ff_onset_m_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged m ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,3); hold on;
for i = 1:14
plot(x,ff_onset_q_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (averaged q ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,4); hold on;
for i = 1:14
plot(x,ff_offset_all_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (averaged all ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,5); hold on;
for i = 1:14
plot(x,ff_offset_m_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (averaged m ROIs)');
ylabel('Fano factor (normalized)');

subplot(2,3,6); hold on;
for i = 1:14
plot(x,ff_offset_q_trialavg_avg(i,:),'color',col(i,:),'linewidth',2);
end
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (averaged q ROIs)');
ylabel('Fano factor (normalized)');

% Plot average fano factor within roi class, compare to mean activity
figure;

subplot(3,2,1); hold on;
plot(x,nanmean(nanmean(ff_onset_all),3),'k')
plot(x,nanmean(nanmean(ff_onset_m),3),'color',[0,0.8,0])
plot(x,nanmean(nanmean(ff_onset_q),3),'r')
ylabel('Fano factor');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset');

subplot(3,2,2); hold on;
plot(x,nanmean(nanmean(ff_offset_all),3),'k')
plot(x,nanmean(nanmean(ff_offset_m),3),'color',[0,0.8,0])
plot(x,nanmean(nanmean(ff_offset_q),3),'r')
ylabel('Fano factor');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset');

subplot(3,2,3); hold on;
plot(x,nanmean(nanmean(ff_onset_all_trialavg),3),'k')
plot(x,nanmean(nanmean(ff_onset_m_trialavg),3),'color',[0,0.8,0])
plot(x,nanmean(nanmean(ff_onset_q_trialavg),3),'r')
ylabel('Fano factor');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset (ROI avg)');

subplot(3,2,4); hold on;
plot(x,nanmean(nanmean(ff_offset_all_trialavg),3),'k')
plot(x,nanmean(nanmean(ff_offset_m_trialavg),3),'color',[0,0.8,0])
plot(x,nanmean(nanmean(ff_offset_q_trialavg),3),'r')
ylabel('Fano factor');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset (ROI avg)');

% Movement onset fano factor: all ROIs, movement rois, quiescent rois
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

subplot(3,2,5); hold on;
plot(x,nanmean(vertcat(avg_act_onset_all{:})),'k');
plot(x,nanmean(vertcat(avg_act_onset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_onset_q{:})),'r');
ylabel('Average \DeltaF/F');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement onset');

subplot(3,2,6); hold on;
plot(x,nanmean(vertcat(avg_act_offset_all{:})),'k');
plot(x,nanmean(vertcat(avg_act_offset_m{:})),'color',[0,0.8,0]);
plot(x,nanmean(vertcat(avg_act_offset_q{:})),'r');
ylabel('Average \DeltaF/F');
legend({'All' 'M' 'Q'},'location','sw')
line([0,0],ylim,'color','k','linestyle','--');
title('Movement offset');



%% Reliability by class

% get reliability by max time
max_reliability_all = cell(8,1);
max_reliability_m = cell(8,1);
max_reliability_q = cell(8,1);
for curr_animal = 1:8;
    for curr_day = 1:14
        
        m_rois = classified_rois(curr_animal).movement(:,curr_day);
        q_rois = classified_rois(curr_animal).quiescent(:,curr_day);
              
        curr_act = ...
            permute(nanmean(analysis(curr_animal).im(curr_day).move_onset_aligned(:,:,:) > 0,1),[3 2 1]);        
        
        max_reliability_all{curr_animal}(curr_day) = nanmean(max(curr_act(:,:),[],2));
        max_reliability_m{curr_animal}(curr_day) = nanmean(max(curr_act(m_rois,:),[],2));
        max_reliability_q{curr_animal}(curr_day) = nanmean(max(curr_act(q_rois,:),[],2));
    end
end

figure; hold on;
plot(nanmean(cell2mat(max_reliability_all)),'k');
plot(nanmean(cell2mat(max_reliability_m)),'g');
plot(nanmean(cell2mat(max_reliability_q)),'r');

ylabel('Average max reliability');
legend({'All','M','Q'});
xlabel('Day');



%% Classification turnover FUTURE
% This compliments the one in the last script 

% Get fraction of currently classified ROIs that will be classified future
m2m = nan(length(data),14);
for curr_animal = 1:length(data)
m2m(curr_animal,:) = ...
    (sum((fliplr(cumsum(fliplr(classified_rois(curr_animal).movement),2)) > 1).* ...
    classified_rois(curr_animal).movement))./ ...
    sum(classified_rois(curr_animal).movement);
end
q2q = nan(length(data),14);
for curr_animal = 1:length(data)
q2q(curr_animal,:) = ...
    (sum((fliplr(cumsum(fliplr(classified_rois(curr_animal).quiescent),2)) > 1).* ...
    classified_rois(curr_animal).quiescent))./ ...
    sum(classified_rois(curr_animal).quiescent);
end
q2m = nan(length(data),14);
for curr_animal = 1:length(data)
q2m(curr_animal,:) = ...
    (sum((fliplr(cumsum(fliplr(classified_rois(curr_animal).movement),2)) >= 1).* ...
    classified_rois(curr_animal).quiescent))./ ...
    sum(classified_rois(curr_animal).quiescent);
end
m2q = nan(length(data),14);
for curr_animal = 1:length(data)
m2q(curr_animal,:) = ...
    (sum((fliplr(cumsum(fliplr(classified_rois(curr_animal).quiescent),2)) >= 1).* ...
    classified_rois(curr_animal).movement))./ ...
    sum(classified_rois(curr_animal).movement);
end
    
    
n_rep = 5000;

m2m_shuff_1 = nan(length(data),14,n_rep);
q2q_shuff_1 = nan(length(data),14,n_rep);
m2q_shuff_1 = nan(length(data),14,n_rep);
q2m_shuff_1 = nan(length(data),14,n_rep);

m2m_shuff_2 = nan(length(data),14,n_rep);
q2q_shuff_2 = nan(length(data),14,n_rep);
m2q_shuff_2 = nan(length(data),14,n_rep);
q2m_shuff_2 = nan(length(data),14,n_rep);

% Get shuffled values of fraction previously classified
for curr_rep = 1:n_rep
    shuff_classified_movement = arrayfun(@(x) shake(classified_rois(x).movement,1),1:length(data),'uni',false);
    shuff_classified_quiescent = arrayfun(@(x) shake(classified_rois(x).quiescent,1),1:length(data),'uni',false);
    
    % Get fraction of currently classified ROIs that were classified previously
    for curr_animal = 1:length(data)
        m2m_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_movement{curr_animal}),2)) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2q_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_quiescent{curr_animal}),2)) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2m_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_movement{curr_animal}),2)) >= 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        m2q_shuff_1(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_quiescent{curr_animal}),2)) >= 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    
    
    shuff_classified_movement = arrayfun(@(x) shake(classified_rois(x).movement,2),1:length(data),'uni',false);
    shuff_classified_quiescent = arrayfun(@(x) shake(classified_rois(x).quiescent,2),1:length(data),'uni',false);
    
    % Get fraction of currently classified ROIs that were classified previously
    for curr_animal = 1:length(data)
        m2m_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_movement{curr_animal}),2)) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2q_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_quiescent{curr_animal}),2)) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        q2m_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_movement{curr_animal}),2)) > 1).* ...
            shuff_classified_quiescent{curr_animal}))./ ...
            sum(shuff_classified_quiescent{curr_animal});
    end
    for curr_animal = 1:length(data)
        m2q_shuff_2(curr_animal,:,curr_rep) = ...
            (sum((fliplr(cumsum(fliplr(shuff_classified_quiescent{curr_animal}),2)) > 1).* ...
            shuff_classified_movement{curr_animal}))./ ...
            sum(shuff_classified_movement{curr_animal});
    end
    disp(curr_rep)
end

% Plot fraction of currently classified ROIs previously classified
figure; 
subplot(2,2,1);hold on;
errorbar(nanmean(m2m), ...
    nanstd(m2m)./sqrt(sum(~isnan(m2m))),'k','linewidth',2)
plot(permute(prctile(nanmean(m2m_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(m2m_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')

xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Movement to movement');
subplot(2,2,2);hold on;
errorbar(nanmean(q2q), ...
    nanstd(q2q)./sqrt(sum(~isnan(q2q))),'k','linewidth',2)
plot(permute(prctile(nanmean(q2q_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(q2q_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Quiescent to quiescent');
subplot(2,2,3);hold on;
errorbar(nanmean(m2q), ...
    nanstd(m2q)./sqrt(sum(~isnan(m2q))),'k','linewidth',2)
plot(permute(prctile(nanmean(m2q_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(m2q_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Movement to quiescent');
subplot(2,2,4);hold on;
errorbar(nanmean(q2m), ...
    nanstd(q2m)./sqrt(sum(~isnan(q2m))),'k','linewidth',2)
plot(permute(prctile(nanmean(q2m_shuff_1,1),[2.5 97.5],3),[2 3 1]),'r')
plot(permute(prctile(nanmean(q2m_shuff_2,1),[2.5 97.5],3),[2 3 1]),'b')
xlabel('Day');
ylabel('Fraction of classified ROIs');
title('Quiescent to movement');
    
    
set(gcf,'Name','% future classified');
    
%% Classification overlap and switching

% Plot the fraction of classified ROIs on one day overlapping with another
m2m_overlap_grid = nan(14,14,length(data));
q2q_overlap_grid = nan(14,14,length(data));
m2q_overlap_grid = nan(14,14,length(data));
q2m_overlap_grid = nan(14,14,length(data));
for curr_animal = 1:length(data)
   for curr_day1 = 1:14;
       for curr_day2 = 1:14;
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


% Plot overlap: of M NOT overlap, how many Q?
m2q_nonm_overlap_grid = nan(14,14,length(data));
for curr_animal = 1:length(data)
   for curr_day1 = 1:14;
       for curr_day2 = 1:14;
           
           curr_nonm = classified_rois(curr_animal).movement(:,curr_day1) & ...
               ~classified_rois(curr_animal).movement(:,curr_day2);
           
           curr_m2q = classified_rois(curr_animal).movement(:,curr_day1) & ...
               classified_rois(curr_animal).quiescent(:,curr_day2);
           
           m2q_nonm_overlap_grid(curr_day1,curr_day2,curr_animal) = sum(curr_m2q)./sum(curr_nonm);
           
       end
   end
end
figure;
imagesc(nanmean(m2q_nonm_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('M2Q/non-M')
colorbar;


% Plot the Q and M/Q ratio with and without the M2Q ROIs
no_m2q_q = nan(length(data),14);
no_m2q_m = nan(length(data),14);
for curr_animal = 1:length(data)
    curr_m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    
    curr_q = classified_rois(curr_animal).quiescent;
    curr_q(curr_m2q,:) = false;
    
    curr_m = classified_rois(curr_animal).movement;
    curr_m(curr_m2q,:) = false;
    
    no_m2q_q(curr_animal,:) = nanmean(curr_q);
    no_m2q_m(curr_animal,:) = nanmean(curr_m);
end

all_q = cell2mat(cellfun(@nanmean,{classified_rois(:).quiescent}','uni',false));
all_m = cell2mat(cellfun(@nanmean,{classified_rois(:).movement}','uni',false));

all_mq_ratio = all_m./all_q;
no_m2q_mq_ratio = all_m./no_m2q_q;
no_m2qever_mq_ratio = no_m2q_m./no_m2q_q;

figure; 

subplot(2,1,1); hold on;
errorbar(nanmean(all_q),nanstd(all_q)./sqrt(sum(~isnan(all_q))),'k','linewidth',2);
errorbar(nanmean(no_m2q_q),nanstd(no_m2q_q)./sqrt(sum(~isnan(no_m2q_q))),'r','linewidth',2);
ylabel('Fraction ROIs')
xlabel('Day')
legend({'All' 'w/o m2q'});
title('Quiescent classified')

subplot(2,1,2); hold on;
errorbar(nanmean(all_mq_ratio),nanstd(all_mq_ratio)./sqrt(sum(~isnan(all_mq_ratio))),'k','linewidth',2);
errorbar(nanmean(no_m2q_mq_ratio),nanstd(no_m2q_mq_ratio)./sqrt(sum(~isnan(no_m2q_mq_ratio))),'r','linewidth',2);
errorbar(nanmean(no_m2qever_mq_ratio),nanstd(no_m2qever_mq_ratio)./sqrt(sum(~isnan(no_m2qever_mq_ratio))),'b','linewidth',2);
ylabel('M/Q ratio')
xlabel('Day')
legend({'All' 'w/o m2q q' 'w/o m2q m/q'});


% Plot fraction of M2Q as a function of M/Q change
mq_ratio_change = max(all_mq_ratio,[],2)./all_mq_ratio(:,1);
m2q = cellfun(@(m,q) any((cumsum(m,2) >= 1).*q,2), ...
     {classified_rois(:).movement},{classified_rois(:).quiescent},'uni',false);
m2q_frac = cellfun(@nanmean,m2q);
figure;plot(mq_ratio_change,m2q_frac,'.k');
    

% Of unique dropped M ROIs, how many ever Q?
dropped_m2q = nan(length(data),1);
for curr_animal = 1:length(data)
    
    early_m = any(classified_rois(curr_animal).movement(:,1:7),2);
    late_m = any(classified_rois(curr_animal).movement(:,8:14),2);
    late_q = any(classified_rois(curr_animal).quiescent(:,8:14),2);
    
    dropped_m2q(curr_animal) = sum(early_m & ~late_m & late_q)/sum(early_m & ~late_m);
    
end


%% Fraction of ROIs / mROIs active per movement vs. frac classified

act_movetime_binned = nan(14,4,length(data));
act_per_movetime = nan(length(data),14);
act_per_move = nan(length(data),14);
movetime = nan(length(data),14);
for curr_animal = 1:length(data)
    
    for curr_session = 1:14;
        
    curr_m = classified_rois(curr_animal).movement(:,curr_session);
    curr_q = classified_rois(curr_animal).quiescent(:,curr_session);
    
    % Parse movements    
    curr_lever_move = analysis(curr_animal).lever(curr_session).lever_move_frames;
    move_frames = find(curr_lever_move);
    move_break = [0;find(diff(move_frames) > 1);length(move_frames)];
    move_durations = diff(move_break);
    move_frames_split = mat2cell(move_frames,move_durations);
    
    % Use only cued-rewarded movements
    reward_frames = cellfun(@(x) x.states.reward(1), ...
        data(curr_animal).bhv(curr_session).bhv_frames( ...
        analysis(curr_animal).lever(curr_session).cued_movement_trials));  
    cr_movements = cellfun(@(x) any(ismember(reward_frames,x)), ...
        move_frames_split);
     % Sanity check: match between reward frames and cr movements
     if length(reward_frames) ~= sum(cr_movements)
         warning on;
         warning(['Animal ' num2str(curr_animal) ' session ' num2str(curr_session) ...
             ': ' num2str(length(reward_frames)) ' rewards, ' ...
             num2str(sum(cr_movements)) ' cr movements']);
     end
    
    move_durations(~cr_movements) = [];
    move_frames_split(~cr_movements) = [];
    
    % Get binary activity during movements
    n_rois = size(data(curr_animal).im(curr_session).roi_trace_thresh,1);
    curr_move_act = nan(n_rois,length(move_frames_split));
    for i = 1:length(move_frames_split)
        curr_move_act(:,i) = ...
            any(data(curr_animal).im(curr_session).roi_trace_thresh(:,move_frames_split{i}),2);
    end
    
    n_rois_move_act = nanmean(curr_move_act,1)';
    
    move_duration_bin_edges = [28:28:5*28];
    [~,bins] = histc(move_durations,move_duration_bin_edges);
    bins_used = unique(bins);
    good_bins = bins_used > 0 & bins_used < length(move_duration_bin_edges);
    n_rois_move_act_binned = grpstats(n_rois_move_act,bins);
    
    act_movetime_binned(curr_session,bins_used(good_bins),curr_animal) = ...
        n_rois_move_act_binned(good_bins);
    act_per_movetime(curr_animal,curr_session) = nanmean(n_rois_move_act./move_durations);
    act_per_move(curr_animal,curr_session) = nanmedian(n_rois_move_act);
    movetime(curr_animal,curr_session) = nanmedian(move_durations);
    
    end
    
end

figure;

subplot(2,3,1);
errorbar(nanmean(movetime/28),nanstd(movetime/28)./ ...
    sqrt(sum(~isnan(movetime))),'k','linewidth',2)
ylabel('Median movement time');
xlabel('Day');
    
subplot(2,3,2);
errorbar(nanmean(act_per_move/28),nanstd(act_per_move/28)./ ...
    sqrt(sum(~isnan(act_per_move))),'k','linewidth',2)
ylabel('Fraction active ROIs per movement');
xlabel('Day');

subplot(2,3,3);
errorbar(nanmean(act_per_movetime/28),nanstd(act_per_movetime/28)./ ...
    sqrt(sum(~isnan(act_per_movetime))),'k','linewidth',2)
ylabel('Fraction active ROIs per movement frame');
xlabel('Day');

subplot(2,3,4);
imagesc(nanmean(act_movetime_binned,3));colormap(hot);
xlabel('Move duration bin');
ylabel('Session');
title('Fraction active ROIs per movement');

subplot(2,3,5);
imagesc(nanmedian(act_movetime_binned,3));colormap(hot);
xlabel('Move duration bin');
ylabel('Session');
title('Fraction active ROIs per movement (median)');

subplot(2,3,6);
act_movetime_binned_norm = bsxfun(@times,act_movetime_binned, ...
    1./act_movetime_binned(:,1));
imagesc(nanmean(act_movetime_binned_norm,3));colormap(hot);
xlabel('Move duration bin');
ylabel('Session');
title('Fraction active ROIs per movement (bin 1 normalized)');



% Plot frac activity per movement against frac M cells
m = cell2mat(cellfun(@(x) nanmean(x,1), {classified_rois(:).movement}','uni',false));
q = cell2mat(cellfun(@(x) nanmean(x,1), {classified_rois(:).quiescent}','uni',false));
mqr = m./q;
col = lines(length(data));

figure; 

subplot(1,2,1);hold on;
for curr_animal = 1:length(data);
   plot(m(curr_animal,:),act_per_movetime(curr_animal,:), ...
       'o','color',col(curr_animal,:), ...
       'MarkerFaceColor',col(curr_animal,:));     
end
legend({'1' '2' '3' '4' '5' '6' '7' '8'})
xlabel('Frac movement ROIs')
ylabel('Frac active ROIs per movement time')

subplot(1,2,2);hold on;
for curr_animal = 1:length(data);
   plot(mqr(curr_animal,:),act_per_movetime(curr_animal,:), ...
       'o','color',col(curr_animal,:), ...
       'MarkerFaceColor',col(curr_animal,:));     
end
legend({'1' '2' '3' '4' '5' '6' '7' '8'})
xlabel('Frac movement ROIs /quiescent ROIs')
ylabel('Frac active ROIs per movement time')


%% Average activity correlation of shared mROIs

avg_onset_act = cell(8,1);
avg_offset_act = cell(8,1);
for curr_animal = 1:length(data);
    for curr_session = 1:14
        
        avg_onset_act{curr_animal}(:,:,curr_session) = ...
            permute(nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned,1),[3 2 1]);
        
        avg_offset_act{curr_animal}(:,:,curr_session) = ...
            permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned,1),[3 2 1]);
        
    end
end

avg_overlap_m_act_corr = nan(14,14,length(data));
avg_combine_m_act_corr = nan(14,14,length(data));
for curr_animal = 1:length(data);
    for curr_day1 = 1:14;
        for curr_day2 = 1:14;
            
            curr_m1 = classified_rois(curr_animal).movement(:,curr_day1);
            curr_m2 = classified_rois(curr_animal).movement(:,curr_day2);
            
            overlap_m = curr_m1 & curr_m2;
            combine_m = curr_m1 | curr_m2;
            
            use_frames = analysis(curr_animal).surrounding_frames > 0 & ...
                analysis(curr_animal).surrounding_frames < 28*1.5;
            
            curr_act_1 = reshape(zscore(avg_onset_act{curr_animal}(overlap_m,use_frames,curr_day1),[],2)',[],1);
            curr_act_2 = reshape(zscore(avg_onset_act{curr_animal}(overlap_m,use_frames,curr_day2),[],2)',[],1);            
            curr_corr = corrcoef(curr_act_1,curr_act_2);           
            avg_overlap_m_act_corr(curr_day1,curr_day2,curr_animal) = curr_corr(2);
            
            curr_act_1 = reshape(zscore(avg_onset_act{curr_animal}(combine_m,use_frames,curr_day1),[],2)',[],1);
            curr_act_2 = reshape(zscore(avg_onset_act{curr_animal}(combine_m,use_frames,curr_day2),[],2)',[],1);            
            curr_corr = corrcoef(curr_act_1,curr_act_2);           
            avg_combine_m_act_corr(curr_day1,curr_day2,curr_animal) = curr_corr(2);
            
        end
    end
end

avg_overlap_m_onoff_act_corr = nan(14,14,length(data));
avg_combine_m_onoff_act_corr = nan(14,14,length(data));
for curr_animal = 1:length(data);
    for curr_day1 = 1:14;
        for curr_day2 = 1:14;
            
            curr_m1 = classified_rois(curr_animal).movement(:,curr_day1);
            curr_m2 = classified_rois(curr_animal).movement(:,curr_day2);
            
            overlap_m = curr_m1 & curr_m2;
            combine_m = curr_m1 | curr_m2;
            
            use_onset_frames = analysis(curr_animal).surrounding_frames > 0 & ...
                analysis(curr_animal).surrounding_frames < 28*1.5;
            
            use_offset_frames = analysis(curr_animal).surrounding_frames < 0 & ...
                analysis(curr_animal).surrounding_frames > -28*1.5;
            
            curr_act_1 = reshape(zscore([avg_onset_act{curr_animal}(overlap_m,use_frames,curr_day1) ...
                avg_offset_act{curr_animal}(overlap_m,use_frames,curr_day1)],[],2)',[],1);
            curr_act_2 = reshape(zscore([avg_onset_act{curr_animal}(overlap_m,use_frames,curr_day2) ...
                avg_offset_act{curr_animal}(overlap_m,use_frames,curr_day2)],[],2)',[],1);           
            curr_corr = corrcoef(curr_act_1,curr_act_2);           
            avg_overlap_m_onoff_act_corr(curr_day1,curr_day2,curr_animal) = curr_corr(2);
            
            curr_act_1 = reshape(zscore([avg_onset_act{curr_animal}(combine_m,use_frames,curr_day1) ...
                avg_offset_act{curr_animal}(combine_m,use_frames,curr_day1)],[],2)',[],1);
            curr_act_2 = reshape(zscore([avg_onset_act{curr_animal}(combine_m,use_frames,curr_day2) ...
                avg_offset_act{curr_animal}(combine_m,use_frames,curr_day2)],[],2)',[],1);           
            curr_corr = corrcoef(curr_act_1,curr_act_2);           
            avg_combine_m_onoff_act_corr(curr_day1,curr_day2,curr_animal) = curr_corr(2);
            
        end
    end
end

figure;

subplot(2,2,1);
imagesc(nanmean(avg_overlap_m_act_corr,3)); colormap(hot);
caxis([0 max(AP_itril(nanmean(avg_overlap_m_act_corr,3),-1))]);
xlabel('Day')
ylabel('Day');
title('Correlation of overlap movement ROI z-scored mean activity')
colorbar;

subplot(2,2,2);
imagesc(nanmean(avg_overlap_m_onoff_act_corr,3)); colormap(hot);
caxis([0 max(AP_itril(nanmean(avg_overlap_m_onoff_act_corr,3),-1))]);
xlabel('Day')
ylabel('Day');
title('Correlation of overlap movement ROI z-scored mean on/off activity')
colorbar;

subplot(2,2,3);
imagesc(nanmean(avg_combine_m_act_corr,3)); colormap(hot);
caxis([0 max(AP_itril(nanmean(avg_overlap_m_act_corr,3),-1))]);
xlabel('Day')
ylabel('Day');
title('Correlation of combined movement ROI z-scored mean activity')
colorbar;

subplot(2,2,4);
imagesc(nanmean(avg_combine_m_onoff_act_corr,3)); colormap(hot);
caxis([0 max(AP_itril(nanmean(avg_overlap_m_onoff_act_corr,3),-1))]);
xlabel('Day')
ylabel('Day');
title('Correlation of combined movement ROI z-scored mean on/off activity')
colorbar;


%% Plot overlap M activity correlation with movement correlation
% (requires last cell and behavior analysis)


col = lines(8);
figure; hold on;
for curr_animal = 1:length(data);
    
    curr_move = AP_itril(movement_corr_grid(:,:,curr_animal),-1);
    curr_act = AP_itril(avg_overlap_m_onoff_act_corr(:,:,curr_animal),-1);
    curr_act = (curr_act - nanmean(curr_act))/nanstd(curr_act);
    
    plot(zscore(curr_move),curr_act,'.','color',col(curr_animal,:));
    
end


% bin and average
n_bins = 10;
zscore_bin = linspace(-2,2,n_bins+1);
bin_centers = zscore_bin(1:end-1) + diff(zscore_bin)/2;
move_act_bins = nan(8,n_bins);
for curr_animal = 1:length(data)
    
    curr_move = zscore(AP_itril(movement_corr_grid(:,:,curr_animal),-1));
    curr_act = AP_itril(avg_overlap_m_onoff_act_corr(:,:,curr_animal),-1);
    curr_act = (curr_act - nanmean(curr_act))/nanstd(curr_act);
    
    [~,bins] = histc(curr_move,zscore_bin);
    bins_used = unique(bins);
    good_bins = bins_used > 0 & bins_used <= n_bins;
    curr_act_binned = grpstats(curr_act,bins,'nanmean');
    
    move_act_bins(curr_animal,bins_used(good_bins)) = curr_act_binned(good_bins);
end

errorbar(bin_centers,nanmean(move_act_bins),nanstd(move_act_bins)./ ...
    sqrt(sum(~isnan(move_act_bins))),'k','linewidth',2);

xlabel('z-scored movement correlation');
ylabel('Overlap mROI z-scored average act corr');



%% Example m2q ROI average activity plots

curr_animal = 3;

curr_m2q = ...
    any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
    classified_rois(curr_animal).quiescent,2);

m_q = classified_rois(curr_animal).movement(:,:) - ...
    classified_rois(curr_animal).quiescent(:,:);

% These ROIs were selected as a handful that are good-looking all days
curr_raw_rois = [103,118,159,165,181,260,307,297,310,362,387];

curr_rois = find(cellfun(@(x) any(ismember(curr_raw_rois,x)),data(curr_animal).roi_group));

figure;
imagesc(m_q(curr_rois,:));
colormap(gray);

figure;
col = jet(14);
for curr_roi_idx = 1:length(curr_rois);
      
    curr_roi = curr_rois(curr_roi_idx);
    
    curr_onset = cellfun(@(x) x(:,:,curr_roi),{analysis(curr_animal).im(:).move_onset_aligned},'uni',false);
    curr_offset = cellfun(@(x) x(:,:,curr_roi),{analysis(curr_animal).im(:).move_offset_aligned},'uni',false);
    
    curr_onset_mean = cell2mat(cellfun(@nanmean,curr_onset,'uni',false)');
    curr_offset_mean = cell2mat(cellfun(@nanmean,curr_offset,'uni',false)');
    
    curr_onoff = [curr_onset_mean(:,90-30:90+10),curr_offset_mean(:,90-10:90+30)];
    
    subplot(4,4,curr_roi_idx);
    imagesc(curr_onset_mean);colormap(hot);
    
    % To plot mean as lines
%     subplot(4,4,curr_roi_idx); hold on;   
%     for i = 1:14
%         plot(curr_onset_mean(i,:),'color',col(i,:));
%     end
    
end

% Plot whole traces of m2q ROI
figure; hold on;
for i = 1:14
    plot(data(curr_animal).im(i).roi_trace_df(236,:)+i*5,'k');
end

%% Activity spacing of m/q ROIs

% Get inter-event-interval fano factor
ff = cell(length(data),1);
for curr_animal = 1:length(data)
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    ff{curr_animal} = nan(n_rois,14);
    for curr_session = 1:14;        
        for curr_roi = 1:n_rois  
            
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session). ...
                roi_trace_thresh(curr_roi,:)>0]) == 1);
                        
            curr_iei = diff(curr_act_starts);
            
            ff{curr_animal}(curr_roi,curr_session) = std(curr_iei)./mean(curr_iei); 
            
        end        
    end
end

% Split IEI ff into m/q
ff_mq = cell(length(data),1);
for curr_animal = 1:length(data)
    
    curr_m = classified_rois(curr_animal).movement;
    curr_q = classified_rois(curr_animal).quiescent;
    curr_ff = ff{curr_animal};
    
    curr_ff_m = curr_ff;
    curr_ff_m(~curr_m) = NaN;
    
    curr_ff_q = curr_ff;
    curr_ff_q(~curr_q) = NaN;
    
    ff_mq{curr_animal} = [nanmean(curr_ff_m,2),nanmean(curr_ff_q,2)];
end

% Get ff change of m2q ROIs
ff_m2q = cell(length(data),1);
for curr_animal = 1:length(data)
    
    curr_m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    
    curr_m = classified_rois(curr_animal).movement(curr_m2q,:);
    curr_q = classified_rois(curr_animal).quiescent(curr_m2q,:);
    curr_ff = ff{curr_animal}(curr_m2q,:);
    
    curr_ff_m = curr_ff;
    curr_ff_m(~curr_m) = NaN;
    
    curr_ff_q = curr_ff;
    curr_ff_q(~curr_q) = NaN;
    
    ff_m2q{curr_animal}= [nanmean(curr_ff_m,2),nanmean(curr_ff_q,2);];
end

figure;hold on
ff_mq_mean = cell2mat(cellfun(@nanmean,ff_mq,'uni',false));
ff_m2q_mean = cell2mat(cellfun(@nanmean,ff_m2q,'uni',false));
errorbar(nanmean(ff_mq_mean),nanstd(ff_mq_mean)./sqrt(sum(~isnan(ff_mq_mean))),'k','linewidth',2);
errorbar(nanmean(ff_m2q_mean),nanstd(ff_m2q_mean)./sqrt(sum(~isnan(ff_m2q_mean))),'r','linewidth',2);

ylabel('IEI fano factor')
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Movement','Quiescent'});
legend({'All ROIs','M2Q ROIs'})


%% Timing of m2q ROI activity

n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
% Movement-binarized lever
back_time = 29*30;
forward_time = 29*30;
act_aligned_move_binary = nan(n_rois,back_time+forward_time+1,14);

curr_class = classified_rois(curr_animal).movement | ...
        classified_rois(curr_animal).quiescent;

for curr_session = 1:14;  
    
    for curr_roi = find(curr_class(:,curr_session))'
        curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_roi,:)>0]) == 1);
        curr_act_starts(curr_act_starts-back_time <= 0 | ...
            curr_act_starts + forward_time > length(analysis(curr_animal).lever(curr_session).lever_move_frames)) = [];
        
        curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_starts(x)-back_time: ...
            curr_act_starts(x)+forward_time,(1:length(curr_act_starts))','uni',false));
        
        % Movement-binarzed lever
        curr_act_move = analysis(curr_animal).lever(curr_session).lever_move_frames(curr_act_move_idx);
        % Movement-binarized onset/offset
        %lever_diff = smooth(diff([0;analysis(curr_animal).lever(curr_session).lever_move_frames]),10);
        %curr_act_move = lever_diff(curr_act_move_idx);
        
        act_aligned_move_binary(curr_roi,:,curr_session) = nanmean(curr_act_move);
    end
    disp(curr_session);
end

% Average surrounding binary movement on m/q sessions
mq_aligned_move_binary = nan(n_rois,back_time+forward_time+1,2);

for curr_roi = find(any(curr_class,2))'
    
    mq_aligned_move_binary(curr_roi,:,1) = ...
        nanmean(act_aligned_move_binary(curr_roi,:, ...
        classified_rois(curr_animal).movement(curr_roi,:)),3);

    mq_aligned_move_binary(curr_roi,:,2) = ...
        nanmean(act_aligned_move_binary(curr_roi,:, ...
        classified_rois(curr_animal).quiescent(curr_roi,:)),3);
end

curr_m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    
m2q_aligned_move_binary = mq_aligned_move_binary(curr_m2q,:,:);
m2m_aligned_move_binary = mq_aligned_move_binary(any(classified_rois(curr_animal).movement,2) & ~ curr_m2q,:,:);
q2q_aligned_move_binary = mq_aligned_move_binary(any(classified_rois(curr_animal).quiescent,2) & ~ curr_m2q,:,:);

figure;hold on
plot([-back_time:forward_time]/28,nanmean(m2q_aligned_move_binary(:,:,1)),'k')
plot([-back_time:forward_time]/28,nanmean(m2q_aligned_move_binary(:,:,2)),'r')
plot([-back_time:forward_time]/28,nanmean(m2m_aligned_move_binary(:,:,1)),'--k')
plot([-back_time:forward_time]/28,nanmean(q2q_aligned_move_binary(:,:,2)),'--r')
legend({'m2q/m' 'm2q/q' 'm2m/m' 'q2q/q'})
ylabel('Movement probability')
xlabel('Time from activity onset')
line([0 0],ylim,'color','b')


%% overlap mROI Trial-by-trial movement-activiy correlation between/across sessions
% (this takes ~ 30 mins to run)

%stages = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]}; 
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
    
    % Define trials to use
    n_sessions = length(analysis(curr_animal).im);
    use_sessions = arrayfun(@(x) ~isempty(analysis(curr_animal).lever(x).rewarded_movement),1:n_sessions);
    
    % Trial-by-trial movement correlation
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time),movement(trials),'uni',false)'), ...
        {analysis(curr_animal).lever(use_sessions).rewarded_movement}, ...
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
    use_cells = classified_rois(curr_animal).movement;
    
    use_activity = cell(n_sessions,1);
    for curr_session = 1:n_sessions
        curr_act = analysis(curr_animal).im(curr_session).move_onset_aligned(:,use_frames,:);
        curr_act(:,:,~use_cells(:,curr_session)) = NaN;
        use_activity{curr_session} = reshape(permute(curr_act,[2,3,1]),[],size(curr_act,1));
    end
    
    activity_corr = cell(n_sessions);
    activity_corr(use_sessions,use_sessions) = ...
        mat2cell(corrcoef(horzcat(use_activity{:}),'rows','pairwise'),n_session_movements,n_session_movements);
    activity_corr(find(eye(size(activity_corr)))) = ...
        cellfun(@(x) AP_itril(x,-1),activity_corr(find(eye(size(activity_corr)))),'uni',false);
    activity_corr(find(~eye(size(activity_corr)))) = ...
        cellfun(@(x) x(:),activity_corr(find(~eye(size(activity_corr)))),'uni',false);
    
    % Bin trials by movement correlation
    corr_bin = linspace(-1,1,20);
    corr_bin_use = [corr_bin(1:end-1) Inf];
    corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];
    
    % Groups stages of learning, bin  
    use_stages = find(cellfun(@(x) any(intersect(x,find(use_sessions))),stages));
    
    [~,stages_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_corr,'uni',false);
        
    curr_stage_mean = cell(length(stages));
    curr_stage_sem = cell(length(stages));
    for curr_stage_1 = use_stages;
        for curr_stage_2 = use_stages
                  
            curr_mean = grpstats(vertcat(activity_corr{curr_stage_1,curr_stage_2}), ...
                vertcat(stages_bin_idx{curr_stage_1,curr_stage_2}),'mean');
            curr_sem = grpstats(vertcat(activity_corr{curr_stage_1,curr_stage_2}), ...
                vertcat(stages_bin_idx{curr_stage_1,curr_stage_2}),'sem');
            
            bins_used = unique(vertcat(stages_bin_idx{curr_stage_1,curr_stage_2}));
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

% Plot act/move corr relative to first , middle, last session
figure;
col = jet(length(stages));

subplot(3,1,1); hold on
for i = 1:length(stages)
    plot(corr_bin_plot,activity_corr_stages_mean_animals{1,i},'color',col(i,:),'linewidth',2)
end
ylabel('Trial-by-trial activity correlation')
xlabel('Trial-by-trial movement correlation')
title('Compared: day 1')

subplot(3,1,2); hold on
for i = 1:length(stages)
    plot(corr_bin_plot,activity_corr_stages_mean_animals{7,i},'color',col(i,:),'linewidth',2)
end
ylabel('Trial-by-trial activity correlation')
xlabel('Trial-by-trial movement correlation')
title('Compared: day 7')

subplot(3,1,3); hold on
for i = 1:length(stages)
    plot(corr_bin_plot,activity_corr_stages_mean_animals{end,i},'color',col(i,:),'linewidth',2)
end
ylabel('Trial-by-trial activity correlation')
xlabel('Trial-by-trial movement correlation')
title('Compared: day 14')



% Estimate slope of positive correlations
use_bins = 5:18;
corr_fit = cellfun(@(x) polyfit(1:length(use_bins),x(use_bins),1), ...
    activity_corr_stages_mean_animals,'uni',false);
corr_slope = cellfun(@(x) x(1),corr_fit);

figure;
imagesc(corr_slope);colormap(hot)
xlabel('Day');
ylabel('Day');
title('Activity/movement correlation slope (overlap mROI only)')


% Plot just the average of the last few bins
max_bins = min(min(cellfun(@(x) find(~isnan(x),1,'last'),activity_corr_stages_mean_animals)));
actmove_corr_mean = cellfun(@(x) nanmean(x(16:max_bins)),activity_corr_stages_mean_animals);
figure;
imagesc(actmove_corr_mean);colormap(hot)
ylabel('Day');
xlabel('Day');
title('Average activity corr in high movement corr bins');

%% Plot principle components in time

curr_animal = 3;
curr_day = 12;

curr_data = data(curr_animal).im(curr_day).roi_trace_df;
for i = 1:size(curr_data,1);
   curr_data(i,:) = smooth(curr_data(i,:),100,'loess');    
end

%[coeff,score,latent] = princomp(zscore(curr_data,[],2)');

% or average or m/q rois:
score = zeros(size(curr_data,2),3);
score(:,1) = nanmean(zscore(curr_data(classified_rois(curr_animal).movement(:,curr_day),:),[],2),1);
score(:,2) = nanmean(zscore(curr_data(classified_rois(curr_animal).quiescent(:,curr_day),:),[],2),1);

curr_move = analysis(curr_animal).lever(curr_day).lever_move_frames;

figure;
axis;
axis vis3d

plot_back = 100;
curr_plot = nan(plot_back,3);

curr_line = plot3(curr_plot(:,1),curr_plot(:,2),curr_plot(:,3),'k');
hold on;

curr_dot = plot3(curr_plot(1,1),curr_plot(1,2),curr_plot(1,3),'ok');

xlim([min(score(:,1)),max(score(:,1))]);
ylim([min(score(:,2)),max(score(:,2))]);
%zlim([min(score(:,3)),max(score(:,3))]);

for curr_frame = 1:length(curr_move)
    
    curr_plot = circshift(curr_plot,[1,0]);
    curr_plot(1,:) = score(curr_frame,1:3);
    
    set(curr_line,'XData',curr_plot(:,1));
    set(curr_line,'YData',curr_plot(:,2));
    set(curr_line,'ZData',curr_plot(:,3));
    
    set(curr_dot,'XData',curr_plot(1,1));
    set(curr_dot,'YData',curr_plot(1,2));
    set(curr_dot,'ZData',curr_plot(1,3));
    
    if curr_move(curr_frame);
        set(curr_line,'color','r');
    else
        set(curr_line,'color','k');
    end
    
    drawnow;
    pause(0.001);
    
end



%% Correlation of activity-triggered movement for mROIs

act_aligned_move_corr_all = nan(14,14,100,length(data));
for curr_animal = 1:length(data)
    
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    % Raw lever trace
    back_time = 3000;
    forward_time = 3000;
    act_aligned_move = nan(n_rois,back_time+forward_time+1,14);
    
    for curr_session = 1:14;
        % downsample lever
        downsample_factor = 100;
        lever_force_resample = resample(data(curr_animal).bhv(curr_session).lever_force,1,downsample_factor);
        for curr_cell = find(classified_rois(curr_animal).movement(:,curr_session))'
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
            curr_act_start_time = round(data(curr_animal).bhv(curr_session).frame_times(curr_act_starts)*(10000/downsample_factor));
            curr_act_start_time(curr_act_start_time-back_time <= 0 | ...
                curr_act_start_time + forward_time > length(lever_force_resample)) = [];
            
            curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_start_time(x)-back_time: ...
                curr_act_start_time(x)+forward_time,(1:length(curr_act_start_time))','uni',false));
            
            % Lever trace
            curr_act_move = lever_force_resample(curr_act_move_idx);
            
            act_aligned_move(curr_cell,:,curr_session) = nanmean(curr_act_move);
        end
    end
        
    %%%%%%%%% CAN BUILD THIS IN LATER
    % TEST STATISTIC: shuffle non-nan ROIs within session
    for curr_session = 1:14
       curr_rois = find(any(~isnan(act_aligned_move(:,:,curr_session)),2)); 
       act_aligned_move(shake(curr_rois),:,curr_session) = ...
           act_aligned_move(curr_rois,:,curr_session);
    end
    %%%%%%%%%
    
    act_aligned_move_corr_movetimes = nan(14,14,100);
    for i = 1:100;
        surround_time = i*10;
        use_movetime = 3000-surround_time:3000+surround_time;
        act_aligned_move_corr = nan(14,14,n_rois);
        for curr_roi = 1:n_rois
            act_aligned_move_corr(:,:,curr_roi) = ...
                corrcoef(permute(act_aligned_move(curr_roi,use_movetime,:),[2,3,1]));
        end
        act_aligned_move_corr_movetimes(:,:,i) = nanmean(act_aligned_move_corr,3);
    end
    
    act_aligned_move_corr_all(:,:,:,curr_animal) = act_aligned_move_corr_movetimes;
    
    disp(curr_animal);
end

act_aligned_move_corr_mean = nanmean(act_aligned_move_corr_all,4);

figure;
imagesc(act_aligned_move_corr_mean(:,:,20));colormap(hot)
ylabel('Day');
xlabel('Day');
title('Activity-triggered movement correlation')



%% Correlation of activity-triggered movement for mROIs (binary movement)

act_aligned_move_corr_all = nan(14,14,29,length(data));
for curr_animal = 1:length(data)
    
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    % Movement-binarized lever
    back_time = 29*30;
    forward_time = 29*30;
    act_aligned_move_binary = nan(n_rois,back_time+forward_time+1,14);
    
    for curr_session = 1:14;
        for curr_cell = find(classified_rois(curr_animal).movement(:,curr_session))'
            curr_act_starts = find(diff([0 data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:)>0]) == 1);
            curr_act_starts(curr_act_starts-back_time <= 0 | ...
                curr_act_starts + forward_time > length(analysis(curr_animal).lever(curr_session).lever_move_frames)) = [];
            
            curr_act_move_idx = cell2mat(arrayfun(@(x) curr_act_starts(x)-back_time: ...
                curr_act_starts(x)+forward_time,(1:length(curr_act_starts))','uni',false));
            
            % Movement-binarzed lever
            curr_act_move = analysis(curr_animal).lever(curr_session).lever_move_frames(curr_act_move_idx);
            % Movement-binarized onset/offset
            %lever_diff = smooth(diff([0;analysis(curr_animal).lever(curr_session).lever_move_frames]),10);
            %curr_act_move = lever_diff(curr_act_move_idx);
            
            act_aligned_move_binary(curr_cell,:,curr_session) = nanmean(curr_act_move);
        end
    end
    
    
    act_aligned_move_corr_movetimes = nan(14,14,29);
    for i = 1:29;
        surround_time = i*1;
        use_movetime = 29*30+1-surround_time:29*30+1+surround_time;
        act_aligned_move_corr = nan(14,14,n_rois);
        for curr_roi = 1:n_rois
            act_aligned_move_corr(:,:,curr_roi) = ...
                corrcoef(permute(act_aligned_move_binary(curr_roi,use_movetime,:),[2,3,1]));
        end
        act_aligned_move_corr_movetimes(:,:,i) = nanmean(act_aligned_move_corr,3);
    end
    
    act_aligned_move_corr_all(:,:,:,curr_animal) = act_aligned_move_corr_movetimes;
    
    disp(curr_animal);
end

act_aligned_move_corr_mean = nanmean(act_aligned_move_corr_all,4);

figure;
imagesc(act_aligned_move_corr_mean(:,:,14));colormap(hot)
ylabel('Day');
xlabel('Day');
title('Activity-triggered binary movement correlation')


%% Cumulative correlation ROIs with movement

% this didn't show anything, didn't bother saving

figure;
for curr_animal = 1:length(data)
    
    n_rois = size(data(curr_animal).im(curr_day).roi_trace_df,1);
    curr_actmovecorr = nan(n_rois,14);
    
    for curr_day = 1:14
        curr_actmovecorr(:,curr_day) = 1 - pdist2(data(curr_animal).im(curr_day).roi_trace_thresh, ...
            analysis(curr_animal).lever(curr_day).lever_move_frames','correlation');
    end
    
    col = jet(14);
    subplot(2,4,curr_animal); hold on;
    for curr_day = 1:14
        plot(sort(curr_actmovecorr(:,curr_day)),(1:n_rois)/n_rois,'color',col(curr_day,:));
    end
    
end



%% Average m/q and m2q vs q2q activity over time

% session x frames x animal x [m,q,q/m2q,q/non-m2q]
avg_onset_act = nan(14,size(analysis(1).im(1).move_onset_aligned,2),length(data),3);
avg_offset_act = nan(14,size(analysis(1).im(1).move_onset_aligned,2),length(data),3);
for curr_animal = 1:length(data);
    
    m2q = ...
    any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
    classified_rois(curr_animal).quiescent,2);

    for curr_session = 1:length(data(curr_animal).im)
        
        m = classified_rois(curr_animal).movement(:,curr_session);
        q = classified_rois(curr_animal).quiescent(:,curr_session);
        q_m2q = q & m2q;
        q_nonm2q = q & ~m2q;
        
        avg_onset_act(curr_session,:,curr_animal,1) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,m),1),[3 2 1]));
        avg_onset_act(curr_session,:,curr_animal,2) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,q),1),[3 2 1]));
        avg_onset_act(curr_session,:,curr_animal,3) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,q_m2q),1),[3 2 1]));
        avg_onset_act(curr_session,:,curr_animal,4) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,q_nonm2q),1),[3 2 1]));
        
        avg_offset_act(curr_session,:,curr_animal,1) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,m),1),[3 2 1]));
        avg_offset_act(curr_session,:,curr_animal,2) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,q),1),[3 2 1]));
        avg_offset_act(curr_session,:,curr_animal,3) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,q_m2q),1),[3 2 1]));
        avg_offset_act(curr_session,:,curr_animal,4) = ...
            nanmean(permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,q_nonm2q),1),[3 2 1]));
        
    end
end

col = jet(14);
figure; 
subplot(2,4,1); hold on;
for i = 1:14
   plot(nanmean(avg_onset_act(i,:,:,1),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Movement ROI - onset')
subplot(2,4,3); hold on;
for i = 1:14
   plot(nanmean(avg_onset_act(i,:,:,2),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Quiescent ROI - onset')
subplot(2,4,5); hold on;
for i = 1:14
   plot(nanmean(avg_onset_act(i,:,:,3),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('M2Q ROI - onset')
subplot(2,4,7); hold on;
for i = 1:14
   plot(nanmean(avg_onset_act(i,:,:,4),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Q2Q ROI - onset')
subplot(2,4,2); hold on;
for i = 1:14
   plot(nanmean(avg_offset_act(i,:,:,1),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Movement ROI - offset')
subplot(2,4,4); hold on;
for i = 1:14
   plot(nanmean(avg_offset_act(i,:,:,2),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Quiescent ROI - offset')
subplot(2,4,6); hold on;
for i = 1:14
   plot(nanmean(avg_offset_act(i,:,:,3),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('M2Q ROI - offset')
subplot(2,4,8); hold on;
for i = 1:14
   plot(nanmean(avg_offset_act(i,:,:,4),3),'color',col(i,:));
end
line([90,90],ylim,'color','k','linestyle','--');
title('Q2Q ROI - offset')

% Compare later stage M2Q vs Q2Q ROIs
use_sessions = 10:14;
m2q_avg_on = nanmean(nanmean(avg_onset_act(use_sessions,:,:,3),3));
q2q_avg_on = nanmean(nanmean(avg_onset_act(use_sessions,:,:,4),3));

m2q_avg_off = nanmean(nanmean(avg_offset_act(use_sessions,:,:,3),3));
q2q_avg_off = nanmean(nanmean(avg_offset_act(use_sessions,:,:,4),3));

m2q_avg_onoff_norm = [m2q_avg_on(1:90+30) m2q_avg_off(90-30:end)]./nanmean(m2q_avg_on(1:30));
q2q_avg_onoff_norm = [q2q_avg_on(1:90+30) q2q_avg_off(90-30:end)]./nanmean(q2q_avg_on(1:30));

figure;
subplot(1,2,1);hold on
plot(q2q_avg_on(1:90+30)./nanmean(q2q_avg_on(1:30)),'k','linewidth',2)
plot(m2q_avg_on(1:90+30)./nanmean(m2q_avg_on(1:30)),'r','linewidth',2)
line([90,90],ylim,'linestyle','--');
title('Movement onset')
subplot(1,2,2);hold on
plot(q2q_avg_off(90-30:end)./nanmean(q2q_avg_on(1:30)),'k','linewidth',2)
plot(m2q_avg_off(90-30:end)./nanmean(m2q_avg_on(1:30)),'r','linewidth',2)
line([30,30],ylim,'linestyle','--');
title('Movement offset')
legend({'Q2Q','M2Q'});



% heatmap of m2q vs q2q ROIs
avg_offset_act_q_rois = cell(14,length(data),2);
for curr_animal = 1:length(data);
    
    m2q = ...
    any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
    classified_rois(curr_animal).quiescent,2);

    for curr_session = 1:14
        
        q = classified_rois(curr_animal).quiescent(:,curr_session);
        q_m2q = q & m2q;
        q_nonm2q = q & ~m2q;
        
        avg_offset_act_q_rois{curr_session,curr_animal,1} = ...
            permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,q_m2q),1),[3 2 1]);
        
        avg_offset_act_q_rois{curr_session,curr_animal,2} = ...
            permute(nanmean(analysis(curr_animal).im(curr_session).move_offset_aligned(:,:,q_nonm2q),1),[3 2 1]);

    end
end






%% SfN poster: additional figures

%%% Example traces

curr_animal = 3;
curr_session = 11;
% m/q/u examples
rois = [39,232,40]; % other optional q = 270
t = 44000:49000;

% Lever / movement
curr_move = analysis(curr_animal).lever(curr_session).lever_move_frames(t);
curr_lever = data(curr_animal).bhv(curr_session).lever_force_plot(t);
curr_lever_move = curr_lever;
curr_lever_move(~curr_move) = NaN;
curr_lever_still = curr_lever;
curr_lever_still(curr_move) = NaN;

% Cue/reward frames
cue_frames = cellfun(@(x) x.states.cue(1),data(curr_animal).bhv(curr_session).bhv_frames);
rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),data(curr_animal).bhv(curr_session).bhv_frames);
reward_frames = cellfun(@(x) x.states.reward(1),data(curr_animal).bhv(curr_session).bhv_frames(rewarded_trials));

cue_frames(~ismember(cue_frames,t)) = [];
reward_frames(~ismember(reward_frames,t)) = [];

cue_frames = cue_frames - t(1) + 1;
reward_frames = reward_frames - t(1) + 1;

% Plot everything
figure; hold on;

plot(curr_move-2,'k');
plot(curr_lever_still,'k');
plot(curr_lever_move,'r');
rois = [39 232 88];
for i = 1:length(rois)
    curr_roi = rois(i);
    plot(data(curr_animal).im(curr_session).roi_trace_df(curr_roi,t) + i*5,'k');
end
for i = 1:length(cue_frames)
    line([cue_frames(i),cue_frames(i)],ylim,'color','g');
end
for i = 1:length(reward_frames)
    line([reward_frames(i),reward_frames(i)],ylim,'color','b');
end


figure;
curr_on_aligned = analysis(curr_animal).im(curr_session).move_onset_aligned(:,1:90+30,rois(1));
curr_off_aligned = analysis(curr_animal).im(curr_session).move_offset_aligned(:,90-30:end,rois(1));
subplot(2,6,1);
imagesc(curr_on_aligned);colormap(hot);
subplot(2,6,2);
imagesc(curr_off_aligned);colormap(hot);
subplot(2,6,7);
plot(nanmean(curr_on_aligned),'k');
subplot(2,6,8);
plot(nanmean(curr_off_aligned),'k');

curr_on_aligned = analysis(curr_animal).im(curr_session).move_onset_aligned(:,1:90+30,rois(2));
curr_off_aligned = analysis(curr_animal).im(curr_session).move_offset_aligned(:,90-30:end,rois(2));
subplot(2,6,3);
imagesc(curr_on_aligned);colormap(hot);
subplot(2,6,4);
imagesc(curr_off_aligned);colormap(hot);
subplot(2,6,9);
plot(nanmean(curr_on_aligned),'k');
subplot(2,6,10);
plot(nanmean(curr_off_aligned),'k');

curr_on_aligned = analysis(curr_animal).im(curr_session).move_onset_aligned(:,1:90+30,rois(3));
curr_off_aligned = analysis(curr_animal).im(curr_session).move_offset_aligned(:,90-30:end,rois(3));
subplot(2,6,5);
imagesc(curr_on_aligned);colormap(hot);
subplot(2,6,6);
imagesc(curr_off_aligned);colormap(hot);
subplot(2,6,11);
plot(nanmean(curr_on_aligned),'k');
subplot(2,6,12);
plot(nanmean(curr_off_aligned),'k');


% Timing heatmap of all mROIs (max normalized)
early_sessions = 1:3;
late_sessions = 12:14;

early_moveact = cell(length(data),1);
for curr_animal = 1:length(data);
    curr_act = {analysis(curr_animal).im(early_sessions).move_onset_aligned};
    for i = 1:length(curr_act)
        curr_act{i}(:,:,~classified_rois(curr_animal).movement(:,early_sessions(i))) = NaN;
    end
    curr_act = permute(nanmean(vertcat(curr_act{:}),1),[3,2,1]);
    curr_act = curr_act(any(classified_rois(curr_animal).movement(:,early_sessions),2),:);
    early_moveact{curr_animal} = curr_act;
end
early_moveact_cat = vertcat(early_moveact{:});
[~,max_idx] = max(early_moveact_cat(:,90:end),[],2);
[~,sort_idx] = sort(max_idx);
early_moveact_cat = early_moveact_cat(sort_idx,:);
early_moveact_cat = bsxfun(@times,early_moveact_cat,1./max(early_moveact_cat,[],2));


late_moveact = cell(length(data),1);
for curr_animal = 1:length(data);
    curr_act = {analysis(curr_animal).im(late_sessions).move_onset_aligned};
    for i = 1:length(curr_act)
        curr_act{i}(:,:,~classified_rois(curr_animal).movement(:,late_sessions(i))) = NaN;
    end
    curr_act = permute(nanmean(vertcat(curr_act{:}),1),[3,2,1]);
    curr_act = curr_act(any(classified_rois(curr_animal).movement(:,late_sessions),2),:);
    late_moveact{curr_animal} = curr_act;
end
late_moveact_cat = vertcat(late_moveact{:});
[~,max_idx] = max(late_moveact_cat(:,90:end),[],2);
[~,sort_idx] = sort(max_idx);
late_moveact_cat = late_moveact_cat(sort_idx,:);
late_moveact_cat = bsxfun(@times,late_moveact_cat,1./max(late_moveact_cat,[],2));

figure;
subplot(1,2,1);
imagesc(early_moveact_cat(:,90-28:end));
colormap(hot);
line([28,28],ylim,'color','w')
title('Days 1-3')
subplot(1,2,2);
imagesc(late_moveact_cat(:,90-28:end));
colormap(hot);
line([28,28],ylim,'color','w')
title('Days 12-14');





%% Visualize / draw maps of M/Q/U ROIs (this is for PC)

load('C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\AP142\AP142_batch_thresh_roi\150803_AP142_001_001_corrected_summed_50.roi','-MAT');
load('C:\Users\Andy\Desktop\roi_info.mat');

roi_map = zeros(512,512,3,14);
roi_masks = cellfun(@(x) poly2mask(x(:,1),x(:,2),512,512),polygon.ROI,'uni',false);
roi_masks = cat(3,roi_masks{:});

for curr_day = 1:14
    
    movement_rois = vertcat(grp{find(classified(:,curr_day) == 1)});
    quiescent_rois = vertcat(grp{find(classified(:,curr_day) == -1)});
    unclassified_rois = vertcat(grp{find(classified(:,curr_day) == 0)});
    
    curr_movement_mask = any(roi_masks(:,:,movement_rois),3);
    curr_quiescent_mask = any(roi_masks(:,:,quiescent_rois),3);
    curr_unclassified_mask = any(roi_masks(:,:,unclassified_rois),3);
    
    curr_map = zeros(512,512,3);
    curr_map(:,:,2) = curr_movement_mask;
    curr_map(:,:,1) = curr_quiescent_mask;
    curr_map(repmat(curr_unclassified_mask,1,1,3)) = 0.5;
    
    roi_map(:,:,:,curr_day) = curr_map;
end

m2q_grps = any((cumsum(classified == 1,2) >= 1).*(classified == -1),2);
m2q_rois = vertcat(grp{m2q_grps});
m2q_mask = any(roi_masks(:,:,m2q_rois),3);









































