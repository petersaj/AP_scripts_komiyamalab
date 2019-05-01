%% Fraction of ROIs that are active, but not classified (update classified)

% Define activity level by that of M/Q ROIs
class_avg = cell(length(data),14,3);
class_onsets = cell(length(data),14,3);
u_active_frac = nan(length(data),14);

for curr_animal = 1:length(data)
    
    classified_rois(curr_animal).unclassified_active = false(size( ...
        classified_rois(curr_animal).movement));
    
    for curr_day = 1:14       
        m = classified_rois(curr_animal).movement(:,curr_day);
        q = classified_rois(curr_animal).quiescent(:,curr_day);
        u = ~m & ~q;
        
        m_avg = nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(m,:),2);
        q_avg = nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(q,:),2);
        c_avg = nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(m|q,:),2);
        u_avg = nanmean(data(curr_animal).im(curr_day).roi_trace_thresh(u,:),2);
        
        m_onsets = sum(diff(data(curr_animal).im(curr_day).roi_trace_thresh(m,:) > 0,[],2) == 1,2);
        q_onsets = sum(diff(data(curr_animal).im(curr_day).roi_trace_thresh(q,:) > 0,[],2) == 1,2);
        c_onsets = sum(diff(data(curr_animal).im(curr_day).roi_trace_thresh(m|q,:) > 0,[],2) == 1,2);
        u_onsets = sum(diff(data(curr_animal).im(curr_day).roi_trace_thresh(u,:) > 0,[],2) == 1,2);
        
        % Define "active unclassified"
        % 1) using minimum
        %u_active = u_avg > min(c_avg) & u_onsets > min(c_onsets);
        % 2) using percentile
        u_active = u_avg > prctile(c_avg,5) & u_onsets > prctile(c_onsets,5);
        
        curr_class_avg = {m_avg,q_avg,u_avg(u_active)};
        curr_class_onsets = {m_onsets,q_onsets,u_onsets(u_active)};
        
        class_avg(curr_animal,curr_day,:) = permute(curr_class_avg,[1,3,2]);
        class_onsets(curr_animal,curr_day,:) = permute(curr_class_onsets,[1,3,2]);
        u_active_frac(curr_animal,curr_day) = sum(u_active)/length(u);
        
        % Update classified ROI structure
        classified_rois(curr_animal).unclassified_active(u,curr_day) = u_active;
    end
    disp(curr_animal);
end

% Plot median activity/onsets for each category, fraction cells
class_avg_median = cellfun(@nanmedian,class_avg);
class_onsets_median = cellfun(@nanmedian,class_onsets);

class_avg_mean = permute(nanmean(class_avg_median),[3,2,1]);
class_onsets_mean = permute(nanmean(class_onsets_median),[3,2,1]);

class_avg_sem = permute(nanstd(class_avg_median)./sqrt(sum(~isnan(class_avg_median),1)),[3,2,1]);
class_onsets_sem = permute(nanstd(class_onsets_median)./sqrt(sum(~isnan(class_onsets_median),1)),[3,2,1]);

figure; 
subplot(1,3,1);
errorbar(class_avg_mean',class_avg_sem','linewidth',2)
ylabel('Average \DeltaF/F');
xlabel('Day');
subplot(1,3,2);
errorbar(class_onsets_mean',class_onsets_sem','linewidth',2)
ylabel('# transients');
xlabel('Day');
subplot(1,3,3); hold on;
col = lines(3);
move_frac = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false));
quiescent_frac = cell2mat(cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false));
errorbar(nanmean(move_frac),nanstd(move_frac)./ ...
    sqrt(sum(~isnan(move_frac))),'color',col(1,:),'linewidth',2)
errorbar(nanmean(quiescent_frac),nanstd(quiescent_frac)./ ...
    sqrt(sum(~isnan(quiescent_frac))),'color',col(2,:),'linewidth',2)
errorbar(nanmean(u_active_frac),nanstd(u_active_frac)./ ...
    sqrt(sum(~isnan(u_active_frac))),'color',col(3,:),'linewidth',2)
ylabel('Fraction of ROIs');
xlabel('Day');
legend({'M','Q','U'});

%% Fraction of unclassified days = unclassified active

m_ua = nan(length(data),3);
q_ua = nan(length(data),3);
m2q_ua = nan(length(data),3);
for curr_animal = 1:length(data)
    m = classified_rois(curr_animal).movement;
    q = classified_rois(curr_animal).quiescent;
    ua = classified_rois(curr_animal).unclassified_active;
    
    curr_m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    
    non_mq = ~m & ~q;
    
    % week 1
    m_ua(curr_animal,1) = nanmedian(sum(ua(any(m,2) & ~curr_m2q,1:7),2)./sum(non_mq(any(m,2) & ~curr_m2q,1:7),2));
    q_ua(curr_animal,1) = nanmedian(sum(ua(any(q,2) & ~curr_m2q,1:7),2)./sum(non_mq(any(q,2) & ~curr_m2q,1:7),2));
    m2q_ua(curr_animal,1) = nanmedian(sum(ua(curr_m2q,1:7),2)./sum(non_mq(curr_m2q,1:7),2));
    
    % week 2
    m_ua(curr_animal,2) = nanmedian(sum(ua(any(m,2) & ~curr_m2q,8:14),2)./sum(non_mq(any(m,2) & ~curr_m2q,8:14),2));
    q_ua(curr_animal,2) = nanmedian(sum(ua(any(q,2) & ~curr_m2q,8:14),2)./sum(non_mq(any(q,2) & ~curr_m2q,8:14),2));
    m2q_ua(curr_animal,2) = nanmedian(sum(ua(curr_m2q,8:14),2)./sum(non_mq(curr_m2q,8:14),2));
    
    % total
    m_ua(curr_animal,3) = nanmedian(sum(ua(any(m,2) & ~curr_m2q,:),2)./sum(non_mq(any(m,2) & ~curr_m2q,:),2));
    q_ua(curr_animal,3) = nanmedian(sum(ua(any(q,2) & ~curr_m2q,:),2)./sum(non_mq(any(q,2) & ~curr_m2q,:),2));
    m2q_ua(curr_animal,3) = nanmedian(sum(ua(curr_m2q,:),2)./sum(non_mq(curr_m2q,:),2));
end

plot_mean = [nanmean(m_ua);nanmean(q_ua);nanmean(m2q_ua)];
plot_sem = [nanstd(m_ua)./sqrt(sum(~isnan(m_ua))); ...
    nanstd(q_ua)./sqrt(sum(~isnan(q_ua))); ...
    nanstd(m2q_ua)./sqrt(sum(~isnan(m2q_ua)))];

figure; hold on;
bar(plot_mean,'FaceColor','w','linewidth',2);
colormap(gray);

errorb(plot_mean,plot_sem,'.k','linewidth',2);

xlim([0,4]);
ylabel('Fraction of unclassified days active');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Movement','Quiescent','M2Q'});
legend({'Week 1','Week 2','Total'});


% Fraction of days unclassified
m_u = nan(length(data),3);
q_u = nan(length(data),3);
m2q_u = nan(length(data),3);
for curr_animal = 1:length(data)
    m = classified_rois(curr_animal).movement;
    q = classified_rois(curr_animal).quiescent;
    ua = classified_rois(curr_animal).unclassified_active;
    
    curr_m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    
    non_mq = ~m & ~q;
    
    % week 1
    m_u(curr_animal,1) = nanmedian(sum(non_mq(any(m,2) & ~curr_m2q,1:7),2)./7);
    q_u(curr_animal,1) = nanmedian(sum(non_mq(any(q,2) & ~curr_m2q,1:7),2)./7);
    m2q_u(curr_animal,1) = nanmedian(sum(non_mq(curr_m2q,1:7),2)./7);
    
    % week 2
    m_u(curr_animal,2) = nanmedian(sum(non_mq(any(m,2) & ~curr_m2q,8:14),2)./7);
    q_u(curr_animal,2) = nanmedian(sum(non_mq(any(q,2) & ~curr_m2q,8:14),2)./7);
    m2q_u(curr_animal,2) = nanmedian(sum(non_mq(curr_m2q,8:14),2)./7);
    
    % total
    m_u(curr_animal,3) = nanmedian(sum(non_mq(any(m,2) & ~curr_m2q,:),2)./14);
    q_u(curr_animal,3) = nanmedian(sum(non_mq(any(q,2) & ~curr_m2q,:),2)./14);
    m2q_u(curr_animal,3) = nanmedian(sum(non_mq(curr_m2q,:),2)./14);
end

plot_mean = [nanmean(m_u);nanmean(q_u);nanmean(m2q_u)];
plot_sem = [nanstd(m_u)./sqrt(sum(~isnan(m_u))); ...
    nanstd(q_u)./sqrt(sum(~isnan(q_u))); ...
    nanstd(m2q_u)./sqrt(sum(~isnan(m2q_u)))];

figure; hold on;
bar(plot_mean,'FaceColor','w','linewidth',2);
colormap(gray);

errorb(plot_mean,plot_sem,'.k','linewidth',2);

xlim([0,4]);
ylabel('Fraction of days unclassified');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Movement','Quiescent','M2Q'});
legend({'Week 1','Week 2','Total'});


%% Turnover from movement/quiescence to unclassified active

m2ua_overlap_grid = nan(14,14,length(data));
q2ua_overlap_grid = nan(14,14,length(data));
ua2m_overlap_grid = nan(14,14,length(data));
ua2q_overlap_grid = nan(14,14,length(data));
for curr_animal = 1:length(data)
   for curr_day1 = 1:14;
       for curr_day2 = 1:14;
           % m2ua
           active_1 = classified_rois(curr_animal).movement(:,curr_day1);           
           active_2 = classified_rois(curr_animal).unclassified_active(:,curr_day2);           
           curr_overlap = +active_1'*+active_2;
           curr_total = sum(active_1);           
           m2ua_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % q2ua
           active_1 = classified_rois(curr_animal).quiescent(:,curr_day1);           
           active_2 = classified_rois(curr_animal).unclassified_active(:,curr_day2);           
           curr_overlap = +active_1'*+active_2;
           curr_total = sum(active_1);           
           q2ua_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % ua2m
           active_1 = classified_rois(curr_animal).unclassified_active(:,curr_day1);           
           active_2 = classified_rois(curr_animal).movement(:,curr_day2);           
           curr_overlap = +active_1'*+active_2;
           curr_total = sum(active_1);           
           ua2m_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
           
           % ua2q
           active_1 = classified_rois(curr_animal).unclassified_active(:,curr_day1);           
           active_2 = classified_rois(curr_animal).quiescent(:,curr_day2);           
           curr_overlap = +active_1'*+active_2;
           curr_total = sum(active_1);           
           ua2q_overlap_grid(curr_day1,curr_day2,curr_animal) = curr_overlap./curr_total;  
       end
   end    
end

figure;
subplot(2,2,1);
imagesc(nanmean(m2ua_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('M2UA')
subplot(2,2,2);
imagesc(nanmean(q2ua_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('Q2UA')
subplot(2,2,3);
imagesc(nanmean(ua2m_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('UA2M')
subplot(2,2,4);
imagesc(nanmean(ua2q_overlap_grid,3));colormap(hot)
ylabel('Denominator ROIs')
xlabel('Overlap ROIs')
title('UA2Q')

%% M2M overlap as stacked bar
% I didn't get much out of this, so didn't put in ppt

m_stage_overlap = nan(3,3,length(data));
m_early_late_overlap = nan(3,8,length(data));
m_middle_late_overlap = nan(3,8,length(data));
for curr_animal = 1:length(data)
    
    m = classified_rois(curr_animal).movement;
    
    m_early = any(m(:,1:3),2);
    m_middle = any(m(:,4:6),2);
    m_late = any(m(:,7:end),2);
    
    m_stage_overlap(:,:,curr_animal) = ...
        [nanmean(m_early & ~m_middle),nanmean(m_early & m_middle),nanmean(~m_early & m_middle); ...
        nanmean(m_middle & ~m_late),nanmean(m_middle & m_late),nanmean(~m_middle & m_late); ...
        nanmean(m_early & ~m_late),nanmean(m_early & m_late),nanmean(~m_early & m_late)];
    
    m_early_late_overlap(:,:,curr_animal) = ...
        [nanmean(repmat(m_early,1,8) & ~m(:,7:end));
        nanmean(repmat(m_early,1,8) & m(:,7:end)); ...
        nanmean(~repmat(m_early,1,8) & m(:,7:end));];
    
    m_middle_late_overlap(:,:,curr_animal) = ...
        [nanmean(repmat(m_middle,1,8) & ~m(:,7:end));
        nanmean(repmat(m_middle,1,8) & m(:,7:end)); ...
        nanmean(~repmat(m_middle,1,8) & m(:,7:end));];
    
end

figure;bar(nanmean(m_stage_overlap,3),'stacked');


%% Correlation structure across ROIs

% All ROIs across days
corr_corr_m = nan(14,14,length(data));
corr_corr_q = nan(14,14,length(data));
corr_corr_m2q = nan(14,14,length(data));
corr_corr_ua = nan(14,14,length(data));

for curr_animal = 1:length(data);    
       
    m2q = ...
        any((cumsum(classified_rois(curr_animal).movement,2) >= 1).* ...
        classified_rois(curr_animal).quiescent,2);
    m = any(classified_rois(curr_animal).movement,2) & ~m2q;
    q = any(classified_rois(curr_animal).quiescent,2) & ~m2q;
    ua = any(classified_rois(curr_animal).unclassified_active,2) & ~m & ~q & ~m2q;
        
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_corr = nan(n_rois,n_rois,14);
    for i = 1:14
        curr_corr(:,:,i) = corrcoef(data(curr_animal).im(i).roi_trace_thresh');
    end   
    
    curr_corr_m = cell2mat(arrayfun(@(x) AP_itril(curr_corr(m,m,x),-1),1:14,'uni',false));
    corr_corr_m(:,:,curr_animal) = corrcoef(curr_corr_m,'rows','pairwise');
    
    curr_corr_q = cell2mat(arrayfun(@(x) AP_itril(curr_corr(q,q,x),-1),1:14,'uni',false));
    corr_corr_q(:,:,curr_animal) = corrcoef(curr_corr_q,'rows','pairwise');
    
    curr_corr_m2q = cell2mat(arrayfun(@(x) AP_itril(curr_corr(m2q,m2q,x),-1),1:14,'uni',false));
    corr_corr_m2q(:,:,curr_animal) = corrcoef(curr_corr_m2q,'rows','pairwise');
    
    curr_corr_ua = cell2mat(arrayfun(@(x) AP_itril(curr_corr(ua,ua,x),-1),1:14,'uni',false));
    corr_corr_ua(:,:,curr_animal) = corrcoef(curr_corr_ua,'rows','pairwise');
end

figure('Name','All ROIs across days');
subplot(2,2,1);
imagesc(nanmean(corr_corr_m,3));colormap(hot);
title('M');
subplot(2,2,2);
imagesc(nanmean(corr_corr_q,3));colormap(hot);
title('Q');
subplot(2,2,3);
imagesc(nanmean(corr_corr_m2q,3));colormap(hot);
title('M2Q');
subplot(2,2,4);
imagesc(nanmean(corr_corr_ua,3));colormap(hot);
title('UA');


% Only ROIs that share classification on pairs of days
corr_corr_m = nan(14,14,length(data));
corr_corr_q = nan(14,14,length(data));
corr_corr_m2q = nan(14,14,length(data));
corr_corr_ua = nan(14,14,length(data));

for curr_animal = 1:length(data);    

    m = classified_rois(curr_animal).movement;
    q = classified_rois(curr_animal).quiescent;
    ua = classified_rois(curr_animal).unclassified_active;
           
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_corr = nan(n_rois,n_rois,14);
    for i = 1:14
        curr_corr(:,:,i) = corrcoef(data(curr_animal).im(i).roi_trace_thresh');
    end   
    
    curr_corr_m = curr_corr;
    curr_corr_q = curr_corr;
    curr_corr_ua = curr_corr;
    for i = 1:14
        curr_corr_m(~m(:,i),~m(:,i),i) = NaN;
        curr_corr_q(~q(:,i),~q(:,i),i) = NaN;
        curr_corr_ua(~ua(:,i),~ua(:,i),i) = NaN;
    end
        
    corr_corr_m(:,:,curr_animal) = corrcoef(cell2mat(arrayfun(@(x) ...
        AP_itril(curr_corr_m(:,:,x),-1),1:14,'uni',false)),'rows','pairwise');
    corr_corr_q(:,:,curr_animal) = corrcoef(cell2mat(arrayfun(@(x) ...
        AP_itril(curr_corr_q(:,:,x),-1),1:14,'uni',false)),'rows','pairwise');
    corr_corr_ua(:,:,curr_animal) = corrcoef(cell2mat(arrayfun(@(x) ...
        AP_itril(curr_corr_ua(:,:,x),-1),1:14,'uni',false)),'rows','pairwise');

end

figure('Name','Shared classified ROIs across days');
subplot(1,3,1);
imagesc(nanmean(corr_corr_m,3));colormap(hot);
title('M');
subplot(1,3,2);
imagesc(nanmean(corr_corr_q,3));colormap(hot);
title('Q');
subplot(1,3,3);
imagesc(nanmean(corr_corr_ua,3));colormap(hot);
title('UA');


% Correlation across classifications on given day comparisons
corr_corr_mq = nan(14,14,length(data));
corr_corr_mua = nan(14,14,length(data));
corr_corr_qua = nan(14,14,length(data));

for curr_animal = 1:length(data);    

    m = classified_rois(curr_animal).movement;
    q = classified_rois(curr_animal).quiescent;
    ua = classified_rois(curr_animal).unclassified_active;
           
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_corr = nan(n_rois,n_rois,14);
    for i = 1:14
        curr_corr(:,:,i) = corrcoef(data(curr_animal).im(i).roi_trace_thresh');
    end   
    
    curr_corr_mq = curr_corr;
    curr_corr_mua = curr_corr;
    curr_corr_qua = curr_corr;
    for i = 1:14
        curr_corr_mq(~m(:,i),:,i) = NaN;
        curr_corr_mq(:,~q(:,i),i) = NaN;
        
        curr_corr_mua(~m(:,i),:,i) = NaN;
        curr_corr_mua(:,~ua(:,i),i) = NaN;
        
        curr_corr_qua(~q(:,i),:,i) = NaN;
        curr_corr_qua(:,~ua(:,i),i) = NaN;
    end
        
    
    corr_corr_mq(:,:,curr_animal) = corrcoef(reshape(curr_corr_mq,[],14),'rows','pairwise');
    corr_corr_mua(:,:,curr_animal) = corrcoef(reshape(curr_corr_mua,[],14),'rows','pairwise');
    corr_corr_qua(:,:,curr_animal) = corrcoef(reshape(curr_corr_qua,[],14),'rows','pairwise');

end

figure('Name','Across classification');
subplot(1,3,1);
imagesc(nanmean(corr_corr_mq,3));colormap(hot);
title('M-Q');
subplot(1,3,2);
imagesc(nanmean(corr_corr_mua,3));colormap(hot);
title('M-UA');
subplot(1,3,3);
imagesc(nanmean(corr_corr_qua,3));colormap(hot);
title('Q-UA');



%% Classifier predict binary move from binary activity

n_split = 5;

correct_move_classify_m = nan(length(data),14);
correct_move_classify_q = nan(length(data),14);
correct_move_classify_ua = nan(length(data),14);

confusion_m = nan(2,2,14,length(data));
confusion_q = nan(2,2,14,length(data));
confusion_ua = nan(2,2,14,length(data));

for curr_animal = 1:length(data);
    for curr_session = 1:14;
        
        m = classified_rois(curr_animal).movement(:,curr_session);
        q = classified_rois(curr_animal).quiescent(:,curr_session);
        ua = classified_rois(curr_animal).unclassified_active(:,curr_session);
        
        n_frames = size(data(curr_animal).im(curr_session).roi_trace_df,2);
        
        n_split_frames = diff(ceil(linspace(0,n_frames,n_split+1)));
        split_frames = mat2cell(randperm(n_frames),1,n_split_frames);
        
        lever_move_predict_m = nan(size(analysis(curr_animal).lever(curr_session).lever_move_frames));
        lever_move_predict_q = nan(size(analysis(curr_animal).lever(curr_session).lever_move_frames));
        lever_move_predict_ua = nan(size(analysis(curr_animal).lever(curr_session).lever_move_frames));
        
        for curr_split = 1:n_split;
            
            curr_train_frames = horzcat(split_frames{setdiff(1:n_split,curr_split)});
            curr_test_frames = split_frames{curr_split};
            
            curr_data = data(curr_animal).im(curr_session).roi_trace_df;
            
            % Movement ROIs
            nb_m = NaiveBayes.fit(curr_data(m,curr_train_frames)', ...
                +analysis(curr_animal).lever(curr_session).lever_move_frames(curr_train_frames)');
            
            lever_move_predict_m(curr_test_frames) = ...
                predict(nb_m,curr_data(m,curr_test_frames)');
            
            % Quiescent ROIs
            nb_q = NaiveBayes.fit(curr_data(q,curr_train_frames)', ...
                +analysis(curr_animal).lever(curr_session).lever_move_frames(curr_train_frames)');
            
            lever_move_predict_q(curr_test_frames) = ...
                predict(nb_q,curr_data(q,curr_test_frames)');
            
            % Unclassified active ROIs
            nb_ua = NaiveBayes.fit(curr_data(ua,curr_train_frames)', ...
                +analysis(curr_animal).lever(curr_session).lever_move_frames(curr_train_frames)');
            
            lever_move_predict_ua(curr_test_frames) = ...
                predict(nb_ua,curr_data(ua,curr_test_frames)');
                       
        end
        
        % Get fraction correct
        correct_move_classify_m(curr_animal,curr_session) = ...
            nanmean(lever_move_predict_m == analysis(curr_animal).lever(curr_session).lever_move_frames);
        
        correct_move_classify_q(curr_animal,curr_session) = ...
            nanmean(lever_move_predict_q == analysis(curr_animal).lever(curr_session).lever_move_frames);
        
        correct_move_classify_ua(curr_animal,curr_session) = ...
            nanmean(lever_move_predict_ua == analysis(curr_animal).lever(curr_session).lever_move_frames);
        
        % Get confusion matricies
        confusion_m(:,:,curr_session,curr_animal) = confusionmat(+analysis(curr_animal). ...
            lever(curr_session).lever_move_frames',lever_move_predict_m);
        confusion_q(:,:,curr_session,curr_animal) = confusionmat(+analysis(curr_animal). ...
            lever(curr_session).lever_move_frames',lever_move_predict_q);
        confusion_ua(:,:,curr_session,curr_animal) = confusionmat(+analysis(curr_animal). ...
            lever(curr_session).lever_move_frames',lever_move_predict_ua);
        
    end
    disp(curr_animal);
end

figure;

% Plot movement classification over days
subplot(1,2,1); hold on;
errorbar(nanmean(correct_move_classify_m), ...
    nanstd(correct_move_classify_m)./sqrt(sum(~isnan(correct_move_classify_m))),'g','linewidth',2);
errorbar(nanmean(correct_move_classify_q), ...
    nanstd(correct_move_classify_q)./sqrt(sum(~isnan(correct_move_classify_q))),'r','linewidth',2);
errorbar(nanmean(correct_move_classify_ua), ...
    nanstd(correct_move_classify_ua)./sqrt(sum(~isnan(correct_move_classify_ua))),'k','linewidth',2);

xlabel('Day');
ylabel('Fraction correct lever prediction')
legend({'M','Q','UA'});

% Plot confusion
confusion_m_norm = nanmean(bsxfun(@times,confusion_m,1./sum(confusion_m,2)),3);
m_misclass_mean = [nanmean(confusion_m_norm(1,2,:),3), ...
    nanmean(confusion_m_norm(2,1,:),3)];
m_misclass_sem = [nanstd(confusion_m_norm(1,2,:),[],3)./sqrt(sum(~isnan(confusion_m_norm(1,2,:)),3)), ...
    nanstd(confusion_m_norm(2,1,:),[],3)./sqrt(sum(~isnan(confusion_m_norm(2,1,:)),3))];

confusion_q_norm = nanmean(bsxfun(@times,confusion_q,1./sum(confusion_q,2)),3);
q_misclass_mean = [nanmean(confusion_q_norm(1,2,:),3), ...
    nanmean(confusion_q_norm(2,1,:),3)];
q_misclass_sem = [nanstd(confusion_q_norm(1,2,:),[],3)./sqrt(sum(~isnan(confusion_q_norm(1,2,:)),3)), ...
    nanstd(confusion_q_norm(2,1,:),[],3)./sqrt(sum(~isnan(confusion_q_norm(2,1,:)),3))];

confusion_ua_norm = nanmean(bsxfun(@times,confusion_ua,1./sum(confusion_ua,2)),3);
ua_misclass_mean = [nanmean(confusion_ua_norm(1,2,:),3), ...
    nanmean(confusion_ua_norm(2,1,:),3)];
ua_misclass_sem = [nanstd(confusion_ua_norm(1,2,:),[],3)./sqrt(sum(~isnan(confusion_ua_norm(1,2,:)),3)), ...
    nanstd(confusion_ua_norm(2,1,:),[],3)./sqrt(sum(~isnan(confusion_ua_norm(2,1,:)),3))];

mean_misclass = [m_misclass_mean;q_misclass_mean;ua_misclass_mean];
sem_misclass = [m_misclass_sem;q_misclass_sem;ua_misclass_sem];

subplot(1,2,2); hold on;
bar(plot_missclass);
errorb(mean_misclass,sem_misclass,'.k','linewidth',2);
colormap(gray)

set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'M','Q','UA'});
ylabel('Fraction movement misclassified');
legend({'Quiescent as Move','Move as Quiescent'});




%% M/Q xcov

% M -> Q
maxlag = 500;
mean_m_q_xcov = nan(14,2*maxlag+1,length(data));

for curr_animal = 1:length(data);
    for curr_session = 1:14;
        
        curr_mean_m = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).movement(:,curr_session),:),1);
        
        curr_mean_q = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).quiescent(:,curr_session),:),1);
        
        curr_mean_ua = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).unclassified_active(:,curr_session),:),1);
        
        mean_m_q_xcov(curr_session,:,curr_animal) = ...
            xcov(curr_mean_m,curr_mean_q,maxlag);
        
    end
end
mean_m_q_xcov_mean = mean(mean_m_q_xcov,3);
t = (-maxlag:maxlag)/28;

figure;
% Plot xcov across days (doesn't seem to change)
col = jet(14);
subplot(2,1,1); hold on
for i = 1:14
    plot(t,mean_m_q_xcov_mean(i,:),'color',col(i,:))
end
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Q lag time (s)');
ylabel('M -> Q xcov')

subplot(2,1,2); hold on;
mean_m_q_xcov_allmean = mean(mean_m_q_xcov_mean,1);
plot(t,mean_m_q_xcov_allmean,'k');
line([0,0],ylim,'color','r','linestyle','--');
xlabel('Q lag time (s)');
ylabel('M -> Q xcov')

% M/Q -> UA
maxlag = 500;
mean_m_ua_xcov = nan(14,2*maxlag+1,length(data));
mean_q_ua_xcov = nan(14,2*maxlag+1,length(data));

for curr_animal = 1:length(data);
    for curr_session = 1:14;
        
        curr_mean_m = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).movement(:,curr_session),:),1);
        
        curr_mean_q = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).quiescent(:,curr_session),:),1);
        
        curr_mean_ua = nanmean(data(curr_animal).im(curr_session). ...
            roi_trace_thresh(classified_rois(curr_animal).unclassified_active(:,curr_session),:),1);
        
        mean_m_ua_xcov(curr_session,:,curr_animal) = ...
            xcov(curr_mean_m,curr_mean_ua,maxlag);
        
        mean_q_ua_xcov(curr_session,:,curr_animal) = ...
            xcov(curr_mean_q,curr_mean_ua,maxlag);
        
    end
end
mean_m_ua_xcov_mean = mean(mean_m_ua_xcov,3);
mean_q_ua_xcov_mean = mean(mean_q_ua_xcov,3);
t = (-maxlag:maxlag)/28;

figure; 
subplot(2,1,1);hold on;
mean_m_ua_xcov_allmean = mean(mean_m_ua_xcov_mean,1);
plot(t,mean_m_ua_xcov_allmean,'k');
line([0,0],ylim,'color','r','linestyle','--');
xlabel('UA lag time (s)');
ylabel('M -> UA xcov')

subplot(2,1,2);hold on;
mean_q_ua_xcov_allmean = mean(mean_q_ua_xcov_mean,1);
plot(t,mean_q_ua_xcov_allmean,'k');
line([0,0],ylim,'color','r','linestyle','--');
xlabel('UA lag time (s)');
ylabel('Q -> UA xcov')


%% Timing / jitter of activity relative to movement onset/offset

roi_timing_allmove = cell(14,length(data));
roi_timing_aligned = cell(14,length(data));
roi_timing_xcorr = cell(14,length(data));

for curr_animal = 1:length(data);
    
    for curr_day = 1:14;
        n_rois = size(data(curr_animal).im(curr_day).roi_trace_df,1);
        
        curr_move_diff = diff(analysis(curr_animal).lever(curr_day).lever_move_frames);
        
        % ALL EVENTS / MOVEMENTS
        % store as rois x event x [mean,std]
        curr_roi_timing_allmove = nan(n_rois,4,2);
        for curr_roi = 1:n_rois
            
            % Get activity onsets
            curr_act_onsets = find(diff([0 data(curr_animal).im(curr_day).roi_trace_thresh(curr_roi,:)>0]) == 1);
            
            % For each activity event: find closest:
            % 1) last move start
            last_move_start = arrayfun(@(x) find([curr_move_diff( ...
                curr_act_onsets(x):-1:1);1] == 1,1),1:length(curr_act_onsets));
            % 2) last move stop
            last_move_stop = arrayfun(@(x) find([curr_move_diff( ...
                curr_act_onsets(x):-1:1);-1] == -1,1),1:length(curr_act_onsets));
            % 3) next move start
            next_move_start = arrayfun(@(x) find([curr_move_diff( ...
                curr_act_onsets(x):end);1] == 1,1),1:length(curr_act_onsets));
            % 4) next move stop
            next_move_stop = arrayfun(@(x) find([curr_move_diff( ...
                curr_act_onsets(x):end);-1] == -1,1),1:length(curr_act_onsets));
            
            curr_timing = [last_move_start;last_move_stop;next_move_start;next_move_stop];
            
            % Store mean/std
            curr_roi_timing_allmove(curr_roi,:,1) = nanmedian(curr_timing,2);
            curr_roi_timing_allmove(curr_roi,:,2) = std(curr_timing,[],2);
            
        end
        roi_timing_allmove{curr_day,curr_animal} = curr_roi_timing_allmove;
        
        % ALIGNED C/R MOVEMENTS
        curr_roi_timing_aligned = nan(n_rois,2,2);
        curr_align_frame = find(analysis(curr_animal).surrounding_frames > 0,1);
        for curr_roi = 1:n_rois
            
            [curr_max,curr_acttime_movestart] = max(analysis(curr_animal).im(curr_day). ...
                move_onset_aligned(:,curr_align_frame:end,curr_roi) > 0,[],2);
            curr_acttime_movestart(curr_max == 0) = [];
            
            [curr_max,curr_acttime_movestop] = max(analysis(curr_animal).im(curr_day). ...
                move_offset_aligned(:,curr_align_frame:end,curr_roi) > 0,[],2);
            curr_acttime_movestop(curr_max == 0) = [];
            
            curr_roi_timing_aligned(curr_roi,1,1) = nanmedian(curr_acttime_movestart);
            curr_roi_timing_aligned(curr_roi,2,1) = nanmedian(curr_acttime_movestop);
            curr_roi_timing_aligned(curr_roi,1,2) = nanstd(curr_acttime_movestart);
            curr_roi_timing_aligned(curr_roi,2,2) = nanstd(curr_acttime_movestop);
            
        end
        roi_timing_aligned{curr_day,curr_animal} = curr_roi_timing_aligned;
        
        % XCORR OF ALL MOVEMENTS
        curr_roi_timing_xcorr = nan(n_rois,2,2);
        max_frames = 90;
        for curr_roi = 1:n_rois
            
            move_onset_xcorr = xcorr(curr_move_diff == 1, ...
                data(curr_animal).im(curr_day).roi_trace_thresh(curr_roi,:),max_frames);
            move_offset_xcorr = xcorr(curr_move_diff == -1, ...
                data(curr_animal).im(curr_day).roi_trace_thresh(curr_roi,:),max_frames);
            
            [curr_onset_max curr_onset_max_time] = max(move_onset_xcorr);
            [curr_offset_max curr_offset_max_time] = max(move_offset_xcorr);
            
            curr_roi_timing_xcorr(curr_roi,1,1) = curr_onset_max_time - max_frames;
            curr_roi_timing_xcorr(curr_roi,2,1) = curr_offset_max_time - max_frames;
            
            curr_roi_timing_xcorr(curr_roi,1,2) = curr_onset_max;
            curr_roi_timing_xcorr(curr_roi,2,2) = curr_offset_max;
            
        end
        roi_timing_xcorr{curr_day,curr_animal} = curr_roi_timing_xcorr;
                
    end
    
    disp(curr_animal);
    
end


%% Get statistics of movement (all movements)

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



%% Get statistics of movement (no movements around reward)

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
        [lever_active,~,~,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated(lever_move,true);
        
        % Get reward times
        first_frame_time = round(data(curr_animal).bhv(curr_session).frame_times(1)* ...
            (data(curr_animal).bhv(curr_session).xsg_sample_rate/10));
        
        rewarded_trials = cellfun(@(x) ~isempty(x.states.reward), ...
            data(curr_animal).bhv(curr_session).bhv_times);

        reward_frames = round(cellfun(@(x) x.states.reward(1), ...
            data(curr_animal).bhv(curr_session).bhv_frames(rewarded_trials)));
        
        reward_times = round(data(curr_animal).bhv(curr_session).frame_times(reward_frames)*1000);
        
        reward_times_offset = reward_times - first_frame_time;
        
        % Split movements/nonmovements
        boundary_times = find(diff([Inf;lever_active;Inf]) ~= 0);
        lever_move_split = mat2cell(lever_move,diff(boundary_times));
        lever_envelope_split = mat2cell(lever_velocity_envelope_smooth,diff(boundary_times));

        lever_move_split_move = cellfun(@any,mat2cell(lever_active,diff(boundary_times)));
        lever_times = mat2cell((1:length(lever_move))',diff(boundary_times));
        
        % Remove movements which are +/- time of reward
        reward_surround = 500; % ms
        reward_surround_times = cell2mat(arrayfun(@(x) x-reward_surround: ...
            x+reward_surround,reward_times_offset','uni',false));
        
        near_reward_movements = cellfun(@(x) ~isempty( ...
            intersect(reward_surround_times, x)), ...
            lever_times) & lever_move_split_move;
        
        lever_move_split(near_reward_movements) = [];
        lever_envelope_split(near_reward_movements) = [];
        lever_move_split_move(near_reward_movements) = [];
        
        % fraction time spent moving (non near-reward movement time)
        lever_active_nonreward = lever_active(vertcat(lever_times{~near_reward_movements}));
        frac_moving = nanmean(lever_active_nonreward);
        
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

figure('Name','Without near-reward movements');

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
















