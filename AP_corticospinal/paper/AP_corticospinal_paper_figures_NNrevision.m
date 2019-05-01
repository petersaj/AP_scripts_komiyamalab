%% Reviewer 1 point 1: fig 5A, fraction classified cells

clearvars -except data analysis classified_rois

m_frac = cellfun(@(x) nanmean(x,1),{classified_rois(:).movement}','uni',false);
m_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),m_frac,'uni',false));

q_frac = cellfun(@(x) nanmean(x,1),{classified_rois(:).quiescent}','uni',false);
q_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),q_frac,'uni',false));

ua_frac = cellfun(@(x) nanmean(x,1),{classified_rois(:).unclassified_active}','uni',false);
ua_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),ua_frac,'uni',false));

s_frac = cellfun(@(m,q,ua) nanmean(m==0 & q==0 & ua==0,1), ...
    {classified_rois(:).movement}', ...
    {classified_rois(:).quiescent}', ...
    {classified_rois(:).unclassified_active}','uni',false);
s_frac_pad = cell2mat(cellfun(@(x) padarray(x(1:min(14,length(x))), ...
    [0,14-min(14,length(x))],NaN,'post'),s_frac,'uni',false));

figure; hold on
errorbar(nanmean(m_frac_pad),nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad))),'k','linewidth',2)
errorbar(nanmean(q_frac_pad),nanstd(q_frac_pad)./sqrt(sum(~isnan(q_frac_pad))),'r','linewidth',2)
errorbar(nanmean(ua_frac_pad),nanstd(ua_frac_pad)./sqrt(sum(~isnan(ua_frac_pad))),'b','linewidth',2)
errorbar(nanmean(s_frac_pad),nanstd(s_frac_pad)./sqrt(sum(~isnan(s_frac_pad))),'m','linewidth',2)
ylabel('Fraction of ROIs');
xlabel('Day');
legend({'Movement-related','Quiescence-related','Unclassified active','Silent'})

% Plot as a cumulatve plot
mean_cat = ...
    [nanmean(m_frac_pad,1);nanmean(q_frac_pad,1);nanmean(ua_frac_pad,1);nanmean(s_frac_pad,1)];
sem_cat = ...
    [nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad))); ...
    nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad))); ...
    nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad))); ...
    nanstd(m_frac_pad)./sqrt(sum(~isnan(m_frac_pad)))];
figure; hold on;
area(mean_cat');
errorbar(cumsum(mean_cat,1)',sem_cat','.k','linewidth',2);
ylim([0 1]);
xlim([1 14]);
legend({'M','Q','UA','S'});
% Statistics

days_1 = 1:6;
days_2 = 6:14;

% Quiescent
p = anova1(q_frac_pad,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(q_frac_pad(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(q_frac_pad(:,days_2),[],1));

% Movement
p = anova1(m_frac_pad,[],'off');

[r,p] = corrcoef(reshape(repmat(days_1,length(data),1),[],1), ...
    reshape(m_frac_pad(:,days_1),[],1));

[r,p] = corrcoef(reshape(repmat(days_2,length(data),1),[],1), ...
    reshape(m_frac_pad(:,days_2),[],1));

% Early quiescent
m_frac_pad_zscore = bsxfun(@rdivide,bsxfun(@minus,m_frac_pad,nanmean(m_frac_pad,2)),nanstd(m_frac_pad,[],2));
q_frac_pad_zscore = bsxfun(@rdivide,bsxfun(@minus,q_frac_pad,nanmean(q_frac_pad,2)),nanstd(q_frac_pad,[],2));

pm = signrank(nanmean(m_frac_pad_zscore(:,1:2),2),nanmean(m_frac_pad_zscore(:,3:4),2));
pq = signrank(nanmean(q_frac_pad_zscore(:,1:2),2),nanmean(q_frac_pad_zscore(:,3:4),2));


%% Reviewer 1 point 1: fig 5D, average day activity of switching cells

clearvars -except data analysis classified_rois

% Get activity of all cells during m/q/total 
avg_movement_act = cell(length(data),1);
avg_quiescent_act = cell(length(data),1);
avg_total_act = cell(length(data),1);

for curr_animal = 1:length(data)
        
    use_cells = true(size(classified_rois(curr_animal).movement,1),1);
        
    curr_avg_movement_act = nan(length(use_cells),14);
    curr_avg_quiescent_act = nan(length(use_cells),14);
    curr_avg_total_act = nan(length(use_cells),14);
    
     for i = 1:14
         
         if length(data(curr_animal).im) < i || ...
                 isempty(data(curr_animal).im(i).roi_trace_thresh)
             continue
         end
         
        curr_avg_movement_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames > 0),2);
        
        curr_avg_quiescent_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames == 0),2);
        
        curr_avg_total_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,:),2);
        
     end
    
    avg_movement_act{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act{curr_animal} = curr_avg_quiescent_act;
    avg_total_act{curr_animal} = curr_avg_total_act;
    
    disp(curr_animal);
    
end

% Get categories
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) ~m & ~q & ~ua,m_pad,q_pad,ua_pad,'uni',false);

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;

all_cells = cellfun(@(m) true(size(m,1),1),m_pad,'uni',false);
m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);
q2m = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);

m2m = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);
q2q = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);

m2u = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) &...
    sum(m(:,8:14),2) < min_classdays,m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2) & ...
    sum(m(:,1:7),2) < min_classdays,m_pad,q_pad,'uni',false);

q2u = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) < min_classdays,m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2) & ...
    sum(q(:,1:7),2) < min_classdays,m_pad,q_pad,'uni',false);

ua2u = cellfun(@(ua,m,q) ...
    sum(ua(:,1:7),2) >= sum(m(:,1:7),2) & ...
    sum(ua(:,1:7),2) >= sum(q(:,1:7),2) & ...
    sum(ua(:,8:end),2) < min_classdays,ua_pad,m_pad,q_pad,'uni',false);
u2ua = cellfun(@(ua,m,q) ...
    sum(ua(:,8:14),2) >= sum(m(:,8:14),2) & ...
    sum(ua(:,8:14),2) >= sum(q(:,8:14),2) & ...
    sum(ua(:,1:7),2) < min_classdays,ua_pad,m_pad,q_pad,'uni',false);

plot_classes = {all_cells,m2m,q2q,m2q,q2m,m2u,u2m,q2u,u2q,ua2u,u2ua};
avg_act = cell(length(plot_classes),1);
for curr_class_idx = 1:length(plot_classes)
    
    curr_class = plot_classes{curr_class_idx};
    
    avg_act{curr_class_idx}{1} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_movement_act,curr_class','uni',false);
    avg_act{curr_class_idx}{2} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_quiescent_act,curr_class','uni',false);
    avg_act{curr_class_idx}{3} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_total_act,curr_class','uni',false);
    
end

avg_act_mean = cellfun(@(x) cellfun(@(x) ...
    nanmean(vertcat(x{:}),1),x,'uni',false),avg_act,'uni',false);
avg_act_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(vertcat(x{:}),[],1)./sqrt(sum(~isnan(vertcat(x{:})),1)), ...
    x,'uni',false),avg_act,'uni',false);

class_labels = {'All','M2M','Q2Q','M2Q','Q2M','M2~M','~M2M','Q2~Q','~Q2Q','UA2~UA','~UA2UA'};
figure;
for curr_class_idx = 1:length(plot_classes)
    subplot(4,3,curr_class_idx); hold on;
    errorbar(avg_act_mean{curr_class_idx}{1},avg_act_sem{curr_class_idx}{1},'r');
    errorbar(avg_act_mean{curr_class_idx}{2},avg_act_sem{curr_class_idx}{2},'b');
    errorbar(avg_act_mean{curr_class_idx}{3},avg_act_sem{curr_class_idx}{3},'k');
    
    xlabel('day')
    ylabel('average df/f')
    title(class_labels{curr_class_idx})
    legend({'total','during movement','during quiescence'})
end


%% Reviewer 1 point 1: fig 5D: classification of switch cells by day

clearvars -except data analysis classified_rois

% Plot cagetorization of m2q cells

% Get categories
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;
m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) >= min_classdays & ...
    sum(q(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
m2u = cellfun(@(m,q) ...
    sum(m(:,1:7),2) >= min_classdays & ...
    sum(m(:,8:end),2) < min_classdays,m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    sum(m(:,1:7),2) < min_classdays & ...
    sum(m(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
q2u = cellfun(@(m,q) ...
    sum(q(:,1:7),2) >= min_classdays & ...
    sum(q(:,8:end),2) < min_classdays,m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    sum(q(:,1:7),2) < min_classdays & ...
    sum(q(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
ua2u = cellfun(@(ua) ...
    sum(ua(:,1:7),2) >= min_classdays & ...
    sum(ua(:,8:end),2) < min_classdays,ua_pad,'uni',false);
u2ua = cellfun(@(ua) ...
    sum(ua(:,1:7),2) < min_classdays & ...
    sum(ua(:,8:end),2) >= min_classdays,ua_pad,'uni',false);

m2q_frac_m = cell2mat(cellfun(@(m,m2q) nanmean(m(m2q,:),1),m_pad,m2q,'uni',false)');
m2q_frac_q = cell2mat(cellfun(@(q,m2q) nanmean(q(m2q,:),1),q_pad,m2q,'uni',false)');
m2q_frac_ua = cell2mat(cellfun(@(ua,m2q) nanmean(ua(m2q,:),1),ua_pad,m2q,'uni',false)');
m2q_frac_s = cell2mat(cellfun(@(m,q,ua,m2q) nanmean(~m(m2q,:)&~q(m2q,:)&~ua(m2q,:),1),m_pad,q_pad,ua_pad,m2q,'uni',false)');

m2u_frac_m = cell2mat(cellfun(@(m,m2u) nanmean(m(m2u,:),1),m_pad,m2u,'uni',false)');
m2u_frac_q = cell2mat(cellfun(@(q,m2u) nanmean(q(m2u,:),1),q_pad,m2u,'uni',false)');
m2u_frac_ua = cell2mat(cellfun(@(ua,m2u) nanmean(ua(m2u,:),1),ua_pad,m2u,'uni',false)');
m2u_frac_s = cell2mat(cellfun(@(m,q,ua,m2u) nanmean(~m(m2u,:)&~q(m2u,:)&~ua(m2u,:),1),m_pad,q_pad,ua_pad,m2u,'uni',false)');

u2m_frac_m = cell2mat(cellfun(@(m,u2m) nanmean(m(u2m,:),1),m_pad,u2m,'uni',false)');
u2m_frac_q = cell2mat(cellfun(@(q,u2m) nanmean(q(u2m,:),1),q_pad,u2m,'uni',false)');
u2m_frac_ua = cell2mat(cellfun(@(ua,u2m) nanmean(ua(u2m,:),1),ua_pad,u2m,'uni',false)');
u2m_frac_s = cell2mat(cellfun(@(m,q,ua,u2m) nanmean(~m(u2m,:)&~q(u2m,:)&~ua(u2m,:),1),m_pad,q_pad,ua_pad,u2m,'uni',false)');

q2u_frac_m = cell2mat(cellfun(@(m,q2u) nanmean(m(q2u,:),1),m_pad,q2u,'uni',false)');
q2u_frac_q = cell2mat(cellfun(@(q,q2u) nanmean(q(q2u,:),1),q_pad,q2u,'uni',false)');
q2u_frac_ua = cell2mat(cellfun(@(ua,q2u) nanmean(ua(q2u,:),1),ua_pad,q2u,'uni',false)');
q2u_frac_s = cell2mat(cellfun(@(m,q,ua,q2u) nanmean(~m(q2u,:)&~q(q2u,:)&~ua(q2u,:),1),m_pad,q_pad,ua_pad,q2u,'uni',false)');

u2q_frac_m = cell2mat(cellfun(@(m,u2q) nanmean(m(u2q,:),1),m_pad,u2q,'uni',false)');
u2q_frac_q = cell2mat(cellfun(@(q,u2q) nanmean(q(u2q,:),1),q_pad,u2q,'uni',false)');
u2q_frac_ua = cell2mat(cellfun(@(ua,u2q) nanmean(ua(u2q,:),1),ua_pad,u2q,'uni',false)');
u2q_frac_s = cell2mat(cellfun(@(m,q,ua,u2q) nanmean(~m(u2q,:)&~q(u2q,:)&~ua(u2q,:),1),m_pad,q_pad,ua_pad,u2q,'uni',false)');

ua2u_frac_m = cell2mat(cellfun(@(m,ua2u) nanmean(m(ua2u,:),1),m_pad,ua2u,'uni',false)');
ua2u_frac_q = cell2mat(cellfun(@(q,ua2u) nanmean(q(ua2u,:),1),q_pad,ua2u,'uni',false)');
ua2u_frac_ua = cell2mat(cellfun(@(ua,ua2u) nanmean(ua(ua2u,:),1),ua_pad,ua2u,'uni',false)');
ua2u_frac_s = cell2mat(cellfun(@(m,q,ua,ua2u) nanmean(~m(ua2u,:)&~q(ua2u,:)&~ua(ua2u,:),1),m_pad,q_pad,ua_pad,ua2u,'uni',false)');

u2ua_frac_m = cell2mat(cellfun(@(m,u2ua) nanmean(m(u2ua,:),1),m_pad,u2ua,'uni',false)');
u2ua_frac_q = cell2mat(cellfun(@(q,u2ua) nanmean(q(u2ua,:),1),q_pad,u2ua,'uni',false)');
u2ua_frac_ua = cell2mat(cellfun(@(ua,u2ua) nanmean(ua(u2ua,:),1),ua_pad,u2ua,'uni',false)');
u2ua_frac_s = cell2mat(cellfun(@(m,q,ua,u2ua) nanmean(~m(u2ua,:)&~q(u2ua,:)&~ua(u2ua,:),1),m_pad,q_pad,ua_pad,u2ua,'uni',false)');

figure; 
subplot(2,4,1);hold on;
errorbar(nanmean(m2q_frac_m,1),nanstd(m2q_frac_m,[],1)./sqrt(sum(~isnan(m2q_frac_m))),'k','linewidth',2);
errorbar(nanmean(m2q_frac_q,1),nanstd(m2q_frac_q,[],1)./sqrt(sum(~isnan(m2q_frac_q))),'r','linewidth',2);
errorbar(nanmean(m2q_frac_ua,1),nanstd(m2q_frac_ua,[],1)./sqrt(sum(~isnan(m2q_frac_ua))),'b','linewidth',2);
errorbar(nanmean(m2q_frac_s,1),nanstd(m2q_frac_s,[],1)./sqrt(sum(~isnan(m2q_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of M2Q');

subplot(2,4,2);hold on;
errorbar(nanmean(m2u_frac_m,1),nanstd(m2u_frac_m,[],1)./sqrt(sum(~isnan(m2u_frac_m))),'k','linewidth',2);
errorbar(nanmean(m2u_frac_q,1),nanstd(m2u_frac_q,[],1)./sqrt(sum(~isnan(m2u_frac_q))),'r','linewidth',2);
errorbar(nanmean(m2u_frac_ua,1),nanstd(m2u_frac_ua,[],1)./sqrt(sum(~isnan(m2u_frac_ua))),'b','linewidth',2);
errorbar(nanmean(m2u_frac_s,1),nanstd(m2u_frac_s,[],1)./sqrt(sum(~isnan(m2u_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of M2~M');

subplot(2,4,3);hold on;
errorbar(nanmean(u2m_frac_m,1),nanstd(u2m_frac_m,[],1)./sqrt(sum(~isnan(u2m_frac_m))),'k','linewidth',2);
errorbar(nanmean(u2m_frac_q,1),nanstd(u2m_frac_q,[],1)./sqrt(sum(~isnan(u2m_frac_q))),'r','linewidth',2);
errorbar(nanmean(u2m_frac_ua,1),nanstd(u2m_frac_ua,[],1)./sqrt(sum(~isnan(u2m_frac_ua))),'b','linewidth',2);
errorbar(nanmean(u2m_frac_s,1),nanstd(u2m_frac_s,[],1)./sqrt(sum(~isnan(u2m_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of ~M2M');

subplot(2,4,4);hold on;
errorbar(nanmean(q2u_frac_m,1),nanstd(q2u_frac_m,[],1)./sqrt(sum(~isnan(q2u_frac_m))),'k','linewidth',2);
errorbar(nanmean(q2u_frac_q,1),nanstd(q2u_frac_q,[],1)./sqrt(sum(~isnan(q2u_frac_q))),'r','linewidth',2);
errorbar(nanmean(q2u_frac_ua,1),nanstd(q2u_frac_ua,[],1)./sqrt(sum(~isnan(q2u_frac_ua))),'b','linewidth',2);
errorbar(nanmean(q2u_frac_s,1),nanstd(q2u_frac_s,[],1)./sqrt(sum(~isnan(q2u_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of Q2~Q');

subplot(2,4,5);hold on;
errorbar(nanmean(u2q_frac_m,1),nanstd(u2q_frac_m,[],1)./sqrt(sum(~isnan(u2q_frac_m))),'k','linewidth',2);
errorbar(nanmean(u2q_frac_q,1),nanstd(u2q_frac_q,[],1)./sqrt(sum(~isnan(u2q_frac_q))),'r','linewidth',2);
errorbar(nanmean(u2q_frac_ua,1),nanstd(u2q_frac_ua,[],1)./sqrt(sum(~isnan(u2q_frac_ua))),'b','linewidth',2);
errorbar(nanmean(u2q_frac_s,1),nanstd(u2q_frac_s,[],1)./sqrt(sum(~isnan(u2q_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of ~Q2Q');

subplot(2,4,6);hold on;
errorbar(nanmean(ua2u_frac_m,1),nanstd(ua2u_frac_m,[],1)./sqrt(sum(~isnan(ua2u_frac_m))),'k','linewidth',2);
errorbar(nanmean(ua2u_frac_q,1),nanstd(ua2u_frac_q,[],1)./sqrt(sum(~isnan(ua2u_frac_q))),'r','linewidth',2);
errorbar(nanmean(ua2u_frac_ua,1),nanstd(ua2u_frac_ua,[],1)./sqrt(sum(~isnan(ua2u_frac_ua))),'b','linewidth',2);
errorbar(nanmean(ua2u_frac_s,1),nanstd(ua2u_frac_s,[],1)./sqrt(sum(~isnan(ua2u_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of UA2~UA');

subplot(2,4,7);hold on;
errorbar(nanmean(u2ua_frac_m,1),nanstd(u2ua_frac_m,[],1)./sqrt(sum(~isnan(u2ua_frac_m))),'k','linewidth',2);
errorbar(nanmean(u2ua_frac_q,1),nanstd(u2ua_frac_q,[],1)./sqrt(sum(~isnan(u2ua_frac_q))),'r','linewidth',2);
errorbar(nanmean(u2ua_frac_ua,1),nanstd(u2ua_frac_ua,[],1)./sqrt(sum(~isnan(u2ua_frac_ua))),'b','linewidth',2);
errorbar(nanmean(u2ua_frac_s,1),nanstd(u2ua_frac_s,[],1)./sqrt(sum(~isnan(u2ua_frac_s))),'g','linewidth',2);
ylabel('Fraction of cells');
xlabel('Days');
title('Classification of ~UA2UA');
legend({'M','Q','UA','S'});

% Get classification between m/q classifications
[~,last_m] = cellfun(@(m,m2q) max(fliplr(m(m2q,:)),[],2),m_pad,m2q,'uni',false);
last_m = cellfun(@(x) 15-x,last_m,'uni',false);

%% Reviewer 1 point 1: fig 5D classification of switch cells by week

clearvars -except data analysis classified_rois

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) ~m & ~q & ~ua,m_pad,q_pad,ua_pad,'uni',false);

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;

all_cells = cellfun(@(m) true(size(m,1),1),m_pad,'uni',false);
m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);
q2m = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);

m2m = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);
q2q = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);

m2u = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) &...
    sum(m(:,8:14),2) < min_classdays,m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2) & ...
    sum(m(:,1:7),2) < min_classdays,m_pad,q_pad,'uni',false);

q2u = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) < min_classdays,m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2) & ...
    sum(q(:,1:7),2) < min_classdays,m_pad,q_pad,'uni',false);

ua2u = cellfun(@(ua,m,q) ...
    sum(ua(:,1:7),2) >= sum(m(:,1:7),2) & ...
    sum(ua(:,1:7),2) >= sum(q(:,1:7),2) & ...
    sum(ua(:,8:end),2) < min_classdays,ua_pad,m_pad,q_pad,'uni',false);
u2ua = cellfun(@(ua,m,q) ...
    sum(ua(:,8:14),2) >= sum(m(:,8:14),2) & ...
    sum(ua(:,8:14),2) >= sum(q(:,8:14),2) & ...
    sum(ua(:,1:7),2) < min_classdays,ua_pad,m_pad,q_pad,'uni',false);

plot_classes = {m2m,q2q,m2q,q2m,m2u,u2m,q2u,u2q,ua2u,u2ua};
class_frac = cell(size(plot_classes));
for curr_class_idx = 1:length(plot_classes);
    
    curr_class = plot_classes{curr_class_idx};
    
    class_frac{curr_class_idx}{1,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),s_pad,curr_class,'uni',false)');

    class_frac{curr_class_idx}{1,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),s_pad,curr_class,'uni',false)');

end

class_frac_mean = cellfun(@(x) cellfun(@(x) nanmean(x,1),x)',class_frac,'uni',false);
class_frac_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(x,[],1)./sqrt(sum(~isnan(x),1)),x)',class_frac,'uni',false);

figure; hold on;
bar(vertcat(class_frac_mean{:}),'stacked')
errorbar(cumsum(vertcat(class_frac_mean{:}),2), ...
    vertcat(class_frac_sem{:}),'.k','linewidth',2);
set(gca,'XTick',1:20);
set(gca,'XTickLabel', ...
    {'M2M 1','M2M 2','Q2Q 1','Q2Q 2', ...
    'M2Q 1','M2Q 2','Q2M 1','Q2M 2', ...
    'M2U 1','M2U 2','U2M 1','U2M 2', ...
    'Q2U 1','Q2U 2','U2Q 1','U2Q 2', ...
    'UA2U 1','UA2U 2','U2UA 1','U2UA 2'});
legend({'M','Q','UA','S'});
ylabel('Frac days in category');



%% Reviewer 1 point 1: fig 5 cumulative category switch (within day)

clearvars -except data analysis classified_rois

% Plot cagetorization of m2q cells

% Get categories
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) ~m & ~q & ~ua,m_pad,q_pad,ua_pad,'uni',false);

plot_classes = {m_pad,q_pad,ua_pad,s_pad};
class_frac = cell(size(plot_classes));
for curr_class_idx = 1:length(plot_classes);
    
    curr_class = plot_classes{curr_class_idx};
    
    class_frac{curr_class_idx}{1} = cell2mat(cellfun(@(cl_alt,cl) arrayfun(@(day) ...
        nanmean(nanmean(cl_alt(cl(:,day)==1,setdiff(1:14,day)),2),1),1:14),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2} = cell2mat(cellfun(@(cl_alt,cl) arrayfun(@(day) ...
        nanmean(nanmean(cl_alt(cl(:,day)==1,setdiff(1:14,day)),2),1),1:14),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3} = cell2mat(cellfun(@(cl_alt,cl) arrayfun(@(day) ...
        nanmean(nanmean(cl_alt(cl(:,day)==1,setdiff(1:14,day)),2),1),1:14),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4} = cell2mat(cellfun(@(cl_alt,cl) arrayfun(@(day) ...
        nanmean(nanmean(cl_alt(cl(:,day)==1,setdiff(1:14,day)),2),1),1:14),s_pad,curr_class,'uni',false)');
    
end

class_frac_mean = cellfun(@(x) cellfun(@(x) nanmean(x,1),x,'uni',false),class_frac,'uni',false);
class_frac_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(x,[],1)./sqrt(sum(~isnan(x),1)),x,'uni',false),class_frac,'uni',false);

figure; 
for curr_class_idx = 1:length(plot_classes)
    subplot(2,2,curr_class_idx); hold on;
    area(vertcat(class_frac_mean{curr_class_idx}{:})');
    errorbar(cumsum(vertcat(class_frac_mean{curr_class_idx}{:}),1)', ...
        vertcat(class_frac_sem{curr_class_idx}{:})','.k');
end
legend({'M','Q','UA','S'})

%% Reviewer 1 point 1: fig 5 cumulative category switch (within week)

clearvars -except data analysis classified_rois

% Plot cagetorization of m2q cells

% Get categories
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) ~m & ~q & ~ua,m_pad,q_pad,ua_pad,'uni',false);

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;
m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) >= min_classdays & ...
    sum(q(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
m2u = cellfun(@(m,q) ...
    sum(m(:,1:7),2) >= min_classdays & ...
    sum(m(:,8:end),2) < min_classdays,m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    sum(m(:,1:7),2) < min_classdays & ...
    sum(m(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
q2u = cellfun(@(m,q) ...
    sum(q(:,1:7),2) >= min_classdays & ...
    sum(q(:,8:end),2) < min_classdays,m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    sum(q(:,1:7),2) < min_classdays & ...
    sum(q(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);
ua2u = cellfun(@(ua) ...
    sum(ua(:,1:7),2) >= min_classdays & ...
    sum(ua(:,8:end),2) < min_classdays,ua_pad,'uni',false);
u2ua = cellfun(@(ua) ...
    sum(ua(:,1:7),2) < min_classdays & ...
    sum(ua(:,8:end),2) >= min_classdays,ua_pad,'uni',false);

plot_classes = {m2u,u2m,q2u,u2q,ua2u,u2ua};
class_frac = cell(size(plot_classes));
for curr_class_idx = 1:length(plot_classes);
    
    curr_class = plot_classes{curr_class_idx};
    
    class_frac{curr_class_idx}{1} = cell2mat(cellfun(@(cl,sw) nanmean(cl(sw,:),1),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2} = cell2mat(cellfun(@(cl,sw) nanmean(cl(sw,:),1),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3} = cell2mat(cellfun(@(cl,sw) nanmean(cl(sw,:),1),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4} = cell2mat(cellfun(@(cl,sw) nanmean(cl(sw,:),1),s_pad,curr_class,'uni',false)');

end

class_frac_mean = cellfun(@(x) cellfun(@(x) nanmean(x,1),x,'uni',false),class_frac,'uni',false);
class_frac_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(x,[],1)./sqrt(sum(~isnan(x),1)),x,'uni',false),class_frac,'uni',false);

figure; 
for curr_class_idx = 1:length(plot_classes)
    subplot(3,2,curr_class_idx); hold on;
    area(vertcat(class_frac_mean{curr_class_idx}{:})');
    errorbar(cumsum(vertcat(class_frac_mean{curr_class_idx}{:}),1)', ...
        vertcat(class_frac_sem{curr_class_idx}{:})','.k');   
end
legend({'M','Q','UA','S'})

%% Reviewer 1 point 1: 5 movement-aligned activity by category

clearvars -except data analysis classified_rois

%%%%% Define categories of cells

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) ~m & ~q & ~ua,m_pad,q_pad,ua_pad,'uni',false);

all_cells = cellfun(@(m) true(size(m,1),1),m_pad,'uni',false);

m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);
q2m = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);

m2m = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2),m_pad,q_pad,'uni',false);
q2q = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2),m_pad,q_pad,'uni',false);

m2u = cellfun(@(m,q) ...
    sum(m(:,1:7),2) > sum(q(:,1:7),2) &...
    sum(m(:,8:14),2) <= sum(q(:,8:14),2),m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    sum(m(:,8:14),2) > sum(q(:,8:14),2) & ...
    sum(m(:,1:7),2) <= sum(q(:,1:7),2),m_pad,q_pad,'uni',false);

q2u = cellfun(@(m,q) ...
    sum(q(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(q(:,8:14),2) <= sum(m(:,8:14),2),m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    sum(q(:,8:14),2) > sum(m(:,8:14),2) & ...
    sum(q(:,1:7),2) < sum(m(:,1:7),2),m_pad,q_pad,'uni',false);

ua2u = cellfun(@(ua,m,q) ...
    sum(ua(:,1:7),2) > sum(m(:,1:7),2) & ...
    sum(ua(:,1:7),2) > sum(q(:,1:7),2) & ...
    sum(ua(:,8:end),2) < sum(m(:,8:end),2) & ...
    sum(ua(:,8:end),2) < sum(q(:,8:end),2),ua_pad,m_pad,q_pad,'uni',false);
u2ua = cellfun(@(ua,m,q) ...
    sum(ua(:,8:14),2) > sum(m(:,8:14),2) & ...
    sum(ua(:,8:14),2) > sum(q(:,8:14),2) & ...
    sum(ua(:,1:7),2) < sum(m(:,1:7),2) & ...
    sum(ua(:,1:7),2) < sum(q(:,1:7),2),ua_pad,m_pad,q_pad,'uni',false);

%%%%% Plot average movement-aligned activity

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
        nonmove_frames = 28*3;
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

plot_classes = {m2m,q2q,m2u,u2m,q2u,u2q,m2q,q2m};
class_labels = {'M2M','Q2Q','M2~M','~M2M','Q2~Q','~Q2Q','M2Q','Q2M'};

move_act_onset_mean = cell(length(plot_classes),1);
move_act_offset_mean = cell(length(plot_classes),1);
for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        curr_onset = move_act_onset_all{curr_animal,curr_day};
        curr_offset = move_act_offset_all{curr_animal,curr_day};
        
        % For control animal
        if isempty(curr_onset)
            continue
        end

        for curr_class_idx = 1:length(plot_classes)
            curr_class = plot_classes{curr_class_idx}{curr_animal};
            move_act_onset_mean{curr_class_idx}{curr_animal,curr_day} = ...
                nanmean(nanmean(curr_onset(:,:,curr_class),1),3);              
            move_act_offset_mean{curr_class_idx}{curr_animal,curr_day} = ...
                nanmean(nanmean(curr_offset(:,:,curr_class),1),3);              
        end   
        
    end
end

group_days = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
figure;
for curr_class_idx = 1:length(plot_classes)
    curr_onset_mean = cell2mat(cellfun(@(x) ...
        nanmean(vertcat(move_act_onset_mean{curr_class_idx}{:,x}),1),group_days,'uni',false)');
    curr_offset_mean = cell2mat(cellfun(@(x) ...
        nanmean(vertcat(move_act_offset_mean{curr_class_idx}{:,x}),1),group_days,'uni',false)');
    subplot(4,2,curr_class_idx); hold on;
    set(gca,'ColorOrder',copper(length(group_days)));
    plot([curr_onset_mean,nan(length(group_days),5),curr_offset_mean]');
    title(class_labels{curr_class_idx});
end
legend(cellfun(@num2str,group_days,'uni',false));


%% Reviewer 1 point 1: fig 5 final: movement/day/week/total timescales

%%%%% Define categories of cells

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

% For the one control animal
m = cellfun(@double,m,'uni',false);
q = cellfun(@double,q,'uni',false);
ua = cellfun(@double,ua,'uni',false);

m_pad = cellfun(@(x) padarray(x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);
s_pad = cellfun(@(m,q,ua) (m==0) & (q==0) & (ua==0),m_pad,q_pad,ua_pad,'uni',false);

all_cells = cellfun(@(m) true(size(m,1),1),m_pad,'uni',false);

m2q = cellfun(@(m,q) ...
    nansum(m(:,1:7),2) > nansum(q(:,1:7),2) & ...
    nansum(q(:,8:14),2) > nansum(m(:,8:14),2),m_pad,q_pad,'uni',false);
q2m = cellfun(@(m,q) ...
    nansum(q(:,1:7),2) > nansum(m(:,1:7),2) & ...
    nansum(m(:,8:14),2) > nansum(q(:,8:14),2),m_pad,q_pad,'uni',false);

m2m = cellfun(@(m,q) ...
    nansum(m(:,1:7),2) > nansum(q(:,1:7),2) & ...
    nansum(m(:,8:14),2) > nansum(q(:,8:14),2),m_pad,q_pad,'uni',false);
q2q = cellfun(@(m,q) ...
    nansum(q(:,1:7),2) > nansum(m(:,1:7),2) & ...
    nansum(q(:,8:14),2) > nansum(m(:,8:14),2),m_pad,q_pad,'uni',false);

m2u = cellfun(@(m,q) ...
    nansum(m(:,1:7),2) > nansum(q(:,1:7),2) &...
    nansum(m(:,8:14),2) <= nansum(q(:,8:14),2),m_pad,q_pad,'uni',false);
u2m = cellfun(@(m,q) ...
    nansum(m(:,8:14),2) > nansum(q(:,8:14),2) & ...
    nansum(m(:,1:7),2) <= nansum(q(:,1:7),2),m_pad,q_pad,'uni',false);

q2u = cellfun(@(m,q) ...
    nansum(q(:,1:7),2) > nansum(m(:,1:7),2) & ...
    nansum(q(:,8:14),2) <= nansum(m(:,8:14),2),m_pad,q_pad,'uni',false);
u2q = cellfun(@(m,q) ...
    nansum(q(:,8:14),2) > nansum(m(:,8:14),2) & ...
    nansum(q(:,1:7),2) <= nansum(m(:,1:7),2),m_pad,q_pad,'uni',false);

u2u = cellfun(@(m,q,s) ...
    nansum(m(:,1:7),2) == nansum(q(:,1:7),2) & ...
    nansum(m(:,8:14),2) == nansum(q(:,8:14),2) & ...
    ~all(s,2),m_pad,q_pad,s_pad,'uni',false);

s2s = cellfun(@(s) ...
    all(s,2),s_pad,'uni',false);

ua2u = cellfun(@(ua,m,q) ...
    nansum(ua(:,1:7),2) > nansum(m(:,1:7),2) & ...
    nansum(ua(:,1:7),2) > nansum(q(:,1:7),2) & ...
    nansum(ua(:,8:end),2) < nansum(m(:,8:end),2) & ...
    nansum(ua(:,8:end),2) < nansum(q(:,8:end),2),ua_pad,m_pad,q_pad,'uni',false);
u2ua = cellfun(@(ua,m,q) ...
    nansum(ua(:,8:14),2) > nansum(m(:,8:14),2) & ...
    nansum(ua(:,8:14),2) > nansum(q(:,8:14),2) & ...
    nansum(ua(:,1:7),2) < nansum(m(:,1:7),2) & ...
    nansum(ua(:,1:7),2) < nansum(q(:,1:7),2),ua_pad,m_pad,q_pad,'uni',false);

%%%%% Plot average movement-aligned activity

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
        nonmove_frames = 28*3;
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

plot_classes = {m2m,q2q,m2u,u2m,q2u,u2q,all_cells};
class_labels = {'M2M','Q2Q','M2~M','~M2M','Q2~Q','~Q2Q','All'};

move_act_onset_mean = cell(length(plot_classes),1);
move_act_offset_mean = cell(length(plot_classes),1);
for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        curr_onset = move_act_onset_all{curr_animal,curr_day};
        curr_offset = move_act_offset_all{curr_animal,curr_day};
        
        % For control animal
        if isempty(curr_onset)
            continue
        end

        for curr_class_idx = 1:length(plot_classes)
            curr_class = plot_classes{curr_class_idx}{curr_animal};
            move_act_onset_mean{curr_class_idx}{curr_animal,curr_day} = ...
                nanmean(nanmean(curr_onset(:,:,curr_class),1),3);              
            move_act_offset_mean{curr_class_idx}{curr_animal,curr_day} = ...
                nanmean(nanmean(curr_offset(:,:,curr_class),1),3);              
        end   
        
    end
end

group_days = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
figure;
for curr_class_idx = 1:length(plot_classes)
    curr_onset_mean = cell2mat(cellfun(@(x) ...
        nanmean(vertcat(move_act_onset_mean{curr_class_idx}{:,x}),1),group_days,'uni',false)');
    curr_offset_mean = cell2mat(cellfun(@(x) ...
        nanmean(vertcat(move_act_offset_mean{curr_class_idx}{:,x}),1),group_days,'uni',false)');
    subplot(4,2,curr_class_idx); hold on;
    set(gca,'ColorOrder',copper(length(group_days)));
    plot([curr_onset_mean,nan(length(group_days),5),curr_offset_mean]');
    title(class_labels{curr_class_idx});
end
legend(cellfun(@num2str,group_days,'uni',false));


%%%%% Plot activity by day

% Get activity of all cells during m/q/total 
avg_movement_act = cell(length(data),1);
avg_quiescent_act = cell(length(data),1);
avg_total_act = cell(length(data),1);

for curr_animal = 1:length(data)
        
    use_cells = true(size(classified_rois(curr_animal).movement,1),1);
        
    curr_avg_movement_act = nan(length(use_cells),14);
    curr_avg_quiescent_act = nan(length(use_cells),14);
    curr_avg_total_act = nan(length(use_cells),14);
    
     for i = 1:14
         
         if length(data(curr_animal).im) < i || ...
                 isempty(data(curr_animal).im(i).roi_trace_thresh)
             continue
         end
         
        curr_avg_movement_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames > 0),2);
        
        curr_avg_quiescent_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames == 0),2);
        
        curr_avg_total_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,:),2);
        
     end
    
    avg_movement_act{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act{curr_animal} = curr_avg_quiescent_act;
    avg_total_act{curr_animal} = curr_avg_total_act;
    
    disp(curr_animal);
    
end

plot_classes = {m2m,q2q,m2u,u2m,q2u,u2q,all_cells};
class_labels = {'M2M','Q2Q','M2~M','~M2M','Q2~Q','~Q2Q','All'};

avg_act = cell(length(plot_classes),1);
for curr_class_idx = 1:length(plot_classes)
    
    curr_class = plot_classes{curr_class_idx};
    
    avg_act{curr_class_idx}{1} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_movement_act,curr_class','uni',false);
    avg_act{curr_class_idx}{2} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_quiescent_act,curr_class','uni',false);
    avg_act{curr_class_idx}{3} = cellfun(@(act,class) nanmean(act(class,:),1), ...
        avg_total_act,curr_class','uni',false);
    
end

avg_act_mean = cellfun(@(x) cellfun(@(x) ...
    nanmean(vertcat(x{:}),1),x,'uni',false),avg_act,'uni',false);
avg_act_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(vertcat(x{:}),[],1)./sqrt(sum(~isnan(vertcat(x{:})),1)), ...
    x,'uni',false),avg_act,'uni',false);

figure;
for curr_class_idx = 1:length(plot_classes)
    subplot(4,2,curr_class_idx); hold on;
    errorbar(avg_act_mean{curr_class_idx}{1},avg_act_sem{curr_class_idx}{1},'r');
    errorbar(avg_act_mean{curr_class_idx}{2},avg_act_sem{curr_class_idx}{2},'b');
    errorbar(avg_act_mean{curr_class_idx}{3},avg_act_sem{curr_class_idx}{3},'k');
    
    xlabel('day')
    ylabel('average df/f')
    title(class_labels{curr_class_idx})
    legend({'during movement','during quiescence','total'})
end

% Stats (IN PROGRESS...)
% all: m stable week 1
curr_data = zscore(vertcat(avg_act{end}{1}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,1:7));
% all: m down week 2
curr_data = zscore(vertcat(avg_act{end}{1}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,8:14));
% all: q stable week 1
curr_data = zscore(vertcat(avg_act{end}{2}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,1:7));
% all: q down week 2
curr_data = zscore(vertcat(avg_act{end}{2}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,8:14));
% m2m: up week 1
curr_data = zscore(vertcat(avg_act{1}{1}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,1:7));
% m2m: down week 2
curr_data = zscore(vertcat(avg_act{1}{1}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,8:14));
% q2q: up week 1
curr_data = zscore(vertcat(avg_act{2}{2}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,1:7));
% q2q: down week 2
curr_data = zscore(vertcat(avg_act{2}{2}{:}),[],2);
[r,p] = corrcoef(repmat(1:7,size(curr_data,1),1),curr_data(:,8:14));

%%%%% Plot classification by week

plot_classes = {m2m,q2q,m2u,u2m,q2u,u2q,all_cells};
class_frac = cell(size(plot_classes));
for curr_class_idx = 1:length(plot_classes);
    
    curr_class = plot_classes{curr_class_idx};
    
    class_frac{curr_class_idx}{1,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4,1} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,1:7),1),2),s_pad,curr_class,'uni',false)');

    class_frac{curr_class_idx}{1,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),m_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{2,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),q_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{3,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),ua_pad,curr_class,'uni',false)');
    class_frac{curr_class_idx}{4,2} = cell2mat(cellfun(@(cl,sw) nanmean(nanmean(cl(sw,8:14),1),2),s_pad,curr_class,'uni',false)');

end

class_frac_mean = cellfun(@(x) cellfun(@(x) nanmean(x,1),x)',class_frac,'uni',false);
class_frac_sem = cellfun(@(x) cellfun(@(x) ...
    nanstd(x,[],1)./sqrt(sum(~isnan(x),1)),x)',class_frac,'uni',false);

figure; hold on;
bar(vertcat(class_frac_mean{:}),'stacked')
errorbar(cumsum(vertcat(class_frac_mean{:}),2), ...
    vertcat(class_frac_sem{:}),'.k','linewidth',2);
set(gca,'XTick',1:16);
set(gca,'XTickLabel', ...
    {'M2M 1','M2M 2','Q2Q 1','Q2Q 2', ...
    'M2U 1','M2U 2','U2M 1','U2M 2', ...
    'Q2U 1','Q2U 2','U2Q 1','U2Q 2', ...
    'All 1','All 2'});
legend({'M','Q','UA','S'});
ylabel('Frac days in category');


%%%%% Plot fraction of cells per category
plot_classes = {m2m,q2q,u2u,s2s,q2u,q2m,u2m,m2u,m2q,u2q,};
class_labels = {'M2M','Q2Q','U2U','S2S','Q2~Q','Q2M','~M2M','M2~M','M2Q','~Q2Q'};

num_cells = cellfun(@length,all_cells);
class_num = cell2mat(cellfun(@(x) cellfun(@sum,x),plot_classes,'uni',false)');
class_frac_mean = nanmean(bsxfun(@rdivide,class_num,num_cells),2);

m2q_category_overlap = mean(cellfun(@(m2u,u2q) sum(m2u&u2q),m2u,u2q)./num_cells);
q2m_category_overlap = mean(cellfun(@(q2u,u2m) sum(q2u&u2m),q2u,u2m)./num_cells);

figure;
pie(class_frac_mean);
legend(class_labels);


%% Reviewer 1 point 1: fig 5c bar plot change
% didn't ask for this, but to make it more intuitive since grids gone
% changed my mind: just change the ylabel

clearvars -except data analysis classified_rois

% Same classification overlap
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};

m2m = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + repmat(sum(m,1),size(m,2),1)),m,'uni',false);
m2m_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),m2m,'uni',false);
m2m_cat = cat(3,m2m_pad{:});

q2q = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + repmat(sum(q,1),size(q,2),1)),q,'uni',false);
q2q_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),q2q,'uni',false);
q2q_cat = cat(3,q2q_pad{:});

% Same classification overlap shuffle
n_rep = 1000;

m2m_cat_shuff = nan(14,14,length(data),n_rep);
q2q_cat_shuff = nan(14,14,length(data),n_rep);

m_use = cellfun(@(x) x(any(x,2),:),m,'uni',false);
q_use = cellfun(@(x) x(any(x,2),:),q,'uni',false);

for i = 1:n_rep
    m_shuff = cellfun(@(x) shake(x,1),m_use,'uni',false);    
    m2m_shuff = cellfun(@(m) (+m'*+m)./(repmat(sum(m,1)',1,size(m,2)) + ...
        repmat(sum(m,1),size(m,2),1)),m_shuff,'uni',false);  
    m2m_shuff_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),m2m_shuff,'uni',false);
    m2m_cat_shuff(:,:,:,i) = cat(3,m2m_shuff_pad{:});
    
    q_shuff = cellfun(@(x) shake(x,1),q_use,'uni',false);    
    q2q_shuff = cellfun(@(q) (+q'*+q)./(repmat(sum(q,1)',1,size(q,2)) + ...
        repmat(sum(q,1),size(q,2),1)),q_shuff,'uni',false);   
    q2q_shuff_pad = cellfun(@(x) padarray(x(1:min(length(x),14),1:min(length(x),14)), ...
    [max(14-length(x),0),max(14-length(x),0)],NaN,'post'),q2q_shuff,'uni',false);
    q2q_cat_shuff(:,:,:,i) = cat(3,q2q_shuff_pad{:});
 
    disp(i);
end

% Real and shuffled overlap within weeks
m2m_overlap = cell2mat(arrayfun(@(x) [nanmean(AP_itril(m2m_cat(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(m2m_cat(8:end,8:end,x),-1))], ...
    transpose(1:size(m2m_cat,3)),'uni',false));
q2q_overlap = cell2mat(arrayfun(@(x) [nanmean(AP_itril(q2q_cat(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(q2q_cat(8:end,8:end,x),-1))], ...
    transpose(1:size(q2q_cat,3)),'uni',false));

m2m_cat_shuff_mean = nanmean(m2m_cat_shuff,4);
m2m_overlap_shuff = cell2mat(arrayfun(@(x) [nanmean(AP_itril(m2m_cat_shuff_mean(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(m2m_cat_shuff_mean(8:14,8:14,x),-1))],transpose(1:size(m2m_cat_shuff_mean,3)),'uni',false));
m2m_overlap_shuff_ci = cell2mat(arrayfun(@(x) [prctile(AP_itril(m2m_cat_shuff_mean(1:7,1:7,x),-1),[2.5,97.5]), ...
    prctile(AP_itril(m2m_cat_shuff_mean(8:14,8:14,x),-1),[2.5,97.5])],transpose(1:size(m2m_cat_shuff_mean,3)),'uni',false));
q2q_cat_shuff_mean = nanmean(q2q_cat_shuff,4);
q2q_overlap_shuff = cell2mat(arrayfun(@(x) [nanmean(AP_itril(q2q_cat_shuff_mean(1:7,1:7,x),-1)), ...
    nanmean(AP_itril(q2q_cat_shuff_mean(8:14,8:14,x),-1))],transpose(1:size(q2q_cat_shuff_mean,3)),'uni',false));
q2q_overlap_shuff_ci = cell2mat(arrayfun(@(x) [prctile(AP_itril(q2q_cat_shuff_mean(1:7,1:7,x),-1),[2.5,97.5]), ...
    prctile(AP_itril(q2q_cat_shuff_mean(8:14,8:14,x),-1),[2.5,97.5])],transpose(1:size(q2q_cat_shuff_mean,3)),'uni',false));

figure; hold on;
bar([nanmean(m2m_overlap),nanmean(q2q_overlap)],'FaceColor','k');
plot([nanmean(m2m_overlap_shuff),nanmean(q2q_overlap_shuff)],'.r');
plot([reshape(nanmean(m2m_overlap_shuff_ci),2,2),reshape(nanmean(q2q_overlap_shuff_ci),2,2)]','--r');


m2m_overlap_norm = m2m_overlap./m2m_overlap_shuff;
q2q_overlap_norm = q2q_overlap./q2q_overlap_shuff;
figure; hold on;
bar([nanmean(m2m_overlap_norm),nanmean(q2q_overlap_norm)],'FaceColor','k');



%% Reviewer 1 point 1: 4A of all cells total and by classification

clearvars -except data analysis classified_rois

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
        nonmove_frames = 28*3;
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
move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_all,'uni',false);

% Activity of all animals and cells, day averaged (like 4A, grand average)

move_act_onset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_onset_all,'uni',false);
move_act_offset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_offset_all,'uni',false);

curr_data_onset_norm_mean = cell(size(data));
curr_data_offset_norm_mean = cell(size(data));
for curr_animal = 1:length(data)
    
    use_days = cellfun(@(x) ~isempty(x),move_act_onset_all_mean(curr_animal,:));
    
    curr_data_onset = cat(3,move_act_onset_all_mean{curr_animal,use_days});
    curr_data_offset = cat(3,move_act_offset_all_mean{curr_animal,use_days});
    
    curr_data_onset_minsub = bsxfun(@minus,curr_data_onset,min(curr_data_onset,[],2));
    curr_data_onset_norm = bsxfun(@times,curr_data_onset_minsub,1./max(curr_data_onset_minsub,[],2));
    
    curr_data_offset_minsub = bsxfun(@minus,curr_data_offset,min(curr_data_offset,[],2));
    curr_data_offset_norm = bsxfun(@times,curr_data_offset_minsub,1./max(curr_data_offset_minsub,[],2));
    
    curr_data_onset_norm_mean{curr_animal} = nanmean(curr_data_onset_norm,3);
    curr_data_offset_norm_mean{curr_animal} = nanmean(curr_data_offset_norm,3); 
    
%     % or to just normalize at the end
%     curr_data_onset_norm_mean{curr_animal} = nanmean(curr_data_onset,3);
%     curr_data_offset_norm_mean{curr_animal} = nanmean(curr_data_offset,3); 
    
end

all_onset = cell2mat(curr_data_onset_norm_mean);
all_offset = cell2mat(curr_data_offset_norm_mean);

all_combined = [all_onset,all_offset];
% m_frames = false(size(all_combined,2),1);
% m_frames(nonmove_frames+1:nonmove_frames+1+move_frames*2) = true;
% m_act = nanmean(all_combined(:,m_frames),2);
% q_act = nanmean(all_combined(:,~m_frames),2);
% mq_ratio = (m_act - q_act)./(q_act + m_act);
% [~,sort_idx] = sort(mq_ratio);
all_combined(isnan(all_combined)) = 0;
[coeff,score,latent] = princomp(all_combined');
[~,sort_idx] = sort(coeff(:,1));
% put empty rows on the bottom of the sort
empty_rows = find(~any(all_combined,2));
sort_idx = [sort_idx(~ismember(sort_idx,empty_rows));sort_idx(ismember(sort_idx,empty_rows))];
figure; colormap(hot);
subplot(1,2,1); imagesc(all_onset(sort_idx,:));
ylabel('Cell number');
title('Movement onset');
subplot(1,2,2); imagesc(all_offset(sort_idx,:));
ylabel('Cell number');
title('Movement offset');

% As above like 4A, but separating between days classified
curr_data_q_onset_norm_mean = cell(size(data));
curr_data_m_onset_norm_mean = cell(size(data));
curr_data_ua_onset_norm_mean = cell(size(data));

curr_data_q_offset_norm_mean = cell(size(data));
curr_data_m_offset_norm_mean = cell(size(data));
curr_data_ua_offset_norm_mean = cell(size(data));
for curr_animal = 1:length(data)
    
    use_days = cellfun(@(x) ~isempty(x),move_act_onset_all_mean(curr_animal,:));
    
    curr_data_onset = cat(3,move_act_onset_all_mean{curr_animal,use_days});
    curr_data_offset = cat(3,move_act_offset_all_mean{curr_animal,use_days});
    
    curr_data_onset_minsub = bsxfun(@minus,curr_data_onset,min(curr_data_onset,[],2));
    curr_data_onset_norm = bsxfun(@times,curr_data_onset_minsub,1./max(curr_data_onset_minsub,[],2));
    
    curr_data_offset_minsub = bsxfun(@minus,curr_data_offset,min(curr_data_offset,[],2));
    curr_data_offset_norm = bsxfun(@times,curr_data_offset_minsub,1./max(curr_data_offset_minsub,[],2));
    
    %     % or to use unnormalized
    %     curr_data_onset_norm = curr_data_onset;
    %     curr_data_offset_norm = curr_data_offset;
    
    % This is for a L2/3 animal
    use_days = min(14,size(classified_rois(curr_animal).movement,2));
    
    q = +classified_rois(curr_animal).quiescent(:,use_days);
    q(~q) = NaN;
    m = +classified_rois(curr_animal).movement(:,use_days);
    m(~m) = NaN;
    ua = +classified_rois(curr_animal).unclassified_active(:,use_days);
    ua(~ua) = NaN;
    
    curr_data_q_onset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_onset_norm, ...
        permute(q,[1,3,2])),3);
    curr_data_m_onset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_onset_norm, ...
        permute(m,[1,3,2])),3);
    curr_data_ua_onset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_onset_norm, ...
        permute(ua,[1,3,2])),3);
    
    curr_data_q_offset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_offset_norm, ...
        permute(q,[1,3,2])),3);
    curr_data_m_offset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_offset_norm, ...
        permute(m,[1,3,2])),3);
    curr_data_ua_offset_norm_mean{curr_animal} = ...
        nanmean(bsxfun(@times,curr_data_offset_norm, ...
        permute(ua,[1,3,2])),3);
    
end
all_q_onset = cell2mat(curr_data_q_onset_norm_mean);
all_m_onset = cell2mat(curr_data_m_onset_norm_mean);
all_ua_onset = cell2mat(curr_data_ua_onset_norm_mean);

all_q_offset = cell2mat(curr_data_q_offset_norm_mean);
all_m_offset = cell2mat(curr_data_m_offset_norm_mean);
all_ua_offset = cell2mat(curr_data_ua_offset_norm_mean);

% only plot filled rows
figure; colormap(hot);
curr_data = all_q_onset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,1); imagesc(curr_data);
title('Q onset');
curr_data = all_q_offset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,2); imagesc(curr_data);
title('Q offset');
curr_data = all_m_onset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,3); imagesc(curr_data);
title('M onset');
curr_data = all_m_offset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,4); imagesc(curr_data);
title('M offset');
curr_data = all_ua_onset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,5); imagesc(curr_data);
title('UA onset');
curr_data = all_ua_offset(sort_idx,:);
curr_data(~any(curr_data,2),:) = [];
subplot(1,6,6); imagesc(curr_data);
title('UA offset');



%% Reviewer 1 point 4: fig 4B, fraction active cells

clearvars -except data analysis classified_rois

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
        nonmove_frames = 28*3;
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

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed > 0;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed > 0;

        
    end
    disp(curr_animal);
end

move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_all,'uni',false);
move_act_m_onset_mean = cell(size(move_act_onset_mean));
move_act_q_onset_mean = cell(size(move_act_onset_mean));

move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_all,'uni',false);
move_act_m_offset_mean = cell(size(move_act_offset_mean));
move_act_q_offset_mean = cell(size(move_act_offset_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        move_act_m_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);       
        
    end
end

% Get mean activity across days for each animal
move_act_mean_onset_animalavg = nan(length(data),total_frames);
move_act_mean_m_onset_animalavg = nan(length(data),total_frames);
move_act_mean_q_onset_animalavg = nan(length(data),total_frames);

move_act_mean_offset_animalavg = nan(length(data),total_frames);
move_act_mean_m_offset_animalavg = nan(length(data),total_frames);
move_act_mean_q_offset_animalavg = nan(length(data),total_frames);

for curr_animal = 1:length(data)
    
    move_act_mean_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_onset_mean{curr_animal,:}),1);
    move_act_mean_m_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_onset_mean{curr_animal,:}),1);
    move_act_mean_q_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_onset_mean{curr_animal,:}),1);
    
    move_act_mean_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_offset_mean{curr_animal,:}),1);
    move_act_mean_m_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_offset_mean{curr_animal,:}),1);
    move_act_mean_q_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_offset_mean{curr_animal,:}),1);
    
end

all_avg_onset = nanmean(move_act_mean_onset_animalavg,1);
m_avg_onset = nanmean(move_act_mean_m_onset_animalavg,1);
q_avg_onset = nanmean(move_act_mean_q_onset_animalavg,1);

all_sem_onset = nanstd(move_act_mean_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_onset_animalavg)));
m_sem_onset = nanstd(move_act_mean_m_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_m_onset_animalavg)));
q_sem_onset = nanstd(move_act_mean_q_onset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_q_onset_animalavg)));

all_avg_offset = nanmean(move_act_mean_offset_animalavg,1);
m_avg_offset = nanmean(move_act_mean_m_offset_animalavg,1);
q_avg_offset = nanmean(move_act_mean_q_offset_animalavg,1);

all_sem_offset = nanstd(move_act_mean_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_offset_animalavg)));
m_sem_offset = nanstd(move_act_mean_m_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_m_offset_animalavg)));
q_sem_offset = nanstd(move_act_mean_q_offset_animalavg,[],1)./sqrt(sum(~isnan(move_act_mean_q_offset_animalavg)));

% Plot all ROIs
df_scalebar = 0.02;
time_scalebar = 14.6; % frames

ymax = max([all_avg_onset+all_sem_onset,all_avg_offset+all_sem_offset]);

figure;
subplot(1,2,1); hold on;
fill([1:length(all_avg_onset),length(all_avg_onset):-1:1], ...
    [all_avg_onset+all_sem_onset,fliplr(all_avg_onset-all_sem_onset)],[0.5,0.5,0.5]);
plot(all_avg_onset,'k','linewidth',2);
ylim([0,ymax]);
line(repmat(nonmove_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Fraction of cells active');
title('Movement onset (all ROIs, all animals, all days)');

subplot(1,2,2); hold on;
fill([1:length(all_avg_offset),length(all_avg_offset):-1:1], ...
    [all_avg_offset+all_sem_offset,fliplr(all_avg_offset-all_sem_offset)],[0.5,0.5,0.5]);
plot(all_avg_offset,'k','linewidth',2);
ylim([0,ymax]);
line(repmat(move_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Fraction of cells active');
title('Movement onset (all ROIs, all animals, all days)');

line([0,0],[0,df_scalebar],'color','k','linewidth',2);
line([0,time_scalebar],[0,0],'color','k','linewidth',2);

% Plot classified ROIs
df_scalebar = 0.05;
time_scalebar = 14.6; % frames

ymax = max([m_avg_onset+m_sem_onset,m_avg_offset+m_sem_offset, ...
    q_avg_onset+q_sem_onset,q_avg_offset+q_sem_offset]);

figure;
subplot(1,2,1); hold on;
fill([1:length(m_avg_onset),length(m_avg_onset):-1:1], ...
    [m_avg_onset+m_sem_onset,fliplr(m_avg_onset-m_sem_onset)],[0.5,0.5,0.5]);
plot(m_avg_onset,'color',[0,0.8,0],'linewidth',2);
fill([1:length(q_avg_onset),length(q_avg_onset):-1:1], ...
    [q_avg_onset+q_sem_onset,fliplr(q_avg_onset-q_sem_onset)],[0.5,0.5,0.5]);
plot(q_avg_onset,'color',[0.8,0,0],'linewidth',2);
ylim([0,ymax]);
line(repmat(nonmove_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Fraction of cells active');
title('Movement onset (all ROIs, all animals, all days)');

subplot(1,2,2); hold on;
fill([1:length(m_avg_offset),length(m_avg_offset):-1:1], ...
    [m_avg_offset+m_sem_offset,fliplr(m_avg_offset-m_sem_offset)],[0.5,0.5,0.5]);
plot(m_avg_offset,'color',[0,0.8,0],'linewidth',2);
fill([1:length(q_avg_offset),length(q_avg_offset):-1:1], ...
    [q_avg_offset+q_sem_offset,fliplr(q_avg_offset-q_sem_offset)],[0.5,0.5,0.5]);
plot(q_avg_offset,'color',[0.8,0,0],'linewidth',2);
ylim([0,ymax]);
line(repmat(move_frames,2,1),ylim,'color','k','linestyle','--');
ylabel('Fraction of cells active');
title('Movement onset (all ROIs, all animals, all days)');

line([0,0],[0,df_scalebar],'color','k','linewidth',2);
line([0,time_scalebar],[0,0],'color','k','linewidth',2);

%% Reviewer 1 point 2: figure 3 lever trajectory examples

clearvars -except data analysis classified_rois

% Cued-rewarded trials

num_animals = length(data);

% Plot the trial-by-trial movement correlation
movement_start_time = 1001; % (defined in prepare_processed);
movement_use_time = 2000; % ms 

% (do in loop in case dispatcher file not saved, was in old cohort)
max_sessions = max(arrayfun(@(x) length(analysis(x).lever),1:num_animals));
movement_use_trials = cell(num_animals,1);
for curr_animal = 1:num_animals    
   for curr_session = 1:length(analysis(curr_animal).lever)
      if ~isempty(analysis(curr_animal).lever(curr_session).rewarded_movement)
          
          curr_move = ...
              horzcat(analysis(curr_animal).lever(curr_session).rewarded_movement_fixtime{ ...
              analysis(curr_animal).lever(curr_session).cued_movement_trials});
          
          movement_use_trials{curr_animal}{curr_session} = ...
              curr_move(1:movement_start_time+movement_use_time,:);          
      end
   end
end

% Histogram of trial correlation to template
template_days = 11:14;
template_move = cellfun(@(x) nanmean(horzcat(x{template_days}),2),movement_use_trials,'uni',false);

template_corr = cellfun(@(template,trials) ...
    cellfun(@(trials) 1-pdist2(template(movement_start_time:end)', ...
    trials(movement_start_time:end,:)','correlation')',trials,'uni',false), ...
    template_move,movement_use_trials,'uni',false);
template_corr = vertcat(template_corr{:});

use_days = {[1:4],[11:14]};
n_bins = 20;
corr_bins = linspace(-1,1,n_bins+1);
corr_bins_plot = corr_bins(1:end-1) + diff(corr_bins);
trial_template_corr_cat = nan(n_bins,length(use_days));
trial_template_corr_animal = nan(n_bins,length(use_days));
trial_template_corr_animal_sem = nan(n_bins,length(use_days));
for curr_use_days = 1:length(use_days)
    
    % all animals together
    curr_corr = vertcat(template_corr{:,use_days{curr_use_days}});
    corr_hist = histc(curr_corr,corr_bins)./length(curr_corr);
    trial_template_corr_cat(:,curr_use_days) = corr_hist(1:end-1);
    
    % mean across animals
    curr_template_corr = nan(n_bins,length(data));
    for curr_animal = 1:length(data)
        curr_corr = vertcat(template_corr{curr_animal,use_days{curr_use_days}});
        corr_hist = histc(curr_corr,corr_bins);
        curr_template_corr(:,curr_animal) = corr_hist(1:end-1)./length(curr_corr);
    end
    trial_template_corr_animal(:,curr_use_days) = nanmean(curr_template_corr,2);
    trial_template_corr_animal_sem(:,curr_use_days) = ...
        nanstd(curr_template_corr,[],2)./sqrt(sum(~isnan(curr_template_corr),2));
end

figure; 
subplot(2,1,1);
hold on; set(gca,'ColorOrder',copper(length(use_days)));
plot(corr_bins_plot,trial_template_corr_cat,'linewidth',2);
xlabel('Correlation with learned movement');
ylabel('Fraction of movements');
title('Concatenated across animals');
legend(cellfun(@num2str,use_days,'uni',false));
subplot(2,1,2);
hold on; set(gca,'ColorOrder',copper(length(use_days)));
for curr_plot = 1:length(use_days)
    fill([corr_bins_plot,fliplr(corr_bins_plot)], ...
        [trial_template_corr_animal(:,curr_plot)+trial_template_corr_animal_sem(:,curr_plot); ...
        flipud(trial_template_corr_animal(:,curr_plot)-trial_template_corr_animal_sem(:,curr_plot))], ...
        [0.5,0.5,0.5]);
end
plot(corr_bins_plot,trial_template_corr_animal,'linewidth',2);
xlabel('Correlation with learned movement');
ylabel('Fraction of movements');
title('Mean across animals');

% Example from mouse (by percentile?)
curr_animal = 3;
curr_corr = cellfun(@(x) vertcat(template_corr{curr_animal,x}),use_days,'uni',false);
n_prctile_bins = 5;
corr_prctile_bins = linspace(0,100,n_prctile_bins+1);
[~,corr_bin_idx] = cellfun(@(x) histc(x,prctile(x,corr_prctile_bins)),curr_corr,'uni',false);

figure;
for curr_days = 1:length(use_days)
    for curr_bin = 1:n_prctile_bins
        subplot(length(use_days),n_prctile_bins,n_prctile_bins*(curr_days-1)+curr_bin); 
        hold on; set(gca,'YDir','reverse');
        
        curr_move = horzcat(movement_use_trials{curr_animal}{use_days{curr_days}});
        curr_move_basesub = bsxfun(@minus,curr_move,nanmean(curr_move(1:movement_start_time-1,:),1));
        
        % random individual trials and mean
        n_plot_trials = 10;
        curr_bin_idx = find(corr_bin_idx{curr_days} == curr_bin);
        plot_trials = curr_bin_idx(randperm(length(curr_bin_idx),n_plot_trials));
        plot(curr_move_basesub(:,plot_trials),'color',[0.5,0.5,0.5]);
        plot(nanmean(curr_move_basesub(:,curr_bin_idx),2),'color','k','linewidth',2);
        line([movement_start_time,movement_start_time],ylim);

        if curr_days == 1
            title(['Percentile bin ',num2str(curr_bin)]);
        end
        if curr_bin == 1;
            ylabel(['Days ' num2str(use_days{curr_days})]);
        end
    end
end

figure;
plot(template_move{curr_animal},'k','linewidth',2);
line([movement_start_time,movement_start_time],ylim);
set(gca,'YDir','reverse');
title('Learned movement');


% % (OR WITH ALL MOVEMENTS - BUT THIS LOOKS BAD)
% 
% movement_start_time = 1001; % (defined in prepare_processed);
% movement_use_time = 2000; % ms 
% 
% % (do in loop in case dispatcher file not saved, was in old cohort)
% max_sessions = max(arrayfun(@(x) length(analysis(x).lever),1:num_animals));
% movement_use_trials = cell(num_animals,1);
% 
% for curr_animal = 1:num_animals
%     for curr_session = 1:length(analysis(curr_animal).lever)
%         
%         % Skip if animal doesn't have current session
%         if length(data(curr_animal).bhv) < curr_session || ...
%                 isempty(data(curr_animal).bhv(curr_session).lever_force)
%             continue
%         end
%         
%         % Get movement durations
%         [lever_active,lever_force_resample] = AP_parseLeverMovement_updated( ...
%             data(curr_animal).bhv(curr_session).lever_force);
%         
%         movement_starts = find(diff([0;lever_active;0]) == 1);
%         movement_stops = find(diff([0;lever_active;0]) == -1);
%         movement_durations = movement_stops - movement_starts;
%         movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
%         pre_movement_iti = [NaN;movement_iti];
%         post_movement_iti = [movement_iti;NaN];
%         
%         % Set criteria for pulling out movements        
%         min_move_time = 2000;
%         max_move_time = Inf;
%         min_pre_move_iti = 1000;
%         min_post_move_iti = 0;
%         use_movements = ...
%             movement_durations > min_move_time & ...
%             movement_durations < max_move_time & ...
%             pre_movement_iti > min_pre_move_iti & ...
%             post_movement_iti > min_post_move_iti;
%         
%         use_movement_start_times = movement_starts(use_movements);
%         
%         movements_idx = bsxfun(@plus,use_movement_start_times, ...
%             -movement_start_time+1:movement_use_time);
%          
%         timed_movements = lever_force_resample(movements_idx');
%         
%         movement_use_trials{curr_animal}{curr_session} = ...
%             timed_movements;
%                 
%     end
%     disp(curr_animal);
% end
% 
% % Histogram of trial correlation to template
% template_move = cellfun(@(x) nanmean(horzcat(x{11:14}),2),movement_use_trials,'uni',false);
% 
% template_corr = cellfun(@(template,trials) ...
%     cellfun(@(trials) 1-pdist2(template(movement_start_time:end)', ...
%     trials(movement_start_time:end,:)','correlation')',trials,'uni',false), ...
%     template_move,movement_use_trials,'uni',false);
% template_corr = vertcat(template_corr{:});
% 
% use_days = {[1:5],[6:10],[11:14]};
% n_bins = 30;
% corr_bins = linspace(-1,1,n_bins+1);
% corr_bins_plot = corr_bins(1:end-1) + diff(corr_bins);
% trial_template_corr_cat = nan(n_bins,length(use_days));
% trial_template_corr_animal = nan(n_bins,length(use_days));
% trial_template_corr_animal_sem = nan(n_bins,length(use_days));
% for curr_use_days = 1:length(use_days)
%     
%     % all animals together
%     curr_corr = vertcat(template_corr{:,use_days{curr_use_days}});
%     corr_hist = histc(curr_corr,corr_bins)./length(curr_corr);
%     trial_template_corr_cat(:,curr_use_days) = corr_hist(1:end-1);
%     
%     % mean across animals
%     curr_template_corr = nan(n_bins,length(data));
%     for curr_animal = 1:length(data)
%         curr_corr = vertcat(template_corr{curr_animal,use_days{curr_use_days}});
%         corr_hist = histc(curr_corr,corr_bins);
%         curr_template_corr(:,curr_animal) = corr_hist(1:end-1)./length(curr_corr);
%     end
%     trial_template_corr_animal(:,curr_use_days) = nanmean(curr_template_corr,2);
%     trial_template_corr_animal_sem(:,curr_use_days) = ...
%         nanstd(curr_template_corr,[],2)./sqrt(sum(~isnan(curr_template_corr),2));
% end
% 
% figure; 
% subplot(2,1,1);
% hold on; set(gca,'ColorOrder',copper(length(use_days)));
% plot(corr_bins_plot,trial_template_corr_cat,'linewidth',2);
% xlabel('Correlation with learned movement');
% ylabel('Fraction of movements');
% title('Concatenated across animals');
% legend(cellfun(@num2str,use_days,'uni',false));
% subplot(2,1,2);
% hold on; set(gca,'ColorOrder',copper(length(use_days)));
% for curr_plot = 1:length(use_days)
%     fill([corr_bins_plot,fliplr(corr_bins_plot)], ...
%         [trial_template_corr_animal(:,curr_plot)+trial_template_corr_animal_sem(:,curr_plot); ...
%         flipud(trial_template_corr_animal(:,curr_plot)-trial_template_corr_animal_sem(:,curr_plot))], ...
%         [0.5,0.5,0.5]);
% end
% plot(corr_bins_plot,trial_template_corr_animal,'linewidth',2);
% xlabel('Correlation with learned movement');
% ylabel('Fraction of movements');
% title('Mean across animals');
% 
% % Example from mouse (by percentile?)
% curr_animal = 3;
% curr_corr = cellfun(@(x) vertcat(template_corr{curr_animal,x}),use_days,'uni',false);
% n_prctile_bins = 5;
% corr_prctile_bins = linspace(0,100,n_prctile_bins+1);
% [~,corr_bin_idx] = cellfun(@(x) histc(x,prctile(x,corr_prctile_bins)),curr_corr,'uni',false);
% 
% figure;
% for curr_days = 1:length(use_days)
%     for curr_bin = 1:n_prctile_bins
%         subplot(length(use_days),n_prctile_bins,n_prctile_bins*(curr_days-1)+curr_bin); 
%         hold on; set(gca,'YDir','reverse');
%         
%         curr_move = horzcat(movement_use_trials{curr_animal}{use_days{curr_days}});
%         curr_move_basesub = bsxfun(@minus,curr_move,nanmean(curr_move(1:movement_start_time-1,:),1));
%         
%         % random individual trials and mean
%         n_plot_trials = 10;
%         curr_bin_idx = find(corr_bin_idx{curr_days} == curr_bin);
%         plot_trials = curr_bin_idx(randperm(length(curr_bin_idx),n_plot_trials));
%         plot(curr_move_basesub(:,plot_trials),'color',[0.5,0.5,0.5]);
%         plot(nanmean(curr_move_basesub(:,curr_bin_idx),2),'color','k','linewidth',2);
%         
%         if curr_days == 1
%             title(['Percentile bin ',num2str(curr_bin)]);
%         end
%         if curr_bin == 1;
%             ylabel(['Days ' num2str(use_days{curr_days})]);
%         end
%     end
% end



%% Reviewer 1 point 6: fig 4E timing changes

clearvars -except data analysis classified_rois

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

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed > 0;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed > 0;

        
    end
    disp(curr_animal);
end

move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_all,'uni',false);
move_act_m_onset_mean = cell(size(move_act_onset_mean));
move_act_q_onset_mean = cell(size(move_act_onset_mean));

move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_all,'uni',false);
move_act_m_offset_mean = cell(size(move_act_offset_mean));
move_act_q_offset_mean = cell(size(move_act_offset_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        move_act_m_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_onset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);
        
        move_act_m_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).movement(:,curr_day) == 1),1),3);
             
        move_act_q_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(move_act_offset_all{curr_animal,curr_day}(:,:, ...
            classified_rois(curr_animal).quiescent(:,curr_day) == 1),1),3);       
        
    end
end

% Get mean activity across days for each animal
move_act_mean_onset_animalavg = nan(length(data),total_frames);
move_act_mean_m_onset_animalavg = nan(length(data),total_frames);
move_act_mean_q_onset_animalavg = nan(length(data),total_frames);

move_act_mean_offset_animalavg = nan(length(data),total_frames);
move_act_mean_m_offset_animalavg = nan(length(data),total_frames);
move_act_mean_q_offset_animalavg = nan(length(data),total_frames);

for curr_animal = 1:length(data)
    
    move_act_mean_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_onset_mean{curr_animal,:}),1);
    move_act_mean_m_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_onset_mean{curr_animal,:}),1);
    move_act_mean_q_onset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_onset_mean{curr_animal,:}),1);
    
    move_act_mean_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_offset_mean{curr_animal,:}),1);
    move_act_mean_m_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_m_offset_mean{curr_animal,:}),1);
    move_act_mean_q_offset_animalavg(curr_animal,:) = nanmean(vertcat(move_act_q_offset_mean{curr_animal,:}),1);
    
end

move_act_mean_q_onset_dayavg = nan(14,total_frames);
for curr_day = 1:14
    move_act_mean_q_onset_dayavg(curr_day,:) = nanmean(vertcat(move_act_q_onset_mean{:,curr_day}),1);
end
figure;hold on
set(gca,'ColorOrder',copper(14));
plot(move_act_mean_q_onset_dayavg','linewidth',2)

move_act_mean_onset_dayavg = nan(14,total_frames);
for curr_day = 1:14
    move_act_mean_onset_dayavg(curr_day,:) = nanmean(vertcat(move_act_onset_mean{:,curr_day}),1);
end
move_act_mean_offset_dayavg = nan(14,total_frames);
for curr_day = 1:14
    move_act_mean_offset_dayavg(curr_day,:) = nanmean(vertcat(move_act_offset_mean{:,curr_day}),1);
end
figure;
subplot(2,2,1);hold on
set(gca,'ColorOrder',copper(14));
plot(move_act_mean_onset_dayavg','linewidth',2)
subplot(2,2,2);hold on
set(gca,'ColorOrder',copper(14));
plot(move_act_mean_offset_dayavg','linewidth',2)
subplot(2,2,3);
plot(nanmean(move_act_mean_onset_dayavg,2),'k','linewidth',2)
ylabel('Average activity');
xlabel('Day');
subplot(2,2,4);
plot(nanmean(move_act_mean_offset_dayavg,2),'k','linewidth',2)
ylabel('Average activity');
xlabel('Day');


move_act_mean_q_offset_dayavg = nan(14,total_frames);
for curr_day = 1:14
    move_act_mean_q_offset_dayavg(curr_day,:) = nanmean(vertcat(move_act_q_offset_mean{:,curr_day}),1);
end





% Activity of all animals and cells, day averaged (like 4A, grand average)

move_act_onset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_onset_all,'uni',false);
move_act_offset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_offset_all,'uni',false);

curr_data_onset_norm_mean = cell(size(data));
curr_data_offset_norm_mean = cell(size(data));
for curr_animal = 1:length(data)
    
    use_days = cellfun(@(x) ~isempty(x),move_act_onset_all_mean(curr_animal,:));
    
    curr_data_onset = cat(3,move_act_onset_all_mean{curr_animal,use_days});
    curr_data_offset = cat(3,move_act_offset_all_mean{curr_animal,use_days});
    
    curr_data_onset_minsub = bsxfun(@minus,curr_data_onset,min(curr_data_onset,[],2));
    curr_data_onset_norm = bsxfun(@times,curr_data_onset_minsub,1./max(curr_data_onset_minsub,[],2));
    
    curr_data_offset_minsub = bsxfun(@minus,curr_data_offset,min(curr_data_offset,[],2));
    curr_data_offset_norm = bsxfun(@times,curr_data_offset_minsub,1./max(curr_data_offset_minsub,[],2));
    
    %curr_data_onset_norm_mean{curr_animal} = nanmean(curr_data_onset_norm,3);
    %curr_data_offset_norm_mean{curr_animal} = nanmean(curr_data_offset_norm,3);
    
    curr_data_onset_norm_mean{curr_animal} = nanmean(curr_data_onset,3);
    curr_data_offset_norm_mean{curr_animal} = nanmean(curr_data_offset,3);
    
end

all_onset = cell2mat(curr_data_onset_norm_mean);
all_offset = cell2mat(curr_data_offset_norm_mean);

all_combined = [all_onset,all_offset];
m_frames = false(size(all_combined,2),1);
m_frames(nonmove_frames+1:nonmove_frames+1+move_frames*2) = true;
m_act = nanmean(all_combined(:,m_frames),2);
q_act = nanmean(all_combined(:,~m_frames),2);
mq_ratio = (m_act - q_act)./(q_act + m_act);

[~,sort_idx] = sort(mq_ratio);

figure; colormap(hot);
subplot(1,2,1); imagesc(all_onset(sort_idx,:));
subplot(1,2,2); imagesc(all_offset(sort_idx,:));



% Activity of all animals and cells, day separate 
use_days = {[1,2],[3,4]};
act_days_norm = cell(length(use_days),1);
for curr_use_days = 1:length(use_days)
    
    curr_days = use_days{curr_use_days};
    
    curr_data_onset = nanmean(cell2mat(permute(move_act_onset_all_mean(:,curr_days),[1,3,2])),3);
    curr_data_onset_minsub = bsxfun(@minus,curr_data_onset,min(curr_data_onset,[],2));
    curr_data_onset_norm = bsxfun(@times,curr_data_onset_minsub,1./max(curr_data_onset_minsub,[],2));    
    
    curr_data_offset = nanmean(cell2mat(permute(move_act_offset_all_mean(:,curr_days),[1,3,2])),3);
    curr_data_offset_minsub = bsxfun(@minus,curr_data_offset,min(curr_data_offset,[],2));
    curr_data_offset_norm = bsxfun(@times,curr_data_offset_minsub,1./max(curr_data_offset_minsub,[],2));
    
    curr_cat = [curr_data_onset_norm,curr_data_offset_norm];  
    
    act_days_norm{curr_use_days} = curr_cat(sort_idx,:);   
    
end

figure; colormap(hot);
for i = 1:length(use_days);
    subplot(1,length(use_days),i);
    imagesc(act_days_norm{i});
end


% Like 4C but over days
use_days = {[1,2],[3,4],[5,6]};
mq_ratio = cell(length(use_days),1);
for curr_use_days = 1:length(use_days)
    
    curr_days = use_days{curr_use_days};
    
    curr_data_onset = nanmean(cell2mat(permute(move_act_onset_all_mean(:,curr_days),[1,3,2])),3);
    curr_data_offset = nanmean(cell2mat(permute(move_act_offset_all_mean(:,curr_days),[1,3,2])),3);
    
    curr_cat = [curr_data_onset,curr_data_offset];
    
    m_frames = false(size(all_combined,2),1);
    m_frames(nonmove_frames+1:nonmove_frames+1+move_frames*2) = true;
    
    m_act = nanmean(curr_cat(:,m_frames),2);
    q_act = nanmean(curr_cat(:,~m_frames),2);
    mq_ratio{curr_use_days} = (m_act - q_act)./(q_act + m_act);
    
end

figure;
for curr_use_days = 1:length(use_days)
    subplot(length(use_days),1,curr_use_days);
    bin_edges = linspace(-1,1,20);
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    hist(mq_ratio{curr_use_days},bin_centers);
end


%% Reviewer 2 point 1: 5A jump in quiescent

clearvars -except data analysis classified_rois

% Get activity of all cells during m/q/total 
avg_movement_act = cell(length(data),1);
avg_quiescent_act = cell(length(data),1);
avg_total_act = cell(length(data),1);

for curr_animal = 1:length(data)
        
    use_cells = true(size(classified_rois(curr_animal).movement,1),1);
        
    curr_avg_movement_act = nan(length(use_cells),14);
    curr_avg_quiescent_act = nan(length(use_cells),14);
    curr_avg_total_act = nan(length(use_cells),14);
    
     for i = 1:14
         
         if length(data(curr_animal).im) < i || ...
                 isempty(data(curr_animal).im(i).roi_trace_thresh)
             continue
         end
         
        curr_avg_movement_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames > 0),2);
        
        curr_avg_quiescent_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,analysis(curr_animal).lever(i).lever_move_frames == 0),2);
        
        curr_avg_total_act(use_cells,i) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(use_cells,:),2);
        
     end
    
    avg_movement_act{curr_animal} = curr_avg_movement_act;
    avg_quiescent_act{curr_animal} = curr_avg_quiescent_act;
    avg_total_act{curr_animal} = curr_avg_total_act;
    
    disp(curr_animal);
    
end

% Get categories
m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);

% Get cumulative activity during quiescence
q_act_cat = vertcat(avg_quiescent_act{:});
q_act_cat_sort = sort(q_act_cat,1,'descend');
q_act_cat_norm = bsxfun(@rdivide,cumsum(q_act_cat_sort,1),sum(q_act_cat,1));

m_act_cat = vertcat(avg_movement_act{:});
m_act_cat_sort = sort(m_act_cat,1,'descend');
m_act_cat_norm = bsxfun(@rdivide,cumsum(m_act_cat_sort,1),sum(m_act_cat,1));

n_shuff = 10000;

q_act_cat_norm_diff_shuff_1234 = nan(size(q_act_cat_norm,1),n_shuff);
q_act_cat_norm_diff_shuff_3456 = nan(size(q_act_cat_norm,1),n_shuff);

m_act_cat_norm_diff_shuff_1234 = nan(size(m_act_cat_norm,1),n_shuff);
m_act_cat_norm_diff_shuff_3456 = nan(size(m_act_cat_norm,1),n_shuff);

for curr_shuff = 1:n_shuff
    
    curr_shuff_data = shake(q_act_cat(:,1:4),2);
    curr_shuff_sort = sort(curr_shuff_data,1,'descend');
    curr_shuff_norm = bsxfun(@rdivide,cumsum(curr_shuff_sort,1),sum(curr_shuff_data,1));
    curr_shuff_diff = nanmean(curr_shuff_norm(:,1:2),2) - ...
        nanmean(curr_shuff_norm(:,3:4),2);
    q_act_cat_norm_diff_shuff_1234(:,curr_shuff) = curr_shuff_diff;   
    
    curr_shuff_data = shake(q_act_cat(:,3:6),2);
    curr_shuff_sort = sort(curr_shuff_data,1,'descend');
    curr_shuff_norm = bsxfun(@rdivide,cumsum(curr_shuff_sort,1),sum(curr_shuff_data,1));
    curr_shuff_diff = nanmean(curr_shuff_norm(:,1:2),2) - ...
        nanmean(curr_shuff_norm(:,3:4),2);
    q_act_cat_norm_diff_shuff_3456(:,curr_shuff) = curr_shuff_diff;
    
    curr_shuff_data = shake(m_act_cat(:,1:4),2);
    curr_shuff_sort = sort(curr_shuff_data,1,'descend');
    curr_shuff_norm = bsxfun(@rdivide,cumsum(curr_shuff_sort,1),sum(curr_shuff_data,1));
    curr_shuff_diff = nanmean(curr_shuff_norm(:,1:2),2) - ...
        nanmean(curr_shuff_norm(:,3:4),2);
    m_act_cat_norm_diff_shuff_1234(:,curr_shuff) = curr_shuff_diff;   
    
    curr_shuff_data = shake(m_act_cat(:,3:6),2);
    curr_shuff_sort = sort(curr_shuff_data,1,'descend');
    curr_shuff_norm = bsxfun(@rdivide,cumsum(curr_shuff_sort,1),sum(curr_shuff_data,1));
    curr_shuff_diff = nanmean(curr_shuff_norm(:,1:2),2) - ...
        nanmean(curr_shuff_norm(:,3:4),2);
    m_act_cat_norm_diff_shuff_3456(:,curr_shuff) = curr_shuff_diff;
    
    disp(curr_shuff);
    
end

q_act_cat_norm_diff_1234_95 = prctile(q_act_cat_norm_diff_shuff_1234,[2.5,97.5],2);
q_act_cat_norm_diff_3456_95 = prctile(q_act_cat_norm_diff_shuff_3456,[2.5,97.5],2);

m_act_cat_norm_diff_1234_95 = prctile(m_act_cat_norm_diff_shuff_1234,[2.5,97.5],2);
m_act_cat_norm_diff_3456_95 = prctile(m_act_cat_norm_diff_shuff_3456,[2.5,97.5],2);

figure;

subplot(2,2,1); hold on;
plot(nanmean(q_act_cat_norm(:,1:2),2),'k','linewidth',2)
plot(nanmean(q_act_cat_norm(:,3:4),2),'r','linewidth',2)
plot(nanmean(q_act_cat_norm(:,5:6),2),'b','linewidth',2)
ylabel('Cumulative proportion of activity');
xlabel('Number of ROIs');
legend({'Days 1:2','Days 3:4','Days 5:6'})
title('Activity during quiescence')

subplot(2,2,2); hold on;
plot(nanmean(q_act_cat_norm(:,1:2),2) - nanmean(q_act_cat_norm(:,3:4),2),'k','linewidth',2);
plot(q_act_cat_norm_diff_1234_95,'--k','linewidth',2);
plot(nanmean(q_act_cat_norm(:,3:4),2) - nanmean(q_act_cat_norm(:,5:6),2),'r','linewidth',2);
plot(q_act_cat_norm_diff_3456_95,'--r','linewidth',2);
 ylabel('Area under curve difference')
xlabel('Number of ROIs')
legend({'Days 1:2 - Days 3:4','CI','CI','Days 1:2 - Days 3:4','CI','CI'})
title('Activity during quiescence')

subplot(2,2,3); hold on;
plot(nanmean(m_act_cat_norm(:,1:2),2),'k','linewidth',2)
plot(nanmean(m_act_cat_norm(:,3:4),2),'r','linewidth',2)
plot(nanmean(m_act_cat_norm(:,5:6),2),'b','linewidth',2)
ylabel('Cumulative proportion of activity');
xlabel('Number of ROIs');
legend({'Days 1:2','Days 3:4','Days 5:6'})
title('Activity during movement')

subplot(2,2,4); hold on;
plot(nanmean(m_act_cat_norm(:,1:2),2) - nanmean(m_act_cat_norm(:,3:4),2),'k','linewidth',2);
plot(m_act_cat_norm_diff_1234_95,'--k','linewidth',2);
plot(nanmean(m_act_cat_norm(:,3:4),2) - nanmean(m_act_cat_norm(:,5:6),2),'r','linewidth',2);
plot(m_act_cat_norm_diff_3456_95,'--r','linewidth',2);
 ylabel('Area under curve difference')
xlabel('Number of ROIs')
legend({'Days 1:2 - Days 3:4','CI','CI','Days 1:2 - Days 3:4','CI','CI'})
title('Activity during movement')

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;

u2q = cellfun(@(m,q) ...
    sum(q(:,1:2),2) < min_classdays & ...
    sum(q(:,3:7),2) >= min_classdays,m_pad,q_pad,'uni',false);

u2q_m_act_avg = cellfun(@(act,class) nanmean(act(class,:),1), ...
    avg_movement_act,u2q','uni',false);
u2q_q_act_avg = cellfun(@(act,class) nanmean(act(class,:),1), ...
    avg_quiescent_act,u2q','uni',false);
u2q_total_act_avg = cellfun(@(act,class) nanmean(act(class,:),1), ...
    avg_total_act,u2q','uni',false);

u2q_total_act_avg_cat = vertcat(u2q_total_act_avg{:});
u2q_m_act_avg_cat = vertcat(u2q_m_act_avg{:});
u2q_q_act_avg_cat = vertcat(u2q_q_act_avg{:});

figure; hold on;
errorbar(nanmean(u2q_total_act_avg_cat),nanstd(u2q_total_act_avg_cat)./sqrt(sum(~isnan(u2q_total_act_avg_cat))),'k')
errorbar(nanmean(u2q_m_act_avg_cat),nanstd(u2q_m_act_avg_cat)./sqrt(sum(~isnan(u2q_m_act_avg_cat))),'r')
errorbar(nanmean(u2q_q_act_avg_cat),nanstd(u2q_q_act_avg_cat)./sqrt(sum(~isnan(u2q_q_act_avg_cat))),'b')
xlabel('day')
ylabel('average df/f')
title('~q2q')
legend({'total','during movement','during quiescence'})

%% Reviewer 2 point 2.1: example fluorescence traces

% Another example
% use_animal = 8;
% use_day = 14;
% use_frames = [37271:40297];
% example_cells_q = [47,93];
% example_cells_m = [45,147];

use_animal = 2;
use_day = 12;
use_frames = [35837:40877];
example_cells_q = [5,29,124];
example_cells_m = [49,9,161];
example_cells_ua = [23,31,67,102];

m = classified_rois(use_animal).movement(:,use_day);
q = classified_rois(use_animal).quiescent(:,use_day);
ua = classified_rois(use_animal).unclassified_active(:,use_day);

curr_data = data(use_animal).im(use_day).roi_trace_df;
curr_data_thresh = data(use_animal).im(use_day).roi_trace_thresh;

curr_move_trace = analysis(use_animal).lever(use_day).lever_move_frames';
curr_lever_trace = data(use_animal).bhv(use_day).lever_force_plot(1:size(curr_data,2))';

figure;

p1 = subplot(3,1,1);hold on;
AP_stackplot([nanmean(curr_data(:,use_frames),1); ...
    nanmean(curr_data(m,use_frames),1); ...
    nanmean(curr_data(q,use_frames),1); ...
    nanmean(curr_data(ua,use_frames),1)]',[],1,false, ...
    [0,0,0;0,1,0;1,0,0;0,0,1])
AP_stackplot([nanmean(curr_data_thresh(:,use_frames),1); ...
    nanmean(curr_data_thresh(m,use_frames),1); ...
    nanmean(curr_data_thresh(q,use_frames),1); ...
    nanmean(curr_data_thresh(ua,use_frames),1)]',[],1,false,'k')

p2 = subplot(3,1,2);hold on;
AP_stackplot(curr_data(...
    [example_cells_m,example_cells_q,example_cells_ua],use_frames)',...
    [],5,false,[repmat([0,1,0],length(example_cells_m),1); ...
    repmat([1,0,0],length(example_cells_q),1); ...
    repmat([0,0,1],length(example_cells_ua),1)]);
AP_stackplot(curr_data_thresh(...
    [example_cells_m,example_cells_q,example_cells_ua],use_frames)',...
    [],5,false,'k');

p3 = subplot(3,1,3);hold on;
plot(curr_lever_trace(use_frames),'b');
plot(curr_move_trace(use_frames),'k');

linkaxes([p1,p2,p3],'x');


    
%% Reviewer 2 point 2.3 (and 1): 4B longer + changes over days


clearvars -except data analysis classified_rois

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
        nonmove_frames = 28*3;
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

% Get activity CONTRIBUTION from cell classes (zero, don't ignore);
move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_all,'uni',false);
move_act_m_onset_mean = cell(size(move_act_onset_mean));
move_act_q_onset_mean = cell(size(move_act_onset_mean));
move_act_ua_onset_mean = cell(size(move_act_onset_mean));

move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_all,'uni',false);
move_act_m_offset_mean = cell(size(move_act_offset_mean));
move_act_q_offset_mean = cell(size(move_act_offset_mean));
move_act_ua_offset_mean = cell(size(move_act_offset_mean));

for curr_animal = 1:length(data)
    for curr_day = 1:min(length(data(curr_animal).im),14)
        
        curr_onset = move_act_onset_all{curr_animal,curr_day};
        curr_offset = move_act_offset_all{curr_animal,curr_day};
        
        % For control animal
        if isempty(curr_onset)
            continue
        end

        curr_onset_class = bsxfun(@times,curr_onset, ...
            permute(classified_rois(curr_animal).movement(:,curr_day),[2,3,1]));
        move_act_m_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_onset_class,1),3);
        
        curr_onset_class = bsxfun(@times,curr_onset, ...
            permute(classified_rois(curr_animal).quiescent(:,curr_day),[2,3,1]));
        move_act_q_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_onset_class,1),3);
        
        curr_onset_class = bsxfun(@times,curr_onset, ...
            permute(classified_rois(curr_animal).unclassified_active(:,curr_day),[2,3,1]));
        move_act_ua_onset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_onset_class,1),3);
        
        curr_offset_class = bsxfun(@times,curr_offset, ...
            permute(classified_rois(curr_animal).movement(:,curr_day),[2,3,1]));
        move_act_m_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_offset_class,1),3);
        
        curr_offset_class = bsxfun(@times,curr_offset, ...
            permute(classified_rois(curr_animal).quiescent(:,curr_day),[2,3,1]));
        move_act_q_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_offset_class,1),3);
        
        curr_offset_class = bsxfun(@times,curr_offset, ...
            permute(classified_rois(curr_animal).unclassified_active(:,curr_day),[2,3,1]));
        move_act_ua_offset_mean{curr_animal,curr_day} = ...
            nanmean(nanmean(curr_offset_class,1),3);       
        
    end
end

%group_days = {[1,2],[3,4],[5,6],[7,8],[9,10],[11,12],[13,14]};
group_days = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};
% Mean activity across days
move_act_onset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_onset_mean{:,x}),1),group_days,'uni',false)');
move_act_offset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_offset_mean{:,x}),1),group_days,'uni',false)');
% Mean activity contribution from classified cells
move_act_m_onset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_m_onset_mean{:,x}),1),group_days,'uni',false)');
move_act_m_offset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_m_offset_mean{:,x}),1),group_days,'uni',false)');
move_act_q_onset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_q_onset_mean{:,x}),1),group_days,'uni',false)');
move_act_q_offset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_q_offset_mean{:,x}),1),group_days,'uni',false)');
move_act_ua_onset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_ua_onset_mean{:,x}),1),group_days,'uni',false)');
move_act_ua_offset_mean_days = cell2mat(cellfun(@(x) ...
    nanmean(vertcat(move_act_ua_offset_mean{:,x}),1),group_days,'uni',false)');

figure; 
subplot(4,2,1); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_onset_mean_days','linewidth',2);
xlim([0 total_frames])
title('All onset')
subplot(4,2,2); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_offset_mean_days','linewidth',2);
xlim([0 total_frames])
title('All offset')
subplot(4,2,3); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_m_onset_mean_days','linewidth',2);
xlim([0 total_frames])
title('M onset')
subplot(4,2,4); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_m_offset_mean_days','linewidth',2);
xlim([0 total_frames])
title('M offset')
subplot(4,2,5); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_q_onset_mean_days','linewidth',2);
xlim([0 total_frames])
title('Q onset')
subplot(4,2,6); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_q_offset_mean_days','linewidth',2);
xlim([0 total_frames])
title('Q offset')
subplot(4,2,7); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_ua_onset_mean_days','linewidth',2);
xlim([0 total_frames])
title('UA onset')
subplot(4,2,8); hold on; set(gca,'ColorOrder',copper(length(group_days)));
plot(move_act_ua_offset_mean_days','linewidth',2);
xlim([0 total_frames])
title('UA offset')


%% Reviewer 2 point 2.4: 5B traces

clearvars -except data analysis classified_rois

% Get categories

m = {classified_rois(:).movement};
q = {classified_rois(:).quiescent};
ua = {classified_rois(:).unclassified_active};

m_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),m,'uni',false);
q_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),q,'uni',false);
ua_pad = cellfun(@(x) padarray(+x(:,1:min(size(x,2),14)),[0,max(14-size(x,2),0)],NaN,'post'),ua,'uni',false);

% Week 1 vs. week 2 classification switch
% NOTE: THERES PROBABLY OVERLAP HERE
min_classdays = 1;
m2q = cellfun(@(m,q) ...
    sum(m(:,1:7),2) >= min_classdays & ...
    sum(q(:,8:end),2) >= min_classdays,m_pad,q_pad,'uni',false);

% Get average onset/offset activity of m2q cells
move_act_onset_m2q = cell(length(data),14);
move_act_offset_m2q = cell(length(data),14);

for curr_animal = 1:length(data);
    for curr_day = 1:min(length(data(curr_animal).im),14);
            
        % Get  activity
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(m2q{curr_animal},:));
        
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

  
        move_act_onset_m2q{curr_animal,curr_day} = move_act_onset_timed > 0;
        move_act_offset_m2q{curr_animal,curr_day} = move_act_offset_timed > 0;

        
    end
    disp(curr_animal);
end

move_act_onset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_onset_m2q,'uni',false);
move_act_offset_mean = cellfun(@(x) nanmean(nanmean(x,1),3),move_act_offset_m2q,'uni',false);

use_days = num2cell(1:14);
move_act_onset_mean_dayavg = nan(length(use_days),total_frames);
move_act_offset_mean_dayavg = nan(length(use_days),total_frames);
for curr_use_days = 1:length(use_days)
    
    curr_days = use_days{curr_use_days};
    
    move_act_onset_mean_dayavg(curr_use_days,:) = ...
        nanmean(vertcat(move_act_onset_mean{:,curr_days}),1);
    move_act_offset_mean_dayavg(curr_use_days,:) = ...
        nanmean(vertcat(move_act_offset_mean{:,curr_days}),1);  
end

figure;

subplot(1,2,1); hold on;
t = [-nonmove_frames:move_frames]/28;
set(gca,'ColorOrder',copper(length(use_days)));
plot(t,move_act_onset_mean_dayavg','linewidth',2)
line([0,0],ylim,'color','b');
title('Movement onset');
xlabel('Time (s)');
ylabel('M2Q average \DeltaF/F');


subplot(1,2,2); hold on;
t = [-move_frames:nonmove_frames]/28;
set(gca,'ColorOrder',copper(length(use_days)));
plot(t,move_act_offset_mean_dayavg','linewidth',2)
line([0,0],ylim,'color','b');
title('Movement offset');
xlabel('Time (s)');
ylabel('M2Q average \DeltaF/F');

%% Reviewer 2 point 2.5: fig 6 examples

%%%% WORK IN PROGRESS

clearvars -except data analysis classified_rois

session_epochs = {[1:4],[11:14]};

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
            use_time = -100:2000; % ms
            use_frames = round(use_time(1)/1000*28):round(use_time(end)/1000*28);
            
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
            edge_movements = ...
                use_movement_start_frames+use_time(1) < 0 | ...
                use_movement_start_times + use_time(end) > ...
                length(lever_force_resample) | ...
                use_movement_start_frames + use_frames(1) < 0 | ...
                (use_movement_start_frames + use_frames(end)) > ...
                size(data(curr_animal).im(curr_session).roi_trace_df,2);
                     
            use_movement_start_times(edge_movements) = [];
            use_movement_start_frames(edge_movements) = [];
            
            movements_idx = bsxfun(@plus,use_movement_start_times,use_time);
            timed_movements = lever_force_resample(movements_idx');
            
            timed_activity = permute(reshape(data(curr_animal).im(curr_session).roi_trace_thresh(:, ...
                bsxfun(@plus,use_movement_start_frames,use_frames)')',length(use_frames),[],...
                size(data(curr_animal).im(curr_session).roi_trace_thresh,1)),[1,3,2]);
            
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

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

% Make template get trial-by-trial activity/movement correlation
corr_bin_use = [-Inf -0.5 0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin_valid = corr_bin_use;
corr_bin_valid(1) = -1;
corr_bin_valid(end) = 1;
corr_bin_plot = [corr_bin_valid(1:end-1) + diff(corr_bin_valid)/2 NaN];

all_act_corr = cell(size(epoch_movement_cat,1),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));

% Make template movement from median of end sessions movements
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,size(epoch_movement_cat,2)),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

n_trials = cellfun(@(x) cellfun(@(x) size(x,3),x),epoch_activity,'uni',false);
movement_template_corr_split = cellfun(@(cor,n_trials) ...
    mat2cell(cor,n_trials),movement_template_corr,n_trials,'uni',false);


n_trials_across = mat2cell(cellfun(@sum,n_trials),ones(length(data),1),2);

epoch_movement_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_movement_cat{x,:})), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_movement_corr_split = cellfun(@(cor,n_trials) ...
    mat2cell(cor,n_trials,n_trials),epoch_movement_corr_all',n_trials_across,'uni',false);

epoch_activity_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_activity_cat_reshape{x,:}),'rows','pairwise'), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
epoch_activity_corr_split = cellfun(@(cor,n_trials) ...
    mat2cell(cor,n_trials,n_trials),epoch_activity_corr_all',n_trials_across,'uni',false);




% plot cell+temporal activity PCs for trials
curr_animal = 2;

act_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(curr_animal,:),'uni',false);
move_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(curr_animal,:),'uni',false);

smooth_kernel = ones(10,1,1)./10;
act_smooth = cellfun(@(x) convn(x,smooth_kernel,'same'),act_cat,'uni',false);

[n_frames,n_cells,n_trials] = size(cat(3,act_smooth{:}));
[coeff,score,latent] = cellfun(@(x) ...
    princomp(reshape(x,[],size(x,3))'),act_smooth,'uni',false);
score_reshape = cellfun(@(x) reshape(x,n_frames,[],n_cells),score,'uni',false);




% plot PCs across all trials like/unlike
curr_animal = 5;
[~,sort_idx] = cellfun(@(x) sort(vertcat(x{:})),movement_template_corr_split(curr_animal,:),'uni',false);

act_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(curr_animal,:),'uni',false);
move_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(curr_animal,:),'uni',false);

smooth_kernel = ones(10,1,1)./10;
act_smooth = cellfun(@(x) convn(x,smooth_kernel,'same'),act_cat,'uni',false);

[n_frames,n_cells,n_trials] = size(cat(3,act_smooth{:}));
[coeff,score,latent] = princomp(reshape(permute(cat(3,act_smooth{:}),[1,3,2]),[],n_cells));
score_reshape = mat2cell(reshape(score,n_frames,n_trials,n_cells), ...
    n_frames,cellfun(@(x) size(x,3),act_smooth),n_cells);

avg_trials = 50;
figure; 

subplot(1,2,1),hold on;
plot3(nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),2),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),3),2),'k')
plot3(nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),2),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),3),2),'r');
view(3);
axis vis3d;

subplot(1,2,2),hold on;
plot3(nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),2),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),3),2),'k')
plot3(nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),2),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),3),2),'r');
view(3);
axis vis3d;



avg_trials = 50;
figure;
subplot(1,2,1); hold on;
plot(nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),2),2),'k')
plot(nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),2),2),'r');

subplot(1,2,2); hold on;
plot(nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),2),2),'k')
plot(nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),2),2),'r');



%%%%% alt: PCA within each epoch of trials
curr_animal = 3;

[~,sort_idx] = cellfun(@(x) sort(vertcat(x{:})),movement_template_corr_split(curr_animal,:),'uni',false);

act_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(curr_animal,:),'uni',false);
move_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(curr_animal,:),'uni',false);

% remove nan cells (L2/3 has some?)
nan_cells = any(any(isnan(cat(3,act_cat{:})),1),3);
act_cat = cellfun(@(x) x(:,~nan_cells,:),act_cat,'uni',false);

smooth_kernel = ones(10,1,1)./10;
act_smooth = cellfun(@(x) convn(x,smooth_kernel,'same'),act_cat,'uni',false);

[n_frames,n_cells,n_trials] = size(cat(3,act_smooth{:}));
[coeff,score,latent] = cellfun(@(x) ...
    princomp(reshape(permute(x,[1,3,2]),[],n_cells)),act_smooth,'uni',false);
score_reshape = cellfun(@(x) reshape(x,n_frames,[],n_cells),score,'uni',false);

avg_trials = 50;
figure; 

subplot(1,2,1),hold on;
plot3(nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),2),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),3),2),'k')
plot3(nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),2),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),3),2),'r');
view(3);
axis vis3d;

subplot(1,2,2),hold on;
plot3(nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),2),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),3),2),'k')
plot3(nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),2),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),3),2),'r');
view(3);
axis vis3d;


avg_trials = 50;
figure;
subplot(1,2,1); hold on;
plot(nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(1:avg_trials),2),2),'k')
plot(nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{1}(:,sort_idx{1}(end-avg_trials+1:end),2),2),'r');

subplot(1,2,2); hold on;
plot(nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(1:avg_trials),2),2),'k')
plot(nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),1),2), ...
    nanmean(score_reshape{2}(:,sort_idx{2}(end-avg_trials+1:end),2),2),'r');


%%% ALT: just plot average activity
curr_animal = 6;
avg_trials = 50;

[~,sort_idx] = cellfun(@(x) sort(vertcat(x{:})),movement_template_corr_split(curr_animal,:),'uni',false);
act_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(curr_animal,:),'uni',false);
move_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(curr_animal,:),'uni',false);
smooth_kernel = ones(10,1,1)./10;
act_smooth = cellfun(@(x) convn(x,smooth_kernel,'same'),act_cat,'uni',false);

early_avg_act_like = nanmean(act_smooth{1}(:,:,sort_idx{1}(end-avg_trials+1:end)),3)';
early_tr_act_like = permute(nanmean(act_smooth{1}(:,:,sort_idx{1}(end-avg_trials+1:end)),1),[2,3,1]);
early_move_like = nanmean(move_cat{1}(:,sort_idx{1}(end-avg_trials+1:end)),2);

early_avg_act_unlike = nanmean(act_smooth{1}(:,:,sort_idx{1}(1:avg_trials)),3)';
early_tr_act_unlike = permute(nanmean(act_smooth{1}(:,:,sort_idx{1}(1:avg_trials)),1),[2,3,1]);
early_move_unlike = nanmean(move_cat{1}(:,sort_idx{1}(1:avg_trials)),2);

late_avg_act_like = nanmean(act_smooth{2}(:,:,sort_idx{2}(end-avg_trials+1:end)),3)';
late_tr_act_like = permute(nanmean(act_smooth{2}(:,:,sort_idx{2}(end-avg_trials+1:end)),1),[2,3,1]);
late_move_like = nanmean(move_cat{2}(:,sort_idx{2}(end-avg_trials+1:end)),2);

late_avg_act_unlike = nanmean(act_smooth{2}(:,:,sort_idx{2}(1:avg_trials)),3)';
late_tr_act_unlike = permute(nanmean(act_smooth{2}(:,:,sort_idx{2}(1:avg_trials)),1),[2,3,1]);
late_move_unlike = nanmean(move_cat{2}(:,sort_idx{2}(1:avg_trials)),2);

figure; colormap(hot);

subplot(3,4,1); 
imagesc(early_avg_act_like);
subplot(3,4,2); 
imagesc(early_avg_act_unlike);
subplot(3,4,3); 
imagesc(late_avg_act_like);
subplot(3,4,4); 
imagesc(late_avg_act_unlike);

subplot(3,4,5); 
imagesc(early_tr_act_like);
subplot(3,4,6); 
imagesc(early_tr_act_unlike);
subplot(3,4,7); 
imagesc(late_tr_act_like);
subplot(3,4,8); 
imagesc(late_tr_act_unlike);

subplot(3,4,9); 
plot(early_move_like,'k');
subplot(3,4,10); 
plot(early_move_unlike,'k');
subplot(3,4,11); 
plot(late_move_like,'k');
subplot(3,4,12); 
plot(late_move_unlike,'k');

%%% ALT: as above, but all cells concatenated

temporal_act = cell(length(data),4);
trial_act = cell(length(data),4);

avg_trials = 50;

for curr_animal = 1:length(data);
    
    [~,sort_idx] = cellfun(@(x) sort(vertcat(x{:})),movement_template_corr_split(curr_animal,:),'uni',false);
    act_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(curr_animal,:),'uni',false);
    move_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(curr_animal,:),'uni',false);
    smooth_kernel = ones(10,1,1)./10;
    act_smooth = cellfun(@(x) convn(x,smooth_kernel,'same'),act_cat,'uni',false);
    
    temporal_act{curr_animal,1} = nanmean(act_smooth{1}(:,:,sort_idx{1}(end-avg_trials+1:end)),3)';
    trial_act{curr_animal,1} = permute(nanmean(act_smooth{1}(:,:,sort_idx{1}(end-avg_trials+1:end)),1),[2,3,1]);
    
    temporal_act{curr_animal,2} = nanmean(act_smooth{1}(:,:,sort_idx{1}(1:avg_trials)),3)';
    trial_act{curr_animal,2} = permute(nanmean(act_smooth{1}(:,:,sort_idx{1}(1:avg_trials)),1),[2,3,1]);
    
    temporal_act{curr_animal,3} = nanmean(act_smooth{2}(:,:,sort_idx{2}(end-avg_trials+1:end)),3)';
    trial_act{curr_animal,3} = permute(nanmean(act_smooth{2}(:,:,sort_idx{2}(end-avg_trials+1:end)),1),[2,3,1]);
    
    temporal_act{curr_animal,4} = nanmean(act_smooth{2}(:,:,sort_idx{2}(1:avg_trials)),3)';
    trial_act{curr_animal,4} = permute(nanmean(act_smooth{2}(:,:,sort_idx{2}(1:avg_trials)),1),[2,3,1]);
    
end

cat_temporal_act = cell(1,4);
cat_trial_act = cell(1,4);
for i = 1:4    
    cat_temporal_act{i} = vertcat(temporal_act{:,i});
    cat_trial_act{i} = vertcat(trial_act{:,i});
end

[~,max_idx] = max(nanmean(cat(3,cat_temporal_act{:}),3),[],2);
[~,sort_idx] = sort(max_idx);

figure; colormap(hot);

subplot(3,4,1); 
imagesc(cat_temporal_act{1}(sort_idx,:));
subplot(3,4,2); 
imagesc(cat_temporal_act{2}(sort_idx,:));
subplot(3,4,3); 
imagesc(cat_temporal_act{3}(sort_idx,:));
subplot(3,4,4); 
imagesc(cat_temporal_act{4}(sort_idx,:));

subplot(3,4,5); 
imagesc(cat_trial_act{1}(sort_idx,:));
subplot(3,4,6); 
imagesc(cat_trial_act{2}(sort_idx,:));
subplot(3,4,7); 
imagesc(cat_trial_act{3}(sort_idx,:));
subplot(3,4,8); 
imagesc(cat_trial_act{4}(sort_idx,:));




%%% ALT ATTEMPT: sort each movement separately by pairs
curr_animal = 4;
figure; hold on;

curr_move = epoch_movement_corr_split{curr_animal}{1,1};
[curr_move_sort,sort_idx] = sort(curr_move,2);

curr_act = epoch_activity_corr_split{curr_animal}{1,1};
curr_act_sort = cell2mat(arrayfun(@(x) curr_act(x,sort_idx(x,:)),1:length(sort_idx),'uni',false)');

plot(nanmean(curr_move_sort(:,1:end-1)),nanmean(curr_act_sort(:,1:end-1)),'.k');

curr_move = epoch_movement_corr_split{curr_animal}{2,2};
[curr_move_sort,sort_idx] = sort(curr_move,2);

curr_act = epoch_activity_corr_split{curr_animal}{2,2};
curr_act_sort = cell2mat(arrayfun(@(x) curr_act(x,sort_idx(x,:)),1:length(sort_idx),'uni',false)');

plot(nanmean(curr_move_sort(:,1:end-1)),nanmean(curr_act_sort(:,1:end-1)),'.r');

%% Reviewer 2 point 2.5: 6B/C alternate

clearvars -except data analysis classified_rois

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

% Make template movement from median of end sessions movements
template_movement = cellfun(@(move) nanmedian(move,2), ...
    epoch_movement_cat(:,end),'uni',false);

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);
movement_template_corr_pair_all = arrayfun(@(x) ...
    bsxfun(@max,vertcat(movement_template_corr{x,:}),vertcat(movement_template_corr{x,:})')./2, ...
    1:size(movement_template_corr,1),'uni',false);
movement_template_corr_pair = cellfun(@(allcorr,act1) ...
    {AP_itril(allcorr(1:size(act1,2),1:size(act1,2)),-1), ...
    AP_itril(allcorr(size(act1,2)+1:end,size(act1,2):end),-1), ...
    reshape(allcorr(1:size(act1,2),size(act1,2)+1:end),[],1)}, ...
    movement_template_corr_pair_all,epoch_activity_cat_reshape(:,1)','uni',false);
movement_template_corr_pair = vertcat(movement_template_corr_pair{:});

% only use trials which are high movement correlated
move_corr_lo = -0.5;
move_corr_hi = 0.5;
lo_move_corr_idx = cellfun(@(x) x <= move_corr_lo,epoch_movement_corr,'uni',false);
hi_move_corr_idx = cellfun(@(x) x >= move_corr_hi,epoch_movement_corr,'uni',false);

epoch_movement_corr_lo = cellfun(@(x,y) x(y),epoch_movement_corr,lo_move_corr_idx,'uni',false);
epoch_activity_corr_lo = cellfun(@(x,y) x(y),epoch_activity_corr,lo_move_corr_idx,'uni',false);
movement_template_corr_pair_lo = cellfun(@(x,y) x(y),movement_template_corr_pair,lo_move_corr_idx,'uni',false);

epoch_movement_corr_hi = cellfun(@(x,y) x(y),epoch_movement_corr,hi_move_corr_idx,'uni',false);
epoch_activity_corr_hi = cellfun(@(x,y) x(y),epoch_activity_corr,hi_move_corr_idx,'uni',false);
movement_template_corr_pair_hi = cellfun(@(x,y) x(y),movement_template_corr_pair,hi_move_corr_idx,'uni',false);

% Bin trials by correlation to template
corr_bin_use = [-Inf 0 Inf];
n_bins = length(corr_bin_use)-1;
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2];

[~,templatecorr_lo_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr_pair_lo,'uni',false);
[~,templatecorr_hi_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr_pair_hi,'uni',false);


figure;

% Low pairwise corr
subplot(1,2,1);

actcorr_bin_mean = cellfun(@(act,movebin) grpstats(act,movebin,'nanmean'), ...
    epoch_activity_corr_lo,templatecorr_lo_bin_idx,'uni',false);
actcorr_bin_sem = cellfun(@(act,movebin) grpstats(act,movebin,'sem'), ...
    epoch_activity_corr_lo,templatecorr_lo_bin_idx,'uni',false);

lo_actcorr_bin_mean_pad = cellfun(@(x) nan(n_bins,1),actcorr_bin_mean,'uni',false);
lo_actcorr_bin_sem_pad = cellfun(@(x) nan(n_bins,1),actcorr_bin_sem,'uni',false);

for i = 1:numel(actcorr_bin_mean)
    curr_used_bins = unique(templatecorr_lo_bin_idx{i});
    lo_actcorr_bin_mean_pad{i}(curr_used_bins) = actcorr_bin_mean{i};
    lo_actcorr_bin_sem_pad{i}(curr_used_bins) = actcorr_bin_sem{i};
end

epoch_combine_mean = nan(n_bins,3);
epoch_combine_sem = nan(n_bins,3);
for i = 1:3
    epoch_combine_mean(:,i) = nanmean(horzcat(lo_actcorr_bin_mean_pad{:,i}),2);
    epoch_combine_sem(:,i) = nanstd(horzcat(lo_actcorr_bin_mean_pad{:,i}),[],2)./ ...
        sqrt(sum(~isnan(horzcat(lo_actcorr_bin_mean_pad{:,i})),2));
end

errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);

legend({'Within early','Within late','Across early/late'});
xlabel('Correlation to template');
ylabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);
xlim([-1,1]);
set(gca,'XTick',-1:0.5:1);
title('Low pairwise correlation');

% High pairwise corr
subplot(1,2,2);

actcorr_bin_mean = cellfun(@(act,movebin) grpstats(act,movebin,'nanmean'), ...
    epoch_activity_corr_hi,templatecorr_hi_bin_idx,'uni',false);
actcorr_bin_sem = cellfun(@(act,movebin) grpstats(act,movebin,'sem'), ...
    epoch_activity_corr_hi,templatecorr_hi_bin_idx,'uni',false);

hi_actcorr_bin_mean_pad = cellfun(@(x) nan(n_bins,1),actcorr_bin_mean,'uni',false);
hi_actcorr_bin_sem_pad = cellfun(@(x) nan(n_bins,1),actcorr_bin_sem,'uni',false);

for i = 1:numel(actcorr_bin_mean)
    curr_used_bins = unique(templatecorr_hi_bin_idx{i});
    hi_actcorr_bin_mean_pad{i}(curr_used_bins) = actcorr_bin_mean{i};
    hi_actcorr_bin_sem_pad{i}(curr_used_bins) = actcorr_bin_sem{i};
end

epoch_combine_mean = nan(n_bins,3);
epoch_combine_sem = nan(n_bins,3);
for i = 1:3
    epoch_combine_mean(:,i) = nanmean(horzcat(hi_actcorr_bin_mean_pad{:,i}),2);
    epoch_combine_sem(:,i) = nanstd(horzcat(hi_actcorr_bin_mean_pad{:,i}),[],2)./ ...
        sqrt(sum(~isnan(horzcat(hi_actcorr_bin_mean_pad{:,i})),2));
end

errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);

legend({'Within early','Within late','Across early/late'});
xlabel('Correlation to template');
ylabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);
xlim([-1,1]);
set(gca,'XTick',-1:0.5:1);
title('High pairwise correlation');

% Statistics
lo_data_early = horzcat(lo_actcorr_bin_mean_pad{:,1});
lo_data_late = horzcat(lo_actcorr_bin_mean_pad{:,2});
lo_data_across = horzcat(lo_actcorr_bin_mean_pad{:,3});

hi_data_early = horzcat(hi_actcorr_bin_mean_pad{:,1});
hi_data_late = horzcat(hi_actcorr_bin_mean_pad{:,2});
hi_data_across = horzcat(hi_actcorr_bin_mean_pad{:,3});

p_lo_early = signrank(lo_data_early(1,:),lo_data_early(2,:));
p_lo_late = signrank(lo_data_late(1,:),lo_data_late(2,:));
p_lo_across = signrank(lo_data_across(1,:),lo_data_across(2,:));

p_hi_early = signrank(hi_data_early(1,:),hi_data_early(2,:));
p_hi_late = signrank(hi_data_late(1,:),hi_data_late(2,:));
p_hi_across = signrank(hi_data_across(1,:),hi_data_across(2,:));


%% Reviewer 3 point 2: example thresholding

use_animal = 7;
use_day = 1;
use_rois = 1:10;

figure; hold on;
AP_stackplot(data(use_animal).im(use_day).roi_trace_df(use_rois,:)',[],10,[],[0.5,0.5,0.5])
AP_stackplot(data(use_animal).im(use_day).roi_trace_thresh(use_rois,:)',[],10,[],'k')


%% Reviewer 3 point 2: 6A soft max norm

clearvars -except data analysis classified_rois

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

% REVIEWER REQUEST: soft-normalize activity by max over days
soft_norm = 0.25;
roi_max_withinday = cellfun(@(x) cellfun(@(x) max(max(x,[],1),[],3),x,'uni',false),epoch_activity,'uni',false);
roi_max_acrossday_epoch = cellfun(@(x) max(vertcat(x{:}),[],1),roi_max_withinday,'uni',false);
roi_max_acrossday_all = cell(size(roi_max_acrossday_epoch,1),length(session_epochs));
for curr_animal = 1:length(data)
    curr_max = max(cat(1,roi_max_acrossday_epoch{curr_animal,:}),[],1);
    for curr_epoch = 1:length(session_epochs)
        roi_max_acrossday_all{curr_animal,curr_epoch} = ...
            repmat({curr_max},1,length(session_epochs{curr_epoch}));
    end
end
epoch_activity_norm = cellfun(@(act,max) cellfun(@(act,max) bsxfun(@rdivide, ...
    act,max+soft_norm),act,max,'uni',false),epoch_activity,roi_max_acrossday_all,'uni',false);

% Get trial-by-trial activity/movement correlation

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity_norm(~any(empty_data,2),:),'uni',false);
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
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use)-1;
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
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
end

figure; hold on;
errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);
legend({'Within early','Within late','Across early/late'});
xlabel('Pairwise movement correlation');
ylabel('Pairwise activity correlation');
legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
    ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);

% Plot histogram of maxes just to show that this worked
max_raw = max(cell2mat(cellfun(@(x) max(max(cat(3,x{:}),[],1),[],3)', ...
    epoch_activity,'uni',false)),[],2);
max_norm = max(cell2mat(cellfun(@(x) max(max(cat(3,x{:}),[],1),[],3)', ...
    epoch_activity_norm,'uni',false)),[],2);
figure;
subplot(2,1,1);
hist(max_raw(max_raw ~= 0),100,'FaceColor','k');
xlabel('Max value')
ylabel('Number of cells')
title('Non-normalized')
subplot(2,1,2);
hist(max_norm(max_norm ~= 0),100,'FaceColor','k');
xlabel('Max value');
ylabel('Number of cells')
title('Normalized');

% New stats (reviewer 3)
early_data = horzcat(actcorr_bin_mean_pad{:,1});
early_data_use = early_data(1:end-1,:);
late_data = horzcat(actcorr_bin_mean_pad{:,2});
late_data_use = late_data(1:end-1,:);
across_data = horzcat(actcorr_bin_mean_pad{:,3});
across_data_use = across_data(1:end-1,:);
% (slope difference)
curr_fit = arrayfun(@(x) polyfit(1:n_bins,early_data_use(:,x)',1),1:length(data),'uni',false);
early_slope = cellfun(@(x) x(1),curr_fit);

curr_fit = arrayfun(@(x) polyfit(1:n_bins,late_data_use(:,x)',1),1:length(data),'uni',false);
late_slope = cellfun(@(x) x(1),curr_fit);

curr_fit = arrayfun(@(x) polyfit(1:n_bins,across_data_use(:,x)',1),1:length(data),'uni',false);
across_slope = cellfun(@(x) x(1),curr_fit);

p_slope_earlylate = signrank(early_slope,late_slope);
p_slope_earlyacross = signrank(early_slope,across_slope);
p_slope_lateacross = signrank(late_slope,across_slope);
% (plot segment difference)
p_negcorr = signrank(reshape(early_data_use(1:3,:),[],1),reshape(late_data_use(1:3,:),[],1));
p_poscorr = signrank(reshape(early_data_use(4:6,:),[],1),reshape(late_data_use(4:6,:),[],1));

%% Reviewer 3 additional point 1: 6a control for time and movement number

n_trials = nan(length(data),14);
for curr_animal = 1:length(data)
    for curr_day = 1:14
        rewarded_trials = sum(cellfun(@(x) ~isempty(x.states.reward), ...
            data(curr_animal).bhv(curr_day).bhv_times));
        n_trials(curr_animal,curr_day) = rewarded_trials;
    end
end

figure;hold on;
plot(n_trials');
xlabel('Day');
ylabel('Number of rewarded trials');


% Get 6A just using the first minimum number of movements
clearvars -except data analysis classified_rois

session_epochs = {[1,2,3,4],[11,12,13,14]};

% Get activity and movements (and time in session)
epoch_movement = cell(length(data),length(session_epochs));
epoch_activity = cell(length(data),length(session_epochs));
movement_time_in_session = cell(length(data),length(session_epochs));

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
            
            %             %%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%
            %             % USE ONLY MOMEMENTS FROM FIRST PART OF SESSION
            %             % (this was adopted at a later stage)
            %             min_frames = min(cellfun(@(x) size(x,2),{data(curr_animal).im.roi_trace_df}));
            %
            %             edge_movements = (use_movement_start_times + (use_time-1) > ...
            %                 length(lever_force_resample)) | ...
            %                 (use_movement_start_frames + (use_frames-1) > ...
            %                 min_frames);
            %
            %             %%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%
            
            movement_time_in_session{curr_animal,curr_epoch}{curr_session_idx} = ...
                use_movement_start_times;
            
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

% %%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%
% % ONLY USE THE MINIMUM NUMBER OF MOVEMENTS PER SESSION
% (now incorporated later)
% 
% min_n_move = min(cell2mat(cellfun(@(x) cellfun(@(x) size(x,2),x),epoch_movement,'uni',false)),[],2);
% epoch_movement = cellfun(@(x,n) cellfun(@(x) x(:,1:n),x,'uni',false),epoch_movement, ...
%     num2cell(repmat(min_n_move,1,2)),'uni',false);
% epoch_activity = cellfun(@(x,n) cellfun(@(x) x(:,:,1:n),x,'uni',false),epoch_activity, ...
%     num2cell(repmat(min_n_move,1,2)),'uni',false);
% 
% %%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%

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
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use)-1;
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),epoch_movement_corr,'uni',false);

%%%%% 3 verions: all movements, minimum time movements, minimum movements

all_movements = cellfun(@(x) true(size(x)),movecorr_bin_idx,'uni',false);

min_time = 1000*arrayfun(@(x) min(cellfun(@(x) size(x,2),{data(x).im.roi_trace_df})),1:length(data))/28;
under_time = arrayfun(@(x) cell2mat(cellfun(@(movetime) vertcat(movetime{:}), ...
    movement_time_in_session(x,:)','uni',false)) <= min_time(x),1:length(data),'uni',false)';
under_time_grid = cellfun(@(x) bsxfun(@and,x,x'),under_time,'uni',false);
min_time_movements = cellfun(@(undertime,moves) ...
    {AP_itril(undertime(1:size(moves,2),1:size(moves,2)),-1), ...
    AP_itril(undertime(size(moves,2)+1:end,size(moves,2):end),-1), ...
    reshape(undertime(1:size(moves,2),size(moves,2)+1:end),[],1)}, ...
    under_time_grid,epoch_movement_cat(:,1),'uni',false);
min_time_movements = vertcat(min_time_movements{:});

min_n_move = min(cell2mat(cellfun(@(x) cellfun(@(x) size(x,2),x),epoch_movement,'uni',false)),[],2);
n_move = cellfun(@(x) cellfun(@(x) [1:size(x,2)]',x,'uni',false),epoch_movement,'uni',false);
under_n_move = arrayfun(@(x) cell2mat(cellfun(@(movetime) vertcat(movetime{:}), ...
    n_move(x,:)','uni',false)) <= min_n_move(x),1:length(data),'uni',false)';
under_n_move_grid = cellfun(@(x) bsxfun(@and,x,x'),under_n_move,'uni',false);
min_num_movements = cellfun(@(undertime,moves) ...
    {AP_itril(undertime(1:size(moves,2),1:size(moves,2)),-1), ...
    AP_itril(undertime(size(moves,2)+1:end,size(moves,2):end),-1), ...
    reshape(undertime(1:size(moves,2),size(moves,2)+1:end),[],1)}, ...
    under_n_move_grid,epoch_movement_cat(:,1),'uni',false);
min_num_movements = vertcat(min_num_movements{:});

use_moves = {all_movements,min_time_movements,min_num_movements};
use_moves_label = {'All','Minimum time','Minimum number'};
% stats
p_slope_earlylate = nan(length(use_moves),1);
p_slope_lateacross = nan(length(use_moves),1);
p_negcorr = nan(length(use_moves),1);
figure;
for curr_use_moves = 1:3;
    
    curr_move_idx = use_moves{curr_use_moves};
    
    % Get mean/sem of activity correlations in each movement correlation bin
    actcorr_bin_mean = cellfun(@(act,movebin,move_idx) grpstats(act(move_idx),movebin(move_idx),'nanmean'), ...
        epoch_activity_corr,movecorr_bin_idx,curr_move_idx,'uni',false);
    actcorr_bin_sem = cellfun(@(act,movebin,move_idx) grpstats(act(move_idx),movebin(move_idx),'sem'), ...
        epoch_activity_corr,movecorr_bin_idx,curr_move_idx,'uni',false);
    
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
    end
    
    subplot(1,3,curr_use_moves); hold on;
    errorbar(repmat(corr_bin_plot',[1,3,1]),epoch_combine_mean,epoch_combine_sem,'linewidth',2);
    legend({'Within early','Within late','Across early/late'});
    xlabel('Pairwise movement correlation');
    ylabel('Pairwise activity correlation');
    legend([cellfun(@(x) ['Within sessions ' num2str(x)],session_epochs,'uni',false), ...
        ['Across sessions ' num2str(session_epochs{1}) ', ' num2str(session_epochs{end})]]);
    title(use_moves_label{curr_use_moves});
    
    % New stats (reviewer 3)
    early_data = horzcat(actcorr_bin_mean_pad{:,1});
    early_data_use = early_data(1:end-1,:);
    late_data = horzcat(actcorr_bin_mean_pad{:,2});
    late_data_use = late_data(1:end-1,:);
    across_data = horzcat(actcorr_bin_mean_pad{:,3});
    across_data_use = across_data(1:end-1,:);
    % (slope difference)
    curr_fit = arrayfun(@(x) polyfit(1:n_bins,early_data_use(:,x)',1),1:length(data),'uni',false);
    early_slope = cellfun(@(x) x(1),curr_fit);
    
    curr_fit = arrayfun(@(x) polyfit(1:n_bins,late_data_use(:,x)',1),1:length(data),'uni',false);
    late_slope = cellfun(@(x) x(1),curr_fit);
    
    curr_fit = arrayfun(@(x) polyfit(1:n_bins,across_data_use(:,x)',1),1:length(data),'uni',false);
    across_slope = cellfun(@(x) x(1),curr_fit);
    
    p_slope_earlylate(curr_use_moves) = signrank(early_slope,late_slope);
    p_slope_earlyacross = signrank(early_slope,across_slope);
    p_slope_lateacross(curr_use_moves) = signrank(late_slope,across_slope);
    % (plot segment difference)
    p_negcorr(curr_use_moves) = signrank(reshape(early_data_use(1:3,:),[],1),reshape(late_data_use(1:3,:),[],1));
    p_poscorr = signrank(reshape(early_data_use(4:6,:),[],1),reshape(late_data_use(4:6,:),[],1));
    
end




%% Reviewer 3 additional point 2: 6a middle, smooth/jump representation

clearvars -except data analysis classified_rois

session_epochs = {[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14]};

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

% Movement/activity correlation within/across all epochs
% (for this on: only use animals with all epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

epoch_movement_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_movement_cat{x,:})), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
n_trials = mat2cell(cellfun(@(x) size(x,2),epoch_activity_cat_reshape), ...
    ones(size(epoch_activity_cat_reshape,1),1),length(session_epochs));
epoch_movement_corr_split = cellfun(@(move_cor,n_trials) ...
    mat2cell(move_cor,n_trials,n_trials),epoch_movement_corr_all,n_trials','uni',false);
epoch_movement_corr_reshape = epoch_movement_corr_split;
for curr_animal = 1:size(epoch_movement_cat,1);
    epoch_movement_corr_reshape{curr_animal}(logical(eye(length(session_epochs)))) = ...
        cellfun(@(x) AP_itril(x,-1), ...
        epoch_movement_corr_split{curr_animal}(logical(eye(length(session_epochs)))),'uni',false);
    epoch_movement_corr_reshape{curr_animal}(~logical(eye(length(session_epochs)))) = ...
        cellfun(@(x) reshape(x,[],1), ...
        epoch_movement_corr_split{curr_animal}(~logical(eye(length(session_epochs)))),'uni',false);  
end
epoch_movement_corr_reshape = cat(3,epoch_movement_corr_reshape{:});

epoch_activity_corr_all = arrayfun(@(x) ...
    corrcoef(horzcat(epoch_activity_cat_reshape{x,:}),'rows','pairwise'), ...
    1:size(epoch_activity_cat_reshape,1),'uni',false);
n_trials = mat2cell(cellfun(@(x) size(x,2),epoch_activity_cat_reshape), ...
    ones(size(epoch_activity_cat_reshape,1),1),length(session_epochs));
epoch_activity_corr_split = cellfun(@(act_cor,n_trials) ...
    mat2cell(act_cor,n_trials,n_trials),epoch_activity_corr_all,n_trials','uni',false);
epoch_activity_corr_reshape = epoch_activity_corr_split;
for curr_animal = 1:size(epoch_movement_cat,1);
    epoch_activity_corr_reshape{curr_animal}(logical(eye(length(session_epochs)))) = ...
        cellfun(@(x) AP_itril(x,-1), ...
        epoch_activity_corr_split{curr_animal}(logical(eye(length(session_epochs)))),'uni',false);
    epoch_activity_corr_reshape{curr_animal}(~logical(eye(length(session_epochs)))) = ...
        cellfun(@(x) reshape(x,[],1), ...
        epoch_activity_corr_split{curr_animal}(~logical(eye(length(session_epochs)))),'uni',false);  
end
epoch_activity_corr_reshape = cat(3,epoch_activity_corr_reshape{:});

% Pick out the comparisons to plot
use_plots = [1,length(session_epochs)];
epoch_movement_corr_plot = permute(epoch_movement_corr_reshape(use_plots,:,:),[3,2,1]);
epoch_activity_corr_plot = permute(epoch_activity_corr_reshape(use_plots,:,:),[3,2,1]);
[~,n_lines,n_plots] = size(epoch_movement_corr_plot);

% Bin trials by movement correlation
corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin = corr_bin_use;
corr_bin(1) = -1;
corr_bin(end) = 1;
corr_bin_plot = [corr_bin(1:end-1) + diff(corr_bin)/2 NaN];

[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),epoch_movement_corr_plot,'uni',false);

% Get mean/sem of activity correlations in each movement correlation bin
actcorr_bin_mean = cellfun(@(act,movebin) grpstats(act,movebin,'nanmean'), ...
    epoch_activity_corr_plot,movecorr_bin_idx,'uni',false);
actcorr_bin_sem = cellfun(@(act,movebin) grpstats(act,movebin,'sem'), ...
    epoch_activity_corr_plot,movecorr_bin_idx,'uni',false);

actcorr_bin_mean_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_mean,'uni',false);
actcorr_bin_sem_pad = cellfun(@(x) nan(length(corr_bin),1),actcorr_bin_sem,'uni',false);

for i = 1:numel(actcorr_bin_mean)
    curr_used_bins = unique(movecorr_bin_idx{i});
    actcorr_bin_mean_pad{i}(curr_used_bins) = actcorr_bin_mean{i};
    actcorr_bin_sem_pad{i}(curr_used_bins) = actcorr_bin_sem{i};
end

epoch_combine_mean = nan(length(corr_bin),n_lines,n_plots);
epoch_combine_sem = nan(length(corr_bin),n_lines,n_plots);
for curr_plot = 1:n_plots
    for curr_line = 1:n_lines
        epoch_combine_mean(:,curr_line,curr_plot) = ...
            nanmean(horzcat(actcorr_bin_mean_pad{:,curr_line,curr_plot}),2);
        epoch_combine_sem(:,curr_line,curr_plot) = ...
            nanstd(horzcat(actcorr_bin_mean_pad{:,curr_line,curr_plot}),[],2)./ ...
            sqrt(sum(~isnan(horzcat(actcorr_bin_mean_pad{:,curr_line,curr_plot})),2));
    end
end

figure;
for curr_plot = 1:n_plots
    subplot(1,n_plots,curr_plot); hold on;
    set(gca,'ColorOrder',copper(n_lines));
    errorbar(repmat(corr_bin_plot',[1,n_lines,1]), ...
        epoch_combine_mean(:,:,curr_plot),epoch_combine_sem(:,:,curr_plot),'linewidth',2);
    xlabel('Pairwise movement correlation');
    ylabel('Pairwise activity correlation');
    legend(cellfun(@(x) ...
        ['Across sessions ' num2str(session_epochs{use_plots(curr_plot)}) ', ' ...
        num2str(x)],session_epochs,'uni',false));
end

% Statistics

% Time within each group
corr_r = nan(3,1);
corr_p = nan(3,1);
anova_p = nan(3,1);
for i = 1:3
    curr_data = horzcat(actcorr_bin_mean_pad{:,i});
    curr_bins = repmat(transpose(1:length(corr_bin)),1,size(curr_data,2));
    
    curr_data = curr_data(1:end-1,:);
    curr_bins = curr_bins(1:end-1,:);
    
    %[r,p] = corrcoef(curr_bins,zscore(curr_data,[],1));
    [r,p] = corrcoef(curr_bins,curr_data);
    corr_r(i) = r(2);
    corr_p(i) = p(2);
    
    anova_p(i) = anova1(zscore(curr_data(3:end,:),[],1)',[],'off');
end

% Across early/late
early_data = horzcat(actcorr_bin_mean_pad{:,1});
early_data_use = early_data(1:end-1,:);
late_data = horzcat(actcorr_bin_mean_pad{:,2});
late_data_use = late_data(1:end-1,:);

p_negcorr = signrank(reshape(early_data_use(1:3,:),[],1),reshape(late_data_use(1:3,:),[],1));
p_poscorr = signrank(reshape(early_data_use(4:6,:),[],1),reshape(late_data_use(4:6,:),[],1));

%% EXTRA: 6A using only cued-rewarded movements


stages = {[1:4],[11:14]};
framerate = 28;

activity_corr_stages_mean = cell(length(stages),length(stages),length(analysis));
activity_corr_stages_sem = cell(length(stages),length(stages),length(analysis));
for curr_animal = 1:length(analysis);
    
    movement_start_time = 1001;
    movement_use_time = 2000;
    
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
    corr_bin_use = [-Inf -0.5:0.25:0.5 Inf];
    n_bins = length(corr_bin_use);
    corr_bin = corr_bin_use;
    corr_bin(1) = -1;
    corr_bin(end) = 1;
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
        
        curr_sem = nanstd(curr_mean,[],1)./sqrt(sum(~isnan(curr_mean)));
        activity_corr_stages_sem_animals{stage_1,stage_2} = curr_sem;
        
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



%% EXTRA: 6B using cued-rewarded movement as template


clearvars -except data analysis classified_rois

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

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

% Make template get trial-by-trial activity/movement correlation
corr_bin_use = [-Inf -0.5 0.5 Inf];
n_bins = length(corr_bin_use);
corr_bin_valid = corr_bin_use;
corr_bin_valid(1) = -1;
corr_bin_valid(end) = 1;
corr_bin_plot = [corr_bin_valid(1:end-1) + diff(corr_bin_valid)/2 NaN];

all_act_corr = cell(size(epoch_movement_cat,1),length(session_epochs),n_bins);
epoch_combine_mean = nan(n_bins,length(session_epochs));
epoch_combine_sem = nan(n_bins,length(session_epochs));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DIFFERENT %%%%%%%%%%%%%%%%%%%

% OLD
% % Make template movement from median of end sessions movements
% template_movement = cellfun(@(move) nanmedian(move,2), ...
%     epoch_movement_cat(:,end),'uni',false);

movement_start_time = 1001;
movement_use_time = 2000;

% Trial-by-trial movement correlation
template_movement = cell(length(data),1);
for curr_animal = 1:length(data)
    
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time-1),movement(trials & ...
        cellfun(@(y) ~isempty(y),movement)),'uni',false)'), ...
        {analysis(curr_animal).lever(end).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(end).cued_movement_trials},'uni',false);   
    
    template_movement{curr_animal} = nanmean(horzcat(use_movement{:}),2);
    
end

template_movement = cellfun(@(x) shake(x,1),template_movement,'uni',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Bin trials by movement correlation
[~,movecorr_bin_idx] = cellfun(@(x) histc(x,corr_bin_use),movement_template_corr,'uni',false);

actcorr_bin_mean = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
actcorr_bin_sem = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
movecorr_bin_mean = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
movecorr_bin_sem = nan(n_bins,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        for curr_bin = 1:n_bins
            
            use_trials =  ...
                find(movecorr_bin_idx{curr_animal,curr_epoch} == curr_bin);
            
%             % Skip if there aren't enough trials
%             if length(use_trials) < 5
%                 continue
%             end
            
            % Get pairwise correlation of activity within bins
            curr_act_corr = AP_itril(corrcoef( ...
                epoch_activity_cat_reshape{curr_animal,curr_epoch}(:,use_trials),'rows','pairwise'),-1);
            actcorr_bin_mean(curr_bin,curr_epoch,curr_animal) = nanmean(curr_act_corr);
            actcorr_bin_sem(curr_bin,curr_epoch,curr_animal) = nanstd(curr_act_corr)./sqrt(sum(~isnan(curr_act_corr)));
            
            % Get pairwise correlation of movement within bins, for use maybe in later sanity check
            curr_move_corr = AP_itril(corrcoef( ...
                epoch_movement_cat{curr_animal,curr_epoch}(:,use_trials),'rows','pairwise'),-1);
            movecorr_bin_mean(curr_bin,curr_epoch,curr_animal) = nanmean(curr_move_corr);
            movecorr_bin_sem(curr_bin,curr_epoch,curr_animal) = nanstd(curr_move_corr)./sqrt(sum(~isnan(curr_move_corr)));
            
            % Store all the activity correlations (total grouping?)
            all_act_corr{curr_animal,curr_epoch,curr_bin} = curr_act_corr;
        end
    end
end

epoch_combine_mean(:,:) = nanmean(actcorr_bin_mean,3);
epoch_combine_sem(:,:) = nanstd(actcorr_bin_mean,[],3)./sqrt(sum(~isnan(actcorr_bin_mean),3));

% few_animals = sum(~isnan(actcorr_bin_mean),3) < 2;
% epoch_combine_mean(few_animals) = NaN;
% epoch_combine_sem(few_animals) = NaN;

figure; 
plot_bins = [1,3];
% Plot pairwise movement correlations
subplot(2,1,1);hold on;
errorbar(repmat(corr_bin_plot(plot_bins)',[1,2]),nanmean(movecorr_bin_mean(plot_bins,:,:),3), ...
    nanmean(movecorr_bin_sem(plot_bins,:,:),3),'linewidth',2);
xlabel('Movement template correlation');
ylabel('Pairwise movement correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));

% Plot pairwise activity correlations
subplot(2,1,2);hold on;
errorbar(repmat(corr_bin_plot(plot_bins)',[1,2]),nanmean(epoch_combine_mean(plot_bins,:,:),3), ...
    nanmean(epoch_combine_sem(plot_bins,:,:),3),'linewidth',2);
xlabel('Movement template correlation');
ylabel('Pairwise activity correlation');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));


%% EXTRA: 6C using cued-rewarded movements as template


clearvars -except data analysis classified_rois

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

% Movement/activity correlation (within epoch 1, within epoch 2, across)
% (for this on: only use animals with both sets of epochs)
epoch_movement_cat = cellfun(@(x) horzcat(x{:}),epoch_movement(~any(empty_data,2),:),'uni',false);
epoch_activity_cat = cellfun(@(x) cat(3,x{:}),epoch_activity(~any(empty_data,2),:),'uni',false);
epoch_activity_cat_reshape = cellfun(@(x) reshape(x,[],size(x,3)),epoch_activity_cat,'uni',false);

% Make template movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DIFFERENT %%%%%%%%%%%%%%%%%%%

% OLD
% % Make template movement from median of end sessions movements
% template_movement = cellfun(@(move) nanmedian(move,2), ...
%     epoch_movement_cat(:,end),'uni',false);

movement_start_time = 1001;
movement_use_time = 2000;

% Trial-by-trial movement correlation
template_movement = cell(length(data),1);
for curr_animal = 1:length(data)
    
    use_movement = cellfun(@(movement,trials) cell2mat(cellfun(@(x) ...
        x(movement_start_time:movement_start_time+movement_use_time-1),movement(trials & ...
        cellfun(@(y) ~isempty(y),movement)),'uni',false)'), ...
        {analysis(curr_animal).lever(end).rewarded_movement_fixtime}, ...
        {analysis(curr_animal).lever(end).cued_movement_trials},'uni',false);
    
    template_movement{curr_animal} = nanmean(horzcat(use_movement{:}),2);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get correlation of all trials with template
movement_template_corr_all = cellfun(@(move,template) ...
    corrcoef([template,move],'rows','pairwise'), ...
    epoch_movement_cat,repmat(template_movement,1,2),'uni',false);
movement_template_corr = cellfun(@(x) x(2:end,1),movement_template_corr_all,'uni',false);

% Get pairwise correlations of all trials
movement_pairwise_corr = cellfun(@(x) corrcoef(x),epoch_movement_cat,'uni',false);

% Get differences in pairwise correlation with template
movement_template_corr_diff = cellfun(@(x) abs(bsxfun(@minus,x,x')), ...
    movement_template_corr,'uni',false);

% Select trials by pairwise correlation AND template correlation
movecorr_cutoff = -0.4;
movetemp_cutoff = 1.2;
ortho_movetemp_cutoff = 0.2;

template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movetempdiff > movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

ortho_template_axis_use = cellfun(@(movecorr,movetempdiff) ...
    tril(movecorr < movecorr_cutoff & ...
    movetempdiff < ortho_movetemp_cutoff,-1), ...
    movement_pairwise_corr,movement_template_corr_diff,'uni',false);

movecorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
actcorr_movedim = nan(2,length(session_epochs),size(epoch_movement_cat,1));
for curr_animal = 1:size(epoch_movement_cat,1)
    for curr_epoch = 1:length(session_epochs);
        
        template_dim_pairs = template_axis_use{curr_animal,curr_epoch};
        ortho_template_dim_pairs = ortho_template_axis_use{curr_animal,curr_epoch};
        
        curr_move_corr_grid = corrcoef( ...
            epoch_movement_cat{curr_animal,curr_epoch},'rows','pairwise');
        curr_act_corr_grid = corrcoef( ...
            epoch_activity_cat_reshape{curr_animal,curr_epoch},'rows','pairwise');      
        
        % Get pairwise movement correlations along template edges
        move_corr_templatedim_mean = nanmean(curr_move_corr_grid(template_dim_pairs));
        
        % Get pairwise movement correlations in selected middle bin trials
        move_corr_orthotemplatedim_mean = nanmean(curr_move_corr_grid(ortho_template_dim_pairs));
        
        % Get pairwise activity correlations along template edges
        act_corr_templatedim_mean = nanmean(curr_act_corr_grid(template_dim_pairs));
        
        % Get pairwise activity correlations in selected middle bin trials
        act_corr_orthotemplatedim_mean = nanmean(curr_act_corr_grid(ortho_template_dim_pairs));
        
        % Store mean pairwise correlations
        movecorr_movedim(:,curr_epoch,curr_animal) = ...
            [move_corr_templatedim_mean;move_corr_orthotemplatedim_mean];
        actcorr_movedim(:,curr_epoch,curr_animal) = ...
            [act_corr_templatedim_mean;act_corr_orthotemplatedim_mean];
        
    end
end

movecorr_movedim_mean = nanmean(movecorr_movedim,3);
movecorr_movedim_sem = nanstd(movecorr_movedim,[],3)./sqrt(sum(~isnan(movecorr_movedim),3));

actcorr_movedim_mean = nanmean(actcorr_movedim,3);
actcorr_movedim_sem = nanstd(actcorr_movedim,[],3)./sqrt(sum(~isnan(actcorr_movedim),3));

figure;

subplot(2,1,1);
errorbar(movecorr_movedim_mean,movecorr_movedim_sem,'linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Template dimension','Ortho-template dimension'});
ylabel('Pairwise movement correlations');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));

subplot(2,1,2);
errorbar(actcorr_movedim_mean,actcorr_movedim_sem,'linewidth',2)
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Template dimension','Ortho-template dimension'});
ylabel('Pairwise activity correlations');
legend(cellfun(@(x) ['Sessions ' num2str(x)],session_epochs,'uni',false));










