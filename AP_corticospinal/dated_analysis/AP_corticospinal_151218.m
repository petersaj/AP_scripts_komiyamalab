%% Autocorrelation of M activity/Q activity/movement over days

maxlag = 100;

m_autocorr_all = cell(length(data),1);
q_autocorr_all = cell(length(data),1);
move_autocorr_all = cell(length(data),1);

for curr_animal = 1:length(data);
    
    m_autocorr = nan(14,maxlag*2+1);
    q_autocorr = nan(14,maxlag*2+1);
    move_autocorr = nan(14,maxlag*2+1);
    
    for curr_day = 1:length(data(curr_animal).im);
        
        % Auto-correlation of calcium event onsets
        
        m = find(classified_rois(curr_animal).movement(:,curr_day))';
        q = find(classified_rois(curr_animal).quiescent(:,curr_day))';
        
        leeway = ones(1,3);
        
        m_onsets = conv2(double(diff(data(curr_animal).im(curr_day). ...
            roi_trace_thresh(m,:) > 0,[],2) == 1),leeway,'same') > 0;
        
        q_onsets = conv2(double(diff(data(curr_animal).im(curr_day). ...
            roi_trace_thresh(q,:) > 0,[],2) == 1),leeway,'same') > 0;
        
        curr_m_autocorr = cell(length(m),1);
        for i = 1:size(m_onsets,1)
            curr_m_autocorr{i} = xcorr(m_onsets(i,:),100,'coeff');
        end
        
        curr_q_autocorr = cell(length(q),1);
        for i = 1:size(q_onsets,1)
            curr_q_autocorr{i} = xcorr(q_onsets(i,:),100,'coeff');
        end
        
        m_autocorr(curr_day,:) = nanmean(vertcat(curr_m_autocorr{:}));
        q_autocorr(curr_day,:) = nanmean(vertcat(curr_q_autocorr{:}));
        
        % Auto-correlation of movement onsets
        move_onsets = conv2(double(diff(analysis(curr_animal).lever(curr_day). ...
            lever_move_frames > 0) == 1),leeway','same') > 0;
        move_autocorr(curr_day,:) = xcorr(move_onsets,100,'coeff');
        
    end
    
    m_autocorr_all{curr_animal} = m_autocorr;
    q_autocorr_all{curr_animal} = q_autocorr;
    move_autocorr_all{curr_animal} = move_autocorr;
    
    disp(curr_animal);
end


% Find the max values and times
m_autocorr_max = nan(14,length(data));
m_autocorr_max_idx = nan(14,length(data));
for curr_animal = 1:length(data);
    for curr_session = 1:length(data(curr_animal).im)
        curr_autocorr = m_autocorr_all{curr_animal};
        curr_autocorr(:,maxlag+1-length(leeway):maxlag+1+length(leeway)) = NaN;
        
        [m_autocorr_max(:,curr_animal), ...
            m_autocorr_max_idx(:,curr_animal)] = ...
            max(curr_autocorr,[],2);
    end
end

% Find the max values and times
q_autocorr_max = nan(14,length(data));
q_autocorr_max_idx = nan(14,length(data));
for curr_animal = 1:length(data);
    for curr_session = 1:length(data(curr_animal).im)
        curr_autocorr = q_autocorr_all{curr_animal};
        curr_autocorr(:,maxlag+1-length(leeway):maxlag+1+length(leeway)) = NaN;
        
        [q_autocorr_max(:,curr_animal), ...
            q_autocorr_max_idx(:,curr_animal)] = ...
            max(curr_autocorr,[],2);
    end
end


% there is periodicity of a little under 1 sec, doesn't change with
% learning, is relatively the same across cell types, isn't directly
% related to movement periodicity?

% did this for L2/3 animals and is present there, I think that this is just
% an artifact of the size of the loess filter during event detection


%% Get distribution of activity during session v during movement

% have other options here to do: classified p-value, correlation

% % Correlate thresholded activity with binary lever movement
% activity_move_corr = cell(length(data),1);
% for curr_animal = 1:length(data);
%     n_days = length(data(curr_animal).im);
%     n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
%     
%     activity_move_corr{curr_animal} = nan(n_rois,n_days);
%     for curr_day = 1:length(data(curr_animal).im)
%         warning off
%         activity_move_corr{curr_animal}(:,curr_day) = ...
%             1-pdist2(data(curr_animal).im(curr_day).roi_trace_thresh, ...
%             analysis(curr_animal).lever(curr_day).lever_move_frames','correlation');
%         warning on
%     end
%     disp(curr_animal);
% end

% Get fraction of activity during session and during movement
activity_frac = cell(length(data),1);
activity_move_frac = cell(length(data),1);

activity_frac_shuff = cell(length(data),1);
activity_move_frac_shuff = cell(length(data),1);

for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    activity_frac{curr_animal} = nan(n_rois,n_days);
    activity_move_frac{curr_animal} = nan(n_rois,n_days);
    
    activity_frac_shuff{curr_animal} = nan(n_rois,n_days);
    activity_move_frac_shuff{curr_animal} = nan(n_rois,n_days);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        activity_frac{curr_animal}(:,curr_day) = ...
            nanmean(curr_act,2);
        
        warning off
        activity_move_frac{curr_animal}(:,curr_day) = ...
            (curr_act* ...
            analysis(curr_animal).lever(curr_day).lever_move_frames) ./ ...
            sum(curr_act,2);
        warning on
        
        % Activity-epoch shuffled data

        boundary_frames = arrayfun(@(x) find(diff([Inf,curr_act(x,:),Inf]) ~= 0),1:size(curr_act,1),'uni',false);
        curr_act_shuff = cell2mat(arrayfun(@(x) cell2mat(shake( ...
            mat2cell(curr_act(x,:),1,diff(boundary_frames{x})))),[1:size(curr_act,1)]','uni',false));  
        
        activity_frac_shuff{curr_animal}(:,curr_day) = ...
            nanmean(curr_act_shuff,2);
        
        warning off
        activity_move_frac_shuff{curr_animal}(:,curr_day) = ...
            (curr_act_shuff* ...
            analysis(curr_animal).lever(curr_day).lever_move_frames) ./ ...
            sum(curr_act_shuff,2);
        warning on
        
    end
    disp(curr_animal);
end

total_active_frac_v_movefrac = [ ...
    cell2mat(cellfun(@(x) x(:),activity_frac,'uni',false)), ...
    cell2mat(cellfun(@(x) x(:),activity_move_frac,'uni',false))];

total_active_frac_v_movefrac_shuff = [ ...
    cell2mat(cellfun(@(x) x(:),activity_frac_shuff,'uni',false)), ...
    cell2mat(cellfun(@(x) x(:),activity_move_frac_shuff,'uni',false))];

% Sanity check: fraction active should be the same in both conditions
if total_active_frac_v_movefrac(:,1) ~= total_active_frac_v_movefrac_shuff(:,1)
    error('Shuffled has different fraction activity')
end

% Plot density plot of early v. late differences
n_bins = 100;

%bin1_edges = linspace(min(total_active_frac_v_movefrac(~isinf(total_active_frac_v_movefrac(:,1)),1)), ...
%    max(total_active_frac_v_movefrac(:,1),[],1),n_bins);
bin1_edges = linspace(0,0.1,n_bins+1);
bin1_centers = bin1_edges(1:end-1) + diff(bin1_edges)/2;
bin1_label = round(bin1_centers*100)/100;
bin1_edges(1) = -Inf;
%bin1_edges(end) = Inf;

bin2_edges = linspace(0,1,n_bins+1);
bin2_centers = bin2_edges(1:end-1) + diff(bin2_edges)/2;
bin2_label = round(bin2_centers*100)/100;
bin2_edges(1) = -Inf;
bin2_edges(end) = Inf;

bin_edges = {bin1_edges,bin2_edges};

figure;

% Real
subplot(1,2,1);axis square;
reliability_density = hist3(total_active_frac_v_movefrac,'Edges',bin_edges);
h = fspecial('gaussian',20,2);
reliability_density_smooth = imfilter(reliability_density,h,'same');
imagesc(reliability_density_smooth);colormap(hot);
set(gca,'YDir','normal')
xlabel('Fraction active during movement')
ylabel('Fraction active during session')
xlim([0.5 n_bins]);
ylim([0.5 n_bins]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'XTick')));
c = colorbar;
ylabel(c,'Number of cells');
title('Real');

% Shuffled
subplot(1,2,2);axis square;
reliability_density_shuff = hist3(total_active_frac_v_movefrac_shuff,'Edges',bin_edges);
h = fspecial('gaussian',20,2);
reliability_density_shuff_smooth = imfilter(reliability_density_shuff,h,'same');
imagesc(reliability_density_shuff_smooth);colormap(hot);
set(gca,'YDir','normal')
xlabel('Fraction active during movement')
ylabel('Fraction active during session')
xlim([0.5 n_bins]);
ylim([0.5 n_bins]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'XTick')));
c = colorbar;
ylabel(c,'Number of cells');
title('Shuffled');

% White to red colormap
% Make blue-to-red colormap
curr_colormap = colormap;
c_vals = size(curr_colormap,1);
wr_colormap = [ones(c_vals,1),repmat(linspace(1,0,c_vals)',1,2)];
colormap(wr_colormap);

% Show difference between real and activity shuffled
figure;axis square;
reliability_density_diff = reliability_density - reliability_density_shuff;
reliability_density_diff_smooth = imfilter(reliability_density_diff,h,'same');
imagesc(reliability_density_diff_smooth);
set(gca,'YDir','normal')
xlabel('Fraction active during movement')
ylabel('Fraction active during session')
xlim([0.5 n_bins]);
ylim([0.5 n_bins]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'XTick')));
c = colorbar;
ylabel(c,'Number of cells');
title('Real - Shuffled');

% Make blue-to-red colormap
curr_colormap = colormap;
c_vals = size(curr_colormap,1);
bwr_colormap = [repmat(linspace(0,1,ceil(c_vals/2))',1,2),ones(ceil(c_vals/2),1); ...
    ones(ceil(c_vals/2),1),repmat(linspace(1,0,ceil(c_vals/2))',1,2)];
colormap(bwr_colormap);


% Plot the measurements ONLY for classified ROIs
total_classified = cell2mat(cellfun(@(x,y) reshape(x(~isnan(x)) | y(~isnan(y)),[],1), ...
    {classified_rois(:).movement},{classified_rois(:).quiescent},'uni',false)');
total_active_frac_v_movefrac_classified = ...
    total_active_frac_v_movefrac(total_classified,:);
total_active_frac_v_movefrac_nonclassified = ...
    total_active_frac_v_movefrac(~total_classified,:);

figure;

subplot(1,2,1);axis square;
reliability_density = hist3(total_active_frac_v_movefrac_classified,'Edges',bin_edges);
h = fspecial('gaussian',20,2);
reliability_density_smooth = imfilter(reliability_density,h,'same');
imagesc(reliability_density_smooth);colormap(hot);
set(gca,'YDir','normal')
xlabel('Fraction active during movement')
ylabel('Fraction active during session')
xlim([0.5 n_bins]);
ylim([0.5 n_bins]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'YTick')));
c = colorbar;
ylabel(c,'Number of cells');
title('Classified');

subplot(1,2,2);axis square;
reliability_density = hist3(total_active_frac_v_movefrac_nonclassified,'Edges',bin_edges);
h = fspecial('gaussian',20,2);
reliability_density_smooth = imfilter(reliability_density,h,'same');
imagesc(reliability_density_smooth);colormap(hot);
set(gca,'YDir','normal')
xlabel('Fraction active during movement')
ylabel('Fraction active during session')
xlim([0.5 n_bins]);
ylim([0.5 n_bins]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'XTick')));
c = colorbar;
ylabel(c,'Number of cells');
title('Nonclassified');

curr_colormap = colormap;
c_vals = size(curr_colormap,1);
wr_colormap = [ones(c_vals,1),repmat(linspace(1,0,c_vals)',1,2)];
colormap(wr_colormap);


% NORMALIZE THIS TO FRAC OF ALL ROIS?


%% Jump in Q cells: is this artifact of sampling, or really jump?

mq_act_u2q = cell(length(data),1);

for curr_animal = 1:length(data);
    curr_quiescentcells = any(classified_rois(curr_animal).quiescent,2);
    
    sessions = {[1,2],[3,4]};
    
    late_notearly_q = ~any(classified_rois(curr_animal).quiescent(:,sessions{1}),2) & ...
        any(classified_rois(curr_animal).quiescent(:,sessions{2}),2);
    
    use_cells = find(late_notearly_q);
    
    curr_mq_act_u2q = cellfun(@(x) nan(length(use_cells),3,length(x)),sessions,'uni',false);
        
    for i = 1:length(sessions)
        for j = 1:length(sessions{i})
            curr_session = sessions{i}(j);
            curr_mq_act_u2q{i}(:,1,j) = nanmean(data(curr_animal).im(curr_session). ...
                roi_trace_thresh(use_cells,analysis(curr_animal).lever(curr_session).lever_move_frames) > 0,2);
            curr_mq_act_u2q{i}(:,2,j) = nanmean(data(curr_animal).im(curr_session). ...
                roi_trace_thresh(use_cells,~analysis(curr_animal).lever(curr_session).lever_move_frames) > 0,2);
            curr_mq_act_u2q{i}(:,3,j) = nanmean(data(curr_animal).im(curr_session). ...
                roi_trace_thresh(use_cells,:) > 0,2);
        end
    end
    
    mq_act_u2q{curr_animal} = cellfun(@(x) nanmean(x,3),curr_mq_act_u2q,'uni',false);
   
end

mq_act_u2q_cat = cellfun(@(x) cell2mat(cellfun(@nanmean,x,'uni',false)'),mq_act_u2q,'uni',false);
mq_act_u2q_fullcat = cat(3,mq_act_u2q_cat{:});
mq_act_u2q_mean = nanmean(mq_act_u2q_fullcat,3);
mq_act_u2q_sem = nanstd(mq_act_u2q_fullcat,[],3)./sqrt(sum(~isnan(mq_act_u2q_fullcat),3));

figure; hold on;
bar(mq_act_u2q_mean'); colormap(gray)
errorb(mq_act_u2q_mean',mq_act_u2q_sem');
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'M','Q','Total'});
title('1:2 non-Q, 1:3 Q');
ylabel('Average activity');
legend({'Days 1-2','Days 3-4'});

% Total quiescent activity
total_act = cell(length(data),1);
for curr_animal = 1:length(data);
        
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_act = nan(n_rois,14,3);
    
    for i = 1:length(data(curr_animal).im)
        curr_act(:,i,1) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,2) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,~analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,3) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,:) > 0,2);
    end
       
    total_act{curr_animal} = permute(nanmean(curr_act,1),[3,2,1]);
   
end

total_act_cat = cat(3,total_act{:});
total_act_norm = bsxfun(@times,total_act_cat,1./total_act_cat(1,:,:));
total_act_mean = nanmean(total_act_norm,3);
total_act_sem = nanstd(total_act_norm,[],3)./sqrt(sum(~isnan(total_act_norm),3));

figure;
errorbar(total_act_mean',total_act_sem')
xlabel('Day');
ylabel('Fraction of activity (normalized to movement activity)')
legend({'Movement','Quiescent','Total'});


%% Q/M ROIs by inter-event interval distribution?

% Get inter-event-interval fano factor
ff = cell(length(data),1);
for curr_animal = 1:length(data)
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    ff{curr_animal} = nan(n_rois,14);
    for curr_session = 1:length(data(curr_animal).im);        
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
    
    ff_mq{curr_animal} = cat(3,curr_ff_m,curr_ff_q);
end

% Total activity
total_act = cell(length(data),1);
for curr_animal = 1:length(data);
        
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_act = nan(n_rois,14,3);
    
    for i = 1:length(data(curr_animal).im)
        curr_act(:,i,1) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,2) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,~analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,3) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,:) > 0,2);
    end
       
    total_act{curr_animal} = curr_act;
   
end

ff_cat = cell2mat(cellfun(@(x) x(:),ff,'uni',false));
ff_mq_cat = cell2mat(cellfun(@(x) reshape(x,[],2),ff_mq,'uni',false));
total_act_cat = cell2mat(cellfun(@(x) reshape(x,[],3,1),total_act,'uni',false));

% Plot density plot of early v. late differences
n_bins = 100;

bin1_edges = linspace(min(ff_cat(~isinf(ff_cat),1)), ...
    max(ff_cat,[],1),n_bins);
bin1_centers = bin1_edges(1:end-1) + diff(bin1_edges)/2;
bin1_label = round(bin1_centers*100)/100;
bin1_edges(1) = -Inf;
%bin1_edges(end) = Inf;

bin2_edges = linspace(0,0.2,n_bins+1);
bin2_centers = bin2_edges(1:end-1) + diff(bin2_edges)/2;
bin2_label = round(bin2_centers*100)/100;
bin2_edges(1) = -Inf;
bin2_edges(end) = Inf;

bin_edges = {bin1_edges,bin2_edges};

figure;

axis square;
reliability_density = hist3([ff_cat,total_act_cat(:,3)],'Edges',bin_edges);
h = fspecial('gaussian',20,2);
reliability_density_smooth = imfilter(reliability_density,h,'same');
imagesc(reliability_density_smooth);colormap(hot);
set(gca,'YDir','normal')
xlabel('Fraction active during session')
ylabel('Fraction active during fano factor')
xlim([0.5 n_bins-1]);
ylim([0.5 n_bins-1]);
set(gca,'XTickLabel',bin2_label(get(gca,'XTick')));
set(gca,'YTickLabel',bin1_label(get(gca,'YTick')));
c = colorbar;
ylabel(c,'Number of cells');


%% Binary population correlation

% Correlation
m_pop_corr = nan(14,14,length(data));
q_pop_corr = nan(14,14,length(data));

for curr_animal = 1:length(data)
    
    curr_m_pop_corr = corrcoef(classified_rois(curr_animal).movement);
    curr_q_pop_corr = corrcoef(classified_rois(curr_animal).quiescent);
    
    m_pop_corr(1:length(curr_m_pop_corr), ...
        1:length(curr_m_pop_corr),curr_animal) = curr_m_pop_corr;
    q_pop_corr(1:length(curr_q_pop_corr),...
        1:length(curr_q_pop_corr),curr_animal) = curr_q_pop_corr;
    
end

% Shuffled p-value
n_shuff = 10000;
m_pop_corr_shuff = nan(14,14,n_shuff);
q_pop_corr_shuff = nan(14,14,n_shuff);
for i = 1:n_shuff
    m_pop_corr_shuff(:,:,i) = nanmean(shake(m_pop_corr,'diag'),3);
    q_pop_corr_shuff(:,:,i) = nanmean(shake(q_pop_corr,'diag'),3);
end

m_pop_corr_shuffdiff = nanmean(m_pop_corr,3) - nanmean(m_pop_corr_shuff,3);
m_pop_corr_rank = permute(tiedrank(permute(cat(3,nanmean(m_pop_corr,3),m_pop_corr_shuff),[3,2,1])),[2,3,1]);
m_pop_corr_p = m_pop_corr_rank(:,:,1)./(n_shuff+1);
   
q_pop_corr_shuffdiff = nanmean(q_pop_corr,3) - nanmean(q_pop_corr_shuff,3);
q_pop_corr_rank = permute(tiedrank(permute(cat(3,nanmean(q_pop_corr,3),q_pop_corr_shuff),[3,2,1])),[2,3,1]);
q_pop_corr_p = q_pop_corr_rank(:,:,1)./(n_shuff+1);

figure;
subplot(1,2,1);
imagesc(nanmean(m_pop_corr,3));
caxis([0 1]);
title('M pop corr')

subplot(1,2,2);
imagesc(nanmean(q_pop_corr,3));
caxis([0 1]);
title('Q pop corr');
colormap(hot);

figure;
subplot(2,2,1);
imagesc(m_pop_corr_p);
caxis([0 1]);
title('M pop corr p')

subplot(2,2,2);
imagesc(q_pop_corr_p);
caxis([0 1]);
title('Q pop corr p');

subplot(2,2,3);
imagesc(m_pop_corr_shuffdiff);
title('M pop corr shuffdiff')

subplot(2,2,4);
imagesc(q_pop_corr_shuffdiff);
title('Q pop corr shuffdiff');

% Make blue-to-red colormap
curr_colormap = colormap;
c_vals = size(curr_colormap,1);
bwr_colormap = [repmat(linspace(0,1,ceil(c_vals/2))',1,2),ones(ceil(c_vals/2),1); ...
    ones(ceil(c_vals/2),1),repmat(linspace(1,0,ceil(c_vals/2))',1,2)];
colormap(bwr_colormap);





% Total activity
total_act = cell(length(data),1);
for curr_animal = 1:length(data);
        
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    curr_act = nan(n_rois,14,3);
    
    for i = 1:length(data(curr_animal).im)
        curr_act(:,i,1) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,2) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,~analysis(curr_animal).lever(i).lever_move_frames) > 0,2);
        curr_act(:,i,3) = nanmean(data(curr_animal).im(i). ...
            roi_trace_thresh(:,:) > 0,2);
    end
       
    total_act{curr_animal} = curr_act;
   
end


total_pop_corr = nan(14,14,length(data));
total_pop_corr_m = nan(14,14,length(data));
total_pop_corr_q = nan(14,14,length(data));

for curr_animal = 1:length(data)
    
    n_days = size(total_act{curr_animal},2);
    
    curr_total_pop_corr = corrcoef(total_act{curr_animal}(:,:,3));
    curr_total_pop_corr_m = corrcoef(total_act{curr_animal}(:,:,1));
    curr_total_pop_corr_q = corrcoef(total_act{curr_animal}(:,:,2));
    
    total_pop_corr(1:n_days, ...
        1:n_days,curr_animal) = curr_total_pop_corr;
    total_pop_corr_m(1:n_days, ...
        1:n_days,curr_animal) = curr_total_pop_corr_m;
    total_pop_corr_q(1:n_days, ...
        1:n_days,curr_animal) = curr_total_pop_corr_q;
    
end


%% Preferred timing of activity relative to movement onset/offset

sig_timing = struct('movement',cell(length(data),1),'quiescence',cell(length(data),1));

for curr_animal = 1:length(data)
    
    disp(['Animal ' num2str(curr_animal)]);
    
    onset_frames = analysis(curr_animal).surrounding_frames >= 0 & ...
        analysis(curr_animal).surrounding_frames <= 45;
    offset_frames = analysis(curr_animal).surrounding_frames >= -45 & ...
        analysis(curr_animal).surrounding_frames <= 0;
    
    total_frames = sum([onset_frames,offset_frames]);
    
    sig_timing(curr_animal).movement = nan(size(data(curr_animal).im(1).roi_trace_df,1), ...
        total_frames,length(data(curr_animal).im));
    
    sig_timing(curr_animal).quiescence = nan(size(data(curr_animal).im(1).roi_trace_df,1), ...
        total_frames,length(data(curr_animal).im));
    
    for curr_day = 1:length(analysis(curr_animal).im);
        
        move_onset_activity = analysis(curr_animal).im(curr_day).move_onset_aligned;
        move_offset_activity = analysis(curr_animal).im(curr_day).move_offset_aligned;
        
        % Move and quiescent activity: concatenate 1.5s after/before onset
        move_activity = [move_onset_activity(:,onset_frames,:), ...
            move_offset_activity(:,offset_frames,:)];
        
        quiescent_activity = [move_offset_activity(:,onset_frames,:), ...
            move_onset_activity(:,offset_frames,:)];
        
        num_shuff = 1000;
        move_activity_shuff_mean = nan(num_shuff,size(move_activity,2),size(move_activity,3));
        quiescent_activity_shuff_mean = nan(num_shuff,size(quiescent_activity,2),size(quiescent_activity,3));
        for curr_shuff = 1:num_shuff
            
            curr_circshift = randi(size(move_activity,2),size(move_activity,1),1);
            curr_circshift_idx = cell2mat(arrayfun(@(x) ...
                [curr_circshift(x):total_frames,1:curr_circshift(x)-1], ...
                transpose(1:length(curr_circshift)),'uni',false));
            
            curr_move_activity_shuff = move_activity;
            curr_quiescent_activity_shuff = quiescent_activity;
            
            for curr_trial = 1:size(move_activity,1)
                
                % circshift function is slower than indexing, don't use
%                 curr_move_activity_shuff(curr_trial,:,:) = circshift( ...
%                     curr_move_activity_shuff(curr_trial,:,:),[0,curr_circshift(curr_trial)]);
%                 
%                 curr_quiescent_activity_shuff(curr_trial,:,:) = circshift( ...
%                     curr_quiescent_activity_shuff(curr_trial,:,:),[0,curr_circshift(curr_trial)]);
                
                curr_move_activity_shuff(curr_trial,:,:) =  ...
                    curr_move_activity_shuff(curr_trial,curr_circshift_idx(curr_trial,:),:);
                
                curr_quiescent_activity_shuff(curr_trial,:,:) =  ...
                    curr_quiescent_activity_shuff(curr_trial,curr_circshift_idx(curr_trial,:),:);
                
            end
            
            move_activity_shuff_mean(curr_shuff,:,:) = nanmean(curr_move_activity_shuff,1);
            quiescent_activity_shuff_mean(curr_shuff,:,:) = nanmean(curr_quiescent_activity_shuff,1);
            
        end
        
        move_activity_mean = permute(nanmean(move_activity,1),[3,2,1]);
        quiescent_activity_mean = permute(nanmean(quiescent_activity,1),[3,2,1]);
        
        move_activity_shuff_prctile = permute(prctile(move_activity_shuff_mean,95,1),[3,2,1]);
        quiescent_activity_shuff_prctile = permute(prctile(quiescent_activity_shuff_mean,95,1),[3,2,1]);
        
        move_activity_sig_timing = move_activity_mean > move_activity_shuff_prctile;
        quiescent_activity_sig_timing = quiescent_activity_mean > quiescent_activity_shuff_prctile;
        
        sig_timing(curr_animal).movement(:,:,curr_day) = move_activity_sig_timing;
        sig_timing(curr_animal).quiescence(:,:,curr_day) = quiescent_activity_sig_timing;
        
        disp(['Day ' num2str(curr_day)]);
    end
    
end

% Fraction of movement/quiescent ROIs with preferred timing
sig_movement_frac_m = nan(length(data),14);
sig_movement_frac_q = nan(length(data),14);
sig_quiescence_frac_m = nan(length(data),14);
sig_quiescence_frac_q = nan(length(data),14);

for curr_animal = 1:length(data)
    
    any_sig_movement = permute(any(sig_timing(curr_animal).movement,2),[1,3,2]);
    any_sig_quiescent = permute(any(sig_timing(curr_animal).quiescence,2),[1,3,2]);
    
    sig_movement_frac_m(curr_animal,1:size(any_sig_movement,2)) = ...
        sum(any_sig_movement & classified_rois(curr_animal).movement)./ ...
        sum(classified_rois(curr_animal).movement);
    
    sig_movement_frac_q(curr_animal,1:size(any_sig_movement,2)) = ...
        sum(any_sig_movement & classified_rois(curr_animal).quiescent)./ ...
        sum(classified_rois(curr_animal).quiescent);
    
    sig_quiescence_frac_m(curr_animal,1:size(any_sig_quiescent,2)) = ...
        sum(any_sig_quiescent & classified_rois(curr_animal).movement)./ ...
        sum(classified_rois(curr_animal).movement);
    
    sig_quiescence_frac_q(curr_animal,1:size(any_sig_quiescent,2)) = ...
        sum(any_sig_quiescent & classified_rois(curr_animal).quiescent)./ ...
        sum(classified_rois(curr_animal).quiescent);
    
end

figure;
subplot(1,2,1); hold on;
errorbar(nanmean(sig_movement_frac_m),nanstd(sig_movement_frac_m)./ ...
    sqrt(sum(~isnan(sig_movement_frac_m))),'r','linewidth',2);
errorbar(nanmean(sig_movement_frac_q),nanstd(sig_movement_frac_q)./ ...
    sqrt(sum(~isnan(sig_movement_frac_q))),'k','linewidth',2);
title('Movement');
ylabel('Fraction classified ROIs sig timing pref');
legend({'M','Q'});

subplot(1,2,2); hold on;
errorbar(nanmean(sig_quiescence_frac_m),nanstd(sig_quiescence_frac_m)./ ...
    sqrt(sum(~isnan(sig_quiescence_frac_m))),'r','linewidth',2);
errorbar(nanmean(sig_quiescence_frac_q),nanstd(sig_quiescence_frac_q)./ ...
    sqrt(sum(~isnan(sig_quiescence_frac_q))),'k','linewidth',2);
title('Quiescence');
ylabel('Fraction classified ROIs sig timing pref');
legend({'M','Q'});


% Fraction of classified ROIs with timing preferred at each frame
sig_timing_movement_m_avg = ...
    cell2mat(permute(cellfun(@(timing,class) cell2mat(arrayfun(@(x) nanmean(timing( ...
    class(:,x),:,x),1),transpose(1:size(class,2)),'uni',false)), ...
    {sig_timing.movement},{classified_rois.movement},'uni',false),[1,3,2]));

sig_timing_movement_q_avg = ...
    cell2mat(permute(cellfun(@(timing,class) cell2mat(arrayfun(@(x) nanmean(timing( ...
    class(:,x),:,x),1),transpose(1:size(class,2)),'uni',false)), ...
    {sig_timing.movement},{classified_rois.quiescent},'uni',false),[1,3,2]));

sig_timing_quiescence_m_avg = ...
    cell2mat(permute(cellfun(@(timing,class) cell2mat(arrayfun(@(x) nanmean(timing( ...
    class(:,x),:,x),1),transpose(1:size(class,2)),'uni',false)), ...
    {sig_timing.quiescence},{classified_rois.movement},'uni',false),[1,3,2]));

sig_timing_quiescence_q_avg = ...
    cell2mat(permute(cellfun(@(timing,class) cell2mat(arrayfun(@(x) nanmean(timing( ...
    class(:,x),:,x),1),transpose(1:size(class,2)),'uni',false)), ...
    {sig_timing.quiescence},{classified_rois.quiescent},'uni',false),[1,3,2]));

figure;
subplot(2,2,1); hold on;
set(gca,'ColorOrder',jet(14));
plot(nanmean(sig_timing_movement_m_avg,3)');
title('Movement, M ROIs');

subplot(2,2,2); hold on;
set(gca,'ColorOrder',jet(14));
plot(nanmean(sig_timing_quiescence_m_avg,3)');
title('Quiescence, M ROIs');

subplot(2,2,3); hold on;
set(gca,'ColorOrder',jet(14));
plot(nanmean(sig_timing_movement_q_avg,3)');
title('Movement, Q ROIs');

subplot(2,2,4); hold on;
set(gca,'ColorOrder',jet(14));
plot(nanmean(sig_timing_quiescence_q_avg,3)');
title('Quiescence, Q ROIs');



% Plot preferred timing of all classified cells during movement/quiescence

sig_timing_movement_m = ...
    cellfun(@(timing,class) arrayfun(@(x) timing( ...
    class(:,x),:,x),transpose(1:size(class,2)),'uni',false), ...
    {sig_timing.movement},{classified_rois.movement},'uni',false);

sig_timing_movement_q = ...
    cellfun(@(timing,class) arrayfun(@(x) timing( ...
    class(:,x),:,x),transpose(1:size(class,2)),'uni',false), ...
    {sig_timing.movement},{classified_rois.quiescent},'uni',false);

sig_timing_quiescence_m = ...
    cellfun(@(timing,class) arrayfun(@(x) timing( ...
    class(:,x),:,x),transpose(1:size(class,2)),'uni',false), ...
    {sig_timing.quiescence},{classified_rois.movement},'uni',false);

sig_timing_quiescence_q = ...
    cellfun(@(timing,class) arrayfun(@(x) timing( ...
    class(:,x),:,x),transpose(1:size(class,2)),'uni',false), ...
    {sig_timing.quiescence},{classified_rois.quiescent},'uni',false);


figure;
for i = 1:4
    switch i
        case 1
            cat_data = vertcat(sig_timing_movement_m{:});
            cat_data = vertcat(cat_data{:});
            time_data = cat_data.*repmat(1:size(cat_data,2),size(cat_data,1),1);
            mean_time_data = sum(time_data,2)./sum(cat_data,2);
            [~,sort_idx] = sort(mean_time_data);
            subplot(1,4,i);
            imagesc(cat_data(sort_idx,:));colormap(gray)
            ylabel('Sorted M ROIs');
            title('Movement, M ROIs')
            line(repmat(ceil(size(cat_data,2)/2),1,2),ylim,'color','r');
        case 2
            cat_data = vertcat(sig_timing_movement_q{:});
            cat_data = vertcat(cat_data{:});
            time_data = cat_data.*repmat(1:size(cat_data,2),size(cat_data,1),1);
            mean_time_data = sum(time_data,2)./sum(cat_data,2);
            [~,sort_idx] = sort(mean_time_data);
            subplot(1,4,i);
            imagesc(cat_data(sort_idx,:));colormap(gray)
            ylabel('Sorted Q ROIs');
            title('Movement, Q ROIs')
            line(repmat(ceil(size(cat_data,2)/2),1,2),ylim,'color','r');
        case 3
            cat_data = vertcat(sig_timing_quiescence_m{:});
            cat_data = vertcat(cat_data{:});
            time_data = cat_data.*repmat(1:size(cat_data,2),size(cat_data,1),1);
            mean_time_data = sum(time_data,2)./sum(cat_data,2);
            [~,sort_idx] = sort(mean_time_data);
            subplot(1,4,i);
            imagesc(cat_data(sort_idx,:));colormap(gray)
            ylabel('Sorted M ROIs');
            title('Quiescence, M ROIs')
            line(repmat(ceil(size(cat_data,2)/2),1,2),ylim,'color','r');
        case 4
            cat_data = vertcat(sig_timing_quiescence_q{:});
            cat_data = vertcat(cat_data{:});
            time_data = cat_data.*repmat(1:size(cat_data,2),size(cat_data,1),1);
            mean_time_data = sum(time_data,2)./sum(cat_data,2);
            [~,sort_idx] = sort(mean_time_data);
            subplot(1,4,i);
            imagesc(cat_data(sort_idx,:));colormap(gray)
            ylabel('Sorted Q ROIs');
            title('Quiescence, Q ROIs')
            line(repmat(ceil(size(cat_data,2)/2),1,2),ylim,'color','r');            
    end   
end


%% Histogram of activity distribution relative to movement
% Meant to be like the 2nd cell, but simplified

% Get fraction of activity during session and during movement
activity_move_frac = cell(length(data),1);
activity_move_frac_shuff = cell(length(data),1);

for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    activity_frac{curr_animal} = nan(n_rois,n_days);
    activity_move_frac{curr_animal} = nan(n_rois,n_days);
    
    activity_frac_shuff{curr_animal} = nan(n_rois,n_days);
    activity_move_frac_shuff{curr_animal} = nan(n_rois,n_days);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        
        warning off
        activity_move_frac{curr_animal}(:,curr_day) = ...
            (curr_act* ...
            analysis(curr_animal).lever(curr_day).lever_move_frames) ./ ...
            sum(curr_act,2);
        warning on
        
        % Activity-epoch shuffled data
        
        % Used to shuffle activity, but might as well do it like classif.
        %         boundary_frames = arrayfun(@(x) find(diff([Inf,curr_act(x,:),Inf]) ~= 0),1:size(curr_act,1),'uni',false);
        %         curr_act_shuff = cell2mat(arrayfun(@(x) cell2mat(shake( ...
        %             mat2cell(curr_act(x,:),1,diff(boundary_frames{x})))),[1:size(curr_act,1)]','uni',false));
        
        % Shuffle movement/quiescent epochs
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
        num_rep = 1;
        shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
        lever_active_shuffle = nan(length(lever_active_frames),num_rep);
        for i = 1:num_rep
            lever_active_shuffle(:,i) = ...
                vertcat(lever_active_frames_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        warning off
        %         activity_move_frac_shuff{curr_animal}(:,curr_day) = ...
        %             (curr_act_shuff* ...
        %             analysis(curr_animal).lever(curr_day).lever_move_frames) ./ ...
        %             sum(curr_act_shuff,2);
        activity_move_frac_shuff{curr_animal}(:,curr_day) = ...
            (curr_act* ...
            lever_active_shuffle) ./ ...
            sum(curr_act,2);
        warning on
        
    end
    disp(curr_animal);
end

% Histogram of real vs. shuffled (all days/ROIs concatenated)
activity_move_frac_cat = cell2mat((cellfun(@(x) x(:),activity_move_frac,'uni',false)));
activity_move_frac_shuff_cat = cell2mat((cellfun(@(x) x(:),activity_move_frac_shuff,'uni',false)));

% Sanity check, have to be more fancy if not full range
if min(activity_move_frac_cat) ~= 0 || ...
        min(activity_move_frac_shuff_cat) ~= 1 || ...
        max(activity_move_frac_cat) ~= 0 || ...
        max(activity_move_frac_shuff_cat) ~= 1 || ...
        error('Not full range')
end

x_plot = linspace(0,1,100);
real_dist = hist(activity_move_frac_cat,100);
shuff_dist = hist(activity_move_frac_shuff_cat,100);

figure; hold on;
plot(x_plot,real_dist,'r');
plot(x_plot,shuff_dist,'k');
legend({'Real','Shuffled'})

ylabel('Number of ROIs');
xlabel('Fraction of activity during movement');

%% Histogram of activity distribution relative to movement/quiescent rate
% Meant to be like the 2nd cell, but simplified

% Get fraction of activity during session and during movement
activity_mq_rate_ratio = cell(length(data),1);
activity_mq_rate_ratio_shuff = cell(length(data),1);
activity_mq_rate_ratio_class = cell(length(data),1);


for curr_animal = 1:length(data);
    
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    activity_mq_rate_ratio{curr_animal} = cell(n_days,1);
    activity_mq_rate_ratio_shuff{curr_animal} = cell(n_days,1);
    activity_mq_rate_ratio_class{curr_animal} = cell(n_days,1);

    
    for curr_day = 1:length(data(curr_animal).im)
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        % Only use "active" ROIs, at least 5 events
        min_events = 5;
        active_rois = sum(diff(curr_act > 0,[],2) == 1,2) > min_events;
        curr_act = curr_act(active_rois,:);
        
        % "Active" classified ROIs
        curr_class = classified_rois(curr_animal).movement(active_rois,curr_day) | ...
            classified_rois(curr_animal).quiescent(active_rois,curr_day);
        
        curr_act_m = (curr_act*lever_active_frames)./sum(lever_active_frames);
        curr_act_q = (curr_act*~lever_active_frames)./sum(~lever_active_frames);
        curr_act_mq = (curr_act_m - curr_act_q)./(curr_act_m + curr_act_q);
        activity_mq_rate_ratio{curr_animal}{curr_day} = curr_act_mq;
        
        activity_mq_rate_ratio_class{curr_animal}{curr_day} = curr_act_mq(curr_class);
        
        % Activity-epoch shuffled data
        
        % Used to shuffle activity, but might as well do it like classif.
        %         boundary_frames = arrayfun(@(x) find(diff([Inf,curr_act(x,:),Inf]) ~= 0),1:size(curr_act,1),'uni',false);
        %         curr_act_shuff = cell2mat(arrayfun(@(x) cell2mat(shake( ...
        %             mat2cell(curr_act(x,:),1,diff(boundary_frames{x})))),[1:size(curr_act,1)]','uni',false));
        
        % Shuffle movement/quiescent epochs
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
        num_rep = 1000;
        shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
        lever_active_shuffle = nan(length(lever_active_frames),num_rep);
        for i = 1:num_rep
            lever_active_shuffle(:,i) = ...
                vertcat(lever_active_frames_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        curr_act_m_shuff = (curr_act*lever_active_shuffle)./sum(lever_active_frames);
        curr_act_q_shuff = (curr_act*~lever_active_shuffle)./sum(~lever_active_frames);
        curr_act_mq_shuff = (curr_act_m_shuff - curr_act_q_shuff)./(curr_act_m_shuff + curr_act_q_shuff);
        activity_mq_rate_ratio_shuff{curr_animal}{curr_day} = curr_act_mq_shuff(:);
     
    end
    disp(curr_animal);
end

% Histogram of real vs. shuffled (all days/ROIs concatenated)
activity_mq_ratio_cat = cell2mat((cellfun(@(x) vertcat(x{:}),activity_mq_rate_ratio,'uni',false)));
activity_mq_ratio_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),activity_mq_rate_ratio_shuff,'uni',false)));
activity_mq_ratio_class_cat = cell2mat((cellfun(@(x) vertcat(x{:}),activity_mq_rate_ratio_class,'uni',false)));

bin_edges = linspace(-1,1,100);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(activity_mq_ratio_cat,bin_centers)./length(activity_mq_ratio_cat);
shuff_dist = hist(activity_mq_ratio_shuff_cat,bin_centers)./length(activity_mq_ratio_shuff_cat);
class_dist = hist(activity_mq_ratio_class_cat,bin_centers)./length(activity_mq_ratio_cat);

figure; hold on;
plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
plot(bin_centers,class_dist,'r','linewidth',2)
legend({'Shuffled','Real','Classified'})

ylabel('Number of ROIs');
xlabel('Rate M-Q/M+Q');



%% Quiescent activity vs. quiescent epoch lever position

q_epoch_act = cell(length(data),1);
q_epoch_lever = cell(length(data),1);

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get quiescent activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_act_q = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        
        % Get lever position/velocity by frame   
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = analysis(curr_animal).lever(curr_day).lever_position_frames;
        lever_velocity_frames = analysis(curr_animal).lever(curr_day).lever_velocity_frames;
        
        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
        lever_split = mat2cell(lever_position_frames,diff(boundary_times));
        
        frames_split = mat2cell(transpose(1:length(lever_active_frames)),diff(boundary_times));
        
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
        
        % Used to look at adjacent movement epochs
        lever_split_q_idx = find(lever_split_q);
        lever_split_q_idx_use = lever_split_q_idx(lever_split_q_idx > 1 & ...
            lever_split_q_idx < length(lever_split_q));
        lever_split_pre_m_idx = lever_split_q_idx_use - 1;
        lever_split_post_m_idx = lever_split_q_idx_use - 1;
        
        split_q_act = cellfun(@(frames) curr_act_q(:,frames), ...
            frames_split(lever_split_q_idx_use),'uni',false);
        split_q_act_mean = cellfun(@(x) nanmean(x(:)),split_q_act);
        
        % Get total median and quiescent epoch median of lever
        median_lever = median(lever_position_frames);
        mean_lever_q = cellfun(@mean,lever_split(lever_split_q_idx_use));
        
        % Store quiescent epoch lever and activity
        q_epoch_act{curr_animal}{curr_day} = split_q_act_mean;
        q_epoch_lever{curr_animal}{curr_day} = mean_lever_q-median_lever;
        
    end
    disp(curr_animal);
end

% Bin lever and get groupstats of activity
q_epoch_lever_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),q_epoch_lever,'uni',false);
q_epoch_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),q_epoch_act,'uni',false);

bin_edges = linspace(-3,3,10);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_edges(end) = Inf;
bin_edges(1) = -Inf;
[~,q_epoch_lever_bin] = cellfun(@(x) cellfun(@(x) histc(x,bin_edges), ...
    x,'uni',false),q_epoch_lever_norm,'uni',false);
q_epoch_act_binmean = cellfun(@(act,bins) cellfun(@(act,bins) ...
    grpstats(act,bins,'nanmean'),act,bins,'uni',false), ...
    q_epoch_act_norm,q_epoch_lever_bin,'uni',false);


q_epoch_act_binmean_cat = cell(length(data),14);
q_epoch_lever_bin_cat = cell(length(data),14);
for curr_animal = 1:length(data)
    q_epoch_act_binmean_cat(curr_animal,1:length(q_epoch_act_binmean{curr_animal})) = ...
        q_epoch_act_binmean{curr_animal};
    
    q_epoch_lever_bin_cat(curr_animal,1:length(q_epoch_lever_bin{curr_animal})) = ...
        q_epoch_lever_bin{curr_animal};
end
use_bins = cellfun(@unique,q_epoch_lever_bin_cat,'uni',false);

q_epoch_act_binmean_pad = cellfun(@(x) nan(length(bin_edges)-1,1), ...
    q_epoch_act_binmean_cat,'uni',false);
for i = 1:numel(q_epoch_act_binmean_cat)
    q_epoch_act_binmean_pad{i}(use_bins{i}) = ...
        q_epoch_act_binmean_cat{i};
end

% Average sessions within animals
q_epoch_act_binmean_animal = nan(length(bin_edges)-1,length(data));
for curr_animal = 1:length(data);
    q_epoch_act_binmean_animal(:,curr_animal) = ....
        nanmean(horzcat(q_epoch_act_binmean_pad{curr_animal,:}),2);
end

% Average animals within sessions
q_epoch_act_binmean_day = nan(length(bin_edges)-1,14);
for curr_day = 1:14;
    q_epoch_act_binmean_day(:,curr_day) = ....
        nanmean(horzcat(q_epoch_act_binmean_pad{:,curr_day}),2);
end

% Plot all animals
figure;
sq = ceil(sqrt(length(data)));
for curr_animal = 1:length(data)
    subplot(sq,sq,curr_animal);hold on;
    set(gca,'ColorOrder',jet(14));
    curr_data = horzcat(q_epoch_act_binmean_pad{curr_animal,:});
    plot(curr_data);
    title(['Animal ' num2str(curr_animal)]);
end

% Plot session / animal average
figure; 
subplot(1,2,1);hold on
plot(bin_centers,q_epoch_act_binmean_animal,'color',[0.5,0.5,0.5]);
plot(bin_centers,nanmean(q_epoch_act_binmean_animal,2),'color','k','linewidth',2);
title('Average across sessions');
xlabel('Normalized lever position');
ylabel('Normalized quiescent activity');

subplot(1,2,2);hold on
set(gca,'ColorOrder',jet(size(q_epoch_act_binmean_day,2)));
plot(bin_centers,q_epoch_act_binmean_day);
plot(bin_centers,nanmean(q_epoch_act_binmean_day,2),'color','k','linewidth',2);
title('Average across animals');
xlabel('Normalized lever position')
ylabel('Normalized quiescent activity');

% Plot average across animals/sessions
q_epoch_act_binmean_pad_cat = horzcat(q_epoch_act_binmean_pad{:});
figure; hold on;
errorbar(bin_centers,nanmean(q_epoch_act_binmean_pad_cat,2), ...
    nanstd(q_epoch_act_binmean_pad_cat,[],2)./ ...
    sqrt(sum(~isnan(q_epoch_act_binmean_pad_cat),2)),'k','linewidth',2);
title('Average across all animals/sessions');
xlabel('Normalized lever position')
ylabel('Normalized quiescent activity');

%% Movement activity vs. movement epoch lever position

m_epoch_act = cell(length(data),1);
m_epoch_lever = cell(length(data),1);

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get quiescent activity
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act_m = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_m,:) > 0);
        
        % Get lever position/velocity by frame   
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = analysis(curr_animal).lever(curr_day).lever_position_frames;
        lever_velocity_frames = analysis(curr_animal).lever(curr_day).lever_velocity_frames;
        
        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
        lever_split = mat2cell(lever_velocity_frames,diff(boundary_times));
        
        frames_split = mat2cell(transpose(1:length(lever_active_frames)),diff(boundary_times));
        
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
        
        % Used to look at adjacent movement epochs
        lever_split_m_idx = find(lever_split_m);
        lever_split_m_idx_use = lever_split_m_idx(lever_split_m_idx > 1 & ...
            lever_split_m_idx < length(lever_split_m));
        lever_split_pre_m_idx = lever_split_m_idx_use - 1;
        lever_split_post_m_idx = lever_split_m_idx_use - 1;
        
        split_m_act = cellfun(@(frames) curr_act_m(:,frames), ...
            frames_split(lever_split_m_idx_use),'uni',false);
        split_m_act_mean = cellfun(@(x) nanmean(x(:)),split_m_act);
        
        % Get total median and quiescent epoch median of lever
        median_lever = median(lever_position_frames);
        mean_lever_m = cellfun(@mean,lever_split(lever_split_m_idx_use));
        
        % Store quiescent epoch lever and activity
        m_epoch_act{curr_animal}{curr_day} = split_m_act_mean;
        m_epoch_lever{curr_animal}{curr_day} = mean_lever_m;
        
    end
    disp(curr_animal);
end

% Bin lever and get groupstats of activity
m_epoch_lever_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),m_epoch_lever,'uni',false);
m_epoch_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),m_epoch_act,'uni',false);

bin_edges = linspace(0,3,10);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_edges(end) = Inf;
bin_edges(1) = -Inf;
[~,m_epoch_lever_bin] = cellfun(@(x) cellfun(@(x) histc(x,bin_edges), ...
    x,'uni',false),m_epoch_lever_norm,'uni',false);
m_epoch_act_binmean = cellfun(@(act,bins) cellfun(@(act,bins) ...
    grpstats(act,bins,'nanmean'),act,bins,'uni',false), ...
    m_epoch_act_norm,m_epoch_lever_bin,'uni',false);


m_epoch_act_binmean_cat = cell(length(data),14);
m_epoch_lever_bin_cat = cell(length(data),14);
for curr_animal = 1:length(data)
    m_epoch_act_binmean_cat(curr_animal,1:length(m_epoch_act_binmean{curr_animal})) = ...
        m_epoch_act_binmean{curr_animal};
    
    m_epoch_lever_bin_cat(curr_animal,1:length(m_epoch_lever_bin{curr_animal})) = ...
        m_epoch_lever_bin{curr_animal};
end
use_bins = cellfun(@unique,m_epoch_lever_bin_cat,'uni',false);

m_epoch_act_binmean_pad = cellfun(@(x) nan(length(bin_edges)-1,1), ...
    m_epoch_act_binmean_cat,'uni',false);
for i = 1:numel(m_epoch_act_binmean_cat)
    m_epoch_act_binmean_pad{i}(use_bins{i}) = ...
        m_epoch_act_binmean_cat{i};
end


% Plot average across animals/sessions
m_epoch_act_binmean_pad_cat = horzcat(m_epoch_act_binmean_pad{:});
figure; hold on;
errorbar(bin_centers,nanmean(m_epoch_act_binmean_pad_cat,2), ...
    nanstd(m_epoch_act_binmean_pad_cat,[],2)./ ...
    sqrt(sum(~isnan(m_epoch_act_binmean_pad_cat),2)),'k','linewidth',2);
title('Average across all animals/sessions');
xlabel('Normalized lever velocity')
ylabel('Normalized movement activity');


%% Quiescent activity vs. running smoothed lever position
% Can modify this to look at whole day or q/m epochs
% qROI activity duirng smoothed trace of entire day: this is an inverted U,
% the more away from median the less quiescent activity (which is expected,
% since far away from median when moving around)

q_act = cell(length(data),1);
smoothed_lever = cell(length(data),1);

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get quiescent activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_act_q = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        
        % Get lever force        
        [lever_active,lever_force_resample,~,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated( ...
            data(curr_animal).bhv(curr_day).lever_force);
        
        frame_times = (data(curr_animal).bhv(curr_day).frame_times( ...
            1:size(curr_act_q,2)))*1000;
        
        imaged_lever_idx = ((1:length(lever_active)) >= ...
            frame_times(1)) & ...
            ((1:length(lever_active)) <= ...
            ceil(frame_times(end)+median(diff(frame_times))));
        
        lever_active_imaged = lever_active(imaged_lever_idx);
        lever_force_resample_imaged = lever_force_resample(imaged_lever_idx);
        smooth_time = 1;
        lever_force_resample_imaged_smoothed = ...
            smooth(lever_force_resample_imaged,smooth_time);
        
        relative_frame_times = ceil(frame_times - frame_times(1));
        lever_time_frames = nan(1,length(lever_force_resample_imaged));
        for curr_frame = 1:length(relative_frame_times);
            if curr_frame ~= length(relative_frame_times)
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    relative_frame_times(curr_frame+1)) = curr_frame;
            else
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    end) = curr_frame;
            end
        end
              
        % Split lever by frames
        frame_boundaries = [0,find(lever_time_frames(2:end) ~= ...
            lever_time_frames(1:end-1)),length(lever_time_frames)];
        lever_frame_split = mat2cell(lever_force_resample_imaged_smoothed, ...
            diff(frame_boundaries)',1);
        
        lever_active_split = mat2cell(lever_active_imaged, ...
            diff(frame_boundaries)',1);
        
        % Get average lever force in each frame
        lever_frame_split_mean = cellfun(@mean,lever_frame_split);
        
        % Get whether movement occurs in each frame
        lever_active_frames = cellfun(@any,lever_active_split);
        
        % Split lever/activity by movement/quiescence
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_times));
        lever_split = mat2cell(lever_frame_split_mean,diff(boundary_times));
        
        mean_act_q_split = mat2cell(nanmean(curr_act_q,1)',diff(boundary_times));
        
        split_q = cellfun(@(x) ~any(x),lever_active_frames_split);
        
        % Only use lever/activity during quiescence
        lever_q = vertcat(lever_split{split_q});
        mean_q_act_q = vertcat(mean_act_q_split{split_q});
               
        % Get median lever of entire day
        median_lever = median(lever_q);
        
        % Store smoothed lever and activity
        q_act{curr_animal}{curr_day} = mean_q_act_q;
        smoothed_lever{curr_animal}{curr_day} = lever_q-median_lever;
        
    end
    disp(curr_animal);
end

% Bin lever and get groupstats of activity
smoothed_lever_norm = cellfun(@(x) cellfun(@(x) x./max(abs(x)), ...
    x,'uni',false),smoothed_lever,'uni',false);
q_act_norm = cellfun(@(x) cellfun(@(x) (x-min(x))./max(x-min(x)), ...
    x,'uni',false),q_act,'uni',false);

bin_edges = linspace(-1,1,10);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_edges(end) = Inf;
[~,smoothed_lever_bin] = cellfun(@(x) cellfun(@(x) histc(x,bin_edges), ...
    x,'uni',false),smoothed_lever_norm,'uni',false);
q_act_binmean = cellfun(@(act,bins) cellfun(@(act,bins) ...
    grpstats(act,bins,'nanmean'),act,bins,'uni',false), ...
    q_act_norm,smoothed_lever_bin,'uni',false);


q_act_binmean_cat = cell(length(data),14);
smoothed_lever_bin_cat = cell(length(data),14);
for curr_animal = 1:length(data)
    q_act_binmean_cat(curr_animal,1:length(q_act_binmean{curr_animal})) = ...
        q_act_binmean{curr_animal};
    
    smoothed_lever_bin_cat(curr_animal,1:length(smoothed_lever_bin{curr_animal})) = ...
        smoothed_lever_bin{curr_animal};
end
use_bins = cellfun(@unique,smoothed_lever_bin_cat,'uni',false);

q_act_binmean_pad = cellfun(@(x) nan(length(bin_edges)-1,1), ...
    q_act_binmean_cat,'uni',false);
for i = 1:numel(q_act_binmean_cat)
    q_act_binmean_pad{i}(use_bins{i}) = ...
        q_act_binmean_cat{i};
end

% Average sessions within animals
q_act_binmean_animal = nan(length(bin_edges)-1,length(data));
for curr_animal = 1:length(data);
    q_act_binmean_animal(:,curr_animal) = ....
        nanmean(horzcat(q_act_binmean_pad{curr_animal,:}),2);
end

% Average animals within sessions
q_act_binmean_day = nan(length(bin_edges)-1,14);
for curr_day = 1:14;
    q_act_binmean_day(:,curr_day) = ....
        nanmean(horzcat(q_act_binmean_pad{:,curr_day}),2);
end

% Plot all animals
figure;
sq = ceil(sqrt(length(data)));
for curr_animal = 1:length(data)
    subplot(sq,sq,curr_animal);hold on;
    set(gca,'ColorOrder',jet(14));
    curr_data = horzcat(q_act_binmean_pad{curr_animal,:});
    plot(curr_data);
    title(['Animal ' num2str(curr_animal)]);
end

% Plot session / animal average
figure; 
subplot(1,2,1);hold on
plot(bin_centers,q_act_binmean_animal,'color',[0.5,0.5,0.5]);
plot(bin_centers,nanmean(q_act_binmean_animal,2),'color','k','linewidth',2);
title('Average across sessions');
xlabel('Normalized lever position');
ylabel('Normalized quiescent activity');

subplot(1,2,2);hold on
set(gca,'ColorOrder',jet(size(q_act_binmean_day,2)));
plot(bin_centers,q_act_binmean_day);
plot(bin_centers,nanmean(q_act_binmean_day,2),'color','k','linewidth',2);
title('Average across animals');
xlabel('Normalized lever position')
ylabel('Normalized quiescent activity');

% Plot average across animals/sessions
q_act_binmean_pad_cat = horzcat(q_act_binmean_pad{:});
figure; hold on;
errorbar(bin_centers,nanmean(q_act_binmean_pad_cat,2), ...
    nanstd(q_act_binmean_pad_cat,[],2)./ ...
    sqrt(sum(~isnan(q_act_binmean_pad_cat),2)),'k','linewidth',2);
title('Average across all animals/sessions');
xlabel('Normalized lever position')
ylabel('Normalized quiescent activity');

%% Quiescent activity vs. BOTH lever position and velocity

q_act = cell(length(data),1);
m_act = cell(length(data),1);
u_act = cell(length(data),1);

lever_position = cell(length(data),1);
lever_velocity = cell(length(data),1);

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get quiescent activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_u = ~curr_q & ~curr_m;
        curr_act_q = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        curr_act_m = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_m,:) > 0);
        curr_act_u = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_u,:) > 0);
        
        % Get lever force        
        [lever_active,lever_force_resample, ...
            lever_force_smooth,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated( ...
            data(curr_animal).bhv(curr_day).lever_force);
        
        frame_times = (data(curr_animal).bhv(curr_day).frame_times( ...
            1:size(curr_act_q,2)))*1000;
        
        imaged_lever_idx = ((1:length(lever_active)) >= ...
            frame_times(1)) & ...
            ((1:length(lever_active)) <= ...
            ceil(frame_times(end)+median(diff(frame_times))));
        
        lever_active_imaged = lever_active(imaged_lever_idx);
        lever_position_resample_imaged = lever_force_resample(imaged_lever_idx);
        lever_velocity_resample_imaged = lever_velocity_envelope_smooth(imaged_lever_idx);
        
        relative_frame_times = ceil(frame_times - frame_times(1));
        lever_time_frames = nan(1,length(lever_position_resample_imaged));
        for curr_frame = 1:length(relative_frame_times);
            if curr_frame ~= length(relative_frame_times)
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    relative_frame_times(curr_frame+1)) = curr_frame;
            else
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    end) = curr_frame;
            end
        end
              
        % Split lever by frames
        frame_boundaries = [0,find(lever_time_frames(2:end) ~= ...
            lever_time_frames(1:end-1)),length(lever_time_frames)];
        lever_position_frame_split = mat2cell(lever_position_resample_imaged, ...
            diff(frame_boundaries)',1);
        lever_velocity_frame_split = mat2cell(lever_velocity_resample_imaged, ...
            diff(frame_boundaries)',1);
        
        lever_active_split = mat2cell(lever_active_imaged, ...
            diff(frame_boundaries)',1);
        
        % Get average lever parameter in each frame
        lever_position_frame = cellfun(@mean,lever_position_frame_split);
        lever_velocity_frame = cellfun(@mean,lever_velocity_frame_split);

        % Get median lever of entire day
        median_lever = median(lever_force_smooth);
        
        % Store smoothed lever and activity
        q_act{curr_animal}{curr_day} = mean(curr_act_q,1);
        m_act{curr_animal}{curr_day} = mean(curr_act_m,1);
        u_act{curr_animal}{curr_day} = mean(curr_act_u,1);
        lever_position{curr_animal}{curr_day} = lever_position_frame-median_lever;
        lever_velocity{curr_animal}{curr_day} = lever_velocity_frame;
    end
    disp(curr_animal);
end

% Bin lever and get groupstats of activity
lever_position_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),lever_position,'uni',false);
lever_velocity_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),lever_velocity,'uni',false);
q_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),q_act,'uni',false);
m_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),m_act,'uni',false);
u_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),u_act,'uni',false);

position_bin_edges = linspace(-3,3,20);
position_bin_centers = position_bin_edges(1:end-1) + diff(position_bin_edges)/2;
position_bin_edges(1) = -Inf;
position_bin_edges(end) = Inf;
[~,lever_position_bins] = cellfun(@(x) cellfun(@(x) histc(x,position_bin_edges), ...
    x,'uni',false),lever_position_norm,'uni',false);

velocity_bin_edges = linspace(0,3,20);
velocity_bin_centers = velocity_bin_edges(1:end-1) + diff(velocity_bin_edges)/2;
velocity_bin_edges(end) = Inf;
[~,lever_velocity_bins] = cellfun(@(x) cellfun(@(x) histc(x,velocity_bin_edges), ...
    x,'uni',false),lever_velocity_norm,'uni',false);

[q_act_binmean,q_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    q_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

[m_act_binmean,m_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    m_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

[u_act_binmean,u_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    u_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

q_act_binmean_idxd = cell(length(data),14);
m_act_binmean_idxd = cell(length(data),14);
u_act_binmean_idxd = cell(length(data),14);
for curr_animal = 1:length(data)
      
    % place bin stats into correct indicies
    curr_bin_cat = cellfun(@(x) cellfun(@str2num,x),q_act_bins{curr_animal},'uni',false);
    curr_bin_idx = cellfun(@(x) sub2ind([length(position_bin_centers), ...
        length(velocity_bin_centers)],x(:,1),x(:,2)),curr_bin_cat,'uni',false);
    
    curr_q_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),q_act_binmean{curr_animal},'uni',false);
    curr_m_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),m_act_binmean{curr_animal},'uni',false);
    curr_u_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),u_act_binmean{curr_animal},'uni',false);
    
    for curr_day = 1:length(curr_q_act_binmean_idxd)
       curr_q_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           q_act_binmean{curr_animal}{curr_day};
       
       curr_m_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           m_act_binmean{curr_animal}{curr_day};
       
       curr_u_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           u_act_binmean{curr_animal}{curr_day};
    end
    
    q_act_binmean_idxd(curr_animal,1:length(q_act_binmean{curr_animal})) = ...
        curr_q_act_binmean_idxd;
    
    m_act_binmean_idxd(curr_animal,1:length(m_act_binmean{curr_animal})) = ...
        curr_m_act_binmean_idxd;
    
    u_act_binmean_idxd(curr_animal,1:length(u_act_binmean{curr_animal})) = ...
        curr_u_act_binmean_idxd;
        
end

% Plot average across all animals/sessions
q_act_binmean_allmean = nanmean(cat(3,q_act_binmean_idxd{:}),3);
m_act_binmean_allmean = nanmean(cat(3,m_act_binmean_idxd{:}),3);
u_act_binmean_allmean = nanmean(cat(3,u_act_binmean_idxd{:}),3);

figure; 

subplot(1,3,1);
imagesc(q_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Quiescent ROI activity');

subplot(1,3,2);
imagesc(m_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Movement ROI activity');

subplot(1,3,3);
imagesc(u_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Unclassified ROI activity');


% Bin along each dimension independently

% Position
[q_act_pbinmean,q_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    q_act_norm,lever_position_bins,'uni',false);
q_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),q_act_pbinmean,'uni',false);
for curr_animal = 1:length(q_act_pbinmean_idxd)
    for curr_session = 1:length(q_act_pbinmean_idxd{curr_animal})
        q_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(q_act_pbins{curr_animal}{curr_session})) = ...
            q_act_pbinmean{curr_animal}{curr_session};
    end
end

[m_act_pbinmean,m_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    m_act_norm,lever_position_bins,'uni',false);
m_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),m_act_pbinmean,'uni',false);
for curr_animal = 1:length(m_act_pbinmean_idxd)
    for curr_session = 1:length(m_act_pbinmean_idxd{curr_animal})
        m_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(m_act_pbins{curr_animal}{curr_session})) = ...
            m_act_pbinmean{curr_animal}{curr_session};
    end
end

[u_act_pbinmean,u_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    u_act_norm,lever_position_bins,'uni',false);
u_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),u_act_pbinmean,'uni',false);
for curr_animal = 1:length(u_act_pbinmean_idxd)
    for curr_session = 1:length(u_act_pbinmean_idxd{curr_animal})
        u_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(u_act_pbins{curr_animal}{curr_session})) = ...
            u_act_pbinmean{curr_animal}{curr_session};
    end
end

% Velocity
[q_act_vbinmean,q_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    q_act_norm,lever_velocity_bins,'uni',false);
q_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),q_act_vbinmean,'uni',false);
for curr_animal = 1:length(q_act_vbinmean_idxd)
    for curr_session = 1:length(q_act_vbinmean_idxd{curr_animal})
        q_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(q_act_vbins{curr_animal}{curr_session})) = ...
            q_act_vbinmean{curr_animal}{curr_session};
    end
end

[m_act_vbinmean,m_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    m_act_norm,lever_velocity_bins,'uni',false);
m_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),m_act_vbinmean,'uni',false);
for curr_animal = 1:length(m_act_vbinmean_idxd)
    for curr_session = 1:length(m_act_vbinmean_idxd{curr_animal})
        m_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(m_act_vbins{curr_animal}{curr_session})) = ...
            m_act_vbinmean{curr_animal}{curr_session};
    end
end

[u_act_vbinmean,u_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    u_act_norm,lever_velocity_bins,'uni',false);
u_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),u_act_vbinmean,'uni',false);
for curr_animal = 1:length(u_act_vbinmean_idxd)
    for curr_session = 1:length(u_act_vbinmean_idxd{curr_animal})
        u_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(u_act_vbins{curr_animal}{curr_session})) = ...
            u_act_vbinmean{curr_animal}{curr_session};
    end
end


figure;
subplot(2,3,1)
plot(position_bin_centers,nanmean(cell2mat(horzcat(q_act_pbinmean_idxd{:})),2),'k');
ylabel('Quiescent ROI activity')
xlabel('Position');

subplot(2,3,2)
plot(position_bin_centers,nanmean(cell2mat(horzcat(m_act_pbinmean_idxd{:})),2),'k');
ylabel('Movement ROI activity')
xlabel('Position');

subplot(2,3,3)
plot(position_bin_centers,nanmean(cell2mat(horzcat(u_act_pbinmean_idxd{:})),2),'k');
ylabel('Unclassified ROI activity')
xlabel('Position');

subplot(2,3,4)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(q_act_vbinmean_idxd{:})),2),'k');
ylabel('Quiescent ROI activity')
xlabel('Velocity');

subplot(2,3,5)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(m_act_vbinmean_idxd{:})),2),'k');
ylabel('Movement ROI activity')
xlabel('Velocity');

subplot(2,3,6)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(u_act_vbinmean_idxd{:})),2),'k');
ylabel('Unclassified ROI activity')
xlabel('Velocity');

%% FOR GENERAL USE: get frame-by-frame position, velocity, speed
% (now a part of normal analysis, so not needed seperately)

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        num_frames = min(size(data(curr_animal).im(curr_day).roi_trace_df,2), ...
            length(data(curr_animal).bhv(curr_day).frame_times));
        
        % Get lever force        
        [lever_active,lever_force_resample, ...
            lever_force_smooth,lever_velocity_envelope_smooth] = ...
            AP_parseLeverMovement_updated( ...
            data(curr_animal).bhv(curr_day).lever_force);
        
        lever_velocity = [0;smooth(diff(lever_force_smooth),5)];
        
        frame_times = (data(curr_animal).bhv(curr_day).frame_times( ...
            1:num_frames))*1000;
        
        imaged_lever_idx = ((1:length(lever_active)) >= ...
            frame_times(1)) & ...
            ((1:length(lever_active)) <= ...
            ceil(frame_times(end)+median(diff(frame_times))));
        
        lever_active_imaged = lever_active(imaged_lever_idx);
        lever_position_resample_imaged = lever_force_resample(imaged_lever_idx);
        lever_velocity_resample_imaged = lever_velocity(imaged_lever_idx);
        lever_speed_resample_imaged = lever_velocity_envelope_smooth(imaged_lever_idx);
        
        relative_frame_times = ceil(frame_times - frame_times(1));
        lever_time_frames = nan(1,length(lever_position_resample_imaged));
        for curr_frame = 1:length(relative_frame_times);
            if curr_frame ~= length(relative_frame_times)
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    relative_frame_times(curr_frame+1)) = curr_frame;
            else
                lever_time_frames((relative_frame_times(curr_frame)+1): ...
                    end) = curr_frame;
            end
        end
              
        % Split lever by frames
        frame_boundaries = [0,find(lever_time_frames(2:end) ~= ...
            lever_time_frames(1:end-1)),length(lever_time_frames)];
        lever_position_frame_split = mat2cell(lever_position_resample_imaged, ...
            diff(frame_boundaries)',1);
        lever_velocity_frame_split = mat2cell(lever_velocity_resample_imaged, ...
            diff(frame_boundaries)',1);
        lever_speed_frame_split = mat2cell(lever_speed_resample_imaged, ...
            diff(frame_boundaries)',1);

        % Get average lever parameter in each frame
        lever_position_frame = cellfun(@mean,lever_position_frame_split);
        lever_velocity_frame = cellfun(@mean,lever_velocity_frame_split);
        lever_speed_frame = cellfun(@mean,lever_speed_frame_split);

        % Get median lever of entire day
        median_lever = median(lever_force_smooth);
        
        % Store lever position, velocity, speed in frames
        analysis(curr_animal).lever(curr_day).lever_position_frames = ...
            lever_position_frame-median_lever;
        analysis(curr_animal).lever(curr_day).lever_velocity_frames = ...
            lever_velocity_frame;
        analysis(curr_animal).lever(curr_day).lever_speed_frames = ...
            lever_speed_frame;
        
    end
    disp(curr_animal);
end

%% Quiescent activity vs. BOTH lever position and velocity (FASTER)

% Get average classified activity
m_act = cellfun(@(act,class) cellfun(@(act,class) ...
    nanmean(act(class,:) > 0,1),act,mat2cell(class,size(class,1), ...
    ones(size(class,2),1)),'uni',false), ...
    cellfun(@(x) {x.roi_trace_thresh},{data.im},'uni',false), ...
    {classified_rois.movement},'uni',false)';

q_act = cellfun(@(act,class) cellfun(@(act,class) ...
    nanmean(act(class,:) > 0,1),act,mat2cell(class,size(class,1), ...
    ones(size(class,2),1)),'uni',false), ...
    cellfun(@(x) {x.roi_trace_thresh},{data.im},'uni',false), ...
    {classified_rois.quiescent},'uni',false)';

u_act = cellfun(@(act,class) cellfun(@(act,class) ...
    nanmean(act(class,:) > 0,1),act,mat2cell(class,size(class,1), ...
    ones(size(class,2),1)),'uni',false), ...
    cellfun(@(x) {x.roi_trace_thresh},{data.im},'uni',false), ...
    {classified_rois.unclassified_active},'uni',false)';

% Get out lever position/velocity in frames
lever_position = cellfun(@(x) {x.lever_position_frames},{analysis.lever},'uni',false)';
lever_velocity = cellfun(@(x) {x.lever_velocity_frames},{analysis.lever},'uni',false)';

% Smooth the lever position (to get an approximated "set point")
lever_position_smooth = cellfun(@(x) cellfun(@(x) smooth(x,100,'loess'), ...
    x,'uni',false),lever_position,'uni',false);

% Bin lever and get groupstats of activity
lever_position_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),lever_position_smooth,'uni',false);
lever_velocity_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),lever_velocity,'uni',false);
q_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),q_act,'uni',false);
m_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),m_act,'uni',false);
u_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),u_act,'uni',false);

position_bin_edges = linspace(-3,3,20);
position_bin_centers = position_bin_edges(1:end-1) + diff(position_bin_edges)/2;
position_bin_edges(1) = -Inf;
position_bin_edges(end) = Inf;
[~,lever_position_bins] = cellfun(@(x) cellfun(@(x) histc(x,position_bin_edges), ...
    x,'uni',false),lever_position_norm,'uni',false);

velocity_bin_edges = linspace(-3,3,20);
velocity_bin_centers = velocity_bin_edges(1:end-1) + diff(velocity_bin_edges)/2;
velocity_bin_edges(1) = -Inf;
velocity_bin_edges(end) = Inf;
[~,lever_velocity_bins] = cellfun(@(x) cellfun(@(x) histc(x,velocity_bin_edges), ...
    x,'uni',false),lever_velocity_norm,'uni',false);

[q_act_binmean,q_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    q_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

[m_act_binmean,m_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    m_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

[u_act_binmean,u_act_bins] = cellfun(@(act,p_bins,v_bins) cellfun(@(act,p_bins,v_bins) ...
    grpstats(act,{p_bins,v_bins},{'nanmean','gname'}),act,p_bins,v_bins,'uni',false), ...
    u_act_norm,lever_position_bins,lever_velocity_bins,'uni',false);

q_act_binmean_idxd = cell(length(data),14);
m_act_binmean_idxd = cell(length(data),14);
u_act_binmean_idxd = cell(length(data),14);
for curr_animal = 1:length(data)
      
    % place bin stats into correct indicies
    curr_bin_cat = cellfun(@(x) cellfun(@str2num,x),q_act_bins{curr_animal},'uni',false);
    curr_bin_idx = cellfun(@(x) sub2ind([length(position_bin_centers), ...
        length(velocity_bin_centers)],x(:,1),x(:,2)),curr_bin_cat,'uni',false);
    
    curr_q_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),q_act_binmean{curr_animal},'uni',false);
    curr_m_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),m_act_binmean{curr_animal},'uni',false);
    curr_u_act_binmean_idxd = cellfun(@(x) nan(length(position_bin_centers), ...
        length(velocity_bin_centers)),u_act_binmean{curr_animal},'uni',false);
    
    for curr_day = 1:length(curr_q_act_binmean_idxd)
       curr_q_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           q_act_binmean{curr_animal}{curr_day};
       
       curr_m_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           m_act_binmean{curr_animal}{curr_day};
       
       curr_u_act_binmean_idxd{curr_day}(curr_bin_idx{curr_day}) = ...
           u_act_binmean{curr_animal}{curr_day};
    end
    
    q_act_binmean_idxd(curr_animal,1:length(q_act_binmean{curr_animal})) = ...
        curr_q_act_binmean_idxd;
    
    m_act_binmean_idxd(curr_animal,1:length(m_act_binmean{curr_animal})) = ...
        curr_m_act_binmean_idxd;
    
    u_act_binmean_idxd(curr_animal,1:length(u_act_binmean{curr_animal})) = ...
        curr_u_act_binmean_idxd;
        
end

% Plot average across all animals/sessions
q_act_binmean_allmean = nanmean(cat(3,q_act_binmean_idxd{:}),3);
m_act_binmean_allmean = nanmean(cat(3,m_act_binmean_idxd{:}),3);
u_act_binmean_allmean = nanmean(cat(3,u_act_binmean_idxd{:}),3);

figure; 

subplot(1,3,1);
imagesc(q_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Quiescent ROI activity');

subplot(1,3,2);
imagesc(m_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Movement ROI activity');

subplot(1,3,3);
imagesc(u_act_binmean_allmean);
colormap(gray);
ylabel('Normalized lever position (stds)');
set(gca,'YTickLabel',round(position_bin_centers(get(gca,'YTick'))*10)/10);
xlabel('Normalized lever velocity (stds)');
set(gca,'XTickLabel',round(velocity_bin_centers(get(gca,'XTick'))*10)/10);
title('Unclassified ROI activity');


% Bin along each dimension independently

% Position
[q_act_pbinmean,q_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    q_act_norm,lever_position_bins,'uni',false);
q_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),q_act_pbinmean,'uni',false);
for curr_animal = 1:length(q_act_pbinmean_idxd)
    for curr_session = 1:length(q_act_pbinmean_idxd{curr_animal})
        q_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(q_act_pbins{curr_animal}{curr_session})) = ...
            q_act_pbinmean{curr_animal}{curr_session};
    end
end

[m_act_pbinmean,m_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    m_act_norm,lever_position_bins,'uni',false);
m_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),m_act_pbinmean,'uni',false);
for curr_animal = 1:length(m_act_pbinmean_idxd)
    for curr_session = 1:length(m_act_pbinmean_idxd{curr_animal})
        m_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(m_act_pbins{curr_animal}{curr_session})) = ...
            m_act_pbinmean{curr_animal}{curr_session};
    end
end

[u_act_pbinmean,u_act_pbins] = cellfun(@(act,p_bins) cellfun(@(act,p_bins) ...
    grpstats(act,{p_bins},{'nanmean','gname'}),act,p_bins,'uni',false), ...
    u_act_norm,lever_position_bins,'uni',false);
u_act_pbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(position_bin_centers),1), ...
    x,'uni',false),u_act_pbinmean,'uni',false);
for curr_animal = 1:length(u_act_pbinmean_idxd)
    for curr_session = 1:length(u_act_pbinmean_idxd{curr_animal})
        u_act_pbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(u_act_pbins{curr_animal}{curr_session})) = ...
            u_act_pbinmean{curr_animal}{curr_session};
    end
end

% Velocity
[q_act_vbinmean,q_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    q_act_norm,lever_velocity_bins,'uni',false);
q_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),q_act_vbinmean,'uni',false);
for curr_animal = 1:length(q_act_vbinmean_idxd)
    for curr_session = 1:length(q_act_vbinmean_idxd{curr_animal})
        q_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(q_act_vbins{curr_animal}{curr_session})) = ...
            q_act_vbinmean{curr_animal}{curr_session};
    end
end

[m_act_vbinmean,m_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    m_act_norm,lever_velocity_bins,'uni',false);
m_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),m_act_vbinmean,'uni',false);
for curr_animal = 1:length(m_act_vbinmean_idxd)
    for curr_session = 1:length(m_act_vbinmean_idxd{curr_animal})
        m_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(m_act_vbins{curr_animal}{curr_session})) = ...
            m_act_vbinmean{curr_animal}{curr_session};
    end
end

[u_act_vbinmean,u_act_vbins] = cellfun(@(act,v_bins) cellfun(@(act,v_bins) ...
    grpstats(act,{v_bins},{'nanmean','gname'}),act,v_bins,'uni',false), ...
    u_act_norm,lever_velocity_bins,'uni',false);
u_act_vbinmean_idxd = cellfun(@(x) cellfun(@(x) nan(length(velocity_bin_centers),1), ...
    x,'uni',false),u_act_vbinmean,'uni',false);
for curr_animal = 1:length(u_act_vbinmean_idxd)
    for curr_session = 1:length(u_act_vbinmean_idxd{curr_animal})
        u_act_vbinmean_idxd{curr_animal}{curr_session}( ...
            str2double(u_act_vbins{curr_animal}{curr_session})) = ...
            u_act_vbinmean{curr_animal}{curr_session};
    end
end


figure;
subplot(2,3,1)
plot(position_bin_centers,nanmean(cell2mat(horzcat(q_act_pbinmean_idxd{:})),2),'k');
ylabel('Quiescent ROI activity')
xlabel('Position');

subplot(2,3,2)
plot(position_bin_centers,nanmean(cell2mat(horzcat(m_act_pbinmean_idxd{:})),2),'k');
ylabel('Movement ROI activity')
xlabel('Position');

subplot(2,3,3)
plot(position_bin_centers,nanmean(cell2mat(horzcat(u_act_pbinmean_idxd{:})),2),'k');
ylabel('Unclassified ROI activity')
xlabel('Position');

subplot(2,3,4)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(q_act_vbinmean_idxd{:})),2),'k');
ylabel('Quiescent ROI activity')
xlabel('Velocity');

subplot(2,3,5)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(m_act_vbinmean_idxd{:})),2),'k');
ylabel('Movement ROI activity')
xlabel('Velocity');

subplot(2,3,6)
plot(velocity_bin_centers,nanmean(cell2mat(horzcat(u_act_vbinmean_idxd{:})),2),'k');
ylabel('Unclassified ROI activity')
xlabel('Velocity');

%% Quiescent activity vs. movement activity surrounding epochs

q_epoch_act = cell(length(data),1);
m_epoch_act = cell(length(data),1);

for curr_animal = 1:length(data)
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_act_q = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_q,:) > 0);
        
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act_m = double(data(curr_animal).im(curr_day).roi_trace_thresh(curr_m,:) > 0);
        
        % Get lever position/velocity by frame   
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;

        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
        
        frames_split = mat2cell(transpose(1:length(lever_active_frames)),diff(boundary_times));
        
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
        
        % Used to look at adjacent movement epochs
        lever_split_q_idx = find(lever_split_q);
        lever_split_q_idx_use = lever_split_q_idx(lever_split_q_idx > 1 & ...
            lever_split_q_idx < length(lever_split_q));
        lever_split_pre_m_idx = lever_split_q_idx_use - 1;
        lever_split_post_m_idx = lever_split_q_idx_use - 1;
        
        split_q_act = cellfun(@(frames) curr_act_q(:,frames), ...
            frames_split(lever_split_q_idx_use),'uni',false);
        split_q_act_mean = cellfun(@(x) nanmean(x(:)),split_q_act);
        
        split_m_act = cellfun(@(frames) curr_act_m(:,frames), ...
            frames_split(lever_split_post_m_idx),'uni',false);
        split_m_act_mean = cellfun(@(x) nanmean(x(:)),split_m_act);       
        
        % Store quiescent epoch lever and activity
        q_epoch_act{curr_animal}{curr_day} = split_q_act_mean;
        m_epoch_act{curr_animal}{curr_day} = split_m_act_mean;
        
    end
    disp(curr_animal);
end

% Bin lever and get groupstats of activity
m_epoch_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),m_epoch_act,'uni',false);
q_epoch_act_norm = cellfun(@(x) cellfun(@(x) x./std(x), ...
    x,'uni',false),q_epoch_act,'uni',false);

bin_edges = linspace(0,3,10);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_edges(end) = Inf;
bin_edges(1) = -Inf;
[~,m_epoch_act_norm_bin] = cellfun(@(x) cellfun(@(x) histc(x,bin_edges), ...
    x,'uni',false),m_epoch_act_norm,'uni',false);
q_epoch_act_binmean = cellfun(@(act,bins) cellfun(@(act,bins) ...
    grpstats(act,bins,'nanmean'),act,bins,'uni',false), ...
    q_epoch_act_norm,m_epoch_act_norm_bin,'uni',false);

q_epoch_act_binmean_cat = cell(length(data),14);
m_epoch_act_bin_cat = cell(length(data),14);
for curr_animal = 1:length(data)
    q_epoch_act_binmean_cat(curr_animal,1:length(q_epoch_act_binmean{curr_animal})) = ...
        q_epoch_act_binmean{curr_animal};
    
    m_epoch_act_bin_cat(curr_animal,1:length(m_epoch_act_norm_bin{curr_animal})) = ...
        m_epoch_act_norm_bin{curr_animal};
end
use_bins = cellfun(@unique,m_epoch_act_bin_cat,'uni',false);

q_epoch_act_binmean_pad = cellfun(@(x) nan(length(bin_edges)-1,1), ...
    q_epoch_act_binmean_cat,'uni',false);
for i = 1:numel(q_epoch_act_binmean_cat)
    q_epoch_act_binmean_pad{i}(use_bins{i}) = ...
        q_epoch_act_binmean_cat{i};
end

% Plot average across animals/sessions
q_epoch_act_binmean_pad_cat = horzcat(q_epoch_act_binmean_pad{:});
figure; hold on;
errorbar(bin_centers,nanmean(q_epoch_act_binmean_pad_cat,2), ...
    nanstd(q_epoch_act_binmean_pad_cat,[],2)./ ...
    sqrt(sum(~isnan(q_epoch_act_binmean_pad_cat),2)),'k','linewidth',2);
title('Average across all animals/sessions');
xlabel('Normalized movement activity')
ylabel('Normalized quiescent activity');


%% Distribution of activity relative to position and velocity/speed

% Get average position and velocity when each ROI is active
roi_position = cell(length(data),1);
roi_position_shuff = cell(length(data),1);

roi_velocity = cell(length(data),1);
roi_velocity_shuff = cell(length(data),1);

roi_speed = cell(length(data),1);
roi_speed_shuff = cell(length(data),1);

for curr_animal = 1:length(data);
    n_days = length(data(curr_animal).im);
    n_rois = size(data(curr_animal).im(1).roi_trace_thresh,1);
    
    roi_position{curr_animal} = cell(n_days,1);
    roi_position_shuff{curr_animal} = cell(n_days,1);

    roi_velocity{curr_animal} = cell(n_days,1);
    roi_velocity_shuff{curr_animal} = cell(n_days,1);
    
    roi_speed{curr_animal} = cell(n_days,1);
    roi_speed_shuff{curr_animal} = cell(n_days,1);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = smooth(analysis(curr_animal). ...
            lever(curr_day).lever_position_frames,100,'loess');
        lever_velocity_frames = analysis(curr_animal).lever(curr_day).lever_velocity_frames;
        lever_speed_frames = analysis(curr_animal).lever(curr_day).lever_speed_frames;
        
        % Normalize lever parameters
        lever_position_norm = lever_position_frames./std(lever_position_frames);
        lever_velocity_norm = lever_velocity_frames./std(lever_velocity_frames);
        lever_speed_norm = lever_speed_frames./std(lever_speed_frames);
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        % Only use "active" ROIs, at least 5 events
        min_events = 5;
        active_rois = sum(diff(curr_act > 0,[],2) == 1,2) > 5;
        curr_act = curr_act(active_rois,:);
        
        curr_roi_position = (curr_act*lever_position_norm)./sum(curr_act,2);
        curr_roi_velocity = (curr_act*lever_velocity_norm)./sum(curr_act,2);
        curr_roi_speed = (curr_act*lever_speed_norm)./sum(curr_act,2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_position_split = mat2cell(lever_position_norm,diff(boundary_frames));
        lever_velocity_split = mat2cell(lever_velocity_norm,diff(boundary_frames));
        lever_speed_split = mat2cell(lever_speed_norm,diff(boundary_frames));
        
        num_rep = 1000;
        shuffle_perms = shake(repmat([1:length(boundary_frames)-1]',1,num_rep),1);
        lever_position_shuffle = nan(length(lever_position_norm),num_rep);
        lever_velocity_shuffle = nan(length(lever_velocity_norm),num_rep);
        lever_speed_shuffle = nan(length(lever_speed_norm),num_rep);
        for i = 1:num_rep
            lever_position_shuffle(:,i) = ...
                vertcat(lever_position_split{shuffle_perms(:,i)});
            lever_velocity_shuffle(:,i) = ...
                vertcat(lever_velocity_split{shuffle_perms(:,i)});
            lever_speed_shuffle(:,i) = ...
                vertcat(lever_speed_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        curr_roi_position_shuff = bsxfun(@times,curr_act*lever_position_shuffle,1./sum(curr_act,2));
        curr_roi_velocity_shuff = bsxfun(@times,curr_act*lever_velocity_shuffle,1./sum(curr_act,2));
        curr_roi_speed_shuff = bsxfun(@times,curr_act*lever_speed_shuffle,1./sum(curr_act,2));
        
        % Store
        roi_position{curr_animal}{curr_day} = curr_roi_position;
        roi_velocity{curr_animal}{curr_day} = curr_roi_velocity;
        roi_speed{curr_animal}{curr_day} = curr_roi_speed;
        
        roi_position_shuff{curr_animal}{curr_day} = curr_roi_position_shuff(:);
        roi_velocity_shuff{curr_animal}{curr_day} = curr_roi_velocity_shuff(:);
        roi_speed_shuff{curr_animal}{curr_day} = curr_roi_speed_shuff(:);
        
    end
    disp(curr_animal);
end

figure; 

% Histogram of real vs. shuffled position (all concatenated)
roi_position_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_position,'uni',false)));
roi_position_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_position_shuff,'uni',false)));

bin_edges = linspace(-2,2,50);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(roi_position_cat,bin_centers)./length(roi_position_cat);
shuff_dist = hist(roi_position_shuff_cat,bin_centers)./length(roi_position_shuff_cat);

subplot(2,3,1); hold on;
plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
legend({'Shuffled','Real'})

ylabel('Number of ROIs');
xlabel('Lever displacement from median (stds)');

% Histogram of real vs. shuffled velocity (all concatenated)
roi_velocity_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_velocity,'uni',false)));
roi_velocity_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_velocity_shuff,'uni',false)));

bin_edges = linspace(-0.5,0.5,50);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(roi_velocity_cat,bin_centers)./length(roi_velocity_cat);
shuff_dist = hist(roi_velocity_shuff_cat,bin_centers)./length(roi_velocity_shuff_cat);

subplot(2,3,2); hold on;
plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
legend({'Shuffled','Real'})

ylabel('Number of ROIs');
xlabel('Lever velocity (stds)');

% Histogram of real vs. shuffled speed (all concatenated)
roi_speed_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_speed,'uni',false)));
roi_speed_shuff_cat = cell2mat((cellfun(@(x) vertcat(x{:}),roi_speed_shuff,'uni',false)));

bin_edges = linspace(0,1,50);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

real_dist = hist(roi_speed_cat,bin_centers)./length(roi_speed_cat);
shuff_dist = hist(roi_speed_shuff_cat,bin_centers)./length(roi_speed_shuff_cat);

subplot(2,3,3); hold on;
plot(bin_centers,shuff_dist,'k');
plot(bin_centers,real_dist,'r');
legend({'Shuffled','Real'})

ylabel('Number of ROIs');
xlabel('Lever speed (stds)');

% 2D histogram of real vs. shuffled position and velocity (all concatenated)
p_bin_edges = linspace(-1,1,20);
p_bin_centers = p_bin_edges(1:end-1) + diff(p_bin_edges)/2;
p_bin_edges(1) = -Inf;
p_bin_edges(end) = Inf;

v_bin_edges = linspace(-0.5,0.5,20);
v_bin_centers = v_bin_edges(1:end-1) + diff(v_bin_edges)/2;
v_bin_edges(1) = -Inf;
v_bin_edges(end) = Inf;

real_dist = hist3([roi_position_cat,roi_velocity_cat],'Edges', ...
    {p_bin_edges,v_bin_edges})./length(roi_velocity_cat);
shuff_dist = hist3([roi_position_shuff_cat,roi_velocity_shuff_cat], ...
    'Edges',{p_bin_edges,v_bin_edges})./length(roi_velocity_shuff_cat);

subplot(2,3,4);
imagesc(real_dist - shuff_dist);
ylabel('Lever displacement from median (stds)');
xlabel('Lever velocity (stds)');
set(gca,'XTick',1:2:length(v_bin_centers));
set(gca,'XTickLabel',round(v_bin_centers(get(gca,'XTick'))*10)/10);
xlim([0.5,length(v_bin_centers)+0.5]);
set(gca,'YTick',1:2:length(p_bin_centers));
set(gca,'YTickLabel',round(p_bin_centers(get(gca,'YTick'))*10)/10);
ylim([0.5,length(p_bin_centers)+0.5]);
colormap(gray);
title('Real - Shuffled');

% 2D histogram of real vs. shuffled position and speed (all concatenated)
p_bin_edges = linspace(-1,1,20);
p_bin_centers = p_bin_edges(1:end-1) + diff(p_bin_edges)/2;
p_bin_edges(1) = -Inf;
p_bin_edges(end) = Inf;

s_bin_edges = linspace(0,1,20);
s_bin_centers = s_bin_edges(1:end-1) + diff(s_bin_edges)/2;
s_bin_edges(1) = -Inf;
s_bin_edges(end) = Inf;

real_dist = hist3([roi_position_cat,roi_speed_cat],'Edges', ...
    {p_bin_edges,s_bin_edges})./length(roi_speed_cat);
shuff_dist = hist3([roi_position_shuff_cat,roi_speed_shuff_cat], ...
    'Edges',{p_bin_edges,s_bin_edges})./length(roi_speed_shuff_cat);

subplot(2,3,5);
imagesc(real_dist - shuff_dist);
ylabel('Lever displacement from median (stds)');
xlabel('Lever speed (stds)');
set(gca,'XTick',1:2:length(s_bin_centers));
set(gca,'XTickLabel',round(s_bin_centers(get(gca,'XTick'))*10)/10);
xlim([0.5,length(s_bin_centers)+0.5]);
set(gca,'YTick',1:2:length(p_bin_centers));
set(gca,'YTickLabel',round(p_bin_centers(get(gca,'YTick'))*10)/10);
ylim([0.5,length(p_bin_centers)+0.5]);
colormap(gray);
title('Real - Shuffled');

%% Classify ROIs according to speed bins

% Speed bins
bin_edges = linspace(0,1,5);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_edges(1) = -Inf;
bin_edges(end) = Inf;

for curr_animal = 1:length(data);
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    n_days = length(data(curr_animal).im);
    
    classified_rois(curr_animal).speed = false( ...
        length(bin_centers),n_rois,n_days);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_speed_frames = analysis(curr_animal).lever(curr_day).lever_speed_frames;
        
        % Normalize lever parameters
        lever_speed_norm = lever_speed_frames./std(lever_speed_frames);
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);
        
        % Bin velocity      
        [~,lever_speed_norm_bin] = ...
            histc(lever_speed_norm,bin_edges);
        
        roi_speed_bin = cell2mat(arrayfun(@(x) grpstats(curr_act(x,:), ...
            lever_speed_norm_bin),1:size(curr_act,1),'uni',false));
        
        roi_speed = (curr_act*lever_speed_norm)./sum(curr_act,2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_speed_norm_bin_split = mat2cell(lever_speed_norm_bin,diff(boundary_frames));
        
        num_rep = 10000;
        shuffle_perms = shake(repmat([1:length(boundary_frames)-1]',1,num_rep),1);
        lever_speed_norm_bin_shuffle = nan(length(lever_speed_norm),num_rep);
        for i = 1:num_rep
            lever_speed_norm_bin_shuffle(:,i) = ...
                vertcat(lever_speed_norm_bin_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        roi_speed_bin_shuff = nan(length(bin_centers),size(curr_act,1),num_rep);
        for curr_bin = 1:length(bin_centers);           
            roi_speed_bin_shuff(curr_bin,:,:) = ...
                permute(bsxfun(@times,curr_act*(lever_speed_norm_bin_shuffle == curr_bin), ...
                1./sum(lever_speed_norm_bin_shuffle == curr_bin,1)),[3,1,2]);           
        end
        
        roi_speed_bin_rank = permute(tiedrank([permute(roi_speed_bin,[3,2,1]); ...
            permute(roi_speed_bin_shuff,[3,2,1])]),[3,2,1]);
        
        roi_speed_bin_p = roi_speed_bin_rank(:,:,1)./(num_rep+1);
        
        
        % Store
        classified_rois(curr_animal).speed(:,:,curr_day) = ...
            roi_speed_bin_p > 0.975;
        
    end
    disp(curr_animal);
end


% Summary figure
figure;

% Fraction of speed-related ROIs
speed_roi_frac = cell2mat(permute(cellfun(@(x) ...
    permute(nanmean(x,2),[1,3,2]),{classified_rois.speed},'uni',false),[1,3,2]));

subplot(1,2,1);
plot(nanmean(speed_roi_frac,3)');
legend({'Speed 1','Speed 2','Speed 3','Speed 4+'});
ylabel('Fraction ROIs');
xlabel('Day');

% Fraction of ROI classification overlap

% quiescent / speed 1
q_s1 = cellfun(@(q,s)  [sum(q & ~s), sum(q & s), sum(~q & s)]/sum(q | s), ...
    cellfun(@(x) x(:),{classified_rois.quiescent},'uni',false), ...
    cellfun(@(x) reshape(permute(x(1,:,:),[2,3,1]),[],1), ...
    {classified_rois.speed},'uni',false),'uni',false);

% quiescent / speed 2-4
q_s24 = cellfun(@(q,s)  [sum(q & ~s), sum(q & s), sum(~q & s)]/sum(q | s), ...
    cellfun(@(x) x(:),{classified_rois.quiescent},'uni',false), ...
    cellfun(@(x) reshape(permute(any(x(2:end,:,:),1),[2,3,1]),[],1), ...
    {classified_rois.speed},'uni',false),'uni',false);

% movement / speed 4
m_s4 = cellfun(@(m,s)  [sum(m & ~s), sum(m & s), sum(~m & s)]/sum(m | s), ...
    cellfun(@(x) x(:),{classified_rois.movement},'uni',false), ...
    cellfun(@(x) reshape(permute(x(4,:,:),[2,3,1]),[],1), ...
    {classified_rois.speed},'uni',false),'uni',false);

% movement / speed 1-3
m_s13 = cellfun(@(m,s)  [sum(m & ~s), sum(m & s), sum(~m & s)]/sum(m | s), ...
    cellfun(@(x) x(:),{classified_rois.movement},'uni',false), ...
    cellfun(@(x) reshape(permute(any(x(1:end-1,:,:),1),[2,3,1]),[],1), ...
    {classified_rois.speed},'uni',false),'uni',false);

% unclassified active / any speed
u_s = cellfun(@(u,s)  [sum(u & ~s), sum(u & s), sum(~u & s)]/sum(u | s), ...
    cellfun(@(x) x(:),{classified_rois.unclassified_active},'uni',false), ...
    cellfun(@(x) reshape(permute(any(x,1),[2,3,1]),[],1), ...
    {classified_rois.speed},'uni',false),'uni',false);

class_overlap_plot = [ ...
    nanmean(vertcat(q_s1{:}),1); ...
    nanmean(vertcat(q_s24{:}),1); ...
    nanmean(vertcat(m_s4{:}),1); ...
    nanmean(vertcat(m_s13{:}),1); ...
    nanmean(vertcat(u_s{:}),1)];

subplot(1,2,2);
bar(class_overlap_plot,'stacked');colormap(gray);
legend({'M/Q','Overlap','Speed'});
set(gca,'XTickLabel',{'Q,S1','Q,S2+','M,S4','M,S3-','U,S'});

%% Classify ROIs by position


for curr_animal = 1:length(data);
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    n_days = length(data(curr_animal).im);
    
    classified_rois(curr_animal).position_up = false(n_rois,n_days);
    classified_rois(curr_animal).position_down = false(n_rois,n_days);
    
    for curr_day = 1:length(data(curr_animal).im)
        
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        lever_position_frames = analysis(curr_animal).lever(curr_day).lever_position_frames;
        
        % Real data
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh > 0);     
        
        roi_position = (curr_act*lever_position_frames)./sum(curr_act,2);
                
        % Activity-epoch shuffled data (movement/quiescent epochs)
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_postion_split = mat2cell(lever_position_frames,diff(boundary_frames));
        
        num_rep = 10000;
        shuffle_perms = shake(repmat([1:length(boundary_frames)-1]',1,num_rep),1);
        lever_position_shuffle = nan(length(lever_position_frames),num_rep);
        for i = 1:num_rep
            lever_position_shuffle(:,i) = ...
                vertcat(lever_postion_split{shuffle_perms(:,i)});
        end
        clear shuffle_perms
        
        roi_position_shuff = bsxfun(@times,curr_act*lever_position_shuffle,1./sum(curr_act,2));
        
        roi_position_rank = tiedrank([roi_position,roi_position_shuff]')';
        roi_position_p = roi_position_rank(:,1)./(num_rep+1);        
        
        % Store
        classified_rois(curr_animal).position_down(:,curr_day) = ...
            roi_position_p > 0.975;
        
        classified_rois(curr_animal).position_up(:,curr_day) = ...
            roi_position_p < 0.025;
        
    end
    disp(curr_animal);
end


%% Quiescent activity across time within session

n_times = 10;

q_act_time = cell(length(data),1);
m_act_time = cell(length(data),1);

for curr_animal = 1:length(data)
    
    q_act_time{curr_animal} = nan(n_times,length(data(curr_animal).im));
    m_act_time{curr_animal} = nan(n_times,length(data(curr_animal).im));
    
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get  activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:) > 0);
        
        % Get lever movement
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);       
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
                
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
                
        split_act = mat2cell(curr_act, ...
            size(curr_act,1),diff(boundary_times));
        
        split_act_q_mean = cell2mat(cellfun(@(x) nanmean(x,1), ...
            split_act(lever_split_q),'uni',false));
        
        split_act_m_mean = cell2mat(cellfun(@(x) nanmean(x,1), ...
            split_act(lever_split_m),'uni',false));
        
        % Split quiescent activity during quiescent epochs into chunks
        time_split_frames_q = linspace(0,length(split_act_q_mean),n_times+1);
        time_split_frames_m = linspace(0,length(split_act_m_mean),n_times+1);
        
        curr_q_act_time = cellfun(@nanmean, ...
            mat2cell(split_act_q_mean',diff(ceil(time_split_frames_q))));       
        curr_m_act_time = cellfun(@nanmean, ...
            mat2cell(split_act_m_mean',diff(ceil(time_split_frames_m))));
        
        curr_q_act_time_norm = curr_q_act_time./curr_q_act_time(1);
        curr_m_act_time_norm = curr_m_act_time./curr_m_act_time(1);
                
        % Store quiescent epoch lever and activity
        q_act_time{curr_animal}(:,curr_day) = curr_q_act_time_norm;
        m_act_time{curr_animal}(:,curr_day) = curr_m_act_time_norm;
        
    end
    disp(curr_animal);
end

q_act_time_cat = horzcat(q_act_time{:});
m_act_time_cat = horzcat(m_act_time{:});

figure; hold on;
errorbar(nanmean(q_act_time_cat,2),nanstd(q_act_time_cat,[],2)./ ...
    sqrt(sum(~isnan(q_act_time_cat),2)),'k');
errorbar(nanmean(m_act_time_cat,2),nanstd(m_act_time_cat,[],2)./ ...
    sqrt(sum(~isnan(m_act_time_cat),2)),'r');

ylabel('Normalized activity')
xlabel('Time bin')
legend({'Quiescence','Movement'})


%% Correlation of M/Q activity during M/Q epochs

q_corr = cell(length(data),1);
m_corr = cell(length(data),1);

for curr_animal = 1:length(data)
    
    q_corr{curr_animal} = nan(length(data(curr_animal).im),1);
    m_corr{curr_animal} = nan(length(data(curr_animal).im),1);
    
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get  activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:) > 0);
        
        % Get lever movement
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);       
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
                
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
                
        split_act = mat2cell(curr_act, ...
            size(curr_act,1),diff(boundary_times));
        
        split_act_q = cell2mat(cellfun(@(x) x(curr_q,:), ...
            split_act(lever_split_q),'uni',false));
        
        split_act_m = cell2mat(cellfun(@(x) x(curr_m,:), ...
            split_act(lever_split_m),'uni',false));
        
        curr_q_corr = nanmean(AP_itril(split_act_q',-1));
        curr_m_corr = nanmean(AP_itril(split_act_m',-1));
               
        % Store quiescent epoch lever and activity
        q_corr{curr_animal}(curr_day) = curr_q_corr;
        m_corr{curr_animal}(curr_day) = curr_m_corr;
        
    end
    disp(curr_animal);
end

q_corr_cat = horzcat(q_corr{:});
m_corr_cat = horzcat(m_corr{:});

figure; hold on;
errorbar(nanmean(q_corr_cat,2),nanstd(q_corr_cat,[],2)./ ...
    sqrt(sum(~isnan(q_corr_cat),2)),'k');
errorbar(nanmean(m_corr_cat,2),nanstd(m_corr_cat,[],2)./ ...
    sqrt(sum(~isnan(m_corr_cat),2)),'r');

ylabel('Correlation within class in corresponding epochs')
xlabel('Day')
legend({'Quiescence','Movement'})


%% Classify ROIs across time within session

classified_rois_split = AP_classify_movement_cells_continuous(data,analysis,4);

q = cellfun(@(x) cat(3,x{:}),{classified_rois_split.quiescent},'uni',false);
q_frac = cellfun(@(x) permute(nanmean(x,1),[2,3,1]),q,'uni',false);

m = cellfun(@(x) cat(3,x{:}),{classified_rois_split.movement},'uni',false);
m_frac = cellfun(@(x) permute(nanmean(x,1),[2,3,1]),m,'uni',false);

% AP152 has 2 days missing, add nans to classified and shift data
AP152 = find(cellfun(@(x) strcmp('AP152',x),{data(:).animal}));

if ~isempty(AP152) && size(q_frac{AP152},2) < 14   
    q_frac{AP152}(:,4:14) = q_frac{AP152}(:,2:12);    
    q_frac{AP152}(:,2:3) = NaN;
    
    m_frac{AP152}(:,4:14) = m_frac{AP152}(:,2:12);    
    m_frac{AP152}(:,2:3) = NaN;
end

q_frac_pad = cellfun(@(x) padarray(x,[0,14-size(x,2)],NaN,'post'),q_frac,'uni',false);
q_frac_padcat = cat(3,q_frac_pad{:});
q_frac_animalmean = nanmean(q_frac_padcat,3);

m_frac_pad = cellfun(@(x) padarray(x,[0,14-size(x,2)],NaN,'post'),m_frac,'uni',false);
m_frac_padcat = cat(3,m_frac_pad{:});
m_frac_animalmean = nanmean(m_frac_padcat,3);

figure; hold on;
set(gca,'ColorOrder', ...
    [cool(size(q_frac_animalmean,1));autumn(size(m_frac_animalmean,1))]);

plot([q_frac_animalmean',m_frac_animalmean']);
xlabel('Day');
ylabel('Fraction classified')
legend({'Q1','Q2','Q3','Q4','M1','M2','M3','M4'});



%% Histogram of Q/M activity during quiescence

n_times = 10;

q_act_q = cell(length(data),1);
q_act_m = cell(length(data),1);

for curr_animal = 1:length(data)
    
    q_act_q{curr_animal} = cell(length(data(curr_animal).im),1);
    q_act_m{curr_animal} = cell(length(data(curr_animal).im),1);
    
    for curr_day = 1:length(data(curr_animal).im);
        
        if isempty(data(curr_animal).im(curr_day).roi_trace_df)
            continue
        end
        
        % Get  activity
        curr_q = logical(classified_rois(curr_animal).quiescent(:,curr_day));
        curr_m = logical(classified_rois(curr_animal).movement(:,curr_day));
        curr_act = double(data(curr_animal).im(curr_day).roi_trace_thresh(:,:) > 0);
        
        % Get lever movement
        lever_active_frames = analysis(curr_animal).lever(curr_day).lever_move_frames;
        
        % Split lever
        boundary_times = find(diff([Inf;lever_active_frames;Inf]) ~= 0);       
        lever_active_split = mat2cell(lever_active_frames,diff(boundary_times));
                
        lever_split_m = cellfun(@(x) any(x),lever_active_split);
        lever_split_q = cellfun(@(x) ~any(x),lever_active_split);
                
        split_act = mat2cell(curr_act, ...
            size(curr_act,1),diff(boundary_times));
        
        act_q_mean = nanmean(cell2mat(split_act(lever_split_q)),2);
        act_m_mean = nanmean(cell2mat(split_act(lever_split_m)),2);
                
        % Store quiescent epoch lever and activity
        q_act_q{curr_animal}{curr_day} = act_q_mean(curr_q);
        q_act_m{curr_animal}{curr_day} = act_q_mean(curr_m);
        
    end
    disp(curr_animal);
end

% this was to look for obvious bimodality in M cells (i.e. are subset of M
% cells high baseline like Q cells?) didn't see anything

%% Fraction movement/quiescent vs position up/down ROIs independently

q = cell2mat(cellfun(@(x) x(:),{classified_rois.quiescent},'uni',false)');
m = cell2mat(cellfun(@(x) x(:),{classified_rois.movement},'uni',false)');
u = cell2mat(cellfun(@(x) x(:),{classified_rois.position_up},'uni',false)');
d = cell2mat(cellfun(@(x) x(:),{classified_rois.position_down},'uni',false)');

qu = q&u;
qd = q&d;
mu = m&u;
md = m&d;

figure;bar(nanmean([qu,qd,mu,md],1));colormap(gray);


%% TO DO: up/down classified? average activity around movement?





