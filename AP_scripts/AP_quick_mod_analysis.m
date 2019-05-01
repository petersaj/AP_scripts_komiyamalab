% This is to apply the quick analysis to all cells over days

clear all

animal = 39;
listing_cell = {...
    '120114' ...
    '120115' ...
    '120116' ...
    '120117' ...
    }'

data = struct;
% save figures?
data.save_figures_flag = 1;

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal data
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = [];
    curr_data = load([curr_loop_folder filesep 'Help_ROI.mat']);
    unstruct(curr_data);
    
    % load behavior data
    curr_data = [];
    curr_data = load([curr_loop_folder filesep 'bhv.mat'], ...
        'frames_reward_all');
    unstruct(curr_data);
    
    data.frames_reward_all{trace_loop_file} = frames_reward_all;
    
    % load ROI polygons
    curr_data = [];
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.roi'],'-mat');
    unstruct(curr_data);
    
    data.polygon{trace_loop_file} = polygon;
    
    %% Step 1: Find calcium transients, the dumb way
    peakstd = [3*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
    peak_back = 5;
    local_max = 20;
    drop_percentage = .5;
    spike_frames = {};
    movie_frames = 1000;
    
    num_rois = size(roi_trace_long,1);
    
    % this is to catch a dumb bug: if ROI is offscreen, just make noise
    for i = 1:num_rois
        if all(roi_trace_long(i,:) == roi_trace_long(i,1))
            roi_trace_long(i,:) = rand(size(roi_trace_long(i,:)));
        end
    end
    
    num_rois = size(roi_trace_long,1);
    parfor i = 1:num_rois
        % get estimate of baseline variance
        try
            norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',2); % there are two hidden curves
        catch me
            norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',1); % there are two hidden curves
        end
        baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
        % I think this is std? not sure
        curr_trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve));
        minpeak = curr_trace_std*peakstd;
        spike_frames{i} = ap_schmittDetect(roi_trace_long(i,:),minpeak,peak_back,local_max,drop_percentage);
        disp([num2str(i) '/' num2str(num_rois)])
    end
    % save event times
    data.spike_frames{trace_loop_file} = spike_frames;
    
    %% Step 2: Find heavily modulated cells
    
    % define how many rewards have to have events around them
    percent_reward = 0.2;
    
    mod_cells_early = [];
    mod_cells_middle = [];
    mod_cells_late = [];
    
    mod_cells_early_percent = [];
    mod_cells_middle_percent = [];
    mod_cells_late_percent = [];
    
    % frames to search for reward (early)
    frame_back = 10;
    frame_forward = 10;
    
    surround_reward_all = zeros(size(roi_trace_long,1),length(frames_reward_all));
    
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        for curr_roi = 1:size(roi_trace_long,1)
            event_detect = [];
            event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
            if any(event_detect)
                surround_reward(curr_roi,curr_reward) = 1;
                surround_reward_all(curr_roi,curr_reward) = 1;
            end
        end
    end
    
    sum_rewards_early = sum(surround_reward,2);
    mod_cells_early = find(sum_rewards_early > percent_reward*length(frames_reward_all));
    mod_cells_early_percent = sum_rewards_early(mod_cells_early)./length(frames_reward_all);
    
    % frames to search for reward (middle)
    frame_back = -10;
    frame_forward = 30;
    
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        for curr_roi = 1:size(roi_trace_long,1)
            event_detect = [];
            event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
            if any(event_detect)
                surround_reward(curr_roi,curr_reward) = 1;
                surround_reward_all(curr_roi,curr_reward) = 2;
            end
        end
    end
    
    sum_rewards_middle = sum(surround_reward,2);
    mod_cells_middle = find(sum_rewards_middle > percent_reward*length(frames_reward_all));
    mod_cells_middle_percent = sum_rewards_middle(mod_cells_middle)./length(frames_reward_all);
    
    % frames to search for reward (late)
    frame_back = -30;
    frame_forward = 50;
    
    % events around reward/left up
    surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
    for curr_reward = 1:length(frames_reward_all)
        curr_frame = frames_reward_all(curr_reward);
        curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
        for curr_roi = 1:size(roi_trace_long,1)
            event_detect = [];
            event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
            if any(event_detect)
                surround_reward(curr_roi,curr_reward) = 1;
                surround_reward_all(curr_roi,curr_reward) = 3;
            end
        end
    end
    
    sum_rewards_late = sum(surround_reward,2);
    mod_cells_late = find(sum_rewards_late > percent_reward*length(frames_reward_all));
    mod_cells_late_percent = sum_rewards_late(mod_cells_late)./length(frames_reward_all);
    
    % if there are duplicates, take which category is higher, first for
    % early/middle, then middle/late, then early/late
    overlap_cells = intersect(mod_cells_early,mod_cells_middle);
    for i = overlap_cells';
        cell_indx1 = find(mod_cells_early == i);
        cell_indx2 = find(mod_cells_middle == i);
        check_more = sum_rewards_early(cell_indx1) > sum_rewards_middle(cell_indx2);
        if check_more == 1;
            mod_cells_middle(cell_indx2) = [];
        elseif check_more == 0;
            mod_cells_early(cell_indx1) = [];
        end
    end
    
    overlap_cells = intersect(mod_cells_middle,mod_cells_late);
    for i = overlap_cells';
        cell_indx1 = find(mod_cells_middle == i);
        cell_indx2 = find(mod_cells_late == i);
        check_more = sum_rewards_middle(cell_indx1) > sum_rewards_late(cell_indx2);
        if check_more == 1;
            mod_cells_late(cell_indx2) = [];
        elseif check_more == 0;
            mod_cells_middle(cell_indx1) = [];
        end
    end
    
    
    overlap_cells = intersect(mod_cells_early,mod_cells_late);
    for i = overlap_cells';
        cell_indx1 = find(mod_cells_early == i);
        cell_indx2 = find(mod_cells_late == i);
        check_more = sum_rewards_early(cell_indx1) > sum_rewards_late(cell_indx2);
        if check_more == 1;
            mod_cells_late(cell_indx2) = [];
        elseif check_more == 0;
            mod_cells_early(cell_indx1) = [];
        end
    end
    
    %% Step 3: Find which of those cells shift relative activity time
    %     shift_cells = [];
    %
    %     for curr_mod_cell = mod_up_cells'
    %
    %         event_time = [];
    %         surround_reward_trace = [];
    %
    %         for curr_reward = 1:length(frames_reward_all)
    %             curr_frame = frames_reward_all(curr_reward);
    %             curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    %
    %             % get all of the traces
    %             if any(curr_frame_search > size(roi_trace_long,2));
    %                 curr_frame_search(curr_frame_search > size(roi_trace_long,2)) = size(roi_trace_long,2);
    %             end
    %
    %             % get only the traces with events
    %             event_detect = [];
    %             event_detect = ismember(spike_frames{curr_mod_cell},curr_frame_search);
    %             if any(event_detect)
    %                 % the event start time is when it was within 1 std of baseline
    %                 event_time = [event_time (min(spike_frames{curr_mod_cell}(event_detect)) - curr_frame)];
    %                 surround_reward_trace = [surround_reward_trace;roi_trace_long(curr_mod_cell,curr_frame_search)];
    %             end
    %         end
    %
    %         % check if there was drift, compare first third to last third
    %         third_trials = floor(length(event_time)/3);
    %         h = ttest2(event_time(1:third_trials),event_time(end-third_trials:end));
    %
    %         if h == 1;
    %             shift_cells = [shift_cells curr_mod_cell];
    %         end
    %     end
    
    %% Step 4: Find the correlations between cells
    
    zoom_factor = 2;
    
    x_microns = (850/512)/zoom_factor;
    y_microns = (850/128)/zoom_factor;
    
    % get centers of all ROIs
    roi_centers = zeros(size(roi_trace_long),2);
    roi_polygon_centers = zeros(size(roi_trace_long),2);
    for i = 1:size(roi_trace_long,1);
        % get center x and y from polygon boundaries
        center_x = mean(polygon.ROI{i}(:,1));
        center_y = mean(polygon.ROI{i}(:,2));
        roi_centers(i,1) = center_x*x_microns;
        roi_centers(i,2) = center_y*y_microns;
        roi_polygon_centers(i,1) = center_x;
        roi_polygon_centers(i,2) = -center_y;
    end
    
    % show distances vs. correlation but using coincident spike index instead
    % look by surround
    roi_dists_all = [];
    coincidence_indx_all = [];
    coincidence_indx_tk_all = [];
    cell_pair_all = [];
    num_coincidence_all = [];
    
    for curr_roi = 1:size(roi_trace_long,1);
        
        % for chosen roi, get distances
        curr_center = roi_centers(curr_roi,:);
        center_diffs = roi_centers(:,1) - curr_center(1);
        center_diffs(:,2) = roi_centers(:,2) - curr_center(2);
        roi_dists = sqrt(center_diffs(:,1).^2 + center_diffs(:,2).^2);
        roi_dists_all = [roi_dists_all; roi_dists(:,1)];
        
        frame_back = 1;
        frame_forward = 1;
        roi_trace_coincident = zeros(size(roi_trace_long));
        for i = 1:length(spike_frames{curr_roi})
            frame_search = spike_frames{curr_roi}(i)-frame_back:spike_frames{curr_roi}(i)+frame_forward;
            for j = 1:size(roi_trace_long,1)
                event_detect = [];
                event_detect = ismember(spike_frames{j},frame_search);
                if any(event_detect)
                    roi_trace_coincident(j,spike_frames{curr_roi}(i)) = 1;
                end
            end
        end
        
        % get coincidence index for each ROI (from Komiyama 2010)
        num_frames = size(roi_trace_long,2);
        num_frames_search = size(roi_trace_long,2)/(frame_back+frame_forward+1);
        coincidence_indx = zeros(size(roi_trace_long,1),1);
        coincidence_indx_tk = zeros(size(roi_trace_long,1),1);
        spike_1 = length(spike_frames{curr_roi});
        spike_1_prob = spike_1/num_frames;
        for i = 1:size(roi_trace_long,1)
            spike_2 = length(spike_frames{i});
            spike_2_prob = spike_2/num_frames;
            % this is what they do in the paper... does this make sense?
            geo_mean_spikes = num_frames*sqrt(spike_1_prob*spike_2_prob);
            co_spike = sum(roi_trace_coincident(i,:));
            coincidence_indx_tk(i) = ((co_spike)-(num_frames*spike_1_prob*spike_2_prob))/(geo_mean_spikes);
            % this makes more intuitive sense to me, co-spikes/all spikes
            coincidence_indx(i) = (co_spike*2)/(spike_1+spike_2);
            cell_pair_all = [cell_pair_all; curr_roi i];
            num_coincidence_all = [num_coincidence_all co_spike];
        end
        coincidence_indx_tk_all = [coincidence_indx_tk_all; coincidence_indx_tk];
        coincidence_indx_all = [coincidence_indx_all;coincidence_indx];
        disp(num2str(curr_roi/size(roi_trace_long,1)))
    end
    
    % get rid of cells with themselves (distance = 0)
    self_compare = find(roi_dists_all == 0);
    roi_dists_all(self_compare)= [];
    coincidence_indx_all(self_compare) = [];
    coincidence_indx_tk_all(self_compare) = [];
    cell_pair_all(self_compare,:) = [];
    num_coincidence_all(self_compare) = [];
    
    num_bins = 30;
    [n,bin] = histc(roi_dists_all,linspace(min(roi_dists_all),max(roi_dists_all),num_bins));
    binMean = [];
    for i = 1:num_bins
        flagBinMembers = (bin == i);
        binMembers = coincidence_indx_all(flagBinMembers);
        binMean(i) = nanmean(binMembers);
    end
    binMean_notnan = find(~isnan(binMean) == 1);
    % don't plot at the moment
    %     figure;plot(binMean_notnan,binMean(binMean_notnan),'k','linewidth',2)
    %     set(gca,'XTick',[1:2:num_bins])
    %     set(gca,'XTickLabel',round(linspace(min(roi_dists_all),max(roi_dists_all),num_bins/2)));
    %     xlabel('Distance (\mum)')
    %     ylabel('Correlation coefficient')
    
    % plot the ROIs
    figure
    hold on;
    set(gca,'color','k');
    set(gcf,'color','k');
    for i = 1:length(polygon.ROI)
        r = patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
        %set(r,'linewidth',length(spike_frames{i})/10+0.01);
        set(r,'edgecolor','w');
    end
    
    for i = 1:length(mod_cells_early)
        curr_color = mod_cells_early_percent(i)*0.5+0.5;
        r = patch(polygon.ROI{mod_cells_early(i)}(:,1),-polygon.ROI{mod_cells_early(i)}(:,2),[0 curr_color 0]);
        %set(r,'linewidth',length(spike_frames{mod_cells_early(i)})/10+0.01);
        set(r,'edgecolor','w');
    end
    
    for i = 1:length(mod_cells_middle)
        curr_color = mod_cells_middle_percent(i)*0.5+0.5;
        r = patch(polygon.ROI{mod_cells_middle(i)}(:,1),-polygon.ROI{mod_cells_middle(i)}(:,2),[curr_color curr_color 0.5]);
        %set(r,'linewidth',length(spike_frames{mod_cells_middle(i)})/10+0.01);
        set(r,'edgecolor','w');
    end
    
    for i = 1:length(mod_cells_late)
        curr_color = mod_cells_late_percent(i)*0.5+0.5
        r = patch(polygon.ROI{mod_cells_late(i)}(:,1),-polygon.ROI{mod_cells_late(i)}(:,2),[curr_color 0.5 0.5]);
        %set(r,'linewidth',length(spike_frames{mod_cells_late(i)})/10+0.01);
        set(r,'edgecolor','w');
    end
    
    ylim([-128 0]);
    xlim([0 512]);
    
    % plot lines on ROIs - thicker = more correlated
    for i = 1:length(coincidence_indx_all)
        cell1 = cell_pair_all(i,1);
        cell2 = cell_pair_all(i,2);
        if coincidence_indx_tk_all(i) > 0;
            x = [roi_polygon_centers(cell1,1) roi_polygon_centers(cell2,1)];
            y = [roi_polygon_centers(cell1,2) roi_polygon_centers(cell2,2)];
            z = [0 0];
            r = patch(x,y,z);
            set(r,'edgealpha',coincidence_indx_tk_all(i));
            set(r,'linewidth',num_coincidence_all(i)/10);
            set(r,'edgecolor','g');
        elseif coincidence_indx_tk_all(i) < 0;
            x = [roi_polygon_centers(cell1,1) roi_polygon_centers(cell2,1)];
            y = [roi_polygon_centers(cell1,2) roi_polygon_centers(cell2,2)];
            z = [0 0];
            r = patch(x,y,z);
            set(r,'edgealpha',0.01);%-coincidence_indx_tk_all(i));
            set(r,'linewidth',num_coincidence_all(i)+1/10);
            set(r,'edgecolor','r');
        end
    end
    
    % write numbers in the centers for number of spikes
    for i = 1:length(polygon.ROI);
        curr_num_spike = length(spike_frames{i});
        roi_numbers_h(i) = text(roi_polygon_centers(i,1),roi_polygon_centers(i,2),num2str(curr_num_spike),'color','red');
    end
    
    if data.save_figures_flag == 1;
        curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\'];
        curr_fig_title = [curr_loop_folder 'AP' num2str(animal) '_' num2str(listing_cell{trace_loop_file}) 'test'];
        set(gcf,'PaperPositionMode','auto');
        print('-dpng',curr_fig_title)
    end
    close all;
    
    %% Step 5: Save the results in the structure
    
    data.mod_cells_early{trace_loop_file} = mod_cells_early;
    data.mod_cells_middle{trace_loop_file} = mod_cells_middle;
    data.mod_cells_late{trace_loop_file} = mod_cells_late;
    
    data.mod_cells_early_percent{trace_loop_file} = mod_cells_early_percent;
    data.mod_cells_middle_percent{trace_loop_file} = mod_cells_middle_percent;
    data.mod_cells_late_percent{trace_loop_file} = mod_cells_late_percent;
    
    data.roi_dists_all{trace_loop_file} = roi_dists_all;
    data.coincidence_indx_all{trace_loop_file} = coincidence_indx_all;
    data.coincidence_indx_tk_all{trace_loop_file} = coincidence_indx_tk_all;
    data.cell_pair_all{trace_loop_file} = cell_pair_all;
    data.num_coincidence_all{trace_loop_file} = num_coincidence_all;
    
    disp(['Finished day ' num2str(trace_loop_file) '/' num2str(length(listing_cell))]);
end

%%
figure
hold on
for i = 1:length(data.mod_cells_early)
    plot(i,data.mod_cells_early{i},'.g')
end
for i = 1:length(data.mod_cells_middle)
    if ~isempty(data.mod_cells_middle{i})
        plot(i,data.mod_cells_middle{i},'.b')
    end
end
for i = 1:length(data.mod_cells_late)
    if ~isempty(data.mod_cells_late{i})
        plot(i,data.mod_cells_late{i},'.r')
    end
end

mod_cells_early_num = cellfun(@length,data.mod_cells_early);
mod_cells_middle_num = cellfun(@length,data.mod_cells_middle);
mod_cells_late_num = cellfun(@length,data.mod_cells_late);
figure; hold on;
plot(mod_cells_early_num,'g')
hold on
plot(mod_cells_middle_num,'b')
plot(mod_cells_late_num,'r')
plot(mod_cells_late_num+mod_cells_early_num+mod_cells_middle_num,'k')

