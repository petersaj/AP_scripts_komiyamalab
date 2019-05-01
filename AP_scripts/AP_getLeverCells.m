%% Categorize cells as responsive to lever activity or not

% get times around lever presses and not around lever presses

frame_surround = 5;
n_sample = 100;


push_time = zeros(1,size(roi_trace_long,2));

frames_all = [frames_left_down_all; frames_left_up_all];
for i = 1:length(frames_left_down_all)
    lever_time = [frames_all(i) - 5:frames_all(i) + 5];
    lever_time(lever_time < 1 | lever_time > size(roi_trace_long,2)) = [];
    push_time(lever_time) = 1;
end

% get random samples and actual counts for events during lever activity

push_time_on = find(push_time == 1);
rand_push_time = zeros(n_sample, size(roi_trace_long,2));
rand_spike_count = zeros(size(roi_trace_long,1),n_sample);

% randomly shuffle lever related times

break_points = find(diff(push_time) == -1)+1;
% define first and last segments
push_segs{1} = push_time(1:break_points(1));
push_segs{length(break_points)+1} = push_time(break_points(end)+1:end);
% define all middle segments
for i = 2:length(break_points);
    push_segs{i} = push_time(break_points(i-1)+1:break_points(i));
end

w = waitbar(0,'Shuffling activity');
for j = 1:n_sample;
    % divide lever active times up and shuffle them
    shuffle_segs = randperm(length(push_segs));
    temp_rand_push_time = [];
    for i = 1:length(shuffle_segs)
        temp_rand_push_time = [temp_rand_push_time push_segs{shuffle_segs(i)}];
    end
    rand_push_time(j,:) = temp_rand_push_time;
    waitbar(j/n_sample, w, 'Shuffling activity');
end
close(w);

w = waitbar(0,'Calculating shuffled events');
for j = 1:n_sample;
    rand_push_time_on = find(rand_push_time(j,:) == 1);
    for i = 1:size(roi_trace_long,1);
        rand_spike_count(i,j) = sum(ismember(spike_frames{i}, rand_push_time_on));
    end
    waitbar(j/n_sample, w, 'Calculating shuffled events');
end
close(w);

lever_spike_count = zeros(size(roi_trace_long,1),1);
for i = 1:size(roi_trace_long,1);
    lever_spike_count(i) = sum(ismember(spike_frames{i}, push_time_on));
end

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [rand_spike_count(curr_roi,:) lever_spike_count(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end  

lever_up_cells = find(p > 0.95);
lever_down_cells = find(p < 0.05);


% show up/down/non-modulated cells

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
end

for i = 1:length(lever_up_cells)
    patch(polygon.ROI{lever_up_cells(i)}(:,1),-polygon.ROI{lever_up_cells(i)}(:,2),[0 1 0]);
end

for i = 1:length(lever_down_cells)
    patch(polygon.ROI{lever_down_cells(i)}(:,1),-polygon.ROI{lever_down_cells(i)}(:,2),[1 0 0]);
end

title('Cells with events correlated with lever activity')

% plot modulated cells during and not during lever activity
figure;imagesc(roi_trace_long(lever_up_cells,push_time == 1));
title('Activity of modulated cells during active lever time');
figure;imagesc(roi_trace_long(lever_up_cells,push_time == 0));
title('Activity of modulated cells during non-active lever time')

%% Get correlation between trace and lever activity

frame_surround = 5;
n_sample = 100;


push_time = zeros(1,size(roi_trace_long,2));

frames_all = [frames_left_down_all; frames_left_up_all];
for i = 1:length(frames_left_down_all)
    lever_time = [frames_all(i) - 5:frames_all(i) + 5];
    lever_time(lever_time < 1 | lever_time > size(roi_trace_long,2)) = [];
    push_time(lever_time) = 1;
end

push_time_on = find(push_time == 1);
rand_push_time = zeros(n_sample, size(roi_trace_long,2));
rand_spike_count = zeros(size(roi_trace_long,1),n_sample);

% randomly shuffle lever related times

break_points = find(diff(push_time) == -1)+1;
% define first and last segments
push_segs{1} = push_time(1:break_points(1));
push_segs{length(break_points)+1} = push_time(break_points(end)+1:end);
% define all middle segments
for i = 2:length(break_points);
    push_segs{i} = push_time(break_points(i-1)+1:break_points(i));
end

w = waitbar(0,'Shuffling activity');
for j = 1:n_sample;
    % divide lever active times up and shuffle them
    shuffle_segs = randperm(length(push_segs));
    temp_rand_push_time = [];
    for i = 1:length(shuffle_segs)
        temp_rand_push_time = [temp_rand_push_time push_segs{shuffle_segs(i)}];
    end
    rand_push_time(j,:) = temp_rand_push_time;
    waitbar(j/n_sample, w, 'Shuffling activity');
end
close(w);

r_rand = zeros(size(roi_trace_long,1),n_sample);
w = waitbar(0,'Calculating random correlation');
for j = 1:n_sample
    for i = 1:size(roi_trace_long,1);
        curr_corrcoef = corrcoef(roi_trace_long(i,:),rand_push_time(j,:));
        r_rand(i,j) = curr_corrcoef(2);
    end
    waitbar(j/n_sample,w,'Calculating random correlation');
end
close(w);

r = zeros(size(roi_trace_long,1),1);
for i = 1:size(roi_trace_long,1);
    curr_corrcoef = corrcoef(roi_trace_long(i,:),push_time);
    r(i) = curr_corrcoef(2);
end

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [r_rand(curr_roi,:) r(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end  

lever_cor_up_cells = find(p > 0.95);
lever_cor_down_cells = find(p < 0.05);

% show up/down/non-modulated cells

p(isnan(p)) = 0.5; % temp fix, if NaN for p, set at 0.5
figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
end

for i = 1:length(lever_cor_up_cells)
    patch(polygon.ROI{lever_cor_up_cells(i)}(:,1),-polygon.ROI{lever_cor_up_cells(i)}(:,2),[0 1 0]);
end

for i = 1:length(lever_cor_down_cells)
    patch(polygon.ROI{lever_cor_down_cells(i)}(:,1),-polygon.ROI{lever_cor_down_cells(i)}(:,2),[1 0 0]);
end

title('Cells with traces correlated with lever activity times')


p(isnan(p)) = 0.5; % temp fix, if NaN for p, set at 0.5
figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1-p(i) p(i) 0]);
end
title('Cells with traces correlated with lever activity times, graded')

