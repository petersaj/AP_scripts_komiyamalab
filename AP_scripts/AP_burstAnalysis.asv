% Burst analysis

% number of frames to do moving window for spike summation
frame_search = 20;
% run through trace and count spikes that occur withing frame_search
num_frames = size(roi_trace_long,2);
frame_scan = frame_search+1:num_frames-(frame_search+1);
num_events = zeros(length(frame_scan),1);
disp('Finding bursts')
tic
parfor i = 1:length(frame_scan);
    curr_frames = frame_scan(i)-frame_search:frame_scan(i)+frame_search;
    for j = 1:length(spike_frames);
        events_in_curr_frames = ismember(spike_frames{j},curr_frames);
        if any(events_in_curr_frames);
            num_events(i) = sum(events_in_curr_frames) + num_events(i);
        end
    end
    disp(['Finding bursts: ' num2str(i) '/' num2str(length(frame_scan))])
end
toc
% find peaks in number of events log
% minimum difference to consider something a peak
min_diff = 25;
[maxtab mintab] = peakdet(num_events,10);
% create cell array for which cells fired during burst
burst_cells = cell(size(maxtab,1),1);
w = waitbar(0,'Finding cells in bursts')
for i = 1:size(maxtab,1);
    curr_frames = maxtab(i,1):maxtab(i,1)+frame_search*2;
    for j = 1:length(spike_frames);
        events_in_curr_frames = ismember(spike_frames{j},curr_frames);
        if any(events_in_curr_frames);
            burst_cells{i} = [burst_cells{i};j];
        end
    end
    waitbar(i/size(maxtab,1),w,'Finding cells in bursts');
end
close(w);

% plot which cells were active during each burst
num_rois = size(roi_trace_long,1)
burst_grid = zeros(num_rois,length(burst_cells));
for i = 1:length(burst_cells)
    burst_grid(burst_cells{i},i) = 1;
end
figure;
imagesc(burst_grid);
colormap(gray);

% correct maxtab for framesearch
maxtab(:,1) = maxtab(:,1) + frame_search;

%% see if bursts correlate with reward


%%
% plot lines for bursts

for i = 1:size(maxtab,1)
    line([maxtab(i,1) maxtab(i,1)],ylim,'linestyle','-','color','r');
end

%% 
% just show cells which participated in burst
roi_trace_burst = zeros(size(roi_trace_long));
for i = 1:length(burst_cells)-1
    roi_trace_burst(1:length(burst_cells{i}),maxtab(i,1)-frame_search:maxtab(i,1)+frame_search) = ...
        roi_trace_long(burst_cells{i},maxtab(i,1)-frame_search:maxtab(i,1)+frame_search);
end
% eliminate space between
zero_column = find(all(roi_trace_burst == 0,1));
roi_trace_burst(:,zero_column) = [];

% order roi_trace_long by which cells participated in each burst
roi_trace_burst = zeros(size(roi_trace_long));
for i = 1:length(burst_cells)-1
    roi_trace_burst(burst_cells{i},maxtab(i,1)-frame_search:maxtab(i,1)+frame_search) = ...
        roi_trace_long(burst_cells{i},maxtab(i,1)-frame_search:maxtab(i,1)+frame_search);
end

% show all, but order for each burst
roi_trace_burst = zeros(size(roi_trace_long));
for i = 1:length(burst_cells)-1
    roi_trace_burst(1:length(burst_cells{i}),maxtab(i,1)-frame_search:maxtab(i+1,1)-1-frame_search) = ...
        roi_trace_long(burst_cells{i},maxtab(i,1)-frame_search:maxtab(i+1,1)-1-frame_search);
    
    non_burst_cells = [1:size(roi_trace_long,1)];;
    non_burst_cells = setdiff(non_burst_cells,burst_cells{i});
    roi_trace_burst(length(burst_cells{i})+1:end,maxtab(i,1)-frame_search:maxtab(i+1,1)-1-frame_search) = ...
        roi_trace_long(non_burst_cells,maxtab(i,1)-frame_search:maxtab(i+1,1)-1-frame_search);
end

%% find cells that bursted together
burst_together = zeros(size(roi_trace_long,1));
for i = 1:size(roi_trace_long,1)
    burst_curr_roi = [];
    burst_curr_roi = burst_grid(i,:) == 1;
    % matrix math yay!
    burst_together(i,:) = (burst_curr_roi*burst_grid')'./sum(burst_curr_roi);
    if sum(burst_curr_roi) == 0;
        burst_together(i,:) = zeros(size(burst_together(i,:)));
    end
end

[temp1 centers] = kmeans(burst_together,3);
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(roi_trace_long));
for i = 1:size(burst_together,1);
    temp_trace(i,:) = roi_trace_long(sort_temp1(i,2),:);
end

% sort cells by similar burst participations
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(burst_grid));
for i = 1:size(burst_grid,1);
    temp_trace(i,:) = burst_grid(sort_temp1(i,2),:);
end



