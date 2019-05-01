%% Find cells to plot

%curr_cells = pyr_unfilled(sort_idx);

curr_cells = [31 67 54 149 6 19 71 107 66 128 157 219];
figure
for i = 1:length(curr_cells)
    subplot(1,length(curr_cells),i);
    curr_cell = curr_cells(i);
    imagesc(vertcat(activity_aligned{curr_cell}.reward{:}));
    colormap(gray);
    axis off
    caxis([0 0.5])
    
    day_markers = cumsum(cellfun(@(x) size(x,1),activity_aligned{curr_cell}.reward));
    for j = 1:length(day_markers)-1
        line(xlim,[day_markers(j) day_markers(j)],'color','r')
    end
end

curr_cells = intersect(pyr_unfilled,find(any(vertcat(classified_cells.move_cells_peak{:}))));
figure;hold off
for i = curr_cells
    imagesc(vertcat(activity_aligned{i}.reward{:}));
    colormap(gray);
    title(i)
    k = waitforbuttonpress;
end

%% Plot means (or aligned traces) of activity over days
% Requires set up from behaviorally-aligned activity

curr_cell = 60;

% Reward-aligned 

day_combine = {[1] [2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
mean_traces = zeros(length(day_combine),401);
sem_traces = zeros(length(day_combine),401);
concat_data = cell(length(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(activity_aligned{curr_cell}.reward{day_combine{i}});
    mean_traces(i,:) = nanmean(curr_data);
    sem_traces(i,:) = nanstd(curr_data)/sqrt(size(curr_data,1));
    
    concat_data{i} = curr_data;
end

figure; set(gcf,'Name',['Reward: ' num2str(curr_cell)])
for i = 1:length(day_combine)
   subplot(1,length(day_combine) + 1,i); hold on;
   plot(mean_traces(i,:),'k','linewidth',2,'linesmoothing','off');
   jbfill(1:size(mean_traces,2),mean_traces(i,:)+sem_traces(i,:), ...
       mean_traces(i,:)-sem_traces(i,:),'k','k',0,0.5);
   xlim([0 size(mean_traces,2)+1])
   ylim([0 max(mean_traces(:))]);
   %imagesc(concat_data{i});colormap(gray);caxis([0 0.5]); 
end

subplot(1,length(day_combine) + 1,length(day_combine) + 1);
line([0 3*bhv.framerate],[0 0],'color','k','linewidth',2)
line([0 0],[0 0.1],'color','k','linewidth',2)
xlim([0 size(mean_traces,2)+1]);
ylim([0 max(mean_traces(:))]);

screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[0 0 ...
    screen_size(3) 200])

% Cue-aligned

day_combine = {[1] [2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
mean_traces = zeros(length(day_combine),401);
sem_traces = zeros(length(day_combine),401);
concat_data = cell(length(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(activity_aligned{curr_cell}.cue{day_combine{i}});
    mean_traces(i,:) = nanmean(curr_data);
    sem_traces(i,:) = nanstd(curr_data)/sqrt(size(curr_data,1));
    
    concat_data{i} = curr_data;
end

figure; set(gcf,'Name',['Cue: ' num2str(curr_cell)])
for i = 1:length(day_combine)
   subplot(1,length(day_combine) + 1,i); hold on;
   plot(mean_traces(i,:),'k','linewidth',2,'linesmoothing','off');
   jbfill(1:size(mean_traces,2),mean_traces(i,:)+sem_traces(i,:), ...
       mean_traces(i,:)-sem_traces(i,:),'k','k',0,0.5);
   xlim([0 size(mean_traces,2)+1])
   ylim([0 max(mean_traces(:))]);
   %imagesc(concat_data{i});colormap(gray);caxis([0 0.5]);
  
end

subplot(1,length(day_combine) + 1,length(day_combine) + 1);
line([0 3*bhv.framerate],[0 0],'color','k','linewidth',2)
line([0 0],[0 0.1],'color','k','linewidth',2)
xlim([0 size(mean_traces,2)+1]);
ylim([0 max(mean_traces(:))]);

screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[0 300 ...
    screen_size(3) 200])

% Cell images

day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
mean_images = nan(41,41,length(day_combine));
for i = 1:length(day_combine)
    curr_data = cat(3,roi_image{curr_cell}{day_combine{i}});
    mean_images(:,:,i) = nanmean(curr_data,3);
end

figure; set(gcf,'Name',[num2str(curr_cell)])
for i = 1:length(day_combine)
    subplot(1,length(day_combine) + 1,i); hold on;
    imagesc(mean_images(:,:,i));
    colormap(gray);
    ylim([1 41]);
    xlim([1 41]);
    axis off
end

subplot(1,length(day_combine) + 1,length(day_combine) + 1);
line([0 0],[0 10],'color','k','linewidth',2)
ylim([1 41]);
xlim([1 41]);

screen_size = get(0,'ScreenSize');
fig_size = get(gcf,'Position');
set(gcf,'Position',[0 600 ...
    screen_size(3) 200])











