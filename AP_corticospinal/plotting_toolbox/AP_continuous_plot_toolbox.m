%% This requires previous loading and processing

%% Plot individual cell with events marked
curr_animal = 2;
curr_session = 13;
curr_cell = 69;

% Plot df/f and resampled lever
figure;hold on;
plot([1:size(data(curr_animal).im(curr_session).roi_trace_df,2)],data(curr_animal).im(curr_session).roi_trace_df(curr_cell,:),'k');
plot([1:size(data(curr_animal).im(curr_session).roi_trace_thresh,2)],data(curr_animal).im(curr_session).roi_trace_thresh(curr_cell,:),'b');
plot([1:length(data(curr_animal).bhv(curr_session).lever_force_plot)],data(curr_animal).bhv(curr_session).lever_force_plot*3-6,'r');

% Get paradigm
catch_paradigm = isfield(data(curr_animal).bhv,'catch_trials');
reward_delay_paradigm = any(cellfun(@(x) isfield(x.states,'reward_delay'),data(curr_animal).bhv(curr_session).bhv_frames));

% Get frames of events
cue_frames = cellfun(@(x) x.states.cue(1),data(curr_animal).bhv(curr_session).bhv_frames);

rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),data(curr_animal).bhv(curr_session).bhv_frames);
reward_frames = cellfun(@(x) x.states.reward(1),data(curr_animal).bhv(curr_session).bhv_frames(rewarded_trials));

reward_delay_trials = cellfun(@(x) isfield(x.states,'reward_delay') && ...
    ~isempty(x.states.reward_delay),data(curr_animal).bhv(curr_session).bhv_frames);
reward_delay_frames = cellfun(@(x) x.states.reward(1),data(curr_animal).bhv(curr_session).bhv_frames(reward_delay_trials));

lick_frames = cell2mat(cellfun(@(x) x.pokes.C(:,1),data(curr_animal).bhv(curr_session).bhv_frames,'uni',false));
 
% plot cues as magenta lines
for i = 1:length(cue_frames)
    temph = line([cue_frames(i) cue_frames(i)],ylim,'color','m');
    if data(curr_animal).bhv(curr_session).catch_trials(i)
       set(temph,'color','b') 
    end
end

% plot rewards as green lines
for i = 1:length(reward_frames)
    line([reward_frames(i) reward_frames(i)],ylim,'color','g');
end

% plot reward delays as dotted green lines
for i = 1:length(reward_delay_frames)
    line([reward_delay_frames(i) reward_delay_frames(i)],ylim,'color','g', ...
        'linestyle','--');
end

% plot licks as black dots
plot(lick_frames,2,'.r');


%% Plot aligned move onset/offset all days for one ROI
curr_animal =3;
curr_roi = 195;

curr_onset = arrayfun(@(x) analysis(curr_animal).im(x).move_onset_aligned(:,:,curr_roi),1:14,'uni',false);
curr_onset_trials = cellfun(@(x) size(x,1),curr_onset);
curr_onset_trials_sum = cumsum(curr_onset_trials);

curr_offset = arrayfun(@(x) analysis(curr_animal).im(x).move_offset_aligned(:,:,curr_roi),1:14,'uni',false);
curr_offset_trials = cellfun(@(x) size(x,1),curr_offset);
curr_offset_trials_sum = cumsum(curr_offset_trials);

figure;

subplot(1,2,1);
imagesc(vertcat(curr_onset{:}));colormap(hot);
for i = 1:13;
    line(xlim,repmat(curr_onset_trials_sum(i),1,2),'color','w')
end

subplot(1,2,2);
imagesc(vertcat(curr_offset{:}));colormap(hot);
for i = 1:13;
    line(xlim,repmat(curr_offset_trials_sum(i),1,2),'color','w')
end



%% Mean/sort day activity (aligned)

curr_animal = 4;

avg_act_onset = nan(size(analysis(curr_animal).im(1).move_onset_aligned,3), ...
    size(analysis(curr_animal).im(1).move_onset_aligned,2),14);
for i = 1:14
    avg_act_onset(:,:,i) = permute(nanmean(analysis(curr_animal).im(i).move_onset_aligned,1),[3 2 1]); 
end

[~,max_idx] = max(nanmean(avg_act_onset,3),[],2);
[~,sort_idx] = sort(max_idx);

AP_image_scroll(avg_act_onset(sort_idx,:,:));colormap(hot)



avg_act_offset = nan(size(analysis(curr_animal).im(1).move_offset_aligned,3), ...
    size(analysis(curr_animal).im(1).move_offset_aligned,2),14);
for i = 1:14
    avg_act_offset(:,:,i) = permute(nanmean(analysis(curr_animal).im(i).move_offset_aligned,1),[3 2 1]); 
end

[~,max_idx] = max(nanmean(avg_act_offset,3),[],2);
[~,sort_idx] = sort(max_idx);

AP_image_scroll(avg_act_offset(sort_idx,:,:));colormap(hot)


%% Plot traces of all given classified ROIs for given day

curr_animal = 1;
curr_session = 3;

curr_class = find(classified_rois(curr_animal).movement(:,curr_session));

figure; hold on
for i = 1:length(curr_class)
    plot(data(curr_animal).im(curr_session).roi_trace_df(curr_class(i),:) + i*5,'k');
end



































