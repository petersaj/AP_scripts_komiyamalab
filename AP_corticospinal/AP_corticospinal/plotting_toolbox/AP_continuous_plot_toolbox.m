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
curr_animal = 4;
curr_roi = 21;

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


%% Heatmap of activity of all movements (on/off stitch) all cells sorted

% Get activity with all movements > n sec

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

  
        move_act_onset_all{curr_animal,curr_day} = move_act_onset_timed;
        move_act_offset_all{curr_animal,curr_day} = move_act_offset_timed;

        
    end
    disp(curr_animal);
end

move_act_onset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_onset_all,'uni',false);
move_act_offset_all_mean = cellfun(@(x) permute(nanmean(x,1),[3,2,1]),move_act_offset_all,'uni',false);

for curr_animal = 1:length(data)   
    figure;
    for curr_day = 1:length(data(curr_animal).im)        
        
        curr_data = [move_act_onset_all_mean{curr_animal,curr_day}, ...
            move_act_offset_all_mean{curr_animal,curr_day}];
       
        m = classified_rois(curr_animal).movement(:,curr_day) == 1;
        q = classified_rois(curr_animal).quiescent(:,curr_day) == 1;
        u = ~m & ~q;
        
        curr_data_m = curr_data(m,:);
        curr_data_q = curr_data(q,:);
        curr_data_u = curr_data(u,:);
        
        % Sort by difference between m/q activity
        m_frames = false(size(curr_data,2),1);
        m_frames(nonmove_frames+1:nonmove_frames+1+move_frames*2) = true;
        
        curr_data_m_diff = ...
            (nanmean(curr_data_m(:,m_frames),2)-nanmean(curr_data_m(:,~m_frames),2))./ ...
            (nanmean(curr_data_m(:,m_frames),2)+nanmean(curr_data_m(:,~m_frames),2));
        
        curr_data_q_diff = ...
            (nanmean(curr_data_q(:,m_frames),2)-nanmean(curr_data_q(:,~m_frames),2))./ ...
            (nanmean(curr_data_q(:,m_frames),2)+nanmean(curr_data_q(:,~m_frames),2));
        
        curr_data_u_diff = ...
            (nanmean(curr_data_u(:,m_frames),2)-nanmean(curr_data_u(:,~m_frames),2))./ ...
            (nanmean(curr_data_u(:,m_frames),2)+nanmean(curr_data_u(:,~m_frames),2));       
        
        [~,sort_idx_m] = sort(curr_data_m_diff);
        [~,sort_idx_q] = sort(curr_data_q_diff);
        [~,sort_idx_u] = sort(curr_data_u_diff);
        
        curr_data_sort = [curr_data_m(sort_idx_m,:); ...
            curr_data_u(sort_idx_u,:); ...
            curr_data_q(sort_idx_q,:)];
        
        curr_data_sort_minsub = bsxfun(@minus,curr_data_sort,min(curr_data_sort,[],2));
        curr_data_sort_norm = bsxfun(@times,curr_data_sort_minsub,1./max(curr_data_sort_minsub,[],2));
        
        subplot(2,7,curr_day);
        imagesc(curr_data_sort_norm);colormap(gray);
        
        num_rois = cumsum([sum(m),sum(u)]);
        line(xlim,[sum(m),sum(m)],'color','r');
        line(xlim,[sum(m)+sum(u),sum(m)+sum(u)],'color','r');
        
        colormap(gray);
        
    end    
end





%% Plot traces of all given classified ROIs for given day

curr_animal = 1;
curr_session = 3;

curr_class = find(classified_rois(curr_animal).movement(:,curr_session));

figure; hold on
for i = 1:length(curr_class)
    plot(data(curr_animal).im(curr_session).roi_trace_df(curr_class(i),:) + i*5,'k');
end



































