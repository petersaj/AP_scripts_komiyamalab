%% Load data
clear all;

animal = 'AP107';

%sessions = {'140805' '140806' '140807','140808','140809','140810','140811'};
sessions = {'140814'};

% change data path structure depending on computer
switch computer
    case 'GLNXA64'
        data_path = '/usr/local/lab/People/Andy/Data/';
    case 'PCWIN'
        data_path = 'Z:\People\Andy\Data\';
end

im = struct;
bhv = struct;
for curr_session = 1:length(sessions)
    session = sessions{curr_session};
    
    im_traces = load([data_path animal filesep 'manual_rois' filesep session '_' animal '_' 'analysis.mat']);
    im(curr_session).roi_trace_long = im_traces.roi_trace_long;
    im(curr_session).roi_trace_df_g = AP_baselineEstimation(im(curr_session).roi_trace_long(:,1:2:end),30);
    im(curr_session).roi_trace_df_r = AP_baselineEstimation(im(curr_session).roi_trace_long(:,2:2:end),30);
    
    xsg_data = AP_read_ephus_dir_raw([data_path animal filesep session filesep 'AP0' animal(3:end)]);
    warning off;
    bhv_filename = dir([data_path animal filesep session filesep 'data_@*']);
    dispatcher_complete = load([data_path animal filesep session filesep bhv_filename.name]);
    warning on;
    bhv_data_all = dispatcher_complete.saved_history.ProtocolsSection_parsed_events;
    
    xsg_samplerate = 10000;
    xsg_frame = xsg_data.channels(:,1);
    xsg_lever = xsg_data.channels(:,2);
    xsg_trial_number = xsg_data.channels(:,3);
    
    all_trial_times = AP_ReadBitCode(xsg_trial_number);
    % convert time s to ms
    all_trial_times(:,1) = all_trial_times(:,1)*1000;
    
    % Match bhv trace to times (only use complete + imaged trials)
    % 1) get rid of trials in xsg that are not in bhv (incomplete trial)
    bhv(curr_session).trial_times = all_trial_times(all_trial_times(:,2) <= length(bhv_data_all),:);
    % 2) get rid of bhv that aren't in xsg (occasionally if false start)
    bhv(curr_session).dispatcher = bhv_data_all(bhv(curr_session).trial_times(:,2));
    
    frames = find(xsg_frame(1:end-1) < 2.5 & xsg_frame(2:end) > 2.5);
    % the end always has a few extra triggers because of manual abort, cut off
    frames = frames(1:size(im(curr_session).roi_trace_df_g,2));
    im(curr_session).frames_ms = frames/10;
    
    [bhv(curr_session).lever_active, bhv(curr_session).lever_force_smooth] = ...
        AP_parseLeverMovement_continuous(xsg_lever);
    
    disp(['Loaded session ' session])
end


%% Get aligned activity (in ms)

for curr_session = 1:length(im)
    
    bhv(curr_session).lick_times = cell2mat(cellfun(@(x,trial_start) (x.pokes.C(:,1) - ...
        x.states.bitcode(1))*1000+trial_start,bhv(curr_session).dispatcher,num2cell(bhv(curr_session).trial_times(:,1)),'uni',false));
    % (too many licks, have to do this differently)
    [offset, lick_frames_idx] = arrayfun(@(x) min(abs(bhv(curr_session).lick_times(x) - im(curr_session).frames_ms)),1:length(bhv(curr_session).lick_times));
    
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),bhv(curr_session).dispatcher);
    bhv(curr_session).rewarded_trials = rewarded_trials;

    cue_timediff = cellfun(@(x) x.states.cue(1) - x.states.bitcode(1),bhv(curr_session).dispatcher(rewarded_trials))*1000;
    reward_timediff = cellfun(@(x) x.states.reward(1) - x.states.bitcode(1),bhv(curr_session).dispatcher(rewarded_trials))*1000;
    punish_timediff = cellfun(@(x) x.states.punish(1) - x.states.bitcode(1),bhv(curr_session).dispatcher(~rewarded_trials))*1000;
    
    [offset, cue_frames_idx] = min(abs(bsxfun(@minus,bhv(curr_session).trial_times(rewarded_trials,1) + cue_timediff,im(curr_session).frames_ms')),[],2);
    [offset, reward_frames_idx] = min(abs(bsxfun(@minus,bhv(curr_session).trial_times(rewarded_trials,1) + reward_timediff,im(curr_session).frames_ms')),[],2);
    [offset, punish_frames_idx] = min(abs(bsxfun(@minus,bhv(curr_session).trial_times(~rewarded_trials,1) + punish_timediff,im(curr_session).frames_ms')),[],2);
    [offset, move_onset_frames_idx] = min(abs(bsxfun(@minus,find(diff([0;bhv(curr_session).lever_active;0]) == 1),im(curr_session).frames_ms')),[],2);
    [offset, move_offset_frames_idx] = min(abs(bsxfun(@minus,find(diff([0;bhv(curr_session).lever_active;0]) == -1),im(curr_session).frames_ms')),[],2);
    
    surround_frames = [30 90];
    % don't use events without full surround before and after
    cue_frames_idx_use = cue_frames_idx-surround_frames(1) > 0 & ...
        cue_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
    im(curr_session).cue_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
        num2cell(cue_frames_idx(cue_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);
    
    reward_frames_idx_use = reward_frames_idx-surround_frames(1) > 0 & ...
        reward_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
    im(curr_session).reward_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
        num2cell(reward_frames_idx(reward_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);
    
    punish_frames_idx_use = punish_frames_idx-surround_frames(1) > 0 & ...
        punish_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
    im(curr_session).punish_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
        num2cell(punish_frames_idx(punish_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);
    
    move_onset_frames_idx_use = move_onset_frames_idx-surround_frames(1) > 0 & ...
        move_onset_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
    im(curr_session).move_onset_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
        num2cell(move_onset_frames_idx(move_onset_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);
    
    % to easily pick out the cued-rewarded movements later
    bhv(curr_session).move_onset_use = move_onset_frames_idx_use;
    
    move_offset_frames_idx_use = move_offset_frames_idx-surround_frames(1) > 0 & ...
        move_offset_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
    im(curr_session).move_offset_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
        num2cell(move_offset_frames_idx(move_offset_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);  
    
    % if reward delay paradigm, grab at the press/reward sound too
    if isfield(bhv(1).dispatcher{1}.states,'reward_delay');
        press_timediff = cellfun(@(x) x.states.reward_delay(1) - x.states.bitcode(1),bhv(curr_session).dispatcher(rewarded_trials))*1000;
        [offset, press_frames_idx] = min(abs(bsxfun(@minus,bhv(curr_session).trial_times(rewarded_trials,1) + press_timediff,im(curr_session).frames_ms')),[],2);
        press_frames_idx_use = press_frames_idx-surround_frames(1) > 0 & ...
            press_frames_idx + surround_frames(2) < size(im(curr_session).roi_trace_df_g,2);
        im(curr_session).press_aligned_df = permute(cell2mat(permute(cellfun(@(x) im(curr_session).roi_trace_df_g(:,x-surround_frames(1):x+surround_frames(2)), ...
            num2cell(press_frames_idx(press_frames_idx_use)),'uni',false),[2 3 1])),[3 2 1]);
    end
    
    
    disp(['Aligned session ' num2str(curr_session)])
end


%% Plot single cell traces
session = 1;
roi = 2;

a = smooth(im(session).roi_trace_df_g(roi,:),30,'loess');

figure;hold on;
plot(im(session).frames_ms,a,'k')
plot(bhv(session).lever_force_smooth+1,'r')
plot(bhv(session).lick_times,1.7,'ob','MarkerSize',5)


%% Get average lever position around activity

curr_session = 7;

% do this by averaging the lever values within that frame time

lever_loess = smooth(bhv(curr_session).lever_force_smooth,30,'loess');
df_loess = zeros(size(im(curr_session).roi_trace_df_g));
for i = 1:size(df_loess,1)
   df_loess(i,:) = smooth(im(curr_session).roi_trace_df_g(i,:),30,'loess'); 
end

avg_lever_frame = arrayfun(@(x) nanmean(lever_loess(round( ...
    im(curr_session).frames_ms(x-1)):round(im(curr_session).frames_ms(x)))),2:size(df_loess,2));

% cut out the extra imaged frame (the first one)
df_loess = df_loess(:,2:end);

% dirty events: get by gad method from before
caEvents = AP_caEvents(im(curr_session).roi_trace_df_g(:,2:end),[],1:size(im(curr_session).roi_trace_df_g,1));

activity_starts = arrayfun(@(x) find(diff(caEvents(x,:) ~= 0,[],2) > 0),1:size(caEvents,1),'uni',false);
lever_surround = 30;
activity_lever = cell(size(df_loess,1),1);
for roi = 1:size(df_loess,1)
    use_activity_starts = activity_starts{roi} - lever_surround > 0 & ...
        activity_starts{roi} + lever_surround < length(avg_lever_frame);
    activity_lever{roi} = cell2mat(cellfun(@(x) avg_lever_frame(x-lever_surround:x+ ...
        lever_surround),num2cell(activity_starts{roi}(use_activity_starts))','uni',false));
end

square_plot = ceil(sqrt(length(activity_lever)));
figure;
for roi = 1:length(activity_lever)
   subplot(square_plot,square_plot,roi);
   
   curr_mean = nanmean(activity_lever{roi},1);
   curr_sem = nanstd(activity_lever{roi},[],1)./sqrt(sum(~isnan(activity_lever{roi}),1));
   plot(curr_mean,'k','linewidth',2)
   jbfill(1:length(curr_mean),curr_mean+curr_sem,curr_mean-curr_sem,[0.5 0.5 0.5]);
   line(repmat(surround_frames(1)+1,2,1),ylim,'color','r','linestyle','--','linewidth',2)
end



%% Plot average traces around movements

% pull out cued-rewarded movements

figure;
cued_rewarded_move_df_mean_all = cell(7,1);
for curr_session = 1:7
    subplot(1,7,curr_session)
    
    cue_timediff = cellfun(@(x) x.states.cue(1) - x.states.bitcode(1), ...
        bhv(curr_session).dispatcher(bhv(curr_session).rewarded_trials))*1000;
    reward_timediff = cellfun(@(x) x.states.reward(1) - x.states.bitcode(1), ...
        bhv(curr_session).dispatcher(bhv(curr_session).rewarded_trials))*1000;
    
    cue_times = round(bhv(curr_session).trial_times(bhv(curr_session).rewarded_trials,1) + cue_timediff);
    reward_times = round(bhv(curr_session).trial_times(bhv(curr_session).rewarded_trials,1) + reward_timediff);
    
    move_lengths = find(diff([0;bhv(curr_session).lever_active;0]) == -1) - find(diff([0;bhv(curr_session).lever_active;0]) == 1);
    move_times = mat2cell(find(bhv(curr_session).lever_active),move_lengths);
    
    cued_rewarded_move = cellfun(@(x) any(intersect(x,reward_times)) & ~any(intersect(x,cue_times)),move_times);
    % only use the ones which were pulled out of df earlier
    cued_rewarded_move_use = cued_rewarded_move(bhv(curr_session).move_onset_use);
    
    cued_rewarded_move_df_mean = permute(nanmean(im(curr_session).move_onset_aligned_df(cued_rewarded_move_use,:,:),1),[3 2 1]);
    
    [asdf max_idx] = max(cued_rewarded_move_df_mean,[],2);
    [asdf sort_idx] = sort(max_idx);
    imagesc(cued_rewarded_move_df_mean(sort_idx,:));colormap(gray);
    line(repmat(surround_frames(1)+1,2,1),ylim,'color','r','linestyle','--','linewidth',2)
    axis off
    
    cued_rewarded_move_df_mean_all{curr_session} = cued_rewarded_move_df_mean;
end









