%% Load raw data

animal = 'AP123';
curr_session = 5;

analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
analysis_files = dir([analysis_path filesep '*.mat']);
num_sessions = length(analysis_files);

curr_file = [analysis_path filesep analysis_files(curr_session).name];
curr_data = load(curr_file);


%% Load and align data

analysis = AP_corticospinal_prepare_processed(animal,sessions);

%% Plot individual cell with events marked

curr_cell = 108;

% Plot df/f and resampled lever
figure;hold on;
plot([1:size(curr_data.im.roi_trace_df,2)],curr_data.im.roi_trace_df(curr_cell,:),'k');
plot([1:length(curr_data.bhv.lever_force_plot)],curr_data.bhv.lever_force_plot,'r');

% Get paradigm
catch_paradigm = isfield(curr_data.bhv,'catch_trials');
reward_delay_paradigm = any(cellfun(@(x) isfield(x.states,'reward_delay'),curr_data.bhv.bhv_frames));

% Get frames of events
cue_frames = cellfun(@(x) x.states.cue(1),curr_data.bhv.bhv_frames);

rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),curr_data.bhv.bhv_frames);
reward_frames = cellfun(@(x) x.states.reward(1),curr_data.bhv.bhv_frames(rewarded_trials));

reward_delay_trials = cellfun(@(x) isfield(x.states,'reward_delay') && ...
    ~isempty(x.states.reward_delay),curr_data.bhv.bhv_frames);
reward_delay_frames = cellfun(@(x) x.states.reward(1),curr_data.bhv.bhv_frames(reward_delay_trials));

lick_frames = cell2mat(cellfun(@(x) x.pokes.C(:,1),curr_data.bhv.bhv_frames,'uni',false));
 
% plot cues as magenta lines
for i = 1:length(cue_frames)
    temph = line([cue_frames(i) cue_frames(i)],ylim,'color','m');
    if curr_data.bhv.catch_trials(i)
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











