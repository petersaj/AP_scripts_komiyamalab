%% Load vel threshold, get False/Total positive with threshold

% Set path where manually scored data is dropped
score_path = '/usr/local/lab/People/Andy/Jun_lickreversal/event_detection/vel_threshold/EventScoreDataDrop';
score_files = dir([score_path filesep '*.eventscores']);

% Load in the raw/summed frame match variable
load([score_path filesep 'numframes_match.mat']);

% Load and concatenate all files;
for curr_file = 1:length(score_files)
    temp_scores = load([score_path filesep score_files(curr_file).name],'-mat');
    if curr_file == 1;
        event_scores = temp_scores.event_scores;
    else
        event_scores = [event_scores;temp_scores.event_scores];
    end
end

% Ignore all events with a score of 0 (unjudged)
event_scores([event_scores.score] == 0) = [];

% Don't use events where raw frames/summed frames ~= sum frames
% 1) find animal by mouse number
event_animal = cellfun(@(x) find(cellfun(@(y) strcmp(x,y),{mice.name})), ...
    {event_scores.animal});
% 2) find session by absolute session
% for some reason this is a really quick for-loop, but an
% incredibly slow arrayfun??
event_session = nan(length(event_scores));
for x = 1:length(event_scores)
    event_session(x) = find(event_scores(x).session == ...
        data_all(event_animal(x)).session);
end
% 3) match values with numframes_match, exclude those
frame_match_events = arrayfun(@(x) ...
    numframes_match{event_animal(x)}(event_session(x)), ...
    1:length(event_scores));

% Ignore all events from session with mismatch sum/raw
event_scores(~frame_match_events) = [];

% Get all animals scored
all_animals = unique({event_scores.animal});

% Set thresholds to draw ROC curve from
thresholds = [4:0.1:10];
thresholds_vel = thresholds;

% Get the actual positive and number of score events
score_positive = [event_scores.score] == 1;
score_negative = [event_scores.score] == -1;

% Initialize true/false positives, seperate animals and gad
true_positive = cell(length(all_animals),2);
false_positive = cell(length(all_animals),2);

% Loop through thresholds, get T/F+ rates
for curr_threshold_idx = 1:length(thresholds)
    
    curr_threshold = thresholds(curr_threshold_idx);
    for curr_animal = 1:length(all_animals)
        
        animal = all_animals{curr_animal};
        use_animal = cellfun(@(x) strcmp(x,animal),{event_scores.animal});
        
        for curr_gad = 1:2
            
            % Seperate ex/inh in each animal
            use_gad = [event_scores.gad] == curr_gad-1;
            
            use_events = use_gad & use_animal;
            
            % Get events rated positive by binary threshold
            thresh_positive = [event_scores(:).std_vel] > curr_threshold;
            
            % Get T/F+ rates
            true_positive{curr_animal,curr_gad}(curr_threshold_idx,1:2) = ...
                [sum(thresh_positive(use_events) & ...
                score_positive(use_events)) sum(score_positive(use_events))];
            
            false_positive{curr_animal,curr_gad}(curr_threshold_idx,1:2) = ...
                [sum(thresh_positive(use_events) & ...
                score_negative(use_events)) sum(score_negative(use_events))];
        end
    end
    disp(curr_threshold_idx / length(thresholds));
end

% Get false/total positive rate
ft_ratio_animal = cellfun(@(t,f) f(:,1)./(f(:,1)+t(:,1)), ...
    true_positive,false_positive,'uni',false);

animal_scope = cellfun(@(x) event_scores(find(cellfun(@(y) ...
    strcmp(x,y),{event_scores.animal}),1)).scope,all_animals);

tp_scope{1,1} = sum(cat(3,true_positive{animal_scope == 1,1}),3);
tp_scope{1,2} = sum(cat(3,true_positive{animal_scope == 1,2}),3);
tp_scope{2,1} = sum(cat(3,true_positive{animal_scope == 2,1}),3);
tp_scope{2,2} = sum(cat(3,true_positive{animal_scope == 2,2}),3);

fp_scope{1,1} = sum(cat(3,false_positive{animal_scope == 1,1}),3);
fp_scope{1,2} = sum(cat(3,false_positive{animal_scope == 1,2}),3);
fp_scope{2,1} = sum(cat(3,false_positive{animal_scope == 2,1}),3);
fp_scope{2,2} = sum(cat(3,false_positive{animal_scope == 2,2}),3);

empty_conds = cellfun(@(x) isempty(x),tp_scope);
ft_ratio_scope_vel = cell(size(tp_scope));
ft_ratio_scope_vel(~empty_conds) = cellfun(@(t,f) f(:,1)./(f(:,1)+t(:,1)), ...
    tp_scope(~empty_conds),fp_scope(~empty_conds),'uni',false);



%% Load df threshold, get False/Total positive with threshold

% Set path where manually scored data is dropped
score_path = '/usr/local/lab/People/Andy/Jun_lickreversal/event_detection/df_threshold/EventScoreDataDrop';
score_files = dir([score_path filesep '*.eventscores']);

% Load in the raw/summed frame match variable
load([score_path filesep 'numframes_match.mat']);

% Load and concatenate all files;
for curr_file = 1:length(score_files)
    temp_scores = load([score_path filesep score_files(curr_file).name],'-mat');
    if curr_file == 1;
        event_scores = temp_scores.event_scores;
    else
        event_scores = [event_scores;temp_scores.event_scores];
    end
end

% Ignore all events with a score of 0 (unjudged)
event_scores([event_scores.score] == 0) = [];

% Don't use events where raw frames/summed frames ~= sum frames
% 1) find animal by mouse number
event_animal = cellfun(@(x) find(cellfun(@(y) strcmp(x,y),{mice.name})), ...
    {event_scores.animal});
% 2) find session by absolute session
% for some reason this is a really quick for-loop, but an
% incredibly slow arrayfun??
event_session = nan(length(event_scores));
for x = 1:length(event_scores)
    event_session(x) = find(event_scores(x).session == ...
        data_all(event_animal(x)).session);
end
% 3) match values with numframes_match, exclude those
frame_match_events = arrayfun(@(x) ...
    numframes_match{event_animal(x)}(event_session(x)), ...
    1:length(event_scores));

% Ignore all events from session with mismatch sum/raw
event_scores(~frame_match_events) = [];

% Get all animals scored
all_animals = unique({event_scores.animal});

% Set thresholds to draw ROC curve from
thresholds = [2:0.1:10];
thresholds_df = thresholds;

% Get the actual positive and number of score events
score_positive = [event_scores.score] == 1;
score_negative = [event_scores.score] == -1;

% Initialize true/false positives, seperate animals and gad
true_positive = cell(length(all_animals),2);
false_positive = cell(length(all_animals),2);

% Loop through thresholds, get T/F+ rates
for curr_threshold_idx = 1:length(thresholds)
    
    curr_threshold = thresholds(curr_threshold_idx);
    for curr_animal = 1:length(all_animals)
        
        animal = all_animals{curr_animal};
        use_animal = cellfun(@(x) strcmp(x,animal),{event_scores.animal});
        
        for curr_gad = 1:2
            
            % Seperate ex/inh in each animal
            use_gad = [event_scores.gad] == curr_gad-1;
            
            use_events = use_gad & use_animal;
            
            % Get events rated positive by binary threshold
            thresh_positive = [event_scores(:).std] > curr_threshold;
            
            % Get T/F+ rates
            true_positive{curr_animal,curr_gad}(curr_threshold_idx,1:2) = ...
                [sum(thresh_positive(use_events) & ...
                score_positive(use_events)) sum(score_positive(use_events))];
            
            false_positive{curr_animal,curr_gad}(curr_threshold_idx,1:2) = ...
                [sum(thresh_positive(use_events) & ...
                score_negative(use_events)) sum(score_negative(use_events))];
        end
    end
    disp(curr_threshold_idx / length(thresholds));
end

% Get false/total positive rate

ft_ratio_animal = cellfun(@(t,f) f(:,1)./(f(:,1)+t(:,1)), ...
    true_positive,false_positive,'uni',false);

animal_scope = cellfun(@(x) event_scores(find(cellfun(@(y) ...
    strcmp(x,y),{event_scores.animal}),1)).scope,all_animals);

tp_scope{1,1} = sum(cat(3,true_positive{animal_scope == 1,1}),3);
tp_scope{1,2} = sum(cat(3,true_positive{animal_scope == 1,2}),3);
tp_scope{2,1} = sum(cat(3,true_positive{animal_scope == 2,1}),3);
tp_scope{2,2} = sum(cat(3,true_positive{animal_scope == 2,2}),3);

fp_scope{1,1} = sum(cat(3,false_positive{animal_scope == 1,1}),3);
fp_scope{1,2} = sum(cat(3,false_positive{animal_scope == 1,2}),3);
fp_scope{2,1} = sum(cat(3,false_positive{animal_scope == 2,1}),3);
fp_scope{2,2} = sum(cat(3,false_positive{animal_scope == 2,2}),3);

ft_ratio_scope_df = cellfun(@(t,f) f(:,1)./(f(:,1)+t(:,1)), ...
    tp_scope,fp_scope,'uni',false);


%% Threshold at similar false/total positive rates, get event rate

% Get df/vel values for give false/total positive rates
ft_rates = [0.1:0.05:0.3];
df_true_ft = nan(size(ft_rates));
ft_thresh_df = nan(size(ft_rates));

vel_true_ft = nan(size(ft_rates));
ft_thresh_vel = nan(size(ft_rates));

for curr_ft = 1:length(ft_rates)
   [~,df_idx] = min(abs(ft_ratio_scope_df{2,1} - ft_rates(curr_ft)));
   df_true_ft(curr_ft) = ft_ratio_scope_df{2,1}(df_idx);
   ft_thresh_df(curr_ft) = thresholds_df(df_idx);
   
   [~,vel_idx] = min(abs(ft_ratio_scope_vel{2,1} - ft_rates(curr_ft)));
   vel_true_ft(curr_ft) = ft_ratio_scope_vel{2,1}(vel_idx);
   ft_thresh_vel(curr_ft) = thresholds_vel(vel_idx);
end


% Loop through given animals/sessions, get event rate for f+/t+ rate
bscope_animals = find([mice.scope] == 2);
thresh_event_rate = cell(length(bscope_animals),1);
for curr_animal_idx = 1:length(bscope_animals)
    curr_animal = bscope_animals(curr_animal_idx);
    
    thresh_event_rate{curr_animal_idx} = cell(length(data_all(curr_animal).im),length(ft_rates));
    % Loop through sessions
    for curr_session = 1:length(data_all(curr_animal).im)      
        
        % Loop through thresholds
        for curr_thresh = 1:length(ft_rates)
            
        curr_concat_data = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});
        
        [caEventIdxs_df] = AP_caEvents_threshtest(curr_concat_data,ft_thresh_df(curr_thresh), ...
            ~data_all(curr_animal).labels.gad);
       
        [caEventIdxs_vel] = AP_caEvents_vel_threshtest(curr_concat_data,ft_thresh_vel(curr_thresh), ...
            ~data_all(curr_animal).labels.gad);        
        
        thresh_event_rate{curr_animal_idx}{curr_session,curr_thresh} = ...
            [cellfun(@length,caEventIdxs_df)/size(curr_concat_data,2) ...
            cellfun(@length,caEventIdxs_vel)/size(curr_concat_data,2)];
        
        end
    end
    disp(curr_animal_idx);
end

thresh_event_rate_diff = cellfun(@(x) cellfun(@(x) diff(x,[],2), ...
    x,'uni',false),thresh_event_rate,'uni',false);

for curr_animal = 1:length(thresh_event_rate_diff)
    for curr_el = 1:numel(thresh_event_rate_diff{curr_animal});
       thresh_event_rate_diff{curr_animal}{curr_el}( ...
           ~any(thresh_event_rate{curr_animal}{curr_el},2)) = NaN;
    end
end

thresh_event_rate_diff_avgsession = cellfun(@(x,y) arrayfun(@(z) ...
    nanmean(horzcat(x{:,z}),2),1:size(x,2),'uni',false), ...
    thresh_event_rate_diff,active_cells,'uni',false);

thresh_event_rate_diff_combine = arrayfun(@(x) ...
    cell2mat(cellfun(@(y) y{x}, ...
    thresh_event_rate_diff_avgsession,'uni',false)), ...
    1:length(ft_rates),'uni',false);

% Plot thresholds given false/true positive rates
figure;
subplot(1,2,1);
h = plotyy(ft_ratio_scope_df{2,1},thresholds_df, ...
    ft_ratio_scope_vel{2,1},thresholds_vel);
xlabel('False / Total positive rate')
ylabel(h(1),'\DeltaF');
ylabel(h(2),'Velocity');

% Plot difference in events between each type
subplot(1,2,2);
boxplot(horzcat(thresh_event_rate_diff_combine{:}),'labels', ...
    cellfun(@num2str,num2cell(ft_rates),'uni',false));
xlabel('False / True positive rate')
ylabel('Velocity event rate - \DeltaF event rate')






