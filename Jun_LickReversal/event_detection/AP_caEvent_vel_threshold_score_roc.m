% NEW! REQUIRES DATA LOADED IN (SAME CONFIGURATION AS WHEN EVENTS SCORED)
% This is to check for frame mismatch between raw and sum

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

%% Plot results

% Get true/false positive rates for animals / scopes
tp_rate_animal = cellfun(@(x) x(:,1)./x(:,2),true_positive,'uni',false);
fp_rate_animal = cellfun(@(x) x(:,1)./x(:,2),false_positive,'uni',false);

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
tp_rate_scope = cell(size(tp_scope));
tp_rate_scope(~empty_conds) = cellfun(@(x) x(:,1)./x(:,2),tp_scope(~empty_conds),'uni',false);
fp_rate_scope = cell(size(fp_scope));
fp_rate_scope(~empty_conds) = cellfun(@(x) x(:,1)./x(:,2),fp_scope(~empty_conds),'uni',false);

% Plot ROC / False positive by threshold

figure;

subplot(1,2,1); hold on;
for curr_animal = 1:length(all_animals)
    switch animal_scope(curr_animal)
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(fp_rate_animal{curr_animal,1}, ...
        tp_rate_animal{curr_animal,1},'color',[0,0.4,0],'linestyle',ls);
    plot(fp_rate_animal{curr_animal,2}, ...
        tp_rate_animal{curr_animal,2},'color',[0.4,0,0],'linestyle',ls);
end

for curr_scope = 2:2
    switch curr_scope
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(fp_rate_scope{curr_scope,1}, ...
        tp_rate_scope{curr_scope,1},'color',[0,0.8,0],'linestyle',ls,'linewidth',3);
    plot(fp_rate_scope{curr_scope,2}, ...
        tp_rate_scope{curr_scope,2},'color',[0.8,0,0],'linestyle',ls,'linewidth',3);
end


xlabel('False positive rate');
ylabel('True positive rate');


subplot(1,2,2); hold on;
for curr_animal = 1:length(all_animals)
    switch animal_scope(curr_animal)
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(thresholds,fp_rate_animal{curr_animal,1},'color',[0,0.4,0],'linestyle',ls);
    plot(thresholds,fp_rate_animal{curr_animal,2},'color',[0.4,0,0],'linestyle',ls);
end

for curr_scope = 2:2
    switch curr_scope
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(thresholds,fp_rate_scope{curr_scope,1},'color',[0,0.8,0],'linestyle',ls,'linewidth',3);
    plot(thresholds,fp_rate_scope{curr_scope,2},'color',[0.8,0,0],'linestyle',ls,'linewidth',3);
end

line(xlim,[0.05 0.05],'color','k','linewidth',3);

xlabel('Threshold');
ylabel('False positive rate');


% Get breakdown of events
figure; hold on;
% Plot events by animal
for curr_animal = 1:length(all_animals)
    animal = all_animals{curr_animal};
    curr_scope = event_scores(find(cellfun(@(x) ...
        strcmp(x,animal),{event_scores.animal}),1)).scope;
    curr_n = sum(cellfun(@(x) strcmp(x,animal),{event_scores.animal}));
    switch curr_scope
        case 1
            col = 'm';
        case 2 
            col = 'b';
    end
    bar(curr_animal,curr_n,'FaceColor',col);
end
ylabel('Number of Events');
xlabel('Animal')



%% Compare early vs. late ROC of grand combined (to get effect of filled)

% Get early/late events
event_sessions = [event_scores.session];
early = event_sessions <= 3;

% Set thresholds to draw ROC curve from
thresholds = [2:0.1:10];

% Get the actual positive and number of score events
score_positive = [event_scores.score] == 1;
score_negative = [event_scores.score] == -1;

% Initialize true/false positives, seperate animals and gad
true_positive = cell(2,1);
false_positive = cell(2,1);

% Loop through thresholds, get T/F+ rates
for curr_threshold_idx = 1:length(thresholds)
    
    curr_threshold = thresholds(curr_threshold_idx);
    
    for earlylate = 1:2
        
        use_events = logical(earlylate-early);
        
        % Get events rated positive by binary threshold
        thresh_positive = [event_scores.std] > curr_threshold;
        
        % Get T/F+ rates
        true_positive{earlylate}(curr_threshold_idx,1:2) = ...
            [sum(thresh_positive(use_events) & ...
            score_positive(use_events)) sum(score_positive(use_events))];
        
        false_positive{earlylate}(curr_threshold_idx,1:2) = ...
            [sum(thresh_positive(use_events) & ...
            score_negative(use_events)) sum(score_negative(use_events))];
        
    end
    
    disp(curr_threshold_idx / length(thresholds));
    
end

tp_rate = cellfun(@(x) x(:,1)./x(:,2),true_positive,'uni',false);
fp_rate = cellfun(@(x) x(:,1)./x(:,2),false_positive,'uni',false);

% Plot 
figure; 

subplot(1,2,1); hold on;
for earlylate = 1:2
    switch earlylate
        case 1
            ls = '-';
        case 2
            ls = '--';
    end
    plot(fp_rate{earlylate}, ...
        tp_rate{earlylate},'color','k','linestyle',ls,'linewidth',2);
end

xlabel('False positive rate');
ylabel('True positive rate');
legend({'Early' 'Late'})

subplot(1,2,2); hold on;
for earlylate = 1:2
    switch earlylate
        case 1
            ls = '-';
        case 2
            ls = '--';
    end
    plot(thresholds,fp_rate_animal{earlylate},'color','k','linestyle',ls,'linewidth',2);
end

xlabel('Threshold');
ylabel('False positive rate');
legend({'Early' 'Late'})


%% Plot the ratio of true to false positives

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
ft_ratio_scope = cell(size(tp_scope));
ft_ratio_scope(~empty_conds) = cellfun(@(t,f) f(:,1)./(f(:,1)+t(:,1)), ...
    tp_scope(~empty_conds),fp_scope(~empty_conds),'uni',false);

% Plot false / total positive by threshold

figure;
hold on;
for curr_animal = 1:length(all_animals)
    switch animal_scope(curr_animal)
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(thresholds,ft_ratio_animal{curr_animal,1},'color',[0,0.4,0],'linestyle',ls);
    plot(thresholds,ft_ratio_animal{curr_animal,2},'color',[0.4,0,0],'linestyle',ls);
end

for curr_scope = 2:2
    switch curr_scope
        case 1
            ls = '--';
        case 2
            ls = '-';
    end
    plot(thresholds,ft_ratio_scope{curr_scope,1},'color',[0,0.8,0],'linestyle',ls,'linewidth',3);
    plot(thresholds,ft_ratio_scope{curr_scope,2},'color',[0.8,0,0],'linestyle',ls,'linewidth',3);
end

xlabel('Threshold');
ylabel('False / Total positive rate');















