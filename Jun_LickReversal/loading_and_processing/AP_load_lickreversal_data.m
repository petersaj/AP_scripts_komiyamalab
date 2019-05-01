function [data_all mice] = AP_load_lickreversal_data(animals,sessions)
% [data_all mice] = AP_load_lickreversal_data(animals);
%
% 'animals' input must be cell array of strings i.e. {'JL072' 'JL073'}
%
% Load in all data from JL lickreversal
%
% - mice: mice data entered in JL_mousedata (of loaded animals, in
%      whichever order they're loaded)
%
% - data_all: all data from all animals (# of animals)
%     - im: structure of imaging data (# of sessions)
%         - framerate
%         - roi_trace_split: raw imaging data
%         - roi_trace_split_bg: raw background imaging data
%         - roi_trace_df: background subtracted, df/f imaging data
%     - bhv: structure of behavior data (# of sessions)
%         - various features of behavior events, timing, and loop assignments
%     - session: the session number for each included session
%         - NOTE: this is important for any animals which were missing a
%         session: should be tracked e.g. as [1 2 4] to include that
%         session 3 happened but doesn't have data
%     - labels: one structure for each label
%         - currently only gad
%         

% Check w/ Jun that this list is comprehensive

% get mice data from mouse data script
mice_all = JL_mousedata;

% get index from 'mice' structure which pertains to which animal
mice_idx = cellfun(@(curr_mouse) ...
    find(arrayfun(@(x) strcmp(mice_all(x).name,curr_mouse),1:length(mice_all))),animals);

% pull out only the mice which were just loaded in 
mice = mice_all(mice_idx);

data_all = struct('im',cell(size(animals)),'bhv',cell(size(animals)), ...
    'session',cell(size(animals)));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    processed_path = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_processedData'];
    processed_dir = dir([processed_path filesep '*processed.mat']);
    processed_filenames = sort(cellfun(@(x) [processed_path filesep x],{processed_dir.name},'uni',false));
    processed_days = cellfun(@(x) x(7:12), {processed_dir.name},'uni',false);

    % get the number of total days. this is to avoid situations where data
    % is missing but the animal WAS trained, e.g. sessions go 1,2,4 and
    % should be saved as 1,2,[],4
    data_path = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_roi'];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    days = cellfun(@(x) x(end-5:end),data_names(3:end),'UniformOutput',false);
    
    % match processed days to image file days
    data_days = cellfun(@(x) find(cellfun(@(y) strcmp(x,y),days)),processed_days);
    
    if nargin == 2
        use_sessions = sessions;
    else
        use_sessions = 1:length(processed_filenames);
    end
    
    % load in current animal's data, organize in master data structure
    for curr_session_idx = 1:length(use_sessions);
        
        curr_session = use_sessions(curr_session_idx);
        
        curr_data = load(processed_filenames{curr_session});
        data_all(curr_animal).im(curr_session_idx) = curr_data.im;
        % ignore the first trial in all analyses
        data_all(curr_animal).bhv(curr_session_idx) = ...
            structfun(@(x) x(2:end),curr_data.bhv,'uni',false);
        data_all(curr_animal).session(curr_session_idx) = ...
        data_days(curr_session);
        disp(['Session ' num2str(curr_session)]);
    end
    disp(['Loaded ' animal]);
end
disp('Finished loading animals');

%%%
% Below was taken out: one animal had ROIs drawn for just an odor lick day,
% but that day was removed, so there should not be any more instances of
% this (all days with data at all should contain discrimination trials)
%
% NOTE: was included again because sometimes back to odor lick at the end
% on purpose (but not included in the analysis of this script currently).
% No way to tell then whether on purpose or not except for timing: if it
% occurs in the first session, error out and display. If otherwise, just
% don't include in analysis.
%%%

% cut out all data with no discrimination trials on that day
for curr_animal = 1:length(mice)
    curr_cut_sessions = false(size(data_all(curr_animal).im));
    for curr_session = 1:length(data_all(curr_animal).im);   
        
        switch mice(curr_animal).scope
            case 1
                loop_frames = 1000;
            case 2
                loop_frames = 4000;
        end
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        if sum(discrimination_trials) == 0;
            if curr_session == 1;
                error([mice(curr_animal).name ': First session not discrimination']);
            end
            curr_cut_sessions(curr_session) = true;
        end
    end
    % cut out the marked sessions
    data_all(curr_animal).im(curr_cut_sessions) = [];
    data_all(curr_animal).bhv(curr_cut_sessions) = [];
    data_all(curr_animal).session(curr_cut_sessions) = [];
end

% Pick out labels (currently only gad, maybe filled in future?) cells for each animal
% adds a "label" field to data_all
for curr_animal = 1:length(mice)
    % get the number of cells: this isn't written anywhere explicitly, just
    % get the maximum number of rows from all sessions first loops
    % (anything imaged should have the same number)
    data_sessions = arrayfun(@(x) ~isempty( ...
        data_all(curr_animal).im(x).roi_trace_df), ...
        1:length(data_all(curr_animal).im));
    
    num_cells = max(cellfun(@(x) size(x,1),cellfun(@(x) x{1},...
        {data_all(curr_animal).im(data_sessions).roi_trace_df},'uni',false)));
    
    curr_gad = false(num_cells,1);
    curr_gad(mice(curr_animal).interneurons) = true;
    
    data_all(curr_animal).labels.gad = curr_gad;
end
disp('Finished getting labels')

