function data = AP_corticospinal_load_all(animals)
% data = AP_corticospinal_load_all(animals);


data = struct('animal',cell(length(animals),1), 'bhv',cell(length(animals),1), ...
    'im',cell(length(animals),1), 'roi_group',cell(length(animals),1));


% Grab anaylsis files (assume dendrites in specific folder)
all_analysis_path = cell(size(animals));
all_analysis_files = cell(size(animals));
if ~ispc
    data_dir = '/usr/local/lab/People/Andy/Data';
else
    data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\analysis_files';
end
for curr_animal = 1:length(animals)
    if ~ispc
        curr_analysis_dir = [data_dir filesep animals{curr_animal} filesep ...
            animals{curr_animal} '_batch_thresh_roi'];
    else
        curr_analysis_dir = [data_dir filesep animals{curr_animal}];
    end
    
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    
    all_analysis_path{curr_animal} = curr_analysis_dir;
    all_analysis_files{curr_animal} = sort({curr_analysis_files.name});   
end

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(['Loading animal ' animal]);
    
    data(curr_animal).animal = animal;
    
    analysis_path = all_analysis_path{curr_animal};
    analysis_files = sort(all_analysis_files{curr_animal});
    num_sessions = length(analysis_files);
    
    
    %% Load ROI data
    
    disp('Loading data...')
    
    for curr_session = 1:num_sessions
        
        curr_file = [analysis_path filesep analysis_files{curr_session}];
        curr_data = load(curr_file);
        
        % Store data
        data(curr_animal).im(curr_session) = curr_data.im;
        % because one day missing bhv data, go through each field of
        % behavior structure and update independently
        bhv_fieldnames = fieldnames(curr_data.bhv);
        for curr_field = 1:length(bhv_fieldnames)
            data(curr_animal).bhv(curr_session).(bhv_fieldnames{curr_field}) = ...
                curr_data.bhv.(bhv_fieldnames{curr_field});
        end
        
    end
    
    %% Get bad ROIs from manual quality control, discard
    
    field_qc_path = '/usr/local/lab/People/Andy/Data/field_qc';
    field_qc_dir = dir([field_qc_path filesep '*.mat']);
    
    curr_animal_field_qc = find(cellfun(@(x) strncmp(animal,x,5),{field_qc_dir.name}));
    if ~isempty(curr_animal_field_qc)
        
        curr_bad_rois = load([field_qc_path filesep ...
            field_qc_dir(curr_animal_field_qc).name]);
        
        for curr_session = 1:length(data(curr_animal).im)
            data(curr_animal).im(curr_session).roi_trace(curr_bad_rois.bad_rois,:) = [];
            data(curr_animal).im(curr_session).roi_trace_bg(curr_bad_rois.bad_rois,:) = [];
            data(curr_animal).im(curr_session).roi_trace_df(curr_bad_rois.bad_rois,:) = [];
        end
    end
    data(curr_animal).roi_idx = find(~curr_bad_rois.bad_rois);
    
    
    %% Find likely same dendrite ROIs / consistent high correlation
    
    % Get the 'correlation' of concatentenated, smoothed velocity
    
    disp('Grouping ROIs...')
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    
    % Loop through days, get L2 normalization value (avoid concat data)
    cat_sum = nan(n_rois,num_sessions);
    for curr_day = 1:num_sessions
        
        curr_smoothpeak = nan(size(data(curr_animal).im(curr_day).roi_trace_df));
        for curr_roi = 1:size(curr_smoothpeak,1)
            curr_smoothpeak(curr_roi,:) = smooth(data(curr_animal).im(curr_day). ...
                roi_trace_df(curr_roi,:),100,'loess');
        end
        cat_sum(:,curr_day) = sum(curr_smoothpeak.^2,2);
    end
    % Loop through days, get total correlation between ROIs
    l2norm = sqrt(sum(cat_sum,2));
    roi_corr = zeros(n_rois,n_rois);
    for curr_day = 1:num_sessions
        
        curr_smoothpeak = nan(size(data(curr_animal).im(curr_day).roi_trace_df));
        for curr_roi = 1:size(curr_smoothpeak,1)
            curr_smoothpeak(curr_roi,:) = smooth(data(curr_animal).im(curr_day). ...
                roi_trace_df(curr_roi,:),100,'loess');
        end
        
        curr_smoothpeak = bsxfun(@times,curr_smoothpeak,1./l2norm);
        roi_corr = roi_corr + curr_smoothpeak*curr_smoothpeak';
    end
    
    roi_corr = tril(roi_corr,-1);
    [d1,d2] = find(roi_corr >= 0.8);
    d_pairs = [d1 d2];
    
    % Group by multi-way combinations
    d_grp = {};
    all_d = unique([d1;d2]);
    for curr_d = all_d'
        
        old_grp = cellfun(@(x) ismember(curr_d,x),d_grp);
        
        if ~any(old_grp)
            curr_grp = length(d_grp) + 1;
            d_grp{curr_grp,1} = [];
        else
            % NOTE: this shouldn't yield more than 1 grp, but sometimes does
            % so do something about this later
            curr_grp = find(old_grp,1);
        end
        
        curr_pair = d_pairs == curr_d;
        paired_d = setdiff(reshape(d_pairs(any(curr_pair,2),:),[],1),curr_d);
        
        d_grp{curr_grp,1} = unique([d_grp{curr_grp};curr_d;paired_d]);
        
    end
    ungrouped = setdiff([1:size(roi_corr,1)]',unique(vertcat(d_grp{:})));
    d_grp = [d_grp;num2cell(ungrouped)];
    
    data(curr_animal).roi_group = d_grp;
    
    
    % Combine same dendrites in loaded data by weighted average using L2
    % norm generated above (more signal = more weight)
    for curr_session = 1:num_sessions       
        data(curr_animal).im(curr_session).roi_trace_df = ...
            cell2mat(arrayfun(@(x) sum(bsxfun(@times,data(curr_animal).im( ...
            curr_session).roi_trace_df(data(curr_animal).roi_group{x},:), ...
            l2norm(data(curr_animal).roi_group{x})),1)/sum(l2norm(data(curr_animal).roi_group{x})), ...
            [1:length(data(curr_animal).roi_group)]','uni',false));       
    end
    
    % Don't keep unused data loaded in to save space
    if isfield(data(curr_animal).im,'roi_trace');
        data(curr_animal).im = rmfield(data(curr_animal).im,'roi_trace');
    end
    if isfield(data(curr_animal).im,'roi_trace_bg');
        data(curr_animal).im = rmfield(data(curr_animal).im,'roi_trace_bg');        
    end
        
    
    %% Threshold data
    disp('Thresholding...')
    for curr_session = 1:num_sessions
        event_thresh = 3;
        method = 2;
        data(curr_animal).im(curr_session).roi_trace_thresh = ...
            AP_caEvents_thresh(data(curr_animal).im(curr_session).roi_trace_df,event_thresh,method);
    end
    disp(['Finished ' animals{curr_animal}]);
end
disp('Finished all');







