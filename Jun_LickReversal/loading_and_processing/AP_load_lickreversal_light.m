function [data_all analysis mice] = AP_load_lickreversal_light(animals)
% [mice,analysis] = AP_load_lickreversal_light(animals);
% 'animals' input is cell array of animal names, i.e. {'JL100','JL101'}
%
% Load in and analyze lickreversal data, retain only thresholded/aligned
% Save memory, prevent matlab closures by only saving things to be directly
% used in analysis



% get mice data from mouse data script
mice_all = JL_mousedata;

% get index from 'mice' structure which pertains to which animal
mice_idx = cellfun(@(curr_mouse) ...
    find(arrayfun(@(x) strcmp(mice_all(x).name,curr_mouse),1:length(mice_all))),animals);

% pull out only the mice which were just loaded in
mice = mice_all(mice_idx);

% set up structures for saving
data_all = struct; % only used for framerate
analysis = struct;

% loop through animals, load and process data, store only processed data
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
    
    % if quality control performed on this animal (post-hoc identification
    % of contaminated ROIs), load in and flag bad cells
    qc_dir = '/usr/local/lab/People/Jun/Data/qualitycontrol_bad';
    qc_filename = [qc_dir filesep animal 'bad.roilabel'];
    if exist(qc_filename,'file')
        load(qc_filename,'-mat','roi_labels');
        bad_rois = find(cellfun(@(x) any(strcmp('bad',x)),roi_labels)');
        use_rois = setdiff(1:length(roi_labels),bad_rois);
    else
        warning(['No quality control file: ' animal]);
    end
      
    % load in current animal's data
    for curr_session = 1:length(processed_filenames);
        
        %% Load data
        curr_data = load(processed_filenames{curr_session});
        % ignore the first trial in all analyses
        curr_data.bhv = ...
            structfun(@(x) x(2:end),curr_data.bhv,'uni',false);
        
        
        %% Threshold data
        switch mice(curr_animal).scope
            case 1
                thresh = 3;
            case 2
                thresh = 3;
        end
        
        curr_concat_data = ...
            horzcat(curr_data.im.roi_trace_df{:});
        
        % exclude flagged bad rois
        if ~exist('use_rois','var')
            use_rois = true(size(curr_concat_data,1),1);
        end
        curr_concat_data = curr_concat_data(use_rois,:);
        
        % threshold events
        curr_caEvents = AP_caEvents_thresh(curr_concat_data,thresh);
        
        % Split the activity trace as the df trace is split
        df_split_sizes = cellfun(@(x) size(x,2), ...
            curr_data.im.roi_trace_df);
        
        curr_caEvents_split = mat2cell(curr_caEvents,size(curr_caEvents,1), ...
            df_split_sizes);
        

        %% User settings and common values
        
        % User: Define time to pull out in sec (0 = odor onset, include baseline)
        epoch_time = [-2 4];
        
        switch mice(curr_animal).scope
            case 1
                loop_frames = 1000;
            case 2
                loop_frames = 4000;
        end
        
        epoch_frames = round(epoch_time * curr_data.im.framerate);
        
        % get trials with 'discrimination' paradigm that were also imaged for the
        % course of the epoch time
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            curr_data.bhv.session_type,curr_data.bhv.bhv_frames) &...
            curr_data.bhv.imaged_trials;
        
        % get number of usable trials and cells
        num_trials = sum(discrimination_trials);
        
        % skip if number of trials is 0
        if num_trials == 0
            continue
        end
        
        
        %% Align data

        
        % pull out the loop in which usable trials occurred
        discrimination_trials_loop = curr_data.bhv.xsg_trials(discrimination_trials);
        
        % get the start time of odor onset and answer period (to nearest frame)
        odor_onsets = cellfun(@(x) round(x.states.apply_odor(1)), ...
            curr_data.bhv.bhv_frames(discrimination_trials));
        answer_onsets = cellfun(@(x) round(x.states.answer(1)), ...
            curr_data.bhv.bhv_frames(discrimination_trials));
        
        % get aligned data faster (vectorized) - the reshape command is
        % there because sometimes df organized vertically, sometimes
        % horizontally (loop df vs concat df)
        cat_data_boundaries = [0;reshape(cumsum(cellfun(@(x) size(x,2), ...
            curr_data.im.roi_trace_df)),[],1)];
        
        grab_odor_frames = cell2mat(cellfun(@(x,y) x+epoch_frames(1)+ ...
            cat_data_boundaries(y):x+epoch_frames(2)+cat_data_boundaries(y), ...
            num2cell(odor_onsets),num2cell(discrimination_trials_loop),'uni',false))';
        
        % grab data from thresholded cat data (not sure what best thresh yet)
        curr_concat_data = horzcat(curr_caEvents_split{:});
        
        odor_aligned_df = permute(reshape(curr_concat_data(:,grab_odor_frames(:)),...
            size(curr_concat_data,1),[],num_trials),[3 2 1]);
        
        
        %% Get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),curr_data.bhv.bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),curr_data.bhv.bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),curr_data.bhv.bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),curr_data.bhv.bhv_frames(discrimination_trials));
        
        condition_trials = ...
            [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
        
        applied_odor_group = curr_data.bhv.applied_odor(discrimination_trials);
        rewarded_odor_group = curr_data.bhv.rewarded_odor(discrimination_trials);
        
        odor_trials = [applied_odor_group rewarded_odor_group];
        
      
        %% Separate data by trial type and epoch
        
        % The alternative way, if everything had data:
        all_reversals = find(diff(odor_trials(:,2)) ~= 0)+1;
        
        % Get the start/ends of epoch (either within reversal or entire day if not)
        epoch_boundaries = ...
            mat2cell([1;reshape([all_reversals-1 all_reversals]',[],1);size(odor_trials,1)], ...
            repmat(2,length(all_reversals)+1,1));
        
        % Get the contingencies and days corresponding to reversals
        epoch_contingencies = cellfun(@(epoch) ...
            odor_trials(epoch(1),2), ...
            epoch_boundaries);
        
        epoch_sessions = repmat(data_days(curr_session),length(epoch_contingencies),1);
        
        % Pull out the activity associated with epochs
        epoch_activity_odor = ...
            cellfun(@(epochtrials) odor_aligned_df(epochtrials(1):epochtrials(2),:,:), ...
            epoch_boundaries,'uni',false);
        
        epoch_activity_CL = ...
            cellfun(@(epochtrials) odor_aligned_df(intersect(find(condition_trials(:,1)), ...
            epochtrials(1):epochtrials(2)),:,:),epoch_boundaries,'uni',false);
        
        epoch_activity_CR = ...
            cellfun(@(epochtrials) odor_aligned_df(intersect(find(condition_trials(:,2)), ...
            epochtrials(1):epochtrials(2)),:,:),epoch_boundaries,'uni',false);
        
        
        %% Get lick data
        
        lick_times = cellfun(@(x) x.pokes.C(:,1) - x.states.apply_odor(1), ...
            curr_data.bhv.bhv_frames(discrimination_trials),'uni',false);
               
        % Make link rate trace
        max_lick_time = 2; % seconds after odor application to analyze
        lick_resolution = 100; % in 1/x seconds (i.e max hz)
              
        num_trials = length(lick_times);
        lick_rate = ...
            nan(num_trials,max_lick_time*lick_resolution);
        
        for curr_trial = 1:num_trials
            
            curr_framerate = curr_data.im.framerate;
            curr_licks = lick_times{curr_trial};
            curr_lick_sec_full = curr_licks/curr_framerate;
            curr_lick_sec = curr_lick_sec_full(curr_lick_sec_full > (1/lick_resolution)/2 & ...
                curr_lick_sec_full < max_lick_time);
            % if there's no licks left after this, skip trial
            if isempty(curr_lick_sec)
                curr_lick_rate = zeros(1,max_lick_time*lick_resolution);
                lick_rate(curr_trial,:) = curr_lick_rate;
                continue
            end
            % Initialize continuous lick rate trace:
            % First lick and before = 0
            % Last lick and after = until next lick in trial (if none, then 0);
            curr_lick_rate = nan(1,max_lick_time*lick_resolution);
            curr_lick_diff = 1./diff(curr_lick_sec);
            
            final_lick = find(curr_lick_sec_full == curr_lick_sec(end));
            if length(curr_lick_sec_full) > final_lick + 1
                last_lick_rate = 1./(curr_lick_sec_full(final_lick + 1) - ...
                    curr_lick_sec(end));
            else
                last_lick_rate = 0;
            end
            
            curr_lick_rate(1:round(curr_lick_sec(1)*lick_resolution)) = 0;
            curr_lick_rate(round(curr_lick_sec(end)*lick_resolution):end) = ...
                last_lick_rate;
            curr_lick_rate(round(curr_lick_sec(1:end-1)*lick_resolution)) = ...
                curr_lick_diff;
            
            % Interpolate lick rate across nans
            curr_lick_rate_interp = curr_lick_rate;
            curr_lick_rate_interp(isnan(curr_lick_rate_interp)) = ...
                interp1(find(~isnan(curr_lick_rate_interp)), ...
                curr_lick_rate_interp(~isnan(curr_lick_rate_interp)), ...
                find(isnan(curr_lick_rate_interp)));
            
            lick_rate(curr_trial,:) = curr_lick_rate_interp;
                    
        end
        
        
        
        %% Store data
               
        % store aligned and partitioned data
        analysis.condition_trials{curr_animal}{curr_session} = condition_trials;
        analysis.odor_trials{curr_animal}{curr_session} = odor_trials;
        analysis.epoch_frames{curr_animal}{curr_session} = epoch_frames;
        analysis.epoch_activity_odor{curr_animal}{curr_session} = epoch_activity_odor;
        analysis.epoch_contingencies{curr_animal}{curr_session} = epoch_contingencies;
        analysis.epoch_sessions{curr_animal}{curr_session} = epoch_sessions;
        analysis.epoch_activity_CL{curr_animal}{curr_session} = epoch_activity_CL;
        analysis.epoch_activity_CR{curr_animal}{curr_session} = epoch_activity_CR;
        
        % store framerate        
        data_all(curr_animal).im(curr_session).framerate = curr_data.im.framerate;
        
        % store cells that were not contaminated (after quality control)
        data_all(curr_animal).use_rois = use_rois;
        
        % store lick information
        analysis.lick_times{curr_animal}{curr_session} = lick_times;
        analysis.lick_rate{curr_animal}{curr_session} = lick_rate;
        
        disp(['Session ' num2str(curr_session) ' animal ' animal]);
    end
   
    clearvars -except data_all analysis mice curr_animal animals
end
disp('Finished loading and processing');



















