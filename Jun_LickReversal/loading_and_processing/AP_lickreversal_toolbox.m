%% Info
%
% The initialization part of this script pulls out relevant information for
% easy access. The important variables are here: 
% 
%   %%% Load in mouse data %%%
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
%
%   %%% Threshold traces %%%
%
% - thresholds traces, adds "roi_trace_df_thresh" field to data_all.im structure
%
%   %%%%%%%%%%%%%%%%%%%%
%   All below are considered 'analysis' and placed in an analysis
%   structure, i.e. all variables are called as (analysis.variable)
%   %%%%%%%%%%%%%%%%%%%%
%
%
%   %%% Grab trial-aligned activity & conditions %%%
% 
% - use_trials: trials which are discrimination, contained entirely by the
%     current loop, and imaged
%
% - epoch_frames: the frames which are pulled out for aligned analysis
%     relative to odor onset = 0 (i.e. negative frames = baseline)
%
% - aligned_df_all: imaging data aligned by trial to odor onset + baseline
%     (# of animals -> # of sessions)
%
% - aligned_df_thresh_all - same as aligned_df_all but thresholded
% 
% - odor_trials_all: odor information, column 1 = applied odor, column 2 = rewarded odor
%     (# of animals -> # of sessions)
%     
% - condition_trials_all: classification of trial as correct/incorrect and 
%     lick/rejection. Columns in order: CL, CR, IL, IR (each row = trial, 
%     which is classified as 1 of the 4 types)
%
%   %%% Parse aligned activity into CL and CR %%%
%
% - epoch_contingencies: cell array (animals) of cells (sessions), the 
%     contingencies (1,2 = rewarded odor) for each epoch
%
% - epoch_sessions: cell array (animals) of cells (sessions), the 
%     sessions (1,2 = rewarded odor) for each epoch
%
% - epoch_activity_CL: aligned CL activity by epoch
%
% - epoch_activity_CR: aligned CR activity by epoch
%
%   %%% Grab licks / lick rate %%%
%
% - lick_rate: licking rate by trial (# of animals -> # of sessions)

%% Load and process data - LIGHT
% If this is run, do not run the others (everything in this step)

animals = {'JL072' 'JL073' 'JL074' 'JL075' 'JL077' 'JL078' 'JL081' ...
    'JL082' 'JL088' 'JL092' 'JL093' 'JL094' 'JL100' 'JL101' ...
    'JL102' 'JL113' 'JL116' 'JL119' 'JL104' 'JL107' 'JL123' ...
    'JL124' 'JL126' 'JL127' 'JL129'};

[data_all,analysis,mice] = AP_load_lickreversal_light(animals);


%% Load in mouse data
% NOTE!!!! THIS DOES NOT CURRENTLY LOAD QUALITY CONTROL LABELS and also
% uses the old threshold

animals = {'JL072' 'JL073' 'JL074' 'JL075' 'JL077' 'JL078' 'JL081' ...
    'JL082' 'JL088' 'JL092' 'JL093' 'JL094' 'JL098' 'JL100' 'JL101' ...
    'JL102' 'JL113' 'JL116'};

[data_all mice] = AP_load_lickreversal_data(animals);

analysis = struct;

clearvars -except mice data_all analysis

%% Threshold Data

for curr_animal = 1:length(data_all)
    
   switch mice(curr_animal).scope
       case 1
           thresh = 4;
       case 2
           thresh = 4;
   end
    
   for curr_session = 1:length(data_all(curr_animal).im)
       
       curr_concat_data = ...
           horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});
       curr_caEvents = AP_caEvents_thresh(curr_concat_data,thresh);
       
       % Split the activity trace as the df trace is split
       df_split_sizes = cellfun(@(x) size(x,2), ...
           data_all(curr_animal).im(curr_session).roi_trace_df);
       
       curr_caEvents_split = mat2cell(curr_caEvents,size(curr_caEvents,1), ...
           df_split_sizes);
       
       data_all(curr_animal).im(curr_session).roi_trace_df_thresh = ...
           curr_caEvents_split;
   end    
end
disp('Finished thresholding data');

%% Grab trial-aligned activity & conditions

% User: Define time to pull out in sec (0 = odor onset, include baseline) 
epoch_time = [-2 10];

% Activity
use_trials = cell(size(mice));
epoch_frames_all = cell(size(mice));
aligned_df_all = cell(size(mice));
for curr_animal = 1:length(mice)
    for curr_session = 1:length(data_all(curr_animal).im);
        
        switch mice(curr_animal).scope
            case 1
                loop_frames = 1000;
            case 2
                loop_frames = 4000;
        end
              
        epoch_frames = round(epoch_time * data_all(curr_animal).im(curr_session).framerate);
        
        % get trials with 'discrimination' paradigm that were also imaged for the
        % course of the epoch time
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        use_trials{curr_animal}{curr_session} = discrimination_trials;
        
        % get number of usable trials and cells
        num_trials = sum(discrimination_trials);
        num_cells = size(data_all(curr_animal).im(curr_session).roi_trace_df{1},1);
        
        % skip if number of trials is 0
        if num_trials == 0
            continue
        end
        
        % pull out the loop in which usable trials occurred
        discrimination_trials_loop = data_all(curr_animal).bhv(curr_session).xsg_trials(discrimination_trials);
        
        % get the start time of odor onset and answer period (to nearest frame)
        odor_onsets = cellfun(@(x) round(x.states.apply_odor(1)), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        answer_onsets = cellfun(@(x) round(x.states.answer(1)), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        % get aligned data faster (vectorized) - the reshape command is
        % there because sometimes df organized vertically, sometimes
        % horizontally (loop df vs concat df)
        cat_data_boundaries = [0;reshape(cumsum(cellfun(@(x) size(x,2), ...
            data_all(curr_animal).im(curr_session).roi_trace_df)),[],1)];
        
        grab_odor_frames = cell2mat(cellfun(@(x,y) x+epoch_frames(1)+ ...
            cat_data_boundaries(y):x+epoch_frames(2)+cat_data_boundaries(y), ...
            num2cell(odor_onsets),num2cell(discrimination_trials_loop),'uni',false))';

        % grab data from thresholded cat data (not sure what best thresh yet)
        curr_concat_data = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df_thresh{:});       
        
        odor_aligned_df = permute(reshape(curr_concat_data(:,grab_odor_frames(:)),...
            size(curr_concat_data,1),[],num_trials),[3 2 1]);
        %odor_aligned_df_thresh = permute(reshape(curr_concat_data_thresh(:,grab_odor_frames(:)),...
        %    size(curr_concat_data_thresh,1),[],num_trials),[3 2 1]);
        
        epoch_frames_all{curr_animal}{curr_session} = epoch_frames;
        aligned_df_all{curr_animal}{curr_session} = odor_aligned_df;  
        %aligned_df_thresh_all{curr_animal}{curr_session} = odor_aligned_df_thresh;  
        
        disp(['Session: ' num2str(curr_session)])
    end
    disp(['Animal: ' mice(curr_animal).name])
end
analysis.use_trials = use_trials;
analysis.epoch_frames = epoch_frames_all;
analysis.aligned_df_all = aligned_df_all;
%analysis.aligned_df_thresh_all = aligned_df_thresh_all;

% Trial conditions
condition_trials_all = cell(1,length(mice));
odor_trials_all = cell(1,length(mice));
for curr_animal = 1:length(mice)
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
        
        
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        condition_trials_all{curr_animal}{curr_session} = ...
            [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
        
        applied_odor_group = data_all(curr_animal).bhv(curr_session).applied_odor(discrimination_trials);
        rewarded_odor_group = data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        
        odor_trials_all{curr_animal}{curr_session} = [applied_odor_group rewarded_odor_group];
    end
end
analysis.condition_trials = condition_trials_all;
analysis.odor_trials = odor_trials_all;

disp('Finished getting trial conditions and aligning activity')

clearvars -except mice data_all analysis

%% Parse aligned activity into CL and CR

% Get reversal trials (reversal = first trial of new contingency)
% (put into loops because missing data doesn't have correct indicies
% all_reversals = cell(size(mice));
% for curr_animal = 1:length(mice)
%     all_reversals{curr_animal} = ...
%         cell(size(analysis.odor_trials_all{curr_animal}));
%     for curr_session = 1:length(analysis.odor_trials_all{curr_animal})
%         if isempty(analysis.odor_trials_all{curr_animal}{curr_session})
%             continue
%         end        
%         all_reversals{curr_animal}{curr_session} = ...
%             find(diff(analysis.odor_trials_all{curr_animal}{curr_session}(:,2)) ~= 0)+1; 
%     end
% end

% The alternative way, if everything had data:
all_reversals = cellfun(@(x) cellfun(@(y) find(diff(y(:,2)) ~= 0)+1,x,'uni',false), ...
   analysis.odor_trials,'uni',false); 

% Get the start/ends of epoch (either within reversal or entire day if not)
epoch_boundaries = cellfun(@(rev,odor) cellfun(@(rev,odor) ...
    mat2cell([1;reshape([rev-1 rev]',[],1);size(odor,1)], ...
    repmat(2,length(rev)+1,1)),rev,odor,'uni',false), ...
    all_reversals,analysis.odor_trials,'uni',false);

% Get the contingencies and days corresponding to reversals 
epoch_contingencies = cellfun(@(epoch,odor) cellfun(@(epoch,odor) cellfun(@(epoch)...
    odor(epoch(1),2),epoch),epoch,odor,'uni',false), ...
    epoch_boundaries,analysis.odor_trials,'uni',false);

epoch_sessions = cellfun(@(epoch) cellfun(@(epoch,sessions) ...
    repmat(sessions,length(epoch),1),epoch,num2cell(1:length(epoch)), ...
    'uni',false),epoch_boundaries,'uni',false);

% Pull out the activity associated with epochs
epoch_activity_odor = cellfun(@(df,epochtrials) ...
    cellfun(@(df,epochtrials) ...
    cellfun(@(epochtrials) df(epochtrials(1):epochtrials(2),:,:), ...
    epochtrials,'uni',false), ...
    df,epochtrials,'uni',false),analysis.aligned_df_all, ...
    epoch_boundaries,'uni',false);

epoch_activity_CL = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,1)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),analysis.aligned_df_all, ...
    analysis.condition_trials,epoch_boundaries,'uni',false);

epoch_activity_CR = cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(df,bhv,epochtrials) ...
    cellfun(@(epochtrials) df(intersect(find(bhv(:,2)), ...
    epochtrials(1):epochtrials(2)),:,:),epochtrials,'uni',false), ...
    df,bhv,epochtrials,'uni',false),analysis.aligned_df_all, ...
    analysis.condition_trials,epoch_boundaries,'uni',false);

analysis.epoch_activity_odor = epoch_activity_odor;
analysis.epoch_contingencies = epoch_contingencies;
analysis.epoch_sessions = epoch_sessions;
analysis.epoch_activity_CL = epoch_activity_CL;
analysis.epoch_activity_CR = epoch_activity_CR;

disp('Finished separating trials')
clearvars -except mice data_all analysis


%% Grab licks / lick rate

lick_times = cell(size(mice));
for curr_animal = 1:length(mice)
    for curr_session = 1:length(data_all(curr_animal).bhv);   
        
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

        lick_times{curr_animal}{curr_session} = cellfun(@(x) x.pokes.C(:,1) - x.states.apply_odor(1), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials),'uni',false);
    end
end
analysis.lick_times = lick_times;

% Make link rate trace
max_lick_time = 2; % seconds after odor application to analyze
lick_resolution = 100; % in 1/x seconds (i.e max hz)
lick_rate = cell(size(data_all));
for curr_animal = 1:length(data_all)
    for curr_session = 1:length(data_all(curr_animal).im)
        num_trials = length(lick_times{curr_animal}{curr_session});
        lick_rate{curr_animal}{curr_session} = ...
            nan(num_trials,max_lick_time*lick_resolution);
        for curr_trial = 1:num_trials
            curr_framerate = data_all(curr_animal).im(curr_session).framerate;
            curr_licks = lick_times{curr_animal}{curr_session}{curr_trial};
            curr_lick_sec_full = curr_licks/curr_framerate;
            curr_lick_sec = curr_lick_sec_full(curr_lick_sec_full > (1/lick_resolution)/2 & ...
                curr_lick_sec_full < max_lick_time);
            % if there's no licks left after this, skip trial
            if isempty(curr_lick_sec)
                curr_lick_rate = zeros(1,max_lick_time*lick_resolution);
                lick_rate{curr_animal}{curr_session}(curr_trial,:) = curr_lick_rate;
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

            lick_rate{curr_animal}{curr_session}(curr_trial,:) = curr_lick_rate_interp;
        end
    end
end
analysis.lick_rate = lick_rate;

disp('Finished getting licks')
clearvars -except mice data_all analysis



%% Plot single cell across sessions

curr_animal = 10;
curr_cell = 54;

%%% For all trials 
curr_sessions = length(analysis.aligned_df_all{curr_animal});
figure;set(gcf,'Name','CL')
for curr_session = 1:curr_sessions
    subplot(ceil(sqrt(curr_sessions)),ceil(sqrt(curr_sessions)),curr_session);
    imagesc(analysis.aligned_df_all{curr_animal}{curr_session}( ...
        analysis.condition_trials_all{curr_animal}{curr_session}(:,1),:,curr_cell));colormap(gray)
    % lines for reversal
    curr_CL = find(analysis.condition_trials_all{curr_animal}{curr_session}(:,1));
    curr_rev = find(diff(analysis.odor_trials_all{curr_animal}{curr_session}(curr_CL,2)) ~= 0);
    for i = 1:length(curr_rev)
       line(xlim,[curr_rev(i) curr_rev(i)],'color','r','linestyle','--'); 
    end
end

figure;set(gcf,'Name','CR')
for curr_session = 1:curr_sessions
    subplot(ceil(sqrt(curr_sessions)),ceil(sqrt(curr_sessions)),curr_session);
    imagesc(analysis.aligned_df_all{curr_animal}{curr_session}( ...
        analysis.condition_trials_all{curr_animal}{curr_session}(:,2),:,curr_cell));colormap(gray)
    % lines for reversal
    curr_CR = find(analysis.condition_trials_all{curr_animal}{curr_session}(:,2));
    curr_rev = find(diff(analysis.odor_trials_all{curr_animal}{curr_session}(curr_CR,2)) ~= 0);
    for i = 1:length(curr_rev)
       line(xlim,[curr_rev(i) curr_rev(i)],'color','r','linestyle','--'); 
    end
end



%% Plot single cell within session

curr_animal = 12;
curr_session = 7;
curr_cell = 5;

%%% Pull out odor onset times

% get trials with 'discrimination' paradigm that were also imaged for the
% course of the epoch time
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

% pull out the loop in which usable trials occurred
discrimination_trials_loop = data_all(curr_animal).bhv(curr_session).xsg_trials(discrimination_trials);

% get the start time of odor onset and answer period (to nearest frame)
odor_onsets = cellfun(@(x) x.states.apply_odor(1), ...
    data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));

% get the lick_times
licks = cellfun(@(x) x.pokes.C(:,1), ...
    data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials),'uni',false);

% correct events by loop
sum_loop_frames = cumsum([0;reshape(cellfun(@(x) size(x,2), ...
    data_all(curr_animal).im(curr_session).roi_trace_df),[],1)]);

odor_onsets_loop = odor_onsets + sum_loop_frames(discrimination_trials_loop);
licks_loop = cellfun(@(licks,loop) licks+loop,licks, ...
    num2cell(sum_loop_frames(discrimination_trials_loop)),'uni',false);

cat_trace = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});

% get thresholded data
caEvents = AP_caEvents_momscope(cat_trace(curr_cell,:),[],1);

% plot everything
figure; hold on

% plot the ROI trace
plot(cat_trace(curr_cell,:),'k');

% plot the thresholded data
plot(caEvents,'r');

% plot licks
plot(vertcat(licks_loop{:}),-0.1,'.b');

% plot odor onsets colored
for i = 1:length(analysis.odor_trials_all{curr_animal}{curr_session});
    switch analysis.odor_trials_all{curr_animal}{curr_session}(i,1)
        case 1
            line([odor_onsets_loop(i) odor_onsets_loop(i)],ylim,'color','b');
        case 2
            line([odor_onsets_loop(i) odor_onsets_loop(i)],ylim,'color','c');
    end
end

% plot reversals
curr_reversals = find(diff(analysis.odor_trials_all{curr_animal}{curr_session}(:,2)) ~= 0);
for i = curr_reversals'
   line([odor_onsets_loop(i) odor_onsets_loop(i)],ylim,'color','r','linestyle','--');    
end

% plot correct trials
curr_correct = find(any(analysis.condition_trials_all{curr_animal}{curr_session}(:,1:2),2));
for i = curr_correct';
    switch diff(analysis.odor_trials_all{curr_animal}{curr_session}(i,:)) == 0
        case true
            plot(odor_onsets_loop(i),0,'ok','MarkerSize',7,'MarkerFaceColor','g'); 
        case false
            plot(odor_onsets_loop(i),0,'ok','MarkerSize',7,'MarkerFaceColor','r'); 
    end
end

% plot incorrect trials
curr_incorrect = find(any(analysis.condition_trials_all{curr_animal}{curr_session}(:,3:4),2));
for i = curr_incorrect';
    switch diff(analysis.odor_trials_all{curr_animal}{curr_session}(i,:)) == 0
        case true
            plot(odor_onsets_loop(i),0,'og','MarkerSize',7,'MarkerFaceColor','w'); 
        case false
            plot(odor_onsets_loop(i),0,'or','MarkerSize',7,'MarkerFaceColor','w'); 
    end
end


%% LIGHT: Plot single cell across sessions

curr_animal = 9;
curr_cell = 131;

%%% For all trials 
curr_sessions = length(analysis.epoch_activity_CL{curr_animal});
figure;set(gcf,'Name','CL')
for curr_session = 1:curr_sessions
    subplot(ceil(sqrt(curr_sessions)),ceil(sqrt(curr_sessions)),curr_session);
    
    curr_activity = vertcat(analysis.epoch_activity_CL{curr_animal}{curr_session}{:});
    imagesc(curr_activity(:,:,curr_cell));
    colormap(gray);
    
    % lines for reversal
    curr_rev = cellfun(@(x) size(x,1),analysis.epoch_activity_CL{curr_animal}{curr_session});
    for i = 1:length(curr_rev)
       line(xlim,[curr_rev(i) curr_rev(i)],'color','r','linestyle','--'); 
    end
end








































