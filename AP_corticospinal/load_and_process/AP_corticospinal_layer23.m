% Load and process layer 2/3 data using updated protocols and methods for
% direct comparison to corticospinal data

% NOTE: CURRENTLY NOTHING IS CHANGED
% Things to change: 
% - AP_parseLeverMovement -> AP_parseLeverMovement_updated
% - classify ROIs by new method



%% PUT IN MEMORY: activity/classification, movements/dispatcher

animals = {'AP71' 'AP72' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

gad_activity_thresh_all = cell(length(animals),15);
gad_class_all = cell(length(animals),15);

pyr_activity_all = cell(length(animals),15);
pyr_class_all = cell(length(animals),15);

pyr_class_rewardmove = cell(length(animals),15);
gad_class_rewardmove = cell(length(animals),15);

movement_trace_all = cell(length(animals),15);
move_epoch_frames = cell(length(animals),15);
rewarded_movements = cell(length(animals),15);
cued_rewarded_movements = cell(length(animals),15);
catch_rewarded_movements = cell(length(animals),15);

lick_times_all = cell(length(animals),15);
lick_trace_all = cell(length(animals),15);

cue_frames_all = cell(length(animals),15);
reward_frames_all = cell(length(animals),15);
cued_reward_frames_all = cell(length(animals),15);
catch_reward_frames_all = cell(length(animals),15);

trial_frames_all = cell(length(animals),15);

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    % load ROIs and get ROI centroids and distances
    roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    roi_dir = dir(roi_path);
    roi_filenames = {roi_dir(:).name};
    roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
    roi_files_sort = sort(roi_filenames(roi_files));
    load([roi_path filesep roi_files_sort{end}],'-MAT');
    num_rois = length(polygon.ROI);
    roi_center_x = zeros(num_rois,1);
    roi_center_y = zeros(num_rois,1);
    for curr_roi  = 1:num_rois
        [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
            polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
    end
    roi_dist = nan(length(roi_center_x));
    for k = 1:length(roi_dist)
        for j = 1:length(roi_dist)
            curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
                (roi_center_y(k) - roi_center_y(j))^2);
            roi_dist(k,j) = curr_dist;
        end
    end
    
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    % load movement classifications for cells
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        'final_class/' animal '_classified_cells_caEvents_10k'])
        
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);

    % Classification based on rewarded moves only
%            
%     clear classified_cells
%     load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
%         'final_class/' animal '_classified_cells_caEvents_rewardedmove_10k'])
%         
%     pyr_class_rewardmove(curr_animal,sessions) = ...
%         cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
%     gad_class_rewardmove(curr_animal,sessions) = ...
%         cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
%         
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);


    if strcmp(animal,'AP98') || strcmp(animal,'AP99') ||strcmp(animal,'AP100')
        days = days(8:end);
        sessions = 1:7;
    end
    
    for curr_day = sessions
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        %%%% Find and store activity
        pyr_activity_all{curr_animal,curr_day} = ...
            AP_caEvents(im.roi_trace_df(pyr_unfilled,:));
        
        gad_activity_thresh_all{curr_animal,curr_day} = ...
        AP_caEvents(im.roi_trace_df(gad_unfilled,:),[],1:sum(gad_unfilled));
       
        %%%% Find and store movements
        
        % get lick frames (in licks per frame);
        lick_trace_all{curr_animal,curr_day} = zeros(1,size(im.roi_trace_df,2));
        lick_trace_all{curr_animal,curr_day}(round( ...
            bhv.lick_frames(~isnan(bhv.lick_frames)))) = 1;
        lick_times_all{curr_animal,curr_day} = bhv.lick_frames(~isnan(bhv.lick_frames));
        
        num_frames = size(im.roi_trace_df,2);
        
        % load split image file to get loop sizes        
        %loop_frames = unique(bhv.add_frames);
        %loop_size = [diff(loop_frames);num_frames-loop_frames(end)];
        framerate = bhv.framerate;
        clear bhv;
        
        %if loading in image anyway, this is probably safer
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        % get behavior files/data
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            continue
        end
        bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(bhv);
        rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv));        
        
        catch_trials = [];
        if isfield(saved_history,'CatchSection_catch_trial')
            catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
            % there's one or more uncompleted trials at the end
            catch_trials = find(catch_trials(1:length(bhv)));
        end
                
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP' repmat('0',1,4-length(animal(3:end))) animal(3:end)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
             
        lever_active_full = cell(length(xsg_filenames),1);
        lever_velocity_full = cell(length(xsg_filenames),1);
        cue_frames_full =  cell(length(xsg_filenames),1);
        reward_frames_full =  cell(length(xsg_filenames),1);
        cued_reward_frames_full = cell(length(xsg_filenames),1);
        catch_reward_frames_full = cell(length(xsg_filenames),1);
        
        trial_frames = cellfun(@(x) nan(1,x),num2cell(loop_size),'uni',false);
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            [lever_active lever_force_smooth lever_force_velocity] = ...
                AP_parseLeverMovement(xsg_data);
            
            % find which movements are (cued) rewarded
            cue_sample = nan(length(curr_trial_list(:,2)),1);
            reward_start_sample = nan(length(curr_trial_list(:,2)),1);
            reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
            cued_movement = nan(length(curr_trial_list(:,2)),1);
            catch_movement = nan(length(curr_trial_list(:,2)),1);
            for curr_trial_idx = 1:length(curr_trial_list(:,2));
                curr_trial = curr_trial_list(curr_trial_idx,2);
                
                % ignore if unrewarded trial
                if ~ismember(curr_trial,rewarded_trials)
                    continue
                end
                
                % get trial start in dispatcher
                curr_bhv_start = bhv{curr_trial}.states.bitcode(1);
                % get cue and reward sample in xsg (downsampled to ms)
                curr_bhv_cue_sample_rel =  (bhv{curr_trial}.states.cue(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                cue_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_cue_sample_rel);
                
                curr_bhv_reward_start_sample_rel =  (bhv{curr_trial}.states.reward(1) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_start_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_start_sample_rel);
                
                curr_bhv_reward_stop_sample_rel =  (bhv{curr_trial}.states.reward(2) - ...
                    curr_bhv_start)*(xsg_samplerate/10);
                reward_stop_sample(curr_trial_idx) = ...
                    round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                    *(xsg_samplerate/10) + ...
                    curr_bhv_reward_stop_sample_rel);
                
                if reward_start_sample(curr_trial_idx) > length(lever_active)
                    continue
                end
                
                rewarded_move_start = find(lever_active(1: ...
                    reward_start_sample(curr_trial_idx)) == 0,1,'last');
                cued_movement(curr_trial_idx) = ...
                    rewarded_move_start > cue_sample(curr_trial_idx) && ...
                    ~any(lever_active(cue_sample(curr_trial_idx)-100: ...
                    cue_sample(curr_trial_idx)));
                                
                if ismember(curr_trial,catch_trials)
                    catch_movement(curr_trial_idx) = true;
                end
                
            end      
            cue_frames = round((cue_sample/1000)*framerate);
            reward_frames = round((reward_start_sample/1000)*framerate);
            cued_reward_frames = reward_frames(cued_movement == 1);
            catch_reward_frames = reward_frames(catch_movement == 1);
            
            cue_frames_full{curr_xsg} = cue_frames( ...
                cue_frames > 0 & ...
                cue_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);
            reward_frames_full{curr_xsg} = reward_frames( ...
                reward_frames > 0 & ...
                reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);            
            cued_reward_frames_full{curr_xsg} = cued_reward_frames( ...
                cued_reward_frames > 0 & ...
                cued_reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);            
            catch_reward_frames_full{curr_xsg} = catch_reward_frames( ...
                catch_reward_frames > 0 & ...
                catch_reward_frames <= loop_size(curr_xsg)) + loop_frames(curr_xsg);
           
            %save the active lever portions in frames
            curr_lever_active_frames = round((find(lever_active)/1000) ...
                *framerate);
            lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
                curr_lever_active_frames <= loop_size(curr_xsg) & ...
                curr_lever_active_frames > 0)) + loop_frames(curr_xsg);
            
            [n d] = rat(loop_size(curr_xsg)/length(lever_force_velocity));
            lever_force_velocity_resample = resample(lever_force_velocity,n,d);
            lever_velocity_full{curr_xsg} = ...
                lever_force_velocity_resample(1:loop_size(curr_xsg));
            
            % save the trial numbers that each frame falls into
            % convert trial start times to frames
            curr_trial_list_frames = curr_trial_list;
            curr_trial_list_frames(:,1) = round( ...
                curr_trial_list_frames(:,1)*framerate);
            % save trials for each frame
            for curr_trial = 1:size(curr_trial_list_frames,1)-1;
                trial_frames{curr_xsg}(curr_trial_list_frames(curr_trial,1) : ...
                    curr_trial_list_frames(curr_trial+1,1)-1) = ...
                    curr_trial_list_frames(curr_trial,2);                
            end
            % fill in the end
            trial_frames{curr_xsg}(curr_trial_list_frames(end,1):end) = ...
                curr_trial_list_frames(end,2);            
            
        end
        trial_frames_all{curr_animal,curr_day} = horzcat(trial_frames{:});
        
        lever_velocity = vertcat(lever_velocity_full{:});
        
        lever_movement = vertcat(lever_active_full{:});
        movement_trace = zeros(1,num_frames);
        movement_trace(lever_movement) = 1;
        
        movement_trace_all{curr_animal,curr_day} = movement_trace;
        move_starts = find(diff([0 movement_trace 0]) == 1);
        move_stops =find(diff([0 movement_trace 0]) == -1);
        curr_move_epoch_frames = arrayfun(@(x) move_starts(x): ...
            move_stops(x),1:length(move_starts),'uni',false);
        
        move_epoch_frames{curr_animal,curr_day} = curr_move_epoch_frames;
        rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(reward_frames_full{:}))),curr_move_epoch_frames);
        cued_rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(cued_reward_frames_full{:}))),curr_move_epoch_frames);
        catch_rewarded_movements{curr_animal,curr_day} = ...
            cellfun(@(x) any(intersect(x,vertcat(catch_reward_frames_full{:}))),curr_move_epoch_frames);
        
        cue_frames_all{curr_animal,curr_day} = vertcat(cue_frames_full{:});
        reward_frames_all{curr_animal,curr_day} = vertcat(reward_frames_full{:});
        cued_reward_frames_all{curr_animal,curr_day} = vertcat(cued_reward_frames_full{:});
        catch_reward_frames_all{curr_animal,curr_day} = vertcat(cued_reward_frames_full{:});
        
    end
  
    disp(num2str(curr_animal));
end


