%% PUT IN MEMORY: activity/classification, movements/dispatcher

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

gad_activity_thresh_all = cell(8,15);
gad_class_all = cell(8,15);

pyr_activity_all = cell(8,15);
pyr_class_all = cell(8,15);

movement_trace_all = cell(8,15);
move_epoch_frames = cell(8,15);
rewarded_movements = cell(8,15);
cued_rewarded_movements = cell(8,15);
catch_rewarded_movements = cell(8,15);

cue_frames_all = cell(8,15);
reward_frames_all = cell(8,15);
cued_reward_frames_all = cell(8,15);
catch_reward_frames_all = cell(8,15);

trial_frames_all = cell(8,15);

for curr_animal = 1:8
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
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_caEvents_5.mat'])
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    
           
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
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
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
                    rewarded_move_start > cue_sample(curr_trial_idx)-50;
                
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
        end
        
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

  
    clc
    curr_animal
end


%% PUT THIS IN ABOVE BEFORE RUNNING AGAIN: get trial numbers for frames

trial_frames_all = cell(8,15);
for curr_animal = 1:8
    animal = animals{curr_animal};
         
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
    
    sessions = find(cellfun(@(x) ~isempty(x),pyr_class_all(curr_animal,:)));
    
    for curr_day = sessions
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name],'bhv');
            
        %%%% Find and store movements
        num_frames = size(pyr_activity_all{curr_animal,curr_day},2);
        
        % load split image file to get loop sizes        
        loop_frames = unique(bhv.add_frames);
        loop_size = [diff(loop_frames);num_frames-loop_frames(end)];
        framerate = bhv.framerate;
        clear bhv;
        
%         %if loading in image anyway, this is probably safer
%         loop_size_split = cellfun(@length,im.roi_trace_split);
%         loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
%         loop_size_sum = cumsum(loop_size_split);
%         loop_frames = [0;loop_size_sum(8:8:end-7)];
        
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
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
             

        % THE POINT: get which trial each frame belongs to
        trial_frames = cellfun(@(x) nan(1,x),num2cell(loop_size),'uni',false);
        
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
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
        curr_day
    end
    curr_animal
end












%% (now fused with above), put in memory movements/dispatcher

movement_trace_all = cell(8,15);
move_epoch_frames = cell(8,15);
rewarded_movements = cell(8,15);
cued_rewarded_movements = cell(8,15);
catch_rewarded_movements = cell(8,15);

cue_frames_all = cell(8,15);
reward_frames_all = cell(8,15);
cued_reward_frames_all = cell(8,15);
catch_reward_frames_all = cell(8,15);
for curr_animal = 1:8
    
    animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
    animal = animals{curr_animal};
    
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
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
    
    for curr_day = sessions
        day = days{curr_day};
        
        num_frames = size(gad_activity_thresh_all{curr_animal,curr_day},2);
        
        % load split image file to get loop sizes
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name],'bhv');
        
        loop_frames = unique(bhv.add_frames);
        loop_size = [diff(loop_frames);num_frames-loop_frames(end)];
        % WARNING: assuming last loop here is 4000 
        framerate = bhv.framerate;
        clear bhv;
        
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
        
        if isfield(saved_history,'CatchSection_catch_trial')
            catch_trials = cell2mat(saved_history.CatchSection_catch_trial);
            % there's one or more uncompleted trials at the end
            catch_trials = catch_trials(1:length(bhv_frames));
        end
                
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        % not necessary with bhv
%         loop_size_split = cellfun(@length,im.roi_trace_split);
%         loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
%         loop_size_sum = cumsum(loop_size_split);
%         loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        
        lever_active_full = cell(length(xsg_filenames),1);
        lever_velocity_full = cell(length(xsg_filenames),1);
        cue_frames_full =  cell(length(xsg_filenames),1);
        reward_frames_full =  cell(length(xsg_filenames),1);
        cued_reward_frames_full = cell(length(xsg_filenames),1);
        catch_reward_frames_full = cell(length(xsg_filenames),1);
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
                    rewarded_move_start > cue_sample(curr_trial_idx)-50;
                
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
        end
        
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
        
        curr_day
        
    end
    
    animal
end









%% get mean population activity over days for classification

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    load([animal '_classified_cells.mat']);
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    %concat_peaks = [peak_matrix_full{curr_animal}{:}];
    %box_filter = ones(1,10000);
    %concat_peaks_spread = conv2(concat_peaks,box_filter,'same');
    concat_peaks_spread = concat_peaks_10000{curr_animal};
    
    % do all the following for pyr and gad
    for curr_cell_type = 1:1;
        
        if curr_cell_type == 1
            cell_type = 'pyr';
            curr_cells = pyr_cells;
        elseif curr_cell_type == 2;
            cell_type = 'gad';
            curr_cells = gad_cells;
        end
        
        concat_peaks_downsamp = concat_peaks_spread(curr_cells,1:1000:end);
        
        % normalize each row by the max
        concat_peaks_downsamp_norm = bsxfun(@times,concat_peaks_downsamp, ...
            1./max(concat_peaks_downsamp,[],2));
        concat_peaks_downsamp_norm(isnan(concat_peaks_downsamp_norm)) = 0;
        
%         concat_peaks_downsamp_norm = bsxfun(@times,concat_peaks_downsamp, ...
%             1./sqrt(sum(concat_peaks_downsamp.^2,2)));
%         concat_peaks_downsamp_norm(isnan(concat_peaks_downsamp_norm)) = 0;
        
        concat_peaks_downsamp = concat_peaks_downsamp_norm;
        
        day_length = cellfun(@(x) size(x,2)/1000, peak_matrix_full{curr_animal});
        day_length_sum = cumsum(day_length);
        
        if curr_animal == 1;
            concat_peaks_downsamp(:,1:day_length_sum(2)) = 0;
        end
        
        % kmeans sort traces
        kidx = kmeans(concat_peaks_downsamp,10);
        [temp sort_idx] = sort(kidx);
%         figure;imagesc(concat_peaks_downsamp(sort_idx,:));colormap(gray);
%         for i = 1:length(day_length_sum);
%             line([day_length_sum(i) day_length_sum(i)],ylim,'color','b','linewidth',2);
%         end
%         title('kmeans sorted')
        
        % condense activity to days
        day_bounds = [1 day_length_sum];
        concat_peaks_daily = zeros(size(concat_peaks_downsamp,1),length(day_length));
        for i = 1:length(day_bounds)-1
            concat_peaks_daily(:,i) = nanmean(concat_peaks_downsamp( ...
                :,day_bounds(i):day_bounds(i+1)),2);
        end

        % get the max of the mean activity in each sorted group
        [temp kmeans_peak_activity] = max(grpstats(concat_peaks_daily,kidx),[],2);
        mean_kmeans = mean(grpstats(concat_peaks_daily,kidx),2);
        
        temp = nanmax([classified_cells.cue_max_xcorr{:}],[],2);
        a = grpstats(temp(pyr_cells),kidx);
        temp2 = nanmax([classified_cells.reward_max_xcorr{:}],[],2);
        b = grpstats(temp2(pyr_cells),kidx);
        scatter3(kmeans_peak_activity,a/max(a),b/max(b),mean_kmeans*100);
    end
    animal
end

%% Get correlation between movement/quiescence over days

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
all_corr = cell(size(animals));
scoresum.moving = cell(size(animals));
scoresum.quiescent = cell(size(animals));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    load([animal '_classified_cells.mat']);
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    
    curr_corr_all = cellfun(@(x,y) corrcoef(x,y),classified_cells.movement_max_xcorr, ...
        classified_cells.quiescent_max_xcorr,'uni',false);
    %     curr_corr_all = cellfun(@(x,y) corrcoef(x,y),classified_cells.cue_max_xcorr, ...
    %         classified_cells.reward_max_xcorr,'uni',false);
    
    curr_corr = cellfun(@(x) x(2), curr_corr_all);
    all_corr{curr_animal} = nan(1,15);
    all_corr{curr_animal}(1:length(curr_corr)) = curr_corr;
    
    scoresum.moving{curr_animal} = nan(1,15);
    scoresum.moving{curr_animal}(1:length(curr_corr)) = ...
        cellfun(@median,classified_cells.movement_max_xcorr);
    
    scoresum.quiescent{curr_animal} = nan(1,15);
    scoresum.quiescent{curr_animal}(1:length(curr_corr)) = ...
        cellfun(@median,classified_cells.quiescent_max_xcorr);
    
    scoresum.movingcells{curr_animal} =nan(1,15);
    scoresum.movingcells{curr_animal}(1:length(curr_corr)) = cellfun(@(x,y) ...
        mean(x > y), classified_cells.movement_max_xcorr, ...
        classified_cells.quiescent_max_xcorr);
    
    scoresum.quiescentcells{curr_animal} =nan(1,15);
    scoresum.quiescentcells{curr_animal}(1:length(curr_corr)) = cellfun(@(x,y) ...
        mean(x < y), classified_cells.movement_max_xcorr, ...
        classified_cells.quiescent_max_xcorr);
        
    
end
% delete AP71 day 1,2
all_corr{1}(1:2) = NaN;
scoresum.moving{1}(1:2) = NaN;
scoresum.quiescent{1}(1:2) = NaN;

%% get the corrcoef matrix of mean lever presses for days

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
cc_grid = nan(15,15,length(animals));
for curr_animal = 1:length(animals)
    
    lever_filename = [animals{curr_animal} 'lever_params.mat'];
    load(lever_filename);
    
    a = cellfun(@(x,y) [x{y}],grab_var.rewarded_lever_force_fixtime ...
        (1:length(grab_var.cued_movement)), ...
        grab_var.cued_movement,'uni',false);
    b = cell2mat(cellfun(@(x) nanmedian(x,2),a,'uni',false));
    b_cc = corrcoef(b);
    
    cc_grid(1:length(b_cc),1:length(b_cc),curr_animal) = b_cc;
    
end



%% get population correlation matrix for binary classified cells

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
movement_corr = nan(15,15,length(animals));
move_cells = nan(length(animals),15);
quiescent_corr = nan(15,15,length(animals));
still_cells = nan(length(animals),15);

movestill_corr = nan(15,15,length(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    load([animal '_classified_cells_move_epoch_shuffle.mat']);
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    curr_cells = setdiff(pyr_cells,filled_cells);
    
    curr_move = vertcat(classified_cells.move_cells_peak{:});%.* ...
        %horzcat(classified_cells.move_deconv_peak{:})';
%     curr_move_norm = bsxfun(@times,curr_move,1./max(curr_move));
%     curr_move_norm(isnan(curr_move_norm)) = 0;
    curr_move_cc = corrcoef(curr_move(:,curr_cells)');
        %(+curr_move(:,curr_cells)*+curr_move(:,curr_cells)')/length(curr_cells);
        %(%*2./ ...
        %(repmat(sum(curr_move(:,curr_cells),2),1,size(curr_move,1)) + ...
        %repmat(sum(curr_move(:,curr_cells),2)',size(curr_move,1),1));
        
        %corrcoef(curr_move(:,curr_cells)');%...
    movement_corr(1:size(curr_move,1),1:size(curr_move,1),curr_animal) = ...
        curr_move_cc;
    move_cells(curr_animal,1:size(curr_move,1)) = mean(curr_move(:,curr_cells),2);
    
    curr_still = vertcat(classified_cells.still_cells_peak{:});%.* ...
        %horzcat(classified_cells.still_deconv_peak{:})';
    curr_still_cc = (+curr_still(:,curr_cells)*+curr_still(:,curr_cells)')/length(curr_cells);
        %(+curr_still(:,curr_cells)*+curr_still(:,curr_cells)')/length(curr_cells);%*2./ ...
        %(repmat(sum(curr_still(:,curr_cells),2),1,size(curr_still,1)) + ...
        %repmat(sum(curr_still(:,curr_cells),2)',size(curr_still,1),1));
        
        %corrcoef(curr_still(:,curr_cells)')
    quiescent_corr(1:size(curr_still,1),1:size(curr_still,1),curr_animal) = ...
        curr_still_cc;
    still_cells(curr_animal,1:size(curr_still,1)) = mean(curr_still(:,curr_cells),2);
    
    curr_movestill_cc = +curr_move(:,curr_cells)*+curr_still(:,curr_cells)';
    movestill_corr(1:size(curr_move,1),1:size(curr_move,1),curr_animal) = ...
        curr_movestill_cc;
    
end
% % Input zeros for AP71 day 1,2
% movement_corr(1:2,:,1) = NaN;
% movement_corr(:,1:2,1) = NaN;
% quiescent_corr(1:2,:,1) = NaN;
% quiescent_corr(:,1:2,1) = NaN;
% movestill_corr(1:2,:,1) = NaN;
% movestill_corr(:,1:2,1) = NaN;
% move_cells(1,1:2) = NaN;
% still_cells(1,1:2) = NaN;

%% get fractions/overlap of classified cells
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

all_dot = nan(15,15,8);
pyr_class_all = cell(length(animals),15);
gad_class_all = cell(length(animals),15);
all_frac_overlap = nan(15,15,length(animals));

all_pyr_reliability = cell(length(animals),15);
all_gad_reliability = cell(length(animals),15);

load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/trial_activity_caEvents.mat');
%figure;
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    %curr_cells = setdiff(pyr_cells,filled_cells);
    clear curr_move curr_still
    curr_move = vertcat(classified_cells.move_cells_peak{:});
    curr_still = vertcat(classified_cells.still_cells_oopsi{:})';
%     subplot(3,3,curr_animal);imagesc(curr_move(:,curr_cells));colormap(gray)
    if curr_animal == 1
        curr_move = [nan(2,size(curr_move,2));curr_move];
        curr_still = [nan(2,size(curr_move,2));curr_move];
    end
    
%     temp = horzcat(classified_cells.move_deconv_peak{:})';
%     a = temp(:,curr_cells);
%     kidx = kmeans(a',5,'EmptyAction','drop');
%     [aseresaf sort_idx] = sort(kidx);
%     subplot(3,3,curr_animal);imagesc(a(:,sort_idx));colormap(gray)
    
%     curr_frac = sum(curr_move(:,curr_cells),2)/length(curr_cells);
%     curr_dot = (+curr_move(:,curr_cells)*+curr_move(:,curr_cells)')/length(curr_cells);
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    pyr_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    gad_class_all(curr_animal,sessions) = ...
        cellfun(@(x) x(gad_unfilled),classified_cells.move_cells_peak(sessions),'uni',false);
    %all_dot(1:length(curr_frac),1:length(curr_frac),curr_animal) = curr_dot;
%     all_dot(1:length(curr_frac),1:length(curr_frac),curr_animal) = ...
%         corrcoef(curr_move');
    
%     curr_overlap = (+curr_move(:,curr_cells)*+curr_move(:,curr_cells)')/length(curr_cells);
%     %curr_overlap(logical(eye(size(curr_overlap)))) = NaN;
%     all_frac_overlap(1:length(curr_frac),1:length(curr_frac),curr_animal) = curr_overlap;
        
    
    
    % pick measure of trial activity
    curr_measure = trial_activity.peak;
    
%     % get/save reliabilities
%     all_pyr_reliability(curr_animal,sessions) = ...
%         cellfun(@(x,y) nanmean(x(pyr_unfilled & y,:),2), ...
%         curr_measure{curr_animal}.cued.moving(sessions),classified_cells.move_cells_peak(sessions),'uni',false);
%     all_gad_reliability(curr_animal,sessions) = ...
%         cellfun(@(x,y) nanmean(x(gad_unfilled & y,:),2), ...
%         curr_measure{curr_animal}.cued.moving(sessions),classified_cells.move_cells_peak(sessions),'uni',false);
%

    % get/save reliabilities
    all_pyr_reliability(curr_animal,sessions) = ...
        cellfun(@(x) nanmean(x(pyr_unfilled,:) > 0,2), ...
        curr_measure{curr_animal}.cued.moving(sessions),'uni',false);
    all_gad_reliability(curr_animal,sessions) = ...
        cellfun(@(x) nanmean(x(gad_unfilled,:) > 0,2), ...
        curr_measure{curr_animal}.cued.moving(sessions),'uni',false);

end



%% get trial-by-trial pyr/gad correlation


animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

daily_pyr_gad_corr = cell(length(animals),15);
daily_pyr_gad_corr_slope = cell(size(animals));

all_pyr_gad_corr = cell(size(animals));
all_pyr_gad_corr_slope = nan(size(animals));

total_concat_pyr_act = [];
total_concat_gad_act = [];

for i = 1:8
    animal = animals{i};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
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
    
    % pick max/peak as measure of activity
    curr_measure = trial_activity.peak{i};
    
    % for AP71 exclude first two days
    if i == 1;
       curr_measure.cued.moving{1} = nan(size(curr_measure.cued.moving{1}));
       curr_measure.cued.moving{2} = nan(size(curr_measure.cued.moving{2}));
    end
    
    pyr_unfilled = setdiff(pyr_cells,filled_cells);
    pyr_act = cellfun(@(x) +[x(pyr_unfilled,:) > 0],curr_measure.cued.moving,'uni',false);
    pyr_act_cat = horzcat(pyr_act{:});
    
    gad_unfilled = setdiff(gad_cells,filled_cells);
    gad_act = cellfun(@(x) +[x(gad_unfilled,:) > 0],curr_measure.cued.moving,'uni',false);
    gad_act_cat = horzcat(gad_act{:});
         
    % prepare bins for fractions of cells active
    binedge = [0:0.01:1-0.01 Inf];
        
    day_start = 1;
    if i == 1;
        day_start = 3;
    end
    
    sessions = day_start:length(pyr_act);
    
    daily_pyr_gad_corr_slope{i} = nan(size(pyr_act));
    for curr_day = day_start:length(pyr_act);
        clear day_frac_activity_means day_frac_activity_sems
        curr_pyr_act = sum(pyr_act{curr_day})/length(pyr_unfilled);
        curr_gad_act = sum(gad_act{curr_day})/length(gad_unfilled);
        
        % get rid of nans
        curr_pyr_act(isnan(curr_pyr_act)) = [];
        curr_gad_act(isnan(curr_gad_act)) = [];
        
        [day_frac_n day_frac_activity_grp] = histc(curr_pyr_act,binedge);
        day_frac_activity_means = grpstats(curr_gad_act,day_frac_activity_grp);
        day_frac_activity_sems = grpstats(curr_gad_act,day_frac_activity_grp,'sem');
        
        day_frac_activity_x_idx = sort(unique(day_frac_activity_grp));
        day_frac_activity_x = binedge(day_frac_activity_x_idx+1);
        daily_pyr_gad_corr{i,curr_day} = nan(length(binedge),3);
        daily_pyr_gad_corr{i,curr_day}(unique(day_frac_activity_grp),1) = day_frac_activity_x;
        daily_pyr_gad_corr{i,curr_day}(unique(day_frac_activity_grp),2) = day_frac_activity_means;
        daily_pyr_gad_corr{i,curr_day}(unique(day_frac_activity_grp),3) = day_frac_activity_sems;
        
        p = polyfit(day_frac_activity_x', ...
            day_frac_activity_means,1);
        
        daily_pyr_gad_corr_slope{i}(curr_day) = p(1);
        
    end
    
    concat_pyr_act = sum(pyr_act_cat)/length(pyr_unfilled);
    concat_gad_act = sum(gad_act_cat)/length(gad_unfilled);  
    
    % get rid of nans
    concat_pyr_act(isnan(concat_pyr_act)) = [];
    concat_gad_act(isnan(concat_gad_act)) = [];
    
    [all_frac_n all_frac_activity_grp] = histc(concat_pyr_act,binedge);
    all_frac_activity_means = grpstats(concat_gad_act,all_frac_activity_grp);
    all_frac_activity_sems = grpstats(concat_gad_act,all_frac_activity_grp,'sem');
    
    all_frac_activity_x_idx = sort(unique(all_frac_activity_grp));
    all_frac_activity_x = binedge(all_frac_activity_x_idx+1);
    
    p = polyfit(all_frac_activity_x', ...
            all_frac_activity_means,1);
    
    all_pyr_gad_corr{i}(:,1) = all_frac_activity_x;
    all_pyr_gad_corr{i}(:,2) = all_frac_activity_means;
    all_pyr_gad_corr{i}(:,3) = all_frac_activity_sems;
    all_pyr_gad_corr_slope(i) = p(1);
    
    total_concat_pyr_act = [total_concat_pyr_act concat_pyr_act];
    total_concat_gad_act = [total_concat_gad_act concat_gad_act];
    
    
    % get distance-activity relationship
    pyr_act_cat(isnan(pyr_act_cat)) = 0;
    gad_act_cat(isnan(gad_act_cat)) = 0;
    
    pyr_gad_dist = roi_dist(pyr_unfilled,gad_unfilled); 
    distance_bins = [0:5:500];
    pyr_gad_dist_act = zeros(length(distance_bins)-1,size(pyr_act_cat,2));
    dist_gad_available = zeros(size(distance_bins));
    for curr_distance = 1:length(distance_bins)
        % find gad cells eligible for this (within distance bin of pyr)
        curr_pyr_gad_dist = pyr_gad_dist < distance_bins(curr_distance);
        
        % pull out trial-by-trial which gad cells are eligible
        curr_gad_idx = cell2mat(arrayfun(@(x) ...
            any(curr_pyr_gad_dist(logical(pyr_act_cat(:,x)),:),1)', ...
            [1:size(pyr_act_cat,2)],'uni',false));
        
        % get activity of eligible cells for each trial
        curr_gad_act = zeros(size(gad_act_cat));
        curr_gad_act(curr_gad_idx) = gad_act_cat(curr_gad_idx);
        
        % save fractions active
        pyr_gad_dist_act(curr_distance,:) = nansum(curr_gad_act)/length(gad_unfilled);
        
        % get what fraction of gad cells are available
        dist_gad_available(curr_distance) = sum(any(curr_pyr_gad_dist))/length(gad_unfilled);
    end
    pyr_gad_dist_diff = bsxfun(@minus,nansum(gad_act_cat)/length(gad_unfilled),pyr_gad_dist_act);
    
end

% get total combined
binedge = [0:0.03:1];
[all_frac_n all_frac_activity_grp] = histc(total_concat_pyr_act,binedge);
all_frac_activity_means = grpstats(total_concat_gad_act,all_frac_activity_grp);
all_frac_activity_sems = grpstats(total_concat_gad_act,all_frac_activity_grp,'sem');

all_frac_activity_x_idx = sort(unique(all_frac_activity_grp));
all_frac_activity_x = binedge(all_frac_activity_x_idx+1);
figure;errorbar(all_frac_activity_x,all_frac_activity_means,all_frac_activity_sems,'k','linewidth',2)

groups = {[1 2] [3 4 5 6 7 8] [9 10 11 12 13 14]};
col = copper(length(groups));
figure; hold on;
for i = 1:length(groups)
    grp = groups{i};
    a = horzcat(daily_pyr_gad_corr{:,grp});
    a = a(:,2:3:end)';
    good_bins = find(any(~isnan(a)));
    a_ci = bootci(1000,@nanmean,a(:,good_bins));
    plot(good_bins,nanmean(a(:,good_bins)),'color',col(i,:));    
    plot(good_bins,a_ci(1,:),'color',col(i,:),'linestyle','--');
    plot(good_bins,a_ci(2,:),'color',col(i,:),'linestyle','--');
end

%% Correlate movement-related cells with other movement-related cells
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

sum_move_day_bins = [0:2:14];

all_sum_corr = cell(length(sum_move_day_bins)-1,15,length(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    unfilled_pyr = setdiff(pyr_cells,filled_cells);
    unfilled_gad = setdiff(gad_cells,filled_cells);
    
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    [n move_day_bin] = histc(sum_move_day,sum_move_day_bins);
    
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
    
    start_day = 1;
    if strmatch(animal,'AP71')
        start_day = 3;
    end
    
    for curr_day = start_day:length(days)
        
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        curr_cells = intersect(find(classified_cells.move_cells_peak{curr_day}),unfilled_pyr);
        
%         curr_peaks = AP_trace_peak(im.roi_trace_df(curr_cells,:));
%         curr_peak_matrix = zeros(size(im.roi_trace_df(curr_cells,:)));
%         for i = 1:length(curr_cells)
%             curr_peak_matrix(i,curr_peaks{i}) = 1;
%         end
        
        curr_peak_matrix = im.roi_concat_oopsi(curr_cells,:);
        
        spread_frames = 5;
        spread_frames_filter = ones(1,spread_frames);
        curr_peak_matrix_spread = conv2(curr_peak_matrix, ...
            spread_frames_filter,'same');
        
        peak_matrix_corrcoef = corrcoef(curr_peak_matrix_spread');
        
        %     kidx = kmeans(curr_peak_matrix_spread,5,'distance','correlation');
        %     [asdf sort_idx] = sort(kidx);
        %     a = corrcoef(curr_peak_matrix_spread(sort_idx,:)');
        %     figure;imagesc(a);colormap(hot)
        
        for i = 1:length(sum_move_day_bins)-1
            curr_bin_cells = move_day_bin(curr_cells) == i;
            tril_matrix = tril(true(nnz(curr_bin_cells)),-1);
            curr_peak_matrix_corrcoef = peak_matrix_corrcoef(curr_bin_cells, ...
                curr_bin_cells);
            all_sum_corr{i,curr_day,curr_animal} = ...
                curr_peak_matrix_corrcoef(tril_matrix);
        end
        
        disp(curr_day)
    end
    disp('end')
    
    
end

cat_data_mean = nan(size(all_sum_corr,1),size(all_sum_corr,2));
cat_data_sem = nan(size(all_sum_corr,1),size(all_sum_corr,2));
for curr_day = 1:size(all_sum_corr,2);
    for curr_bin = 1:size(all_sum_corr,1);
        curr_data = vertcat(all_sum_corr{curr_bin,curr_day,:});
        cat_data_mean(curr_bin,curr_day) = nanmean(curr_data);
        cat_data_sem(curr_bin,curr_day) = nanstd(curr_data)/sqrt(length(curr_data));
    end
end


%% Get reliability of cells dependent on classification days
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

all_reliability = cell(length(animals),15);
all_sum_move = cell(length(animals),15);

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);
% gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%    gad_activity_thresh_all,move_epoch_frames,'uni',false);
% pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%    pyr_activity_all,move_epoch_frames,'uni',false);

all_reliability = cellfun(@(x) nanmean(x,2),pyr_move_active,'uni',false);
all_sum_move = cell(8,1);
all_reliability_mean = cell(8,1);
for i = 1:8
    curr_move = vertcat(pyr_class_all{i,:})';
    curr_reliability = horzcat(all_reliability{i,:});
    all_sum_move{i} = nanmean(curr_move,2);
    
    move_reliability = curr_reliability;
    move_reliability(~curr_move) = NaN;
    all_reliability_mean{i} = nanmean(move_reliability,2);
end

% concatenate multiple animals
cat_reliability = cell(1,15);
cat_sum_move = cell(1,15);
for i = 1:15
   cat_reliability{i} = vertcat(all_reliability{:,i}); 
   %cat_sum_move{i} = vertcat(all_sum_move{:,i});   
end

% get groupstats for concatenated animals
a = nan(15,15);
b = nan(15,15);
for i = 1:15
    curr_grp = ceil(cat_sum_move{i}/1);
    curr_days = length(unique(curr_grp));
    a(i,1:curr_days) = grpstats(cat_reliability{i},curr_grp);
    b(i,1:curr_days) = grpstats(cat_reliability{i},curr_grp,'sem');
end

% get groupstats for individual animals
a = nan(15,15,8);
b = nan(15,15,8);
for j = 1:8
    for i = 1:15
        curr_grp = ceil(all_sum_move{j,i}/1);
        curr_days = unique(curr_grp);
        if length(curr_days) == 0
            continue
        end
        a(i,curr_days,j) = grpstats(all_reliability{j,i},curr_grp);
        b(i,curr_days,j) = grpstats(all_reliability{j,i},curr_grp,'sem');
    end
end

% combine days
cat_cat_reliability = cell(1,7);
cat_cat_sum_move = cell(1,7);
for i = 1:7
    cat_cat_reliability{i} = vertcat(cat_reliability{(i-1)*2+1:(i*2)});
    cat_cat_sum_move{i} = vertcat(cat_sum_move{(i-1)*2+1:(i*2)});
end

% split up into days 1:2, 3:14, combine move sum days, errorbar
a = vertcat(cat_reliability{1:2});
a_sum = vertcat(cat_sum_move{1:2});

b = vertcat(cat_reliability{3:14});
b_sum = vertcat(cat_sum_move{3:14});

a(a_sum == 15) = [];
a_sum(a_sum == 15) = [];
b(b_sum == 15) = [];
b_sum(b_sum == 15) = [];

a_grp = ceil(a_sum/4);
a_mean = grpstats(a,a_grp);
a_sem = grpstats(a,a_grp,'sem');

b_grp = ceil(b_sum/4);
b_mean = grpstats(b,b_grp);
b_sem = grpstats(b,b_grp,'sem');

% get mean/sem reliability over days
a = cellfun(@mean,cat_reliability);
b = cellfun(@(x) std(x)/length(x),cat_reliability);
figure
for i = 1:8
   subplot(3,3,i)
   a = cellfun(@mean,all_reliability(i,:));
   b = cellfun(@(x) std(x)/length(x),all_reliability(i,:));
   errorbar(a,b,'k')
end

% scatter plot all reliabilities for move cells in a session, color by days
% classified as movement-related - by animal
figure;
for j = 1:8
    subplot(3,3,j); hold on
    for i = 1:15
        scatter(i+0.5+0.3*rand(size(all_reliability{j,i})),all_reliability{j,i}, ...
            20,all_sum_move{j,i},'filled')
    end
end



%% Dig into reliability: why dip on day 3, cell increase v decrease reliability, etc
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

sum_move_day_all = [];
under_7_all = cell(8,15);
over_7_all = cell(8,15);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    unfilled_pyr = setdiff(pyr_cells,filled_cells);
    unfilled_gad = setdiff(gad_cells,filled_cells);

    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    
    sum_move_day_all = [sum_move_day_all sum_move_day(unfilled_pyr)];

    curr_reliability = cellfun(@(x) nanmean(x,2), ...
            cell_trial_peak_full{curr_animal}.cued.moving,'uni',false);
    
    use_day = 1:4;
    if max(use_day) > length(curr_reliability)
        continue
    end
    
    use_move_cells = any(vertcat(classified_cells.move_cells_peak{use_day}),1);
    use_reliability = horzcat(curr_reliability{use_day});
    
    under_7 = sum_move_day(unfilled_pyr) > 1 & ...
        sum_move_day(unfilled_pyr) <= 7 & ...
        any(use_reliability(unfilled_pyr,:) > 0.2,2)'; %& ...
        %use_move_cells(unfilled_pyr); 
    over_7 = sum_move_day(unfilled_pyr) > 7 & ...
        any(use_reliability(unfilled_pyr,:) > 0.2,2)'; %& ...
        %use_move_cells(unfilled_pyr); 
    
%     use_day = 3;
%     under_7 = sum_move_day(unfilled_pyr) > 1 & ...
%         sum_move_day(unfilled_pyr) <= 7 & ...
%         curr_reliability{use_day}(unfilled_pyr)' < 0.2; 
%     over_7 = sum_move_day(unfilled_pyr) > 7 & ...
%         curr_reliability{use_day}(unfilled_pyr)' < 0.2; 
        
    under_7_all(curr_animal,1:length(curr_reliability)) = ...
        cellfun(@(x) x(unfilled_pyr(under_7)),curr_reliability,'uni',false);
    over_7_all(curr_animal,1:length(curr_reliability)) = ...
        cellfun(@(x) x(unfilled_pyr(over_7)),curr_reliability,'uni',false);
end

% concatenate multiple animals
cat_under_7 = cell(1,15);
cat_over_7 = cell(1,15);
for i = 1:15
   cat_under_7{i} = vertcat(under_7_all{:,i}); 
   cat_over_7{i} = vertcat(over_7_all{:,i}); 
end
figure; hold on
errorbar(cellfun(@nanmean,cat_under_7), ...
    cellfun(@(x) std(x)/sqrt(length(x)),cat_under_7),'r','linewidth',2);
errorbar(cellfun(@nanmean,cat_over_7), ...
    cellfun(@(x) std(x)/sqrt(length(x)),cat_over_7),'k','linewidth',2);

% concatenate multiple animals
cat_reliability = cell(1,15);
cat_sum_move = cell(1,15);
for i = 1:15
   cat_reliability{i} = vertcat(all_reliability{:,i}); 
   cat_sum_move{i} = vertcat(all_sum_move{:,i});   
end

% scatter plot all reliabilities for move cells in a session, color by days
% classified as movement-related
figure;hold on
for i = 1:15
    scatter(i+0.3*(rand(size(cat_reliability{i}))-0.5),cat_reliability{i}, ...
        20,cat_sum_move{i},'filled')
end

% show reliability trajectories for cells classified on a particular day
start_reliability = cell(size(animals));
end_reliability = cell(size(animals));
sum_move_cells = cell(size(animals));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    unfilled_pyr = setdiff(pyr_cells,filled_cells);
    unfilled_gad = setdiff(gad_cells,filled_cells);
    
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    vertcat_classified = vertcat(classified_cells.move_cells_peak{:});
    
    curr_cells = intersect(unfilled_pyr,find(sum_move_day > 2)); 
    
    start_reliability{curr_animal} = nan(size(curr_cells));
    end_reliability{curr_animal} = nan(size(curr_cells));
    
    for i = 1:length(curr_cells)
       start_day = find(vertcat_classified(:,curr_cells(i)),1);
       end_day = find(vertcat_classified(:,curr_cells(i)),1,'last');
       
       if curr_animal == 1;
           start_day = start_day + 2;
           end_day = end_day + 2;
       end
       
       start_reliability{curr_animal}(i) = ...
           nanmean(cell_trial_peak_full{curr_animal}.cued.moving{start_day} ...
           (curr_cells(i),:),2);
       end_reliability{curr_animal}(i) = ...
           nanmean(cell_trial_peak_full{curr_animal}.cued.moving{end_day} ...
           (curr_cells(i),:),2);
    end
    sum_move_cells{curr_animal} = sum_move_day(curr_cells);
end


figure; hold on
for i = 1:8
    scatter(start_reliability{i},end_reliability{i}, ...
        50,sum_move_cells{i},'filled');
end




% make # days classified vs. day classified plot
all_move_days = cell(size(animals));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels))';
    
    unfilled_pyr = setdiff(pyr_cells,filled_cells);
    unfilled_gad = setdiff(gad_cells,filled_cells);
    
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    vertcat_classified = vertcat(classified_cells.move_cells_peak{:});
    
    all_move_days{curr_animal} = arrayfun(@(x) ...
        find(vertcat_classified(:,x)), ...
        unfilled_pyr,'uni',false);
end

figure; hold on
for j = 1:8
    for i = 1:length(all_move_days{j})
        if isempty(all_move_days{j}{i})
            continue
        end
        
        y = all_move_days{j}{i};
        x = repmat(length(y),[length(y) 1]);
        
        plot(x+0.3*(rand(size(x))-0.5),y+0.3*(rand(size(x))-0.5),'.k');
    end
end
xlabel('# days classified');ylabel('day classified');


move_cell_grid = zeros(15,15,8);
for j = 1:8
    for i = 1:length(all_move_days{j})
        if isempty(all_move_days{j}{i})
            continue
        end
        
        y = all_move_days{j}{i};
        x = repmat(length(y),[length(y) 1]);
        
        move_cell_grid(y,x,j) = move_cell_grid(y,x,j) + 1;
    end
end

move_grid_split = cell2mat(reshape(arrayfun(@(x) ...
    [sum(move_cell_grid(1:14,1:7,x),2) ...
    sum(move_cell_grid(1:14,8:14,x),2)],1:8,'uni',false),[1,1,8]));

move_grid_split_norm = bsxfun(@times,move_grid_split, ...
    1./sum(move_grid_split));

a = sum(move_cell_grid,3);
a_norm = bsxfun(@times,a,1./sum(a,1));
figure;imagesc(a_norm);colormap(hot);
xlabel('# days classified');ylabel('day classified');

a_sum = [sum(a(1:14,1:7),2) sum(a(1:14,8:14),2)];
a_sum_norm = bsxfun(@times,a_sum,1./sum(a_sum));



%% Get trial correlations over time
% Mean correlation/day correlation by max df

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

% load trial-by-trial activity
load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/cell_trial_all_activity_rewmove.mat');

all_mean_correlation = nan(15,15,8);
all_day_correlation = nan(8,15);
all_day_dot = nan(8,15);
concat_day_correlation = cell(8,15);
all_day_correlation_full = cell(8,15);
all_cell_omit_corr_diff = cell(8,15);
all_sum_move_day = cell(8,1);
all_measure_cells = cell(8,15);

all_move_cells = cell(8,15);

all_activity_transient = nan(8,15);
all_activity_persistent = nan(8,15);

all_activity_transient_std = nan(8,15);
all_activity_persistent_std = nan(8,15);

all_cells = cell(8,1);

all_reliability = cell(8,15);

for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';

    % get usable cells
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    all_cells{curr_animal} = pyr_unfilled;
    
    % get days classified
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_caEvents_5.mat'])
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    all_sum_move_day{curr_animal} = sum_move_day(pyr_unfilled);
    
    % pick measure of trial activity
    curr_measure = trial_activity.peak;
% 
%     curr_measure_cells = cellfun(@(x,y) max(...
%         cat(3,x(pyr_unfilled,:),y(pyr_unfilled,:)),[],3), ...
%         curr_measure{curr_animal}.cued.moving, ...
%         curr_measure{curr_animal}.cued.quiescent,'uni',false);
    
    if strmatch(animal,'AP71')
        classified_cells.move_cells_peak{1} = false(1,length(roi_labels));
        classified_cells.move_cells_peak{2} = false(1,length(roi_labels));
    end
    
    persistent_cells = sum_move_day > 7;
    transient_cells = sum_move_day > 1 & sum_move_day <= 7;
    other_cells = ~persistent_cells & ~transient_cells;
    
    %curr_measure_cells = cellfun(@(x,y) x(pyr_unfilled & y,:), ...
    %    curr_measure{curr_animal}.cued.moving,classified_cells.move_cells_peak,'uni',false);
    curr_measure_cells = cellfun(@(x) x(pyr_unfilled,:), ...
        curr_measure{curr_animal}.cued.moving,'uni',false);
    
    all_move_cells(curr_animal,1:length(classified_cells.move_cells_peak)) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak,'uni',false);
    
    % get/save reliabilities
    all_reliability(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = ...
        cellfun(@(x) nanmean(x,2), curr_measure_cells,'uni',false);
    if strmatch(animal,'AP71')
        all_reliability{1,1} = [];
        all_reliability{1,2} = [];
    end
    
    for i = 1:length(curr_measure_cells)
        nan_trials = all(isnan(curr_measure_cells{i}));
        curr_measure_cells{i}(:,nan_trials) = [];
        curr_measure_cells{i}(isnan(curr_measure_cells{i})) = 0;
    end
    
    all_measure_cells(curr_animal,1:length(curr_measure_cells)) = curr_measure_cells;
    
%     all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
%         nanmean(nansum(x(pyr_unfilled,:)/sum(pyr_unfilled))), ...
%         curr_measure{curr_animal}.cued.moving);
    % I'm kind of butchering this variable to get stuff out temporarily
    all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nansum(any(x(pyr_unfilled&transient_cells,:),2))/sum(pyr_unfilled&transient_cells), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanmean(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    all_activity_transient_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & transient_cells,:)/sum(pyr_unfilled & transient_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    curr_measure_cells_corrcoef = cellfun(@corrcoef, ...
        curr_measure_cells,'uni',false);
    curr_measure_cells_corrcoef_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_corrcoef); 
    
    curr_measure_bin = cellfun(@(x) round(x - 0.5), ...
        curr_measure_cells,'uni',false);   
    curr_measure_cells_dot = cellfun(@(x) (x'*x)/size(x,1), ...
        curr_measure_bin,'uni',false);
    curr_measure_cells_dot_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_dot);
    
    start_day = 1;
    if curr_animal == 1;
        start_day = 3;
    end
    
    all_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef_mean)) = ...
        curr_measure_cells_corrcoef_mean(start_day:end);
    all_day_dot(curr_animal, ...
        start_day:length(curr_measure_cells_dot_mean)) = ...
        curr_measure_cells_dot_mean(start_day:end);
    
%     curr_measure_cells_mean = cell2mat(cellfun(@(x) nanmean(x,2), ...
%         curr_measure_cells,'uni',false));
%     curr_mean_corrcoef = corrcoef(curr_measure_cells_mean);
%     all_mean_correlation(start_day:size(curr_measure_cells_mean,2), ...
%         start_day:size(curr_measure_cells_mean,2),curr_animal) = ...
%         curr_mean_corrcoef(start_day:end,start_day:end);
    
    % correlate mean of trials with all trials
    curr_measure_cells_trialmean_corrcoef = cellfun(@(x) ...
        corrcoef([nanmean(x,2) x]), ...
        curr_measure_cells,'uni',false);
%     concat_day_correlation(curr_animal, ...
%         start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
%         x(2:end,1), curr_measure_cells_trialmean_corrcoef(start_day:end),'uni',false);
    concat_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
        x(tril(true(size(x)),-1)), curr_measure_cells_corrcoef(start_day:end),'uni',false);

    all_day_correlation_full(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef)) = ...
        curr_measure_cells_corrcoef(start_day:end);
    
    % get contribution of each cell to mean correlation of the day
    cell_omit_corr = cellfun(@(x) jackknife(@(y) corrcoef(y),x),...
        curr_measure_cells,'uni',false);
    tril_indx = cellfun(@(x) find(tril(true(size(x,2)),-1)),curr_measure_cells,'uni',false);    
    cell_omit_corr_diff = cellfun(@(x,y,z) z-nanmean(x(:,y),2), ...
        cell_omit_corr,tril_indx,num2cell(curr_measure_cells_corrcoef_mean),'uni',false);
 
    all_cell_omit_corr_diff(curr_animal, ...
        start_day:length(curr_measure_cells)) = cell_omit_corr_diff(start_day:end);
end

%day_combine = num2cell(1:14);
day_combine = {[1 2],[3 4],[5 6],[7:10],[11:14]};
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = temp(:,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'color',[0 0.8 0],'linewidth',2);

% treat each trial-trial comparison as an n
trial_correlation = cell(1,15);
for i = 1:15
    trial_correlation{i} = vertcat(concat_day_correlation{:,i});
end
figure;errorbar(cellfun(@nanmean,trial_correlation), ...
    cellfun(@(x) nanstd(x)./sqrt(length(x)),trial_correlation),'k','linewidth',2);

% compare omitted corr diff with sum move days
concat_mean_omit = [];
concat_sum_move_day = horzcat(all_sum_move_day{:})';
for i = 1:8
   curr_mean_omit = max(horzcat(all_cell_omit_corr_diff{i,7:14}),[],2);
   concat_mean_omit = [concat_mean_omit;curr_mean_omit];
end
figure;hold on;
plot(0.3*(rand(size(concat_sum_move_day))-0.5)+concat_sum_move_day,concat_mean_omit,'.k');
binary_move_day = concat_sum_move_day;
binary_move_day(binary_move_day > 1 & binary_move_day <= 7) = 2;
binary_move_day(binary_move_day > 7) = 3;
a = grpstats(concat_mean_omit,binary_move_day);
b = grpstats(concat_mean_omit,binary_move_day,'sem');
figure;errorbar([2:14],a(3:end-1),b(3:end-1),'k','linewidth',2);

all_transient_cells = cellfun(@(x) x > 1 & x <= 7,all_sum_move_day,'uni',false);
all_persistent_cells = cellfun(@(x) x > 7,all_sum_move_day,'uni',false);
transient_contribution = cell(1,15);
persistent_contribution = cell(1,15);
for i = 1:15;
    curr_data = vertcat(all_cell_omit_corr_diff{:,i});
    curr_transient = horzcat(all_transient_cells{cellfun(@(x) ~isempty(x),all_cell_omit_corr_diff(:,i))})';
    curr_persistent = horzcat(all_persistent_cells{cellfun(@(x) ~isempty(x),all_cell_omit_corr_diff(:,i))})';

    transient_contribution{i} = curr_data(curr_transient);
    persistent_contribution{i} = curr_data(curr_persistent);
end

cat_reliability = cell(1,14);
for i = 1:14
   cat_reliability{i} = vertcat(all_reliability{:,i});    
end

% concat_concat_day_correlation = cell(1,15);
% for i = 1:15
%     concat_concat_day_correlation{i} = vertcat(concat_day_correlation{:,i});
% end
% concat_day_correlation_mean = cellfun(@nanmean,concat_concat_day_correlation);
% concat_day_correlation_sem = cellfun(@(x) std(x)/sqrt(length(x)),concat_concat_day_correlation);
% figure;errorbar(concat_day_correlation_mean,concat_day_correlation_sem,'b','linewidth',2);

% for j = 1:8
%    a = figure; set(a,'Name',num2str(j));
%    for i = 1:length(all_day_correlation_full{j});
%       n = length(all_day_correlation_full{j}); 
%       subplot(ceil(sqrt(n)),ceil(sqrt(n)),i);
%       
%       x = all_day_correlation_full{j}{i};
% %       x(isnan(x)) = 0;
% %       %x(all(~x),:) = [];
% %       
% %       kidx = kmeans(x,2);
% %       [asdf sort_idx] = sort(kidx);
%       
%       imagesc(x);
%       caxis([0 1]);
%    end  
% end
figure;
for i = 1:8
    for j = 1:15
        subplot(8,15,j+15*(i-1))
        hist(all_reliability{i,j},20);
    end
end

bin_edges = 0:0.1:1;
a = cellfun(@(x) histc(x,bin_edges)/length(x),all_reliability,'uni',false);
figure;
colors = jet(15);
for i = 1:8
    subplot(3,3,i);hold on;
    for j = 1:15
    	plot(cumsum(a{i,j}),'color',colors(j,:));
    end
end

figure; hold on
temp_data = nan(length(bin_edges),14);
for i = 1:14
    curr_data = nanmean(cumsum(horzcat(a{:,i})),2);
    plot(curr_data,'color',colors(i,:));    
    
    temp_data(:,i) = curr_data;
end

figure; 
colors = jet(14);
for i = 1:8
    subplot(3,3,i);hold on;
   for j = 1:14;
       plot(all_reliability{i,j},all_cell_omit_corr_diff{i,j},'.','color',colors(j,:));
    end
end

% find pairs of negative correlation trials
[i_c j_c] = cellfun(@(x) find(x < 0), all_day_correlation_full,'uni',false);
ij_unique = cellfun(@(x,y) ...
    unique(sort([x y],2),'rows'),i_c,j_c,'uni',false);
ij_dot = cellfun(@(x,y) diag(x(:,y(:,1))'*x(:,y(:,2))),all_measure_cells,ij_unique);

ij_interleave = cellfun(@(x,y) ...
    reshape(unique(sort([x y],2),'rows')',[],1),i_c,j_c,'uni',false);
paired_neg_trials = cellfun(@(x,y) x(:,y),all_measure_cells,ij_interleave,'uni',false);
paired_neg_dot = cellfun(@(x) diag(x(:,1:2:end)'*x(:,2:2:end)), ...
    paired_neg_trials,'uni',false);
paired_neg_std = cellfun(@(x) [nanstd(x(:,1:2:end)).*nanstd(x(:,2:2:end))]', ...
    paired_neg_trials,'uni',false);

bin_edge = [{[0:0.1:0.5 1]} {[0:0.1:0.5 1]}];
paired_neg_mean = cellfun(@(x) [nanmean(x(:,1:2:end));nanmean(x(:,2:2:end))]', ...
    paired_neg_trials,'uni',false);
paired_neg_mean_hist = cellfun(@(x,y) hist3(x,'edges',bin_edge)/numel(y), ...
    paired_neg_mean,all_day_correlation_full,'uni',false);



% find and plot groups of negative-correlation trials
h_corr = figure;
h_act = figure;
h_hist = figure;
all_reliability_hist = cell(8,15);
temp_animal = 5;
for i = 1:length(all_day_correlation_full{temp_animal});
    if temp_animal == 1 && (i == 1 || i == 2)
        continue
    end
    b = all_day_correlation_full{temp_animal}{i};
    b = b < 0;
    %b(b > 0) = 0;
    %(isnan(b)) = 0;
    try
        kidx = kmeans(b,2);
        [asdf sort_idx] = sort(kidx);
    catch me
        kidx = ones(length(b),1);
        sort_idx = 1:length(b);
    end
    
    % figure out which kmeans group is the low-correlation one
    [asdf pos_corr_kidx] = min(grpstats(nanmean(b,2),kidx));
    pos_kidx = find(kidx == pos_corr_kidx);
    neg_kidx = find(kidx ~= pos_corr_kidx);
    
    figure(h_corr);
    subplot(4,4,i);
    imagesc(b(sort_idx,sort_idx));colormap(gray);
    colormap(summer);
    
    % !!TEMP! plot only high corr group, sort by cells active
    
    [asdf sort_trials] = sortrows([kidx nanmean(all_measure_cells{temp_animal,i})'],[1 2]);
    [asdf sort_cells] = sort(nanmean(all_measure_cells{temp_animal,i}(:,pos_kidx),2));
    figure(h_act);
    subplot(4,4,i);
    %imagesc(all_measure_cells{temp_animal,i}(:,sort_idx));
    imagesc(all_measure_cells{temp_animal,i}(sort_cells,sort_trials));
    colormap(gray); caxis([0 1])
    line([sum(kidx == 1)+0.5 sum(kidx == 1)+0.5],ylim,'color','r');
    
    figure(h_hist);
    subplot(4,4,i);
    bin_edge = 0:0.1:1;
    [n_temp bins] = histc(nanmean(all_measure_cells{temp_animal,i},2),bin_edge);
    bar(bin_edge,n_temp,'k');
%     [f xi] = ksdensity(all_day_correlation_full{temp_animal}{i}(:)); 
%     plot(xi,f,'k');
    xlim([-0.1 1.3])
        
end

all_reliability_hist = cell(8,15);
for animal = 1:8
    for i = find(cellfun(@(x) ~isempty(x),all_measure_cells(animal,:)));
    if animal == 1 && (i == 1 || i == 2)
        continue
    end
    bin_edge = [0 0.001:0.1:1 1];
    [n_temp bins] = histc(nanmean(all_measure_cells{animal,i},2),bin_edge);   
    all_reliability_hist{animal,i} = n_temp/size(all_measure_cells{animal,i},1);   
    end
end
cat_hist = nan(12,14);
for i = 1:14
    cat_hist(:,i) = nanmean(horzcat(all_reliability_hist{:,i}),2);
end

% pull out trials which aren't negative-correlation driving
pos_concat_day_correlation = cell(size(concat_day_correlation));
pos_reliability = cell(size(concat_day_correlation));

frac_active_normal = nan(size(concat_day_correlation));
frac_active_abnormal = nan(size(concat_day_correlation));
frac_active_all = cell(size(concat_day_correlation));

for temp_animal = 1:8;
    for i = 1:length(all_day_correlation_full{temp_animal});
        if temp_animal == 1 && (i == 1 || i == 2)
            continue
        end
        b = all_day_correlation_full{temp_animal}{i};
        b = b < 0;
        %b(b > 0) = 0;
        %(isnan(b)) = 0;
        try
            kidx = kmeans(b,2);
            [asdf sort_idx] = sort(kidx);
        catch me
            kidx = ones(length(b),1);
            sort_idx = 1:length(b);
        end
        [asdf pos_corr_kidx] = min(grpstats(nanmean(b,2),kidx));
        use_kidx = kidx == pos_corr_kidx;
        corr_sort = all_day_correlation_full{temp_animal}{i}(use_kidx,use_kidx);
        pos_concat_day_correlation{temp_animal,i} = ...
            corr_sort(tril(true(size(corr_sort)),-1));
        
        pos_reliability{temp_animal,i} = ...
            nanmean(all_measure_cells{temp_animal,i}(:,use_kidx),2);
        
        frac_active_normal(temp_animal,i) = ...
            nanmean(nanmean(all_measure_cells{temp_animal,i}(:,use_kidx)));
        frac_active_abnormal(temp_animal,i) = ...
            nanmean(nanmean(all_measure_cells{temp_animal,i}(:,~use_kidx)));
        
        frac_active_all{temp_animal,i} = ...
            nanmean(all_measure_cells{temp_animal,i}(:,~use_kidx));
    end
end

cat_corr = cell(1,14);
cat_corr_pos = cell(1,14);
figure
for i = 1:14
   cat_corr{i} = vertcat(concat_day_correlation{:,i});    
   cat_corr_pos{i} = vertcat(pos_concat_day_correlation{:,i});
   bin_edge = -1:0.05:1;
   [n_corr bins] = histc(cat_corr{i},bin_edge);
   [n_corr_pos bins] = histc(cat_corr_pos{i},bin_edge);
   
   subplot(4,4,i); hold on;
   bar(bin_edge,n_corr,'k');
   bar(bin_edge,n_corr_pos,'r');
   xlim([-1 1])
end



figure; hold on
colors = jet(14);
for i = 1:14
   [f xi] = ksdensity(cat_reliability{i}); 
   plot(xi,f,'color',colors(i,:));
end

temp = cell(1,14);
%figure;
for i = 1:14;
    temp{i} = horzcat(frac_active_all{:,i});
    subplot(4,4,i); hold on
    bin_edge = 0:0.02:1;
    [n_temp bins] = histc(temp{i},bin_edge);
    bar(bin_edge,n_temp,'r');
    xlim([-0.1 0.5]);
end

%% Get reaction times of trials which give many negative correlations

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

% load trial-by-trial activity
load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/cell_trial_all_activity_rewmove.mat');

all_mean_correlation = nan(15,15,8);
all_day_correlation = nan(8,15);
all_day_dot = nan(8,15);
concat_day_correlation = cell(8,15);
all_day_correlation_full = cell(8,1);
all_cell_omit_corr_diff = cell(8,15);
all_sum_move_day = cell(8,1);
all_measure_cells = cell(8,15);

all_activity_transient = nan(8,15);
all_activity_persistent = nan(8,15);

all_activity_transient_std = nan(8,15);
all_activity_persistent_std = nan(8,15);

all_reliability = cell(8,15);

for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';

    % get usable cells
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    % get days classified
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    all_sum_move_day{curr_animal} = sum_move_day(pyr_unfilled);
    
    % pick measure of trial activity
    curr_measure = trial_activity.peak;
% 
%     curr_measure_cells = cellfun(@(x,y) max(...
%         cat(3,x(pyr_unfilled,:),y(pyr_unfilled,:)),[],3), ...
%         curr_measure{curr_animal}.cued.moving, ...
%         curr_measure{curr_animal}.cued.quiescent,'uni',false);
    
    if strmatch(animal,'AP71')
        classified_cells.move_cells_peak{1} = false(1,length(roi_labels));
        classified_cells.move_cells_peak{2} = false(1,length(roi_labels));
    end
    
    persistent_cells = sum_move_day > 7;
    transient_cells = sum_move_day > 1 & sum_move_day <= 7;
    other_cells = ~persistent_cells & ~transient_cells;
    
    % get/save reliabilities
     all_reliability(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = ...
         cellfun(@(x,y) nanmean(x(pyr_unfilled & y,:),2), ...
         curr_measure{curr_animal}.cued.moving,classified_cells.move_cells_peak,'uni',false);
%     all_reliability(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = ...
%         cellfun(@(x) nanmean(x(pyr_unfilled,:),2), ...
%         curr_measure{curr_animal}.cued.moving,'uni',false);
    if strmatch(animal,'AP71')
        all_reliability{1,1} = [];
        all_reliability{1,2} = [];
    end

    %curr_measure_cells = cellfun(@(x,y) x(pyr_unfilled & y,:), ...
    %    curr_measure{curr_animal}.cued.moving,classified_cells.move_cells_peak,'uni',false);
    curr_measure_cells = cellfun(@(x) x(pyr_unfilled,:), ...
        curr_measure{curr_animal}.cued.moving,'uni',false);
    
%     for i = 1:length(curr_measure_cells)
%         nan_trials = all(isnan(curr_measure_cells{i}));
%         curr_measure_cells{i}(:,nan_trials) = [];
%         curr_measure_cells{i}(isnan(curr_measure_cells{i})) = 0;
%     end
    
    all_measure_cells(curr_animal,1:length(curr_measure_cells)) = curr_measure_cells;
    
%     all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
%         nanmean(nansum(x(pyr_unfilled,:)/sum(pyr_unfilled))), ...
%         curr_measure{curr_animal}.cued.moving);
    % I'm kind of butchering this variable to get stuff out temporarily
    all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nansum(any(x(pyr_unfilled&transient_cells,:),2))/sum(pyr_unfilled&transient_cells), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanmean(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    all_activity_transient_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & transient_cells,:)/sum(pyr_unfilled & transient_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    curr_measure_cells_corrcoef = cellfun(@corrcoef, ...
        curr_measure_cells,'uni',false);
    curr_measure_cells_corrcoef_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_corrcoef); 
    
    curr_measure_bin = cellfun(@(x) round(x - 0.5), ...
        curr_measure_cells,'uni',false);   
    curr_measure_cells_dot = cellfun(@(x) (x'*x)/size(x,1), ...
        curr_measure_bin,'uni',false);
    curr_measure_cells_dot_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_dot);
    
    start_day = 1;
    if curr_animal == 1;
        start_day = 3;
    end
    
    all_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef_mean)) = ...
        curr_measure_cells_corrcoef_mean(start_day:end);
    all_day_dot(curr_animal, ...
        start_day:length(curr_measure_cells_dot_mean)) = ...
        curr_measure_cells_dot_mean(start_day:end);
    
%     curr_measure_cells_mean = cell2mat(cellfun(@(x) nanmean(x,2), ...
%         curr_measure_cells,'uni',false));
%     curr_mean_corrcoef = corrcoef(curr_measure_cells_mean);
%     all_mean_correlation(start_day:size(curr_measure_cells_mean,2), ...
%         start_day:size(curr_measure_cells_mean,2),curr_animal) = ...
%         curr_mean_corrcoef(start_day:end,start_day:end);
    
    % correlate mean of trials with all trials
    curr_measure_cells_trialmean_corrcoef = cellfun(@(x) ...
        corrcoef([nanmean(x,2) x]), ...
        curr_measure_cells,'uni',false);
%     concat_day_correlation(curr_animal, ...
%         start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
%         x(2:end,1), curr_measure_cells_trialmean_corrcoef(start_day:end),'uni',false);
    concat_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
        x(tril(true(size(x)),-1)), curr_measure_cells_corrcoef(start_day:end),'uni',false);

    all_day_correlation_full{curr_animal} = curr_measure_cells_corrcoef;
    
%     % get contribution of each cell to mean correlation of the day
%     cell_omit_corr = cellfun(@(x) jackknife(@(y) corrcoef(y),x),...
%         curr_measure_cells,'uni',false);
%     tril_indx = cellfun(@(x) find(tril(true(size(x,2)),-1)),curr_measure_cells,'uni',false);    
%     cell_omit_corr_diff = cellfun(@(x,y,z) z-nanmean(x(:,y),2), ...
%         cell_omit_corr,tril_indx,num2cell(curr_measure_cells_corrcoef_mean),'uni',false);
%  
%     all_cell_omit_corr_diff(curr_animal, ...
%         start_day:length(curr_measure_cells)) = cell_omit_corr_diff(start_day:end);
end

all_reward_reaction = cell(8,15);
all_pos_kidx = cell(8,15);
all_trial_activity = cell(8,15);
for curr_animal = 1:8;
    load(['/usr/local/lab/People/Andy/Data/lever_data' filesep ...
        animals{curr_animal} 'lever_params.mat']);
    for i = 1:length(all_day_correlation_full{curr_animal});
        if curr_animal == 1 && (i == 1 || i == 2)
            continue
        end
        b = all_day_correlation_full{curr_animal}{i};
        b = b < 0;

        try
            kidx = kmeans(b,2);
            [asdf sort_idx] = sort(kidx);
        catch me
            kidx = ones(length(b),1);
            sort_idx = 1:length(b);
        end
        
        [asdf pos_corr_kidx] = min(grpstats(nanmean(b,2),kidx));
        use_kidx = kidx == pos_corr_kidx;
        
        all_pos_idx{curr_animal,i} = use_kidx;
        all_reward_reaction{curr_animal,i} = grab_var.reward_reaction{i}( ...
            trial_type.cued_movement{curr_animal}{i} & ...
            ~trial_type.catch_trials{curr_animal}{i});
        all_trial_activity{curr_animal,i} = nanmean(all_measure_cells{curr_animal,i});
    end
end

% get reaction times from groups
pos_corr_reaction = cellfun(@(x,y) x(y),all_reward_reaction,all_pos_idx,'uni',false);
neg_corr_reaction = cellfun(@(x,y) x(~y),all_reward_reaction,all_pos_idx,'uni',false);

cat_reaction = cell(15,1);
cat_activity = cell(15,1);
for i = 1:14
   cat_reaction{i} = vertcat(all_reward_reaction{:,i});
   cat_activity{i} = horzcat(all_trial_activity{:,i});
end


%% Noise/trial correlations within persistent and transient cells

load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/cell_trial_all_activity_rewmove.mat');

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
all_under_7_corr = cell(8,15);
all_over_7_corr = cell(8,15);
for curr_animal = 2:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    unfilled_pyr = pyr_cells & ~filled_cells;
    unfilled_gad = gad_cells & ~filled_cells;

    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    under_7 = sum_move_day <= 7;
    over_7 = sum_move_day > 7;
    
    % get reliability, frac of trials with activity
    curr_reliability = cellfun(@(x) nanmean(x,2)', ...
            cell_trial_peak_full{curr_animal}.cued.moving,'uni',false);
    
    % pick out persistent/transient cells classified on each day
    move_under_7 = cellfun(@(x,y) x & unfilled_pyr & under_7 & y < 0.2, ...
       classified_cells.move_cells_peak,curr_reliability,'uni',false);

    move_over_7 = cellfun(@(x,y) x & unfilled_pyr & over_7 & y < 0.2, ...
       classified_cells.move_cells_peak,curr_reliability,'uni',false);
   
   % get trial activity for all days, all cells
    trial_activity = cellfun(@(x,y) [x(:,~all(isnan(x))) y(:,~all(isnan(y)))], ...
        cell_trial_peak_full{curr_animal}.cued.moving, ...
        cell_trial_peak_full{curr_animal}.uncued.moving,'uni',false);
    
    % fix nan, ugghh have to do this with loop?
    for i = 1:length(trial_activity);
        trial_activity{i}(isnan(trial_activity{i})) = 0;
    end
    
    % get average noise correlation for persistent/transient cells
    corr_under_7 = cellfun(@(x,y) corrcoef(x(y,:)'), ...
        trial_activity,move_under_7,'uni',false);
    corr_over_7 = cellfun(@(x,y) corrcoef(x(y,:)'), ...
        trial_activity,move_over_7,'uni',false);
    
    tril_under_7 = cellfun(@(x) tril(true(size(x)),-1),corr_under_7,'uni',false);
    tril_over_7 = cellfun(@(x) tril(true(size(x)),-1),corr_over_7,'uni',false);
    
    all_under_7_corr(curr_animal,1:length(trial_activity)) = cellfun(@(x,y) x(y), ...
        corr_under_7,tril_under_7,'uni',false);
    all_over_7_corr(curr_animal,1:length(trial_activity)) = cellfun(@(x,y) x(y), ...
        corr_over_7,tril_over_7,'uni',false);
    
end

all_under_7_corr_cat = cell(1,15);
all_over_7_corr_cat = cell(1,15);
for i = 1:15
   all_under_7_corr_cat{i} = vertcat(all_under_7_corr{:,i});
   all_over_7_corr_cat{i} = vertcat(all_over_7_corr{:,i});
end
figure;hold on;
errorbar(cellfun(@nanmean,all_under_7_corr_cat), ...
    cellfun(@(x) nanstd(x)/nnz(~isnan(x)),all_under_7_corr_cat),'r','linewidth',2);
errorbar(cellfun(@nanmean,all_over_7_corr_cat), ...
    cellfun(@(x) nanstd(x)/nnz(~isnan(x)),all_over_7_corr_cat),'k','linewidth',2);


a = cellfun(@nanmean,all_under_7_corr);
b = cellfun(@nanmean,all_over_7_corr);

% plot number of over/under 7 for each day
num_over_7 = cellfun(@length,all_over_7_corr);
num_under_7 = cellfun(@length,all_under_7_corr);
figure;hold on;
plot(sum(num_under_7),'r','linewidth',2);
plot(sum(num_over_7),'k','linewidth',2);


%% plot classified cell percentage over time


load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/cell_trial_all_activity_rewmove.mat');

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
all_classified_gad = nan(8,15);
all_classified_pyr = nan(8,15);
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    all_classified_pyr(curr_animal,sessions) = ...
        cellfun(@(x) nanmean(x(pyr_unfilled)),classified_cells.move_cells_peak(sessions));
    all_classified_gad(curr_animal,sessions) = ...
        cellfun(@(x) nanmean(x(gad_unfilled)),classified_cells.move_cells_peak(sessions));
end



%% FIGURE 2 - pyramidal cell 


use_animals = [1 2 4 5 6 7 8];

% plot movement-cell classification over time
b = cellfun(@nanmean,pyr_class_all);
b(1,1:2) = NaN;
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = horzcat(pyr_class_all{use_animals,day_combine{i}});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data),'color','k','linewidth',2);
set(gca,'XTick',1:length(day_combine));
set(gca,'XTickLabel',{'1-2' '3-4' '5-6' '7-8' '9-10' '11-12' '13-14'});
xlabel('Session')
ylabel('Fraction of cells classified as movement-related')

% plot reliability vs. persistence 
pyr_move_active = cellfun(@(x,y,r) cell2mat(cellfun(@(z) ...
    any(x(:,z),2),y(r),'uni',false)), ...
   pyr_activity_all,move_epoch_frames,rewarded_movements,'uni',false);
pyr_reliability = cellfun(@(x,y) nanmean(x,2), ...
    pyr_move_active,'uni',false);
bin_edges = [0 0.0001 0.1:0.1:0.9 Inf];
pyr_class_all{1,1} = [];
pyr_class_all{1,2} = [];
reliability_binned = nan(12,8);
for i = 1:8
   curr_reliability = horzcat(pyr_reliability{i,:});
   curr_idx = logical(vertcat(pyr_class_all{i,:})');
   curr_mean_reliability = ...
       arrayfun(@(x) nanmean(curr_reliability(x,curr_idx(x,:))),1:size(curr_reliability,1));
   curr_sum_days = nanmean(curr_idx,2);
   [n bins] = histc(curr_sum_days,bin_edges);
   curr_bins = unique(bins);
   curr_binmean = grpstats(curr_mean_reliability,bins);
   reliability_binned(curr_bins(curr_bins ~= 1),i) = curr_binmean(curr_bins ~= 1);
end
figure;errorbar(nanmean(reliability_binned(:,use_animals)'), ...
   nanstd(reliability_binned(:,use_animals)') ...
   ./sqrt(sum(~isnan(reliability_binned(:,use_animals)'))), ...
    'color','k','linewidth',2);
set(gca,'XTick',2:length(bin_edges)-1);
set(gca,'XTickLabel',{'10' '20' '30' '40' '50' '60' '70' '80' '90' '100'});
xlabel('Persistence (%)')
ylabel('Reliability')
title('Mean reliability of cells by persistence')

% plot mean persistance of movement-related cells over days
mean_persistence = nan(8,14);
for i = 1:8
   sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
   curr_move = vertcat(pyr_class_all{i,:})';
   curr_move_sum = bsxfun(@times,curr_move,nanmean(curr_move,2));
   curr_move_sum(curr_move_sum == 0) = NaN;

   mean_persistence(i,sessions) = nanmean(curr_move_sum);
end
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = mean_persistence(:,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','k','linewidth',2);
set(gca,'XTick',1:length(day_combine));
set(gca,'XTickLabel',{'1-2' '3-4' '5-6' '7-8' '9-10' '11-12' '13-14'});
xlabel('Session')
ylabel('Mean persistence of movement cells')


% get daily overlap
overlap = nan(15,15,8);
new_overlap_frac = nan(8,14);
old_overlap_frac = nan(8,14);
new_frac = nan(8,15);
old_frac = nan(8,15);
for i = 1:8
   sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(i,:)));
   curr_move = vertcat(pyr_class_all{i,:})';
   overlap(sessions,sessions,i) = +curr_move'*+curr_move;
   
   % what frac of new/old move cells are carried between days?
   new_move = diff(curr_move,[],2) > 0;
   new_overlap = +new_move'*curr_move(:,2:end);
   a = diag(new_overlap); b = diag(new_overlap,1);
   new_overlap_frac(i,sessions(2:end-1)) = b./a(1:end-1);
   new_frac(i,sessions(2:end)) = sum(new_move)./sum(curr_move(:,2:end));
   
   old_move = diff(curr_move,[],2) == 0 & curr_move(:,2:end) == 1;
   old_overlap = +old_move'*curr_move(:,2:end);
   a = diag(old_overlap); b = diag(old_overlap,1);
   old_overlap_frac(i,sessions(2:end-1)) = b./a(1:end-1);
   old_frac(i,sessions(2:end)) = sum(old_move)./sum(curr_move(:,2:end));
   
   curr_persistence = sum(curr_move,2)/size(curr_move,2); 
end
overlap_diag = cell2mat(arrayfun(@(x) diag(overlap(:,:,x)),1:8,'uni',false));
overlap_diag_next = cell2mat(arrayfun(@(x) diag(overlap(:,:,x),-1),1:8,'uni',false));
total_overlap = [nan(1,size(overlap_diag,2)); ...
    overlap_diag_next(1:13,:)./overlap_diag(1:13,:)]';
overlap_fraction = [overlap_diag_next(1:13,:)./overlap_diag(2:14,:)]';
new_cells_fraction = [(overlap_diag(2:14,:)-overlap_diag_next(1:13,:)) ./ ...
    overlap_diag(2:14,:)]';
%day_combine = {[2 3] [4 5] [6 7] [8 9] [10 11] [12 13] [14]};
day_combine = num2cell(2:14);
day_combine_data = cell(size(day_combine));
figure; hold on
for i = 1:length(day_combine)
    curr_data = old_frac(use_animals,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','b','linewidth',2);
for i = 1:length(day_combine)
    curr_data = new_frac(use_animals,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','r','linewidth',2);
set(gca,'XTick',1:length(day_combine));
set(gca,'XTickLabel',cellfun(@num2str,day_combine,'uni',false));
xlabel('Session')
legend({'Frac old' 'Frac new'})
ylabel('Composition of current day')

%day_combine = {[2 3] [4 5] [6 7] [8 9] [10 11] [12 13]};
day_combine = num2cell(2:14);
day_combine_data = cell(size(day_combine));
figure; hold on
for i = 1:length(day_combine)
    curr_data = total_overlap(use_animals,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','k','linewidth',2,'linestyle','--');
for i = 1:length(day_combine)
    curr_data = old_overlap_frac(use_animals,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','b','linewidth',2);
for i = 1:length(day_combine)
    curr_data = new_overlap_frac(use_animals,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'color','r','linewidth',2);
set(gca,'XTick',1:length(day_combine));
set(gca,'XTickLabel',cellfun(@num2str,day_combine,'uni',false));
xlabel('Session')
legend({'Old overlap' 'Frac old' 'Frac new'})
ylabel('Overlap with previous day')


%% Active pyr population correlation

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,'uni',false);


pyr_active_mean = cellfun(@(x) nanmean(x,2),pyr_move_active,'uni',false);

pop_corr = nan(15,15,8);
for i = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_active_mean(i,:)));
    %curr_data = horzcat(pyr_active_mean{i,:});
    curr_data = vertcat(pyr_class_all{i,:})';
    pop_corr(sessions,sessions,i) = corrcoef(curr_data);   
end

use_animals = [1 2 4 5 6 7 8];
pop_corr_mean = nanmean(pop_corr(1:14,1:14,use_animals),3);
figure;imagesc(pop_corr_mean);colormap(hot);




%% Pyr: persistence, reliability trends
high_reliability = nan(8,15);
low_reliability = nan(8,15);

mean_reliability = cell(8,1);
for i = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), all_move_cells(i,:)));
    
    curr_move = vertcat(all_move_cells{i,:})';
    
    curr_persistence = repmat(horzcat(all_persistence{i,:}),1,size(curr_move,2));
    curr_reliability = horzcat(all_reliability{i,:});
    
    high_persistence_reliability = nan(size(curr_move));
    high_persistence_reliability(curr_move & curr_persistence >= 0.6) = ...
        curr_reliability(curr_move & curr_persistence >= 0.6);
    high_reliability(i,sessions) = nanmean(high_persistence_reliability);
    
    low_persistence_reliability = nan(size(curr_move));
    low_persistence_reliability(curr_move & curr_persistence < 0.6) = ...
        curr_reliability(curr_move & curr_persistence < 0.6);
    low_reliability(i,sessions) = nanmean(low_persistence_reliability);
    
    curr_idx = logical(curr_move);
    curr_mean_reliability = ...
       arrayfun(@(x) nanmean(curr_reliability(x,curr_idx(x,:))),1:size(curr_reliability,1));
    mean_reliability{i} = curr_mean_reliability';
end



%% get cell activity to reward relationship 
% requires peak_matrix_full and the trial reward one

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

% load trial-by-trial activity
load('/usr/local/lab/People/Andy/Data/Analysis/trial_activity/cell_trial_all_activity_rewmove.mat');

all_mean_correlation = nan(15,15,8);
all_day_correlation = nan(8,15);
all_day_dot = nan(8,15);
concat_day_correlation = cell(8,15);
all_day_correlation_full = cell(8,15);
all_cell_omit_corr_diff = cell(8,15);
all_sum_move_day = cell(8,1);
all_measure_cells = cell(8,15);
all_timing = cell(8,15);

all_move_cells = cell(8,15);

all_activity_transient = nan(8,15);
all_activity_persistent = nan(8,15);

all_activity_transient_std = nan(8,15);
all_activity_persistent_std = nan(8,15);

all_cells = cell(8,1);

all_reliability = cell(8,15);

for curr_animal = 1:8
    animal = animals{curr_animal};
       
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';

    % get usable cells
    pyr_unfilled = pyr_cells & ~filled_cells;
    gad_unfilled = gad_cells & ~filled_cells;
    
    all_cells{curr_animal} = pyr_unfilled;
    
    % get days classified
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    sum_move_day = sum(vertcat(classified_cells.move_cells_peak{:}));
    all_sum_move_day{curr_animal} = sum_move_day(pyr_unfilled);
    
    % pick measure of trial activity
    curr_measure = trial_activity.peak;
% 
%     curr_measure_cells = cellfun(@(x,y) max(...
%         cat(3,x(pyr_unfilled,:),y(pyr_unfilled,:)),[],3), ...
%         curr_measure{curr_animal}.cued.moving, ...
%         curr_measure{curr_animal}.cued.quiescent,'uni',false);
    
    if strmatch(animal,'AP71')
        classified_cells.move_cells_peak{1} = false(1,length(roi_labels));
        classified_cells.move_cells_peak{2} = false(1,length(roi_labels));
    end
    
    persistent_cells = sum_move_day > 7;
    transient_cells = sum_move_day > 1 & sum_move_day <= 7;
    other_cells = ~persistent_cells & ~transient_cells;
    
    %curr_measure_cells = cellfun(@(x,y) x(pyr_unfilled & y,:), ...
    %    curr_measure{curr_animal}.cued.moving,classified_cells.move_cells_peak,'uni',false);
    %curr_measure_cells = cellfun(@(x) x(pyr_unfilled,:), ...
    %    curr_measure{curr_animal}.cued.moving,'uni',false);
    curr_measure_cells = cellfun(@(x,y) [x(pyr_unfilled,:) y(pyr_unfilled,:)], ...
        curr_measure{curr_animal}.cued.moving, ...
        curr_measure{curr_animal}.uncued.moving,'uni',false);
    
    sessions = find(cellfun(@(x) ~isempty(x),curr_measure_cells));
    if curr_animal == 1
        sessions = 3:14;
    end
    all_timing(curr_animal,sessions) = cellfun(@(x) x(pyr_unfilled), ...
        classified_cells.move_timing_peak(sessions),'uni',false);
    
    all_move_cells(curr_animal,1:length(classified_cells.move_cells_peak)) = ...
        cellfun(@(x) x(pyr_unfilled),classified_cells.move_cells_peak,'uni',false);
    
    % get/save reliabilities
    all_reliability(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = ...
        cellfun(@(x) nanmean(x,2), curr_measure_cells,'uni',false);
    if strmatch(animal,'AP71')
        all_reliability{1,1} = [];
        all_reliability{1,2} = [];
    end
    
    for i = 1:length(curr_measure_cells)
        nan_trials = all(isnan(curr_measure_cells{i}));
        curr_measure_cells{i}(:,nan_trials) = [];
        curr_measure_cells{i}(isnan(curr_measure_cells{i})) = 0;
    end
    
    all_measure_cells(curr_animal,1:length(curr_measure_cells)) = curr_measure_cells;
    
%     all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
%         nanmean(nansum(x(pyr_unfilled,:)/sum(pyr_unfilled))), ...
%         curr_measure{curr_animal}.cued.moving);
    % I'm kind of butchering this variable to get stuff out temporarily
    all_activity_transient(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nansum(any(x(pyr_unfilled&transient_cells,:),2))/sum(pyr_unfilled&transient_cells), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanmean(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    all_activity_transient_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & transient_cells,:)/sum(pyr_unfilled & transient_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    all_activity_persistent_std(curr_animal,1:size(curr_measure{curr_animal}.cued.moving,2)) = cellfun(@(x) ...
        nanstd(nansum(x(pyr_unfilled & persistent_cells,:)/sum(pyr_unfilled & persistent_cells))), ...
        curr_measure{curr_animal}.cued.moving);
    
    curr_measure_cells_corrcoef = cellfun(@corrcoef, ...
        curr_measure_cells,'uni',false);
    curr_measure_cells_corrcoef_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_corrcoef); 
    
    curr_measure_bin = cellfun(@(x) round(x - 0.5), ...
        curr_measure_cells,'uni',false);   
    curr_measure_cells_dot = cellfun(@(x) (x'*x)/size(x,1), ...
        curr_measure_bin,'uni',false);
    curr_measure_cells_dot_mean = cellfun(@(x) ...
        nanmean(x(tril(true(size(x)),-1))), curr_measure_cells_dot);
    
    start_day = 1;
    if curr_animal == 1;
        start_day = 3;
    end
    
    all_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef_mean)) = ...
        curr_measure_cells_corrcoef_mean(start_day:end);
    all_day_dot(curr_animal, ...
        start_day:length(curr_measure_cells_dot_mean)) = ...
        curr_measure_cells_dot_mean(start_day:end);
    
%     curr_measure_cells_mean = cell2mat(cellfun(@(x) nanmean(x,2), ...
%         curr_measure_cells,'uni',false));
%     curr_mean_corrcoef = corrcoef(curr_measure_cells_mean);
%     all_mean_correlation(start_day:size(curr_measure_cells_mean,2), ...
%         start_day:size(curr_measure_cells_mean,2),curr_animal) = ...
%         curr_mean_corrcoef(start_day:end,start_day:end);
    
    % correlate mean of trials with all trials
    curr_measure_cells_trialmean_corrcoef = cellfun(@(x) ...
        corrcoef([nanmean(x,2) x]), ...
        curr_measure_cells,'uni',false);
%     concat_day_correlation(curr_animal, ...
%         start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
%         x(2:end,1), curr_measure_cells_trialmean_corrcoef(start_day:end),'uni',false);
    concat_day_correlation(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef)) = cellfun(@(x) ...
        x(tril(true(size(x)),-1)), curr_measure_cells_corrcoef(start_day:end),'uni',false);

    all_day_correlation_full(curr_animal, ...
        start_day:length(curr_measure_cells_corrcoef)) = ...
        curr_measure_cells_corrcoef(start_day:end);
%     
%     % get contribution of each cell to mean correlation of the day
%     cell_omit_corr = cellfun(@(x) jackknife(@(y) corrcoef(y),x),...
%         curr_measure_cells,'uni',false);
%     tril_indx = cellfun(@(x) find(tril(true(size(x,2)),-1)),curr_measure_cells,'uni',false);    
%     cell_omit_corr_diff = cellfun(@(x,y,z) z-nanmean(x(:,y),2), ...
%         cell_omit_corr,tril_indx,num2cell(curr_measure_cells_corrcoef_mean),'uni',false);
%  
%     all_cell_omit_corr_diff(curr_animal, ...
%         start_day:length(curr_measure_cells)) = cell_omit_corr_diff(start_day:end);
end



dropped_cells = cell(8,14);
retained_cells = cell(8,14);
for i = 1:8
   sessions = find(cellfun(@(x) ~isempty(x), all_move_cells(i,:)));
   curr_reliability = horzcat(all_reliability{i,:});
   curr_move = vertcat(all_move_cells{i,:})';

   dropped_cells(i,sessions(1:end-1)) = mat2cell(diff(curr_move,[],2) < 0, size(curr_move,1), ...
       ones(size(curr_move,2)-1,1));
   retained_cells(i,sessions(1:end-1)) = mat2cell(diff(curr_move,[],2) == 0 & ...
       curr_move(:,2:end) == 1, size(curr_move,1), ...
       ones(size(curr_move,2)-1,1));
end


% get fraction of activity associated with rewarded movements
% dropped_rewarded = nan(8,15);
% retained_rewarded = nan(8,15);
% for i = 1:8
%    sessions = find(cellfun(@(x) ~isempty(x), all_move_cells(i,:)));
%    curr_reliability = horzcat(all_reliability{i,:});
%    curr_move = vertcat(all_move_cells{i,:})';
% 
%    dropped_cells = mat2cell(diff(curr_move,[],2) < 0, size(curr_move,1), ...
%        ones(size(curr_move,2)-1,1));
%    retained_cells = mat2cell(diff(curr_move,[],2) == 0 & ...
%        curr_move(:,2:end) == 1, size(curr_move,1), ...
%        ones(size(curr_move,2)-1,1));
%    
%    dropped_cell_activity = cellfun(@(x,y) x(y,:), ...
%        all_measure_cells(i,sessions(1:end-1)),dropped_cells,'uni',false);
%    retained_cell_activity = cellfun(@(x,y) x(y,:), ...
%        all_measure_cells(i,sessions(1:end-1)),retained_cells,'uni',false);
%    
%    % get total number of spikes per cell
%    rewarded_activity = cellfun(@(x,y) nansum(x,2)./nansum(y(all_cells{i},:),2), ...
%        all_measure_cells(i,sessions),peak_matrix_full{i},'uni',false);
%    
%    dropped_rewarded(i,sessions(1:end-1)) = cellfun(@(x,y) nanmean(x(y)), ...
%        rewarded_activity(sessions(1:end-1)),dropped_cells);
%    retained_rewarded(i,sessions(1:end-1)) = cellfun(@(x,y) nanmean(x(y)), ...
%        rewarded_activity(sessions(1:end-1)),retained_cells);
%    
% end





%% Gad cell

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

low_kidx_cc = nan(8,15);
high_kidx_cc = nan(8,15);

raw_low_kidx_cc = nan(8,15);
raw_high_kidx_cc = nan(8,15);

gad_grouped_cc = cell(8,1);
gad_concat_cc = cell(8,1);
gad_high_kidx = cell(8,1);

gad_activity_all = cell(8,15);
gad_activity_thresh_all = cell(8,15);
pyr_activity_all = cell(8,15);

for curr_animal = 1:8
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
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_caEvents_5.mat'])
        
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
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
    
    gad_activity = cell(size(days));
    gad_activity_thresh = cell(size(days));
    pyr_activity = cell(size(days));
    for curr_day = sessions
        day = num2str(days{curr_day});
        
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        gad_activity{curr_day} = im.roi_trace_df(gad_unfilled,:);
        
        % get and save pyramidal cell activity (peaks)
%         [peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df(pyr_unfilled,:));     
%         peak_matrix = zeros(size(im.roi_trace_df));
%         pyr_unfilled_idx = find(pyr_unfilled);
%         for i = 1:length(pyr_unfilled_idx)
%             curr_pyr = pyr_unfilled_idx(i);
%             peak_matrix(curr_pyr,peak_frames{i}) = 1;
%         end
        pyr_activity{curr_day} = AP_caEvents(im.roi_trace_df(pyr_unfilled,:));
        
        gad_activity_thresh{curr_day} = ...
        AP_caEvents(im.roi_trace_df(gad_unfilled,:),[],1:sum(gad_unfilled));
        
    end
    pyr_activity_all(curr_animal,sessions) = pyr_activity(sessions);
    gad_activity_all(curr_animal,sessions) = gad_activity(sessions);
    gad_activity_thresh_all(curr_animal,sessions) = gad_activity_thresh(sessions);
    % threshold gad activity
%     gad_activity_thresh = gad_activity;
%     for curr_day = sessions
%         gad_activity{curr_day}(isnan(gad_activity{curr_day})) = 0;
%         temp_smooth = cell2mat(arrayfun(@(x) smooth(gad_activity{curr_day} ...
%             (x,:),30,'loess')',1:size(gad_activity{curr_day},1),'uni',false)');
%         noise_est = mean(abs(gad_activity{curr_day} - temp_smooth),2);
%         noise_thresh = noise_est*3;
%         below_thresh = bsxfun(@lt,temp_smooth,noise_thresh);
%         gad_activity_thresh{curr_day}(below_thresh) = 0;
%     end
    
    
           
    concat_gad_activity = horzcat(gad_activity_thresh{:});
    concat_gad_activity(isnan(concat_gad_activity)) = 0;
    concat_gad_activity_cc = corrcoef(concat_gad_activity');
    kidx = kmeans(corrcoef(concat_gad_activity_cc),2);
    [asdf low_kidx] = min( ...
    [nanmean(AP_itril(concat_gad_activity_cc(kidx == 1,kidx == 1),-1)) ...
        nanmean(AP_itril(concat_gad_activity_cc(kidx == 2,kidx == 2),-1))]);
    [asdf high_kidx] = max( ...
    [nanmean(AP_itril(concat_gad_activity_cc(kidx == 1,kidx == 1),-1)) ...
        nanmean(AP_itril(concat_gad_activity_cc(kidx == 2,kidx == 2),-1))]);
    
    low_kidx_cc(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril( ...
        corrcoef(x(kidx == low_kidx,:)'),-1)), ...
        gad_activity_thresh(sessions));
    high_kidx_cc(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril( ...
        corrcoef(x(kidx == high_kidx,:)'),-1)), ...
        gad_activity_thresh(sessions));
    
    raw_low_kidx_cc(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril( ...
        corrcoef(x(kidx == low_kidx,:)'),-1)), ...
        gad_activity(sessions));
    raw_high_kidx_cc(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril( ...
        corrcoef(x(kidx == high_kidx,:)'),-1)), ...
        gad_activity(sessions));
    
    gad_grouped_cc{curr_animal}(sessions) = cellfun(@(x) ...
        corrcoef(x([find(kidx == low_kidx);find(kidx == high_kidx)],:)'), ...
        gad_activity_thresh(sessions),'uni',false);
    
    gad_concat_cc{curr_animal} = concat_gad_activity_cc( ...
        [find(kidx == low_kidx);find(kidx == high_kidx)], ...
        [find(kidx == low_kidx);find(kidx == high_kidx)]);
    gad_high_kidx{curr_animal} = kidx == high_kidx;
    
    curr_animal
end


day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = a(:,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'k','linewidth',2);


% get average activity / above threshold for high/low
high_avg_gad_act = nan(8,15);
high_avg_gad_abovethresh = nan(8,15);
low_avg_gad_act = nan(8,15);
low_avg_gad_abovethresh = nan(8,15);
for i = 1:8
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(i,:)));
    
    gad_nan_act = gad_activity_thresh_all(i,sessions);
    for j = 1:length(gad_nan_act);
        gad_nan_act{j}(gad_nan_act{j} == 0) = NaN;
    end
    
    high_avg_gad_act(i,sessions) = cellfun(@(x) ...
        nanmean(nanmean(x(gad_high_kidx{i},:))),gad_nan_act);
    high_avg_gad_abovethresh(i,sessions) = cellfun(@(x) ...
        nanmean(nanmean(x(gad_high_kidx{i},:)>0)),gad_activity_thresh_all(i,sessions));
    
    low_avg_gad_act(i,sessions) = cellfun(@(x) ...
        nanmean(nanmean(x(~gad_high_kidx{i},:))),gad_nan_act);
    low_avg_gad_abovethresh(i,sessions) = cellfun(@(x) ...
        nanmean(nanmean(x(~gad_high_kidx{i},:)>0)),gad_activity_thresh_all(i,sessions));
end

day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = low_avg_gad_act(:,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),day_combine_data), ...
    'k','linewidth',2);



% shuffle to get correlations by chance
n_rep = 10;
real_corr = nan(8,15);
shuff_corr = nan(8,15);
real_overlap = nan(8,15);
real_overlap_full = cell(8,15);
for i = 1:8
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(i,:)));    
    curr_gad_act = cellfun(@(x) x(gad_high_kidx{i},:), ...
        gad_activity_thresh_all(i,sessions),'uni',false);
    
    %real_corr(i,sessions) = cellfun(@(x) nanmean(AP_itril(corrcoef(x'),-1)),curr_gad_act);
    real_corr(i,sessions) = cellfun(@(x) nanmean(AP_itril(corrcoef([x>0]'),-1)),curr_gad_act);
    
    curr_overlap = cellfun(@(x) +[x>0]*+[x>0]',curr_gad_act,'uni',false);
    real_overlap(i,sessions) = cellfun(@(x) nanmean(AP_itril( ...
        x./sqrt([diag(x)*diag(x)']),-1)),curr_overlap);
    real_overlap_full(i,sessions) = curr_overlap;
%     shuff_corr_rep = nan(n_rep,length(sessions));
%     for rep = 1:n_rep
%         shuff_corr_rep(rep,:) = cellfun(@(x) ...
%             nanmean(AP_itril(corrcoef(shake(x,2)'),-1)),curr_gad_act);
%         rep
%     end
%     shuff_corr(i,sessions) = nanmean(shuff_corr_rep);
    i
end

gad_trial_corr = nan(8,15);
for curr_animal = 1:8
    animal = animals{curr_animal};
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    gad_unfilled = find(gad_cells & ~filled_cells);
    pyr_unfilled = find(pyr_cells & ~ filled_cells);
    gad_trial_corr(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril(corrcoef(x( ...
       gad_unfilled(gad_high_kidx{curr_animal}),:)'),-1)), ...
       trial_activity.df{curr_animal}.cued.moving(sessions));
%     gad_trial_corr(curr_animal,sessions) = cellfun(@(x) nanmean(AP_itril(corrcoef(x( ...
%         pyr_unfilled,:)'),-1)), ...
%         trial_activity.peak{curr_animal}.cued.moving(sessions));
end
    
    
curr_gad_activity = gad_activity_thresh_all(4,:);
curr_pyr_activity = pyr_activity_all(4,:);
    
    
spread_filter = ones(1,10);
pyr_spread = cellfun(@(x) conv2(x,spread_filter),pyr_activity_all,'uni',false);
pyr_spread_overlap = cellfun(@(x) nanmean(AP_itril( ...
    [x*x']./sqrt([diag([x*x'])*diag([x*x'])']),-1)),pyr_spread);
    
    

    
%% Gad cells - show low active = preference for higher active trials

reliability_popavg = cell(8,15);
reliability_popavg_shuff = cell(8,15);
distance_corr = cell(8,15);

for curr_animal = 1:8
    
    animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
    animal = animals{curr_animal};
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
        animal '_classified_cells_move_epoch_shuffle.mat'])
    
    % load ROIs and get distances
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
    
    sessions = find(cellfun(@(x) ~isempty(x),classified_cells.move_cells_peak));
    
    
    %figure; hold on;
    col = copper(max(sessions));
    
    
    for day = sessions
        
        reliability_popavg{curr_animal,day} = nan(6,1);
        reliability_popavg_shuff{curr_animal,day} = nan(6,1);
        
        high_gad = gad_unfilled(gad_high_kidx{curr_animal});
        %a = trial_activity.peak{curr_animal}.cued.moving{day}(high_gad,:);
        a = trial_activity.peak{curr_animal}.cued.moving{day}(pyr_unfilled,:);
        a(:,all(isnan(a))) = [];
        trial_mean = nanmean(a,1);
        t = (a*mean(a)')./sum(a,2);
        k = arrayfun(@(x) mean(trial_mean(logical(a(x,:)))),1:size(a,1));
        [asdf sort_idx] = sort(t);
        
        % shuffle the data to get chance preference for each cell
        n_rep = 200;
        rep_t = nan(n_rep,size(a,1));
        for rep = 1:n_rep
            shake_a = shake(a,2);
            trial_mean = mean(shake_a,1);
            rep_t(rep,:) = (shake_a*mean(shake_a)')./sum(shake_a,2);
            %rep_k(rep,:) = arrayfun(@(x) mean(trial_mean( ...
            %    randperm(length(trial_mean),sum(shake_a(x,:))))),1:size(shake_a,1));
            rep
        end
        rep_t_mean = nanmean(rep_t);
        
        % bin cells by fraction trials active/reliability
        bin_edges = [0:0.2:0.8 Inf];
        [asdf reliability_bins] = histc(nanmean(a,2),bin_edges);
        %norm_t = (t-nanmean(a(:)))./(max(nanmean(a,1))-nanmean(a(:)));
        norm_t = t./nanmean(a(:));
        t_grouped = grpstats(norm_t,reliability_bins);
        reliability_popavg{curr_animal,day}(unique(reliability_bins)) = t_grouped;
        
        norm_t_shuff = [rep_t_mean./nanmean(a(:))]';
        t_shuff_grouped = grpstats(norm_t_shuff,reliability_bins);
        reliability_popavg_shuff{curr_animal,day}(unique(reliability_bins)) = ...
            t_shuff_grouped;
        
        %plot(unique(reliability_bins),t_grouped,'col',col(day,:));
        
        %plot(sum(a,2)/size(a,2),t./nanmean(a(:)),'.','color',col(day,:));%'.k');
        %plot(sum(a,2)/size(a,2),nanmean(rep_k)/nanmean(a(:)),'.r');
        
        %     plot(sum(a,2)/size(a,2), ...
        %         (t-nanmean(a(:)))./(max(nanmean(a,1))-nanmean(a(:))), ...
        %         '.','color',col(day,:));%'.k');
        %     plot(sum(a,2)/size(a,2), ...
        %         (nanmean(rep_k)-nanmean(a(:)))./(max(nanmean(a,1))-nanmean(a(:))),'.r');
        %ylabel('Preferred pop activity');
        %xlabel('Fraction trials active');
        
        
        % get pyr-gad distance correlation relationship
        distance_bin_edges = [0:20:480 Inf];
        distance_corr{curr_animal,day} = nan(length(distance_bin_edges),1);
        all_act = trial_activity.peak{curr_animal}.cued.moving{day};
        all_act(:,all(isnan(all_act))) = [];
        all_cc = corrcoef(all_act');
        curr_cc = all_cc(pyr_unfilled,high_gad);
        curr_distance = roi_dist(pyr_unfilled,high_gad);
        [n distance_bins] = histc(curr_distance,distance_bin_edges);
        distance_corr{curr_animal,day}(unique(distance_bins)) = ...
            grpstats(curr_cc(:),distance_bins(:),'nanmean');
    end
    curr_animal
end

a = horzcat(reliability_popavg{:})';
b = horzcat(reliability_popavg_shuff{:})';
figure; hold on;
errorbar(nanmean(a),nanstd(a)./sqrt(sum(~isnan(a))),'k','linewidth',2)
errorbar(nanmean(b),nanstd(b)./sqrt(sum(~isnan(b))),'--r','linewidth',2)
legend({'Real','Shuffle'})
xlim([0 6])
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'0-20' '20-40' '40-60' '60-80' '80-100'});
xlabel('Reliability')
ylabel('Fraction active cells/mean')


day_groups = {[1 2] [3 5] [6 10] [11 14]};
col = copper(length(day_groups));
figure; hold on;
for i = 1:length(day_groups)
    a = horzcat(reliability_popavg{:,day_groups{i}})';
    b = horzcat(reliability_popavg_shuff{:,day_groups{i}})';
    errorbar(nanmean(a,1),nanstd(a)./sqrt(sum(~isnan(a),1)),'color',col(i,:),'linewidth',2)
    errorbar(nanmean(b),nanstd(b)./sqrt(sum(~isnan(b),1)),'--','color',col(i,:),'linewidth',2)   
end




% % TEMP: plot activity preference of cells
% % load ROIs and get ROI centroids and distances
% a = trial_activity.peak{curr_animal}.cued.moving{day}(pyr_unfilled,:);
% a(:,all(isnan(a))) = [];
% t_pyr = ((a*mean(a)')./sum(a,2))/nanmean(a(:));
% 
% gad_high = gad_unfilled(gad_high_kidx{curr_animal});
% b = trial_activity.peak{curr_animal}.cued.moving{day}(gad_high,:);
% b(:,all(isnan(b))) = [];
% t_gad = ((b*mean(b)')./sum(b,2))/nanmean(b(:));
% 
% roi_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
% roi_dir = dir(roi_path);
% roi_filenames = {roi_dir(:).name};
% roi_files = cellfun(@(x) any(strfind(x,'.roi')),roi_filenames);
% roi_files_sort = sort(roi_filenames(roi_files));
% load([roi_path filesep roi_files_sort{end}],'-MAT');
% num_rois = length(polygon.ROI);
% roi_center_x = zeros(num_rois,1);
% roi_center_y = zeros(num_rois,1);
% for curr_roi  = 1:num_rois
%     [area roi_center_x(curr_roi) roi_center_y(curr_roi)] = ...
%         polycenter(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2));
% end
% roi_dist = nan(length(roi_center_x));
% for k = 1:length(roi_dist)
%     for j = 1:length(roi_dist)
%         curr_dist = sqrt((roi_center_x(k) - roi_center_x(j))^2 + ...
%             (roi_center_y(k) - roi_center_y(j))^2);
%         roi_dist(k,j) = curr_dist;
%     end
% end
% 
% figure;
% roi_handle = [];
% for i = 1:length(polygon.ROI)
%     roi_handle(i) = patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0 0 0]);
% end
% 
% t_pyr = (t_pyr-nanmean(a(:))) / (max(t_pyr)-nanmean(a(:)));
% t_gad = (t_gad-nanmean(a(:))) / (max(t_gad)-nanmean(a(:)));
% for i = 1:length(pyr_unfilled)
%     if isnan(t_pyr(i)) || sign(t_pyr(i)) == -1
%         continue
%     end
%     set(roi_handle(pyr_unfilled(i)),'FaceColor',[0 t_pyr(i) 0]);
% end
% for i = 1:length(gad_high)
%     if isnan(t_gad(i)) || sign(t_gad(i)) == -1
%         continue
%     end
%     set(roi_handle(gad_high(i)),'FaceColor',[t_gad(i) 0 0]);
% end
% ylim([-512 0]);
% xlim([0 512]);

%% show gad firing hierarchy over whole session: EDIT! FIXED BELOW!
hierarchy_score = nan(8,15);
for curr_animal = 1:8
    animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
    animal = animals{curr_animal};
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
for curr_session = sessions

% % load ROI labels and identify cell type
%     analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
%     roilabel_name = [animal '_roilabels.roilabel'];
%     load([analysis_path filesep roilabel_name],'-MAT')
%     cells = 1:length(roi_labels);
%     gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
%     pyr_cells = ~gad_cells;
%     filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
%     pyr_unfilled = find(pyr_cells & ~filled_cells);
%     gad_unfilled = find(gad_cells & ~filled_cells);    
    
spread_filter = ones([1 60]);
%spread_filter = ones([1 1]);

x = gad_activity_thresh_all{curr_animal,curr_session}(gad_high_kidx{curr_animal},:);
%x = pyr_activity_all{curr_animal,curr_session};
%x = trial_activity.peak{curr_animal}.cued.moving{day}(gad_unfilled,:);
%x = gad_activity_thresh_all{curr_animal,curr_session};
y = +[conv2(x,spread_filter,'same') > 0];

[asdf sort_rows] = sort(nanmean(y,2));
y_id = bsxfun(@times,y(sort_rows,:)>0,[size(y,1):-1:1]');

% a = sum(y_id > 0);
% b = nanmean(y_id);
% temp_cc = corrcoef(a,b);
% hier_cc = temp_cc(2);
% 
% reps = 200;
% hier_cc_shuff = nan(reps,1);
% for rep = 1:reps
%     id_shake = repmat(randperm(size(y,1))',1,size(y,2));
%     y_id_shake = y.*id_shake;
%     
%     a = sum(y_id_shake > 0);
%     b = nanmean(y_id_shake);
%     temp_cc = corrcoef(a,b);
%     hier_cc_shuff(rep) = temp_cc(2);
%     
%     rep
% 
% end

% different approach: find coincident activity events, get hierarchy score
r = nansum(y);
maxtab = peakdet(r,3);
coactive_events = zeros(size(y,1),size(maxtab,1));
for i = 1:size(maxtab,1);
    coactive_events(find(y(:,maxtab(i,1))),i) = 1;
end

%coactive_events = y;

[asdf sort_rows] = sort(nanmean(coactive_events,2));
[asdf sort_rows2] = sort(nanmean(coactive_events,1),'descend');
coactive_events = coactive_events(sort_rows,sort_rows2);

% get how off from a perfect hierarchy it is
coactive_perfect_hierarchy = zeros(size(coactive_events));
for i = 1:size(coactive_events,1)
   coactive_perfect_hierarchy(i,1:sum(coactive_events(i,:))) = 1;
end

coactive_hierarchy_score = coactive_events(:)'*coactive_perfect_hierarchy(:);

coactive_perfect_hierarchy_score = coactive_perfect_hierarchy(:)'*coactive_perfect_hierarchy(:);

reps = 100;
coactive_hierarchy_shuff_score = nan(reps,1);
coactive_perfect_hierarchy_shuff_score = nan(reps,1);
for rep = 1:reps
    coactive_shake = shake(coactive_events,2);
    [asdf sort_rows2] = sort(nanmean(coactive_shake,1),'descend');
    
    coactive_shake_sort = coactive_shake(:,sort_rows2);
    
    coactive_hierarchy_shuff_score(rep) = ...
        coactive_shake_sort(:)'*coactive_perfect_hierarchy(:);
    
    rep
end

hierarchy_score(curr_animal,curr_session) = ...
    (coactive_hierarchy_score-nanmean(coactive_hierarchy_shuff_score)) ...
    /coactive_perfect_hierarchy_score;
end
end



%% Store all movement traces (also tried correlation v lag for gad cells)
movement_trace_all = cell(8,15);
for curr_animal = 1:8
    
    animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
    animal = animals{curr_animal};
    
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
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
    
    for curr_day = sessions
        day = days{curr_day};
        
        % load split image file to get loop sizes
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name],'im');
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        lever_active_full = cell(length(xsg_filenames),1);
        lever_velocity_full = cell(length(xsg_filenames),1);
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
            curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
            % if trials are not consecutive, delete the offending trial
            consec_trials = [0;diff(curr_trial_list(:,2))];
            curr_trial_list(consec_trials ~= 1,:) = [];
            
            [lever_active lever_force_smooth lever_force_velocity] = ...
                AP_parseLeverMovement(xsg_data);
            
            
            
            %save the active lever portions in frames
            curr_lever_active_frames = round((find(lever_active)/1000) ...
                *bhv.framerate);
            lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
                curr_lever_active_frames <= loop_size(curr_xsg) & ...
                curr_lever_active_frames > 0)) + loop_frames(curr_xsg);
            
            [n d] = rat(loop_size(curr_xsg)/length(lever_force_velocity));
            lever_force_velocity_resample = resample(lever_force_velocity,n,d);
            lever_velocity_full{curr_xsg} = ...
                lever_force_velocity_resample(1:loop_size(curr_xsg));
        end
        
        lever_velocity = vertcat(lever_velocity_full{:});
        
        lever_movement = vertcat(lever_active_full{:});
        movement_trace = zeros(1,size(im.roi_trace_df,2));
        movement_trace(lever_movement) = 1;
        
        movement_trace_all{curr_animal,curr_day} = movement_trace;
        
        curr_day
        
    end
    
    animal
end

figure;
all_binned_lagcorr = cell(8,15);
for curr_animal = 1:8
subplot(3,3,curr_animal);hold on;
%curr_animal = 3;
sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
clear a b c d
a = gad_activity_thresh_all(curr_animal,sessions);
b(sessions) = cellfun(@(x) +[x(gad_high_kidx{curr_animal},:)>0],a,'uni',false);
c(sessions) = cellfun(@(x,y) corrcoef([x;y]'),b(sessions),movement_trace_all(curr_animal,sessions),'uni',false);
d(sessions) = cellfun(@(x) x(end,1:end-1),c(sessions),'uni',false);
% this was for checking out movement-gad activity lag
%figure; hold on;
col = jet(max(sessions));
for curr_session = sessions
    rr = nan(sum(gad_high_kidx{curr_animal}),1);
    jj = nan(size(rr));
    curr_pyr_activity = +[sum(pyr_activity_all{curr_animal,curr_session}) > 0];
    curr_pyr_activity_spread = +[smooth(curr_pyr_activity,20) > 0];
    for i = 1:length(rr)
        [r kk] = max(xcov(movement_trace_all{curr_animal,curr_session},b{curr_session}(i,:),'coeff'));
        %[r kk] = max(xcov(curr_pyr_activity_spread,b{curr_session}(i,:)));
        rr(i) = size(b{curr_session},2) - kk;
        jj(i) = r;%/sum(b{curr_session}(i,:));
    end
    %figure;plot(rr(d{curr_session}>0.1),jj(d{curr_session}>0.1),'.k')  
    %plot(rr,jj,'.','color',col(curr_session,:));xlim([-10 100]);
    
    % cut out cells with less than 0.1 correlation with movement
    lag_bin_edges = [0:10:100];
    movecorr_gad = d{curr_session} > 0.1;
    [n lag_bins] = histc(rr(movecorr_gad),lag_bin_edges);
    lag_bin_mean = grpstats(jj(movecorr_gad),lag_bins);
    lag_bin_sem = grpstats(jj(movecorr_gad),lag_bins,'sem');
    unique_lag_bins = unique(lag_bins);
    
    errorbar(unique_lag_bins(unique_lag_bins ~= 0), ...
        lag_bin_mean(unique(lag_bins) ~= 0), ...
        lag_bin_sem(unique(lag_bins) ~= 0),'color',col(curr_session,:));
    
    all_binned_lagcorr{curr_animal,curr_session} = nan(length(lag_bin_edges),1);
    all_binned_lagcorr{curr_animal,curr_session}(unique_lag_bins(...
        unique_lag_bins ~= 0)) = lag_bin_mean(unique_lag_bins ~= 0);
    
end
end
% plot stacked lag-corrected gad activity
figure; hold on
for i = 1:length(ff)
    plot(b{curr_session}(ff(i),rr(ff(i)):end-rr(ff(i)))+(1.5*i),'k')
end
plot(movement_trace_all{4,curr_session}+(1.5*(i+1)),'r');


% find the best smoothing value for pyr activity
pyr_smooth_curves = cell(size(gad_activity_thresh_all));
for curr_animal = 1:8
    %curr_animal = 1
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
    col = jet(length(sessions));
    for curr_session = sessions
        %curr_session = 1;
        blarg = nan(100,1);
        blarg_xcorr = nan(100,1);
        curr_pyr_activity = +[sum(pyr_activity_all{curr_animal,curr_session}) > 0];
        for i = 1:100
            temp = corrcoef(+[smooth(curr_pyr_activity,i) > 0],movement_trace_all{curr_animal,curr_session});
            blarg(i) = temp(2);
            
            [bb mm] = max(xcov(movement_trace_all{curr_animal,curr_session},+[smooth(curr_pyr_activity,i) > 0]));
            blarg_xcorr(i) = bb;%/sum(+[smooth(curr_pyr_activity,i) > 0]);
        end
        %plot(blarg/max(blarg),'k')
        pyr_smooth_curves{curr_animal,curr_session} = blarg/max(blarg);
        %plot(blarg_xcorr/max(blarg_xcorr),'r')
    end
    curr_animal
end





%% Gad/pyr hierarchy via activity during all movement times

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,'uni',false);
% pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,'uni',false);

% this is to use only the movement-classified ones
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);

gad_move_sort = cellfun(@(x) sortrows([x nanmean(x,2)],-(size(x,2)+1)), ...
    gad_move_active,'uni',false);

gad_hierarchy_cells = cellfun(@(x) arrayfun(@(y) ...
    sum(x(1:sum(x(:,y)),y)),1:size(x,2)-1),gad_move_sort,'uni',false);

% get the chance of having active cells be the most active ones
gad_hierarchy_chance = cellfun(@(x,y) arrayfun(@(z) ...
    (nchoosek(sum(x(:,z)),y(z)) * ...
    nchoosek(size(x,1)-sum(x(:,z)),sum(x(:,z))-y(z))) / ...
      nchoosek(size(x,1),sum(x(:,z))),1:length(y)), ...
      gad_move_active,gad_hierarchy_cells,'uni',false);
  
gad_hierarchy_chance_mean = cellfun(@nanmean,gad_hierarchy_chance);

day_combine = {[1 2],[3 4],[5 6],[7 8],[9 10],[11 12],[13 14]};
day_combine_data = cell(size(day_combine));
for i = 1:length(day_combine)
    curr_data = gad_hierarchy_chance_mean(:,day_combine{i});
    day_combine_data{i} = curr_data(:);
end
figure;errorbar(cellfun(@nanmean,day_combine_data), ...
    cellfun(@(x) nanstd(x)./sqrt(length(x)),day_combine_data),'k','linewidth',2);



% how to find chance of having given number of correct hierarchy cells in a
% trial:
%
% (nchoosek(hierarchy_group,hierarchy_cells) * ...
%     nchoosek(num_cells-hierarchy_group,active_cells-hierarchy_cells)) / ...
%       nchoosek(num_cells,active_cells);



%% ''

gad_num_active = cellfun(@(x) sum(x(:)),gad_move_sort);

gad_hierarchy_cells = cellfun(@(x) sum(arrayfun(@(y) ...
    sum(x(1:sum(x(:,y)),y)),1:size(x,2)-1)),gad_move_sort);

% and try also getting chance by shuffling active cells on a trial
rep = 100;
all_hierarchy_cells_shuffle = cell(size(gad_move_shuffle));
for i = 1:rep
    gad_move_shuffle = cellfun(@(x) shake(x,1),gad_move_active,'uni',false);
    
    gad_move_shuffle_sort = cellfun(@(x) sortrows([x nanmean(x,2)],-(size(x,2)+1)), ...
    gad_move_shuffle,'uni',false);

    gad_hierarchy_cells_shuffle = cellfun(@(x) sum(arrayfun(@(y) ...
    sum(x(1:sum(x(:,y)),y)),1:size(x,2)-1)),gad_move_shuffle_sort,'uni',false);

    all_hierarchy_cells_shuffle = cellfun(@(x,y) [x;y],all_hierarchy_cells_shuffle, ...
        gad_hierarchy_cells_shuffle,'uni',false);
i
end

r = cellfun(@nanmean,all_hierarchy_cells_shuffle);
j = (gad_hierarchy_cells-r)./(gad_num_active - r);
figure;errorbar(nanmean(j),nanstd(j)./sqrt(sum(~isnan(j))),'k')


% this was to get correlation of ranks, i.e. same hierarchy order, but I'm
% not sure how to execute this and just take daily movecells into acct
% usem = cellfun(@(x) ~isempty(x),pyr_class_all);
% pyr_movecell_mean = cell(size(pyr_class_all));
% pyr_move_mean = cellfun(@(x) nanmean(x,2),pyr_move_active,'uni',false);
% pyr_movecell_mean(usem) = cellfun(@(x,y) x.*y',pyr_move_mean(usem),pyr_class_all(usem),'uni',false);
% pyr_movecell_rank = cellfun(@(x) tiedrank(x),pyr_movecell_mean,'uni',false);
% pyr_movecell_rank_cc = nan(15,15,8);
% for i = 1:8
%    curr_data = horzcat(pyr_movecell_rank{i,:});
%    curr_cc = corrcoef(curr_data);
%    pyr_movecell_rank_cc(usem(i,:),usem(i,:),i) = curr_cc;
% end
% figure;imagesc(nanmean(pyr_movecell_rank_cc(1:14,1:14,:),3));colormap(hot)
% title('pyr rank correlation')
% diag_data_1 = arrayfun(@(x) diag(pyr_movecell_rank_cc(:,:,x),-1),...
%     1:size(pyr_movecell_rank_cc,3),'uni',false); 
% diag_data_2 = arrayfun(@(x) diag(pyr_movecell_rank_cc(:,:,x),-2),...
%     1:size(pyr_movecell_rank_cc,3),'uni',false); 
% diag_data_3 = arrayfun(@(x) diag(pyr_movecell_rank_cc(:,:,x),-3),...
%     1:size(pyr_movecell_rank_cc,3),'uni',false); 
% diag_data_combine = cellfun(@(x,y,z) [x(1:12)';y(1:12)';z(1:12)'],diag_data_1, ...
%     diag_data_2,diag_data_3,'uni',false)';
% diag_data_concat = vertcat(diag_data_combine{:});
% figure;errorbar(nanmean(diag_data_concat), ...
%     nanstd(diag_data_concat)./sqrt(sum(~isnan(diag_data_concat))),'k');
% title('pyr rank correlation (day+1:day+3)')
% 
% 
% usem = cellfun(@(x) ~isempty(x),gad_class_all);
% gad_movecell_mean = cell(size(gad_class_all));
% gad_move_mean = cellfun(@(x) nanmean(x,2),gad_move_active,'uni',false);
% gad_movecell_mean(usem) = cellfun(@(x,y) x.*y',gad_move_mean(usem),gad_class_all(usem),'uni',false);
% gad_movecell_rank = cellfun(@(x) tiedrank(x),gad_movecell_mean,'uni',false);
% gad_movecell_rank_cc = nan(15,15,8);
% for i = 1:8
%    curr_data = horzcat(gad_movecell_rank{i,:});
%    curr_cc = corrcoef(curr_data);
%    gad_movecell_rank_cc(usem(i,:),usem(i,:),i) = curr_cc;
% end
% figure;imagesc(nanmean(gad_movecell_rank_cc(1:14,1:14,:),3));colormap(hot)
% title('gad rank correlation')
% diag_data_1 = arrayfun(@(x) diag(gad_movecell_rank_cc(:,:,x),-1),...
%     1:size(gad_movecell_rank_cc,3),'uni',false); 
% diag_data_2 = arrayfun(@(x) diag(gad_movecell_rank_cc(:,:,x),-2),...
%     1:size(gad_movecell_rank_cc,3),'uni',false); 
% diag_data_3 = arrayfun(@(x) diag(gad_movecell_rank_cc(:,:,x),-3),...
%     1:size(gad_movecell_rank_cc,3),'uni',false); 
% diag_data_combine = cellfun(@(x,y,z) [x(1:12)';y(1:12)';z(1:12)'],diag_data_1, ...
%     diag_data_2,diag_data_3,'uni',false)';
% diag_data_concat = vertcat(diag_data_combine{:});
% figure;errorbar(nanmean(diag_data_concat), ...
%     nanstd(diag_data_concat)./sqrt(sum(~isnan(diag_data_concat))),'k');
% title('gad rank correlation (day+1:day+3)')


%% HIERARCHY: MOST UP-TO DATE
% DO THIS: 
% shuffle x
% hierarchy score by mean of cell hierarchy scores

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,'uni',false);
% pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,'uni',false);

% this is to use only the movement-classified ones
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);

%%%% FOR GAD
gad_move_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    gad_move_active,'uni',false);

gad_hierarchy_events = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',gad_move_sort,'uni',false);

gad_active_events = cellfun(@(x) sum(x(1:end-1,:),2),gad_move_sort,'uni',false);


% and try also getting chance by shuffling active cells on a trial
rep = 100;
all_hierarchy_events_shuffle = cell(size(gad_move_active));
for i = 1:rep
    active_move_shuffle = cellfun(@(x) shake(x,2),gad_move_active,'uni',false);
    
    active_move_shuffle_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    active_move_shuffle,'uni',false);

    hierarchy_events_shuffle = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',active_move_shuffle_sort,'uni',false);

    all_hierarchy_events_shuffle = cellfun(@(x,y) [x y],all_hierarchy_events_shuffle, ...
        hierarchy_events_shuffle,'uni',false);
i
end

mean_hierarchy_shuffle = cellfun(@(x) nanmean(x,2),all_hierarchy_events_shuffle,'uni',false);
% hierarchy score by mean of all individual cells hierarchy scores
hierarchy_score = cellfun(@(x,y,z) nanmean((x-y)./(z - y)), ...
    gad_hierarchy_events,mean_hierarchy_shuffle,gad_active_events);

% hierarchy score by looking at total hierarchy events
hierarchy_score_sum = cellfun(@(x,y,z) (sum(x)-sum(y))/(sum(z) - sum(y)), ...
    gad_hierarchy_events,mean_hierarchy_shuffle,gad_active_events);

figure; hold on
errorbar(nanmean(hierarchy_score), ...
    nanstd(hierarchy_score)./sqrt(sum(~isnan(hierarchy_score))),'k')
errorbar(nanmean(hierarchy_score_sum), ...
    nanstd(hierarchy_score_sum)./sqrt(sum(~isnan(hierarchy_score_sum))),'r')
ylabel('Hierarchy Score')
xlabel('Session')
legend({'Mean' 'Sum'});
title('Gad hierarchy score')

%%%% FOR PYR

pyr_move_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    pyr_move_active,'uni',false);

pyr_hierarchy_events = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',pyr_move_sort,'uni',false);

pyr_active_events = cellfun(@(x) sum(x(1:end-1,:),2),pyr_move_sort,'uni',false);


% and try also getting chance by shuffling active cells on a trial
rep = 100;
all_hierarchy_events_shuffle = cell(size(pyr_move_active));
for i = 1:rep
    active_move_shuffle = cellfun(@(x) shake(x,2),pyr_move_active,'uni',false);
    
    active_move_shuffle_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    active_move_shuffle,'uni',false);

    hierarchy_events_shuffle = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',active_move_shuffle_sort,'uni',false);

    all_hierarchy_events_shuffle = cellfun(@(x,y) [x y],all_hierarchy_events_shuffle, ...
        hierarchy_events_shuffle,'uni',false);
i
end

mean_hierarchy_shuffle = cellfun(@(x) nanmean(x,2),all_hierarchy_events_shuffle,'uni',false);
% hierarchy score by mean of all individual cells hierarchy scores
hierarchy_score = cellfun(@(x,y,z) nanmean((x-y)./(z - y)), ...
    pyr_hierarchy_events,mean_hierarchy_shuffle,pyr_active_events);

% hierarchy score by looking at total hierarchy events
hierarchy_score_sum = cellfun(@(x,y,z) (sum(x)-sum(y))/(sum(z) - sum(y)), ...
    pyr_hierarchy_events,mean_hierarchy_shuffle,pyr_active_events);

figure; hold on
errorbar(nanmean(hierarchy_score), ...
    nanstd(hierarchy_score)./sqrt(sum(~isnan(hierarchy_score))),'k')
errorbar(nanmean(hierarchy_score_sum), ...
    nanstd(hierarchy_score_sum)./sqrt(sum(~isnan(hierarchy_score_sum))),'r')
ylabel('Hierarchy Score')
xlabel('Session')
legend({'Mean' 'Sum'});
title('Pyr hierarchy score')

% how to find chance of having given number of correct hierarchy cells in a
% trial:
%
% (nchoosek(active_events,hierarchy_events) * ...
%     nchoosek(num_events-active_events,active_events-hierarchy_events)) / ...
%       nchoosek(num_events,active_events);


  

% how to find chance of having given number of correct hierarchy cells in a
% trial:
%
% (nchoosek(hierarchy_group,hierarchy_cells) * ...
%     nchoosek(num_cells-hierarchy_group,active_cells-hierarchy_cells)) / ...
%       nchoosek(num_cells,active_cells);



% BY PROBABILITY INSTEAD OF SHUFFLING

% %%%% WARNING! THIS GIVES VERY LOW PROBS... I THINK SOMETHING WRONG
% % get the chance of having active cells be the most active ones
% gad_hierarchy_chance = cellfun(@(x,y) arrayfun(@(z) ...
%     (nchoosek(sum(x(z,:)),y(z)) * ...
%     nchoosek(size(x,2)-sum(x(z,:)),sum(x(z,:))-y(z))) / ...
%       nchoosek(size(x,2),sum(x(z,:))),1:length(y)), ...
%       gad_move_active,gad_hierarchy_events,'uni',false);
%   





% Compare data to simulated data: does correlation produce hierarchy?
use_animalday = cellfun(@(x) ~isempty(x),gad_move_active);
gad_move_simulated = cell(size(use_animalday));
for i = find(use_animalday)'
    curr_data = gad_move_active{i};
    exclude_data = nanmean(curr_data,2) == 0 | nanmean(curr_data,2) == 1;
    curr_data = curr_data(~exclude_data,:);
    
    mean_corrcoef = nanmean(AP_itril(corrcoef(curr_data'),-1))/2;
    mean_corrcoef_matrix = ones(size(curr_data,1))*mean_corrcoef;
    mean_corrcoef_matrix = ones(size(curr_data,1));
    mean_corrcoef_matrix(eye(size(mean_corrcoef_matrix))>0) = 1;
    std_cov = (std(curr_data,[],2)*std(curr_data,[],2)');
    mean_cov_matrix = mean_corrcoef_matrix.*std_cov;
    
    activity_mean = repmat(nanmean(curr_data(:)),size(curr_data,1),1)/4;
    
    warning off
    %sim_data = sampleDichGauss01( ...
    %    nanmean(curr_data,2),cov(curr_data'),size(curr_data,2));
    sim_data = sampleDichGauss01( ...
        activity_mean,cov(curr_data'),size(curr_data,2));
    warning on
    sim_data_cat = [sim_data;gad_move_active{i}(exclude_data,:)];
    gad_move_simulated{i} = sim_data_cat;
    i
end

gad_move_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    gad_move_simulated,'uni',false);

gad_hierarchy_events = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',gad_move_sort,'uni',false);

gad_active_events = cellfun(@(x) sum(x(1:end-1,:),2),gad_move_sort,'uni',false);


% and try also getting chance by shuffling active cells on a trial
rep = 100;
all_hierarchy_events_shuffle = cell(size(gad_move_simulated));
for i = 1:rep
    active_move_shuffle = cellfun(@(x) shake(x,2),gad_move_simulated,'uni',false);
    
    active_move_shuffle_sort = cellfun(@(x) sortrows([x' nanmean(x',2)],-(size(x,1)+1))', ...
    active_move_shuffle,'uni',false);

    hierarchy_events_shuffle = cellfun(@(x) arrayfun(@(y) ...
    sum(x(y,1:sum(x(y,:)))),1:size(x,1)-1)',active_move_shuffle_sort,'uni',false);

    all_hierarchy_events_shuffle = cellfun(@(x,y) [x y],all_hierarchy_events_shuffle, ...
        hierarchy_events_shuffle,'uni',false);
i
end

mean_hierarchy_shuffle = cellfun(@(x) nanmean(x,2),all_hierarchy_events_shuffle,'uni',false);
% hierarchy score by mean of all individual cells hierarchy scores
hierarchy_score_sim = cellfun(@(x,y,z) nanmean((x-y)./(z - y)), ...
    gad_hierarchy_events,mean_hierarchy_shuffle,gad_active_events);





%% gad-pyr distance relationship (just needs correlation - pairwise)

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
subset_corr = cell(size(all_corr));
reliable_frac = nan(size(all_corr));
gadradius_maps = cell(size(all_corr));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % load movement classification
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    
    sessions = find(cellfun(@(x) ~isempty(x),all_corr(curr_animal,:)));
    
    %figure;
    for curr_day = sessions
                
        
        
        curr_reliabilities = nanmean(trial_activity.peak{curr_animal}.cued.moving{curr_day},2);
        
        gad_unfilled_move = gad_unfilled(...
            classified_cells.move_cells_peak{curr_day}(gad_unfilled));
        pyr_unfilled_move = pyr_cells(...
            classified_cells.move_cells_peak{curr_day}(pyr_unfilled));
        
        %%% TEMPORARY
        curr_gad = gad_unfilled_move(curr_reliabilities(gad_unfilled_move) > 0.9 ...
            & curr_reliabilities(gad_unfilled_move) <= 1);
        reliable_frac(curr_animal,curr_day) = length(curr_gad);
        %%%
        
        gad_move = gad_cells(...
            classified_cells.move_cells_peak{curr_day}(gad_cells));
        pyr_move = pyr_cells(...
            classified_cells.move_cells_peak{curr_day}(pyr_cells));
        
        radius_bin_edges = 1:50:501;
        
        %curr_gad = gad_unfilled_move(curr_reliabilities(gad_unfilled_move) > 0.9 ...
        %    & curr_reliabilities(gad_unfilled_move) <= 1);
        curr_gad = gad_unfilled_move;
        %curr_gad = pyr_move;
        curr_pyr = pyr_unfilled_move;

        
        if isempty(curr_pyr) || isempty(curr_gad)
            continue
        end
        
        num_bins = length(radius_bin_edges);
        
        gadpyr_corr = all_corr{curr_animal,curr_day}(curr_gad,curr_pyr);
        gadpyr_dist = roi_dist(curr_gad,curr_pyr);
        
        [n,gadpyr_bin] = histc(gadpyr_dist(:),radius_bin_edges);
        %gadpyr_binMean = grpstats(gadpyr_corr(:),gadpyr_bin);
        %gadpyr_binSem = grpstats(gadpyr_corr(:),gadpyr_bin,'sem');
        
        % count each cell as an n
        gadpyr_binmembers = cell(num_bins,1);
        for i = 1:num_bins
            curr_binmembers = gadpyr_corr(gadpyr_bin == i);
            gadpyr_binmembers{i} = vertcat(curr_binmembers(:));
        end
                
        subset_corr{curr_animal,curr_day} = gadpyr_binmembers;
        
        % for plotting gad cells by reliability
        
        %         subplot(4,4,curr_day); hold on
        %         for i = curr_gad
        %             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2), ...
        %                 [curr_reliabilities(i) curr_reliabilities(i) curr_reliabilities(i)]);
        %             ylim([-512 0]);
        %             xlim([0 512]);
        %         end
        
        % For plotting gad cells with mean of pyr correlations
        
%          if curr_animal == 2 && curr_day == 3;
%              keyboard
%          end
%         subplot(4,4,curr_day); hold on
%         for i = curr_gad
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2), ...
%                 [curr_reliabilities(i) 0 0]);
%         end
%         for i = 1:length(curr_pyr)
%             pyr = curr_pyr(i);
%             max_corr = max(nanmean(gadpyr_corr,1));
%             min_corr = min(nanmean(gadpyr_corr,1));
%             curr_corr = (nanmean(gadpyr_corr(:,i)) - min_corr) / ...
%                 (max_corr - min_corr);
%             if isnan(curr_corr)
%                 patch(polygon.ROI{pyr}(:,1),-polygon.ROI{pyr}(:,2), ...
%                     [1 1 1]);
%                 continue
%             end
%             patch(polygon.ROI{pyr}(:,1),-polygon.ROI{pyr}(:,2), ...
%                 [0 curr_corr 0]);
%         end
%         ylim([-512 0]);
%         xlim([0 512]);
        

%         % TEMPORARY: this is to visually identify high gad-pyr corr dist
%         if curr_animal == 6 && curr_day == 9;
%             keyboard
%         end
%         asdf = 0;
%         if asdf
%             figure; hold on
%             for i = curr_gad
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2), ...
%                     [curr_reliabilities(i) 0 0]);
%             end
%             ylim([-512 0]);
%             xlim([0 512]);
% 
%             for gad = 1:length(curr_gad)
%             figure; hold on
%             patch(polygon.ROI{curr_gad(gad)}(:,1),-polygon.ROI{curr_gad(gad)}(:,2), ...
%                     [1 curr_reliabilities(i) 0]);
%             for i = 1:length(curr_pyr)
%                 pyr = curr_pyr(i);
%                 %for i = setdiff(pyr_cells,filled);
%                 temp_corr = (gadpyr_corr(gad,i)-min(gadpyr_corr(gad,:)))/...
%                     (max(gadpyr_corr(gad,:))-min(gadpyr_corr(gad,:)));
%                 if isnan(temp_corr)
%                     patch(polygon.ROI{pyr}(:,1),-polygon.ROI{pyr}(:,2), ...
%                     [1 1 1]);
%                 continue
%                 end
%                 patch(polygon.ROI{pyr}(:,1),-polygon.ROI{pyr}(:,2), ...
%                     [0 temp_corr 0]);
%                 %[max_idx(i)/max(max_idx) max_peaks(i)/max(max_peaks) 0]);
%             end
%             ylim([-512 0]);
%             xlim([0 512]);
%             end
%         end
        
        
%         % Show space map of pyr cells correlated with gad cell
%         subplot(4,4,curr_day); hold on
%         space_matrix = nan(512,512);
%         spread_filter = fspecial('gaussian',512,30);
%         for i = 1:length(pyr_move)
%             pyr = pyr_move(i);
%             curr_space_matrix = zeros(512,512);
%             curr_x = round(roi_center_x(pyr));
%             curr_y = round(roi_center_y(pyr));
%             if curr_x < 1
%                 curr_x = 1;
%             end
%             if curr_y < 1
%                 curr_y = 1;
%             end          
%             if curr_x > 512
%                 curr_x = 512;
%             end
%             if curr_y > 512
%                 curr_y = 512;
%             end     
%             %curr_space_matrix(curr_y,curr_x) = nanmean(gadpyr_corr(:,i));
%             curr_space_matrix(curr_y,curr_x) = curr_reliabilities(pyr);
%             curr_space_matrix_spread = conv2(curr_space_matrix,...
%                 spread_filter,'same');
%             
%             space_matrix = nansum(cat(3,space_matrix, ...
%                 curr_space_matrix_spread),3);
%         end
%         imagesc(space_matrix);colormap(gray)
%         ylim([0 512])
%         xlim([0 512])
        


%         % Show space map of pyramidal move cells
%         curr_img = subplot(4,4,curr_day);
%         spread_filter = fspecial('gaussian',512,50);
%         
%         pyr_space_matrix = nan(512,512);
%         for i = 1:length(curr_pyr);
%             pyr = curr_pyr(i);
%             curr_space_matrix = zeros(512,512);
%             curr_x = round(roi_center_x(pyr));
%             curr_y = round(roi_center_y(pyr));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             %curr_space_matrix(curr_y,curr_x) = nanmean(gadpyr_corr(:,i));
%             curr_space_matrix(curr_y,curr_x) = curr_reliabilities(pyr);
%             curr_space_matrix_spread = conv2(curr_space_matrix,...
%                 spread_filter,'same');          
%             pyr_space_matrix = nansum(cat(3,pyr_space_matrix, ...
%                 curr_space_matrix_spread),3);
%         end
%         pyr_space_matrix_norm = (pyr_space_matrix - min(pyr_space_matrix(:))) ./ ...
%             (max(pyr_space_matrix(:)) - min(pyr_space_matrix(:)));
%         
%         gad_space_matrix = nan(512,512);
%         for gad = gad_unfilled_move;
%             curr_space_matrix = zeros(512,512);
%             curr_x = round(roi_center_x(gad));
%             curr_y = round(roi_center_y(gad));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             curr_space_matrix(curr_y,curr_x) = curr_reliabilities(gad);
%             curr_space_matrix_spread = conv2(curr_space_matrix,...
%                 spread_filter,'same');          
%             gad_space_matrix = nansum(cat(3,gad_space_matrix, ...
%                 curr_space_matrix_spread),3);
%         end
%         gad_space_matrix_norm = (gad_space_matrix - min(gad_space_matrix(:))) ./ ...
%             (max(gad_space_matrix(:)) - min(gad_space_matrix(:)));
%         
%         
%         pyrgad_space_matrix = gad_space_matrix_norm;
%         pyrgad_space_matrix(:,:,2) = pyr_space_matrix_norm;
%         pyrgad_space_matrix(:,:,3) = zeros(size(gad_space_matrix));
%         
%         imagesc(pyrgad_space_matrix);colormap(gray)
%         ylim([0 512])
%         xlim([0 512])
%         set(curr_img,'YDir','reverse')


%         if curr_animal == 6 && curr_day == 9;
%             keyboard
%         end
      %  curr_img = subplot(4,4,curr_day);
        %spread_filter = fspecial('disk',20);
        %spread_filter = fspecial('gaussian',512,50);
        
        gadpyr_space_matrix = nan(1024,1024);
        for i = 1:length(curr_gad);
            gad = curr_gad(i);
            gad_x = round(roi_center_x(gad));
            gad_y = round(roi_center_y(gad));
            for j = 1:length(curr_pyr);
                pyr = curr_pyr(j);
                curr_x = round(roi_center_x(pyr)) - gad_x + 512;
                curr_y = round(roi_center_y(pyr)) - gad_y + 512;
                if curr_x < 1;curr_x = 1;end
                if curr_y < 1;curr_y = 1;end
                if curr_x > 1024;curr_x = 1024;end
                if curr_y > 1024;curr_y = 1024;end
                if isnan(gadpyr_corr(i,j))
                    continue
                end
                gadpyr_space_matrix(curr_y,curr_x) = nansum(...
                    [gadpyr_space_matrix(curr_y,curr_x) gadpyr_corr(i,j)]);
            end
        end
        %curr_space_matrix_spread = conv2(curr_gadpyr_space_matrix,...
        %        spread_filter,'same');
        %gadpyr_space_matrix_norm = (gadpyr_space_matrix - min(gadpyr_space_matrix(:))) ./ ...
        %    (max(gadpyr_space_matrix(:)) - min(gadpyr_space_matrix(:)));
        gadpyr_space_matrix_norm = gadpyr_space_matrix;
        
%         imagesc(gadpyr_space_matrix_norm);colormap(gray)
%         ylim([0 1024])
%         xlim([0 1024])
%         set(curr_img,'YDir','normal')

        gadradius_maps{curr_animal,curr_day} = gadpyr_space_matrix_norm;
        
        curr_day
    end
    
    
    curr_animal
end

%For condition 2 (pairwise)
use_cells = cellfun(@(x) ~isempty(x),subset_corr);
radius_concat = cell(num_bins,1);
for i = 1:num_bins
   curr_data = cellfun(@(x) x{i},subset_corr(use_cells),'uni',false); 
   radius_concat{i} = vertcat(curr_data{:});
end
data_mean = cellfun(@nanmean,radius_concat);
data_sem = cellfun(@(x) nanstd(x)/sum(~isnan(x)),radius_concat);
figure;errorbar(data_mean,data_sem,'k','linewidth',2);

figure; hold on
use_cells = cellfun(@(x) ~isempty(x),subset_corr);
temp = cell(size(use_cells));
temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
l = horzcat(temp{:})';
errorbar(nanmean(l),nanstd(l)./sqrt(sum(~isnan(l))),'k','linewidth',2)
col = jet(8);
for i = 1:8
    use_cells = cellfun(@(x) ~isempty(x),subset_corr);
    use_cells(i,:) = false;
    temp = cell(size(use_cells));
    temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
    l = horzcat(temp{:})';
    errorbar(nanmean(l),nanstd(l)./sqrt(sum(~isnan(l))),'color',col(i,:),'linewidth',2)
end
title('exclude animals')

figure; hold on
use_cells = cellfun(@(x) ~isempty(x),subset_corr);
temp = cell(size(use_cells));
temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
l = horzcat(temp{:})';
errorbar(nanmean(l),nanstd(l)./sqrt(sum(~isnan(l))),'k','linewidth',2)
col = jet(15);
for i = 1:14
    use_cells = cellfun(@(x) ~isempty(x),subset_corr);
    use_cells(:,i) = false;
    temp = cell(size(use_cells));
    temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
    l = horzcat(temp{:})';
    errorbar(nanmean(l),nanstd(l)./sqrt(sum(~isnan(l))),'color',col(i,:),'linewidth',2)
end
title('exclude days')

figure; hold on
use_cells = cellfun(@(x) ~isempty(x),subset_corr);
temp = cell(size(use_cells));
temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
l = horzcat(temp{:})';
errorbar(nanmedian(l),nanstd(l)./sqrt(sum(~isnan(l))),'k','linewidth',2)
col = jet(15);
for i = 1:14
    use_cells = cellfun(@(x) ~isempty(x),subset_corr);
    null_days = setdiff(1:15,i);
    use_cells(:,null_days) = false;
    temp = cell(size(use_cells));
    temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
    l = horzcat(temp{:})';
    subplot(4,4,i)
    errorbar(nanmean(l),nanstd(l)./sqrt(sum(~isnan(l))),'color',col(i,:),'linewidth',2)
end
title('individual days')




use_cells = cellfun(@(x) ~isempty(x),subset_corr);
temp = cell(size(use_cells));
temp(use_cells) = cellfun(@(x) cellfun(@nanmean,x),subset_corr(use_cells),'uni',false);
figure; hold on
day_combine = {[1 2 3] [4 5 6] [7 8 9] [10 11 12] [13 14]};
col = jet(length(day_combine));
for i = 1:length(day_combine)
   curr_data = horzcat(temp{:,day_combine{i}})'; 
   errorbar(nanmean(curr_data),nanstd(curr_data)./ ...
       sqrt(sum(~isnan(curr_data))),'color',col(i,:),'linewidth',2) 
end





sum_map = nan(1024,1024);
num_map = zeros(1024,1024);
cat_maps = cat(3,gadradius_maps{:});
num_map = sum(~isnan(cat_maps),3);
sum_map = nansum(cat_maps,3);
sum_map(num_map == 0) = NaN;
    


%day_combine = {[1 2 3] [4 5 6] [7 8 9] [10 11 12] [13 14]};
day_combine = num2cell(1:14);
all_sum_maps = cell(length(day_combine),1);
all_num_maps = cell(length(day_combine),1);
for i = 1:length(day_combine);
    curr_days = day_combine{i};
    sum_map = nan(1024,1024);
    num_map = zeros(1024,1024);
    cat_maps = cat(3,gadradius_maps{:,curr_days});
    num_map = sum(~isnan(cat_maps),3);
    sum_map = nansum(cat_maps,3);
    sum_map(num_map == 0) = NaN;
    
    all_sum_maps{i} = sum_map;
    all_num_maps{i} = num_map;
    i
end
for i = 1:length(all_sum_maps)
    j = all_sum_maps{i};
    k = all_num_maps{i};
    %j = nanmean(sum_map,3);
    %k = nanmean(num_map,3);
    j(isnan(j)) = 0;
    r = conv2(j,spread_filter,'same');
    r2 = conv2(k,spread_filter,'same');
    m = r./r2;
    figure;
    imagesc(m);colormap(gray);
end


[rr cc] = meshgrid(1:1024);
dists = radius_bin_edges;
dist_mean = nan(length(dists)-1,8);
for j = 1:8
    for i = 1:length(dists)-1
        c = sqrt((rr - 512).^2+(cc-512).^2) >= dists(i) & ...
            sqrt((rr - 512).^2+(cc-512).^2) < dists(i+1);
        dist_mean(i,j) = nansum(all_sum_maps{j}(c))/nansum(all_num_maps{j}(c));
    end
end
figure;plot(dist_mean);

figure;imagesc(sum_map);colormap(gray);
set(gca,'XTick',[0:128:1024])
set(gca,'XTickLabel',[-512:128:512])
set(gca,'YTick',[0:128:1024])
set(gca,'YTickLabel',[-512:128:512])

[rr cc] = meshgrid(1:1024);
dists = radius_bin_edges;
dist_mean = nan(length(dists)-1,1);
for i = 1:length(dists)-1
    c = sqrt((rr - 512).^2+(cc-512).^2) >= dists(i) & ...
        sqrt((rr - 512).^2+(cc-512).^2) < dists(i+1);   
    %dist_mean(i) = nansum(sum_map(c))/nansum(num_map(c));
    dist_mean(i) = nanmean(m(c));
end
figure;plot(dist_mean);




% averaging it this way doesn't make it look better
%day_combine = {[1 2 3] [4 5 6] [7 8 9] [10 11 12] [13 14]};
day_combine = num2cell(1:14);
for i = 1:length(day_combine);
    curr_days = day_combine{i};
    
    sum_map = nan(1024,1024);
    num_map = zeros(1024,1024);
    cat_maps = cat(3,gadradius_maps{:,curr_days});
    num_map = sum(~isnan(cat_maps),3);
    sum_map = nansum(cat_maps,3);
    sum_map(num_map == 0) = NaN;
    
%     nonan_maps = cat_maps;
%     nonan_maps(isnan(nonan_maps)) = 0;
%     cat_maps_smooth = arrayfun(@(x) ...
%     conv2(nonan_maps(:,:,x),spread_filter,'same'),1:size(cat_maps,3),'uni',false);
%     cat_maps_smooth = cat(3,cat_maps_smooth{:});
%     
%     num_maps_smooth = arrayfun(@(x) ...
%     conv2(+[nonan_maps(:,:,x) > 0],spread_filter,'same'),1:size(cat_maps,3),'uni',false);
%     num_maps_smooth = cat(3,num_maps_smooth{:});
    
    
    j = sum_map;
    k = num_map;
    %j = nanmean(sum_map,3);
    %k = nanmean(num_map,3);
    j(isnan(j)) = 0;
    r = conv2(j,spread_filter,'same');
    r2 = conv2(k,spread_filter,'same');
    m = r./r2;
    
    figure;imagesc(m);colormap(gray);
   % caxis([0 0.2]);
end



%% Gad-pyr distance assorted ways (requires activity in memory)

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
radius_corr_all = cell(size(gad_activity_thresh_all));
num_pyrcells_all = cell(size(gad_activity_thresh_all));
num_pyrcorr_all = cell(size(gad_activity_thresh_all));
mean_radius_all = cell(size(gad_activity_thresh_all));

dist_weight_all = cell(size(gad_activity_thresh_all));
field_lscorr_all = cell(size(gad_activity_thresh_all));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % load movement classification
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);
    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    
    sessions = find(cellfun(@(x) ~isempty(x),gad_activity_thresh_all(curr_animal,:)));
    
    
    for curr_day = sessions
        
%         if curr_animal == 5 && curr_day == 3;
%             keyboard
%         end
        
        gad_unfilled_move = classified_cells.move_cells_peak{curr_day}(gad_unfilled);
        pyr_unfilled_move = classified_cells.move_cells_peak{curr_day}(pyr_unfilled);
        
        gad_unfilled_move_idx = gad_unfilled(gad_unfilled_move);
        pyr_unfilled_move_idx = pyr_unfilled(pyr_unfilled_move);
        
        curr_gad_act = gad_activity_thresh_all{curr_animal,curr_day}(gad_unfilled_move,:);

        % pick unfilled move-related or all cells
        curr_dist = roi_dist(gad_unfilled_move_idx,pyr_unfilled_move_idx);
        curr_pyr_act = pyr_activity_all{curr_animal,curr_day}(pyr_unfilled_move,:);
        %curr_dist = roi_dist(gad_unfilled_move_idx,pyr_unfilled);
        %curr_pyr_act = pyr_activity_all{curr_animal,curr_day};
        
        % Find correlation between gad cells and summed activity within
        % ring around gad cells at different radii
        
        % Spread pyramidal activity out to get better corr measure
        spread_filter = ones(1,30);
        pyr_spread = conv2(curr_pyr_act,spread_filter);
        pyr_spread = pyr_spread(:,1:size(curr_pyr_act,2));
        pyr_spread(isnan(pyr_spread)) = 0;
 
        
%         gad_activity_celld = mat2cell(curr_gad_act,...
%             ones(size(curr_gad_act,1),1),size(curr_gad_act,2));
%         radius_bin_edges = 1:50:600;
%         %radius_bin_edges = [1 120 180 230 290 370 600];
%         radius_corr = nan(length(radius_bin_edges)-1,sum(gad_unfilled_move));
%         num_pyrcells = nan(length(radius_bin_edges)-1,sum(gad_unfilled_move));
%         for radius_idx = 1:length(radius_bin_edges) -1;
%             curr_gad_rad_idx = curr_dist >= radius_bin_edges(radius_idx) ...
%                 & curr_dist < radius_bin_edges(radius_idx+1);
%             curr_pyr_act = arrayfun(@(x) ...
%                 sum(pyr_spread(curr_gad_rad_idx ...
%                 (x,:),:),1),1:sum(gad_unfilled_move),'uni',false)';
%             curr_radius_corr = cellfun(@(x,y) corrcoef(x,y), ...
%                 gad_activity_celld,curr_pyr_act,'uni',false);
%             radius_corr(radius_idx,:) = cellfun(@(x) x(2),curr_radius_corr);
%             
%             num_pyrcells(radius_idx,:) = sum(curr_gad_rad_idx,2);          
%         end
        

% 
%         % Max correlation of least squares pyr linsum
%         gad_activity_celld = mat2cell(curr_gad_act,...
%             ones(size(curr_gad_act,1),1),size(curr_gad_act,2));
%         radius_bin_edges = 1:50:600;
%         %radius_bin_edges = [1 120 180 230 290 370 600];
%         radius_corr = nan(length(radius_bin_edges)-1,sum(gad_unfilled_move));
%         num_pyrcells = nan(length(radius_bin_edges)-1,sum(gad_unfilled_move));
%         for radius_idx = 1:length(radius_bin_edges) -1;
%             % get usable pyr cells
%             
%             % ring
%             curr_gad_rad_idx = curr_dist >= radius_bin_edges(radius_idx) ...
%                 & curr_dist < radius_bin_edges(radius_idx+1);
%             
%             % circle
%             %curr_gad_rad_idx = curr_dist < radius_bin_edges(radius_idx+1);
%             % get individual correlations
%             
%             curr_corr = corrcoef([curr_gad_act;pyr_spread]');
%             curr_corr_subset = curr_corr(1:size(curr_gad_act,1),...
%                 size(curr_gad_act,1)+1:end);
%             curr_idx_corr = curr_gad_rad_idx.*curr_corr_subset;
%             
%             warning off
%             pyr_weights = arrayfun(@(x) ...
%                 pyr_spread(curr_gad_rad_idx(x,:),:)'\curr_gad_act(x,:)', ...
%                 1:sum(gad_unfilled_move),'uni',false)';
%             warning on 
%             
%             curr_pyr_linsum = arrayfun(@(x) ...
%                 pyr_spread(curr_gad_rad_idx(x,:),:)'*pyr_weights{x}, ...
%                 1:sum(gad_unfilled_move),'uni',false)';
%             
%             curr_radius_corr = cellfun(@(x,y) corrcoef(x,y), ...
%                 gad_activity_celld,curr_pyr_linsum,'uni',false);
%             
%             radius_corr(radius_idx,:) = cellfun(@(x) x(2),curr_radius_corr);
%             
%             num_pyrcells(radius_idx,:) = sum(curr_gad_rad_idx,2);
%         end
%         radius_corr_all{curr_animal,curr_day} = radius_corr;
%         num_pyrcells_all{curr_animal,curr_day} = num_pyrcells;
%         mean_radius_all{curr_animal,curr_day} = mean_radius;

        
%         % Max correlation of least squares pyr linsum: rings by num cells
%         gad_activity_celld = mat2cell(curr_gad_act,...
%             ones(size(curr_gad_act,1),1),size(curr_gad_act,2));
%         radius_cells = 10;
%         radius_cell_groups = floor(size(curr_dist,2)/radius_cells);
%         radius_corr = nan(radius_cell_groups,sum(gad_unfilled_move));
%         num_pyrcells = nan(radius_cell_groups,sum(gad_unfilled_move));
%         mean_radius = nan(radius_cell_groups,sum(gad_unfilled_move));
%         for radius_idx = 1:radius_cell_groups;
%             % get usable pyr cells
%             [asdf curr_dist_rank] = sort(curr_dist,2);
%             
%             % ring
%             curr_gad_rad_idx = curr_dist_rank > radius_cells*(radius_idx-1) ...
%                 &  curr_dist_rank <= radius_cells*(radius_idx);
%             
%             curr_corr = corrcoef([curr_gad_act;pyr_spread]');
%             curr_corr_subset = curr_corr(1:size(curr_gad_act,1),...
%                 size(curr_gad_act,1)+1:end);
%             curr_idx_corr = curr_gad_rad_idx.*curr_corr_subset;           
%             
%             pyr_weights = arrayfun(@(x) ...
%                 pyr_spread(curr_gad_rad_idx(x,:),:)'\curr_gad_act(x,:)', ...
%                 1:sum(gad_unfilled_move),'uni',false)';
%             
%             curr_pyr_linsum = arrayfun(@(x) ...
%                 pyr_spread(curr_gad_rad_idx(x,:),:)'*pyr_weights{x}, ...
%                 1:sum(gad_unfilled_move),'uni',false)';
%             
%             curr_radius_corr = cellfun(@(x,y) corrcoef(x,y), ...
%                 gad_activity_celld,curr_pyr_linsum,'uni',false);
%             
%             radius_corr(radius_idx,:) = cellfun(@(x) x(2),curr_radius_corr);          
%             num_pyrcells(radius_idx,:) = sum(curr_gad_rad_idx,2);
%             
%             curr_usedist = curr_dist;
%             curr_usedist(~curr_gad_rad_idx) = nan;
%             mean_radius(radius_idx,:) = nanmean(curr_usedist,2);
%             
%         end      
%         radius_corr_all{curr_animal,curr_day} = radius_corr;
%         num_pyrcells_all{curr_animal,curr_day} = num_pyrcells;
%         mean_radius_all{curr_animal,curr_day} = mean_radius;


        % whole field: LS weight vs. distance      
        %pyr_spread = +[pyr_spread > 0];
        %curr_gad_act = +[curr_gad_act > 0];
        gad_activity_celld = mat2cell(curr_gad_act,...
            ones(size(curr_gad_act,1),1),size(curr_gad_act,2));
        dist_weight = nan(size(curr_dist,1),size(curr_dist,2),2);
        field_lscorr = nan(size(curr_gad_act,1),1);
        
        pyr_weights = pyr_spread'\curr_gad_act';
        
        curr_pyr_linsum = transpose(pyr_spread'*pyr_weights);
        curr_pyr_linsum_celld = mat2cell(curr_pyr_linsum,...
            ones(size(curr_pyr_linsum,1),1),size(curr_pyr_linsum,2));
        
        curr_radius_corr = cellfun(@(x,y) corrcoef(x,y), ...
            gad_activity_celld,curr_pyr_linsum_celld,'uni',false);
        
        dist_weight(:,:,1) = curr_dist;
        dist_weight(:,:,2) = pyr_weights';
        
        field_lscorr = cellfun(@(x) x(2),curr_radius_corr);
        
        dist_weight_all{curr_animal,curr_day} = dist_weight;
        field_lscorr_all{curr_animal,curr_day} = field_lscorr;%/size(curr_pyr_act,1);
        
        
     
% whoops: didn't need this
%         [xx yy] = meshgrid(1:512);
%         dist_mesh = sqrt((xx - 256).^2+(yy-256).^2);
%         radius_bin_edges = 1:50:501;
%         dist_mean = nan(length(radius_bin_edges)-1,1);
%         for i = 1:length(radius_bin_edges)-1
%             curr_radius = dist_mesh >= radius_bin_edges(i) & ...
%                 dist_mesh < radius_bin_edges(i+1);
%         end     
        
        
%         %for checking cells within radii from gad cells by eye
%         curr_gad = 3;
%         curr_x = roi_center_x(curr_gad);
%         curr_y = roi_center_y(curr_gad);
%         figure; hold on
%         patch(polygon.ROI{curr_gad}(:,1),-polygon.ROI{curr_gad}(:,2), ...
%                    [1 0 0]);
%         for i = pyr_unfilled(pyr_unfilled_move)
%              patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2), ...
%                    [0 1 0]);
%         end
%         for i = setdiff(1:length(roi_labels),[curr_gad pyr_unfilled(pyr_unfilled_move)]);
%              patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1 1 1]);
%         end
%         for i = radius_bin_edges
%             t = 0:0.01:2*pi;
%             x = cos(t)*i + curr_x;
%             y = sin(t)*i - curr_y;
%             plot(x,y,'k');
%         end
%         xlim([0 512])
%         ylim([-512 0])

        curr_day
    end
    

    curr_animal
end

% correlation plot
a = cellfun(@(x) nanmean(x,2),radius_corr_all,'uni',false);
a2 = horzcat(a{:})';
figure;errorbar(nanmedian(a2),nanstd(a2)./sum(~isnan(a2)),'k','linewidth',2);
set(gca,'XTick',1:length(radius_bin_edges)-1)
set(gca,'XTickLabel',radius_bin_edges(2:end))
xlabel('Distance')
ylabel('Correlation')

% number cell plot
z = cellfun(@(x) sum(x,2),num_pyrcells_all,'uni',false);
a = horzcat(z{:})';
a(a == 0) = nan;
figure;errorbar(nanmean(a),nanstd(a)./sum(~isnan(a)),'k','linewidth',2);
set(gca,'XTick',1:length(radius_bin_edges)-1)
set(gca,'XTickLabel',radius_bin_edges(2:end))
xlabel('Distance')
ylabel('Avg num cells')

% correlation plot (by set number of cells)
a = cellfun(@(x) x(:),radius_corr_all,'uni',false);
b = cellfun(@(x) x(:),mean_radius_all,'uni',false);
a1 = vertcat(a{:});
b1 = vertcat(b{:});

radius_bin_edges = [1:50:600];
[n,dist_bin] = histc(b1,radius_bin_edges);
dist_binMean = grpstats(a1,dist_bin);
dist_binSem = grpstats(a1,dist_bin,'sem');
use_bins = unique(dist_bin);
x = radius_bin_edges(use_bins(use_bins ~= 0)+1);
figure;errorbar(x,dist_binMean(use_bins ~= 0), ...
   dist_binSem(use_bins ~= 0),'k','linewidth',2);


% LS weight by distance
use_cells = cellfun(@(x) ~isempty(x),dist_weight_all);
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = {[1:14]};
col = summer(length(day_combine));
figure; hold on;
for curr_day = 1:length(day_combine)
    days = day_combine{curr_day};
    curr_use_cells = false(size(use_cells));
    curr_use_cells(:,days) = use_cells(:,days);
    
    cat_dist_weight = cell(size(curr_use_cells));
    cat_dist_weight(curr_use_cells) = cellfun(@(x) reshape(x,numel(x(:,:,1)),2), ...
        dist_weight_all(curr_use_cells),'uni',false);
    all_dist_weight = vertcat(cat_dist_weight{:});
    radius_bin_edges = [1:50:600];
    [n,dist_weight_bin] = histc(all_dist_weight(:,1),radius_bin_edges);
    dist_weight_binMean = grpstats(all_dist_weight(:,2),dist_weight_bin);
    dist_weight_binSem = grpstats(all_dist_weight(:,2),dist_weight_bin,'sem');
    
    use_bins = unique(dist_weight_bin);
    x = radius_bin_edges(use_bins(use_bins ~= 0)+1);
    %subplot(4,4,curr_day);
    errorbar(x,dist_weight_binMean(use_bins ~= 0), ...
        dist_weight_binSem(use_bins ~= 0),'k','linewidth',2,'color',col(curr_day,:));

end

a = cellfun(@nanmean,field_lscorr_all);
%day_combine = {[1] [2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
lscorr_mean = nan(size(day_combine));
lscorr_sem = nan(size(day_combine));
for curr_day = 1:length(day_combine)
    days = day_combine{curr_day};
    curr_data = reshape(a(:,days),1,[]);
    lscorr_mean(curr_day) = nanmean(curr_data);
    lscorr_sem(curr_day) = nanstd(curr_data)/sum(~isnan(curr_data));
end
figure;errorbar(lscorr_mean,lscorr_sem,'k','linewidth',2)




%% Hierarchy - distance by activity

% get distance matricies
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
gad_dist_matrix = cell(length(animal),1);
pyr_dist_matrix = cell(length(animal),1);
gadpyr_dist_matrix = cell(length(animal),1);
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
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    gad_dist_matrix{curr_animal} = roi_dist(gad_unfilled,gad_unfilled);
    pyr_dist_matrix{curr_animal} = roi_dist(pyr_unfilled,pyr_unfilled);
    gadpyr_dist_matrix{curr_animal} = roi_dist(gad_unfilled,pyr_unfilled);
    
end

use_cells = cellfun(@(x) ~isempty(x),movement_trace_all);
gad_dist_all = cell(size(use_cells));
temp_dist = repmat(gad_dist_matrix,[1 15]);
gad_dist_all(use_cells) = temp_dist(use_cells);

pyr_dist_all = cell(size(use_cells));
temp_dist = repmat(pyr_dist_matrix,[1 15]);
pyr_dist_all(use_cells) = temp_dist(use_cells);

gadpyr_dist_all = cell(size(use_cells));
temp_dist = repmat(gadpyr_dist_matrix,[1 15]);
gadpyr_dist_all(use_cells) = temp_dist(use_cells);

% get subsets of distances for movement cells only
gad_move_dist_all = cell(size(gad_dist_all));
gad_move_dist_all(use_cells) = cellfun(@(x,y) x(y,y), ...
    gad_dist_all(use_cells),gad_class_all(use_cells),'uni',false);

pyr_move_dist_all = cell(size(pyr_dist_all));
pyr_move_dist_all(use_cells) = cellfun(@(x,y) x(y,y), ...
    pyr_dist_all(use_cells),pyr_class_all(use_cells),'uni',false);

gadpyr_move_dist_all = cell(size(pyr_dist_all));
gadpyr_move_dist_all(use_cells) = cellfun(@(x,y,z) x(y,z), ...
    gadpyr_dist_all(use_cells),gad_class_all(use_cells), ...
    pyr_class_all(use_cells),'uni',false);

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get activity during movements

% this is to use all cells
% gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,'uni',false);
% pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,'uni',false);

% this is to use only the movement-classified ones
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);

%%% GAD
% get the geometric mean of activity pairs, normalize from min to max
gad_move_active_meanpair = cellfun(@(x) sqrt(...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))))* ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))))'), ...
    gad_move_active,'uni',false);

gad_move_active_meanpair_itril = ...
    cellfun(@(x) AP_itril(x,-1),gad_move_active_meanpair,'uni',false);
gad_move_dist_all_itril = cellfun(@(x) AP_itril(x,-1),gad_move_dist_all,'uni',false);

gad_cat_meanpair = vertcat(gad_move_active_meanpair_itril{:});
gad_cat_dist = vertcat(gad_move_dist_all_itril{:});

dist_bin_edges = [0:0.1:0.9 Inf];
[n,dist_bins] = histc(gad_cat_meanpair(:,1),dist_bin_edges);
dist_binMean = grpstats(gad_cat_dist,dist_bins);
dist_binSem = grpstats(gad_cat_dist,dist_bins,'sem');
figure;errorbar(dist_binMean,dist_binSem,'k','linewidth',2);
line(xlim,[mean(gad_cat_dist) mean(gad_cat_dist)],'color','r', ...
    'linestyle','--','linewidth',2)
set(gca,'XTick',1:length(dist_bin_edges)-1)
set(gca,'XTickLabel',dist_bin_edges(2:end))
ylabel('Distance')
xlabel('Geometric mean of normalized activity (max range for active during movements)')
title('Gad-gad distance by activity')

%%% PYR
% get the geometric mean of activity pairs, normalize from min to max
pyr_move_active_meanpair = cellfun(@(x) sqrt(...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))))* ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))))'), ...
    pyr_move_active,'uni',false);

pyr_move_active_meanpair_itril = ...
    cellfun(@(x) AP_itril(x,-1),pyr_move_active_meanpair,'uni',false);
pyr_move_dist_all_itril = cellfun(@(x) AP_itril(x,-1),pyr_move_dist_all,'uni',false);

pyr_cat_meanpair = vertcat(pyr_move_active_meanpair_itril{:});
pyr_cat_dist = vertcat(pyr_move_dist_all_itril{:});

dist_bin_edges = [0:0.1:0.9 Inf];
[n,dist_bins] = histc(pyr_cat_meanpair(:,1),dist_bin_edges);
dist_binMean = grpstats(pyr_cat_dist,dist_bins);
dist_binSem = grpstats(pyr_cat_dist,dist_bins,'sem');
figure;errorbar(dist_binMean,dist_binSem,'k','linewidth',2);
line(xlim,[mean(pyr_cat_dist) mean(pyr_cat_dist)],'color','r', ...
    'linestyle','--','linewidth',2)
set(gca,'XTick',1:length(dist_bin_edges)-1)
set(gca,'XTickLabel',dist_bin_edges(2:end))
ylabel('Distance')
xlabel('Geometric mean of normalized activity (max range for active during movements)')
title('Pyr-pyr distance by activity')

%%% GAD-PYR
% get the geometric mean of activity pairs, normalize from min to max
gad_move_active_meannorm = cellfun(@(x) ...
    (nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))), ...
    gad_move_active,'uni',false);
pyr_move_active_meannorm = cellfun(@(x) ...
    (nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2))), ...
    pyr_move_active,'uni',false);

gadpyr_move_active_meanpair = cellfun(@(x,y) x*y', ...
    gad_move_active_meannorm, pyr_move_active_meannorm,'uni',false);

gadpyr_move_active_meanpair_itril = ...
    cellfun(@(x) x(:),gadpyr_move_active_meanpair,'uni',false);
gadpyr_move_dist_all_itril = cellfun(@(x) x(:),gadpyr_move_dist_all,'uni',false);

gadpyr_cat_meanpair = vertcat(gadpyr_move_active_meanpair_itril{:});
gadpyr_cat_dist = vertcat(gadpyr_move_dist_all_itril{:});

dist_bin_edges = [0:0.1:0.9 Inf];
[n,dist_bins] = histc(gadpyr_cat_meanpair(:,1),dist_bin_edges);
dist_binMean = grpstats(gadpyr_cat_dist,dist_bins);
dist_binSem = grpstats(gadpyr_cat_dist,dist_bins,'sem');
figure;errorbar(dist_binMean,dist_binSem,'k','linewidth',2);
line(xlim,[mean(gadpyr_cat_dist) mean(gadpyr_cat_dist)],'color','r', ...
    'linestyle','--','linewidth',2)
set(gca,'XTick',1:length(dist_bin_edges)-1)
set(gca,'XTickLabel',dist_bin_edges(2:end))
ylabel('Distance')
xlabel('Geometric mean of normalized activity (max range for active during movements)')
title('Gad-pyr distance by activity')


%% Hierarchy - plot cells color coded by activity

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% get normalized mean activity during movement
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
gad_move_active_meanpair = cellfun(@(x) ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2)))), ...
    gad_move_active,'uni',false);

pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);
pyr_move_active_meanpair = cellfun(@(x) ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2)))), ...
    pyr_move_active,'uni',false);

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
gad_dist_matrix = cell(length(animal),1);
pyr_dist_matrix = cell(length(animal),1);
gadpyr_dist_matrix = cell(length(animal),1);
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % load movement classification
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);

    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    sessions = find(cellfun(@(x) ~isempty(x),movement_trace_all(curr_animal,:)));
    
    figure;
    set(gcf,'Name',animal);
    for curr_day = sessions
        
        gad_unfilled_move = classified_cells.move_cells_peak{curr_day}(gad_unfilled);
        pyr_unfilled_move = classified_cells.move_cells_peak{curr_day}(pyr_unfilled);
        
        subplot(4,4,curr_day);
        
        circle_radius = 20;
        circle_sub = sqrt(2*circle_radius^2);
        
        % Show space map of gad colored move / gad unfilled
        %use_cells = gad_unfilled(gad_unfilled_move);
        use_cells = gad_unfilled(gad_unfilled_move);
        for i = 1:length(use_cells)
            curr_cell = use_cells(i);
            curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
            curr_x = round(roi_center_x(curr_cell));
            curr_y = -round(roi_center_y(curr_cell));
            
            circle_corner = [curr_x-circle_sub curr_y+circle_sub];
            rectangle('Position',[circle_corner circle_radius*2 circle_radius*2], ...
                'Curvature',[1 1],'FaceColor',[0 curr_reliability 0])
        end
        use_cells = gad_unfilled(~gad_unfilled_move);
        for i = 1:length(use_cells)
            curr_cell = use_cells(i);
            curr_x = round(roi_center_x(curr_cell));
            curr_y = -round(roi_center_y(curr_cell));
            
            circle_corner = [curr_x-circle_sub curr_y+circle_sub];
            rectangle('Position',[circle_corner circle_radius*2 circle_radius*2], ...
                'Curvature',[1 1])
        end
        
        xlim([1 512]);
        ylim([-512 -1]);
        axis off
        
         % Didn't need this approach: only useful for blurring
        
%         % Show space map of pyramidal move cells
%         curr_img = subplot(4,4,curr_day);
%         spread_filter = fspecial('disk',20);
%         spread_filter(spread_filter > 0) = 1;
%         
%         use_cells = gad_unfilled(gad_unfilled_move);
%         gadmove_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
%             curr_x = round(roi_center_x(curr_cell));
%             curr_y = round(roi_center_y(curr_cell));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             gadmove_activity_map(curr_y,curr_x) = curr_reliability;
%         end
%         gadmove_activity_map_conv = conv2(gadmove_activity_map,spread_filter,'same');
%         
%         use_cells = gad_unfilled(~gad_unfilled_move);
%         gadnonmove_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_x = round(roi_center_x(curr_cell));
%             curr_y = round(roi_center_y(curr_cell));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             gadnonmove_activity_map(curr_y,curr_x) = 1;
%         end
%         gadnonmove_activity_map_conv = conv2(gadnonmove_activity_map,spread_filter,'same');
%         gadnonmove_activity_map_conv(gadnonmove_activity_map_conv > 0) = 0.2;
%         
%         gad_map = gadmove_activity_map_conv;
%         gad_map(:,:,2) = zeros(size(gad_map));
%         gad_map(:,:,3) = gadnonmove_activity_map_conv;
%         
%         % make background white
%         bg = repmat(~any(gad_map,3),[1 1 3]);
%         gad_map(bg) = 1;      
%         use_cells = setdiff(1:length(roi_labels),gad_unfilled);
%         pyr_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_x = round(roi_center_x(curr_cell));
%             curr_y = round(roi_center_y(curr_cell));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             pyr_activity_map(curr_y,curr_x) = 1;
%         end
%         pyr_activity_map_conv = conv2(pyr_activity_map,spread_filter,'same');
        

        
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
%             patch(polygon.ROI{curr_cell}(:,1),-polygon.ROI{curr_cell}(:,2), ...
%                 [curr_reliability 0 0]);
%         end
%         for i = gad_unfilled(~gad_unfilled_move);
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0 0 1]);
%         end
%         for i = setdiff(1:length(roi_labels),gad_unfilled)
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1 1 1]);
%         end
%         xlim([1 512]);
%         ylim([-512 -1]);
        
    end
end

%% Hierarchy - spatial entropy

% get normalized mean activity during movement
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
gad_move_active_meanpair = cellfun(@(x) ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2)))), ...
    gad_move_active,'uni',false);

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
jitter_dist = 50:50:500;
field_energy_p = nan(8,15,length(jitter_dist));
field_entropy_p = nan(8,15,length(jitter_dist));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % load movement classification
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);

    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    sessions = find(cellfun(@(x) ~isempty(x),movement_trace_all(curr_animal,:)));
    

    for curr_day = sessions
        
        gad_unfilled_move = classified_cells.move_cells_peak{curr_day}(gad_unfilled);
        pyr_unfilled_move = classified_cells.move_cells_peak{curr_day}(pyr_unfilled);
        
        % Energy from Hira et al 2013
        % get the energy from the paper
        use_cells = gad_unfilled(gad_unfilled_move);
        curr_roi_dist = roi_dist(use_cells,use_cells);
        roi_dist_unique = AP_itril(curr_roi_dist,-1);
        Jnorm = sum(1./roi_dist_unique);       
        FiFj_all = gad_move_active_meanpair{curr_animal,curr_day} * ...
            gad_move_active_meanpair{curr_animal,curr_day}';
        FiFj = AP_itril(FiFj_all,-1);
        Jij = (1./roi_dist_unique)/Jnorm;
        Eg = -FiFj'*Jij;
        Eg_hat = -mean(FiFj);        
        delta_Eg = Eg - Eg_hat;
        
        % shuffle distances to get chance
        use_cells = gad_unfilled(gad_unfilled_move);
        use_cells_centers_sub = [round(roi_center_y(use_cells)) ...
            round(roi_center_x(use_cells))];
        use_cells_centers_sub(use_cells_centers_sub < 1) = 1;
        use_cells_centers_sub(use_cells_centers_sub > 512) = 512;
        
        delta_Eg_p = nan(length(jitter_dist),1);
        for curr_jitter = 1:length(jitter_dist)
            
            num_shuff = 10000;
            
            [jitter_x jitter_y] = pol2cart(2*pi*rand( ...
                length(use_cells),num_shuff),jitter_dist(curr_jitter));
            use_cells_centers_sub_jitter_y = ...
                repmat(use_cells_centers_sub(:,2),1,num_shuff)+jitter_y;
            use_cells_centers_sub_jitter_x = ...
                repmat(use_cells_centers_sub(:,1),1,num_shuff)+jitter_x;
            roi_dist_jitter = cell2mat(arrayfun(@(x) pdist( ...
                [use_cells_centers_sub_jitter_y(:,x) ...
                use_cells_centers_sub_jitter_x(:,x)])',1:num_shuff,'uni',false))';
                      
            %curr_perm = meshgrid(1:length(roi_dist_unique),1:num_shuff);
            %curr_perm = shake(curr_perm,2);
            %roi_dist_shuff = roi_dist_unique(curr_perm);
            
            Jnorm = sum(1./roi_dist_jitter,2);
            FiFj_all = gad_move_active_meanpair{curr_animal,curr_day} * ...
                gad_move_active_meanpair{curr_animal,curr_day}';
            FiFj = AP_itril(FiFj_all,-1);
            Jij = bsxfun(@times,1./roi_dist_jitter,1./Jnorm);
            Eg = -FiFj'*Jij';
            Eg_hat = -mean(FiFj);
            delta_Eg_shuff = Eg - Eg_hat;
            delta_Eg_rank = tiedrank([delta_Eg_shuff delta_Eg]);
            delta_Eg_p(curr_jitter) = delta_Eg_rank(end)/num_shuff;
        end
        field_energy_p(curr_animal,curr_day,:) = delta_Eg_p;

        
        
        
        
        
%         % Entropy from Vergasolla et al. 2007
%         spread_filter = fspecial('gaussian',800,50);
%         use_cells = gad_unfilled(gad_unfilled_move);
%         use_cells_centers_sub = [round(roi_center_y(use_cells)) ...
%             round(roi_center_x(use_cells))];
%         use_cells_centers_sub(use_cells_centers_sub < 1) = 1;
%         use_cells_centers_sub(use_cells_centers_sub > 512) = 512;
%         use_cells_centers = sub2ind([512 512],use_cells_centers_sub(:,1), ...
%             use_cells_centers_sub(:,2));
%         gadmove_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
%             gadmove_activity_map(use_cells_centers(i)) = curr_reliability;
%         end
%         gadmove_activity_map_conv = conv2(gadmove_activity_map,spread_filter,'same');
%         gadmove_activity_map_conv = gadmove_activity_map_conv/ ...
%                 sum(gadmove_activity_map_conv(:));
%         entropy_map = gadmove_activity_map_conv*log(gadmove_activity_map_conv);
%         sum_entropy = -sum(entropy_map(:));
%         % shuffle centers to get chance
%         sum_entropy_p = nan(length(jitter_dist),1);
%         for curr_jitter = 1:length(jitter_dist)
%             num_shuff = 1000;
%             sum_entropy_shuff = nan(num_shuff,1);
%             for curr_shuff = 1:num_shuff
%                 %use_cells_centers_shuff =  ...
%                 %    use_cells_centers(randperm(length(use_cells_centers)));
%                 % jitter centers by distance
%                 [jitter_x jitter_y] = pol2cart(2*pi*rand( ...
%                     length(use_cells_centers),1),jitter_dist(curr_jitter));                
%                 use_cells_centers_sub_jitter = ...
%                     [use_cells_centers_sub(:,1)+round(jitter_x), ...
%                     use_cells_centers_sub(:,2)+round(jitter_y)];
%                 use_cells_centers_sub_jitter(use_cells_centers_sub_jitter < 1) = 1;
%                 use_cells_centers_sub_jitter(use_cells_centers_sub_jitter > 512) = 512;
%                 use_cells_centers_jitter = sub2ind([512 512],use_cells_centers_sub_jitter(:,1), ...
%                     use_cells_centers_sub_jitter(:,2));
%                 
%                 gadmove_activity_map = zeros(512,512);
%                 for i = 1:length(use_cells)
%                     curr_cell = use_cells(i);
%                     curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
%                     gadmove_activity_map(use_cells_centers_jitter(i)) = curr_reliability;
%                 end
%                 gadmove_activity_map_conv = conv2(gadmove_activity_map,spread_filter,'same');
%                 gadmove_activity_map_conv = gadmove_activity_map_conv/ ...
%                     sum(gadmove_activity_map_conv(:));
%                 entropy_map = gadmove_activity_map_conv.*log(gadmove_activity_map_conv);
%                 sum_entropy_shuff(curr_shuff) = -sum(entropy_map(:));
%             end
%             sum_entropy_rank = tiedrank([sum_entropy_shuff;sum_entropy]);
%             sum_entropy_p(curr_jitter) = sum_entropy_rank(end)/num_shuff;
%         end
%         field_entropy_p(curr_animal,curr_day,:) = sum_entropy_p;
        
        curr_day
%         
%         use_cells = gad_unfilled(~gad_unfilled_move);
%         gadnonmove_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_x = round(roi_center_x(curr_cell));
%             curr_y = round(roi_center_y(curr_cell));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             gadnonmove_activity_map(curr_y,curr_x) = 1;
%         end
%         gadnonmove_activity_map_conv = conv2(gadnonmove_activity_map,spread_filter,'same');
%         gadnonmove_activity_map_conv(gadnonmove_activity_map_conv > 0) = 0.2;
%         
%         gad_map = gadmove_activity_map_conv;
%         gad_map(:,:,2) = zeros(size(gad_map));
%         gad_map(:,:,3) = gadnonmove_activity_map_conv;
%         
%         % make background white
%         bg = repmat(~any(gad_map,3),[1 1 3]);
%         gad_map(bg) = 1;      
%         use_cells = setdiff(1:length(roi_labels),gad_unfilled);
%         pyr_activity_map = zeros(512,512);
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_x = round(roi_center_x(curr_cell));
%             curr_y = round(roi_center_y(curr_cell));
%             if curr_x < 1;curr_x = 1;end
%             if curr_y < 1;curr_y = 1;end
%             if curr_x > 512;curr_x = 512;end
%             if curr_y > 512;curr_y = 512;end
%             pyr_activity_map(curr_y,curr_x) = 1;
%         end
%         pyr_activity_map_conv = conv2(pyr_activity_map,spread_filter,'same');
        

        
%         for i = 1:length(use_cells)
%             curr_cell = use_cells(i);
%             curr_reliability = gad_move_active_meanpair{curr_animal,curr_day}(i);
%             patch(polygon.ROI{curr_cell}(:,1),-polygon.ROI{curr_cell}(:,2), ...
%                 [curr_reliability 0 0]);
%         end
%         for i = gad_unfilled(~gad_unfilled_move);
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0 0 1]);
%         end
%         for i = setdiff(1:length(roi_labels),gad_unfilled)
%             patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1 1 1]);
%         end
%         xlim([1 512]);
%         ylim([-512 -1]);
        
    end

    curr_animal

end







%% TEMPORARY - look at spatial distributions of cells with each movement
activity_map_corr_p = nan(8,15);
activity_map_meancorr_p = nan(8,15);


% get normalized mean activity during movement
gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);

gad_move_active_df = cellfun(@(x,y,w) cell2mat(cellfun(@(z) max(x(w,z),[],2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active_df = cellfun(@(x,y,w) cell2mat(cellfun(@(z) max(x(w,z),[],2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % load movement classification
    classification_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    load([classification_path filesep animal '_classified_cells_caEvents_5.mat']);

    
    % load ROI labels and identify cell type
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
    roilabel_name = [animal '_roilabels.roilabel'];
    load([analysis_path filesep roilabel_name],'-MAT')
    cells = 1:length(roi_labels);
    gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels)';
    pyr_cells = ~gad_cells;
    filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels)';
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    sessions = find(cellfun(@(x) ~isempty(x),movement_trace_all(curr_animal,:)));
    

    for curr_day = sessions
        
        gad_unfilled_move = classified_cells.move_cells_peak{curr_day}(gad_unfilled);
        pyr_unfilled_move = classified_cells.move_cells_peak{curr_day}(pyr_unfilled);
        
        imscale = 10;
        imedge = floor(512/imscale);
        spread_filter = fspecial('gaussian',800/imscale,100/imscale);
        
        gad_use_cells = gad_unfilled(gad_unfilled_move);
        gad_use_cells_centers_sub = round([roi_center_y(gad_use_cells) ...
            roi_center_x(gad_use_cells)]/imscale);
        gad_use_cells_centers_sub(gad_use_cells_centers_sub < 1) = 1;
        gad_use_cells_centers_sub(gad_use_cells_centers_sub > imedge) = imedge;
        gad_use_cells_centers = sub2ind([imedge imedge],gad_use_cells_centers_sub(:,1), ...
            gad_use_cells_centers_sub(:,2));
        
        pyr_use_cells = pyr_unfilled(pyr_unfilled_move);
        pyr_use_cells_centers_sub = round([roi_center_y(pyr_use_cells) ...
            roi_center_x(pyr_use_cells)]/imscale);
        pyr_use_cells_centers_sub(pyr_use_cells_centers_sub < 1) = 1;
        pyr_use_cells_centers_sub(pyr_use_cells_centers_sub > imedge) = imedge;
        pyr_use_cells_centers = sub2ind([imedge imedge],pyr_use_cells_centers_sub(:,1), ...
            pyr_use_cells_centers_sub(:,2));
        
        gad_move_matrix = gad_move_active_df{curr_animal,curr_day};
        pyr_move_matrix = pyr_move_active_df{curr_animal,curr_day};

        [asdf activity_sort] = sort(nanmean(gad_move_matrix > 0));
        
        %gadmove_activity_map = zeros(512,512,size(move_matrix,2));
        activity_map = zeros(imedge,imedge,3,size(gad_move_matrix,2));
        gad_activity_map = zeros(imedge^2,size(gad_move_matrix,2));
        pyr_activity_map = zeros(imedge^2,size(gad_move_matrix,2));
        for curr_move = 1:size(gad_move_matrix,2)
            gad_map = zeros(imedge,imedge);
            for i = 1:length(gad_use_cells)                
                curr_cell = gad_use_cells(i);
                curr_df = gad_move_matrix(i,curr_move);
                gad_map(gad_use_cells_centers(i)) = curr_df;
            end
            gad_conv = conv2(gad_map,spread_filter,'same');
            gad_activity_map(:,curr_move) = gad_conv(:);
            activity_map(:,:,1,activity_sort==curr_move) = gad_conv/max(gad_conv(:));
            
            pyr_map = zeros(imedge,imedge);
            for i = 1:length(pyr_use_cells)                
                curr_cell = pyr_use_cells(i);
                curr_df = pyr_move_matrix(i,curr_move);
                pyr_map(pyr_use_cells_centers(i)) = curr_df;
            end
            pyr_conv = conv2(pyr_map,spread_filter,'same');
            pyr_activity_map(:,curr_move) = pyr_conv(:);
            activity_map(:,:,2,activity_sort==curr_move) = pyr_conv/max(pyr_conv(:));
        end
        
        temp_corr = corrcoef(gad_activity_map,pyr_activity_map);
        activity_map_corr = temp_corr(2);
        
        temp_corr = arrayfun(@(x) corrcoef( ...
            gad_activity_map(:,x),pyr_activity_map(:,x)), ...
            1:size(gad_activity_map,2),'uni',false);
        activity_map_meancorr = nanmean(cellfun(@(x) x(2), temp_corr));
        
        num_shuff = 100;
        activity_map_corr_shuff = nan(num_shuff,1);       
        activity_map_meancorr_shuff = nan(num_shuff,1);
        for curr_shuff = 1:num_shuff
            gad_use_cells_centers_shuff = ...
                gad_use_cells_centers(randperm(length(gad_use_cells_centers)));
            pyr_use_cells_centers_shuff = ...
                pyr_use_cells_centers(randperm(length(pyr_use_cells_centers)));
            
            gad_activity_map = zeros(imedge^2,size(gad_move_matrix,2));
            pyr_activity_map = zeros(imedge^2,size(gad_move_matrix,2));
            for curr_move = 1:length(activity_sort);
                gad_map = zeros(imedge,imedge);
                for i = 1:length(gad_use_cells)
                    curr_cell = gad_use_cells(i);
                    curr_df = gad_move_matrix(i,curr_move);
                    gad_map(gad_use_cells_centers_shuff(i)) = curr_df;
                end
                gad_conv = conv2(gad_map,spread_filter,'same');
                gad_activity_map(:,curr_move) = gad_conv(:);
                
                pyr_map = zeros(imedge,imedge);
                for i = 1:length(pyr_use_cells)
                    curr_cell = pyr_use_cells(i);
                    curr_df = pyr_move_matrix(i,curr_move);
                    pyr_map(pyr_use_cells_centers_shuff(i)) = curr_df;
                end
                pyr_conv = conv2(pyr_map,spread_filter,'same');
                pyr_activity_map(:,curr_move) = pyr_conv(:);
            end          
            disp(curr_shuff)
            temp_corr = corrcoef(gad_activity_map,pyr_activity_map);
            activity_map_corr_shuff(curr_shuff) = temp_corr(2);
        
            temp_corr = arrayfun(@(x) corrcoef( ...
                gad_activity_map(:,x),pyr_activity_map(:,x)), ...
                1:size(gad_activity_map,2),'uni',false);
            activity_map_meancorr_shuff(curr_shuff) = nanmean(cellfun(@(x) x(2), temp_corr));
                
        end
        activity_map_corr_rank = tiedrank([activity_map_corr_shuff;...
            activity_map_corr]);
        activity_map_corr_p(curr_animal,curr_day) = ...
            activity_map_corr_rank(end)/num_shuff;
        
        activity_map_meancorr_rank = tiedrank([activity_map_meancorr_shuff;...
            activity_map_meancorr]);
        activity_map_meancorr_p(curr_animal,curr_day) = ...
            activity_map_meancorr_rank(end)/num_shuff;
        
        curr_day
    end

    curr_animal
end



%% Gad-Pyr relatedness weights by distance


animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
gad_dist_matrix = cell(length(animal),1);
pyr_dist_matrix = cell(length(animal),1);
gadpyr_dist_matrix = cell(length(animal),1);
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
    pyr_unfilled = find(pyr_cells & ~filled_cells);
    gad_unfilled = find(gad_cells & ~filled_cells);
    
    gad_cells = find(gad_cells);
    pyr_cells = find(pyr_cells);
    
    % load ROIs and get distances
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
    
    gad_dist_matrix{curr_animal} = roi_dist(gad_unfilled,gad_unfilled);
    pyr_dist_matrix{curr_animal} = roi_dist(pyr_unfilled,pyr_unfilled);
    gadpyr_dist_matrix{curr_animal} = roi_dist(gad_unfilled,pyr_unfilled);
    
end

use_cells = cellfun(@(x) ~isempty(x),movement_trace_all);
gad_dist_all = cell(size(use_cells));
temp_dist = repmat(gad_dist_matrix,[1 15]);
gad_dist_all(use_cells) = temp_dist(use_cells);

pyr_dist_all = cell(size(use_cells));
temp_dist = repmat(pyr_dist_matrix,[1 15]);
pyr_dist_all(use_cells) = temp_dist(use_cells);

gadpyr_dist_all = cell(size(use_cells));
temp_dist = repmat(gadpyr_dist_matrix,[1 15]);
gadpyr_dist_all(use_cells) = temp_dist(use_cells);

% get normalized mean activity during movement
% gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
% pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);
gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,'uni',false);

gad_move_active_reliability = cellfun(@(x) ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2)))), ...
    gad_move_active,'uni',false);
pyr_move_active_reliability = cellfun(@(x) ...
    ((nansum(x,2)-min(nansum(x,2)))/(max(nansum(x,2))-min(nansum(x,2)))), ...
    pyr_move_active,'uni',false);

gad_move_active_df = cellfun(@(x,y,w) cell2mat(cellfun(@(z) max(x(w,z),[],2),y,'uni',false)), ...
    gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
pyr_move_active_df = cellfun(@(x,y,w) cell2mat(cellfun(@(z) max(x(w,z),[],2),y,'uni',false)), ...
    pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);


weights = cell(8,1);
for i = 1:8
    curr_pyr_r = horzcat(pyr_move_active_reliability{i,:})';
    curr_gad_r = horzcat(gad_move_active_reliability{i,:})';
    curr_gadpyr_dist = gadpyr_dist_matrix{i};
    
    weights{i} = transpose(curr_pyr_r\curr_gad_r);   
end
a = cellfun(@(x) x(:),gadpyr_dist_matrix,'uni',false);
b = cellfun(@(x) x(:),weights,'uni',false);
a = vertcat(a{:});
b = vertcat(b{:});

radius_bin_edges = [1:50:600];
[n,dist_weight_bin] = histc(a,radius_bin_edges);
dist_weight_binMean = grpstats(b,dist_weight_bin);
dist_weight_binSem = grpstats(b,dist_weight_bin,'sem');

% do this but daily: with movements seperated

animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
gad_dist_matrix = cell(length(animal),1);
pyr_dist_matrix = cell(length(animal),1);
gadpyr_dist_matrix = cell(length(animal),1);






%% Move-onset/reward modulated cells


% animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
% for i = 1:8
%     animal = animals{i};
%     
%     % load ROI labels and identify cell type
%     analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];    
%     roilabel_name = [animal '_roilabels.roilabel'];
%     load([analysis_path filesep roilabel_name],'-MAT')
%     cells = 1:length(roi_labels);
%     gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels);
%     pyr_cells = ~gad_cells;
%     filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels);
%     
%     pyr_unfilled = pyr_cells & ~filled_cells;
%     gad_unfilled = gad_cells & ~filled_cells;
% end
% 
% 
% 
% caEvent_matrix = AP_caEvents(im.roi_trace_df,pyr_unfilled,gad_unfilled);
% movement_trace_all
% move_epoch_frames
% rewarded_movements
% cued_rewarded_movements
% 
% cue_frames_all
% reward_frames_all
% cued_reward_frames_all

reward_cells_all = cell(size(movement_trace_all));
move_onset_cells_all = cell(size(movement_trace_all));

animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';

for curr_animalday = animaldays
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_reward_frames_mult = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_reward_frames = cellfun(@(x) x(1),curr_reward_frames_mult);
    
    curr_move_onset_frames = cellfun(@(x) x(1),curr_rewardmoves);
    
    %curr_reward_frames = reward_frames_all{curr_animalday};
    %curr_move_onset_frames = cellfun(@(x) x(1),move_epoch_frames{curr_animalday});
    
    spread_frames = 1;
    spread_filter = ones(1,spread_frames);
    pyr_spread = conv2(pyr_activity_all{curr_animalday},spread_filter,'same');
    
    num_shift = 1000;
    % this is given by the minimum and maximum values for having seperate
    % distributions given spread (function of spike uncertainty and leeway for
    % being overlapping)
    % I think that reasoning is wrong, but not sure...
    %frame_shift = [-spread_frames:-spread_frames/2 spread_frames/2:spread_frames];
    %frame_shift = [-spread_frames*2:-spread_frames spread_frames:spread_frames*2];
    % This was changed: instead of spreading, just use start to peak
    frame_shift = [-10 10];
    frames_back = 20;
    frames_forward = 60;
    
    reward_align_p_all = nan(size(pyr_spread,1),1);
    reward_align_mean = nan(size(pyr_spread,1),1);
    reward_align_95 = nan(size(pyr_spread,1),1);
    
    move_onset_align_p_all = nan(size(pyr_spread,1),1);
    move_onset_align_mean = nan(size(pyr_spread,1),1);
    move_onset_align_95 = nan(size(pyr_spread,1),1);
    for curr_cell = find(pyr_class_all{curr_animalday});
        
        reward_align = nan(length(curr_reward_frames),frames_back+frames_forward+1);
        for i = 1:length(curr_reward_frames);
            curr_reward = round(curr_reward_frames(i));
            roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
            if any(roi_reward_frames<=0) || any(roi_reward_frames>size(pyr_spread,2));
                continue
            end
            reward_align(i,:) = pyr_spread(curr_cell,roi_reward_frames);
        end
        reward_align_mean(curr_cell) = max(nanmean(reward_align));
        
        [m n] = size(reward_align);
        reward_align_buffer = [reward_align(:,end-(abs(min(frame_shift))-1):end) ...
            reward_align reward_align(:,1:max(frame_shift))];
        reward_align_shift_mean = nan(num_shift,1);
        for curr_shift = 1:num_shift
            randshift = frame_shift(randi(length(frame_shift),m,1));
            reward_align_shift = nan(m,n);
            for i = 1:m
                reward_align_shift(i,:) = ...
                    reward_align_buffer(i,abs(min(frame_shift)) + randshift(i)+1: ...
                    n + max(frame_shift) + randshift(i));
            end
            reward_align_shift_mean(curr_shift) = max(nanmean(reward_align_shift));
        end
        reward_align_rank = tiedrank([reward_align_shift_mean;reward_align_mean(curr_cell)]);
        reward_align_p = reward_align_rank(end)/num_shift;
        
        reward_align_95(curr_cell) = prctile(reward_align_shift_mean,95);
        reward_align_p_all(curr_cell) = reward_align_p;
        
        
        move_onset_align = nan(length(curr_move_onset_frames),frames_back+frames_forward+1);
        for i = 1:length(curr_move_onset_frames);
            curr_move_onset = round(curr_move_onset_frames(i));
            roi_move_onset_frames = curr_move_onset-frames_back:curr_move_onset+frames_forward;
            if any(roi_move_onset_frames<=0) || any(roi_move_onset_frames>size(pyr_spread,2));
                continue
            end
            move_onset_align(i,:) = pyr_spread(curr_cell,roi_move_onset_frames);
        end
        move_onset_align_mean(curr_cell) = max(nanmean(move_onset_align));
        
        [m n] = size(move_onset_align);
        move_onset_align_buffer = [move_onset_align(:,end-(abs(min(frame_shift))-1):end) ...
            move_onset_align move_onset_align(:,1:max(frame_shift))];
        move_onset_align_shift_mean = nan(num_shift,1);
        for curr_shift = 1:num_shift
            randshift = frame_shift(randi(length(frame_shift),m,1));
            move_onset_align_shift = nan(m,n);
            for i = 1:m
                move_onset_align_shift(i,:) = ...
                    move_onset_align_buffer(i,abs(min(frame_shift)) + randshift(i)+1: ...
                    n + max(frame_shift) + randshift(i));
            end
            move_onset_align_shift_mean(curr_shift) = max(nanmean(move_onset_align_shift));
        end
        move_onset_align_rank = tiedrank([move_onset_align_shift_mean;move_onset_align_mean(curr_cell)]);
        move_onset_align_p = move_onset_align_rank(end)/num_shift;
        
        move_onset_align_95(curr_cell) = prctile(move_onset_align_shift_mean,95);
        move_onset_align_p_all(curr_cell) = move_onset_align_p;
        
        %curr_cell
        
    end
  
    reward_cells = reward_align_p_all > 0.95;
    move_onset_cells = move_onset_align_p_all > 0.95;
    
    % resolve dual classified cells with shuffled 95 prctile
    dual_class_cells = reward_cells & move_onset_cells;   
    reward_higher = reward_align_95 > move_onset_align_95;
    move_onset_higher = move_onset_align_95 > reward_align_95;
    
    reward_cells(dual_class_cells & move_onset_higher) = false;
    move_onset_cells(dual_class_cells & reward_higher) = false;
    
    reward_cells_all{curr_animalday} = reward_cells;
    move_onset_cells_all{curr_animalday} = move_onset_cells;
    
    % To plot aligned activity: 
 
%     reward_cells_idx = find(reward_cells);
%     move_onset_cells_idx = find(move_onset_cells);
%     
%     hr = figure;
%     set(hr,'Name','Reward-aligned')
%     n_subs = ceil(sqrt(length(reward_cells_idx)));
%     reward_align_all = cell(size(reward_cells_idx));
%     for curr_cells = 1:length(reward_cells_idx)
%         curr_cell = reward_cells_idx(curr_cells);
%     
%         % plot around reward
%         frames_back = 30;
%         frames_forward = 90;
%         reward_align = nan(length(curr_reward_frames),frames_back+frames_forward+1);
%         for i = 1:length(curr_reward_frames);
%             curr_reward = round(curr_reward_frames(i));
%             roi_reward_frames = curr_reward-frames_back:curr_reward+frames_forward;
%             if any(roi_reward_frames<=0) || any(roi_reward_frames>size(pyr_spread,2));
%                 continue
%             end
%             reward_align(i,:) = pyr_spread(curr_cell,roi_reward_frames);
%         end
%         reward_align_all{curr_cells} = reward_align;
%         figure(hr);
%         subplot(n_subs,n_subs,curr_cells);
%         imagesc(reward_align);
%         colormap(gray);
%     end
%     
%     hc = figure;
%     set(hc,'Name','Move onset-aligned')
%     n_subs = ceil(sqrt(length(move_onset_cells_idx)));
%     move_onset_align_all = cell(size(move_onset_cells_idx));
%     for curr_cells = 1:length(move_onset_cells_idx)
%         curr_cell = move_onset_cells_idx(curr_cells);
%         % plot around move_onset
%         move_onset_align = nan(length(curr_move_onset_frames),frames_back+frames_forward+1);
%         for i = 1:length(curr_move_onset_frames);
%             curr_move_onset = round(curr_move_onset_frames(i));
%             roi_move_onset_frames = curr_move_onset-frames_back:curr_move_onset+frames_forward;
%             if any(roi_move_onset_frames<=0) || any(roi_move_onset_frames>size(pyr_spread,2));
%                 continue
%             end
%             move_onset_align(i,:) = pyr_spread(curr_cell,roi_move_onset_frames);
%         end
%         move_onset_align_all{curr_cells} = move_onset_align;
%         figure(hc);
%         subplot(n_subs,n_subs,curr_cells);
%         imagesc(move_onset_align);
%         colormap(gray);
%     end
%     
    disp(curr_animalday / max(animaldays))
    
end

%% Plots from move-onest/reward cells

movement_cells = cellfun(@nanmean,pyr_class_all);
movement_cells = movement_cells(:,1:14);
bhv_cells = cellfun(@(x,y) nanmean(x|y),move_onset_cells_all,reward_cells_all);
bhv_cells = bhv_cells(:,1:14);

move_onset_cells = cellfun(@nanmean,move_onset_cells_all);
reward_cells = cellfun(@nanmean,reward_cells_all);

bhv_movement_frac = bhv_cells./movement_cells;

day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = bhv_movement_frac(:,day_combine{i});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
figure;errorbar(data_mean,data_sem,'k','linewidth',2);
title('Fraction of movement-cells which are bhv-correlated');


day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = move_onset_cells(:,day_combine{i});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
figure; hold on
errorbar(data_mean,data_sem,'r','linewidth',2);
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = reward_cells(:,day_combine{i});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
errorbar(data_mean,data_sem,'b','linewidth',2);
legend({'Move-onset' 'Reward'});
title('Fraction of bhv-correlated cells');


%% Find reward-related activity compared to catch trials
rewardcell_mean_reward_activity = nan(size(catch_rewarded_movements));
rewardcell_mean_catch_activity = nan(size(catch_rewarded_movements));

moveonsetcell_mean_reward_activity = nan(size(catch_rewarded_movements));
moveonsetcell_mean_catch_activity = nan(size(catch_rewarded_movements));

animaldays_catch = find(cellfun(@sum,catch_rewarded_movements) > 0)';
for curr_animalday = animaldays_catch
    % get mean df/f for reward cells reward v catch trials
    curr_rewardmoves_real = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday} & ...
        ~catch_rewarded_movements{curr_animalday});
    
    curr_rewardmoves_catch = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday} & ...
        catch_rewarded_movements{curr_animalday});
        
    curr_pyrmoveact_real = pyr_activity_all{curr_animalday} ...
        (reward_cells_all{curr_animalday},horzcat(curr_rewardmoves_real{:}));
    
    curr_pyrmoveact_catch = pyr_activity_all{curr_animalday} ...
        (reward_cells_all{curr_animalday},horzcat(curr_rewardmoves_catch{:}));
    
    rewardcell_mean_reward_activity(curr_animalday) = nanmean(curr_pyrmoveact_real(:));
    rewardcell_mean_catch_activity(curr_animalday) = nanmean(curr_pyrmoveact_catch(:));
    
    
    % get mean df/f for move onset cells reward v catch trials
    curr_rewardmoves_real = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday} & ...
        ~catch_rewarded_movements{curr_animalday});
    
    curr_rewardmoves_catch = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday} & ...
        catch_rewarded_movements{curr_animalday});    
    
    curr_pyrmoveact_real = pyr_activity_all{curr_animalday} ...
        (move_onset_cells_all{curr_animalday},horzcat(curr_rewardmoves_real{:}));
    
    curr_pyrmoveact_catch = pyr_activity_all{curr_animalday} ...
        (move_onset_cells_all{curr_animalday},horzcat(curr_rewardmoves_catch{:}));
    
    moveonsetcell_mean_reward_activity(curr_animalday) = nanmean(curr_pyrmoveact_real(:));
    moveonsetcell_mean_catch_activity(curr_animalday) = nanmean(curr_pyrmoveact_catch(:));

end



%% Get time preferences for cells based on rewarded movements
% movement_trace_all
% move_epoch_frames
% rewarded_movements
% cued_rewarded_movements
% 
% cue_frames_all
% reward_frames_all
% cued_reward_frames_all
%%%%%%%% FIX WITH LINSPACE (contained in Similarity by correlation/distance
%%%%%%%% section)
movement_time_all = cell(size(movement_trace_all));
movement_rescaled_time_all = cell(size(movement_trace_all));
mean_activity_all = cell(size(movement_trace_all));
movement_time_std_all = cell(size(movement_trace_all));

pyrmoveact_timerescale_full_all = cell(size(movement_trace_all));
pyrmoveact_timerescale_all = cell(size(movement_trace_all));
pyrmoveact_timerescale_mean_all = cell(size(movement_trace_all));
animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';

for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}) - x(1) + 1,curr_rewardmoves,'uni',false); 
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    curr_timeframes = cellfun(@(x) linspace(0,1,length(x)), ...
        curr_rewardmoves,'uni',false);
    curr_rescaled_timeframes = cellfun(@(x,y) [(1:y(1))/(y(1)*2) ...
        0.5+((1:length(x)-y(1))/((length(x)-y(1))*2))], ...
        curr_rewardmoves,curr_rewardframes,'uni',false);
    
    % add leeway to movements to get accurate distributions
    leeway = 1;
    leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
        round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
    curr_rescaled_timeframes_leeway = cellfun(@(x,y,z) ...
        [linspace(-leeway,0.5,y+z(1)) linspace(0.5,1+leeway,length(x)-y+z(2))], ...
        curr_rewardmoves,curr_rewardframes,leeway_frames,'uni',false);
    curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
        curr_rewardmoves,leeway_frames,'uni',false);
    eliminate_movements = cellfun(@(x) ...
        any(x < 1 | x > length(movetime_trace)),curr_rewardmoves_leeway);
    curr_rewardmoves_leeway(eliminate_movements) = [];
    curr_rescaled_timeframes_leeway(eliminate_movements) = [];
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    curr_pyrclass = pyr_class_all{curr_animalday};
    %curr_pyrmoveact = pyr_activity_all{curr_animalday} ...
    %    (curr_pyrclass,horzcat(curr_rewardmoves{:}));
    curr_pyrmoveact_split = cellfun(@(x) pyr_activity_onset(:,x), ...
        curr_rewardmoves_leeway,'uni',false);
    curr_pyrmoveact = horzcat(curr_pyrmoveact_split{:}) > 0;
   
    %curr_pyrtime = curr_pyrmoveact.* ...
    %    repmat(horzcat(curr_timeframes{:}),size(curr_pyrmoveact,1),1);   
    %curr_pyrtime_mean = sum(curr_pyrtime,2)./sum(curr_pyrmoveact,2);
    
    curr_rescaled_pyrtime = curr_pyrmoveact.* ...
        repmat(horzcat(curr_rescaled_timeframes_leeway{:}),size(curr_pyrmoveact,1),1);   
    curr_rescaled_pyrtime_mean = sum(curr_rescaled_pyrtime,2)./sum(curr_pyrmoveact,2);
    
    %movement_time_all{curr_animalday} = curr_pyrtime_mean;
    movement_rescaled_time_all{curr_animalday} = curr_rescaled_pyrtime_mean;
    mean_activity_all{curr_animalday} = sum(curr_pyrmoveact,2)./ ...
        length(curr_rewardmoves_leeway);
    
    movement_time_std_all{curr_animalday} = arrayfun(@(x) ...
        std(curr_rescaled_pyrtime(x,curr_rescaled_pyrtime(x,:)~=0)), ...
        1:size(curr_rescaled_pyrtime,1));
     
    curr_pyrmoveact_split = cellfun(@(x) pyr_activity_all ...
        {curr_animalday}(:,x),curr_rewardmoves_leeway,'uni',false);
    %curr_pyrmoveact_split = cellfun(@(x) pyr_activity_onset(:,x), ...
    %    curr_rewardmoves_leeway,'uni',false);
    
    bin_frac = 0.05;
    timefrac_bins = [-leeway:bin_frac:1+leeway-bin_frac Inf];
    [n curr_timeframes_bins] = cellfun(@(x) histc(x,timefrac_bins), ...
        curr_timeframes,'uni',false);
    [n curr_rescaled_timeframes_bins] = cellfun(@(x) histc(x,timefrac_bins), ...
        curr_rescaled_timeframes_leeway,'uni',false);
    
    full_bins = cellfun(@(x) length(unique(x)) == ...
        length(timefrac_bins)-1, curr_rescaled_timeframes_bins);
    
    curr_pyrmoveact_timerescale = cellfun(@(x,y) ...
        grpstats(x(curr_pyrclass,:)',y)', ...
        curr_pyrmoveact_split(full_bins), ...
        curr_rescaled_timeframes_bins(full_bins),'uni',false);
    
    cat_pyrmoveact_timerescale = permute(cat(3,curr_pyrmoveact_timerescale{:}),[3 2 1]);
    pyrmoveact_timerescale_full_all{curr_animalday} = curr_pyrmoveact_timerescale;
    %pyrmoveact_timerescale_all{curr_animalday} = vertcat(curr_pyrmoveact_timerescale{:});
    pyrmoveact_timerescale_all{curr_animalday} = cat_pyrmoveact_timerescale;
    pyrmoveact_timerescale_mean_all{curr_animalday} = permute(nanmean(cat_pyrmoveact_timerescale,1),[3 2 1]);

        
    %fh = fspecial('gaussian',[1 20],2);    
    %j3 = conv2(cat(1,curr_pyrmoveact_timerescale{:}),fh,'same');
    %pyrmoveact_timerescale_mean_all{curr_animalday} = nanmean(j3);
    
    curr_animalday/max(animaldays)
    
%     spread_frames = 1;
%     spread_filter = ones(1,spread_frames);
%     pyr_spread = conv2(pyr_activity_all{curr_animalday}(curr_pyrclass,:),spread_filter,'same');
%     curr_move_onset_frames = cellfun(@(x) x(1), curr_rewardmoves);
%     curr_reward_frames = reward_frames_all{curr_animalday};
%     frames_back = 30;
%     frames_forward = 90;
%     hc = figure;
%     set(hc,'Name','Aligned')
%     n_subs = ceil(sqrt(length(curr_pyrtime_mean)));
%     [asdf sort_idx] = sort(curr_pyrtime_mean);
%     aligned_mean = nan(length(curr_pyrtime_mean),frames_back+frames_forward+1);
%     for curr_cell = 1:length(curr_pyrtime_mean)
%         % plot around move_onset
%         if curr_pyrtime_mean(curr_cell) < 2
%             curr_align = round(curr_move_onset_frames);
%         else
%             curr_align = round(curr_reward_frames);
%         end
%         activity_align = nan(length(curr_align),frames_back+frames_forward+1);
%         for i = 1:length(curr_align);          
%             align_frames = curr_align(i)-frames_back:curr_align(i)+frames_forward;
%             if any(align_frames<=0) || any(align_frames>size(pyr_spread,2));
%                 continue
%             end
%             activity_align(i,:) = pyr_spread(curr_cell,align_frames);
%         end
%         figure(hc);
%         subplot(n_subs,n_subs,curr_cell);%find(sort_idx==curr_cell));
%         imagesc(activity_align);
%         %imagesc(cat_pyrmoveact_timerescale(:,:,curr_cell));
%         colormap(gray);
%         
%         aligned_mean(find(sort_idx == curr_cell),:) = nanmean(activity_align);
%     end
    
end



figure
day_combine = {[2 3] [4 5] [6 7] [8 9] [10 11] [12 13 14]};
for i = 1:length(day_combine)
    subplot(length(day_combine),1,i)
    
    curr_data = vertcat(movement_rescaled_time_all{:,day_combine{i}});
    bin_edges = 0:0.05:1;
    [n bins] = histc(curr_data,bin_edges);
    
    curr_activity = vertcat(mean_activity_all{:,day_combine{i}});
    binned_activity = grpstats(curr_activity,bins);
    binned_activity_sem = grpstats(curr_activity,bins,'sem');
    
    use_bins = unique(bins);
    errorbar(use_bins,binned_activity,binned_activity_sem,'k','linewidth',2);
    
    %bar(use_bins(use_bins ~= 0),n(use_bins ~= 0),'k');
    xlim([0 length(bin_edges)-1]);
    set(gca,'XTick',2:2:length(bin_edges)-1);
    set(gca,'XTickLabel',bin_edges(3:2:end));
end

round_timefrac_bins = round(timefrac_bins*100)/100;
figure
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
for i = 1:length(day_combine)
    subplot(1,length(day_combine),i)
    
    curr_timing = vertcat(movement_rescaled_time_all{:,day_combine{i}});
    curr_data = vertcat(pyrmoveact_timerescale_all{:,day_combine{i}});
    curr_data_norm = bsxfun(@times,curr_data,1./max(curr_data,[],2));
    [asdf sort_idx] = sort(vertcat(pc1_scores{day_combine{i}}));
    imagesc(curr_data_norm(sort_idx,:));colormap(gray);
    %caxis([0 0.01]);
    
    line([find(round_timefrac_bins == 0) find(round_timefrac_bins == 0)],ylim,'color','r')
    line([find(round_timefrac_bins == 0.5) find(round_timefrac_bins == 0.5)],ylim,'color','r')
    line([find(round_timefrac_bins == 1) find(round_timefrac_bins == 1)],ylim,'color','r')
    
end

figure
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
bb = nan(length(day_combine),60);
for i = 1:length(day_combine)
    subplot(1,length(day_combine),i)
    
    curr_data = vertcat(pyrmoveact_timerescale_mean_all{:,day_combine{i}});
    curr_data_norm = bsxfun(@times,curr_data,1./max(curr_data,[],2));
    [cm cmi] = max(curr_data_norm,[],2);
    [asdf sort_idx] = sort(cmi);
    imagesc(curr_data_norm(sort_idx,:));%colormap(gray);
    
    bb(i,:) = nanmean(curr_data);
end

figure;
use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
for i = 1:length(day_combine)
    subplot(length(day_combine),2,(i*2)-1)
    bin_edges = -leeway:0.05:1+leeway;
    curr_time = vertcat(movement_rescaled_time_all{use_animals,day_combine{i}});
    [n bins] = histc(curr_time,bin_edges);
    use_bins = unique(bins);
    bar(bin_edges,n,'k');
    xlim([min(bin_edges) max(bin_edges)]);
    
    subplot(length(day_combine),2,2*i)
    curr_act = vertcat(mean_activity_all{use_animals,day_combine{i}});
    use_bins = unique(bins); 
    binned_act = grpstats(curr_act,bins);
    binned_act_sem = grpstats(curr_act,bins,'sem');
    errorbar(use_bins(use_bins ~= 0),binned_act(use_bins ~= 0), ...
        binned_act_sem(use_bins ~= 0),'k','linewidth',2);
    xlim([0 length(bin_edges)]);
end


% Plot std
movement_time_std_movecells = cellfun(@(x,y) x(y),...
    movement_time_std_all,pyr_class_all,'uni',false);
use_animals = [1 2 4 5 6 7 8];
std_cat = cell(1,14);
for i = 1:14
    std_cat{i} = horzcat(movement_time_std_movecells{use_animals,i});
end
figure;
errorbar(cellfun(@nanmedian,std_cat),cellfun(@(x) ...
    nanstd(x)/sqrt(sum(~isnan(x))),std_cat),'.k','linewidth',2)

%% Get time preferences for cells based on all movements
% movement_trace_all
% move_epoch_frames
% rewarded_movements
% cued_rewarded_movements
% 
% cue_frames_all
% reward_frames_all
% cued_reward_frames_all

movement_time_all = cell(size(movement_trace_all));
mean_activity_all = cell(size(movement_trace_all));

animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';

pyrmoveact_all = cell(size(movement_trace_all));


for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    % TEMPORARY: take curr moves start +/- 15 frames
    curr_moves = cellfun(@(x) x(1)-15:x(1)+15,curr_moves,'uni',false);
    elim_moves = cellfun(@(x) any(x < 1 | x > length(movetime_trace)),curr_moves);
    curr_moves(elim_moves) = [];
    
    curr_timeframes = cellfun(@(x) (1:length(x))/length(x), ...
        curr_moves,'uni',false);
  
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    curr_pyrclass = pyr_class_all{curr_animalday};
    %curr_pyrmoveact = pyr_activity_all{curr_animalday} ...
    %    (curr_pyrclass,horzcat(curr_moves{:}));
    curr_pyrmoveact = pyr_activity_onset ...
        (curr_pyrclass,horzcat(curr_moves{:})) > 0;
    
    curr_pyrtime = curr_pyrmoveact.* ...
        repmat(horzcat(curr_timeframes{:}),size(curr_pyrmoveact,1),1);   
    curr_pyrtime_mean = sum(curr_pyrtime,2)./sum(curr_pyrmoveact,2);
    
    %curr_pyrmoveact_split = cellfun(@(x) pyr_activity_all ...
    %    {curr_animalday}(:,x),curr_moves,'uni',false);
    curr_pyrmoveact_split = cellfun(@(x) pyr_activity_onset(curr_pyrclass,x), ...
        curr_moves,'uni',false);
    
    % get timed activity by first activity within movement
    [active start_idx] = cellfun(@(x) max(x>0,[],2),curr_pyrmoveact_split,'uni',false);
    move_times = cellfun(@(x,y,z) x(y),curr_timeframes,start_idx,'uni',false);
    
    cat_move_times = vertcat(move_times{:})';
    cat_active = horzcat(active{:});
    
    cat_move_times(~cat_active) = NaN;
    
    movement_time_all{curr_animalday} = nanmean(cat_move_times,2);
    mean_activity_all{curr_animalday} = nanmean(cat_active,2);
    
    
%     % store activity 
%     curr_pyrmoveact_split = cellfun(@(x) pyr_activity_all{curr_animalday}(:,x), ...
%        curr_moves,'uni',false);
%     pyrmoveact_all{curr_animalday} = permute(cat(3,curr_pyrmoveact_split{:}),[3,2,1]);

    curr_animalday/max(animaldays)

end


figure;
use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
for i = 1:length(day_combine)
    subplot(length(day_combine),2,(i*2)-1)
    bin_edges = 0:0.05:1;
    curr_time = vertcat(movement_time_all{use_animals,day_combine{i}});
    [n bins] = histc(curr_time,bin_edges);
    use_bins = unique(bins);
    bar(bin_edges,n,'k');
    xlim([0 max(bin_edges)]);
    
    subplot(length(day_combine),2,2*i)
    bin_edges = 0:0.05:1;
    curr_act = vertcat(mean_activity_all{use_animals,day_combine{i}});
    use_bins = unique(bins); 
    binned_act = grpstats(curr_act,bins);
    binned_act_sem = grpstats(curr_act,bins,'sem');
    errorbar(use_bins(use_bins ~= 0),binned_act(use_bins ~= 0), ...
        binned_act_sem(use_bins ~= 0),'k','linewidth',2);
    xlim([0 length(bin_edges)]);
end



figure;hold on;
use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
col = jet(length(day_combine));
for i = 1:length(day_combine)
    bin_edges = 0:0.05:1;
    curr_time = vertcat(movement_time_all{use_animals,1:14});
    [n bins] = histc(curr_time,bin_edges);
    use_bins = unique(bins);

    bin_edges = 0:0.05:1;
    curr_act = vertcat(mean_activity_all{use_animals,1:14});
    use_bins = unique(bins); 
    binned_act = grpstats(curr_act,bins);
    binned_act_sem = grpstats(curr_act,bins,'sem');
    errorbar(use_bins(use_bins ~= 0),binned_act(use_bins ~= 0), ...
        binned_act_sem(use_bins ~= 0),'color',col(i,:),'linewidth',2);
end




% TEMPORARY: to look at movement onset
pyrmoveact_mean = cellfun(@(x) nanmean(x,3),pyrmoveact_all,'uni',false);
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
sig_act = nan(length(day_combine),1);
for i = 1:length(day_combine)
    curr_data = vertcat(pyrmoveact_mean{:,day_combine{i}});
    curr_premove = curr_data(:,1:5);
    curr_prctile = prctile(bootstrp(10000,@nanmean,curr_premove(:)),99);
    over_thresh = +(nanmean(curr_data) > curr_prctile);
    sig_act(i) = find(conv(over_thresh,[1 1],'valid') == 2,1);    
    i
end


%% Get aligned activity +/- frames
% movement_trace_all
% move_epoch_frames
% rewarded_movements
% cued_rewarded_movements
% 
% cue_frames_all
% reward_frames_all
% cued_reward_frames_all

animaldays = find(cellfun(@(x) ~isempty(x),movement_trace_all))';

moveonset_aligned = cell(size(movement_trace_all));
reward_aligned = cell(size(movement_trace_all));
moveoffset_aligned = cell(size(movement_trace_all));

for curr_animalday = animaldays
    movetime_trace = nan(size(movement_trace_all{curr_animalday}));
    curr_moves = move_epoch_frames{curr_animalday};
    
    curr_rewardmoves = move_epoch_frames{curr_animalday}( ...
        rewarded_movements{curr_animalday});
    
    curr_moveonsets = cellfun(@(x) x(1),curr_rewardmoves,'uni',false);
    curr_moveoffsets = cellfun(@(x) x(end),curr_rewardmoves,'uni',false);
    
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{curr_animalday}),curr_rewardmoves,'uni',false);
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
    
    % Movement onset/offset reward +/- frames
    fs_b = 60;fs_f = 200;
    curr_moveonset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveonsets,'uni',false)';
    curr_reward_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_rewardframes,'uni',false)';
    curr_moveoffset_surround = cellfun(@(x) x-fs_b:x+fs_f,curr_moveoffsets,'uni',false)';
       
    elim_moves = cellfun(@(x,y,z) any([x y z]  < 1 | [x y z] > ...
        length(movetime_trace)),curr_moveonset_surround, ...
        curr_reward_surround,curr_moveoffset_surround);
    curr_moveonset_surround(elim_moves) = [];
    curr_reward_surround(elim_moves) = [];
    curr_moveoffset_surround(elim_moves) = [];
     
    % to get just activity onsets (ignore sustained, to peak)
    event_onsets = [false(size(pyr_activity_all{curr_animalday},1),1) ...
        diff(pyr_activity_all{curr_animalday}>0,[],2) == 1];
    pyr_activity_onset = zeros(size(pyr_activity_all{curr_animalday}));
    pyr_activity_onset(event_onsets) = pyr_activity_all{curr_animalday}(event_onsets);
    
    curr_pyrclass = pyr_class_all{curr_animalday};
      
    use_activity = pyr_activity_all{curr_animalday};%(curr_pyrclass,:);
    
    % store activity 
    curr_moveonset_aligned = cellfun(@(x) use_activity(:,x), ...
       curr_moveonset_surround,'uni',false);
    moveonset_aligned{curr_animalday} = nanmean(cat(3,curr_moveonset_aligned{:}),3);

    curr_reward_aligned = cellfun(@(x) use_activity(:,x), ...
       curr_reward_surround,'uni',false);
    reward_aligned{curr_animalday} = nanmean(cat(3,curr_reward_aligned{:}),3);
    
    curr_moveoffset_aligned = cellfun(@(x) use_activity(:,x), ...
       curr_moveoffset_surround,'uni',false);
    moveoffset_aligned{curr_animalday} = nanmean(cat(3,curr_moveoffset_aligned{:}),3);
    
    curr_animalday/max(animaldays)

end

% Make heatmap of all activity aligned to movement onset
use_animals = [1 2 4 5 6 7 8];
figure;
for i = 1:14
    subplot(1,14,i)
    curr_act = cellfun(@(x,y) x(y,:),moveonset_aligned(use_animals,i), ...
        pyr_class_all(use_animals,i),'uni',false);
    cat_data = vertcat(curr_act{:});
    cat_data_norm = bsxfun(@times,cat_data,1./max(cat_data,[],2));
    [c ci] = max(cat_data,[],2);
    [asdf sort_idx] = sort(ci);
    imagesc(cat_data_norm(sort_idx,:));colormap(hot);
    line([60 60],ylim,'color','b','linewidth',3)
    line([144 144],ylim,'color','b','linewidth',3)
    caxis([0 1]);
    axis off
end

%% Plot pyr:gad balance by movement

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
% pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);
gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   pyr_activity_all,move_epoch_frames,'uni',false);

use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);
gc = horzcat(g{use_animals,:});

bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
figure;
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement-based activity')


%% Movement/task related cells: history of classification

task_prior_movement = nan(8,15);
task_prior_task = nan(8,15);
for curr_animal = 1:8
    curr_reward = horzcat(reward_cells_all{curr_animal,:});
    curr_move_onset = horzcat(move_onset_cells_all{curr_animal,:});
    curr_task = curr_reward | curr_move_onset;    
    curr_movement = vertcat(pyr_class_all{curr_animal,:})';
    
    curr_task_diff = diff(curr_task,[],2);
    cur_movement_diff = diff(curr_movement,[],2);
    
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    task_prior_movement(curr_animal,sessions(2:end)) = ...
        sum(curr_task(:,2:end) & curr_movement(:,1:end-1))./sum(curr_task(:,2:end));
    task_prior_task(curr_animal,sessions(2:end)) = ...
        sum(curr_task(:,2:end) & curr_task(:,1:end-1))./sum(curr_task(:,2:end));

    
end

figure; hold on;
day_combine = {[2 3] [4 5] [6 7] [8 9] [10 11] [12 13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = task_prior_movement(:,day_combine{i});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
errorbar(data_mean,data_sem,'k','linewidth',2);
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = task_prior_task(:,day_combine{i});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
errorbar(data_mean,data_sem,'r','linewidth',2);
legend({'Previously movement-correlated','Previously task-aligned'});
title('Task-aligned classification history');
xlabel('Session')
ylabel('Fraction of task-aligned cells')
xlim([0 length(day_combine)+1])
set(gca,'XTick',1:length(day_combine))
set(gca,'XTickLabel',cellfun(@num2str,day_combine,'uni',false));



%% Task-aligned timing sequence/time preference std

task_aligned_all = cellfun(@(x,y,z) x(z)|y(z), ...
    reward_cells_all,move_onset_cells_all,pyr_class_all,'uni',false);
task_aligned_movement = cellfun(@(x,y) x(y), ...
    movement_time_all,task_aligned_all,'uni',false);
task_aligned_traces = cellfun(@(x,y) x(y,:), ...
    pyrmoveact_timerescale_all,task_aligned_all,'uni',false);

reward_aligned_all = cellfun(@(x,z) x(z), ...
    reward_cells_all,pyr_class_all,'uni',false);
reward_aligned_movement = cellfun(@(x,y) x(y), ...
    movement_time_all,reward_aligned_all,'uni',false);
reward_aligned_traces = cellfun(@(x,y) x(y,:), ...
    pyrmoveact_timerescale_all,reward_aligned_all,'uni',false);

move_onset_aligned_all = cellfun(@(x,z) x(z), ...
    move_onset_cells_all,pyr_class_all,'uni',false);
move_onset_aligned_movement = cellfun(@(x,y) x(y), ...
    movement_time_all,move_onset_aligned_all,'uni',false);
move_onset_aligned_traces = cellfun(@(x,y) x(y,:), ...
    pyrmoveact_timerescale_all,move_onset_aligned_all,'uni',false);

% Break them up by <,> 0.5 movement mean
after_reward_cells_all = cellfun(@(x) x >= 0.5,task_aligned_movement,'uni',false);
before_reward_cells_all = cellfun(@(x) x < 0.5,task_aligned_movement,'uni',false);

after_reward_aligned_movement = cellfun(@(x,y) x(y), ...
    task_aligned_movement,after_reward_cells_all,'uni',false);
after_reward_aligned_traces = cellfun(@(x,y) x(y,:), ...
    task_aligned_traces,after_reward_cells_all,'uni',false);

before_reward_aligned_movement = cellfun(@(x,y) x(y), ...
    task_aligned_movement,before_reward_cells_all,'uni',false);
before_reward_aligned_traces = cellfun(@(x,y) x(y,:), ...
    task_aligned_traces,before_reward_cells_all,'uni',false);


m1 = cellfun(@nanstd,reward_aligned_movement,'uni',false);
m2 = cellfun(@nanstd,move_onset_aligned_movement,'uni',false);


figure; hold on;
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(m1{:,day_combine{i}});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
errorbar(data_mean,data_sem,'k','linewidth',2);
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(m2{:,day_combine{i}});
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
errorbar(data_mean,data_sem,'r','linewidth',2);

%% Rank sequence and sequence correlation over time



use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
pyr_move_active = cell(size(use_animaldays));
pyr_move_time = cell(size(use_animaldays));

pyr_trial_rank_corr = cell(size(use_animaldays));
pyr_all_rank_corr = cell(size(use_animaldays));
pyr_time_corr = cell(size(use_animaldays));

num_cells = cell(size(use_animaldays));
all_reliability = cell(size(use_animaldays));
pyr_move_active_relative = cell(size(use_animaldays));
for i = find(use_animaldays)'
    
    curr_animal = mod(i,8);
    if curr_animal == 0;
       curr_animal = 8; 
    end
    %use_cells = pyr_class_all{i} & ...
    %    (move_onset_cells_all{i}' &...
    %    ~reward_cells_all{i}');
    %use_cells = pyr_class_all{i} & ...
    %    ~any(horzcat(move_onset_cells_all{curr_animal,:}),2)' &...
    %    any(horzcat(reward_cells_all{curr_animal,:}),2)';
    %use_cells = pyr_class_all{i} & ...
    %    any(horzcat(reward_cells_all{curr_animal,:}),2)';
    use_cells = pyr_class_all{i};
    use_cells_idx = find(use_cells);
    
    curr_rewardmoves = move_epoch_frames{i}( ...
        rewarded_movements{i});   
    curr_rewardframes = cellfun(@(x) intersect(x, ...
        reward_frames_all{i}) - x(1) + 1,curr_rewardmoves,'uni',false); 
    curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);    
     % add leeway to movements
    leeway = 0.25;
    %leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
    %    round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
    leeway_frames = cellfun(@(x,y) [10 ...
        10],curr_rewardmoves,curr_rewardframes,'uni',false);
    curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
        curr_rewardmoves,leeway_frames,'uni',false);
    eliminate_movements = cellfun(@(x) ...
        any(x < 1 | x > length(movement_trace_all{i})),curr_rewardmoves_leeway);
    curr_rewardmoves_leeway(eliminate_movements) = [];
    
    [curr_activity curr_time] = cellfun(@(x) ...
        max(pyr_activity_all{i}(use_cells,x)>0,[],2), ...
        curr_rewardmoves_leeway,'uni',false);
    
    curr_movement_lengths = cellfun(@length,...
        move_epoch_frames{i}(rewarded_movements{i}));
    
    pyr_move_active{i} = cell2mat(curr_activity);
    pyr_move_time{i} = cell2mat(curr_time);
    
    
    % try 3 different templates: 1) rank in trial, 2) rank in total,
    % 3) relative timing in movement
    [asdf trial_rank_template] = cellfun(@(x) ...
        sort(movement_rescaled_time_all{i}(use_cells_idx(x))), ...
        curr_activity,'uni',false);
    
    all_rank = tiedrank(movement_rescaled_time_all{i}(use_cells));
    all_rank_template = cellfun(@(x) all_rank(x),curr_activity,'uni',false);
    
    time_template = cellfun(@(x) movement_rescaled_time_all{i} ...
        (use_cells_idx(x)),curr_activity,'uni',false);
    
    % get actual rank
    [asdf curr_sort_idx] = cellfun(@(x,y) sort(y(x)), ...
        curr_activity,curr_time,'uni',false);
      
    trial_rank_corr = cellfun(@(x,y) corrcoef(x,y),trial_rank_template,curr_sort_idx,'uni',false);
    all_rank_corr = cellfun(@(x,y) corrcoef(x,y),all_rank_template,curr_sort_idx,'uni',false);
    time_corr = cellfun(@(x,y) corrcoef(x,y),time_template,curr_sort_idx,'uni',false);
    
    sort_corr_use = cellfun(@(x) ~any(isnan(x(:))),trial_rank_corr);   
    
    pyr_trial_rank_corr{i} = cellfun(@(x) x(2),trial_rank_corr(sort_corr_use));
    pyr_all_rank_corr{i} = cellfun(@(x) x(2),all_rank_corr(sort_corr_use));
    pyr_time_corr{i} = cellfun(@(x) x(2),time_corr(sort_corr_use));
    
    num_cells{i} = cellfun(@length,curr_sort_idx(sort_corr_use));
    all_reliability{i} = nanmean(pyr_move_active{i}(:,sort_corr_use),2);
    
    pyr_move_active_relative{i} = nanmean(pyr_move_active{i} ...
        (:,sort_corr_use))./curr_movement_lengths(sort_corr_use);
end

pyr_corr_mean = cellfun(@nanmean,pyr_trial_rank_corr);
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
use_animals = [1 2 4 5 6 7 8];
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    %curr_data = horzcat(pyr_trial_rank_corr{use_animals,day_combine{i}});
    curr_data = vertcat(pyr_corr_mean(use_animals,day_combine{i}));
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
figure;
errorbar(data_mean,data_sem,'--r','linewidth',2);


% visualize changes in sorting
figure; 
for curr_animal = 1:8;
subplot(1,8,curr_animal);hold on;
curr_sessions = find(cellfun(@(x) ...
    ~isempty(x),pyrmoveact_timerescale_all(curr_animal,:)));
total_rows = sum(cellfun(@(x) size(x,1),pyrmoveact_timerescale_all(curr_animal,:)));
row_idx = total_rows;
xlim([0 61]);ylim([0 total_rows+1])
for session = curr_sessions
    %subplot(max(curr_sessions),1,session); hold on
    use_cells_all = pyr_class_all{curr_animal,session}';
    use_cells = use_cells_all...
        (pyr_class_all{curr_animal,session});
    [asdf sort_idx] = sort(movement_rescaled_time_all ...
        {curr_animal,session}(use_cells_all));
    curr_data = pyrmoveact_timerescale_all{curr_animal,session}(:,:,use_cells);
    curr_data_sort = curr_data(:,:,sort_idx);
    col = jet(size(curr_data,3));
    for j = 1:size(curr_data,1)
        for i = 1:size(curr_data,3)
            %              plot(find(curr_data_sort(j,:,i)),row_idx* ...
            %                  ones(length(find(curr_data_sort(j,:,i))),1),'.', ...
            %                  'MarkerSize',5,'color',col(i,:));
            plot(find(curr_data_sort(j,:,i),1),row_idx* ...
                 ones(length(find(curr_data_sort(j,:,i),1)),1),'.', ...
                'MarkerSize',5,'color',col(i,:));
            %             curr_active = find(curr_data_sort(j,:,i));
            %             for k = 1:length(curr_active)
            %                line([curr_active(k) curr_active(k)], ...
            %                    [row_idx row_idx+0.5],'color',col(i,:),'linewidth',2);
            %             end
        end
        row_idx = row_idx - 1;
    end
    line(xlim,[row_idx-0.5 row_idx - 0.5],'color','k','linewidth',2);
end
end

% visualize changes in sorting by matrix
figure;
for curr_animal = 1:8
curr_sessions = find(cellfun(@(x) ...
    ~isempty(x),pyrmoveact_timerescale_all(curr_animal,:)));
total_rows = sum(cellfun(@(x) size(x,1),pyrmoveact_timerescale_all(curr_animal,:)));
active_plot = nan(total_rows,60,3);
row_idx = 1;
for session = curr_sessions
    %subplot(max(curr_sessions),1,session); hold on
    [asdf sort_idx] = sort(movement_rescaled_time_all ...
        {curr_animal,session}(pyr_class_all{curr_animal,session}));
    curr_data = pyrmoveact_timerescale_all{curr_animal,session};
    curr_data_sort = curr_data(:,:,sort_idx);
    sort_col = reshape(jet(size(curr_data,3)),[],1,3);
    for j = 1:size(curr_data,1)
        for i = 1:size(curr_data,3)
            curr_active = find(curr_data_sort(j,:,i));
            for k = 1:length(curr_active)
                active_plot(row_idx,curr_active(k),:) = ...
                    nanmean([active_plot(row_idx,curr_active(k),:)  ...
                    sort_col(i,:,:)]);
            end
%             active_plot(row_idx,:) = ...
%                 active_plot(row_idx,:) + ...
%                 (curr_data_sort(j,:,i) > 0)*sort_col(i);
        end
        row_idx = row_idx + 1;
    end
end
subplot(1,8,curr_animal);
imagesc(active_plot);
end


%% Get activity during rewarded / other movements

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   pyr_activity_all,move_epoch_frames,'uni',false);


dropped_all = cell(8,15);
retained_all = cell(8,15);
for curr_animal = 1:8
 
    curr_movement = vertcat(pyr_class_all{curr_animal,:})';
    
    dropped = curr_movement(:,1:end-1) & ~curr_movement(:,2:end);
    retained = curr_movement(:,1:end-1) & curr_movement(:,2:end);
    
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    dropped_all(curr_animal,sessions(2:end)) = ...
        mat2cell(dropped,size(dropped,1),ones(1,size(dropped,2)));
    
    retained_all(curr_animal,sessions(2:end)) = ...
        mat2cell(retained,size(retained,1),ones(1,size(retained,2)));

end

pyr_cuedrewarded_move_active = cellfun(@(x,y) nanmean(x(:,~y),2), ...
    pyr_move_active,cued_rewarded_movements,'uni',false);
pyr_rewarded_move_active = cellfun(@(x,y) nanmean(x(:,y),2), ...
    pyr_move_active,rewarded_movements,'uni',false);
pyr_nonrewarded_move_active = cellfun(@(x,y) nanmean(x(:,~y),2), ...
    pyr_move_active,rewarded_movements,'uni',false);


r_r_mean = cellfun(@(x,y) nanmean(x(y)), pyr_rewarded_move_active,retained_all);
d_r_mean = cellfun(@(x,y) nanmean(x(y)), pyr_rewarded_move_active,dropped_all);

r_n_mean = cellfun(@(x,y) nanmean(x(y)), pyr_nonrewarded_move_active,retained_all);
d_n_mean = cellfun(@(x,y) nanmean(x(y)), pyr_nonrewarded_move_active,dropped_all);

r_c_mean = cellfun(@(x,y) nanmean(x(y)), pyr_cuedrewarded_move_active,retained_all);
d_c_mean = cellfun(@(x,y) nanmean(x(y)), pyr_cuedrewarded_move_active,dropped_all);


%% Get gad-pyr correlations v reconstruction correlation


gad_moveact = cellfun(@(x,y) x(y,:),gad_activity_thresh_all,gad_class_all,'uni',false);
pyr_moveact = cellfun(@(x,y) x(y,:),pyr_activity_all,pyr_class_all,'uni',false);



all_corrcoef = cellfun(@(x,y) corrcoef([x' y']),pyr_moveact,gad_moveact,'uni',false);
all_corrcoef_pyrgad = cellfun(@(x,y) x(1:size(y,1),size(y,1)+1:end), ...
    all_corrcoef,pyr_moveact,'uni',false);
all_corrcoef_mean = cellfun(@(x) nanmean(x(:)),all_corrcoef_pyrgad);

all_pyr_weights = cellfun(@(x,y) y'\x',gad_moveact,pyr_moveact,'uni',false);
all_pyr_reconstruct = cellfun(@(x,y) x'*y,pyr_moveact,all_pyr_weights,'uni',false);
all_reconstruct_corr = cellfun(@(x,y) corrcoef([x' y]), ...
    gad_moveact,all_pyr_reconstruct,'uni',false);
all_reconstruct_corr_gadpyr = cellfun(@(x,y) x(1:size(y,1),size(y,1)+1:end), ...
    all_reconstruct_corr,gad_moveact,'uni',false);
all_reconstruct_corr_diag = cellfun(@(x) nanmean(diag(x)), ...
    all_reconstruct_corr_gadpyr);

figure;
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(all_corrcoef_mean(:,day_combine{i}));
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
subplot(2,1,1);
errorbar(data_mean,data_sem,'k','linewidth',2);
title('Mean move gad-pyr pairwise correlation')
set(gca,'XTick',[1:7])
set(gca,'XTickLabel',cellfun(@num2str,day_combine,'uni',false));
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(all_reconstruct_corr_diag(:,day_combine{i}));
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
subplot(2,1,2);
errorbar(data_mean,data_sem,'k','linewidth',2);
title('Mean gad-reconstructed pyr correlation')
set(gca,'XTick',[1:7])
set(gca,'XTickLabel',cellfun(@num2str,day_combine,'uni',false));







%% Similarity by PCA - eh

pc_dist = nan(8,15);
pc_corr = cell(8,15);
for curr_animal = 1:8;
sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
for curr_day = sessions;
% get coeffs from mean
a = pyrmoveact_timerescale_mean_all{curr_animal,curr_day};
a_mean_norm = bsxfun(@minus,a,nanmean(a,2));
a_norm = bsxfun(@times,a_mean_norm,1./std(a_mean_norm,[],2));
a_norm(isnan(a_norm)) = 0;
[coeffs scores latents] = princomp(a_norm');

% try blurring average or cleaner trajectory
%a_norm_blur = conv2(a_norm,fh','same');
%[coeffs scores latents] = princomp(a_norm_blur);

m = cat(3,pyrmoveact_timerescale_full_all{curr_animal,curr_day}{:});
m_mean_norm = bsxfun(@minus,m,nanmean(m,2));
m_std_norm = bsxfun(@times,m_mean_norm,1./nanstd(m_mean_norm,[],2));
m_std_norm(isnan(m_std_norm)) = 0;
m_movement_reshape = reshape(permute(m_std_norm,[2 3 1]),size(m_std_norm,2),[]);
fh = fspecial('gaussian',[20 1],2);fh = fh/max(fh);
m_blur = conv2(m_movement_reshape,fh,'same');
m_norm_blur = reshape(m_blur,[],size(m,1));

X = m_norm_blur*scores(:,1:3);
X1_reshape = reshape(X(:,1),size(m,2),[]);
X2_reshape = reshape(X(:,2),size(m,2),[]);
X3_reshape = reshape(X(:,3),size(m,2),[]);

X1_reshape_minsub = bsxfun(@minus,X1_reshape,min(X1_reshape));
X2_reshape_minsub = bsxfun(@minus,X2_reshape,min(X2_reshape));
X3_reshape_minsub = bsxfun(@minus,X3_reshape,min(X3_reshape));

X1_reshape_norm = bsxfun(@times,X1_reshape_minsub,1./max(X1_reshape_minsub));
X2_reshape_norm = bsxfun(@times,X2_reshape_minsub,1./max(X2_reshape_minsub));
X3_reshape_norm = bsxfun(@times,X3_reshape_minsub,1./max(X3_reshape_minsub));

corrcoef_grid = corrcoef([X1_reshape_norm coeffs(:,1);X2_reshape_norm ...
    coeffs(:,2); X3_reshape_norm coeffs(:,3)]);
% corrcoef_grid = corrcoef([X1_reshape_norm(21:41,:) coeffs(21:41,1);X2_reshape_norm(21:41,:) ...
%     coeffs(21:41,2); X3_reshape_norm(21:41,:) coeffs(21:41,3)]);
trial_corrcoef = corrcoef_grid(end,1:end-1);
pc_corr{curr_animal,curr_day} = trial_corrcoef;
% figure; hold on
% for i = 1:size(X1_reshape,2)
%     plot3(X1_reshape(:,i),X2_reshape(:,i),X3_reshape(:,i),'k');
% end

a1 = mean(X1_reshape,2)/max(mean(X1_reshape,2));
a2 = mean(X2_reshape,2)/max(mean(X2_reshape,2));
a3 = mean(X3_reshape,2)/max(mean(X3_reshape,2));
r = bsxfun(@times,coeffs(:,1:3),1./max(coeffs(:,1:3)));

% figure; hold on;
% plot3(a1,a2,a3,'k')
% plot3(r(:,1),r(:,2),r(:,3),'b')
% 
% plot3(a1(21:41),a2(21:41),a3(21:41),'k','linewidth',3)
% plot3(r(21:41,1),r(21:41,2),r(21:41,3),'b','linewidth',3)

% this is to clean up the template a little
a_norm_blur = conv2(a_norm,fh','same');
[coeffs scores latents] = princomp(a_norm_blur);
r = bsxfun(@times,coeffs(:,1:3),1./max(coeffs(:,1:3)));

pc_dist(curr_animal,curr_day) = nanmean(sqrt((a1(21:41)-r(21:41,1)).^2+ ...
     (a2(21:41)-r(21:41,2)).^2+(a3(21:41)-r(21:41,3)).^2));
%curr_corr = corrcoef([a1;a2;a2],reshape(r,length(r(:)),1));
%pc_dist(curr_animal,curr_day) = curr_corr(2);
curr_day
end
end

pc_corr_mean = cellfun(@nanmean,pc_corr);
use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
for i = 1:length(day_combine)
    curr_data = horzcat(pc_corr{use_animals,day_combine{i}});
    %curr_data = horzcat(pc_corr_mean(use_animals,day_combine{i}));
    data_mean(i) = nanmean(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
end
figure;
errorbar(data_mean,data_sem,'k','linewidth',2);






%% Similarity by correlation/distance
activity_distance = cell(8,15);
trial_corr = cell(8,15);
trial_corr_zeromean = cell(8,15);
template_origin_dist = cell(8,15);
activity_full_zeromean_all = cell(8,15);
activity_template_zeromean_all = cell(8,15);

for curr_animal = 1:8;
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_day = sessions;
        
        resample_frames = 60;
        
        %%%% Get activity
        movetime_trace = nan(size(movement_trace_all{curr_animal,curr_day}));
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_day}( ...
            rewarded_movements{curr_animal,curr_day});
        
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_day}) - x(1) + 1,curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        
        % add leeway to movements to get accurate distributions
        leeway = 1;
        leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
            round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false); 
        
        curr_prereward_leeway = cellfun(@(x,y,z) (x(1)-y(1)):(x(1)+z-1), ...
            curr_rewardmoves,leeway_frames,curr_rewardframes,'uni',false);
        curr_postreward_leeway = cellfun(@(x,y,z) (x(1)+z):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,curr_rewardframes,'uni',false);
        
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,'uni',false);
        
        eliminate_movements = cellfun(@(x,y) ...
            any([x y] < 1 | [x y] > length(movetime_trace)), ...
            curr_prereward_leeway,curr_postreward_leeway);
        curr_prereward_leeway(eliminate_movements) = [];
        curr_postreward_leeway(eliminate_movements) = [];
        curr_rewardmoves_leeway(eliminate_movements) = [];
       
        curr_pyrclass = pyr_class_all{curr_animal,curr_day};
     
        curr_pyrmoveact_split = cellfun(@(x) pyr_activity_all ...
            {curr_animal,curr_day}(:,x),curr_rewardmoves_leeway,'uni',false);
        
        % Get activity during movements
        curr_pyrmoveact_prereward_split = cellfun(@(x) pyr_activity_all ...
            {curr_animal,curr_day}(:,x)',curr_prereward_leeway,'uni',false);
        curr_pyrmoveact_postreward_split = cellfun(@(x) pyr_activity_all ...
            {curr_animal,curr_day}(:,x)',curr_postreward_leeway,'uni',false);
        
        % Pad movement activity for resampling
        curr_pyrmoveact_prereward_pad = cellfun(@(x) [repmat(x(1,:),100,1); ...
            x;repmat(x(end,:),100,1)],curr_pyrmoveact_prereward_split,'uni',false);
        curr_pyrmoveact_postreward_pad = cellfun(@(x) [repmat(x(1,:),100,1); ...
            x;repmat(x(end,:),100,1)],curr_pyrmoveact_postreward_split,'uni',false);
        
        % Resample padded movement activity
        curr_pyrmoveact_prereward_resample_pad = cellfun(@(x,y) resample(x, ...
            round(resample_frames/2),size(y,1)), ...
            curr_pyrmoveact_prereward_pad,curr_pyrmoveact_prereward_split,'uni',false);
        curr_pyrmoveact_postreward_resample_pad = cellfun(@(x,y) resample(x, ...
            round(resample_frames/2),size(y,1)), ...
            curr_pyrmoveact_postreward_pad,curr_pyrmoveact_postreward_split,'uni',false);
        
        % Take out center of resampled padded activity
        curr_pyrmoveact_prereward_resample = cellfun(@(x) ...
            x(round((size(x,1)-round(resample_frames/2))/2)+1 : ...
            round((size(x,1)-round(resample_frames/2))/2)+ ...
            round(resample_frames/2),:),curr_pyrmoveact_prereward_resample_pad,'uni',false);
         curr_pyrmoveact_postreward_resample = cellfun(@(x) ...
            x(round((size(x,1)-round(resample_frames/2))/2)+1 : ...
            round((size(x,1)-round(resample_frames/2))/2)+ ...
            round(resample_frames/2),:),curr_pyrmoveact_postreward_resample_pad,'uni',false);
        
         % TRY SHUFFLE: pull out movement time and shuffle
         curr_pyrmoveact_prereward_resample = cellfun(@(x) ...
             x(16:end,:),curr_pyrmoveact_prereward_resample,'uni',false);
         curr_pyrmoveact_postreward_resample = cellfun(@(x) ...
             x(1:15,:),curr_pyrmoveact_postreward_resample,'uni',false);
         
        
        cat_pyrmoveact_resample = nanmean([cat(3,curr_pyrmoveact_prereward_resample{:}); ...
            cat(3,curr_pyrmoveact_postreward_resample{:})]);
        
        activity_full = cat_pyrmoveact_resample;
        activity_full(isnan(activity_full)) = 0;
        activity_template = nanmean(activity_full,3);
        
        activity_full_zeromean = bsxfun(@minus,activity_full, ...
            mean(activity_full,1));
        activity_template_zeromean = nanmean(activity_full_zeromean,3);
        
        trial_corr_grid = corrcoef([activity_template(:) ...
            reshape(activity_full,[],size(activity_full,3))]);
        trial_corr_zeromean_grid = corrcoef([activity_template_zeromean(:) ...
            reshape(activity_full_zeromean,[],size(activity_full_zeromean,3))]);

        activity_distance{curr_animal,curr_day} = permute(sqrt(nansum((activity_full - ...
            repmat(activity_template,[1 1 size(activity_full,3)])).^2,2)),[1 3 2]);

        trial_corr{curr_animal,curr_day} = trial_corr_grid(1,2:end);   
        trial_corr_zeromean{curr_animal,curr_day} = trial_corr_zeromean_grid(1,2:end); 
        
        template_origin_dist{curr_animal,curr_day} = ...
            sqrt(sum(activity_template.^2,2));
        
        activity_full_zeromean_all{curr_animal,curr_day} = activity_full_zeromean;
        activity_template_zeromean_all{curr_animal,curr_day} = activity_template_zeromean;
    end
    curr_animal
end

pop_mean = cellfun(@nanmean,trial_corr_zeromean);
use_animals = [1 2 4 5 6 7 8];
%day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
day_combine = num2cell(1:14);
data_mean = nan(size(day_combine));
data_sem = nan(size(day_combine));
data_cat = cell(1,length(day_combine));
for i = 1:length(day_combine)
    curr_data = horzcat(trial_corr_zeromean{use_animals,day_combine{i}});
    %curr_data = horzcat(pop_mean(use_animals,day_combine{i}));
    data_mean(i) = nanmedian(curr_data(:));
    data_sem(i) = nanstd(curr_data(:))/sqrt(sum(~isnan(curr_data(:))));
    data_cat{i} = curr_data;
end
figure;
errorbar(data_mean,data_sem,'k','linewidth',2);
title('Trial correlation with average activity')


% use_animals = [1 2 4 5 6 7 8];
% day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
% col = jet(length(day_combine));
% figure; hold on;
% for i = 1:length(day_combine)
%     curr_data = horzcat(temp{use_animals,day_combine{i}});
%     %curr_data = horzcat(pop_mean(use_animals,day_combine{i}));
%     plot(nanmean(curr_data,2),'color',col(i,:));
% end
% 
% use_animals = [1 2 4 5 6 7 8];
% day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
% col = jet(length(day_combine));
% figure; hold on;
% for i = 1:length(day_combine)
%     curr_data = horzcat(template_origin_dist{use_animals,day_combine{i}});
%     %curr_data = horzcat(pop_mean(use_animals,day_combine{i}));
%     plot(nanmean(curr_data,2),'color',col(i,:));
% end

template_max = cellfun(@max,template_origin_dist,'uni',false);
activity_dist_norm = cellfun(@(x,y) x/y,activity_distance,template_max,'uni',false);
use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
col = jet(length(day_combine));
figure; hold on;
for i = 1:length(day_combine)
    curr_data = horzcat(activity_dist_norm{use_animals,day_combine{i}});
    %curr_data = horzcat(pop_mean(use_animals,day_combine{i}));
    plot(nanmean(curr_data,2),'color',col(i,:));
end
title('Normalized distance from template')

% Create correlation grid between days 
trial_correlation_days = cell(15,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    curr_templates_block = horzcat(activity_template_zeromean_all{curr_animal,:});
    curr_templates = reshape(curr_templates_block,[],length(sessions));
    curr_corr_grids = cellfun(@(x) corrcoef([curr_templates ...
        reshape(x,[],size(x,3))]),activity_full_zeromean_all( ...
        curr_animal,sessions),'uni',false);
    curr_corr_days = cellfun(@(x) x(1:length(sessions), ...
        length(sessions)+1:end),curr_corr_grids,'uni',false);
    curr_corr_split = cellfun(@(x) mat2cell(x, ...
        (ones(length(sessions),1)),size(x,2)),curr_corr_days,'uni',false);
    
    trial_correlation_days(sessions,sessions,curr_animal) = ...
        horzcat(curr_corr_split{:});  
    
    %mean_curr_corr_days = cell2mat(cellfun(@(x) nanmean(x,2),curr_corr_days,'uni',false));
    %trial_correlation_days(sessions,sessions,curr_animal) = mean_curr_corr_days;
end
use_animals = [1 2 4 5 6 7 8];
trial_correlation_cat = cell(14,14);
for i = 1:14
    for j = 1:14
        trial_correlation_cat{i,j} = ...
            horzcat(trial_correlation_days{i,j,use_animals});
    end
end
trial_corr_pop = cellfun(@nanmedian,trial_correlation_cat);
figure;imagesc(trial_corr_pop);colormap(hot)
title('Mean trial correlation with session template')






%% Gad correlation


use_animals = [1 2 4 5 6 7 8];

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) x(:,z),y,'uni',false)), ...
   gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) x(:,z),y,'uni',false)), ...
   pyr_activity_all,move_epoch_frames,'uni',false);

gad_corr = cellfun(@(x,y) corrcoef(x(y,:)'), ...
   gad_move_active,gad_class_all,'uni',false);
gad_pairmean_activity = cellfun(@(x,y) sqrt(nanmean(x(y,:)>0,2) * ...
    nanmean(x(y,:)>0,2)'),gad_move_active,gad_class_all,'uni',false);

use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
use_animaldays(setdiff(1:size(use_animaldays,1),use_animals),:) = false;

all_gad_corr = cellfun(@(x) AP_itril(x,-1), ...
    gad_corr(use_animaldays),'uni',false);
all_gad_pairmean_activity = cellfun(@(x) AP_itril(x,-1), ...
    gad_pairmean_activity(use_animaldays),'uni',false);

bin_edges = [0:0.1:0.9 Inf];
[n bins] = histc(vertcat(all_gad_pairmean_activity{:}),bin_edges);
binned_all_gad_corr = grpstats(vertcat(all_gad_corr{:}),bins);






%% Massimo question: get activity for trials with similar movement

trial_template_corr = cell(8,15);
trial_template_dot = cell(8,15);
trial_template_trials = cell(8,15);
for curr_animal = 1:8
   
    % Make template: average of last 3 days for animal
    curr_class = vertcat(pyr_class_all{curr_animal,:});
    curr_template = nanmean(curr_class(end-2:end,:));
    
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    for curr_session = sessions
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 10;
        %leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
        %    round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        cat_activity = horzcat(curr_activity{:});
        curr_template_corr = corrcoef([curr_template' cat_activity]);

        trial_template_corr{curr_animal,curr_session} = ...
            curr_template_corr(1,2:end);
        
        trial_template_trials{curr_animal,curr_session} = ...
            curr_trials;
        
        trial_template_dot{curr_animal,curr_session} = ...
            (curr_template*cat_activity)/sum(curr_template);
        
    end


end

use_animals = [1 2 4 5 6 7 8];
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
cat_data = cell(length(day_combine),1);
for i = 1:length(day_combine)
    cat_data{i} = horzcat(trial_template_corr{use_animals,day_combine{i}});    
end
figure;
for i = 1:length(day_combine)
    subplot(length(day_combine),1,i);
    hist(cat_data{i});
    xlim([0 1]);
end



%% Instead of above: go from movement to activity

activity_corr = cell(8,15);
cell_corr_lin = cell(8,15);
movement_activity = cell(8,15);
used_trials = cell(8,15);
% make activity template for each animal, mean of last 3 days
for curr_animal = 1:8
   
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    for curr_session = sessions
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        %leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
        %    round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
%         % get rid of frames longer than cutoff amount
%         frame_cutoff = 100;
%         long_moves = cellfun(@(x) length(x) > frame_cutoff,curr_rewardmoves_leeway);
%         curr_rewardmoves_leeway_cutoff = curr_rewardmoves_leeway;
%         curr_rewardmoves_leeway_cutoff(long_moves) = cellfun(@(x) x(1:frame_cutoff), ...
%             curr_rewardmoves_leeway(long_moves),'uni',false);
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        % get trials which have both lever correlation and activity
        use_moves = ismember(curr_trials,lever_trials_used{curr_animal,curr_session});
        use_curr_trials = curr_trials(use_moves);
        use_curr_activity = curr_activity(:,use_moves);
    end
    
end
for curr_animal = 1:8
   
    % Make template: average of last 3 days for animal
    curr_class = vertcat(pyr_class_all{curr_animal,:});
    curr_template = nanmean(curr_class(end-2:end,:));
   
    
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    for curr_session = sessions
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        %leeway_frames = cellfun(@(x,y) [round(y*leeway) ...
        %    round((length(x)-y)*leeway)],curr_rewardmoves,curr_rewardframes,'uni',false);
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
%         % get rid of frames longer than cutoff amount
%         frame_cutoff = 100;
%         long_moves = cellfun(@(x) length(x) > frame_cutoff,curr_rewardmoves_leeway);
%         curr_rewardmoves_leeway_cutoff = curr_rewardmoves_leeway;
%         curr_rewardmoves_leeway_cutoff(long_moves) = cellfun(@(x) x(1:frame_cutoff), ...
%             curr_rewardmoves_leeway(long_moves),'uni',false);
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        % get trials which have both lever correlation and activity
        use_moves = ismember(curr_trials,lever_trials_used{curr_animal,curr_session});
        use_curr_trials = curr_trials(use_moves);
        use_curr_activity = curr_activity(:,use_moves);
        
%         % check trials match
%         if any(use_curr_trials ~= lever_corr_trials{curr_animal,curr_session});
%             error('Trial numbers don''t match');
%         end

        cat_activity = horzcat(use_curr_activity{:});
        curr_template_corr_grid = corrcoef([act_template_1{curr_animal} cat_activity]);
        curr_template_corr = curr_template_corr_grid(1,2:end);
        activity_corr{curr_animal,curr_session} = curr_template_corr;
        
%         nan_trial = isnan(lever_corr_1{curr_animal,curr_session});
%         cell_corr_lin{curr_animal,curr_session} ...
%             = cat_activity(curr_template ~= 0,~nan_trial)'\ ...
%             lever_corr{curr_animal,curr_session}(~nan_trial)';
        
        
        movement_activity{curr_animal,curr_session} = cat_activity;
        used_trials{curr_animal,curr_session} = ismember( ...
            lever_corr_trials{curr_animal,curr_session},use_curr_trials);
    end
    curr_animal
end



use_animals = [1 2 4 5 6 7 8];
up_clusters = [2 3 NaN 2 2 3 1 3];

cat_dot = horzcat(activity_corr{use_animals,1:5});
template_cluster_cc
cat_cluster_cc = cell(8,15);
for i = 1:8
    sessions = find(cellfun(@(x) ~isempty(x),activity_corr(i,:)));
    cat_cluster_cc(curr_animal,sessions) = cellfun(@(x) ...
        x(up_clusters(i),:),template_cluster_cc,'uni',false);
end
cat_lever_corr = horzcat(cat_cluster_cc{use_animals,1:5});

corr_bins = [-1:0.25:1-0.25 Inf];
[n bins] = histc(cat_lever_corr,corr_bins);
lever_corr_bin_mean = grpstats(cat_dot(bins ~= 0),bins(bins ~= 0));
lever_corr_bin_sem = grpstats(cat_dot(bins ~= 0),bins(bins ~= 0),'sem');

figure;
errorbar(dot_bins(unique(bins(bins ~= 0))), ...
    lever_corr_bin_mean,lever_corr_bin_sem,'r','linewidth',2)


% get trial-trial correlation by binned lever correlation
corr_bins = [-1:0.25:1-0.25 Inf];
act_corr_binned = cell(8,15,length(corr_bins)-1,2);
lever_corr_bins = cell(8,15,length(corr_bins)-1,2);
for curr_animal = 1:8
    for curr_lever = 1:2
        if curr_lever == 1
            curr_lever_corr = lever_corr_1;
        elseif curr_lever == 2
            curr_lever_corr = lever_corr_2;
        end
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        for curr_session = sessions;
            curr_cat_act = horzcat(movement_activity{curr_animal,curr_session});
            lever_corr_use = cellfun(@(x,y) x(y),curr_lever_corr(curr_animal,curr_session), ...
                used_trials(curr_animal,curr_session),'uni',false);
            lever_corr_cat = horzcat(lever_corr_use{:});
            [n bins] = histc(lever_corr_cat,corr_bins);
            
            curr_lever_use = cellfun(@(x,y) x(y,:),lever_use(curr_animal,curr_session), ...
                used_trials(curr_animal,curr_session),'uni',false);
            curr_lever_use_cat = vertcat(curr_lever_use{:});
            
            for i = 1:length(corr_bins) - 1
                if sum(bins == i) == 0;
                    continue
                end
                %curr_corr_grid = +curr_cat_act(:,bins == i)'*+curr_cat_act(:,bins == i);
                curr_corr_grid = corrcoef(curr_cat_act(:,bins == i));
                act_corr_binned{curr_animal,curr_session,i,curr_lever} = ...
                    AP_itril(curr_corr_grid,-1);
                
                
                curr_lever_corr_grid = corrcoef(curr_lever_use_cat(bins == i,:)');
                lever_corr_bins{curr_animal,curr_session,i,curr_lever} = ...
                    AP_itril(curr_lever_corr_grid,-1);
            end
        end       
    end    
end

use_animals = [1 2 4 5 6 7 8];
cat_act_corr_binned = cell(14,2);
for i = 1:2
    for j = 1:14
    cat_act_corr_binned{j,i} = vertcat(act_corr_binned{use_animals,j,i});
    end
end

m = cellfun(@nanmean,act_corr_binned);
%m = cellfun(@nanmean,lever_corr_bins);
figure;
for i = use_animals
    subplot(3,3,i);
    plot(m(i,:,1),'k');
    hold on
    plot(m(i,:,2),'r');
end

figure;
col = jet(size(m,3));
for i = use_animals
    subplot(3,3,i);hold on;
    for j = 1:size(m,3);
        plot(m(i,:,j,1),'color',col(j,:));
    end    
end


% get trial-trial correlation by binned lever correlation
lever_act_corr = nan(8,15,3);
for curr_animal = 1:8
    for curr_lever = 1:2
        if curr_lever == 1
            curr_lever_corr = lever_corr_1;
        elseif curr_lever == 2
            curr_lever_corr = lever_corr_2;
        end
        curr_act = movement_activity(curr_animal,:);
        curr_act_corr = cellfun(@(x) nanmean(AP_itril(corrcoef(x),-1)),curr_act);
        
        lever_corr_use = cellfun(@(x,y) x(y),curr_lever_corr(curr_animal,:), ...
            used_trials(curr_animal,:),'uni',false);
        curr_lever_corr = cellfun(@nanmean,lever_corr_use);
        
        lever_act_corr(curr_animal,:,1) = curr_act_corr;
        lever_act_corr(curr_animal,:,1+curr_lever) = curr_lever_corr;
    end
end

figure;
for i = use_animals
    subplot(3,3,i);hold on
    plotyy(1:15,lever_act_corr(i,:,2),1:15,lever_act_corr(i,:,1));
    plot(lever_act_corr(i,:,3),'r');
end
figure; hold on
col = jet(8);
for i = use_animals
    plot(lever_act_corr(i,:,1),lever_act_corr(i,:,2),'o','color',col(i,:));
end


% new template: based on average activity from trials for lever template
act_template_1 = cell(8,1);
act_template_2 = cell(8,1);
act_template_1_corr = cell(8,15);
act_template_2_corr = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_act = movement_activity(curr_animal,sessions(end-2:end));
    
    curr_trials_1 = cellfun(@(x,y) ismember(find(x),y), ...
        used_trials(curr_animal,sessions(end-2:end))', ...
        template_trials{curr_animal,1},'uni',false);    
    act_template_1_all = cellfun(@(x,y) x(:,y), ...
        movement_activity(curr_animal,sessions(end-2:end))',curr_trials_1,'uni',false);
    curr_act_template_1 = nanmean(horzcat(act_template_1_all{:}),2);    
    act_template_1_corr_grid = cellfun(@(x) corrcoef([curr_act_template_1 x]), ...
        movement_activity(curr_animal,sessions),'uni',false);    
    act_template_1_corr(curr_animal,sessions) = ...
        cellfun(@(x) x(1,2:end),act_template_1_corr_grid,'uni',false);
    act_template_1{curr_animal} = curr_act_template_1;
    
    curr_trials_2 = cellfun(@(x,y) ismember(find(x),y), ...
        used_trials(curr_animal,sessions(end-2:end))', ...
        template_trials{curr_animal,2},'uni',false);    
    act_template_2_all = cellfun(@(x,y) x(:,y), ...
        movement_activity(curr_animal,sessions(end-2:end))',curr_trials_2,'uni',false);
    curr_act_template_2 = nanmean(horzcat(act_template_2_all{:}),2);    
    act_template_2_corr_grid = cellfun(@(x) corrcoef([curr_act_template_2 x]), ...
        movement_activity(curr_animal,sessions),'uni',false);    
    act_template_2_corr(curr_animal,sessions) = ...
        cellfun(@(x) x(1,2:end),act_template_2_corr_grid,'uni',false);   
    act_template_2{curr_animal} = curr_act_template_2;
end

a1 = cellfun(@nanmean,act_template_1_corr);
a2 = cellfun(@nanmean,act_template_2_corr);
figure
for i = use_animals
    subplot(3,3,i);hold on;
    plot(a1(i,:),'k');
    plot(a2(i,:),'r');
end


%% template correlation: binary trial activity

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    curr_animal
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);


subsets = 5;
corr_bins = [-1:0.1:1-0.1 Inf];
m = nan(subsets,length(corr_bins)-1);
s = nan(subsets,length(corr_bins)-1);
for curr_subset = 1:subsets
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_col = cell(8,1);

    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        %use_sessions = sessions(end-6:end-3);
        use_sessions = intersect(1:5,sessions);
        
        curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
            setdiff(1:size(x,2),y)), ...
            template_cluster_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        

        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;

    end
    
    
    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    curr_subset
end
figure;errorbar(corr_bins(1:end-1),nanmean(m,1),nanmean(s,1),'k','linewidth',2);

%---------------------------------
% another control from above: do same but with template made randomly
% for up template: look at trials with high lever correlation for activity

% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step

random_template_trials = cell(8,15);

up_clusters = [2 3 NaN 3 1 2 2 2];
activity_corr = cell(8,15);
act_template_trials = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_template_trials = cellfun(@(x) randsample(length(x), ...
        sum(x == up_clusters(curr_animal))), ...
        template_trials_crop(curr_animal,sessions),'uni',false);
    random_template_trials(curr_animal,sessions) = curr_template_trials;
    act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
        round(length(x)/2))),curr_template_trials,'uni',false);
    movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
        movement_activity_crop(curr_animal,sessions), ...
        act_template_trials(curr_animal,sessions),'uni',false);
    activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
    for curr_session = sessions
        cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
        curr_template_corr_grid = corrcoef([activity_template cat_activity]);
        curr_template_corr = curr_template_corr_grid(1,2:end);
        activity_corr{curr_animal,curr_session} = curr_template_corr;
    end
end

random_lever_template = cell(8,1);
random_lever_template_cc = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_lever_use = cellfun(@(x,y) x(y,:),lever_used_crop(curr_animal,:), ...
        random_template_trials(curr_animal,:),'uni',false);
    random_lever_template{curr_animal} = nanmean(vertcat(curr_lever_use{:}));
    curr_lever_cc_grid = cellfun(@(x) corrcoef([random_lever_template{curr_animal}' ...
        x']),lever_used_crop(curr_animal,:),'uni',false);
    curr_lever_cc = cellfun(@(x) x(1,2:end),curr_lever_cc_grid(sessions),'uni',false);
    random_lever_template_cc(curr_animal,sessions) = curr_lever_cc;
end

use_animals = [1 2 4 5 6 7 8];
total_lever_cc = cell(8,1);
total_act_cc = cell(8,1);
total_col = cell(8,1);
%figure; hold on
for curr_animal = use_animals;
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    use_sessions = sessions(end-2:end);
    
    curr_lever_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
        random_lever_template_cc(curr_animal,use_sessions), ...
        act_template_trials(curr_animal,use_sessions),'uni',false);
    cat_lever_cc = horzcat(curr_lever_cc{:});
    
    curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
        activity_corr(curr_animal,use_sessions), ...
        act_template_trials(curr_animal,use_sessions),'uni',false);
    cat_act_cc = horzcat(curr_act_cc{:});
    
    col = cellfun(@(x,y) ones(length(x),1)*y,curr_lever_cc, ...
        num2cell(1:length(curr_lever_cc)),'uni',false);
    x = 1:length(use_sessions);
    %scatter(cat_act_cc,cat_lever_cc,10,x);
    total_lever_cc{curr_animal} = cat_lever_cc;
    total_act_cc{curr_animal} = cat_act_cc;
    total_col{curr_animal} = vertcat(col{:});
end
corr_bins = [-1:0.1:1-0.1 Inf];
total_lever_cc_cat = horzcat(total_lever_cc{:});  
total_act_cc_cat = horzcat(total_act_cc{:});
[n bins] = histc(total_lever_cc_cat,corr_bins);
m = grpstats(total_act_cc_cat,bins);
s = grpstats(total_act_cc_cat,bins,'sem');
figure;errorbar(m,s,'--k','linewidth',2);



%% TK NYSCF figure

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x:x+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);


subsets = 100;
m_template = nan(subsets,4);
s_template = nan(subsets,4);
m_pairwise = nan(subsets,6);
s_pairwise = nan(subsets,6);
for curr_subset = 1:100
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,2);
    total_act_cc = cell(8,2);
    total_act = cell(8,2);
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        
        for earlylate = 1:2
            if earlylate == 1
                use_sessions = sessions(1:5);
            elseif earlylate == 2
                use_sessions = sessions(end-2:end);
            end
            
            curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
                setdiff(1:size(x,2),y)), ...
                template_cluster_cc_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_lever_cc = horzcat(curr_lever_cc{:});
            
            curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
                activity_corr(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_act_cc = horzcat(curr_act_cc{:});
            
            curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
                movement_activity_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            curr_act_cat = horzcat(curr_act{:});
            
            total_act{curr_animal,earlylate} = curr_act_cat;
            total_lever_cc{curr_animal,earlylate} = cat_lever_cc;
            total_act_cc{curr_animal,earlylate} = cat_act_cc;
        end
    end
    
    total_act_pairwise = cellfun(@corrcoef,total_act,'uni',false);
    total_act_pairwise_cross = cell(8,1);
    for i = 1:8
        curr_cat_act = horzcat(total_act{i,:});
        curr_pairwise_cc = corrcoef(curr_cat_act);
        total_act_pairwise_cross{i} = curr_pairwise_cc( ...
            1:size(total_act{i,1},2),size(total_act{i,1},2)+1:end);
    end
    
    total_lever_cc_cat_early = horzcat(total_lever_cc{:,1});
    total_act_cc_cat_early = horzcat(total_act_cc{:,1});
    
    total_lever_cc_cat_late = horzcat(total_lever_cc{:,2});
    total_act_cc_cat_late = horzcat(total_act_cc{:,2});
    
    dissimilar_animal = cellfun(@(x) x > -0.2 & ...
        x < 0.2,total_lever_cc,'uni',false);
    similar_animal = cellfun(@(x) x > 0.6, ...
        total_lever_cc,'uni',false);
  
    early_dissimilar = horzcat(dissimilar_animal{:,1});
    early_similar = horzcat(similar_animal{:,1});
    late_dissimilar = horzcat(dissimilar_animal{:,2});
    late_similar = horzcat(similar_animal{:,2});
    
    % get pairwise activity comparison
    
    dissimilar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,dissimilar_animal,'uni',false);
    similar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,similar_animal,'uni',false);
    simdis_pairwise = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise,similar_animal,dissimilar_animal,'uni',false);
    simsim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    dissim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    
    dissimilar_pairwise_early = vertcat(dissimilar_pairwise{:,1});
    similar_pairwise_early = vertcat(similar_pairwise{:,1});
    dissimilar_pairwise_late = vertcat(dissimilar_pairwise{:,2});
    similar_pairwise_late = vertcat(similar_pairwise{:,2});
    simdis_pairwise_early = vertcat(simdis_pairwise{:,1});
    simsim_pairwise_earlylate = vertcat(simsim_pairwise_cross{:});
    dissim_pairwise_earlylate = vertcat(dissim_pairwise_cross{:});
    
    % the pairs of within-comparison
    dissimilar_pairwise_early_mean = nanmean(dissimilar_pairwise_early);
    similar_pairwise_early_mean = nanmean(similar_pairwise_early);
    dissimilar_pairwise_late_mean = nanmean(dissimilar_pairwise_late);
    similar_pairwise_late_mean = nanmean(similar_pairwise_late);
    
    dissimilar_pairwise_early_sem = nanstd(dissimilar_pairwise_early)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_early)));
    similar_pairwise_early_sem = nanstd(similar_pairwise_early)/ ...
        sqrt(sum(~isnan(similar_pairwise_early)));
    dissimilar_pairwise_late_sem = nanstd(dissimilar_pairwise_late)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_late)));
    similar_pairwise_late_sem = nanstd(similar_pairwise_late)/ ...
        sqrt(sum(~isnan(similar_pairwise_late)));
    
    % the pairs of cross-comparison
    simdis_pairwise_early_mean = nanmean(simdis_pairwise_early);
    simsim_pairwise_earlylate_mean = nanmean(simsim_pairwise_earlylate);
    dissim_pairwise_earlylate_mean = nanmean(dissim_pairwise_earlylate);
    
    simdis_pairwise_early_sem = nanstd(simdis_pairwise_early)/ ...
        sqrt(sum(~isnan(simdis_pairwise_early)));
    simsim_pairwise_earlylate_sem = nanstd(simsim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(simsim_pairwise_earlylate)));
    dissim_pairwise_earlylate_sem = nanstd(dissim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(dissim_pairwise_earlylate)));
    
    m_pairwise(curr_subset,:) = [dissimilar_pairwise_early_mean similar_pairwise_early_mean ...
        simdis_pairwise_early_mean dissim_pairwise_earlylate_mean ...
        simsim_pairwise_earlylate_mean similar_pairwise_late_mean];
    s_pairwise(curr_subset,:) = [dissimilar_pairwise_early_sem similar_pairwise_early_sem ...
        simdis_pairwise_early_sem dissim_pairwise_earlylate_sem ...
        simsim_pairwise_earlylate_sem similar_pairwise_late_sem];
    
    % get activity comparison to template means
    
    early_dissimilar_cc_mean = nanmean(total_act_cc_cat_early(early_dissimilar));
    early_similar_cc_mean = nanmean(total_act_cc_cat_early(early_similar));
    late_dissimilar_cc_mean = nanmean(total_act_cc_cat_late(late_dissimilar));
    late_similar_cc_mean = nanmean(total_act_cc_cat_late(late_similar));
    
    early_dissimilar_cc_sem = nanstd(total_act_cc_cat_early(early_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_dissimilar))));
    early_similar_cc_sem = nanstd(total_act_cc_cat_early(early_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_similar))));
    late_dissimilar_cc_sem = nanstd(total_act_cc_cat_late(late_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_dissimilar))));
    late_similar_cc_sem = nanstd(total_act_cc_cat_late(late_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_similar))));
    
    m_template(curr_subset,:) = [early_dissimilar_cc_mean early_similar_cc_mean ...
        late_dissimilar_cc_mean late_similar_cc_mean];
    s_template(curr_subset,:) = [early_dissimilar_cc_sem early_similar_cc_sem ...
        late_dissimilar_cc_sem late_similar_cc_sem];
    

    
    disp(curr_subset);
end

m2 = nanmean(m_template);
s2 = nanmean(s_template);
figure; hold on
bar(m2([1 2 4]),'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2([1 2 4]),s2([1 2 4]),'.k','linewidth',2)
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early dissimilar lever' 'Early similar lever' 'Late similar lever'});
xlabel('Lever correlation')
ylabel('Template activity correlation')
ylim([0.15 0.5]);

m2 = nanmean(m_pairwise);
s2 = nanmean(s_pairwise);
figure; hold on
bar(m2,'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2,s2,'.k','linewidth',2)
set(gca,'XTick',1:6);
set(gca,'XTickLabel',{'Early dissimilar' 'Early similar' ...
    'Early dissimilar/similar' 'Late/early dissimilar' ...
    'Late/early similar' 'Late similar'});
xlabel('Lever correlation')
ylabel('Pairwise activity correlation')
ylim([0.15 0.30])
set(gca,'YTick',0.15:0.05:0.3)

%% template corrleation: include temporal information 

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);


subsets = 100;
m_template = nan(subsets,4);
s_template = nan(subsets,4);
m_pairwise = nan(subsets,6);
s_pairwise = nan(subsets,6);
for curr_subset = 1:100
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,2);
    total_act_cc = cell(8,2);
    total_act = cell(8,2);
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        
        for earlylate = 1:2
            if earlylate == 1
                use_sessions = sessions(1:5);
            elseif earlylate == 2
                use_sessions = sessions(end-2:end);
            end
            
            curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
                setdiff(1:size(x,2),y)), ...
                template_cluster_cc_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_lever_cc = horzcat(curr_lever_cc{:});
            
            curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
                activity_corr(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_act_cc = horzcat(curr_act_cc{:});
            
            curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
                movement_activity_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            curr_act_cat = horzcat(curr_act{:});
            
            total_act{curr_animal,earlylate} = curr_act_cat;
            total_lever_cc{curr_animal,earlylate} = cat_lever_cc;
            total_act_cc{curr_animal,earlylate} = cat_act_cc;
        end
    end
    
    total_act_pairwise = cellfun(@corrcoef,total_act,'uni',false);
    total_act_pairwise_cross = cell(8,1);
    for i = 1:8
        curr_cat_act = horzcat(total_act{i,:});
        curr_pairwise_cc = corrcoef(curr_cat_act);
        total_act_pairwise_cross{i} = curr_pairwise_cc( ...
            1:size(total_act{i,1},2),size(total_act{i,1},2)+1:end);
    end
    
    total_lever_cc_cat_early = horzcat(total_lever_cc{:,1});
    total_act_cc_cat_early = horzcat(total_act_cc{:,1});
    
    total_lever_cc_cat_late = horzcat(total_lever_cc{:,2});
    total_act_cc_cat_late = horzcat(total_act_cc{:,2});
    
    dissimilar_animal = cellfun(@(x) x > -0.2 & ...
        x < 0.2,total_lever_cc,'uni',false);
    similar_animal = cellfun(@(x) x > 0.6, ...
        total_lever_cc,'uni',false);
  
    early_dissimilar = horzcat(dissimilar_animal{:,1});
    early_similar = horzcat(similar_animal{:,1});
    late_dissimilar = horzcat(dissimilar_animal{:,2});
    late_similar = horzcat(similar_animal{:,2});
    
    % get pairwise activity comparison
    
    dissimilar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,dissimilar_animal,'uni',false);
    similar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,similar_animal,'uni',false);
    simdis_pairwise = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise,similar_animal,dissimilar_animal,'uni',false);
    simsim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    dissim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    
    dissimilar_pairwise_early = vertcat(dissimilar_pairwise{:,1});
    similar_pairwise_early = vertcat(similar_pairwise{:,1});
    dissimilar_pairwise_late = vertcat(dissimilar_pairwise{:,2});
    similar_pairwise_late = vertcat(similar_pairwise{:,2});
    simdis_pairwise_early = vertcat(simdis_pairwise{:,1});
    simsim_pairwise_earlylate = vertcat(simsim_pairwise_cross{:});
    dissim_pairwise_earlylate = vertcat(dissim_pairwise_cross{:});
    
    % the pairs of within-comparison
    dissimilar_pairwise_early_mean = nanmean(dissimilar_pairwise_early);
    similar_pairwise_early_mean = nanmean(similar_pairwise_early);
    dissimilar_pairwise_late_mean = nanmean(dissimilar_pairwise_late);
    similar_pairwise_late_mean = nanmean(similar_pairwise_late);
    
    dissimilar_pairwise_early_sem = nanstd(dissimilar_pairwise_early)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_early)));
    similar_pairwise_early_sem = nanstd(similar_pairwise_early)/ ...
        sqrt(sum(~isnan(similar_pairwise_early)));
    dissimilar_pairwise_late_sem = nanstd(dissimilar_pairwise_late)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_late)));
    similar_pairwise_late_sem = nanstd(similar_pairwise_late)/ ...
        sqrt(sum(~isnan(similar_pairwise_late)));
    
    % the pairs of cross-comparison
    simdis_pairwise_early_mean = nanmean(simdis_pairwise_early);
    simsim_pairwise_earlylate_mean = nanmean(simsim_pairwise_earlylate);
    dissim_pairwise_earlylate_mean = nanmean(dissim_pairwise_earlylate);
    
    simdis_pairwise_early_sem = nanstd(simdis_pairwise_early)/ ...
        sqrt(sum(~isnan(simdis_pairwise_early)));
    simsim_pairwise_earlylate_sem = nanstd(simsim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(simsim_pairwise_earlylate)));
    dissim_pairwise_earlylate_sem = nanstd(dissim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(dissim_pairwise_earlylate)));
    
    m_pairwise(curr_subset,:) = [dissimilar_pairwise_early_mean similar_pairwise_early_mean ...
        simdis_pairwise_early_mean dissim_pairwise_earlylate_mean ...
        simsim_pairwise_earlylate_mean similar_pairwise_late_mean];
    s_pairwise(curr_subset,:) = [dissimilar_pairwise_early_sem similar_pairwise_early_sem ...
        simdis_pairwise_early_sem dissim_pairwise_earlylate_sem ...
        simsim_pairwise_earlylate_sem similar_pairwise_late_sem];
    
    % get activity comparison to template means
    
    early_dissimilar_cc_mean = nanmean(total_act_cc_cat_early(early_dissimilar));
    early_similar_cc_mean = nanmean(total_act_cc_cat_early(early_similar));
    late_dissimilar_cc_mean = nanmean(total_act_cc_cat_late(late_dissimilar));
    late_similar_cc_mean = nanmean(total_act_cc_cat_late(late_similar));
    
    early_dissimilar_cc_sem = nanstd(total_act_cc_cat_early(early_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_dissimilar))));
    early_similar_cc_sem = nanstd(total_act_cc_cat_early(early_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_similar))));
    late_dissimilar_cc_sem = nanstd(total_act_cc_cat_late(late_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_dissimilar))));
    late_similar_cc_sem = nanstd(total_act_cc_cat_late(late_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_similar))));
    
    m_template(curr_subset,:) = [early_dissimilar_cc_mean early_similar_cc_mean ...
        late_dissimilar_cc_mean late_similar_cc_mean];
    s_template(curr_subset,:) = [early_dissimilar_cc_sem early_similar_cc_sem ...
        late_dissimilar_cc_sem late_similar_cc_sem];
    

    
    disp(curr_subset);
end

m2 = nanmean(m_template);
s2 = nanmean(s_template);
figure; hold on
bar(m2([1 2 4]),'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2([1 2 4]),s2([1 2 4]),'.k','linewidth',2)
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early dissimilar lever' 'Early similar lever' 'Late similar lever'});
xlabel('Lever correlation')
ylabel('Template activity correlation')


m2 = nanmean(m_pairwise);
s2 = nanmean(s_pairwise);
figure; hold on
bar(m2,'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2,s2,'.k','linewidth',2)
set(gca,'XTick',1:6);
set(gca,'XTickLabel',{'Early dissimilar' 'Early similar' ...
    'Early dissimilar/similar' 'Late/early dissimilar' ...
    'Late/early similar' 'Late similar'});
xlabel('Lever correlation')
ylabel('Pairwise activity correlation')


% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step
subsets = 5;
corr_bins = [-1:0.1:1-0.1 Inf];
m = nan(subsets,length(corr_bins)-1);
s = nan(subsets,length(corr_bins)-1);
for curr_subset = 1:subsets
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_col = cell(8,1); 

    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        %use_sessions = sessions(1:end);
        use_sessions = intersect(8:14,sessions);
        
        curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
            setdiff(1:size(x,2),y)), ...
            template_cluster_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        

        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;

    end
    
    
    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    curr_subset
end
figure;errorbar(corr_bins(1:end-1),nanmean(m,1),nanmean(s,1),'k','linewidth',2);

%% template corrleation: include temporal information TIME-ALIGNED/SHIFTED

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
use_animals = [1 2 4 5 6 7 8];
for curr_animal = use_animals
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x:x+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
     
        % align activity based on lever alignment

        % get lever trials/offsets corresponding to activity trials
        framerate = 28;
        curr_offsets = arrayfun(@(x) ...
            round(align_shift{curr_animal,curr_session}( ...
            lever_trials_used{curr_animal,curr_session} == curr_trials(x)) / ...
            (framerate/100)), ...
            1:length(curr_trials),'uni',false);
        curr_activity_aligned = ...
            arrayfun(@(x) reshape(circshift(reshape(cat_activity(:,x), ...
            move_time_framespan+2*leeway+1,[]), ...
            [curr_offsets{x} 0]),[],1), ...
            1:size(cat_activity,2),'uni',false);
        cat_activity_aligned = horzcat(curr_activity_aligned{:});
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity_aligned,gauss_blur,'same');            
        movement_activity{curr_animal,curr_session} = cat_activity_blur;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);
lever_aligned_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_aligned,lever_trials_used,lever_act_trials,'uni',false);
lever_aligned_cc_crop = cellfun(@(x,y,z) x(ismember(y,z)), ...
    lever_aligned_cc,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);

subsets = 100;
m_template = nan(subsets,4);
s_template = nan(subsets,4);
m_pairwise = nan(subsets,6);
s_pairwise = nan(subsets,6);
for curr_subset = 1:100
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,2);
    total_act_cc = cell(8,2);
    total_act = cell(8,2);
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        
        for earlylate = 1:2
            if earlylate == 1
                use_sessions = sessions(1:5);
            elseif earlylate == 2
                use_sessions = sessions(end-2:end);
            end
            
            curr_lever_cc = cellfun(@(x,y) x(...
                setdiff(1:size(x,2),y)), ...
                lever_aligned_cc_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_lever_cc = horzcat(curr_lever_cc{:});
            
            curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
                activity_corr(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            cat_act_cc = horzcat(curr_act_cc{:});
            
            curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
                movement_activity_crop(curr_animal,use_sessions), ...
                act_template_trials(curr_animal,use_sessions),'uni',false);
            curr_act_cat = horzcat(curr_act{:});
            
            total_act{curr_animal,earlylate} = curr_act_cat;
            total_lever_cc{curr_animal,earlylate} = cat_lever_cc;
            total_act_cc{curr_animal,earlylate} = cat_act_cc;
        end
    end
    
    total_act_pairwise = cellfun(@corrcoef,total_act,'uni',false);
    total_act_pairwise_cross = cell(8,1);
    for i = 1:8
        curr_cat_act = horzcat(total_act{i,:});
        curr_pairwise_cc = corrcoef(curr_cat_act);
        total_act_pairwise_cross{i} = curr_pairwise_cc( ...
            1:size(total_act{i,1},2),size(total_act{i,1},2)+1:end);
    end
    
    total_lever_cc_cat_early = horzcat(total_lever_cc{:,1});
    total_act_cc_cat_early = horzcat(total_act_cc{:,1});
    
    total_lever_cc_cat_late = horzcat(total_lever_cc{:,2});
    total_act_cc_cat_late = horzcat(total_act_cc{:,2});
    
    dissimilar_animal = cellfun(@(x) x > -0.2 & ...
        x < 0.2,total_lever_cc,'uni',false);
    similar_animal = cellfun(@(x) x > 0.6, ...
        total_lever_cc,'uni',false);
  
    early_dissimilar = horzcat(dissimilar_animal{:,1});
    early_similar = horzcat(similar_animal{:,1});
    late_dissimilar = horzcat(dissimilar_animal{:,2});
    late_similar = horzcat(similar_animal{:,2});
    
    % get pairwise activity comparison
    
    dissimilar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,dissimilar_animal,'uni',false);
    similar_pairwise = cellfun(@(x,y) AP_itril(x(y,y),-1), ...
        total_act_pairwise,similar_animal,'uni',false);
    simdis_pairwise = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise,similar_animal,dissimilar_animal,'uni',false);
    simsim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    dissim_pairwise_cross = cellfun(@(x,y,z)  reshape(x(y,z),[],1), ...
        total_act_pairwise_cross,similar_animal(:,1),similar_animal(:,2),'uni',false);
    
    dissimilar_pairwise_early = vertcat(dissimilar_pairwise{:,1});
    similar_pairwise_early = vertcat(similar_pairwise{:,1});
    dissimilar_pairwise_late = vertcat(dissimilar_pairwise{:,2});
    similar_pairwise_late = vertcat(similar_pairwise{:,2});
    simdis_pairwise_early = vertcat(simdis_pairwise{:,1});
    simsim_pairwise_earlylate = vertcat(simsim_pairwise_cross{:});
    dissim_pairwise_earlylate = vertcat(dissim_pairwise_cross{:});
    
    % the pairs of within-comparison
    dissimilar_pairwise_early_mean = nanmean(dissimilar_pairwise_early);
    similar_pairwise_early_mean = nanmean(similar_pairwise_early);
    dissimilar_pairwise_late_mean = nanmean(dissimilar_pairwise_late);
    similar_pairwise_late_mean = nanmean(similar_pairwise_late);
    
    dissimilar_pairwise_early_sem = nanstd(dissimilar_pairwise_early)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_early)));
    similar_pairwise_early_sem = nanstd(similar_pairwise_early)/ ...
        sqrt(sum(~isnan(similar_pairwise_early)));
    dissimilar_pairwise_late_sem = nanstd(dissimilar_pairwise_late)/ ...
        sqrt(sum(~isnan(dissimilar_pairwise_late)));
    similar_pairwise_late_sem = nanstd(similar_pairwise_late)/ ...
        sqrt(sum(~isnan(similar_pairwise_late)));
    
    % the pairs of cross-comparison
    simdis_pairwise_early_mean = nanmean(simdis_pairwise_early);
    simsim_pairwise_earlylate_mean = nanmean(simsim_pairwise_earlylate);
    dissim_pairwise_earlylate_mean = nanmean(dissim_pairwise_earlylate);
    
    simdis_pairwise_early_sem = nanstd(simdis_pairwise_early)/ ...
        sqrt(sum(~isnan(simdis_pairwise_early)));
    simsim_pairwise_earlylate_sem = nanstd(simsim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(simsim_pairwise_earlylate)));
    dissim_pairwise_earlylate_sem = nanstd(dissim_pairwise_earlylate)/ ...
        sqrt(sum(~isnan(dissim_pairwise_earlylate)));
    
    m_pairwise(curr_subset,:) = [dissimilar_pairwise_early_mean similar_pairwise_early_mean ...
        simdis_pairwise_early_mean dissim_pairwise_earlylate_mean ...
        simsim_pairwise_earlylate_mean similar_pairwise_late_mean];
    s_pairwise(curr_subset,:) = [dissimilar_pairwise_early_sem similar_pairwise_early_sem ...
        simdis_pairwise_early_sem dissim_pairwise_earlylate_sem ...
        simsim_pairwise_earlylate_sem similar_pairwise_late_sem];
    
    % get activity comparison to template means
    
    early_dissimilar_cc_mean = nanmean(total_act_cc_cat_early(early_dissimilar));
    early_similar_cc_mean = nanmean(total_act_cc_cat_early(early_similar));
    late_dissimilar_cc_mean = nanmean(total_act_cc_cat_late(late_dissimilar));
    late_similar_cc_mean = nanmean(total_act_cc_cat_late(late_similar));
    
    early_dissimilar_cc_sem = nanstd(total_act_cc_cat_early(early_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_dissimilar))));
    early_similar_cc_sem = nanstd(total_act_cc_cat_early(early_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_early(early_similar))));
    late_dissimilar_cc_sem = nanstd(total_act_cc_cat_late(late_dissimilar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_dissimilar))));
    late_similar_cc_sem = nanstd(total_act_cc_cat_late(late_similar)) ...
        /sqrt(sum(~isnan(total_act_cc_cat_late(late_similar))));
    
    m_template(curr_subset,:) = [early_dissimilar_cc_mean early_similar_cc_mean ...
        late_dissimilar_cc_mean late_similar_cc_mean];
    s_template(curr_subset,:) = [early_dissimilar_cc_sem early_similar_cc_sem ...
        late_dissimilar_cc_sem late_similar_cc_sem];
    

    
    disp(curr_subset);
end

m2 = nanmean(m_template);
s2 = nanmean(s_template);
figure; hold on
bar(m2([1 2 4]),'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2([1 2 4]),s2([1 2 4]),'.k','linewidth',2)
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Early dissimilar lever' 'Early similar lever' 'Late similar lever'});
xlabel('Lever correlation')
ylabel('Template activity correlation')


m2 = nanmean(m_pairwise);
s2 = nanmean(s_pairwise);
figure; hold on
bar(m2,'linewidth',2,'FaceColor','w','BarWidth',0.5);
errorbar(m2,s2,'.k','linewidth',2)
set(gca,'XTick',1:6);
set(gca,'XTickLabel',{'Early dissimilar' 'Early similar' ...
    'Early dissimilar/similar' 'Late/early dissimilar' ...
    'Late/early similar' 'Late similar'});
xlabel('Lever correlation')
ylabel('Pairwise activity correlation')




% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step
subsets = 100;
corr_bins = [-1:0.1:1-0.1 Inf];
m = nan(subsets,length(corr_bins)-1);
s = nan(subsets,length(corr_bins)-1);
for curr_subset = 1:100
    
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            activity_corr{curr_animal,curr_session} = curr_template_corr;
        end
    end
    
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_col = cell(8,1);
    %figure; hold on
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = sessions(end-2:end);
        
        curr_lever_cc = cellfun(@(x,y) x(...
            setdiff(1:size(x,2),y)), ...
            lever_aligned_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        
        col = cellfun(@(x,y) ones(length(x),1)*y,curr_lever_cc, ...
            num2cell(1:length(curr_lever_cc)),'uni',false);
        x = 1:length(use_sessions);
        %scatter(cat_act_cc,cat_lever_cc,10,x);
        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;
        total_col{curr_animal} = vertcat(col{:});
    end
    

    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    curr_subset
end
figure;errorbar(corr_bins(1:end-1),nanmean(m),nanmean(s),'k','linewidth',2);

%% template correlation: include temporal informaiton but align to act template


% try aligning activity: is it the same active cells, just different
% timing?


% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);



% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step
subsets = 5;
corr_bins = [-1:0.1:1-0.1 Inf];
m = nan(subsets,length(corr_bins)-1);
s = nan(subsets,length(corr_bins)-1);
for curr_subset = 1:subsets
    
    % create activity template based on trials used to create lever template
    % use random subset of template trials at end so trials aren't compared to
    % themselves in the next step
    up_clusters = [2 3 NaN 3 1 2 2 2];
    activity_corr = cell(8,15);
    act_template_trials = cell(8,15);
    activity_template_all = cell(8,1);
    for curr_animal = 1:8
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
            template_trials_crop(curr_animal,sessions),'uni',false);
        act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
            round(length(x)/2))),curr_template_trials,'uni',false);
        movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
            movement_activity_crop(curr_animal,sessions), ...
            act_template_trials(curr_animal,sessions),'uni',false);
        activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
        activity_template_all{curr_animal} = activity_template;
        for curr_session = sessions
            cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
            curr_template_corr_grid = corrcoef([activity_template cat_activity]);
            curr_template_corr = curr_template_corr_grid(1,2:end);
            
            curr_activity_split = mat2cell(cat_activity,size(cat_activity,1), ...
                ones(1,size(cat_activity,2)));
            curr_activity_zerocorr = cellfun(@(x) corrcoef( ...
                activity_template(activity_template ~= 0 | x ~= 0), ...
                x(activity_template ~= 0 | x ~= 0)),curr_activity_split,'uni',false);
            curr_template_zerocorr = cellfun(@(x) x(2), curr_activity_zerocorr);
            
            curr_template_corr_dist = sqrt(sum(bsxfun(@minus, ...
                cat_activity,activity_template).^2,1));
            
            activity_corr{curr_animal,curr_session} = curr_template_zerocorr;
            %activity_corr{curr_animal,curr_session} = curr_template_corr;
            %activity_corr{curr_animal,curr_session} = curr_template_corr_dist;
        end
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_col = cell(8,1); 

    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        %use_sessions = sessions(1:end);
        use_sessions = intersect(8:14,sessions);
        
        curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
            setdiff(1:size(x,2),y)), ...
            template_cluster_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        
%         % align activity for each cell/trial to template
%         curr_act = horzcat(movement_activity_crop{curr_animal,use_sessions});
%         curr_frames = move_time_framespan+2*leeway+1;
%         curr_act_reshape = reshape(curr_act,curr_frames,[],size(curr_act,2));
%         
%         curr_template_reshape = reshape(activity_template_all{curr_animal},curr_frames,[]);
%         
%         curr_act_align = nan(size(curr_act_reshape));
%         
%         for curr_cell = 1:size(curr_act_reshape,2);
%             for curr_trial = 1:size(curr_act_reshape,3)
%                 act_xcov = xcov(curr_act_reshape(:,curr_cell,curr_trial), ...
%                     curr_template_reshape(:,curr_cell));
%                 [c ci] = max(act_xcov);
%                 align_shift = -(ci - size(curr_template_reshape,1));
%                 curr_act_align(:,curr_cell,curr_trial) = ...
%                     circshift(curr_act_reshape(:,curr_cell,curr_trial), ...
%                     [align_shift 0]);
%             end
%             curr_cell
%         end
%         
%         cat_act_align = reshape(curr_act_align,[],size(curr_act,2));
%         curr_template_corr_grid = corrcoef( ...
%             [activity_template_all{curr_animal} cat_act_align]);
%         curr_template_corr = curr_template_corr_grid(1,2:end);
%         curr_template_corr_split = mat2cell(curr_template_corr, ...
%             1,cellfun(@(x) size(x,2),movement_activity_crop(curr_animal,use_sessions)));
%         curr_act_align_cc = cellfun(@(x,y) x(setdiff(1:length(x),y)), ...
%             curr_template_corr_split, ...
%             act_template_trials(curr_animal,use_sessions),'uni',false);
%         curr_act_align_cc_cat = horzcat(curr_act_align_cc{:});
        
        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;
        %total_act_cc{curr_animal} = curr_act_align_cc_cat;
    end
    
    
    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    curr_subset
end
figure;errorbar(corr_bins(1:end-1),nanmean(m,1),nanmean(s,1),'k','linewidth',2);


% TK summary of trial/lever activity

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
movement_frames_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x:x+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);
        
        movement_frames_used{curr_animal,curr_session} = ...
            vertcat(curr_rewardmoves_leeway{:});

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
        
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);

movement_frames_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    movement_frames_used,act_trials_used,lever_act_trials,'uni',false);

 


% This was the trial example figure TK asked for

curr_animal = 1;

sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
use_sessions = intersect(8:14,sessions);

% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step
up_clusters = [1 1 NaN 1 2 2 3 2];
activity_corr = cell(8,15);
act_template_trials = cell(8,15);

sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
    template_trials_crop(curr_animal,sessions),'uni',false);
act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
    round(length(x)/2))),curr_template_trials,'uni',false);
movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
    movement_activity_crop(curr_animal,sessions), ...
    act_template_trials(curr_animal,sessions),'uni',false);
activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
for curr_session = sessions
    cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
    curr_template_corr_grid = corrcoef([activity_template cat_activity]);
    curr_template_corr = curr_template_corr_grid(1,2:end);
    activity_corr{curr_animal,curr_session} = curr_template_corr;
end

curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
    setdiff(1:size(x,2),y)), ...
    template_cluster_cc_crop(curr_animal,use_sessions), ...
    act_template_trials(curr_animal,use_sessions),'uni',false);
cat_lever_cc = horzcat(curr_lever_cc{:});

curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
    activity_corr(curr_animal,use_sessions), ...
    act_template_trials(curr_animal,use_sessions),'uni',false);
cat_act_cc = horzcat(curr_act_cc{:});

temp_act_all = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
            movement_activity_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
temp_act = horzcat(temp_act_all{:}) > 0;
temp_lever_all = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
            lever_used_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
temp_lever = vertcat(temp_lever_all{:});

[asdf sort_lever] = sort(cat_lever_cc,'descend');
% pick which cells to use
%use_cells = find(any(vertcat(pyr_class_all{curr_animal,:}),1));
curr_act_cellact = any(any(reshape(temp_act(:,sort_lever(1:35)), ...
    95,[],35,1),3));
use_cells = find(curr_act_cellact);
[c ci] = max(reshape(activity_template,95,[]),[],1);
ci(c == 0) = nan;
[asdf sort_cells] = sort(ci(use_cells));



num_active = permute(sum(any(reshape(temp_act,95,[],size(temp_act,2)) > 0),2),[3 2 1]);
figure;scatter3(cat_lever_cc,cat_act_cc,1:length(cat_lever_cc),num_active*2,'k');
grid off;view(2);
xlabel('Lever template corr')
ylabel('Act template corr')




% plot the top 30 lever/activity correlations
[asdf sort_act] = sort(cat_act_cc,'descend');
sort_act(isnan(cat_act_cc(sort_act))) = [];
act_trials = sort_act(1:35);
figure;
for i = 1:length(act_trials)
    subplot(ceil(sqrt(length(act_trials))),ceil(sqrt(length(act_trials))),i);
    curr_act = reshape(temp_act(:,act_trials(i)),95,[])';
    imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
    ylim([0.5 length(use_cells)+0.5]);
    xlim([0.5 size(curr_act,2)+0.5]);
    axis off
end
subplot(ceil(sqrt(length(act_trials))),ceil(sqrt(length(act_trials))),length(act_trials)+1);
curr_act = reshape(activity_template,95,[])';
imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
ylim([0.5 length(use_cells)+0.5]);
xlim([0.5 size(curr_act,2)+0.5]);
axis off

lever_trials = sort_lever(1:36);
figure;
for i = 1:length(lever_trials)
    subplot(ceil(sqrt(length(lever_trials))),ceil(sqrt(length(lever_trials))),i);
    hold on
    plot(temp_lever(lever_trials(i),:),'r','linewidth',2)
    plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:),'--k','linewidth',2)
    axis off
end

% plot the top 30 lever+activity correlations
lever_trials = sort_lever(1:35);
subplot_side = ceil(sqrt(length(lever_trials)));
figure
for i = 1:length(lever_trials)
    subplot(subplot_side,subplot_side*2,i*2-1);hold on
    plot(temp_lever(lever_trials(i),:),'r','linewidth',2)
    plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:),'--k','linewidth',2)
    axis off

    subplot(subplot_side,subplot_side*2,i*2)
    curr_act = reshape(temp_act(:,lever_trials(i)),95,[])';
    imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
    ylim([0.5 length(use_cells)+0.5]);
    xlim([0.5 size(curr_act,2)+0.5]);
    axis off
end
subplot(subplot_side,subplot_side*2,(length(lever_trials)+1)*2)
curr_act = reshape(activity_template,95,[])';
imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
ylim([0.5 length(use_cells)+0.5]);
xlim([0.5 size(curr_act,2)+0.5]);
ylim([0.5 length(use_cells)+0.5]);
xlim([0.5 size(curr_act,2)+0.5]);
axis off


binary_act = +permute(any(reshape(temp_act,95,[],size(temp_act,2)),1),[2 3 1]);
overlap_cells = binary_act(:,lever_trials)'*binary_act(:,lever_trials);

% sort that in some way, come up with trials with > 4 similar cells

lever_trials = sort_lever(1:35);
lever_trials = lever_trials([4 14 16]);
coactive_cells = find(all(binary_act(:,lever_trials),2));

subplot_side = ceil(sqrt(length(lever_trials)));
figure
for i = 1:length(lever_trials)
    subplot(subplot_side,subplot_side*2,i*2-1);hold on
    plot(temp_lever(lever_trials(i),:),'r','linewidth',2)
    plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:),'--k','linewidth',2)
    axis off

    subplot(subplot_side,subplot_side*2,i*2)
    curr_act = reshape(temp_act(:,lever_trials(i)),95,[])';
    imagesc(curr_act(coactive_cells,:));colormap(gray);caxis([0 0.5]);
    ylim([0.5 length(coactive_cells)+0.5]);
    xlim([0.5 size(curr_act,2)+0.5]);
    axis off
end

% find those frames for trials used
trials_split = mat2cell(1:length(cat_lever_cc),1,cellfun(@length,curr_lever_cc));
trials_split_used = cellfun(@(x) ismember(x,lever_trials),trials_split,'uni',false);

movement_frames_used_crop_nontemplate = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
            movement_frames_used_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        
leveract_plot_frames = cell(15,1);
leveract_plot_frames(use_sessions) = cellfun(@(x,y) x(y,:), ...
    movement_frames_used_crop_nontemplate,trials_split_used,'uni',false);



% trials = [144 116 138 109 90];
% figure
% for i = 1:length(trials)
%     subplot(2,length(trials),i);hold on
%     plot(temp_lever(trials(i),:),'r','linewidth',2)
%     plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:),'--k','linewidth',2)
%     axis off
% 
%     subplot(2,length(trials),i+length(trials));
%     curr_act = reshape(temp_act(:,trials(i)),95,[])';
%     imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
%     ylim([0.5 length(use_cells)+0.5]);
%     xlim([0.5 size(curr_act,2)+0.5]);
%     axis off
% end
% 
% figure;
% max_plot = 10;
% for i = 1:max_plot
%     subplot(2,max_plot,i);hold on
%     plot(temp_lever(sort_lever(i),:),'r','linewidth',2)
%     plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:),'--k','linewidth',2)
%     axis off
% 
%     subplot(2,max_plot,i+max_plot);
%     curr_act = reshape(temp_act(:,sort_lever(i)),95,[])';
%     imagesc(curr_act(use_cells(sort_cells),:));colormap(gray);caxis([0 0.5]);
%     ylim([0.5 length(use_cells)+0.5]);
%     xlim([0.5 size(curr_act,2)+0.5]);
%     axis off
% end







%% k-means lever corr vs activity corr (like above, new)

activity_corr = cell(8,15);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
   
    % Make template: average of last 3 days for animal
    curr_class = vertcat(pyr_class_all{curr_animal,:});
    curr_template = nanmean(curr_class(end-2:end,:));
      
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    
    for curr_session = sessions
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        curr_template_corr_grid = corrcoef([act_template_1{curr_animal} cat_activity]);
        curr_template_corr = curr_template_corr_grid(1,2:end);
        activity_corr{curr_animal,curr_session} = curr_template_corr;

        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
    end
    curr_animal
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);

% plot the lever/activity correlation over days
use_animals = [1 2 4 5 6 7 8];
up_clusters = [2 3 NaN 2 2 3 1 3];
down_clusters = [3 1 NaN 3 2 3 2 2];
cc_thresh = [0.6 0.6 NaN 0.6 0.6 0.7 0.6 0.6];
noncc_thresh_low = -0.1;
noncc_thresh_high = 0.1;
corr_fig = figure;
lever_act_fig = figure;
for curr_animal = use_animals
    curr_template = up_clusters(curr_animal);
    
    template_cluster_cc_crop_cat = horzcat(template_cluster_cc_crop{curr_animal,:});
    trial_num = cumsum(cellfun(@(x) size(x,2),template_cluster_cc_crop(curr_animal,:)));
    figure(corr_fig); subplot(3,3,curr_animal);
    plot(template_cluster_cc_crop_cat(curr_template,:),'.k');
    for i = 1:length(trial_num)
        line([trial_num(i) trial_num(i)],ylim);
    end
    title(['Animal ' num2str(curr_animal) ' lever correlation, template ' ...
        num2str(curr_template)]);
    
    curr_cc_thresh = cc_thresh(curr_animal);
    line(xlim,[curr_cc_thresh curr_cc_thresh],'color','r');
    sessions = cellfun(@(x) ~isempty(x),template_cluster_cc_crop(curr_animal,:));
    cc_thresh_trials = cellfun(@(x) x(curr_template,:) > ...
        curr_cc_thresh,template_cluster_cc_crop(curr_animal,sessions),'uni',false);
    lever_thresh_corr = cellfun(@(x,y) x(curr_template,y), ...
        template_cluster_cc_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    lever_thresh_act_corr = cellfun(@(x,y) AP_itril(corrcoef(+x(:,y)),-1), ...
        movement_activity_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    
    m_lever = cellfun(@nanmean,lever_thresh_corr);
    s_lever = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),lever_thresh_corr);
    m_act = cellfun(@nanmean,lever_thresh_act_corr);
    s_act = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),lever_thresh_act_corr);
    figure(lever_act_fig);subplot(3,3,curr_animal); hold on
    errorbar(m_lever,s_lever,'r','linewidth',2)
    errorbar(m_act,s_act,'k','linewidth',2)
    title(['Animal ' num2str(curr_animal) ' lever correlation, template ' ...
        num2str(curr_template)]);
    legend({'Lever','Activity'});
    ylabel('Mean correlation')
    xlabel('Session')
    
end
% plot the summary data for lever/template corr

template_cluster_cc_crop_mean = cellfun(@(x) nanmean(x,2), ...
    template_cluster_cc_crop,'uni',false);

up_template_cluster_cc = cell(14,1);
for curr_day = 1:14
    curr_animals = intersect(use_animals,find(cellfun(@(x) ~isempty(x), ...
        template_cluster_cc_crop(:,curr_day))));
    curr_template_cc = cellfun(@(x,y) x(y,:), ...
        template_cluster_cc_crop(curr_animals,curr_day), ...
        num2cell(up_clusters(curr_animals))','uni',false);
    up_template_cluster_cc{curr_day} = horzcat(curr_template_cc{:});
end

down_template_cluster_cc = cell(14,1);
for curr_day = 1:14
    curr_animals = intersect(use_animals,find(cellfun(@(x) ~isempty(x), ...
        template_cluster_cc_crop(:,curr_day))));
    curr_template_cc = cellfun(@(x,y) x(y,:), ...
        template_cluster_cc_crop(curr_animals,curr_day), ...
        num2cell(down_clusters(curr_animals))','uni',false);
    down_template_cluster_cc{curr_day} = horzcat(curr_template_cc{:});
end

m_up = cellfun(@nanmean,up_template_cluster_cc);
s_up = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),up_template_cluster_cc);
m_down = cellfun(@nanmean,down_template_cluster_cc);
s_down = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),down_template_cluster_cc);
figure; hold on;
errorbar(m_up,s_up,'k','linewidth',2)
errorbar(m_down,s_down,'r','linewidth',2)
title('Lever template correlation summary')
xlabel('Session')
ylabel('Mean correlation')
legend({'Up templates' 'Down templates'});

% plot summary data for activity corr
lever_thresh_corr = cell(8,15);
lever_thresh_act_corr = cell(8,15);

for curr_animal = use_animals;
    curr_template = up_clusters(curr_animal);
   
    curr_cc_thresh = cc_thresh(curr_animal);
    sessions = cellfun(@(x) ~isempty(x),template_cluster_cc_crop(curr_animal,:));
    cc_thresh_trials = cellfun(@(x) x(curr_template,:) > ...
        curr_cc_thresh,template_cluster_cc_crop(curr_animal,sessions),'uni',false);
    lever_thresh_corr(curr_animal,sessions) = cellfun(@(x,y) x(curr_template,y), ...
        template_cluster_cc_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    lever_thresh_act_corr(curr_animal,sessions) = cellfun(@(x,y) AP_itril(corrcoef(+x(:,y)),-1),...
        movement_activity_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    
end
cat_act_data = cell(14,1);
cat_lever_data = cell(14,1);
for i = 1:14
   cat_act_data{i} = vertcat(lever_thresh_act_corr{:,i}); 
   cat_lever_data{i} = horzcat(lever_thresh_corr{:,i})';
end
m_act = cellfun(@nanmean,cat_act_data);
s_act = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),cat_act_data);
m_lever = cellfun(@nanmean,cat_lever_data);
s_lever = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),cat_lever_data);
figure; hold on;
errorbar(m_act,s_act,'-k','linewidth',2)
errorbar(m_lever,s_lever,'-r','linewidth',2)

lever_thresh_corr = cell(8,15);
lever_thresh_act_corr = cell(8,15);

for curr_animal = use_animals;
    curr_template = up_clusters(curr_animal);
   
    curr_cc_thresh = low_cc_thresh;   
    sessions = cellfun(@(x) ~isempty(x),template_cluster_cc_crop(curr_animal,:));
%     cc_thresh_trials = cellfun(@(x) x(curr_template,:) < ...
%         curr_cc_thresh,template_cluster_cc_crop(curr_animal,sessions),'uni',false);
    cc_thresh_trials = cellfun(@(x) x(curr_template,:) < ...
         noncc_thresh_high & x(curr_template,:) > noncc_thresh_low, ...
         template_cluster_cc_crop(curr_animal,sessions),'uni',false);
    lever_thresh_corr(curr_animal,sessions) = cellfun(@(x,y) x(curr_template,y), ...
        template_cluster_cc_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    lever_thresh_act_corr(curr_animal,sessions) = cellfun(@(x,y) AP_itril(corrcoef(+x(:,y)),-1),...
        movement_activity_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
    
end

% Use random trials instead of down trials for comparison
% for curr_animal = use_animals;
%     curr_template = up_clusters(curr_animal);
%    
%     curr_cc_thresh = cc_thresh;
%     sessions = cellfun(@(x) ~isempty(x),template_cluster_cc_crop(curr_animal,:));
%     num_cc_thresh_trials = cellfun(@(x) sum(x(curr_template,:) > ...
%         curr_cc_thresh),template_cluster_cc_crop(curr_animal,sessions),'uni',false);
%     cc_thresh_trials = cellfun(@(x,y) randsample(size(x,2),y), ...
%         template_cluster_cc_crop(curr_animal,sessions),num_cc_thresh_trials,'uni',false);
%     
%     lever_thresh_corr(curr_animal,sessions) = cellfun(@(x,y) x(curr_template,y), ...
%         template_cluster_cc_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
%     lever_thresh_act_corr(curr_animal,sessions) = cellfun(@(x,y) AP_itril(corrcoef(+x(:,y)),-1),...
%         movement_activity_crop(curr_animal,sessions),cc_thresh_trials,'uni',false);
%     
% end
cat_act_data = cell(14,1);
cat_lever_data = cell(14,1);
for i = 1:14
   cat_act_data{i} = vertcat(lever_thresh_act_corr{:,i}); 
   cat_lever_data{i} = horzcat(lever_thresh_corr{:,i})';
end
m_act = cellfun(@nanmean,cat_act_data);
s_act = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),cat_act_data);
m_lever = cellfun(@nanmean,cat_lever_data);
s_lever = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),cat_lever_data);
errorbar(m_act,s_act,'--k','linewidth',2)
errorbar(m_lever,s_lever,'--r','linewidth',2)

title('Activity with lever template correlation summary')
xlabel('Session')
ylabel('Mean correlation')
legend({'Activity up' 'Lever up' 'Activity down' 'Lever down'});


%% Updated timing preference (based on seconds, not relative)

time_pref = cell(8,15);

for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        movetime_trace = nan(size(movement_trace_all{curr_animal,curr_session}));
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            pyr_activity_all{curr_animal,curr_session}(:,x), ...
            curr_rewardmoves_leeway,'uni',false);
        
        curr_activity_cat = cat(3,curr_activity{:});
        [c ci] = max(curr_activity_cat,[],2);
        curr_onset_timing = permute((c > 0).*ci,[1,3,2]);
        curr_onset_timing(~curr_onset_timing) = NaN;
        %time_pref{curr_animal,curr_session} = nanmean(curr_onset_timing,2);
        
        %time_pref{curr_animal,curr_session} = ...
        %    curr_onset_timing(pyr_class_all{curr_animal,curr_session},:);
        %temp_timing = curr_onset_timing(pyr_class_all{curr_animal,curr_session},:);
        %time_pref{curr_animal,curr_session} = temp_timing(~isnan(temp_timing));  
        
        % to get grid of all onsets       
        curr_onset_act = nan(size(curr_activity_cat));
        [x_grid y_grid z_grid] = meshgrid(1:size(c,2),1:size(c,1),1:size(c,3));
        onset_idx = sub2ind(size(curr_onset_act),y_grid(:),ci(:),z_grid(:));
        curr_onset_act(onset_idx(c(:) ~= 0)) = 1;
        time_pref{curr_animal,curr_session} = ...
            curr_onset_act(pyr_class_all{curr_animal,curr_session},:,:);       
                
    end
    disp(curr_animal);
end

% plot mean 
time_pref_class = cellfun(@(x,y) x(y),time_pref,pyr_class_all,'uni',false);
time_pref_cat = cell(14,1);
use_animals = [1 2 4 5 6 7 8];
for i = 1:14
    time_pref_cat{i} = vertcat(time_pref_class{use_animals,i});    
end
m = cellfun(@nanmean,time_pref_cat)./28;
s = cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),time_pref_cat)./28;
figure;errorbar(m,s,'k','linewidth',2);

% plot cumulative sum distribution
figure; hold on
day_combine = {[1 2] [3 4] [5 6] [7 8] [9 10] [11 12] [13 14]};
col = jet(length(day_combine));
for i = 1:length(day_combine)
    curr_data = vertcat(time_pref_cat{day_combine{i}});
    plot(sort(curr_data)./28,linspace(0,1,length(curr_data)),'color',col(i,:));
end
ylabel('Fraction of cells')
xlabel('Average timing (s)')
title('Cumulative distribution of activity onsets')
legend(cellfun(@num2str,day_combine,'uni',false))


% find changes within timing from when it first appeared
pref_change = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_pref = horzcat(time_pref{curr_animal,:});
    curr_class = vertcat(pyr_class_all{curr_animal,:})';
    [c ci] = max(curr_class,[],2);
    idx = sub2ind(size(curr_class),1:size(curr_class,1),ci');
    
    curr_pref_change = bsxfun(@minus,curr_pref,curr_pref(idx)');
    
    pref_change(curr_animal,sessions) = mat2cell(curr_pref_change, ...
        size(curr_pref_change,1),ones(1,size(curr_pref_change,2)));
end
use_animals = [1 2 4 5 6 7 8];
cat_pref_change = cell(1,14);
for i = 1:14
   cat_pref_change{i} = vertcat(pref_change{use_animals,i});
end


% using activity onsets, find distribution
use_animals = [1 2 4 5 6 7 8];
day_combine = num2cell(1:14);
cat_act_onset = cell(length(day_combine),1);
col = [linspace(0,1,length(day_combine))' zeros(length(day_combine),1) ...
    zeros(length(day_combine),1)];
figure; hold on;
for i = 1:length(day_combine)
    curr_data = cellfun(@(x) permute(x,[2 3 1]), ...
        time_pref(use_animals,day_combine{i}),'uni',false);
    curr_data_reshape = cellfun(@(x) reshape(x,size(x,1),[]),curr_data,'uni',false);    
    cat_act_onset{i} = horzcat(curr_data_reshape{:}); 
    
    total_onset = nansum(cat_act_onset{i}');
    plot(linspace(0,length(total_onset)/28,length(total_onset))-5/28, ...
        cumsum(total_onset)/sum(total_onset),'color',col(i,:),'linewidth',2)
end
xlabel('Time (s)')
ylabel('Fraction of activity onsets')
legend({'Sessions 1-3' 'Sessions 4-6' 'Sessions 7-9' 'Sessions 10-12' 'Sessions 13-14'})

%% Make figure for act/lever corr with mean examples at bottom

% get activity with lever 
activity_template = cell(8,1);
movement_activity = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        curr_activity = cellfun(@(x) ...
            reshape(pyr_activity_all{curr_animal,curr_session}(:,x)',[],1), ...
            curr_rewardmoves_leeway,'uni',false);

        cat_activity = horzcat(curr_activity{:});
        % sometimes have nans, make those zero
        cat_activity(isnan(cat_activity)) = 0;
        
        % to blur activity
        gauss_blur = fspecial('gaussian',[100,1],10);
        cat_activity_blur = conv2(cat_activity,gauss_blur,'same');
            
        movement_activity{curr_animal,curr_session} = cat_activity;
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
    end
    disp(curr_animal)
    
end

% for up template: look at trials with high lever correlation for activity
lever_act_trials = cellfun(@(x,y) intersect(x,y), ...
    act_trials_used,lever_trials_used,'uni',false);
% crop activity/lever to have same trials
movement_activity_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    movement_activity,act_trials_used,lever_act_trials,'uni',false);
template_cluster_cc_crop = cellfun(@(x,y,z) x(:,ismember(y,z)), ...
    template_cluster_cc,lever_trials_used,lever_act_trials,'uni',false);
lever_used_crop = cellfun(@(x,y,z) x(ismember(y,z),:), ...
    lever_used,lever_trials_used,lever_act_trials,'uni',false);

crop_lever_trials = cellfun(@(x,y) ismember(x,y), ...
    lever_trials_used,lever_act_trials,'uni',false);
template_trials_crop = cellfun(@(x,y) x(y(1:length(x))), ...
    template_trials,crop_lever_trials,'uni',false);


% create activity template based on trials used to create lever template
% use random subset of template trials at end so trials aren't compared to
% themselves in the next step
up_clusters = [1 1 NaN 1 2 2 3 2];
activity_corr = cell(8,15);
act_template_trials = cell(8,15);
activity_template_all = cell(8,1);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    curr_template_trials = cellfun(@(x) find(x == up_clusters(curr_animal)), ...
        template_trials_crop(curr_animal,sessions),'uni',false);
    act_template_trials(curr_animal,sessions) = cellfun(@(x) x(randsample(length(x), ...
        round(length(x)/2))),curr_template_trials,'uni',false);
    movement_activity_lever_template = cellfun(@(x,y) x(:,y), ...
        movement_activity_crop(curr_animal,sessions), ...
        act_template_trials(curr_animal,sessions),'uni',false);
    activity_template = nanmean(horzcat(movement_activity_lever_template{:}),2);
    activity_template_all{curr_animal} = activity_template;
    for curr_session = sessions
        cat_activity = horzcat(movement_activity_crop{curr_animal,curr_session});
        curr_template_corr_grid = corrcoef([activity_template cat_activity]);
        curr_template_corr = curr_template_corr_grid(1,2:end);
        activity_corr{curr_animal,curr_session} = curr_template_corr;
    end
end

for earlylate = 1:2
    
    if earlylate == 1
        curr_sessions = 1:3;
    elseif earlylate == 2
        curr_sessions = 8:14;
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,1);
    total_act_cc = cell(8,1);
    total_act = cell(8,1);
    
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        %use_sessions = sessions(1:end);
        use_sessions = intersect(curr_sessions,sessions);
        
        curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
            setdiff(1:size(x,2),y)), ...
            template_cluster_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        
        curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
            movement_activity_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        curr_act_cat = horzcat(curr_act{:});
        
        
        total_lever_cc{curr_animal} = cat_lever_cc;
        total_act_cc{curr_animal} = cat_act_cc;
        total_act{curr_animal} = curr_act_cat;
        
    end
    
    
    total_lever_cc_cat = horzcat(total_lever_cc{:});
    total_act_cc_cat = horzcat(total_act_cc{:});
    
    corr_bins = [-1:0.1:1-0.1 Inf];
    [n bins] = histc(total_lever_cc_cat,corr_bins);
    m(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins);
    s(curr_subset,unique(bins)) = grpstats(total_act_cc_cat,bins,'sem');
    
    curr_subset
    
    figure;errorbar(corr_bins(1:end-1),nanmean(m,1),nanmean(s,1),'k','linewidth',2);
    xlabel('Lever template correlation')
    ylabel('Activity template correlation')
    title(['Sessions ' num2str(curr_sessions(1)) ':' ...
        num2str(curr_sessions(end))]);
    
    % plot the average activity of cells in bins
    
    %act_corr_bins = [-1:0.25:1-0.25 Inf];
    %act_corr_bins = [-1 -0.5 -0.25 0 0.25 0.5 Inf];
    act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
    
    [n act_bins] = histc(total_lever_cc_cat,act_corr_bins);
    
    bin_split = mat2cell(act_bins,1,cellfun(@length,total_act_cc));
    act_grp_all = cell(8,1);
    for curr_animal = use_animals
        act_grp = nan(size(total_act{curr_animal},1),length(act_corr_bins));
        act_grp(:,unique(bin_split{curr_animal})) = ...
            grpstats([total_act{curr_animal}]',bin_split{curr_animal})';
        act_grp_reshape = reshape(act_grp,move_time_framespan+2*leeway+1, ...
            [],length(act_corr_bins));
        act_grp_all{curr_animal} = act_grp_reshape;
    end
    
    act_grp_cat = cat(2,act_grp_all{:});
    
    activity_template_cat = vertcat(activity_template_all{use_animals});
    activity_template_reshape = reshape(activity_template_cat, ...
        move_time_framespan+2*leeway+1,[])';
    [c ci] = max(activity_template_reshape,[],2);
    use_cells = find(c > 0.01);
    [asdf sort_idx] = sort(ci(use_cells));
    
    activity_template_reshape_norm = bsxfun(@times,activity_template_reshape,1./c);
    
    figure;
    for i = 1:length(act_corr_bins)-1
        subplot(1,length(act_corr_bins),i);
        
        imagesc(act_grp_cat(:,use_cells(sort_idx),i)');
        
        %curr_act_grp_cat_norm = bsxfun(@times,act_grp_cat(:,:,i)',1./c);
        %imagesc(curr_act_grp_cat_norm(use_cells(sort_idx),:));
        
        colormap(hot);
        caxis([0 0.2]);
        axis off
        title([num2str(act_corr_bins(i)) ...
            ':' num2str(act_corr_bins(i+1))]);
    end
    subplot(1,length(act_corr_bins),length(act_corr_bins));
    imagesc(activity_template_reshape(use_cells(sort_idx),:));
    caxis([0 0.2]);
    colormap(hot);
    axis off
    title('Template')
    colormap(hot);
    
    set(gcf,'Name',['Sessions ' num2str(curr_sessions(1)) ':' ...
        num2str(curr_sessions(end))]);
end

% plot the average lever trace for each used animal in each bin

for earlylate = 1:2
    
    if earlylate == 1
        curr_sessions = 1:3;
        figure;
        col = [0 0 1];
    elseif earlylate == 2
        curr_sessions = 8:14;
        col = [1 0 0];
    end
    
    use_animals = [1 2 4 5 6 7 8];
    total_lever_cc = cell(8,1);
    
    for curr_animal = use_animals;
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        %use_sessions = sessions(1:end);
        use_sessions = intersect(curr_sessions,sessions);
        
        curr_lever_cc = cellfun(@(x,y) x(up_clusters(curr_animal), ...
            setdiff(1:size(x,2),y)), ...
            template_cluster_cc_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_lever_cc = horzcat(curr_lever_cc{:});
        
        curr_act_cc = cellfun(@(x,y) x(setdiff(1:size(x,2),y)), ...
            activity_corr(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        cat_act_cc = horzcat(curr_act_cc{:});
        
        curr_act = cellfun(@(x,y) x(:,setdiff(1:size(x,2),y)), ...
            movement_activity_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        curr_act_cat = horzcat(curr_act{:});
        
        total_lever_cc{curr_animal} = cat_lever_cc;
    end
    total_lever_cc_cat = horzcat(total_lever_cc{:});
       
    %act_corr_bins = [-1:0.25:1-0.25 Inf];
    %act_corr_bins = [-1 -0.5 -0.25 0 0.25 0.5 Inf];
    act_corr_bins = [-1 -0.5 -0.1667 0.1667 0.5 Inf];
    
    [n act_bins] = histc(total_lever_cc_cat,act_corr_bins);
    
    bin_split = mat2cell(act_bins,1,cellfun(@length,total_lever_cc));
    
    lever_used_bin = nan(length(act_corr_bins)-1,302,8);
    lever_used_bin_sem = nan(length(act_corr_bins)-1,302,8);
    for curr_animal = use_animals
        
        sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
        use_sessions = intersect(curr_sessions,sessions);
        
        lever_used_crop_nontemplate = cellfun(@(x,y) x(setdiff(1:size(x,1),y),:), ...
            lever_used_crop(curr_animal,use_sessions), ...
            act_template_trials(curr_animal,use_sessions),'uni',false);
        lever_cat = vertcat(lever_used_crop_nontemplate{:});
        
        lever_used_bin(unique(bin_split{curr_animal}),:,curr_animal) = ...
            grpstats(lever_cat,bin_split{curr_animal});
        lever_used_bin_sem(unique(bin_split{curr_animal}),:,curr_animal) = ...
            grpstats(lever_cat,bin_split{curr_animal},'sem');
        
    end
    
    for curr_animal = use_animals
        for i = 1:length(act_corr_bins)-1
            subplot(length(act_corr_bins),8,8*(i-1)+curr_animal);hold on;
            plot(lever_used_bin(i,:,curr_animal),'color',col,'linewidth',2);
            
            if earlylate == 1
                plot(lever_k_templates{curr_animal}(up_clusters(curr_animal),:), ...
                    '--k','linewidth',2)
            end
            
            
            if i == 1
                title(curr_animal)
            end
            
            if curr_animal == 1
                ylabel(num2str(act_corr_bins(i)));
            end
            
            %    plot(lever_used_bin(i,:,curr_animal) + ...
            %        lever_used_bin_sem(i,:,curr_animal),'r','linewidth',2);
            %    plot(lever_used_bin(i,:,curr_animal) - ...
            %        lever_used_bin_sem(i,:,curr_animal),'r','linewidth',2);
            axis off
        end
    end
end







%% updated pyr:gad balance by movement

% separate movement traces into movement epochs
move_starts = cellfun(@(x) find(diff([0 x 0]) == 1),movement_trace_all,'uni',false);
move_stops = cellfun(@(x) find(diff([0 x 0]) == -1),movement_trace_all,'uni',false);

move_epoch_frames = cellfun(@(x,y) ...
    arrayfun(@(z) x(z):y(z),1:length(x),'uni',false), ...
    move_starts,move_stops,'uni',false);

% gad_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     gad_activity_thresh_all,move_epoch_frames,gad_class_all,'uni',false);
% pyr_move_active = cellfun(@(x,y,w) cell2mat(cellfun(@(z) any(x(w,z),2),y,'uni',false)), ...
%     pyr_activity_all,move_epoch_frames,pyr_class_all,'uni',false);
gad_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   gad_activity_thresh_all,move_epoch_frames,'uni',false);
pyr_move_active = cellfun(@(x,y) cell2mat(cellfun(@(z) any(x(:,z),2),y,'uni',false)), ...
   pyr_activity_all,move_epoch_frames,'uni',false);

use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);
gc = horzcat(g{use_animals,:});

bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
figure;
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement-based activity')


% get pyr:gad ratio for restricted movement lengths
movement_lengths = cellfun(@(x) cellfun(@length,x),move_epoch_frames,'uni',false);
movement_lengths_cat = horzcat(movement_lengths{use_animals,:})/28;

movement_bins = [0:0.5:5.5];
col = jet(length(movement_bins)-1);

figure; hold on
for i = 1:length(movement_bins)-1      
    curr_bin_idx = movement_lengths_cat > movement_bins(i) & ...
        movement_lengths_cat <= movement_bins(i+1);
    
    curr_move = movement_lengths_cat(curr_bin_idx);
    curr_pc = pc(curr_bin_idx);
    curr_gc = gc(curr_bin_idx);
    
    bin_edges = [0:0.02:0.25];
    %bin_edges = [0:0.1:1-0.1 Inf];
    [n bins] = histc(curr_pc,bin_edges);
    gc_mean = grpstats(curr_gc,bins);
    gc_sem = grpstats(curr_gc,bins,'sem');
    use_bins = unique(bins) ~= 0;
    use_bin_edges = bin_edges(unique(bins(bins ~= 0)));
    errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins), ...
        'color',col(i,:),'linewidth',2);  
end
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title(['Movement-based activity, time ' ...
    num2str(movement_bins(i)) ':' num2str(movement_bins(i+1))])
ylim([0 0.8])
xlim([0 0.25])



%%%%%%% Do above but with movements being 3s and only rewarded movements



% get activity with lever 
pyr_move_active = cell(8,15);
gad_move_active = cell(8,15);
movement_lengths = cell(8,15);
act_trials_used = cell(8,15);
for curr_animal = 1:8
    sessions = find(cellfun(@(x) ~isempty(x), pyr_class_all(curr_animal,:)));
    for curr_session = sessions
        
        % look at 3 seconds following movement onset (to go with lever)
        move_time_framespan = 28*3;
        %move_time_frames = cellfun(@(x) x(1):x(1)+move_time_framespan, ...
        %    move_epoch_frames{curr_animal,curr_session},'uni',false);
        move_time_frames = cellfun(@(x) x, ...
            move_epoch_frames{curr_animal,curr_session},'uni',false);
        
        curr_rewardmoves = move_epoch_frames{curr_animal,curr_session}( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardmoves_timed = move_time_frames( ...
            rewarded_movements{curr_animal,curr_session});
        curr_rewardframes = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) - x(1) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_rewardframes = cellfun(@(x) x(1),curr_rewardframes,'uni',false);
        % get corresponding trials
        curr_rewardframes_abs = cellfun(@(x) intersect(x, ...
            reward_frames_all{curr_animal,curr_session}) + 1, ...
            curr_rewardmoves,'uni',false);
        curr_trials = cellfun(@(x) trial_frames_all{curr_animal,curr_session}(x(1)), ...
            curr_rewardframes_abs);
        
        % add leeway to movements (fixed frames)
        leeway = 5;
        leeway_frames = cellfun(@(x,y) [leeway ...
            leeway],curr_rewardmoves_timed,curr_rewardframes,'uni',false);
        curr_rewardmoves_leeway = cellfun(@(x,y) (x(1)-y(1)):(x(end)+y(2)), ...
            curr_rewardmoves_timed,leeway_frames,'uni',false);
        eliminate_movements = cellfun(@(x) ...
            any(x < 1 | x > length(movement_trace_all ...
            {curr_animal,curr_session})),curr_rewardmoves_leeway);
        curr_rewardmoves_leeway(eliminate_movements) = [];
        curr_trials(eliminate_movements) = [];
        
        [curr_pyr_activity curr_time] = cellfun(@(x) ...
            max(pyr_activity_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);
        
        [curr_gad_activity curr_time] = cellfun(@(x) ...
            max(gad_activity_thresh_all{curr_animal,curr_session}(:,x)>0,[],2), ...
            curr_rewardmoves_leeway,'uni',false);

        pyr_move_active{curr_animal,curr_session} = horzcat(curr_pyr_activity{:});
        gad_move_active{curr_animal,curr_session} = horzcat(curr_gad_activity{:});
        act_trials_used{curr_animal,curr_session} = curr_trials;
        
        movement_lengths{curr_animal,curr_session} = ...
            cellfun(@length,curr_rewardmoves_leeway);
        
    end

    curr_animal
    
end


% get pyr:gad ratio for movements

use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);
gc = horzcat(g{use_animals,:});


bin_edges = [0:0.02:0.25];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
figure;
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement-based activity')


% get pyr:gad ratio for restricted movement lengths
use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x,1),pyr_move_active(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
g(use_animaldays) = cellfun(@(x) nanmean(x,1),gad_move_active(use_animaldays),'uni',false);
gc = horzcat(g{use_animals,:});

movement_lengths_cat = horzcat(movement_lengths{use_animals,:})/28;

movement_bins = [0:0.5:5.5];
col = jet(length(movement_bins)-1);

figure; hold on
for i = 1:length(movement_bins)-1      
    curr_bin_idx = movement_lengths_cat > movement_bins(i) & ...
        movement_lengths_cat <= movement_bins(i+1);
    
    curr_move = movement_lengths_cat(curr_bin_idx);
    curr_pc = pc(curr_bin_idx);
    curr_gc = gc(curr_bin_idx);
    
    bin_edges = [0:0.02:0.85];
    %bin_edges = [0:0.1:1-0.1 Inf];
    [n bins] = histc(curr_pc,bin_edges);
    gc_mean = grpstats(curr_gc,bins);
    gc_sem = grpstats(curr_gc,bins,'sem');
    use_bins = unique(bins) ~= 0;
    use_bin_edges = bin_edges(unique(bins(bins ~= 0)));
    if sum(curr_bin_idx) == 0
        continue
    end
    %errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins), ...
    %    'color',col(i,:),'linewidth',2);  
    
    plot(curr_pc+0.005*randn(size(curr_move)),curr_gc+0.005*rand(size(curr_move)),'.','color',col(i,:));
end
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title(['Movement-based activity, time ' ...
    num2str(movement_bins(i)) ':' num2str(movement_bins(i+1))])
ylim([0 0.8])
xlim([0 0.25])



% Get continuous pyr-gad ratio

use_animals = [1 2 4 5 6 7 8];
p = cell(size(movement_trace_all));
g = cell(size(movement_trace_all));
use_animaldays = cellfun(@(x) ~isempty(x),movement_trace_all);
p(use_animaldays) = cellfun(@(x) nanmean(x>0,1), ...
    pyr_activity_all(use_animaldays),'uni',false);
pc = horzcat(p{use_animals,:});
g(use_animaldays) = cellfun(@(x) nanmean(x>0,1), ...
    gad_activity_thresh_all(use_animaldays),'uni',false);
gc = horzcat(g{use_animals,:});

bin_edges = [0:0.005:0.08];
%bin_edges = [0:0.1:1-0.1 Inf];
[n bins] = histc(pc,bin_edges);
figure;
gc_mean = grpstats(gc,bins);
gc_sem = grpstats(gc,bins,'sem');
use_bins = unique(bins) ~= 0;
use_bin_edges = bin_edges(unique(bins(bins ~= 0))+1);
errorbar(use_bin_edges,gc_mean(use_bins),gc_sem(use_bins),'k','linewidth',2);
xlabel('Fraction of pyramidal move cells')
ylabel('Fraction of gad move cells')
title('Movement-based activity')




