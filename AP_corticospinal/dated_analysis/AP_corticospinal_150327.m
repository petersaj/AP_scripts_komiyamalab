%% Correlation between ROIs (identify same cells, or co-active groups)

d_grp_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    dendrite_corr = [];
    
    analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
    analysis_files = dir([analysis_path filesep '*.mat']);
    num_sessions = length(analysis_files);
    
    for curr_session = 1:num_sessions;
        
        curr_file = [analysis_path filesep analysis_files(curr_session).name];
        curr_data = load(curr_file);
        
        % Smooth data
        roi_trace_df_smooth = nan(size(curr_data.im.roi_trace_df));
        for curr_roi = 1:size(roi_trace_df_smooth,1);
            roi_trace_df_smooth(curr_roi,:) = smooth(curr_data.im.roi_trace_df(curr_roi,:),30,'loess');
        end
        
        curr_corr = corrcoef(roi_trace_df_smooth');
        dendrite_corr(:,:,curr_session) = curr_corr;
        
    end
    
    median_dendrite_corr = tril(median(dendrite_corr,3),-1);
    
    [d1,d2] = find(median_dendrite_corr >= 0.8);
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
    ungrouped = setdiff([1:size(dendrite_corr,1)]',unique(vertcat(d_grp{:})));
    d_grp = [d_grp;num2cell(ungrouped)];
    
    d_grp_all{curr_animal} = d_grp;
    
    disp(['Grouped animal ' animal]);
    
end

%% Apply dendrite grouping to move-onset aligned analysis
% requires analysis structure from AP_corticospinal_prepare_processed

for curr_animal = 1:length(analysis)
   for curr_session = 1:length(analysis(curr_animal).im)
       
       if isempty(analysis(curr_animal).im(curr_session).move_onset_aligned)
           continue
       end
       
        curr_grp_act = cell2mat(cellfun(@(x) ...
            nanmean(analysis(curr_animal).im(curr_session).move_onset_aligned(:,:,x),3), ...
            permute(d_grp_all{curr_animal},[2 3 1]),'uni',false));
    
        analysis(curr_animal).im(curr_session).move_onset_aligned = ...
            curr_grp_act;
        
   end
end


%% Get activity-triggered lever press

animal = 'AP123';

analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
analysis_files = dir([analysis_path filesep '*.mat']);
num_sessions = length(analysis_files);

avg_activity_lever = cell(num_sessions,1);
active_dendrite = cell(num_sessions,1);

for curr_session = 1:num_sessions
    
    curr_file = [analysis_path filesep analysis_files(curr_session).name];
    curr_data = load(curr_file);
    
    % Threshold data
    event_thresh = 3;
    curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
    n_rois = size(curr_df_thresh,1);
    
    % Make matricies for pulling out frames
    surrounding_time = [-300 300]; % (ms)
    surrounding_samples = (surrounding_time/1000)*curr_data.bhv.xsg_sample_rate;
    
    avg_activity_lever{curr_session} = ...
        nan(size(curr_df_thresh,1),surrounding_samples(2)-surrounding_samples(1)+1);
    
    active_dendrite{curr_session} = nan(n_rois,1);
    
    for curr_roi = 1:n_rois;
        
        curr_event_onset_frames = find(diff(curr_df_thresh(curr_roi,:) > 0) == 1);
        
        curr_event_onset_samples = ...
            round(curr_data.bhv.frame_times(curr_event_onset_frames)*curr_data.bhv.xsg_sample_rate);
        
        pull_samples = cell2mat(arrayfun(@(x) curr_event_onset_samples(x) + ...
            surrounding_samples(1) : curr_event_onset_samples(x) + ...
            surrounding_samples(2),[1:length(curr_event_onset_samples)]','uni',false));
        
        curr_activity_lever = curr_data.bhv.lever_force(pull_samples);
        
        if ~isempty(curr_activity_lever)
            avg_activity_lever{curr_session}(curr_roi,:) = nanmean(curr_activity_lever,1);
        end
        
        active_dendrite{curr_session}(curr_roi) = length(curr_event_onset_frames)/size(curr_df_thresh,2);
        
    end
    
    disp(['Session ' num2str(curr_session)]);
end


%% Identify avg activity around movements and correlation to movement

animal = 'AP116';

analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
%analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_roi_fixed_bg'];
analysis_files = dir([analysis_path filesep '*.mat']);
n_sessions = length(analysis_files);

move_act_corr = cell(n_sessions,1);
shuff_move_act_corr = cell(n_sessions,1);

for curr_session = 1:n_sessions
    
    curr_file = [analysis_path filesep analysis_files(curr_session).name];
    curr_data = load(curr_file);
    
    event_thresh = 3;
    curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
    
    [lever_active,lever_force_resample] = ...
        AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
    
    %     % Get activity surrounding movements
    %
    %     surround_frames = [-90 90];
    %
    %
    %     curr_move_onsets_samples = find(diff(lever_active) == 1);
    %     [~,curr_move_onsets_frames] = ...
    %         arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
    %         curr_move_onsets_samples(x))), 1:length(curr_move_onsets_samples));
    %
    %     pull_move_onset_frames = cell2mat(arrayfun(@(x) curr_move_onsets_frames(x) + ...
    %         surround_frames(1) : curr_move_onsets_frames(x) + ...
    %         surround_frames(2),[1:length(curr_move_onsets_frames)]','uni',false));
    %     pull_move_onset_frames(any(pull_move_onset_frames < 1 | ...
    %         pull_move_onset_frames > size(curr_data.im.roi_trace_df,2),2),:) = [];
    %
    %     move_onset_activity = ...
    %         reshape(curr_data.im.roi_trace_df(:,pull_move_onset_frames')', ...
    %         [],size(pull_move_onset_frames,1),size(curr_data.im.roi_trace_df,1));
    %
    %     curr_move_offsets_samples = find(diff(lever_active) == -1);
    %     [~,curr_move_offsets_frames] = ...
    %         arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
    %         curr_move_onsets_samples(x))), 1:length(curr_move_offsets_samples));
    %
    %     pull_move_offset_frames = cell2mat(arrayfun(@(x) curr_move_offsets_frames(x) + ...
    %         surround_frames(1) : curr_move_offsets_frames(x) + ...
    %         surround_frames(2),[1:length(curr_move_offsets_frames)]','uni',false));
    %     pull_move_offset_frames(any(pull_move_offset_frames < 1 | ...
    %         pull_move_offset_frames > size(curr_data.im.roi_trace_df,2),2),:) = [];
    %
    %     move_offset_activity = ...
    %         reshape(curr_data.im.roi_trace_df(:,pull_move_offset_frames')', ...
    %         [],size(pull_move_offset_frames,1),size(curr_data.im.roi_trace_df,1));
    %
    %     %%%% CURRENTLY HERE
    %     % looking at how activity relates to movements: looked not correlated,
    %     % or even anti-correlated at first glance?
    %     % could also just look at move onset aligned from analysis...
    
    
    % Analyze general activity-movement correlation
    n_frames = size(curr_data.im.roi_trace_df,2);
    movement_frames = lever_active(round(curr_data.bhv.frame_times(1:n_frames)*1000))';
    
    % Shuffle movement trace (currently move/still epochs all intact)
    move_epochs = [1 diff([movement_frames]) 1];
    move_epochs_length = diff(find(move_epochs ~= 0));
    movement_frames_epochs = mat2cell(movement_frames',move_epochs_length);
    
    num_rep = 1000;
    perms = shake(repmat(transpose(1:length(movement_frames_epochs)),...
        1,num_rep),1);
    
    movement_frames_epochs_shuffle = nan(length(movement_frames),num_rep);
    for i = 1:num_rep
        movement_frames_epochs_shuffle(:,i) = ...
            vertcat(movement_frames_epochs{perms(:,i)});
    end
    
    % Get correlation of each cell with movement / shuffled movement
    move_act_corr_grid = ...
        corrcoef([movement_frames' movement_frames_epochs_shuffle curr_df_thresh']);
    
    curr_move_act_corr = move_act_corr_grid(num_rep+2:end,1);
    curr_shuff_move_act_corr = move_act_corr_grid(num_rep+2:end,2:1+num_rep);
    
    move_act_corr{curr_session} = curr_move_act_corr;
    shuff_move_act_corr{curr_session} = prctile(sort(curr_shuff_move_act_corr),[2.5 97.5],2);
    
    disp(['Session ' num2str(curr_session)]);
end

col = jet(n_sessions);
figure; hold on;
for curr_session = 1:n_sessions
    use_cells = ~isnan(move_act_corr{curr_session});
    plot(sort(move_act_corr{curr_session}(use_cells)),(1:sum(use_cells))/sum(use_cells), ...
        'color',col(curr_session,:));
end

figure;
for curr_session = 1:n_sessions
    subplot(4,4,curr_session); hold on;
    use_cells = ~isnan(move_act_corr{curr_session});
    plot(shuff_move_act_corr{curr_session},'--r');
    plot(sort(move_act_corr{curr_session}(use_cells)),'k');
end

figure;hold on;
for curr_session = 1:n_sessions
    plot(move_act_corr{curr_session}(sort_idx),'.','color',col(curr_session,:));
end



%% Look for "blocks" of active/coactive dendrites (didn't see anything?)

animal = 'AP123';

analysis_path = ['/usr/local/lab/People/Andy/Data' filesep animal filesep animal '_batch_thresh_roi'];
analysis_files = dir([analysis_path filesep '*.mat']);
n_sessions = length(analysis_files);

for curr_session = 1:n_sessions
    
    curr_file = [analysis_path filesep analysis_files(curr_session).name];
    curr_data = load(curr_file);
    
    event_thresh = 3;
    curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
    
    df_thresh_binary = curr_df_thresh > 0;
    
    activity_spread = 30;
    activity_conv = ones(1,activity_spread);
    
    df_thresh_spread = conv2(+df_thresh_binary,activity_conv,'same') > 0;
    




end



%% Slide 3) Correlation of population activity over days

% Population correlation (of average temporal activity) over days
pop_corr = nan(14,14,length(analysis));
for curr_animal = 1:length(analysis)
    curr_act = cell(1,length(analysis(curr_animal).im));
    for j = 1:length(analysis(curr_animal).im)
        curr_session_act = permute(nanmean(analysis(curr_animal).im(j).move_onset_aligned,1),[2,3,1]);
        
        if ~isempty(curr_session_act)
            
            % group dendrites (by avg at the moment)
            %curr_session_act_grp = cell2mat(cellfun(@(x) nanmean(curr_session_act(:,x),2), ...
            %    d_grp_all{curr_animal}','uni',false));
            
            curr_act{j} = reshape(curr_session_act(30:90,:),[],1);
            
        end
    end
    use_sessions = find(cellfun(@(x) ~isempty(x),curr_act));
    pop_corr(use_sessions,use_sessions,curr_animal) = ...
        corrcoef(horzcat(curr_act{:}),'rows','complete');
end


% % Average population activity over days
% avg_act = nan(length(analysis),14);
% for curr_animal = 1:length(analysis)
%     for j = 1:length(analysis(curr_animal).im)
%         curr_session_act = permute(nanmean(analysis(curr_animal).im(j).move_onset_aligned,1),[2,3,1]);
%         
%         if ~isempty(curr_session_act)
%             
%             % group dendrites (by avg at the moment)
%             %curr_session_act_grp = cell2mat(cellfun(@(x) nanmean(curr_session_act(:,x),2), ...
%             %    d_grp_all{curr_animal}','uni',false));
%             
%             avg_act(curr_animal,j) = nanmean(curr_session_act(:));
%             
%         end
%     end
% end




%% Slide 6) Make average movie surrounding lever press
animals = {'AP140' 'AP141' 'AP142' 'AP145' 'AP146' 'AP147' 'AP148' 'AP150'};

for curr_animal = 1:length(animals);
    
animal = animals{curr_animal};

%curr_session = 10;
for curr_session = 7;
    
    clearvars -except animal curr_session animals curr_animal
    
    
    data_dir = '/usr/local/lab/People/Andy/Data';
    curr_analysis_dir = [data_dir filesep animal filesep ...
        animal '_batch_thresh_roi'];
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
    curr_data = load(curr_file);
    
  
    num_trials = length(curr_data.bhv.bhv_times);
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),curr_data.bhv.bhv_times);
    
    [lever_active,lever_force_resample] = ...
        AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
    
    movement_starts = find(diff([0;lever_active;0]) == 1);
    movement_stops = find(diff([0;lever_active;0]) == -1);
    movement_durations = movement_stops - movement_starts;
    movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
    pre_movement_iti = [NaN;movement_iti];
    post_movement_iti = [movement_iti;NaN];
        % TO USE ANY MOVEMENTS THAT FULFILL TIME CRITERIA
    min_move_time = 1000;
    max_move_time = 4000;
    min_pre_move_iti = 1000;
    min_post_move_iti = 1000;
    use_movements = find(movement_durations > min_move_time & ...
        movement_durations < max_move_time & ...
        pre_movement_iti > min_pre_move_iti & ...
        post_movement_iti > min_post_move_iti);
    
    [~,use_movement_start_frames] = ...
        arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
        movement_starts(use_movements(x)))), 1:length(use_movements));
    
    [~,use_movement_stop_frames] = ...
        arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
        movement_stops(use_movements(x)))), 1:length(use_movements));
    
    % TO USE CUED REWARDED MOVEMENTS
    cr_move_onset_frames = nan(length(intersect(find(rewarded_trials)',curr_data.bhv.trial_list(:,2)')),1);
    cr_move_offset_frames = nan(length(intersect(find(rewarded_trials)',curr_data.bhv.trial_list(:,2)')),1);   
    curr_idx = 1;
    for curr_trial = intersect(find(rewarded_trials)',curr_data.bhv.trial_list(:,2)');
        
        % Get the offset between xsg and bhv
        curr_list_idx = find(curr_data.bhv.trial_list(:,2) == curr_trial);
        curr_xsg_bhv_offset = curr_data.bhv.trial_list(curr_list_idx,1) - ...
            curr_data.bhv.bhv_times{curr_trial}.states.bitcode(1);
        
        % Get the reward in xsg time (in ms)
        if isfield(curr_data.bhv.bhv_times{curr_trial}.states,'reward_delay')
            curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward_delay(1) + ...
                curr_xsg_bhv_offset)*1000);
        else
            curr_reward_time = round((curr_data.bhv.bhv_times{curr_trial}.states.reward(1) + ...
                curr_xsg_bhv_offset)*1000);
        end
        
        % Get the boundaries of the movement that overlaps with reward
        curr_rewarded_movement_onset = curr_reward_time - find( ...
            ~lever_active(curr_reward_time:-1:1),1) + 2;
        curr_rewarded_movement_offset = curr_reward_time + find( ...
            ~lever_active(curr_reward_time:end),1) - 2;
        
        [~,curr_move_onset_frames(curr_trial)] = ...
            min(abs(curr_data.bhv.frame_times*1000 - ...
            curr_rewarded_movement_onset));
        
        [~,curr_move_offset_frames(curr_trial)] = ...
            min(abs(curr_data.bhv.frame_times*1000 - ...
            curr_rewarded_movement_offset));
        
        cr_move_onset_frames(curr_idx) = curr_move_onset_frames(curr_trial);
        cr_move_offset_frames(curr_idx) = curr_move_offset_frames(curr_trial);
        
        curr_idx = curr_idx+1;
    end
    
    use_frames = use_movement_start_frames;
    
    % load frames around movement
    pathsep = strfind(curr_file,'/');
    curr_day = curr_file(pathsep(end)+1:pathsep(end)+6);
    data_path = [data_dir filesep animal filesep curr_day];
    tiff_dir = dir([data_path filesep '*.tif']);
    tiff_files = cellfun(@(x) [data_path filesep x], sort({tiff_dir.name}),'uni',false);
          
    frames_premoveonset = 30;
    frames_postmoveonset = 15;
    
    frames_premoveoffset = 15;
    frames_postmoveoffset = 30;
    
    move_start_surround_im = zeros(512,512,frames_premoveonset + frames_postmoveonset + 1);
    move_stop_surround_im = zeros(512,512,frames_premoveoffset + frames_postmoveoffset + 1);

    %%% movement onset
    move_onset_files = ceil(use_movement_start_frames/800);
    move_onset_file_frames = rem(use_movement_start_frames,800);
 
    num_ims = zeros(size(move_start_surround_im,3),1);
    for i = 1:length(use_movement_start_frames)
        
        curr_file = tiff_files{move_onset_files(i)};
        curr_frames_back = move_onset_file_frames(i) - max(1,move_onset_file_frames(i) - frames_premoveonset);
        curr_frames_forward = min(800,move_onset_file_frames(i) + frames_postmoveonset) - move_onset_file_frames(i);
        
        imageinfo=imfinfo(curr_file,'tiff');
        
        for curr_relative_frame = -curr_frames_back:curr_frames_forward;
            
            curr_load_frame = move_onset_file_frames(i) + curr_relative_frame;
            
            curr_im = zeros(512,'int16');
            curr_im = imread(curr_file,'tiff',curr_load_frame,'Info',imageinfo);
            
            curr_surround_frame = frames_premoveonset + 1 + curr_relative_frame;
            move_start_surround_im(:,:,curr_surround_frame) = ...
                move_start_surround_im(:,:,curr_surround_frame) + double(curr_im);
            
            num_ims(curr_surround_frame) = num_ims(curr_surround_frame) + 1;
        end
        
        disp(i/length(use_movement_start_frames));
        
    end
    
    move_start_surround_im_avg = zeros(size(move_start_surround_im));
    for i = 1:size(move_start_surround_im_avg,3)
        move_start_surround_im_avg(:,:,i) = move_start_surround_im(:,:,i)./num_ims(i);
    end
    
    %%% movement offset 
    move_offset_files = ceil(use_movement_stop_frames/800);
    move_offset_file_frames = rem(use_movement_stop_frames,800);
 
    num_ims = zeros(size(move_stop_surround_im,3),1);
    for i = 1:length(use_movement_stop_frames)
        
        curr_file = tiff_files{move_offset_files(i)};
        curr_frames_back = move_offset_file_frames(i) - max(1,move_offset_file_frames(i) - frames_premoveoffset);
        curr_frames_forward = min(800,move_offset_file_frames(i) + frames_postmoveoffset) - move_offset_file_frames(i);
        
        imageinfo=imfinfo(curr_file,'tiff');
        
        for curr_relative_frame = -curr_frames_back:curr_frames_forward;
            
            curr_load_frame = move_offset_file_frames(i) + curr_relative_frame;
            
            curr_im = zeros(512,'int16');
            curr_im = imread(curr_file,'tiff',curr_load_frame,'Info',imageinfo);
            
            curr_surround_frame = frames_premoveoffset + 1 + curr_relative_frame;
            move_stop_surround_im(:,:,curr_surround_frame) = ...
                move_stop_surround_im(:,:,curr_surround_frame) + double(curr_im);
            
            num_ims(curr_surround_frame) = num_ims(curr_surround_frame) + 1;
        end
        
        disp(i/length(use_movement_stop_frames));
        
    end
    
    move_stop_surround_im_avg = zeros(size(move_stop_surround_im));
    for i = 1:size(move_stop_surround_im_avg,3)
        move_stop_surround_im_avg(:,:,i) = move_stop_surround_im(:,:,i)./num_ims(i);
    end
    
    %%% combine
    move_surround_im_avg = cat(3,move_start_surround_im_avg, ...
        move_stop_surround_im_avg);
    
    
    % normalize image, prepare for uint16 tiff viewing
    baseline_img = nanmedian(move_surround_im_avg(:,:,1:frames_premoveonset-20),3);
    
    smooth_filter = fspecial('gaussian',20,2);
    
    move_surround_im_df = cell2mat(arrayfun(@(x) imfilter((move_surround_im_avg(:,:,x) ...
        - baseline_img),smooth_filter), permute(1:size(move_surround_im_avg,3),[1 3 2]),'uni',false));
    
    % make lowest value nonnegative for tiff purposes
    move_surround_im_df = (move_surround_im_df - min(move_surround_im_df(:)))/ ...
        (max(move_surround_im_df(:)) - min(move_surround_im_df(:)))*1000;
    
    save_path = '/usr/local/lab/People/Andy/Corticospinal/dated_analyses/150327/move_surround_videos';
    save_filename = [save_path filesep animal '_s' num2str(curr_session)];
    for i = 1:size(move_surround_im_avg,3)
        imwrite(uint16(move_surround_im_avg(:,:,i)), ...
            [save_filename '.tif'], ...
            'writemode','append','Compression','none');
        
        imwrite(uint16(move_surround_im_df(:,:,i)), ...
            [save_filename 'df.tif'], ...
            'writemode','append','Compression','none');
    end
    
    
end

disp('done');

end

%% Slide 7) Get activity around selected movements and plot corrcoef

% NOTE! THIS NEEDS TO BE NARROWED TO ONLY IMAGED MOVEMENTS (probably)




binned_act_corr_sessiondiff = cell(length(animals),1);
binned_act_corr_session = cell(length(animals),1);

%animal = 'AP118';
for curr_animal = 1:length(animals);
    animal = animals{curr_animal};
%for curr_animal = 1:length(data)
    movements = cell(14,1);
    activity = cell(14,1);
    
    %n_sessions = length(data(curr_animal).im);
    %for curr_session = 1:n_sessions;
        
        % Use this if data not loaded in
            data_dir = '/usr/local/lab/People/Andy/Data';
            curr_analysis_dir = [data_dir filesep animal filesep ...
                animal '_batch_thresh_roi'];
            curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
            n_sessions = length(curr_analysis_files);
            for curr_session = 1:n_sessions
            curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
            curr_data = load(curr_file);
            
            % Threshold data
            event_thresh = 3;
            if exist('d_grp_all','var')
                % If ROIs grouped, then average together and threshold
                curr_df_grp = cell2mat(cellfun(@(x) nanmean(curr_data.im.roi_trace_df(x,:),1), ...
                    d_grp_all{curr_animal},'uni',false));
                curr_df_thresh = AP_caEvents_thresh(curr_df_grp,event_thresh);
            else
                % If not, just threshold all
                curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
            end
        
            
        
%         % Use this if data loaded in
%         curr_data.im = data(curr_animal).im(curr_session);
%         curr_data.bhv = data(curr_animal).bhv(curr_session);
%         curr_df_thresh = data(curr_animal).im(curr_session).roi_trace_thresh;
        
        
        [lever_active,lever_force_resample] = ...
            AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        movement_starts = find(diff([0;lever_active;0]) == 1);
        movement_stops = find(diff([0;lever_active;0]) == -1);
        movement_durations = movement_stops - movement_starts;
        movement_iti = movement_starts(2:end) - movement_stops(1:end-1);
        pre_movement_iti = [NaN;movement_iti];
        post_movement_iti = [movement_iti;NaN];
        
        % Set criteria for pulling out movements
        min_move_time = 200;
        max_move_time = 3000;
        min_pre_move_iti = 1000;
        min_post_move_iti = 1000;
        use_movements = find(movement_durations > min_move_time & ...
            movement_durations < max_move_time & ...
            pre_movement_iti > min_pre_move_iti & ...
            post_movement_iti > min_post_move_iti);
        criteria_movements = movement_starts(use_movements);
        
        [~,use_movement_start_frames] = ...
            arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
            criteria_movements(x))), 1:length(criteria_movements));
        
        
        frames_back = 0;
        frames_forward = round(max_move_time./curr_data.bhv.framerate);
        valid_frames = (use_movement_start_frames - frames_back > 0) & ...
            (use_movement_start_frames + frames_forward < size(curr_data.im.roi_trace_df,2));
        
        
        activity{curr_session} = cell2mat(cellfun(@(x) curr_df_thresh(:, ...
            x-frames_back:x+frames_forward)',permute(num2cell( ...
            use_movement_start_frames(valid_frames)),[1 3 2]),'uni',false));
        
        movements{curr_session} = cell2mat(cellfun(@(x) lever_force_resample(x:x+max_move_time), ...
            num2cell(criteria_movements(valid_frames))','uni',false));
        
        disp(curr_session)
    end
    
    % Plot all act vs. move corr
    % All temporal activity
    activity_reshape = cellfun(@(x) reshape(x,[],size(x,3)),activity,'uni',false);
    % Binary active / not active
    %activity_reshape = cellfun(@(x) permute(any(x,1),[2 3 1]),activity,'uni',false);
    
    movements_use = cellfun(@(x) x,movements,'uni',false);
    
    act_corr = AP_itril(corrcoef(horzcat(activity_reshape{:}),'rows','complete'),-1);
    move_corr = AP_itril(corrcoef(horzcat(movements_use{:}),'rows','complete'),-1);
    
    corr_bins = linspace(-1,1,30);
    corr_bins(end) = Inf;
    
    [~,bin_idx] = histc(move_corr,corr_bins);
    move_corr_binned = grpstats(act_corr,bin_idx,'mean');
    
    figure;plot(move_corr_binned,'k','linewidth',2);
    
    session_grp = num2cell(1:n_sessions);%{[1:3] [4:6] [7:9] [10:12] [13:14]}%
    
    % Plot by session separation
    sessions = cell2mat(arrayfun(@(x) repmat(x,size(movements_use{x},2),1),[1:length(movements_use)]','uni',false));
    sessions_tril_1 = AP_itril(repmat(sessions,1,length(sessions)),-1);
    sessions_tril_2 = AP_itril(repmat(sessions,1,length(sessions))',-1);
    session_diff = AP_itril(abs(repmat(sessions,1,length(sessions)) - repmat(sessions,1,length(sessions))'),-1);
    
    col = jet(length(session_grp));
    figure; hold on;
    binned_act_corr_sessiondiff{curr_animal} = nan(length(corr_bins)-1,14);
    for i = 1:length(session_grp)
        
        [~,bin_idx] = histc(move_corr(ismember(session_diff,session_grp{i}-1)),corr_bins);
        act_corr_binned = grpstats(act_corr(ismember(session_diff,session_grp{i}-1)),bin_idx,'mean');
        
        use_bins = unique(bin_idx);
        plot(use_bins,act_corr_binned,'color',col(i,:),'linewidth',2);
        
        binned_act_corr_sessiondiff{curr_animal}(use_bins,i) = act_corr_binned;
        
    end
    
    % Plot by session
    figure; hold on;
    binned_act_corr_session{curr_animal} = nan(length(corr_bins)-1,14);
    for i = 1:length(session_grp)
        
        curr_act_corr = AP_itril(corrcoef(horzcat(activity_reshape{session_grp{i}}),'rows','complete'),-1);
        curr_move_corr = AP_itril(corrcoef(horzcat(movements_use{session_grp{i}}),'rows','complete'),-1);
        
        [~,bin_idx] = histc(curr_move_corr,corr_bins);
        act_corr_binned = grpstats(curr_act_corr,bin_idx,'mean');
        
        use_bins = unique(bin_idx);
        plot(use_bins,act_corr_binned,'color',col(i,:),'linewidth',2);
        
        binned_act_corr_session{curr_animal}(use_bins,i) = act_corr_binned;
        
    end
    
end

% can also try this two ways: average the correlation of trials, or average
% the trials together and take the correlation




% This was to look at some day separation but restricted to certain time
% period (i.e. is the change the same amount every day, or more in the
% beginning than the end, etc?) I don't know yet how to check whether this
% is just some calcium imaging artifact, which I'm guessing it is...
% col = jet(length(session_grp));
% figure; hold on;
% for i = 1:length(session_grp)
%     
%     use_pair = sessions_tril_1 > 10 & sessions_tril_2 > 10 & ismember(session_diff,session_grp{i}-1);
%     
%     [~,bin_idx] = histc(move_corr(use_pair),corr_bins);
%     act_corr_binned = grpstats(act_corr(use_pair),bin_idx,'mean');
%     
%     use_bins = unique(bin_idx);
%     plot(use_bins,act_corr_binned,'color',col(i,:),'linewidth',2);
%     
% end


% plot when all animals loaded in
corr_bins_plot = corr_bins(1:end-1) + nanmedian(diff(corr_bins(1:end-1)));
use_bins = 4:length(corr_bins_plot) - 5;

all_plot = cat(3,binned_act_corr_sessiondiff{:});
col = jet(14);
figure; hold on;
for i = 1:size(all_plot,2);

    plot(corr_bins_plot(use_bins),nanmean(all_plot(use_bins,i,:),3),'color',col(i,:),'linewidth',2);
    
end

all_plot = cat(3,binned_act_corr_session{:});
figure; hold on;
for i = 1:size(all_plot,2);

    plot(corr_bins_plot(use_bins),nanmean(all_plot(use_bins,i,:),3),'color',col(i,:),'linewidth',2);
    
end

% early within session v late within session
all_plot = cat(3,binned_act_corr_session{:});
early_within = permute(nanmean(all_plot(:,1:4,:),2),[1 3 2]);
late_within = permute(nanmean(all_plot(:,11:14,:),2),[1 3 2]);

figure; hold on

%errorbar(nanmean(early_within,2),nanstd(early_within,[],2)./sqrt(sum(~isnan(early_within),2)),'k','linewidth',2)
%errorbar(nanmean(late_within,2),nanstd(late_within,[],2)./sqrt(sum(~isnan(late_within),2)),'r','linewidth',2)

early_within_mean = nanmean(early_within,2);
early_within_sem = nanstd(early_within,[],2)./sqrt(sum(~isnan(early_within),2));
late_within_mean = nanmean(late_within,2);
late_within_sem = nanstd(late_within,[],2)./sqrt(sum(~isnan(late_within),2));

plot(corr_bins_plot(use_bins),early_within_mean(use_bins),'k','linewidth',2)
plot(corr_bins_plot(use_bins),late_within_mean(use_bins),'r','linewidth',2)

jbfill(corr_bins_plot(use_bins),early_within_mean(use_bins)'-early_within_sem(use_bins)', ...
    early_within_mean(use_bins)'+early_within_sem(use_bins)','k','k');
jbfill(corr_bins_plot(use_bins),late_within_mean(use_bins)'-late_within_sem(use_bins)', ...
    late_within_mean(use_bins)'+late_within_sem(use_bins)','r','r');

legend({'S1-4','S11-14'})
ylabel('Activity correlation')
xlabel('Movement correlation')


%% Get activity around selected quiescence

animal = 'AP124';

quiescents = cell(14,1);
activity = cell(14,1);
for curr_session = 1:14;    
    
    data_dir = '/usr/local/lab/People/Andy/Data';
    curr_analysis_dir = [data_dir filesep animal filesep ...
        animal '_batch_thresh_roi'];
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
    curr_data = load(curr_file);
    
    
    num_trials = length(curr_data.bhv.bhv_times);
    rewarded_trials = cellfun(@(x) ~isempty(x.states.reward),curr_data.bhv.bhv_times);
    
    [lever_active,lever_force_resample] = ...
        AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
    
    quiescent_starts = find(diff([1;lever_active;1]) == -1);
    quiescent_stops = find(diff([1;lever_active;1]) == 1);
    quiescent_durations = quiescent_stops - quiescent_starts;
    quiescent_iti = quiescent_starts(2:end) - quiescent_stops(1:end-1);
    pre_quiescent_iti = [NaN;quiescent_iti];
    post_quiescent_iti = [quiescent_iti;NaN];
    
    % Set criteria for pulling out movements
    min_quiescent_time = 1500;
    max_quiescent_time = 2500;
    min_pre_quiescent_iti = 1000;
    min_post_quiescent_iti = 1000;
    use_quiescent = find(quiescent_durations > min_quiescent_time & ...
        quiescent_durations < max_quiescent_time & ...
        pre_quiescent_iti > min_pre_quiescent_iti & ...
        post_quiescent_iti > min_post_quiescent_iti);
    criteria_quiescent = quiescent_starts(use_quiescent);
    
    [~,use_quiescent_start_frames] = ...
        arrayfun(@(x) min(abs(curr_data.bhv.frame_times*1000 - ...
        criteria_quiescent(x))), 1:length(criteria_quiescent));
    
    
    frames_back = 30;
    frames_forward = 90;
    valid_frames = (use_quiescent_start_frames - frames_back > 0) & ...
        (use_quiescent_start_frames + frames_forward < size(curr_data.im.roi_trace_df,2));
      
    activity{curr_session} = cell2mat(cellfun(@(x) curr_data.im.roi_trace_df(:, ...
        x-frames_back:x+frames_forward)',permute(num2cell( ...
        use_quiescent_start_frames(valid_frames)),[1 3 2]),'uni',false));
    
    quiescents{curr_session} = cell2mat(cellfun(@(x) lever_force_resample(x-1000:x+3000), ...
        num2cell(criteria_quiescent(valid_frames))','uni',false));
    
    disp(curr_session)
end

% Plot activity by level of lever (NOTE: this could just vary by day if the
% lever is in a different place)
quiescents_cat = horzcat(quiescents{:});
q_level = nanmean(quiescents_cat(1001:2500,:),1);
activity_reshape = cellfun(@(x) reshape(x(31:60,:,:),[],size(x,3)),activity,'uni',false);

bin_edges = linspace(min(q_level),max(q_level),5);
[n bins] = histc(q_level,bin_edges);
act_grp = grpstats(horzcat(activity_reshape{:})',bins);

% Check for correlation of q_level with session (if there is one, bad)
sessions = cell2mat(cellfun(@(x) repmat(x,size(quiescents{x},2),1), ...
    num2cell([1:length(quiescents)]'),'uni',false));
[r,p] = corrcoef(sessions,q_level); % it looks like there is one



quiescents_cat = horzcat(quiescents{:});
q_cat_level = nanmean(quiescents_cat(1001:2500,:),1);
bin_edges = linspace(min(q_cat_level),max(q_cat_level),8);
figure; hold on
col = jet(14);
binned_act = nan(length(bin_edges),14);
for curr_session = 1:14
    q_level = nanmean(quiescents{curr_session}(1001:2500,:),1);
    if length(q_level) < 5;
        continue
    end
    
    [n bins] = histc(q_level,bin_edges);
    act_reshape = reshape(activity{curr_session}(31:60,:,:),[],size(activity{curr_session},3));
    act_grp = grpstats(act_reshape',bins);
    
    binned_act(unique(bins),curr_session) = zscore(nanmean(act_grp,2));
    
    plot(unique(bins),zscore(nanmean(act_grp,2)),'color',col(curr_session,:));
end




% do something like this during movement too? highest/lowest position?




%% Slide 8) Get up/down pattern stuff

% Create vector with average activity of all cells for select days
use_sessions = 1:4;

curr_act = cell(length(analysis),1);
for curr_animal = 1:length(analysis)   
    use_sessions_animal = intersect(use_sessions,1:length(analysis(curr_animal).im));
    curr_act{curr_animal} = ...
        permute(nanmean(vertcat(analysis(curr_animal).im(use_sessions_animal).move_onset_aligned),1),[3 2 1]);   
end
curr_act_cat = vertcat(curr_act{:});

figure; hold on;

t = analysis(1).surrounding_frames./28;

start_frames = 1:61;
up_col = autumn(length(start_frames));
down_col = cool(length(start_frames));

max_up = nan(length(start_frames),1);
min_down = nan(length(start_frames),1);
for start_frame = start_frames;
    
    curr_idx = find(start_frames == start_frame);
    
    pattern = [zeros(15,1);ones(15,1)];
    corr_grid = corrcoef([pattern curr_act_cat(:,start_frame:start_frame+length(pattern)-1)']);
    pattern_corr = corr_grid(1,2:end);
    
    mean_up = nanmean(curr_act_cat(pattern_corr > 0,:));
    mean_down = nanmean(curr_act_cat(pattern_corr < 0,:));
    
    diff_up(curr_idx) = max(mean_up) - min(mean_up);
    diff_down(curr_idx) = max(mean_down) - min(mean_down);
    
    % Plot every nth
    if mod(start_frame,3) == 0;
        subplot(2,1,1);hold on;
        plot(t,mean_up,'color',up_col(curr_idx,:),'linewidth',2);
        subplot(2,1,2);hold on;
        plot(t,mean_down,'color',down_col(curr_idx,:),'linewidth',2);
    end
    
end
subplot(2,1,1);
line([0 0],ylim,'color','k','linewidth',2,'linestyle','--');
subplot(2,1,2);
line([0 0],ylim,'color','k','linewidth',2,'linestyle','--');


% Find the greatest difference for both
[~,max_idx] = max(diff_up + diff_down);
max_start_frame = start_frames(max_idx);

pattern = [zeros(15,1);ones(15,1)];
[corr_grid corr_p]= corrcoef([pattern curr_act_cat(:,max_start_frame:max_start_frame+length(pattern)-1)']);
pattern_corr = corr_grid(1,2:end);
pattern_corr_p = corr_p(1,2:end);

mean_up = nanmean(curr_act_cat(pattern_corr > 0,:));
mean_down = nanmean(curr_act_cat(pattern_corr < 0,:));

subplot(2,1,1);
plot(t,mean_up,'k','linewidth',2)
ylabel('Average activity')
xlabel('Time (s)')
subplot(2,1,2);
plot(t,mean_down,'k','linewidth',2);
ylabel('Average activity')
xlabel('Time (s)')

% Histogram everything
bin_edges = linspace(-1,1,20);
bin_edges(end) = Inf;
n_all = histc(pattern_corr,bin_edges);
n_sig = histc(pattern_corr(pattern_corr_p < 0.05),bin_edges);
figure; hold on;
bar(bin_edges(1:end-1),n_all(1:end-1)/length(pattern_corr),'FaceColor','k')
bar(bin_edges(1:end-1),n_sig(1:end-1)/length(pattern_corr),'FaceColor','r')
line([0 0],ylim,'color','k','linewidth',2,'linestyle','--');
ylabel('Frequency (fraction)')
xlabel('Down/Up correlation')


% (as above, but do up/down seperately)
% Find the greatest difference for both
% [~,max_up_idx] = max(diff_up);
% [~,max_down_idx] = max(diff_down);
% 
% pattern = [zeros(15,1);ones(15,1)];
% 
% max_up_start_frame = start_frames(max_up_idx);
% corr_grid = corrcoef([pattern b(:,max_up_start_frame:max_up_start_frame+length(pattern)-1)']);
% pattern_corr = corr_grid(1,2:end);
% mean_up = nanmean(b(pattern_corr > 0,:));
% 
% max_down_start_frame = start_frames(max_down_idx);
% corr_grid = corrcoef([pattern b(:,max_up_start_frame:max_up_start_frame+length(pattern)-1)']);
% pattern_corr = corr_grid(1,2:end);
% mean_down = nanmean(b(pattern_corr < 0,:));
% 
% subplot(2,1,1);
% plot(mean_up,'k','linewidth',2)
% subplot(2,1,2);
% plot(mean_down,'k','linewidth',2);

% Get every cell, every session
curr_act = cell(length(analysis),1);
pattern_corr = cell(length(analysis),1);
pattern_corr_p = cell(length(analysis),1);

pattern = [zeros(15,1);ones(15,1)];

for curr_animal = 1:length(analysis)  
    
    curr_act{curr_animal} = ...
        arrayfun(@(x) permute(nanmean( ...
        analysis(curr_animal).im(x).move_onset_aligned,1),[3 2 1]), ...
        1:length(analysis(curr_animal).im),'uni',false);   
    
    use_sessions = cellfun(@(x) ~isempty(x),curr_act{curr_animal});
    
    [curr_corr_grid curr_corr_p] = cellfun(@(x) ...
        corrcoef([pattern x(:,max_start_frame:max_start_frame+ ...
        length(pattern)-1)']),curr_act{curr_animal}(use_sessions),'uni',false);
    
    pattern_corr{curr_animal}(use_sessions,:) = cell2mat(cellfun(@(x) x(1,2:end),curr_corr_grid','uni',false));
    pattern_corr_p{curr_animal}(use_sessions,:) = cell2mat(cellfun(@(x) x(1,2:end),curr_corr_p','uni',false));
    
end

% Plot fraction active, up, down
frac_active_rois = cell2mat(cellfun(@(x,y) [nanmean(y < 0.05,2);nan(14-size(x,1),1)], ...
    pattern_corr',pattern_corr_p','uni',false));
frac_up_rois = cell2mat(cellfun(@(x,y) [nanmean(x > 0 & y < 0.05,2);nan(14-size(x,1),1)], ...
    pattern_corr',pattern_corr_p','uni',false));
frac_down_rois = cell2mat(cellfun(@(x,y) [nanmean(x < 0 & y < 0.05,2);nan(14-size(x,1),1)], ...
    pattern_corr',pattern_corr_p','uni',false));

figure; hold on;
errorbar(nanmean(frac_active_rois,2),nanstd(frac_active_rois,[],2)./sqrt(sum(~isnan(frac_active_rois),2)),'k','linewidth',2);
errorbar(nanmean(frac_up_rois,2),nanstd(frac_up_rois,[],2)./sqrt(sum(~isnan(frac_up_rois),2)),'g','linewidth',2);
errorbar(nanmean(frac_down_rois,2),nanstd(frac_down_rois,[],2)./sqrt(sum(~isnan(frac_down_rois),2)),'r','linewidth',2);
legend({'All' 'Up' 'Down'})
ylabel('Fraction of ROIs')
xlabel('Session')

% Plot turnover of up/down rois
up_rois = cellfun(@(x,y) x > 0 & y < 0.05,pattern_corr,pattern_corr_p,'uni',false);
down_rois = cellfun(@(x,y) x < 0 & y < 0.05,pattern_corr,pattern_corr_p,'uni',false);

% For day before
% up_history_frac = cellfun(@(x,y) [ ...
%     (sum(x(2:end,:) & x(1:end-1,:),2))./sum(x(2:end,:),2) ...
%     (sum(x(2:end,:) & y(1:end-1,:),2))./sum(x(2:end,:),2); nan(13-(size(x,1)-1),2)], ...
%     up_rois,down_rois,'uni',false);
% 
% down_history_frac = cellfun(@(x,y) [ ...
%     (sum(x(2:end,:) & x(1:end-1,:),2))./sum(x(2:end,:),2) ...
%     (sum(x(2:end,:) & y(1:end-1,:),2))./sum(x(2:end,:),2); nan(13-(size(x,1)-1),2)], ...
%     down_rois,up_rois,'uni',false);

% For ever before
up_history_frac = cellfun(@(x,y) [ ...
    (sum(x(2:end,:) & cumsum(x(1:end-1,:)) > 0,2))./sum(x(2:end,:),2) ...
    (sum(x(2:end,:) & cumsum(y(1:end-1,:)) > 0,2))./sum(x(2:end,:),2); nan(13-(size(x,1)-1),2)], ...
    up_rois,down_rois,'uni',false);

down_history_frac = cellfun(@(x,y) [ ...
    (sum(x(2:end,:) & cumsum(x(1:end-1,:)) > 0,2))./sum(x(2:end,:),2) ...
    (sum(x(2:end,:) & cumsum(y(1:end-1,:)) > 0,2))./sum(x(2:end,:),2); nan(13-(size(x,1)-1),2)], ...
    down_rois,up_rois,'uni',false);


up_history_cat = cat(3,up_history_frac{:});
up_history_mean = nanmean(up_history_cat,3);
up_history_sem = nanstd(up_history_cat,[],3)./sqrt(sum(~isnan(up_history_cat),3));

down_history_cat = cat(3,down_history_frac{:});
down_history_mean = nanmean(down_history_cat,3);
down_history_sem = nanstd(down_history_cat,[],3)./sqrt(sum(~isnan(down_history_cat),3));

figure;
subplot(2,1,1); hold on;
errorb(up_history_mean,up_history_sem)
colormap(gray);
ylabel('Fraction of up rois')
xlabel('Session')

subplot(2,1,2); hold on;
errorb(down_history_mean,down_history_sem)
colormap(gray);
ylabel('Fraction of down rois')
xlabel('Session')
legend({'Previously same' 'Previously opposite'})


%% PCA data and get correlation over time

act_session_corr = nan(14,14,length(analysis));
pca_session_corr = nan(14,14,length(analysis));
for curr_animal = 1:length(analysis)
    
    use_sessions = cellfun(@(x) ~isempty(x), {analysis(curr_animal).im(:).move_onset_aligned});
    
    temp_act = cellfun(@(act) zscore(permute(nanmean(act),[3 2 1]),[],2), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    temp_act_cat = horzcat(temp_act{:});
    temp_act_cat(isnan(temp_act_cat)) = 0;
    
    [coeff score latent] = princomp(temp_act_cat');
    act_size = cellfun(@(x) size(x,2),temp_act);
    score_sessions = mat2cell(score,act_size);
    
    score_cat = cell2mat(cellfun(@(x) reshape(x(:,1:3),[],1),score_sessions','uni',false));
    score_cat_corr = corrcoef(score_cat);   
    pca_session_corr(use_sessions,use_sessions,curr_animal) = score_cat_corr;
    
    act_reshape = cell2mat(cellfun(@(x) reshape(x',[],1),temp_act,'uni',false));
    act_corr = corrcoef(act_reshape);
    act_session_corr(use_sessions,use_sessions,curr_animal) = act_corr;
end

figure;imagesc(nanmean(act_session_corr,3));colormap(hot);
figure;imagesc(nanmean(pca_session_corr,3));colormap(hot);
% NOTE: pop corr from upper cell is better b/c deals w nans seperately


%% Slide 9) Get differences /correlation of differences between days

pop_corr = nan(14,14,length(analysis));
pop_diff_corr = nan(13,13,length(analysis));
pop_diff_corr_allcomp = nan(13,13,length(analysis));

for curr_animal = 1:length(analysis);
    
    use_sessions = find(cellfun(@(x) ~isempty(x), {analysis(curr_animal).im(:).move_onset_aligned}));
    
    % All temporal activity
    %temp_act = cellfun(@(act) reshape(permute(nanmean(act),[2 3 1]),[],1), ...
    %    {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    % Mean binary activity (% active trials)
    temp_act = cellfun(@(act) permute(nanmean(any(act,2),1),[3 2 1]), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
    
    temp_act_cat = horzcat(temp_act{:});

   %%%% TEMPORARY SHUFFLING
   
%    %temp_act_cat = shake(temp_act_cat,2);
%    %td = diff(temp_act_cat,[],2);
%    %td = medfilt1(td,5,[],2);
%    %td = td.*sign(randn(size(td)));
%    %td = shake(td,2);
%    %td = td(:,randperm(size(td,2)));
%    %td = td(randperm(size(td,1)),:);
%    temp_act_cat = temp_act{1};
%    for i = 1:size(td,2)
%        temp_act_cat(:,i+1) = temp_act_cat(:,i) + td(:,i);
%    end
    
   %%%%
    
    pop_corr(use_sessions,use_sessions,curr_animal) = corrcoef(temp_act_cat);
    
    temp_act_diff_all = cell(size(temp_act_cat,2));
    for pre = 1:size(temp_act_cat,2)
        for post = 1:size(temp_act_cat,2)
            temp_act_diff_all{pre,post} = temp_act_cat(:,post) - ...
                temp_act_cat(:,pre);
        end
    end
        
    diff_inout = cell(size(temp_act_cat,2),1);
    for curr_session = 1:size(temp_act_cat,2)   
        diff_in = horzcat(temp_act_diff_all{:,curr_session});
        diff_out = horzcat(temp_act_diff_all{curr_session,:});
        diff_inout{curr_session} = [diff_in diff_out];        
    end
    diff_inout_allcorr = mat2cell(corrcoef(horzcat(diff_inout{:}),'rows','complete'), ...
        repmat(size(temp_act_cat,2)*2,size(temp_act_cat,2),1), ...
        repmat(size(temp_act_cat,2)*2,size(temp_act_cat,2),1));
    diff_inout_corr_mean = cellfun(@(x) nanmean(AP_itril(x,-size(temp_act_cat,2)+1)), ...
        diff_inout_allcorr);
    
    pop_diff_corr_allcomp(use_sessions,use_sessions,curr_animal) = diff_inout_corr_mean;
    
    
    temp_act_diff = diff(temp_act_cat,[],2);
    act_diff_corr = corrcoef(temp_act_diff);
     
    use_sessions_diff = use_sessions(1:end-1);
    pop_diff_corr(use_sessions_diff,use_sessions_diff,curr_animal) = act_diff_corr;
end

figure;

subplot(1,3,1);
imagesc(nanmean(pop_corr,3));
colormap(gray);
caxis([0 1]);
xlabel('Session')
ylabel('Session')
title('Pop corr')

subplot(1,3,2);imagesc(nanmean(pop_diff_corr,3));
colormap(gray);
caxis([-1 1]);
xlabel('Session')
ylabel('Session')
title('Pop diff corr')

subplot(1,3,3);imagesc(nanmean(pop_diff_corr_allcomp,3));
colormap(gray);
caxis([-1 1]);
ylabel('Forward difference (session)')
xlabel('Backward difference (session)')
title('Pop diff corr allcomp');



%% Slide 11) Up/Down rois by reliability + chi-square

premove_frames = 1:20;
move_frames = 31:51;
for curr_animal = 1:length(analysis)   
    use_sessions_animal = intersect(use_sessions,1:length(analysis(curr_animal).im));
    
    timing_reliability = cellfun(@(x) [ ...
        permute(nanmean(any(x(:,premove_frames,:) > 0,2),1),[3 2 1]), ...
        permute(nanmean(any(x(:,move_frames,:) > 0,2),1),[3 2 1])], ...
        {analysis(curr_animal).im(use_sessions_animal).move_onset_aligned},'uni',false);
    
    comparison_chi_sq = cellfun(@(x) AP_chisquare_matrix( ...
        permute(any(x(:,premove_frames,:) > 0,2),[1 3 2]), ...
        permute(any(x(:,move_frames,:) > 0,2),[1 3 2])), ...
        {analysis(curr_animal).im(use_sessions_animal).move_onset_aligned},'uni',false);

    
end

% Group sessions, look at changes
session_grp = {[1:3],[12:14]};
premove_frames = 1:20;
move_frames = 31:51;
timing_reliability = cell(length(animals),length(session_grp));
comparison_chi_sq = cell(length(animals),length(session_grp));
for curr_animal = 1:length(analysis)   
    
    use_sessions = find(cellfun(@(x) ~isempty(x), {analysis(curr_animal).im(:).move_onset_aligned}));
    curr_sessions_grp = cellfun(@(x) intersect(use_sessions,x),session_grp,'uni',false);
    if any(cellfun(@isempty,curr_sessions_grp)) 
        continue
    end
    
    timing_reliability(curr_animal,:) = cellfun(@(x) [ ...
        permute(nanmean(any(x(:,premove_frames,:) > 0,2),1),[3 2 1]), ...
        permute(nanmean(any(x(:,move_frames,:) > 0,2),1),[3 2 1])], ...
        cellfun(@(x) vertcat(analysis(curr_animal).im(x).move_onset_aligned), ...
        curr_sessions_grp,'uni',false),'uni',false);

    comparison_chi_sq(curr_animal,:) = cellfun(@(x) AP_chisquare_matrix( ...
        permute(any(x(:,premove_frames,:) > 0,2),[1 3 2]), ...
        permute(any(x(:,move_frames,:) > 0,2),[1 3 2])), ...
        cellfun(@(x) vertcat(analysis(curr_animal).im(x).move_onset_aligned), ...
        curr_sessions_grp,'uni',false),'uni',false);
    
end

reliability_diff = cell2mat(cellfun(@(x) diff(x,[],2),timing_reliability,'uni',false));
sig_rois = cell2mat(cellfun(@(x) x' < 0.05,comparison_chi_sq,'uni',false));

% Plot density plot of early v. late differences
bin_edges = linspace(-1,1,200);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
bin_label = round(bin_centers*100)/100;
bin_edges(1) = -Inf;
bin_edges(end) = Inf;
reliability_density = hist3(fliplr(reliability_diff(any(sig_rois,2),:)),'Edges',repmat({bin_edges},2,1));
h = fspecial('gaussian',10,2);
reliability_density_smooth = imfilter(reliability_density,h,'same');
figure;imagesc(reliability_density_smooth);colormap(hot);
line(repmat(length(bin_edges)/2,2,1),ylim,'color','w','linewidth',2)
line(xlim,repmat(length(bin_edges)/2,2,1),'color','w','linewidth',2)
xlabel('Early reliability difference');
ylabel('Late reliability difference')

lim_center = max(abs(length(bin_centers)/2 - [find(any(reliability_density_smooth,1),1); ...
    find(any(reliability_density_smooth,2),1);
    find(any(reliability_density_smooth,1),1,'last'); ...
    find(any(reliability_density_smooth,2),1,'last')]));
lim = [floor(length(bin_centers)/2-lim_center), ceil(length(bin_centers)/2+lim_center)];
xlim(lim); ylim(lim);
set(gca,'XTickLabel',bin_label(str2num(get(gca,'XTickLabel'))));
set(gca,'YTickLabel',bin_label(str2num(get(gca,'YTickLabel'))));
set(gca,'YDir','normal');



%% Slide 11) Make movie of reliability changes / density heatmap (as above)

savename = '3session_ngap_early_bothsig_comparison_regline';
writerObj = VideoWriter(['/usr/local/lab/People/Andy/Corticospinal/dated_analyses/150327/after_retreat/reliability_density_' savename '.avi']);
writerObj.FrameRate = 2;
open(writerObj);
% Group sessions, look at changes
for i = 1:8
    
    session_grp = {[1:3] [4:6]+i-1};
    premove_frames = 1:20;
    move_frames = 31:51;
    timing_reliability = cell(length(animals),length(session_grp));
    comparison_chi_sq = cell(length(animals),length(session_grp));
    session_comparison_chi_sq = cell(length(animals),1);
    for curr_animal = 1:length(analysis)
        
        use_sessions = find(cellfun(@(x) ~isempty(x), {analysis(curr_animal).im(:).move_onset_aligned}));
        curr_sessions_grp = cellfun(@(x) intersect(use_sessions,x),session_grp,'uni',false);
        if any(cellfun(@isempty,curr_sessions_grp))
            continue
        end
        
        timing_reliability(curr_animal,:) = cellfun(@(x) [ ...
            permute(nanmean(any(x(:,premove_frames,:) > 0,2),1),[3 2 1]), ...
            permute(nanmean(any(x(:,move_frames,:) > 0,2),1),[3 2 1])], ...
            cellfun(@(x) vertcat(analysis(curr_animal).im(x).move_onset_aligned), ...
            curr_sessions_grp,'uni',false),'uni',false);
        
        comparison_chi_sq(curr_animal,:) = cellfun(@(x) AP_chisquare_matrix( ...
            permute(any(x(:,premove_frames,:) > 0,2),[1 3 2]), ...
            permute(any(x(:,move_frames,:) > 0,2),[1 3 2])), ...
            cellfun(@(x) vertcat(analysis(curr_animal).im(x).move_onset_aligned), ...
            curr_sessions_grp,'uni',false),'uni',false);
        
        session_comparison_chi_sq{curr_animal} = cellfun(@(x) AP_chisquare_matrix( ...
            permute(any(x(:,premove_frames,:) > 0,2),[1 3 2]), ...
            permute(any(x(:,move_frames,:) > 0,2),[1 3 2])), ...
            {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false);
        
    end
    
    reliability_diff = cell2mat(cellfun(@(x) diff(x,[],2),timing_reliability,'uni',false));
    sig_rois = cell2mat(cellfun(@(x) x' < 0.05,comparison_chi_sq,'uni',false));

    use_animals = cellfun(@(x) ~isempty(x),session_comparison_chi_sq);
    ever_sig_rois = cell2mat(cellfun(@(x) any(vertcat(x{:}) < 0.05,1)', ...
        session_comparison_chi_sq(use_animals),'uni',false));
    
    % Plot density plot of early v. late differences
    bin_edges = linspace(-1,1,200);
    bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
    bin_label = round(bin_centers*100)/100;
    bin_edges(1) = -Inf;
    bin_edges(end) = Inf;
    reliability_density = hist3(fliplr(reliability_diff(all(sig_rois,2),:)),'Edges',repmat({bin_edges},2,1));
    %reliability_density = hist3(fliplr(reliability_diff(ever_sig_rois,:)),'Edges',repmat({bin_edges},2,1));
    h = fspecial('gaussian',10,2);
    reliability_density_smooth = imfilter(reliability_density,h,'same');
    figure;imagesc(reliability_density_smooth);colormap(hot);
    line(repmat(length(bin_edges)/2,2,1),ylim,'color','w','linewidth',2)
    line(xlim,repmat(length(bin_edges)/2,2,1),'color','w','linewidth',2)
    xlabel(['Early reliability difference - ' num2str(session_grp{1})]);
    ylabel(['Late reliability difference - ' num2str(session_grp{2})]);
    
    % lim_center = max(abs(length(bin_centers)/2 - [find(any(reliability_density_smooth,1),1); ...
    %     find(any(reliability_density_smooth,2),1);
    %     find(any(reliability_density_smooth,1),1,'last'); ...
    %     find(any(reliability_density_smooth,2),1,'last')]));
    % lim = [floor(length(bin_centers)/2-lim_center), ceil(length(bin_centers)/2+lim_center)];
    
    lim = [40 160];
    xlim(lim); ylim(lim);
    set(gca,'XTickLabel',bin_label(str2num(get(gca,'XTickLabel'))));
    set(gca,'YTickLabel',bin_label(str2num(get(gca,'YTickLabel'))));
    set(gca,'YDir','normal');
    
    % Fit points with line
    [a,b] = polyfit(reliability_diff(any(sig_rois,2),1),reliability_diff(any(sig_rois,2),2),1);
    y = (bin_centers(xlim)*a(1) + a(2))*(diff(ylim)/diff(bin_centers(ylim)))+mean(ylim);
    line(xlim,y,'color','w','linewidth',2,'linestyle','--')
     
    F = getframe(gcf);
    writeVideo(writerObj,F);
    
end

close(writerObj);
disp('done')


%% Slide 12) Get move/premove significant cells by session

premove_frames = 1:20;
move_frames = 31:61;
reliability_diff_session = cell(length(animals),1);
comparison_chi_sq_session = cell(length(animals),1);
for curr_animal = 1:length(analysis)
    
    use_sessions = find(cellfun(@(x) ~isempty(x), {analysis(curr_animal).im(:).move_onset_aligned}));
    n_rois = length(d_grp_all{curr_animal});
    
    reliability_diff_session{curr_animal} = nan(n_rois,14);
    comparison_chi_sq_session{curr_animal} = nan(n_rois,14);
    
    reliability_diff_session{curr_animal}(:,use_sessions) = cell2mat(cellfun(@(x) diff([ ...
        permute(nanmean(any(x(:,premove_frames,:) > 0,2),1),[3 2 1]), ...
        permute(nanmean(any(x(:,move_frames,:) > 0,2),1),[3 2 1])],[],2), ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false));

    comparison_chi_sq_session{curr_animal}(:,use_sessions) = cell2mat(cellfun(@(x) AP_chisquare_matrix( ...
        permute(any(x(:,premove_frames,:) > 0,2),[1 3 2]), ...
        permute(any(x(:,move_frames,:) > 0,2),[1 3 2]))', ...
        {analysis(curr_animal).im(use_sessions).move_onset_aligned},'uni',false));
    
end

% Plot relative fraction of significant cells
all_mean_sig = cell2mat(cellfun(@(x) nanmean(x < 0.05,1),comparison_chi_sq_session,'uni',false));
all_mean_movesig = cell2mat(cellfun(@(x,y) nanmean(x < 0.05 & y > 0,1), ...
    comparison_chi_sq_session,reliability_diff_session,'uni',false));
all_mean_premovesig = cell2mat(cellfun(@(x,y) nanmean(x < 0.05 & y < 0,1), ...
    comparison_chi_sq_session,reliability_diff_session,'uni',false));

figure; hold on;
errorbar(nanmean(all_mean_sig),nanstd(all_mean_sig)./sqrt(sum(~isnan(all_mean_sig))),'k','linewidth',2)
errorbar(nanmean(all_mean_movesig),nanstd(all_mean_movesig)./sqrt(sum(~isnan(all_mean_movesig))),'g','linewidth',2)
errorbar(nanmean(all_mean_premovesig),nanstd(all_mean_premovesig)./sqrt(sum(~isnan(all_mean_premovesig))),'r','linewidth',2)

ylabel('Fraction of ROIs');
xlabel('Session');
legend({'All active' 'Movement' 'Pre-movement'})
 
 
% Get ROIs which switch between down and up (I doubt these are real)
double_sig_cells = cellfun(@(x,y) any(x < 0.05 & y > 0,2) & ...
    any(x < 0.5 & y < 0,2), comparison_chi_sq_session, ...
    reliability_diff_session,'uni',false);


%% Find significantly modulated ROIs the old way (more robust?)

% NOTE: First run cell to find dendrite groups
classified_rois = struct('movement',cell(length(animals),1),'quiescent',cell(length(animals),1));
for curr_animal = 1:length(animals);
    
    % clean workspace
    clearvars -except animals d_grp_all classified_rois curr_animal
    
    animal = animals{curr_animal};
 
    % Use this if data not loaded in
    data_dir = '/usr/local/lab/People/Andy/Data';
    curr_analysis_dir = [data_dir filesep animal filesep ...
        animal '_batch_thresh_roi'];
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    n_sessions = length(curr_analysis_files);
    
    classified_rois(curr_animal).movement = nan(length(d_grp_all{curr_animal}),n_sessions);
    classified_rois(curr_animal).quiescent = nan(length(d_grp_all{curr_animal}),n_sessions);
    
    for curr_session = 1:n_sessions
        curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
        curr_data = load(curr_file);
        
        % Threshold data
        event_thresh = 3;
        if exist('d_grp_all','var')
            % If ROIs grouped, then average together and threshold
            curr_df_grp = cell2mat(cellfun(@(x) nanmean(curr_data.im.roi_trace_df(x,:),1), ...
                d_grp_all{curr_animal},'uni',false));
            curr_df_thresh = AP_caEvents_thresh(curr_df_grp,event_thresh);
        else
            % If not, just threshold all
            curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
        end
        
        % Get periods of lever movement/quiescence
        [lever_active,lever_force_resample] = ...
            AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        % Get movement/quiescence for each frame
        % (in at least one instance, the number of imaged frames actually
        % exceeds the number of frame triggers? - AP132 s11)
        num_frames = min(size(curr_data.im.roi_trace_df,2),length(curr_data.bhv.frame_times));
        lever_active_frames = lever_active(round(curr_data.bhv.frame_times(1:num_frames)*1000));
        
        % Split active frames by movement blocks / quiescent frames
        % (to split up by movement blocks / quiescent blocks)
        %lever_active_frames_nan = +lever_active_frames;
        %lever_active_frames_nan(~lever_active_frames) = NaN;
        boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
        lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
        
        % Get activity during movement/quiescence
        movement_activity = +(curr_df_thresh(:,1:num_frames) > 0)*lever_active_frames;

        % Get shuffled activity distribution
        num_rep = 10000;
        shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
        lever_active_shuffle = nan(length(lever_active_frames),num_rep);
        for i = 1:num_rep
           lever_active_shuffle(:,i) = ...
               vertcat(lever_active_frames_split{shuffle_perms(:,i)});
        end
        
        shuffle_movement_activity = +(curr_df_thresh(:,1:num_frames) > 0)*lever_active_shuffle;

        movement_rank = tiedrank([movement_activity shuffle_movement_activity]')';
        movement_p = movement_rank(:,1)/(num_rep+1);
        movement_cells = movement_p > 0.975;
        quiescent_cells = movement_p < 0.025;
       
        classified_rois(curr_animal).movement(:,curr_session) = movement_cells;
        classified_rois(curr_animal).quiescent(:,curr_session) = quiescent_cells;
        
    end 
    
    disp(animal);
end





%% Get average activity per ROI during movement/quiescence

avg_activity = struct('movement',cell(length(animals),1),'quiescent',cell(length(animals),1));
for curr_animal = 1:length(animals);
    
    animal = animals{curr_animal};
 
    % Use this if data not loaded in
    data_dir = '/usr/local/lab/People/Andy/Data';
    curr_analysis_dir = [data_dir filesep animal filesep ...
        animal '_batch_thresh_roi'];
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    n_sessions = length(curr_analysis_files);
    
    avg_activity(curr_animal).movement = nan(length(d_grp_all{curr_animal}),n_sessions);
    avg_activity(curr_animal).quiescent = nan(length(d_grp_all{curr_animal}),n_sessions);
    
    for curr_session = 1:n_sessions
        curr_file = [curr_analysis_dir filesep curr_analysis_files(curr_session).name];
        curr_data = load(curr_file);
        
        % Threshold data
        event_thresh = 3;
        if exist('d_grp_all','var')
            % If ROIs grouped, then average together and threshold
            curr_df_grp = cell2mat(cellfun(@(x) nanmean(curr_data.im.roi_trace_df(x,:),1), ...
                d_grp_all{curr_animal},'uni',false));
            curr_df_thresh = AP_caEvents_thresh(curr_df_grp,event_thresh);
        else
            % If not, just threshold all
            curr_df_thresh = AP_caEvents_thresh(curr_data.im.roi_trace_df,event_thresh);
        end
        
        % Get periods of lever movement/quiescence
        [lever_active,lever_force_resample] = ...
            AP_parseLeverMovement_updated(curr_data.bhv.lever_force);
        
        % Get movement/quiescence for each frame
        % (in at least one instance, the number of imaged frames actually
        % exceeds the number of frame triggers? - AP132 s11)
        num_frames = min(size(curr_data.im.roi_trace_df,2),length(curr_data.bhv.frame_times));
        lever_active_frames = lever_active(round(curr_data.bhv.frame_times(1:num_frames)*1000));
        
        % Get activity during movement/quiescence
        movement_activity = +(curr_df_thresh(:,1:num_frames) > 0)*lever_active_frames;
        quiescent_activity = +(curr_df_thresh(:,1:num_frames) > 0)*~lever_active_frames;
        
        avg_activity(curr_animal).movement(:,curr_session) = movement_activity/(sum(lever_active_frames));
        avg_activity(curr_animal).quiescent(:,curr_session) = quiescent_activity/(sum(~lever_active_frames));
        
    end
end

















