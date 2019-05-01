%% Classify cells by activity relating to behavior - movement epoch shuffle (peak+oopsi)
clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

for curr_animal = 1:length(animals)
    clearvars -except curr_animal animals
    
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
    
    classified_cells = {};
    
    for curr_day = 1:length(days)
        
        % Don't bother if AP71 day 1/2
        if strcmp(animal,'AP71') && (curr_day == 1 || curr_day == 2);
            continue
        end
        
        %%%% Load/Initialize
        clearvars -except animal days curr_day classified_cells animal animals
        
        day = num2str(days{curr_day});
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            error('Problem with behavior file')
        end
        raw_bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(raw_bhv);
        
        % load labels file
        label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
            animal '_roi_template'];
        labels_name = [animal '_roilabels.roilabel'];
        load([label_path filesep labels_name],'-MAT');
        
        % get cell labels
        cells = 1:size(im.roi_trace_df,1);
        gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
        pyr_cells = cells(~ismember(cells,gad_cells));
        contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels));
        
        %%%% Process       
        
        spread_frames = 5;
        spread_frames_filter = ones(1,spread_frames);
        
        % get binary trace lick, also get lick epochs
        lick_trace = zeros(1,size(im.roi_trace_df,2));
        lick_trace(round(bhv.lick_frames( ...
            ~isnan(bhv.lick_frames) & ...
            bhv.lick_frames<size(im.roi_trace_df,2) ...
            & bhv.lick_frames > 0))) = 1;
        lick_trace = +(conv(lick_trace,spread_frames_filter,'same') > 0);
        %         lick_stopstart = diff(lick_trace);
        %         lick_epochs = [0 find(lick_stopstart ~= 0) length(lick_trace)];
        %         lick_epoch_size = diff(lick_epochs);
        %         lick_epoch_shuffle = mat2cell(lick_trace',lick_epoch_size,1);
        lick_stopstart = diff([0 lick_trace 0]);
        lick_starts = find(lick_stopstart == 1);
        lick_stops = find(lick_stopstart == -1);
        lick_move_epochs = cellfun(@(x,y) ones(y-x,1), ...
            num2cell(lick_starts),num2cell(lick_stops),'uni',false);
        lick_epoch_shuffle = [lick_move_epochs ...
            num2cell(zeros(1,length(lick_trace) - ...
            sum(lick_stops-lick_starts)))];
        
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
        
        lever_quiescent = setdiff(1:size(im.roi_trace_df,2),lever_movement);
        quiescent_trace = zeros(1,size(im.roi_trace_df,2));
        quiescent_trace(lever_quiescent) = 1;
        
        % NOTE: if they're shuffled like this, tend to get bimodal
        % distributions of probability for cells that fire rarely
        
        % prepare movements to be shuffled by movement/quiecsent epochs
        %         lever_stopstart = diff(movement_trace);
        %         lever_epochs = [0 find(lever_stopstart ~= 0) length(movement_trace)];
        %         lever_epoch_size = diff(lever_epochs);
        
        %         lever_epoch_shuffle = mat2cell(movement_trace',lever_epoch_size,1);
        %         lever_velocity_epoch_shuffle = mat2cell(lever_velocity,lever_epoch_size,1);
        
        % this way of shufling probably doesn't make sense, but hopefully
        % better than the other way
        
        lever_stopstart = diff([0 movement_trace 0]);
        lever_starts = find(lever_stopstart == 1);
        lever_stops = find(lever_stopstart == -1);
        lever_move_epoch_lengths = lever_stops - lever_starts;
        lever_still_epoch_lengths = [lever_starts-[1 lever_stops(1:end-1)] ...
            1+length(movement_trace)-lever_stops(end)];
        
        % create cell array for splitting
        lever_still_epoch_split = cellfun(@(x) ones(1,x), ...
            num2cell(lever_still_epoch_lengths),'uni',false);
        
        lever_epoch_size = cell(1,length(lever_move_epoch_lengths) + ...
            length(lever_still_epoch_lengths));
        lever_epoch_size(1:2:end) = lever_still_epoch_split;
        lever_epoch_size(2:2:end) = num2cell(lever_move_epoch_lengths);
        lever_epoch_size = horzcat(lever_epoch_size{:});
        
        lever_epoch_shuffle = mat2cell(movement_trace',lever_epoch_size,1);
        lever_vel_epoch_shuffle = mat2cell(lever_velocity,lever_epoch_size,1);
        
        % Make timing traces: from 0 to 1 for each epoch
        lever_epoch_timing = cellfun(@(x) ...
            [0+1/length(x):1/length(x):1],lever_epoch_shuffle,'uni',false);
        
        % Define activity: Two ways: peak detection, thresholded oopsi
        % get activity, spread out
        
        % If using oopsi: filter oopsi results by thresholded traces
        for activity_measure = 1:2
            
            switch activity_measure
                case 1
                    %[peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);
                    %
                    %peak_matrix = zeros(size(im.roi_trace_df));
                    %for i = pyr_cells
                    %    peak_matrix(i,peak_frames{i}) = 1;
                    %end
%                     %for gad cells, use threshold
%                     for i = gad_cells
%                         if any(isnan(im.roi_trace_df(i,:)));
%                             continue
%                         end
%                         temp_smooth = smooth(im.roi_trace_df(i,:),30,'loess')';
%                         noise_est = mean(abs(im.roi_trace_df(i,:) - temp_smooth));
%                         noise_thresh = noise_est*3;
%                         peak_matrix(i,(im.roi_trace_df(i,:) ...
%                             > noise_thresh)) = im.roi_trace_df(i, ...
%                             (im.roi_trace_df(i,:) > noise_thresh));
%                     end

%                     peak_matrix = AP_make_peak_matrix(im.roi_trace_df,pyr_cells,gad_cells);
%                     % binarize pyramidal peaks
%                     peak_matrix(pyr_cells,:) = +[peak_matrix(pyr_cells,:) > 0];
                      peak_matrix = AP_caEvents(im.roi_trace_df,pyr_cells,gad_cells);

                   
                case 2
                    peak_matrix = im.roi_concat_oopsi;
                    
                    for curr_cell = 1:size(peak_matrix,1)
                        if ~ismember(curr_cell,gad_cells)
                            threshfac = 3.5;%how many times sd
                        else
                            threshfac = 2;
                        end
                        abovethresh_concat=[];
                        belowzero = im.roi_trace_df(curr_cell,im.roi_trace_df(curr_cell,:)<0);
                        thresh = std([belowzero,-belowzero]);
                        thresh = thresh*threshfac;
                        abovethresh = im.roi_trace_df(curr_cell,:)>thresh;
                        % has to be above thresh twice within 5 frames
                        abovethresh2 = abovethresh + [abovethresh(2:end),0] ...
                            + [abovethresh(3:end),0,0] +...
                            [abovethresh(4:end),0,0,0] + [abovethresh(5:end),0,0,0,0];
                        abovethresh2 = abovethresh2>=2;
                        peak_matrix(curr_cell,~abovethresh2) = 0;
                    end
                    
            end
            % get the actual values for movement/lick/timing traces
            move_deconv = peak_matrix*movement_trace';
            still_deconv = peak_matrix*(1-movement_trace)';
            vel_deconv = peak_matrix*lever_velocity;
            move_timing = (peak_matrix*[lever_epoch_timing{:}]')...
                ./nansum(peak_matrix,2);
            lick_deconv = peak_matrix*lick_trace';
            
            % this is to find significantly modulated cells for all movement
            % Differently from above! Shuffle the movement/quiescent epochs
            num_rois = size(im.roi_trace_df,1);            
            num_rep = 1000;
            move_deconv_shuffle = zeros(num_rois,num_rep);
            perms = zeros(length(lever_epoch_shuffle),num_rep);
            for i = 1:num_rep;
                perms(:,i) = randperm(length(lever_epoch_shuffle));
            end;
            lever_epoch_perms = lever_epoch_shuffle(perms);
            lever_epoch_perms_mat = zeros(length(movement_trace),num_rep);
            % this is stupid, but matlab has no better cell2mat solution
            for i = 1:num_rep;
                lever_epoch_perms_mat(:,i) = cell2mat(lever_epoch_perms(:,i));
            end
            move_deconv_shuffle = peak_matrix*lever_epoch_perms_mat;
            
            move_deconv_prctile = prctile(move_deconv_shuffle',[0.5 99.5]);
            move_cells = move_deconv' > move_deconv_prctile(2,:);
            still_cells = move_deconv' < move_deconv_prctile(1,:);
            
            % this is to find significantly modulated cells for all movement
            % Differently from above! Shuffle the movement/quiescent epochs
            num_rois = size(im.roi_trace_df,1);
            num_rep = 1000;
            move_deconv_shuffle = zeros(num_rois,num_rep);
            perms = zeros(length(lever_vel_epoch_shuffle),num_rep);
            for i = 1:num_rep;
                perms(:,i) = randperm(length(lever_vel_epoch_shuffle));
            end;
            lever_epoch_perms = lever_vel_epoch_shuffle(perms);
            lever_epoch_perms_mat = zeros(length(movement_trace),num_rep);
            % this is stupid, but matlab has no better cell2mat solution
            for i = 1:num_rep;
                lever_epoch_perms_mat(:,i) = cell2mat(lever_epoch_perms(:,i));
            end
            move_deconv_shuffle = peak_matrix*lever_epoch_perms_mat;
            
            move_deconv_prctile = prctile(move_deconv_shuffle',[0.5 99.5]);
            move_vel_cells = vel_deconv' > move_deconv_prctile(2,:);
            still_vel_cells = vel_deconv' < move_deconv_prctile(1,:);
            
            % same as above, except for licking
            num_rois = size(im.roi_trace_df,1);
            still_deconv = peak_matrix*(1-lick_trace)';
            num_rep = 1000;
            lick_deconv_shuffle = zeros(num_rois,num_rep);
            perms = zeros(length(lick_epoch_shuffle),num_rep);
            for i = 1:num_rep;
                perms(:,i) = randperm(length(lick_epoch_shuffle));
            end;
            lick_epoch_perms = lick_epoch_shuffle(perms);
            lick_epoch_perms_mat = zeros(length(lick_trace),num_rep);
            % this is stupid, but matlab has no better cell2mat solution
            for i = 1:num_rep;
                lick_epoch_perms_mat(:,i) = cell2mat(lick_epoch_perms(:,i));
            end
            lick_deconv_shuffle = peak_matrix*lick_epoch_perms_mat;
            
            lick_deconv_prctile = prctile(lick_deconv_shuffle',99.5);
            lick_cells = lick_deconv' > lick_deconv_prctile;
            
            % the old way: this is slower for some reason? maybe because of
            % less matrix operations?
            %         for i = 1:num_rep
            %             tic
            %             curr_perm = randperm(length(lever_epoch_shuffle));
            %             curr_movetrace = vertcat(lever_epoch_shuffle{curr_perm});
            %             tic
            %             move_deconv_shuffle(:,i) = ...
            %                 peak_matrix*curr_movetrace;
            %         end
            
            
            switch activity_measure
                
                case 1
                    % save
                    classified_cells.move_cells_peak{curr_day} = move_cells;
                    classified_cells.still_cells_peak{curr_day} = still_cells;
                    classified_cells.move_deconv_peak{curr_day} = move_deconv/sum(movement_trace);
                    classified_cells.still_deconv_peak{curr_day} = still_deconv/sum(1-movement_trace);
                    
                    classified_cells.move_vel_cells_peak{curr_day} = move_vel_cells;
                    classified_cells.still_vel_cells_peak{curr_day} = still_vel_cells;
                    classified_cells.vel_deconv_peak{curr_day} = vel_deconv/sum(lever_velocity);
                    
                    classified_cells.move_timing_peak{curr_day} = move_timing;
                    classified_cells.lick_cells_peak{curr_day} = lick_cells;
                    classified_cells.lick_deconv_peak{curr_day} = lick_deconv/sum(lick_trace);
                    
                case 2
                    % save
                    classified_cells.move_cells_oopsi{curr_day} = move_cells;
                    classified_cells.still_cells_oopsi{curr_day} = still_cells;
                    classified_cells.move_deconv_oopsi{curr_day} = move_deconv/sum(movement_trace);
                    classified_cells.still_deconv_oopsi{curr_day} = still_deconv/sum(1-movement_trace);
                    
                    classified_cells.move_vel_cells_oopsi{curr_day} = move_vel_cells;
                    classified_cells.still_vel_cells_oopsi{curr_day} = still_vel_cells;
                    classified_cells.vel_deconv_oopsi{curr_day} = vel_deconv/sum(lever_velocity);
                    
                    classified_cells.move_timing_oopsi{curr_day} = move_timing;
                    classified_cells.lick_cells_oopsi{curr_day} = lick_cells;
                    classified_cells.lick_deconv_oopsi{curr_day} = lick_deconv/sum(lick_trace);
            end
            
            disp(['Finished activity measure: ' num2str(activity_measure)]);
            
        end
        
        disp(['Finished day: ' day]);
        
    end
    disp(['Finished all ' animal]);
    save_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification';
    save([save_path filesep animal ...
        '_classified_cells_caEvents_5'],'classified_cells');
    
end

%% Compare classified cells by method
clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
class_compare = cell(size(animals));
for curr_animal = 1:length(animals)
    clearvars -except curr_animal animals class_compare
    
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
    
    use_days = 1:length(classified_cells.move_cells_peak);
    if strcmp(animal,'AP71');
        use_days = 3:length(classified_cells.move_cells_peak);
    end
    
    class_compare{curr_animal} = zeros(length(use_days), ...
        length(curr_cells),3);
    
    class_compare{curr_animal}(:,:,1) = ...
        cell2mat(cellfun(@(x) x(curr_cells),...
        classified_cells.move_cells_peak(use_days),'uni',0)');
    %vertcat(classified_cells.move_cells{:});%.* ...
    %horzcat(classified_cells.move_deconv{:})';
    
    class_compare{curr_animal}(:,:,2) = ...
        cell2mat(cellfun(@(x) x(curr_cells),...
        classified_cells.move_cells_oopsi(use_days),'uni',0)');
    
    %     load([animal '_classified_cells_epoch_shuffle_oopsi_thresh.mat']);
    %     class_compare{curr_animal}(:,:,2) = ...
    %         cell2mat(cellfun(@(x) x(curr_cells),...
    %         classified_cells.move_cells(use_days),'uni',0)');
    %         %vertcat(classified_cells.move_cells{:});%.* ...
    %         %horzcat(classified_cells.move_deconv{:})';
    %
    % %     load([animal '_classified_cells_epoch_shuffle.mat']);
    % %     class_compare{curr_animal}(:,:,3) = ...
    % %         cell2mat(cellfun(@(x) x(curr_cells),...
    % %         classified_cells.move_cells(use_days),'uni',0)');
    % %         %vertcat(classified_cells.move_cells{:});%.* ...
    % %         %horzcat(classified_cells.move_deconv{:})';
    %
    %     load([animal '_classified_cells.mat']);
    %     class_compare{curr_animal}(:,:,3) = ...
    %         cell2mat(cellfun(@(x) x(curr_cells),...
    %         classified_cells.move_cells_vel(use_days),'uni',0)');
    %         %vertcat(classified_cells.move_cells{:});%.* ...
    %         %horzcat(classified_cells.move_deconv{:})';
end

%% Check cells: plot cells of classified type with lever movement

animal = 'AP76';
curr_day = 10; 
% class: move,still,move vel, still vel (vel = based on lever velocity
% instead of binary active/inactive portions)
curr_class = 'move';

load(['/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/' ...
    animal '_classified_cells_epoch_shuffle_oopsi.mat']);

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

day = num2str(days{curr_day});
analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
analysis_name = [day '_' animal '_bgROI_analysis.mat'];
load([analysis_path filesep analysis_name]);

% load labels file
label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
    animal '_roi_template'];
labels_name = [animal '_roilabels.roilabel'];
load([label_path filesep labels_name],'-MAT');

% get cell labels
cells = 1:size(im.roi_trace_df,1);
gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));
filled_cells = find(cellfun(@(x) any(strcmp('filled',x)),roi_labels));

switch curr_class
    case 'move'
        cell_class = classified_cells.move_cells{curr_day};
    case 'still'
        cell_class = classified_cells.still_cells{curr_day};
    case 'move vel'
        cell_class = classified_cells.move_cells_vel{curr_day};
    case 'still vel'
        cell_class = classified_cells.still_cells_vel{curr_day};   
end

cell_class_idx = setdiff(find(cell_class),filled_cells);

[peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df(cell_class_idx,:));

figure; 
p1= subplot(4,1,1);
plot(bhv.lever_force_resample,'k'); 
ylim([min(bhv.lever_force_resample),max(bhv.lever_force_resample)]);
p2 = subplot(4,1,2:4); hold on;
for i = 1:length(cell_class_idx)
    if ismember(cell_class_idx(i),pyr_cells)
        curr_color = 'k';
    else
        curr_color = 'b';
    end
    plot(im.roi_trace_df(cell_class_idx(i),:)+i,curr_color);
    if ~isempty(peak_frames{i});
        plot(peak_frames{i},i,'.m','MarkerSize',20);
    end
    plot(im.roi_concat_oopsi(cell_class_idx(i),:)+i,'r')
end

linkaxes([p1 p2],'x');

% peak_matrix = zeros(size(im.roi_trace_df));
% for i = 1:length(peak_frames)
%     peak_matrix(i,peak_frames{i}) = 1;
% end





%% Plot results from cell classification
classified_cells_fixed = classified_cells;
for i = 1:length(classified_cells.cue_cells)
    % find double-classified cells
    double_cells = intersect(classified_cells.cue_cells{i}, ...
        classified_cells.reward_cells{i});
    % for each double cell, find which bhv is used more
    for curr_cell = double_cells
        more_cue = sum(classified_cells.cued_movement_aligned_active{i}(curr_cell,:)>0) ...
            > sum(classified_cells.reward_aligned_active{i}(curr_cell,:)>0);
        if more_cue
            classified_cells_fixed.reward_cells{i}( ...
                classified_cells_fixed.reward_cells{i} == curr_cell) = [];
        else
            classified_cells_fixed.cue_cells{i}( ...
                classified_cells_fixed.cue_cells{i} == curr_cell) = [];
        end
    end
end

% plot number of behavior-related cells
figure; hold on; set(gcf,'Name',animal)
plot(cellfun(@length,classified_cells_fixed.cue_cells),'k');
plot(cellfun(@length,classified_cells_fixed.reward_cells),'r');
title('Number of cells')

% plot robustness, daily and stable cells
figure; hold on; set(gcf,'Name',animal)
cue_robustness_mean = ...
    cellfun(@(x,y) mean(mean(x(y,:) > 0,2)),classified_cells.cued_movement_aligned_active, ...
    classified_cells_fixed.cue_cells);
cue_robustness_sem = ...
    cellfun(@(x,y) std(mean(x(y,:) > 0,2))/sqrt(length(y)),classified_cells.cued_movement_aligned_active, ...
    classified_cells_fixed.cue_cells);
reward_robustness_mean = ...
    cellfun(@(x,y) mean(mean(x(y,:) > 0,2)),classified_cells.reward_aligned_active, ...
    classified_cells_fixed.reward_cells);
reward_robustness_sem = ...
    cellfun(@(x,y) std(mean(x(y,:) > 0,2))/sqrt(length(y)),classified_cells.cued_movement_aligned_active, ...
    classified_cells_fixed.cue_cells);
subplot(2,1,1); hold on;
errorbar(cue_robustness_mean,cue_robustness_sem,'k');
errorbar(reward_robustness_mean,reward_robustness_sem,'r');
title('Robustness of response (daily cells)')

cue_cell_matrix = [];
reward_cell_matrix = [];
for i = 1:length(classified_cells_fixed.cue_cells)
    cue_cell_matrix(classified_cells_fixed.cue_cells{i},i) = 1;
    reward_cell_matrix(classified_cells_fixed.reward_cells{i},i) = 1;
end
stable_cue_cells = find(mean(cue_cell_matrix,2) > 0.33);
stable_reward_cells = find(mean(reward_cell_matrix,2) > 0.33);

cue_stable_robustness_mean = ...
    cellfun(@(x) mean(mean(x(stable_cue_cells,:) > 0,2)), ...
    classified_cells.cued_movement_aligned_active);
cue_stable_robustness_sem = ...
    cellfun(@(x) std(mean(x(stable_cue_cells,:) > 0,2))/sqrt(length(stable_cue_cells)), ...
    classified_cells.cued_movement_aligned_active);
reward_stable_robustness_mean = ...
    cellfun(@(x) mean(mean(x(stable_reward_cells,:) > 0,2)), ...
    classified_cells.reward_aligned_active);
reward_stable_robustness_sem = ...
    cellfun(@(x) std(mean(x(stable_reward_cells,:) > 0,2))/sqrt(length(stable_reward_cells)), ...
    classified_cells.reward_aligned_active);
subplot(2,1,2); hold on;
errorbar(cue_stable_robustness_mean,cue_stable_robustness_sem,'k');
errorbar(reward_stable_robustness_mean,reward_stable_robustness_sem,'r');
title('Robustness of response (stable cells)')

% "plot noise correlations"
figure; hold on; set(gcf,'Name',animal)
cue_noise_corr = cellfun(@(x,y) corrcoef(x(y,:)'), ...
    classified_cells.cued_movement_aligned_active, ...
    classified_cells_fixed.cue_cells,'UniformOutput',false);
reward_noise_corr = cellfun(@(x,y) corrcoef(x(y,:)'), ...
    classified_cells.reward_aligned_active, ...
    classified_cells_fixed.reward_cells,'UniformOutput',false);

cue_noise_corr_mean = cellfun(@(x) nanmean(x(:)), cue_noise_corr);
cue_noise_corr_sem = cellfun(@(x) nanstd(x(:))/sqrt(length(x)), cue_noise_corr);
reward_noise_corr_mean = cellfun(@(x) nanmean(x(:)), reward_noise_corr);
reward_noise_corr_sem = cellfun(@(x) nanstd(x(:))/sqrt(length(x)), reward_noise_corr);

subplot(2,1,1); hold on
errorbar(cue_noise_corr_mean,cue_noise_corr_sem,'k')
errorbar(reward_noise_corr_mean,reward_noise_corr_sem,'r')
title('Noise correlations (daily)')

cue_stable_noise_corr = cellfun(@(x) corrcoef(x(stable_cue_cells,:)'), ...
    classified_cells.cued_movement_aligned_active, ...
    'UniformOutput',false);
reward_stable_noise_corr = cellfun(@(x) corrcoef(x(stable_reward_cells,:)'), ...
    classified_cells.reward_aligned_active, ...
    'UniformOutput',false);

cue_stable_noise_corr_mean = cellfun(@(x) nanmean(x(:)), cue_stable_noise_corr);
cue_stable_noise_corr_sem = cellfun(@(x) nanstd(x(:))/sqrt(length(x)), cue_stable_noise_corr);
reward_stable_noise_corr_mean = cellfun(@(x) nanmean(x(:)), reward_stable_noise_corr);
reward_stable_noise_corr_sem = cellfun(@(x) nanstd(x(:))/sqrt(length(x)), reward_stable_noise_corr);

subplot(2,1,2); hold on;
errorbar(cue_stable_noise_corr_mean,cue_stable_noise_corr_sem,'k')
errorbar(reward_stable_noise_corr_mean,reward_stable_noise_corr_sem,'r')
title('Noise correlations (stable)')

% plot similarity of responsive cells
figure; hold on; set(gcf,'Name',animal)
same_cells_cue = nan(length(classified_cells_fixed.cue_cells),1);
same_cells_reward = nan(length(classified_cells_fixed.cue_cells),1);
diff_cells_cue = nan(length(classified_cells_fixed.cue_cells),1);
diff_cells_reward = nan(length(classified_cells_fixed.cue_cells),1);
for i = 2:length(classified_cells_fixed.cue_cells)
    same_cells_cue(i) = length(intersect(classified_cells_fixed.cue_cells{i}, ...
        classified_cells_fixed.cue_cells{i-1}))/length(classified_cells_fixed.cue_cells{i});
    same_cells_reward(i) = length(intersect(classified_cells_fixed.reward_cells{i}, ...
        classified_cells_fixed.reward_cells{i-1}))/length(classified_cells_fixed.reward_cells{i});
    
    diff_cells_cue(i) = length(setdiff(classified_cells_fixed.cue_cells{i}, ...
        classified_cells_fixed.cue_cells{i-1}))/length(classified_cells_fixed.cue_cells{i});
    diff_cells_reward(i) = length(setdiff(classified_cells_fixed.reward_cells{i}, ...
        classified_cells_fixed.reward_cells{i-1}))/length(classified_cells_fixed.reward_cells{i});
end
plot(same_cells_cue,'k');
plot(same_cells_reward,'r');
title('Cell similarity')

% plot number of cued-rewarded movements
figure; set(gcf,'Name',animal)
plot(cellfun(@(x) size(x,2),classified_cells.cued_movement_aligned_active),'k');
title('Number of cued-rewarded presses')


% plot number of cells classified as opposite on the day before
reward_crossclass = [];
cue_crossclass = [];
for i = 2:length(classified_cells_fixed.cue_cells)
   reward_crossclass(i) = length(intersect(classified_cells_fixed.reward_cells{i}, ...
        classified_cells_fixed.cue_cells{i-1}));
   cue_crossclass(i) = length(intersect(classified_cells_fixed.cue_cells{i}, ...
        classified_cells_fixed.reward_cells{i-1}));
end
figure; hold on
plot(cue_crossclass,'k')
plot(reward_crossclass,'r')
title('Cells which cross-classify between days')

% plot how many "stable" cells were found on a given day
num_stable_cue_cells = cellfun(@(x) length(intersect(stable_cue_cells,x)), ...
    classified_cells_fixed.cue_cells);
num_stable_reward_cells = cellfun(@(x) length(intersect(stable_reward_cells,x)), ...
    classified_cells_fixed.reward_cells);
figure; hold on;
plot(num_stable_cue_cells,'k')
plot(num_stable_reward_cells,'r')
title('Number of stable cells classified on each day')

% plot robustness of stable cells which were also classified on that day
cue_dstable_robustness = ...
    cell2mat(cellfun(@(x,y) [mean(mean(x(intersect(stable_cue_cells,y),:) > 0,2)); ...
    std(mean(x(intersect(stable_cue_cells,y),:) > 0,2))/sqrt(length(intersect(stable_cue_cells,y)))], ...
    classified_cells.cued_movement_aligned_active,classified_cells_fixed.cue_cells,'UniformOutput',false));

reward_dstable_robustness = ...
    cell2mat(cellfun(@(x,y) [mean(mean(x(intersect(stable_reward_cells,y),:) > 0,2)); ...
    std(mean(x(intersect(stable_reward_cells,y),:) > 0,2))/sqrt(length(intersect(stable_reward_cells,y)))], ...
    classified_cells.reward_aligned_active,classified_cells_fixed.reward_cells,'UniformOutput',false));

figure; hold on;
errorbar(cue_dstable_robustness(1,:),cue_dstable_robustness(2,:),'k');
errorbar(reward_dstable_robustness(1,:),reward_dstable_robustness(2,:),'r');
title('Robustness of response (stable cells intersect daily)')




%% playground

% plot all animals together
clear all
cue_dstable_robustness_all = cell(6,20);
reward_dstable_robustness_all = cell(6,20);
for curr_animal = [1:6 8 9];
    animal = ['AP7' num2str(curr_animal)];
    load([animal '_classified_cells.mat']);
    
    classified_cells_fixed = classified_cells;
    for i = 1:length(classified_cells.cue_cells)
        % find double-classified cells
        double_cells = intersect(classified_cells.cue_cells{i}, ...
            classified_cells.reward_cells{i});
        % for each double cell, find which bhv is used more
        for curr_cell = double_cells
            more_cue = sum(classified_cells.cued_movement_aligned_active{i}(curr_cell,:)>0) ...
                > sum(classified_cells.reward_aligned_active{i}(curr_cell,:)>0);
            if more_cue
                classified_cells_fixed.reward_cells{i}( ...
                    classified_cells_fixed.reward_cells{i} == curr_cell) = [];
            else
                classified_cells_fixed.cue_cells{i}( ...
                    classified_cells_fixed.cue_cells{i} == curr_cell) = [];
            end
        end
    end
    cue_cell_matrix = [];
    reward_cell_matrix = [];
    for i = 1:length(classified_cells_fixed.cue_cells)
        cue_cell_matrix(classified_cells_fixed.cue_cells{i},i) = 1;
        reward_cell_matrix(classified_cells_fixed.reward_cells{i},i) = 1;
    end
    stable_cue_cells = find(mean(cue_cell_matrix,2) > 0.33);
    stable_reward_cells = find(mean(reward_cell_matrix,2) > 0.33);
      
    % intersect(daily, stable)
%     cue_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x,y) mean(x(intersect(stable_cue_cells,y),:) > 0,2), ...
%     classified_cells.cued_movement_aligned_active,classified_cells_fixed.cue_cells,'UniformOutput',false);
%     
%     reward_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x,y) mean(x(intersect(stable_reward_cells,y),:) > 0,2), ...
%     classified_cells.reward_aligned_active,classified_cells_fixed.reward_cells,'UniformOutput',false);
%     
    % daily
    cue_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
    cellfun(@(x,y) mean(x(y,:) > 0,2), ...
    classified_cells.cued_movement_aligned_active,classified_cells_fixed.cue_cells,'UniformOutput',false);
    
    reward_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
    cellfun(@(x,y) mean(x(y,:) > 0,2), ...
    classified_cells.reward_aligned_active,classified_cells_fixed.reward_cells,'UniformOutput',false);

%     % stable
%     cue_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x) mean(x(stable_cue_cells,:) > 0,2), ...
%     classified_cells.cued_movement_aligned_active,'UniformOutput',false);
%     
%     reward_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x) mean(x(stable_reward_cells,:) > 0,2), ...
%     classified_cells.reward_aligned_active,'UniformOutput',false);

    % number of classified cells
%     max_cue = max(cellfun(@length,classified_cells.cue_cells));
%     max_reward = max(cellfun(@length,classified_cells.reward_cells));
%     
%     cue_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x) 100*length(x)/size(classified_cells.reward_aligned_active{1},1), ...
%     classified_cells_fixed.cue_cells,'UniformOutput',false);
%     
%     reward_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) = ...
%     cellfun(@(x) 100*length(x)/size(classified_cells.reward_aligned_active{1},1), ...
%     classified_cells_fixed.reward_cells,'UniformOutput',false);


%     % noise correlation
%     
%     cue_dstable_temp = ...
%     cellfun(@(x,y) corrcoef(x(intersect(stable_cue_cells,y),:)'), ...
%     classified_cells.cued_movement_aligned_active,classified_cells_fixed.cue_cells, ...
%     'UniformOutput',false);
% 
%     cue_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) ...
%         = cellfun(@(x) x(tril(true(length(x)),-1)),cue_dstable_temp,'UniformOutput',false);
% 
%     reward_dstable_temp = ...
%     cellfun(@(x,y) corrcoef(x(intersect(stable_reward_cells,y),:)'), ...
%     classified_cells.reward_aligned_active,classified_cells_fixed.reward_cells, ...
%     'UniformOutput',false);
% 
%     reward_dstable_robustness_all(curr_animal,1:length(classified_cells.cue_cells)) ...
%         = cellfun(@(x) x(tril(true(length(x)),-1)),reward_dstable_temp,'UniformOutput',false);

end



plot_cue = [];
plot_reward = [];
plot_all = [];
plot_cue_ci = [];
plot_reward_ci = [];
for i = 1:14%size(cue_dstable_robustness_all,2)
   curr_cue_mat = cell2mat(cue_dstable_robustness_all(:,i));
   curr_reward_mat= cell2mat(reward_dstable_robustness_all(:,i));
   
   robustness_all = [curr_cue_mat;curr_reward_mat];
   
   plot_cue(i,:) = [nanmean(curr_cue_mat)];
   plot_reward(i,:) = [nanmean(curr_reward_mat)];  
   
   plot_cue_ci(i,:) = bootci(1000,@nanmean,curr_cue_mat);
   plot_reward_ci(i,:) = bootci(1000,@nanmean,curr_reward_mat);
   
   
   %plot_cue(i,:) = [nanmean(curr_cue_mat) nanstd(curr_cue_mat)/sqrt(length(curr_cue_mat))];
   %plot_reward(i,:) = [nanmean(curr_reward_mat) nanstd(curr_reward_mat)/sqrt(length(curr_reward_mat))];    
   
   %plot_all(i,:) = [nanmean(robustness_all) nanstd(robustness_all)/sqrt(length(robustness_all))];
   
end
figure; hold on;
errorbar(1:length(plot_cue),plot_cue(:,1),plot_cue_ci(:,1),plot_cue_ci(:,2),'k','linewidth',2)
errorbar(1:length(plot_reward),plot_reward(:,1),plot_reward_ci(:,1),plot_reward_ci(:,2),'r','linewidth',2)
%errorbar(plot_all(:,1),plot_all(:,2),'k','linewidth',2)
xlim([0 sum(~isnan(plot_cue(:,1)))+1])
title('Combined animal cue/reward daily stable robustness')


%% The stuff above is old and probably crap - now it's not

%% Classify cells as rewarded-movement related (as opposed to quiescent)

clear all

animal = 'AP75';

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

classified_cells.rewarded_movement = cell(length(days),1);
classified_cells.movement = cell(length(days),1);
classified_cells.quiescent = cell(length(days),1);

for curr_day = 1:length(days)
   
    %%%% Load/Initialize
    clearvars -except animal days curr_day classified_cells
    
    day = num2str(days{curr_day});
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_name = [day '_' animal '_bgROI_analysis.mat'];
    load([analysis_path filesep analysis_name]);
      
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
          
    % get behavior file
    dir_currfolder = dir([data_path filesep day]);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
    
    % get xsg files
    xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
    xsg_dir = dir([xsg_folder filesep '*.xsg']);
    xsg_filenames = {xsg_dir.name};
    xsg_filenames = sort(xsg_filenames);
    xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
        'UniformOutput',false);
    
    % get behavior data and xsg sample rate
    try
        % ignore if something's wrong with datafile (usually, >1 of them)
        warning off
        load([data_path filesep day filesep bhv_filename],'-MAT');
        warning on
    catch me
        error('Problem with behavior file')
    end
    raw_bhv = saved_history.ProtocolsSection_parsed_events;
    num_trials = length(raw_bhv);
    
    % load labels file
    label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
        animal '_roi_template'];
    labels_name = [animal '_roilabels.roilabel'];
    load([label_path filesep labels_name],'-MAT');
    
    % get cell labels
    cells = 1:size(im.roi_trace_df,1);
    gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
    pyr_cells = cells(~ismember(cells,gad_cells));
    contaminated = find(cellfun(@(x) any(strcmp('contaminated',x)),roi_labels));
           
    % Break up baseline corrected trace into loops
    num_loops = length(im.roi_trace_split)/8;
    file_lengths = cellfun(@length,im.roi_trace_split);
    file_lengths_loops = reshape(file_lengths,8,num_loops);
    loop_lengths = sum(file_lengths_loops);
    roi_trace_df_split = mat2cell(im.roi_trace_df,size(im.roi_trace_df,1), ...
        loop_lengths);    

    xsg_sample = load(xsg_fullfilenames{1},'-MAT');
    xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
    
    rewarded_trials = find(cellfun(@(x) ~isempty(x.states.reward), bhv.bhv_frames));
    
    rewarded_movement_frames = {};
    quiescent_frames = {};
    movement_frames = {};
    
    % go through all xsg files, find movement times
    for curr_xsg = 1:length(xsg_fullfilenames)
        xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
        
        curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
        % if trials are not consecutive, delete the offending trial
        consec_trials = [0;diff(curr_trial_list(:,2))];
        curr_trial_list(consec_trials ~= 1,:) = [];
        
        [lever_active lever_force_smooth] = AP_parseLeverMovement(xsg_data);
                
        % get cue/reward times in current loop
        cue_sample = nan(length(curr_trial_list(:,2)),1);
        reward_start_sample = nan(length(curr_trial_list(:,2)),1);
        reward_stop_sample = nan(length(curr_trial_list(:,2)),1);
        for curr_trial_idx = 1:length(curr_trial_list(:,2));
            curr_trial = curr_trial_list(curr_trial_idx,2);
            
            % ignore if unrewarded trial
            if ~ismember(curr_trial,rewarded_trials)
                continue
            end
            
            % get trial start in dispatcher
            curr_bhv_start = raw_bhv{curr_trial}.states.bitcode(1);
            % get cue and reward sample in xsg (downsampled to ms)
            curr_bhv_cue_sample_rel =  (raw_bhv{curr_trial}.states.cue(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            cue_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_cue_sample_rel);
            
            curr_bhv_reward_start_sample_rel =  (raw_bhv{curr_trial}.states.reward(1) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_start_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_start_sample_rel);
            
            curr_bhv_reward_stop_sample_rel =  (raw_bhv{curr_trial}.states.reward(2) - ...
                curr_bhv_start)*(xsg_samplerate/10);
            reward_stop_sample(curr_trial_idx) = ...
                round(curr_trial_list(curr_trial == curr_trial_list(:,2),1)...
                *(xsg_samplerate/10) + ...
                curr_bhv_reward_stop_sample_rel);
        end       
        
        % find rewarded movement epochs
        rewarded_movement_samples = cell(length(reward_start_sample),1);
        for i = 1:length(reward_start_sample)
            if isnan(cue_sample(i)) || isnan(reward_stop_sample(i)) || ...
                    reward_stop_sample(i) > length(lever_active);
                continue
            end
            curr_move_start = find(lever_active(1:reward_start_sample(i)) ...
                == 0,1,'last');
            curr_move_end = find(lever_active(reward_start_sample(i):end) ...
                == 0,1) + reward_start_sample(i);
            rewarded_movement_samples{i} = curr_move_start:curr_move_end;
        end
        
        % find movement/quiescent samples
        quiescent_samples = find(~lever_active);
        movement_samples = find(lever_active);
              
        % convert rewarded movement epochs and quiescent periods to frames
        rewarded_movement_frames{curr_xsg} = ...
            cellfun(@(x) unique(round((x/(1000)*bhv.framerate))),rewarded_movement_samples, ...
            'UniformOutput',false);        
        quiescent_frames{curr_xsg} = unique(floor((quiescent_samples/1000)*bhv.framerate));         
        movement_frames{curr_xsg} = unique(floor((movement_samples/1000)*bhv.framerate));
    end
    
    % prepare frames for concatenating: add on loop lengths
    loop_lengths_add = num2cell([0 cumsum(loop_lengths(1:end-1))]);
    
    rewarded_movement_frames_concat = cellfun(@(x,y) ...
        horzcat(x{:}) + y,rewarded_movement_frames,loop_lengths_add, ...
        'UniformOutput',false);
    rewarded_movement_frames_concat = horzcat(rewarded_movement_frames_concat{:});
    
    quiescent_frames_concat = cellfun(@(x,y) ...
        horzcat(x) + y,quiescent_frames,loop_lengths_add, ...
        'UniformOutput',false);
    quiescent_frames_concat = vertcat(quiescent_frames_concat{:});
    quiescent_frames_concat(quiescent_frames_concat ...
        > size(im.roi_trace_df,2) | quiescent_frames_concat < 1) = [];
    
    movement_frames_concat = cellfun(@(x,y) ...
        horzcat(x) + y,movement_frames,loop_lengths_add, ...
        'UniformOutput',false);
    movement_frames_concat = vertcat(movement_frames_concat{:});
    movement_frames_concat(movement_frames_concat ...
        > size(im.roi_trace_df,2) | movement_frames_concat < 1) = [];
    
    % get peaks in trace
    [peak_frames peak_amplitudes] = AP_trace_peak(im.roi_trace_df);
    peak_matrix = zeros(size(im.roi_trace_df));
    for i = 1:size(im.roi_trace_df,1)
        if ~isempty(peak_frames{i})
            peak_matrix(i,peak_frames{i}) = 1;%peak_amplitudes{i};      
        end
    end
    
    % Ignore anything within the first 10 frames of a loop
    loop_junctions = cumsum(loop_lengths);
    loop_start_frames = repmat(loop_junctions(1:end-1),10,1) + ...
        repmat([1:10]',1,length(loop_junctions)-1);
    peak_matrix(:,loop_start_frames) = 0;
    
    % get peaks for classification
    rewarded_movement_peaks = ...
        sum(peak_matrix(:,rewarded_movement_frames_concat),2);
    movement_peaks = ...
        sum(peak_matrix(:,movement_frames_concat),2);
    quiescent_peaks = ...
        sum(peak_matrix(:,quiescent_frames_concat),2);

    % DO THE CLASSIFICATION
    
    % rewarded movement: vs. shuffled points in quiescent times
    peak_matrix_quiescent = peak_matrix(:,quiescent_frames_concat);
    num_rewarded_movement_frames = ...
        length(rewarded_movement_frames_concat);
    num_shuff = 1000;
    quiescent_peaks_shuff = nan(size(peak_matrix,1),num_shuff);
    for shuff = 1:num_shuff
        curr_shuff = randi(size(peak_matrix_quiescent,2), ...
            [num_rewarded_movement_frames 1]);
        quiescent_peaks_shuff(:,shuff) = sum(peak_matrix_quiescent(:,curr_shuff),2);
    end
    
    quiescent_peaks_ci = prctile(quiescent_peaks_shuff',99);
    classified_cells.rewarded_movement_cells{curr_day} =  ...
        find(rewarded_movement_peaks > quiescent_peaks_ci');
    classified_cells.rewarded_movement_peaks{curr_day} = rewarded_movement_peaks;
    
    % NOTE: not sure if this makes sense because quiescent frames ~
    % movement frames, probably should bootstrap both and compare
    
    % movement: vs. shuffled points in quiescent times
    peak_matrix_quiescent = peak_matrix(:,quiescent_frames_concat);
    num_movement_frames = ...
        length(movement_frames_concat);
    num_shuff = 1000;
    quiescent_peaks_shuff = nan(size(peak_matrix,1),num_shuff);
    for shuff = 1:num_shuff
        curr_shuff = randi(size(peak_matrix_quiescent,2), ...
            [num_movement_frames 1]);
        quiescent_peaks_shuff(:,shuff) = sum(peak_matrix_quiescent(:,curr_shuff),2);
    end
    
    quiescent_peaks_ci = prctile(quiescent_peaks_shuff',99);
    classified_cells.movement_cells{curr_day} =  ...
        find(movement_peaks > quiescent_peaks_ci');
    classified_cells.movement_peaks{curr_day} = movement_peaks;
    
    % quiescent: vs. shuffled points in movement times
    peak_matrix_movement = peak_matrix(:,movement_frames_concat);
    num_quiescent_frames = ...
        length(quiescent_frames_concat);
    num_shuff = 1000;
    movement_peaks_shuff = nan(size(peak_matrix,1),num_shuff);
    for shuff = 1:num_shuff
        curr_shuff = randi(size(peak_matrix_movement,2), ...
            [num_quiescent_frames 1]);
        movement_peaks_shuff(:,shuff) = sum(peak_matrix_movement(:,curr_shuff),2);
    end
    
    movement_peaks_ci = prctile(movement_peaks_shuff',99);
    classified_cells.quiescent_cells{curr_day} =  ...
        find(quiescent_peaks > movement_peaks_ci');
    classified_cells.quiescent_peaks{curr_day} = quiescent_peaks;
    
    disp(['Finished ' num2str(curr_day)]);
end

% 
% 
% a = zeros(size(im.roi_trace_df,1),length(days));
% for i = 1:length(days)
%     a(classified_cells.rewarded_movement_cells{i},i) = 1;
% end
% 
% a = a(pyr_cells,:);
% 
% kidx = kmeans(a(:,:),20,'EmptyAction','drop');
% [temp sort_idx] = sort(kidx);
% figure;imagesc(a(sort_idx,:));colormap(gray);










%% FOR PAPER: classify cells using AP_caEvents for ALL movements

clear all
%animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};
animals = {'AP103'};
for curr_animal = 1:length(animals)
    clearvars -except curr_animal animals
    
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
    
    classified_cells = {};
    
    for curr_day = 1:length(days);
        
        % Don't bother if AP71 day 1/2
        if strcmp(animal,'AP71') && (curr_day == 1 || curr_day == 2);
            continue
        end
        
        %%%% Load/Initialize
        clearvars -except animal days curr_day classified_cells animal animals
        
        day = num2str(days{curr_day});
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP' repmat('0',1,4-length(animal(3:end))) animal(3:end)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
%         % get behavior data and xsg sample rate
%         try
%             % ignore if something's wrong with datafile (usually, >1 of them)
%             warning off
%             load([data_path filesep day filesep bhv_filename],'-MAT');
%             warning on
%         catch me
%             warning('No behavior file')
%         end
%         
        % load labels file
        label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
            animal '_roi_template'];
        labels_name = [animal '_roilabels.roilabel'];
        load([label_path filesep labels_name],'-MAT');
        
        % get cell labels
        cells = 1:length(roi_labels);
        gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
        pyr_cells = cells(~ismember(cells,gad_cells));
        
        %%%% Process
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        lever_active_full = cell(length(xsg_filenames),1);
        for curr_xsg = 1:length(xsg_fullfilenames)
            xsg_data = load(xsg_fullfilenames{curr_xsg},'-MAT');
            
%             curr_trial_list = AP_ReadBitCode(xsg_fullfilenames{curr_xsg});
%             % if trials are not consecutive, delete the offending trial
%             consec_trials = [0;diff(curr_trial_list(:,2))];
%             curr_trial_list(consec_trials ~= 1,:) = [];
            
            [lever_active lever_force_smooth lever_force_velocity] = ...
                AP_parseLeverMovement(xsg_data);
            
            %save the active lever portions in frames
            curr_lever_active_frames = round((find(lever_active)/1000) ...
                *bhv.framerate);
            lever_active_full{curr_xsg} = unique(curr_lever_active_frames( ...
                curr_lever_active_frames <= loop_size(curr_xsg) & ...
                curr_lever_active_frames > 0)) + loop_frames(curr_xsg);
        end
        
        lever_movement = vertcat(lever_active_full{:});
        movement_trace = zeros(1,size(im.roi_trace_df,2));
        movement_trace(lever_movement) = 1;
        
        % Add leeway to movements to catch pre-move and just post-move cells
        movement_leeway = 5; % in frames
        movement_leeway_filter = ones(1,2*movement_leeway+1);
        movement_trace = +(conv(movement_trace,movement_leeway_filter,'same') > 0);
                  
        % Find and parse conserved movement epochs 
        lever_stopstart = diff([0 movement_trace 0]);
        lever_starts = find(lever_stopstart == 1);
        lever_stops = find(lever_stopstart == -1);
        lever_move_epoch_lengths = lever_stops - lever_starts;
        lever_still_epoch_lengths = [lever_starts-[1 lever_stops(1:end-1)] ...
            1+length(movement_trace)-lever_stops(end)];
        
        % create cell array for splitting
        lever_still_epoch_split = cellfun(@(x) ones(1,x), ...
            num2cell(lever_still_epoch_lengths),'uni',false);
        
        lever_epoch_size = cell(1,length(lever_move_epoch_lengths) + ...
            length(lever_still_epoch_lengths));
        lever_epoch_size(1:2:end) = lever_still_epoch_split;
        lever_epoch_size(2:2:end) = num2cell(lever_move_epoch_lengths);
        lever_epoch_size = horzcat(lever_epoch_size{:});
        
        lever_epoch_shuffle = mat2cell(movement_trace',lever_epoch_size,1);
      
        % Find activity
        peak_matrix = AP_caEvents(im.roi_trace_df,pyr_cells,gad_cells);
        
        % get the actual values for movement/lick/timing traces
        move_deconv = peak_matrix*movement_trace';
        still_deconv = peak_matrix*(1-movement_trace)';
        

        % Shuffle conserved movements through time to find activity which
        % is significantly restricted to movement
        num_rois = size(im.roi_trace_df,1);
        num_rep = 10000;
        perms = shake(repmat(transpose(1:length(lever_epoch_shuffle)),...
            1,num_rep),1);        

        % fucking dumb - matlab cell2mat breaks on this matrix for unknown
        % reasons, so putting in loop
        lever_epoch_perms_mat = nan(length(movement_trace),num_rep);
        for i = 1:num_rep
           lever_epoch_perms_mat(:,i) = ...
               vertcat(lever_epoch_shuffle{perms(:,i)});
        end
        
        clear move_deconv_shuffle
        move_deconv_shuffle = peak_matrix*lever_epoch_perms_mat;

        move_rank = tiedrank([move_deconv_shuffle move_deconv]')';
        move_p = move_rank(:,end)/(num_rep+1);
        move_cells = move_p > 0.995;
        still_cells = move_p < 0.005;
  
        % save
        classified_cells.move_cells_peak{curr_day} = move_cells;
        classified_cells.still_cells_peak{curr_day} = still_cells;
        classified_cells.move_deconv_peak{curr_day} = move_deconv/sum(movement_trace);
        classified_cells.still_deconv_peak{curr_day} = still_deconv/sum(1-movement_trace);
        classified_cells.move_p{curr_day} = move_p;       
        
        disp(['Finished day: ' day]);
        
    end
    disp(['Finished all ' animal]);
    save_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/final_class';
    save([save_path filesep animal ...
        '_classified_cells_caEvents_10k'],'classified_cells');
end

%% FOR PAPER: classify cells using AP_caEvents for REWARDED movements

clear all
animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

for curr_animal = 1:length(animals)
    clearvars -except curr_animal animals
    
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
    
    classified_cells = {};
    
    for curr_day = 1:length(days)
        
        % Don't bother if AP71 day 1/2
        if strcmp(animal,'AP71') && (curr_day == 1 || curr_day == 2);
            continue
        end
        
        %%%% Load/Initialize
        clearvars -except animal days curr_day classified_cells animal animals
        
        day = num2str(days{curr_day});
        analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
        analysis_name = [day '_' animal '_bgROI_analysis.mat'];
        load([analysis_path filesep analysis_name]);
        
        data_path = ['/usr/local/lab/People/Andy/Data/' animal];
        
        
        % get behavior file
        dir_currfolder = dir([data_path filesep day]);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = cell2mat((dir_filenames(bhv_file == 1)));
        
        % get xsg files
        xsg_folder = [data_path filesep day filesep 'AP00' animal(3:4)];
        xsg_dir = dir([xsg_folder filesep '*.xsg']);
        xsg_filenames = {xsg_dir.name};
        xsg_filenames = sort(xsg_filenames);
        xsg_fullfilenames = cellfun(@(x) [xsg_folder filesep x],xsg_filenames, ...
            'UniformOutput',false);
        
        % get behavior data and xsg sample rate
        try
            % ignore if something's wrong with datafile (usually, >1 of them)
            warning off
            load([data_path filesep day filesep bhv_filename],'-MAT');
            warning on
        catch me
            error('Problem with behavior file')
        end
        raw_bhv = saved_history.ProtocolsSection_parsed_events;
        num_trials = length(raw_bhv);
        
        % load labels file
        label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
            animal '_roi_template'];
        labels_name = [animal '_roilabels.roilabel'];
        load([label_path filesep labels_name],'-MAT');
        
        % get cell labels
        cells = 1:size(im.roi_trace_df,1);
        gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
        pyr_cells = cells(~ismember(cells,gad_cells));
        
        %%%% Process
        
        % get binary traces for movement, quiescence
        xsg_sample = load(xsg_fullfilenames{1},'-MAT');
        xsg_samplerate = xsg_sample.header.acquirer.acquirer.sampleRate;
        
        loop_size_split = cellfun(@length,im.roi_trace_split);
        loop_size = sum(reshape(loop_size_split,[8 length(loop_size_split)/8]));
        loop_size_sum = cumsum(loop_size_split);
        loop_frames = [0;loop_size_sum(8:8:end-7)];
        
        lever_active_full = cell(length(xsg_filenames),1);
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
        end
        
        lever_movement = vertcat(lever_active_full{:});
        movement_trace = zeros(1,size(im.roi_trace_df,2));
        movement_trace(lever_movement) = 1;
        
        % Add leeway to movements to catch pre-move and just post-move cells
        movement_leeway = 5; % in frames
        movement_leeway_filter = ones(1,2*movement_leeway+1);
        movement_trace = +(conv(movement_trace,movement_leeway_filter,'same') > 0);
                  
        % Find and parse conserved movement epochs 
        lever_stopstart = diff([0 movement_trace 0]);
        lever_starts = find(lever_stopstart == 1);
        lever_stops = find(lever_stopstart == -1);
        lever_move_epoch_lengths = lever_stops - lever_starts;
        lever_still_epoch_lengths = [lever_starts-[1 lever_stops(1:end-1)] ...
            1+length(movement_trace)-lever_stops(end)];
        
        % create cell array for splitting
        lever_still_epoch_split = cellfun(@(x) ones(1,x), ...
            num2cell(lever_still_epoch_lengths),'uni',false);
        
        lever_epoch_size = cell(1,length(lever_move_epoch_lengths) + ...
            length(lever_still_epoch_lengths));
        lever_epoch_size(1:2:end) = lever_still_epoch_split;
        lever_epoch_size(2:2:end) = num2cell(lever_move_epoch_lengths);
        lever_epoch_size = horzcat(lever_epoch_size{:});
        
        lever_epoch_shuffle = mat2cell(movement_trace',lever_epoch_size,1);
        lever_epoch_shuffle_frames = mat2cell(transpose(1:length(movement_trace)), ...
            lever_epoch_size,1);
        
        % mark movements that don't overlap with reward
        nonrewarded_movements = cellfun(@(x,y) any(x) && ~any(ismember( ...
            y,round(bhv.reward_frames))),lever_epoch_shuffle, ...
            lever_epoch_shuffle_frames);
        
        lever_epoch_shuffle_use = lever_epoch_shuffle(~nonrewarded_movements);
        
        use_frames = vertcat( ...
            lever_epoch_shuffle_frames{~nonrewarded_movements});
      
        % Find activity
        peak_matrix = AP_caEvents(im.roi_trace_df,pyr_cells,gad_cells);
        peak_matrix_use = peak_matrix(:,use_frames);
        
        movement_trace_use = movement_trace(use_frames);
        
        % get the actual values for movement/lick/timing traces
        move_deconv = peak_matrix_use*movement_trace_use';
        still_deconv = peak_matrix_use*(1-movement_trace_use)';
        
        % Shuffle conserved movements through time to find activity which
        % is significantly restricted to movement
        num_rois = size(im.roi_trace_df,1);
        num_rep = 10000;
        perms = shake(repmat(transpose(1:length(lever_epoch_shuffle_use)),...
            1,num_rep),1); 
        
        % fucking dumb - matlab cell2mat breaks on this matrix for unknown
        % reasons, so putting in loop
        lever_epoch_perms_mat = nan(length(movement_trace_use),num_rep);
        for i = 1:num_rep
           lever_epoch_perms_mat(:,i) = ...
               vertcat(lever_epoch_shuffle_use{perms(:,i)});
        end
        
        clear move_deconv_shuffle
        move_deconv_shuffle = peak_matrix_use*lever_epoch_perms_mat;

        move_rank = tiedrank([move_deconv_shuffle move_deconv]')';
        move_p = move_rank(:,end)/(num_rep+1);
        move_cells = move_p > 0.995;
        still_cells = move_p < 0.005;
  
        % save
        classified_cells.move_cells_peak{curr_day} = move_cells;
        classified_cells.still_cells_peak{curr_day} = still_cells;
        classified_cells.move_deconv_peak{curr_day} = move_deconv/sum(movement_trace_use);
        classified_cells.still_deconv_peak{curr_day} = still_deconv/sum(1-movement_trace_use);
        classified_cells.move_p{curr_day} = move_p;       
        
        disp(['Finished day: ' day]);
        
    end
    disp(['Finished all ' animal]);
    save_path = '/usr/local/lab/People/Andy/Data/Analysis/bhv_classification/final_class';
    save([save_path filesep animal ...
        '_classified_cells_caEvents_rewardedmove_10k'],'classified_cells');
end
























