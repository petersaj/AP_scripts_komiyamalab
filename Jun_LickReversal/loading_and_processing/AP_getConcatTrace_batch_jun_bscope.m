function AP_getConcatTrace_batch_jun_bscope(animal,roi_path)
% AP_getConcatTrace_batch(animal,days,roi_path)
% days = the days of the experiment. This is separate from roi_path because
% occasionally ROIs won't be drawn, in that case - mark as missed day by
% having empty data

% set flags to overwrite (1) or not (0)
im_flag = 0;
bhv_flag = 0;

disp('SPECIFIC TO JUN')

% get roi files and corresponding days
roi_dir = dir(roi_path);
roi_dir_filenames = {roi_dir.name};
roi_files_indx = cellfun(@(x) ~isempty(strfind(x,'.roi')), roi_dir_filenames);
roi_filenames = roi_dir_filenames(roi_files_indx);
roi_filenames = sort(roi_filenames);

days = cellfun(@(x) x(1:6), roi_filenames,'uni',false);

for curr_day = 1:length(days);
    
    %% Set up savename, if it already exists, skip processing
    processed_flag = 0;
    savename = ['/usr/local/lab/People/Jun/Data/'  ...
        animal filesep animal '_processedData' ...
        filesep animal '_' days{curr_day} ...
        '_processed.mat'];
    if exist(savename,'file')
        disp(['Already processed: ' savename]);
        load(savename);
        processed_flag = 1;
    end
    
    % Within this day - check / only use where tiff matches xsg number
    tiff_path = ['/usr/local/lab/People/Jun/Data/' animal filesep ...
        animal '_roi' filesep animal '_' days{curr_day}];  
    dir_currfolder = dir([tiff_path filesep '*.tif']);
    tiff_filenames = sort({dir_currfolder.name});
    tiff_underscore_idx = strfind(tiff_filenames,'_');
    tiff_loops = cellfun(@(x,y) str2num(x(y(end-2)+1:y(end-1)-1)), ...
        tiff_filenames,tiff_underscore_idx)-3; %!! -3 because always starts at 4
   
    xsg_path = ['/usr/local/lab/People/Jun/Data/' animal filesep ...
        animal '_xsg' filesep animal '_' days{curr_day}];
    xsg_filename = dir([xsg_path filesep '*.xsg']);
    xsg_fullfilenames = sort(cellfun(@(x) [xsg_path filesep x], ...
        {xsg_filename.name},'uni',false));
    xsg_loops = cellfun(@(x) str2num(x(end-7:end-4)),xsg_fullfilenames);
    
    
    tiff_loops_use = ismember(tiff_loops,xsg_loops);
    xsg_loops_use = ismember(xsg_loops,tiff_loops);
    
    
    %% Get image data
    if ~processed_flag || (processed_flag && im_flag)
        % initialize the structure to save the data
        im = struct;
        
        tiff_path = ['/usr/local/lab/People/Jun/Data/' animal filesep ...
            animal '_roi' filesep animal '_' days{curr_day}];
        
        % assume all tiff files in directory are motion corrected data
        dir_currfolder = dir([tiff_path filesep '*.tif']);
        tiff_filenames = sort({dir_currfolder.name});
        disp(vertcat(tiff_filenames{:}));
        
        % append directory to tiff filename
        tiff_corrected_filenames = cellfun(@(x) [tiff_path filesep x], ...
            tiff_filenames,'UniformOutput',false);
        
        % load ROIs
        load([roi_path filesep roi_filenames{curr_day}],'-MAT')
        
        disp('Getting and concatenating traces');
        % convert tiff names to cell, in case there's only one file
        if ~iscell(tiff_corrected_filenames);
            tiff_corrected_filenames = {tiff_corrected_filenames};
        end
        
        % only use tiff loops with matching xsg loops
        tiff_corrected_filenames(~tiff_loops_use) = [];
        
        % get traces from current tiff file
        imageinfo=imfinfo(tiff_corrected_filenames{1},'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        % get framerate (pick first image file, read header
        disp('Getting framerate from first file')
        img_info = imfinfo(tiff_corrected_filenames{1});
        img_info = img_info(1).ImageDescription;
        [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
        framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x, ...
            'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
        im.framerate = str2num(img_value{framerate_indx});
        
        % create masks from ROIs
        for n_polygon = 1:length(polygon.ROI);
            if ~isfield(polygon,'autosort');
                temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
                cellMask(:,n_polygon) = temp_mask(:);
            elseif isfield(polygon,'autosort');
                temp_mask = polygon.autosort(:,:,n_polygon);
                cellMask(:,n_polygon) = temp_mask(:);
            end
        end
        
        % Create mask for background ROIs, if available
        if isfield(polygon,'bgROI')
            for n_polygon = 1:length(polygon.bgROI);
                cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
            end
        end
        
        % get the traces for the ROIs
        % concatenate files that are part of the same loop
        tiff_loops_unique = unique(tiff_loops(tiff_loops_use));
        num_loops = length(tiff_loops_unique);
        
        roi_trace_long_split = cell(num_loops,1);
        roi_trace_long_split_bg = cell(num_loops,1);
        
        for curr_loop = 1:length(tiff_loops_unique)
            loop_num = tiff_loops_unique(curr_loop);
            curr_tiff_loop_filenames = tiff_corrected_filenames( ...
                tiff_loops == loop_num);
            
            curr_loop_trace = cell(length(curr_tiff_loop_filenames),1);
            curr_loop_trace_bg = cell(length(curr_tiff_loop_filenames),1);
            for i = 1:length(curr_tiff_loop_filenames);
                img_filename = [curr_tiff_loop_filenames{i}];
                
                % get traces from current tiff file
                imageinfo=imfinfo(img_filename,'tiff');
                numframes=length(imageinfo);
                M=imageinfo(1).Width;
                N=imageinfo(1).Height;
                
                im_raw = [];
                im_raw = zeros(N*M,numframes);
                
                tic
                disp('Loading file....')
                for loadframe = 1:numframes
                    im_temp = [];
                    im_temp = imread(img_filename,'tiff',loadframe);
                    % image: columns are frames
                    im_raw(:,loadframe) = im_temp(:);
                end
                disp('Done')
                toc
                
                disp(['Getting activity from file ' num2str(i) '/' num2str(length(curr_tiff_loop_filenames))]);
                
                traceVal = zeros(size(cellMask,2),numframes);               
                traceVal_bg = zeros(size(cellMask,2),numframes);
                
                if ~isfield(polygon,'autosort');
                    
                    % treat zeros as missing values
                    %im_raw(im_raw == 0) = NaN;
                    % Changed this: the new scanimage has 0s (and negs)
                    
                    for n_polygon = 1:size(cellMask,2)
                        % sum of pixel brighness / number of pixels, avg fluor in ROI
                        traceVal(n_polygon,:) = nanmean(im_raw(cellMask(:,n_polygon),:));
                        traceVal_bg(n_polygon,:) = nanmean(im_raw(cellMask_bg(:,n_polygon),:));
                    end
                    
                elseif isfield(polygon,'autosort');
                    % weighted fluorescence
                    traceVal = (cellMask'*im)./repmat(sum(cellMask',2),1,size(im,2));
                end
                
                curr_loop_trace{i} = traceVal;
                curr_loop_trace_bg{i} = traceVal_bg;
                
            end
            
            % save current file's traces
            roi_trace_long_split{curr_loop} = horzcat(curr_loop_trace{:});
            roi_trace_long_split_bg{curr_loop} = horzcat(curr_loop_trace_bg{:});
                                   
            clear curr_loop_trace curr_loop_trace_bg
        end
        
        %% DF, not sure what's good here
        
        % NOT SURE YET WHETHER INDIVIDUAL DF OK, OR SHOULD DO ALL CONCAT
        disp('Calculating df/f: using concatenated trace');
        % there are spikes at the loops (where the single is bigger, then
        % smaller than the concat), but they're only ~1% df/f, I thought concat
        % would be better but now maybe single, because maybe water added
        % between loops, and maybe edge taken to be activity? faster to df the
        % concat trace, but don't want loop related artifacts.
        
        % I changed my mind here: I think the loop artifacts are pretty
        % minorly changed by concatenation, and probably should be ignored
        % later anyway. On the other hand: if lots of activity during a
        % loop and no baseline is found, then a loop can be set to NaN,
        % which is much worse
        
        %         % get df/f values for traces (in individual loops)
        %         roi_trace_df = cell(size(roi_trace_long_split_bg));
        %         for curr_loop = 1:length(roi_trace_long_split_bg);
        %             [curr_bg_df curr_bg_baseline] = ...
        %                 AP_baselineEstimation(roi_trace_long_split_bg{curr_loop},im.framerate);
        %             curr_bgsubtract = roi_trace_long_split{curr_loop} - ...
        %                 (roi_trace_long_split_bg{curr_loop} - curr_bg_baseline);
        %
        %             [roi_trace_df{curr_loop} curr_baseline] = ...
        %                 AP_baselineEstimation(curr_bgsubtract,im.framerate);
        %         end
        
        % get df/f values for traces (from concatenated trace)
        [curr_bg_df curr_bg_baseline] = ...
            AP_baselineEstimation(horzcat(roi_trace_long_split_bg{:}),im.framerate);
        curr_bgsubtract = horzcat(roi_trace_long_split{:}) - ...
            (horzcat(roi_trace_long_split_bg{:}) - curr_bg_baseline);
        
        [roi_trace_df curr_baseline] = ...
            AP_baselineEstimation(curr_bgsubtract,im.framerate);
        
        % split roi_trace_df into loop sizes
        loop_sizes = cellfun(@(x) size(x,2),roi_trace_long_split_bg);
        roi_trace_df_split = mat2cell(roi_trace_df,size(roi_trace_df,1),loop_sizes);
        
        % save data into image-related structure
        im.roi_trace_split = roi_trace_long_split;
        im.roi_trace_split_bg = roi_trace_long_split_bg;
        im.roi_trace_df = roi_trace_df_split;
        
        % clear to save memory
        im_raw = [];
        
    end
    
    %% Get behavior data
    
    if ~processed_flag || (processed_flag && bhv_flag)
        
        bhv = struct;
        
        % dispatcher behavior filenames
        bhv_path = ['/usr/local/lab/People/Jun/Data/' animal filesep ...
            animal '_bhv'];
        bhv_filename = dir([bhv_path filesep ...
            'data_@JunOdor_experimenter_' animal '_' days{curr_day} '*.mat']);
        if length(bhv_filename) ~= 1
            error([animal ' ' days{curr_day} ': behavior files ~ = 1']);
        end
        bhv_fullfilename = [bhv_path filesep bhv_filename.name];
        
        % get applied/rewarded odors from dispatcher
        warning off
        curr_bhv = load(bhv_fullfilename);
        warning on
        num_completed_trials = length(curr_bhv.saved_history.ProtocolsSection_parsed_events);
        
        bhv.applied_odor = vertcat(curr_bhv.saved_history.OdorChoiceSection_AppliedOdor{1:num_completed_trials});
        bhv.rewarded_odor = vertcat(curr_bhv.saved_history.OdorChoiceSection_RewardedOdor{1:num_completed_trials});
        
        % get relevant state time settings
        bhv.odor_before_lick_time = vertcat(curr_bhv.saved_history.TimesSection_odor_before_lick_time{1:num_completed_trials});
        bhv.answer_time = vertcat(curr_bhv.saved_history.TimesSection_answer_time{1:num_completed_trials});
        bhv.water_time = vertcat(curr_bhv.saved_history.TimesSection_water_time{1:num_completed_trials});
        bhv.iti_step = vertcat(curr_bhv.saved_history.TimesSection_iti_step{1:num_completed_trials});
        bhv.iti_jitter = vertcat(curr_bhv.saved_history.TimesSection_iti_jitter{1:num_completed_trials});
        bhv.extra_iti = vertcat(curr_bhv.saved_history.TimesSection_extra_iti{1:num_completed_trials});
        bhv.iti = vertcat(curr_bhv.saved_history.TimesSection_iti{1:num_completed_trials});
        bhv.session_type = curr_bhv.saved_history.SessionTypeSection_SessionType(1:num_completed_trials);
        
        % xsg filenames
        xsg_path = ['/usr/local/lab/People/Jun/Data/' animal filesep ...
            animal '_xsg' filesep animal '_' days{curr_day}];
        xsg_filename = dir([xsg_path filesep '*.xsg']);
        xsg_fullfilenames = sort(cellfun(@(x) [xsg_path filesep x], ...
            {xsg_filename.name},'uni',false));
        
        % only use xsg loops with matching tiff loops
        xsg_fullfilenames(~xsg_loops_use) = [];
        
        [bhv.bhv_frames, bhv.imaged_trials, bhv.xsg_trials] ...
            = AP_getBehavior_dispatcher(im.framerate, bhv_fullfilename, xsg_fullfilenames);
        
    end
    
    %% Save
   
    % save roi_trace_long and roi_trace_long_split in the ROI path
    save(savename,'im','bhv');

    
end