%% load in image file

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','off');

img_filename = [tiff_path tiff_filename];

% get traces from current tiff file
imageinfo=imfinfo(img_filename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

try
    matlabpool
catch me
end

im = zeros(N,M,numframes);
im_temp = zeros(N,M);
disp('Loading file....')
parfor loadframe = 1:numframes
    im_temp = double(imread(img_filename,'tiff',loadframe));
    im(:,:,loadframe) = im_temp;
end
disp('Done')

%% Get trace from only ROI by loading in whole movies and multithreading

[roi_file roi_path] = uigetfile('*','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
w = waitbar(0,'Getting and concatenating traces');
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % Create mask for selected ROI, get trace
    for n_polygon = 1:length(polygon.ROI);
        if ~isfield(polygon,'autosort');
            temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            cellMask(:,n_polygon) = temp_mask(:);
        elseif isfield(polygon,'autosort');
            temp_mask = polygon.autosort(:,:,n_polygon);
            cellMask(:,n_polygon) = temp_mask(:);
        end
    end
    
    im = [];
    im = zeros(N*M,numframes);
    
    try
        matlabpool
    catch me
    end
    
    tic
    disp('Loading file....')
    parfor loadframe = 1:numframes
        im_temp = double(imread(img_filename,'tiff',loadframe));
        % image: columns are frames
        im(:,loadframe) = im_temp(:);
    end
    disp('Done')
    toc
    
    roi_trace = zeros(length(polygon.ROI),numframes);
    w_curr = waitbar(0,'Getting activity from current file');
        traceVal = zeros(length(polygon.ROI),numframes);
        for n_polygon = 1:length(polygon.ROI)
            % sum of pixel brighness / number of pixels, avg fluor in ROI
            traceVal(n_polygon,:) = sum(im(cellMask(:,n_polygon),:)./sum(sum(cellMask(:,n_polygon))));
            waitbar(n_polygon/length(polygon.ROI),w_curr,'Getting activity from current file');
        end
        roi_trace = traceVal;
        close(w_curr);
        
        % create df/f for roi_trace
        roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        try
            norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two gaussians
        catch me
            try
                norm_fit_obj = gmdistribution.fit(cellTrace',1); % there is one gaussian
            catch me
                continue
            end
        end
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f
        cellTrace = cellTrace./baseline_firing - 1;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear roi_trace
    waitbar(i/length(tiff_filename),w,'Getting and concatenating traces');
end
close(w);

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);
% clear to save memory
im = [];


%% create df/f for roi_trace
roi_trace_df = zeros(size(roi_trace));
for i = 1:size(roi_trace,1)
    cellTrace = roi_trace(i,:);
    % don't bother if they're all NaNs
    if ~any(isfinite(cellTrace));
        continue
    end
    try
        norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
    catch me
        norm_fit_obj = gmdistribution.fit(cellTrace',1); % one hidden curve if error
    end
    dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
    x = min(cellTrace):dx:max(cellTrace);
    x = x(1:end-1); % it adds an extra, not sure why but can fix later
    norm_fit = pdf(norm_fit_obj,x');
    baseline_firing = x(norm_fit == max(norm_fit));
    
    % change to df/f make peakin relative to df
    cellTrace = (cellTrace-baseline_firing)./baseline_firing;
    roi_trace_df(i,:) = cellTrace;
end

%% Get distances, correlations by distance between rois

% get centers of all ROIs
roi_centers = zeros(size(roi_trace_long),1);
for i = 1:size(roi_trace_long,1);
    % get center x and y from polygon boundaries
    center_x = mean(polygon.ROI{i}(:,1));
    center_y = mean(polygon.ROI{i}(:,2))*4; % correct for resolution
    roi_centers(i,1) = center_x;
    roi_centers(i,2) = center_y;
end

% loop through ROIs, sort by distance and get corrcoefs
roi_corrcoef = zeros(size(roi_trace_long,1));
for j = 1:size(roi_trace_long,1);
    % for chosen roi, get sort distances and plot in order
    curr_roi = j;
    curr_center = roi_centers(curr_roi,:);
    center_diffs = roi_centers(:,1) - curr_center(1);
    center_diffs(:,2) = roi_centers(:,2) - curr_center(2);
    roi_dists = sqrt(center_diffs(:,1).^2 + center_diffs(:,2).^2);
    roi_dists(:,2) = [1:size(roi_dists,1)];
    roi_dists_sort = sortrows(roi_dists);
    
    
    roi_trace_sort = (roi_trace_long(roi_dists_sort(:,2),:));
    for i = 1:size(roi_trace_sort,1)
        curr_corrcoef = corrcoef(roi_trace_sort(1,:),roi_trace_sort(i,:));
        curr_roi_corrcoef(i) = curr_corrcoef(2);
    end
    
    roi_corrcoef(curr_roi,:) = curr_roi_corrcoef;
end
figure;
imagesc(roi_corrcoef);

%% Compare peaks of traces

minpeak = 0.1; % min increase (std) to be considered active
roi_peak = {};
for i = 1:size(roi_trace_df,1)
    [max_ca min_ca] = peakdet(roi_trace_df(i,:),minpeak);
    roi_peak{i} = max_ca;
end

%% go through entire concat trace, get all calcium transients

peakstd = [1.5]; % definitions for max trigger, in stds
roi_peak_max = {};
roi_peak_min = {};
num_rois = size(roi_trace_long,1);
for i = 1:num_rois
    % get estimate of baseline variance
    try
        norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',2); % there are two hidden curves
    catch me
        norm_fit_obj = gmdistribution.fit(roi_trace_long(i,:)',1); % there are two hidden curves
    end
    baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
    % I think this is std? not sure
    curr_trace_std = sqrt(norm_fit_obj.Sigma(baseline_curve));
    minpeak = curr_trace_std*peakstd;
    [max_ca min_ca] = ap_peakdet(roi_trace_long(i,:),minpeak);
    % eliminate peaks at file junctions because they're mistimed
    if ~isempty(max_ca) && ~isempty(min_ca)
        file_stitch = [];
        file_stitch = find(mod(max_ca(:,1)-1,600) == 0);
        max_ca(file_stitch,:) = [];
        roi_peak_max{i} = [max_ca];
        roi_peak_min{i} = [min_ca];
    end
    % eliminate peaks that are negative
    roi_peak_max{i}(roi_peak_max{i}(:,2) < 0,:) = [];
end

% don't know if this made sense, not using now
% peaks_all = [];
% for j = 1:num_rois;
%     curr_peaks = [];
%     if ~isempty(roi_peak_max{j}) && ~isempty(roi_peak_min{j});
%         % if starts on a high, start with second
%         if roi_peak_max{j}(1,1) - roi_peak_min{j}(1,1) < 0;
%             roi_peak_max{j} = roi_peak_max{j}(2:end,:);
%             if isempty(roi_peak_max{j});
%                 continue
%             end
%         end
%         % if ends on a low, cut out the last low
%         if roi_peak_min{j}(end,1) - roi_peak_max{j}(end,1) > 0;
%             roi_peak_min{j} = roi_peak_min{j}(1:end-1,:);
%             if isempty(roi_peak_min{j});
%                 continue
%             end
%         end
%         curr_peaks = roi_peak_max{j}(:,2) - roi_peak_min{j}(:,2);
%         peaks_all = [peaks_all; curr_peaks];
%     end
% end


%% Concatenatenate traces of .roi already made

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','on');
for i = 1:length(tiff_filename);
    
    load([tiff_path tiff_filename],'-MAT')
    % create df/f for roi_trace
    roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f make peakin relative to df
        cellTrace = cellTrace./baseline_firing;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear polygon roi_mask roi_trace
end

num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_acq = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_acq);

%% Concatenatenate traces from only ROI

auto_align_flag = 0; % do you want to auto align ROI for each file?
if auto_align_flag == 1;
    [align_file align_path] = uigetfile('*.tif','Pick image to align to');
    img_align = double(imread([align_path align_file],'tiff'));
end

[roi_file roi_path] = uigetfile('*','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
w = waitbar(0,'Getting and concatenating traces');
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    if auto_align_flag == 0;
        % Create mask for selected ROI, get trace
        for n_polygon = 1:length(polygon.ROI);
            if ~isfield(polygon,'autosort');
                cellMask(:,:,n_polygon) = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            elseif isfield(polygon,'autosort');
                cellMask(:,:,n_polygon) = polygon.autosort(:,:,n_polygon);
            end
        end
        
        % if auto align selected, make masks from aligned ROIs
    elseif auto_align_flag == 1;
        aligned_ROI = cell(size(polygon.ROI));
        % get current file average
        img_avg = zeros(N,M);
        w_avg = waitbar(0,'Creating average of current file');
        for frame = 1:numframes;
            img_temp(:,:) = double(imread(img_filename,'tiff',frame));
            img_avg(:,:) = img_avg + img_temp./numframes;
            waitbar(frame/numframes,w_avg,'Creating average of current file');
        end
        close(w_avg);
        
        % DON'T HARDCODE THIS IN THE FUTURE
        disp 'xy_ratio hardcoded as 4'
        xy_ratio = 4;
        
        % get estimate for dx dy
        cc = normxcorr2(img_align,img_avg);
        [max_cc, imax] = max(abs(cc(:)));
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        
        % maximize correlation through transformation
        starting = [xpeak ypeak 0]; % dx dy theta
        options = optimset('Algorithm','interior-point');
        lb = [-100 -100/xy_ratio -pi/2];
        ub = [100 100/xy_ratio pi/2];
        % estimate transform of orig -> curr
        display 'Finding best fit...'
        estimates = fmincon(@(x) AP_affineAlign(x,img_align,img_avg,xy_ratio),starting,[],[],[],[],lb,ub,[],options);
        display 'Done'
        % get estimates from best guess
        dx = estimates(1);
        dy = estimates(2);
        theta = estimates(3);
        
        % transform polygons according to image match
        for poly_align = 1:length(polygon.ROI)
            
            udata = [1 N] - (N+1)/2;% input image x
            vdata = [1 M] - (M+1)/2;% input image y
            
            xdata = udata;% output image x
            ydata = vdata;% output image y
            
            % define affine transform matrix
            A = ...
                [cos(theta) sin(theta) 0;
                -sin(theta) cos(theta) 0;
                dx          dy          1];
            
            tform_roi = maketform('affine',A);
            
            % transform polygon, fix for xy_ratio
            aligned_ROI{poly_align} = tformfwd(tform_roi,[polygon.ROI{poly_align}(:,1)...
                polygon.ROI{poly_align}(:,2).*xy_ratio]);%, ...
            % 'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
            % correct final ROIs for xy_ratio
            aligned_ROI{poly_align}(:,2) = aligned_ROI{poly_align}(:,2)./xy_ratio;
            
            cellMask(:,:,poly_align) = poly2mask(aligned_ROI{poly_align}(:,1)',aligned_ROI{poly_align}(:,2)',N,M);
        end
    end
    
    roi_trace = zeros(length(polygon.ROI),numframes);
    w_curr = waitbar(0,'Getting activity from current file');
    for frame=1:numframes
        im_s=double(imread([tiff_path tiff_filename{i}],frame));
        traceVal = [];
        for n_polygon = 1:length(polygon.ROI)
            normalactivityMap = [];
            % sum of pixel brighness / number of pixels, avg fluor in ROI
            normactivityMap = (im_s.*cellMask(:,:,n_polygon))./sum(sum(cellMask(:,:,n_polygon)));
            traceVal(n_polygon,1) = sum(normactivityMap(:));
        end
        roi_trace(:,frame) = traceVal;
        waitbar(frame/numframes,w_curr,'Getting activity from current file');
    end
    close(w_curr);
    
    % create df/f for roi_trace
    roi_trace_df = zeros(size(roi_trace));
    for curr_roi = 1:size(roi_trace,1)
        cellTrace = roi_trace(curr_roi,:);
        if ~any(isfinite(cellTrace));
            continue
        end
        try
            norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two gaussians
        catch me
            try
                norm_fit_obj = gmdistribution.fit(cellTrace',1); % there is one gaussian
            catch me
                continue
            end
        end
        dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
        x = min(cellTrace):dx:max(cellTrace);
        x = x(1:end-1); % it adds an extra, not sure why but can fix later
        norm_fit = pdf(norm_fit_obj,x');
        baseline_firing = x(norm_fit == max(norm_fit));
        
        % change to df/f
        cellTrace = cellTrace./baseline_firing - 1;
        roi_trace_df(curr_roi,:) = cellTrace;
    end
    
    roi_trace_all(:,:,i) = roi_trace_df;
    clear roi_trace
    waitbar(i/length(tiff_filename),w,'Getting and concatenating traces');
    
end
close(w);

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);

%% Get and concatenate events
n_channels = 1; % true at the moment for all motion corrected files
n_frames = 1000; %get this later

[bhv_filename bhv_path] = uigetfile('*','Behavior file');
bhv = load([bhv_path bhv_filename],'-MAT');

saved_history = [];
saved_history = bhv.saved_history; % this in for now, somewhere below fix

[img_filename_all img_path]= uigetfile('*.tif','Choose TIFF files','Multiselect','on');
[acq_filename_all acq_path] = uigetfile('*.xsg','Choose XSG files','Multiselect','on');

w = waitbar(0, 'Getting and concatenating events');

left_down_all = [];
left_up_all = [];
reward_all = [];
unrewarded_all = [];
cue_onset_all = [];
frames_left_down_all = [];
frames_left_up_all = [];
frames_reward_all = [];
frames_unrewarded_all = [];
frames_cue_onset_all = [];
timeDiff = [];
reward_trial = [];
unrewarded_trial = [];

for i = 1:length(img_filename_all);
    % get current img filename
    img_filename = [img_path img_filename_all{i}];
    
    % get current acq filename
    acq_filename = [acq_path acq_filename_all{i}];
    
    % this is hard coded until the header problem is fixed
    % % %     % pull out time for each for each frame from header
    % % %     img_info = imfinfo([img_pathname img_filename]);
    % % %     msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
    % % %     msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
    % % %     msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
    % % %     lines = img_info(1).Height;
    % % %     fpms = msPerLine*lines;
    
    msPerLine = 1.24;
    lines = 128;
    fpms = msPerLine*lines;
    
    % Get trials in s since started
    TrialNumList = TK_ReadBitCode([acq_filename]);
    close(gcf)
    % put trial time in ms
    TrialNumList(2,:) = TrialNumList(2,:)*1000;
    
    left_down = [];
    left_up = [];
    reward = [];
    unrewarded = [];
    cue_onset = [];
    
    for trial = TrialNumList(1,:)
        if trial > length(bhv.saved_history.RewardsSection_LastTrialEvents);
            continue
        end
        bhv_current = bhv.saved_history.RewardsSection_LastTrialEvents{trial,1};
        % fixes possible bug of putting in old timepoints at the end
        LastGoodTrial = find(diff(bhv_current(:,3)) < 0,1);
        if ~isempty(LastGoodTrial)
            bhv_current = bhv_current(1:LastGoodTrial,:);
        end
        % get time lag of behavior to acq (t_bhv + timeDiff = t_acq) in ms
        trialStart_bhv = bhv_current(find(bhv_current(:,1) == 101 & bhv_current(:,2) == 0,1,'last'),3)*1000;
        lag_time = TrialNumList(2,find(TrialNumList(1,:) == trial)) - trialStart_bhv;
        timeDiff(trial,1) = lag_time;
        % make a second column for which acquisition you're in to frame shift
        timeDiff(trial,2) = i*ones(length(lag_time),1);
        % time for all lever up and down (left down = 5, left up = 6) in ms
        %left_down = [left_down; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff(trial)];
        %left_up = [left_up; bhv_current(find(bhv_current(:,2) == 6),3)*1000 + timeDiff(trial)];
        reward = [reward; bhv_current(find(bhv_current(:,1) == 44,1),3)*1000 + timeDiff(trial,1)];
        unrewarded = [unrewarded; bhv_current(find(bhv_current(:,1) == 45,1),3)*1000 + timeDiff(trial,1)];
        if ~isempty(find(bhv_current(:,1) == 44,1))
            reward_trial = [reward_trial trial];
        end
        if ~isempty(find(bhv_current(:,1) == 45,1))
            unrewarded_trial = [unrewarded_trial trial];
        end        
        cue_onset = [cue_onset; bhv_current(find(bhv_current(:,1) == 41,1),3)*1000 + timeDiff(trial,1)];
    end
    
    % concatenate all events, multiply by img number to get total time
    
    %left_down_all = [left_down_all; left_down];
    %left_up_all = [left_up_all; left_up];
    unrewarded_all = [unrewarded_all; unrewarded];
    reward_all = [reward_all; reward];
    cue_onset_all = [cue_onset_all; cue_onset];
    
    %frames_left_down = ceil(left_down(~isnan(left_down))/fpms)';
    %frames_left_up = ceil(left_up(~isnan(left_up))/fpms)';
    frames_reward = ceil(reward(~isnan(reward))/fpms)';
    frames_unrewarded = ceil(unrewarded(~isnan(unrewarded))/fpms)';
    frames_cue_onset = ceil(cue_onset(~isnan(cue_onset))/fpms)';
    
    % round up for number of channels
    %frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
    %frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
    frames_cue_onset = (1-mod(frames_cue_onset,n_channels))+frames_cue_onset;
    frames_reward = (1-mod(frames_reward,n_channels))+frames_reward;
    frame_unrewarded = (1-mod(frames_unrewarded,n_channels))+frames_unrewarded;
    
    % Get frame for each event
    %frames_left_down_all = [frames_left_down_all frames_left_down+n_frames*(i-1)];
    %frames_left_up_all = [frames_left_up_all frames_left_up+n_frames*(i-1)];
    frames_reward_all = [frames_reward_all frames_reward+n_frames*(i-1)];
    frames_unrewarded_all = [frames_unrewarded_all frames_unrewarded+n_frames*(i-1)];    
    frames_cue_onset_all = [frames_cue_onset_all frames_cue_onset+n_frames*(i-1)];
    
    waitbar(i/length(img_filename_all),w,'Getting and concatenating events');
    
    
end

% get all lever timing at once to get rid of artifacts correctly
AP_getLeftLeverTimesImg
% correct lever times for time delay, keep old ones
left_down_acq = left_down.*1000 + timeDiff(trial_down,1);
trial_down_acq = trial_down;
states_down_acq = states_down;
left_up_acq = left_up.*1000 + timeDiff(trial_up,1);
trial_up_acq = trial_up;
states_up_acq = states_up;
% get rid of lever times when trial not recorded
trial_notRecorded = find(timeDiff(:,1) == 0);
elim_left_down = ismember(trial_down,trial_notRecorded);
elim_left_up = ismember(trial_up,trial_notRecorded);

left_down_acq(elim_left_down) = [];
trial_down_acq(elim_left_down) = [];
states_down_acq(elim_left_down) = [];
left_up_acq(elim_left_up) = [];
trial_up_acq(elim_left_up) = [];
states_up_acq(elim_left_up) = [];
% pick out the left_down that were in the right place (state 41)
left_down_air = left_down_acq(states_down_acq == 41);

% convert lever times to frames
frames_left_down = ceil(left_down_acq(~isnan(left_down_acq))/fpms)';
frames_left_up = ceil(left_up_acq(~isnan(left_up_acq))/fpms)';
% round frames based on channels
frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
% account for concatenation
% get amount of frames to shift based on trial number
frame_shift_left_down = [];
frame_shift_left_up = [];
frame_shift_left_down = (timeDiff(trial_down_acq,2)-1)*n_frames;
frame_shift_left_up = (timeDiff(trial_up_acq,2)-1)*n_frames;
frames_left_down_all = [frames_left_down_all frames_left_down'+frame_shift_left_down];
frames_left_up_all = [frames_left_up_all frames_left_up'+frame_shift_left_up];

% get the frames where left lever down during correct time
frames_left_down_air_all = frames_left_down_all(states_down_acq == 41);

close(w)

%% Make traces for plots, etc

mod_cells = mod_up_cells;
%mod_cells = 54;

% hist_flag = 1 if you want to see histograms and raw trace with events
hist_flag = 0;
% reward_trace_flag = 1 if you want to see all traces plotted for each rew
reward_trace_flag = 0;

frame_sample = 30;
reward_trace_concat = zeros(length(frames_reward_all),frame_sample*2+1,length(mod_cells));

%for curr_mod = 1:size(roi_trace_long,1)
for curr_mod = 1:length(mod_cells);
    
    curr_roi = mod_cells(curr_mod);
    %curr_roi = curr_mod;
    
    % get traces time-locked to events
    reward_trace = zeros(length(frames_reward_all),frame_sample*2+1);
    for reward_num = 1:length(frames_reward_all);
        if frames_reward_all(reward_num)-frame_sample > 1 && ...
                frames_reward_all(reward_num)+frame_sample < size(roi_trace_long,2);
            reward_trace(reward_num,:) = roi_trace_long...
                (curr_roi,frames_reward_all(reward_num)-frame_sample:frames_reward_all(reward_num)+frame_sample);
        else
            reward_trace(reward_num,:) = NaN;
        end
    end
    
    left_down_trace = [];
    for left_down_num = 1:length(frames_left_down_all);
        if frames_left_down_all(left_down_num)-frame_sample > 1 && ...
                frames_left_down_all(left_down_num)+frame_sample < size(roi_trace_long,2);
            left_down_trace(left_down_num,:) = roi_trace_long...
                (curr_roi,frames_left_down_all(left_down_num)-frame_sample:frames_left_down_all(left_down_num)+frame_sample);
        end
    end
    
    left_up_trace = [];
    for left_up_num = 1:length(frames_left_up_all);
        if frames_left_up_all(left_up_num)-frame_sample > 1 && ...
                frames_left_up_all(left_up_num)+frame_sample < size(roi_trace_long,2);
            left_up_trace(left_up_num,:) = roi_trace_long...
                (curr_roi,frames_left_up_all(left_up_num)-frame_sample:frames_left_up_all(left_up_num)+frame_sample);
        end
    end
    
    % get the left downs and up that count (down during 41, up during 43)
    left_down_cued = find(states_down_acq == 41);
    left_up_cued = find(states_up_acq == 43);
    
    left_down_cued_trace = [];
    for left_down_cued_num = left_down_cued';
        if frames_left_down_all(left_down_cued_num)-frame_sample > 1 && ...
                frames_left_down_all(left_down_cued_num)+frame_sample < size(roi_trace_long,2);
            left_down_cued_trace(find(left_down_cued_num == left_down_cued),:) = roi_trace_long...
                (curr_roi,frames_left_down_all(left_down_cued_num)-frame_sample:frames_left_down_all(left_down_cued_num)+frame_sample);
        end
    end
    
    left_up_cued_trace = [];
    for left_up_cued_num = left_up_cued';
        if frames_left_up_all(left_up_cued_num)-frame_sample > 1 && ...
                frames_left_up_all(left_up_cued_num)+frame_sample < size(roi_trace_long,2);
            left_up_cued_trace(find(left_up_cued_num == left_up_cued),:) = roi_trace_long...
                (curr_roi,frames_left_up_all(left_up_cued_num)-frame_sample:frames_left_up_all(left_up_cued_num)+frame_sample);
        end
    end
    
    
    % get time differences for events around peaks
    %peak_frames = roi_peak_max{curr_roi}(:,1);
    peak_frames = spike_frames{curr_roi};
    frame_search = 5;
    left_down_peak_diff = [];
    left_down_peak_event = [];
    left_up_peak_diff = [];
    left_up_peak_event = [];
    reward_peak_diff = [];
    reward_peak_event = [];
    
    for curr_peak = 1:length(peak_frames)
        left_down_peak = [];
        left_up_peak = [];
        reward_peak = [];
        
        frame_range = peak_frames(curr_peak)-frame_search:peak_frames(curr_peak);
        
        left_down_peak = ismember(frames_left_down_all,frame_range);
        left_down_peak_diff = [left_down_peak_diff; frames_left_down_all(left_down_peak)-peak_frames(curr_peak)];
        if sum(left_down_peak)~=0
            left_down_peak_event(curr_peak) = 1;
        end
        left_up_peak = ismember(frames_left_up_all,frame_range);
        left_up_peak_diff = [left_up_peak_diff; frames_left_up_all(left_up_peak)-peak_frames(curr_peak)];
        if sum(left_up_peak)~=0
            left_up_peak_event(curr_peak) = 1;
        end
        reward_peak = ismember(frames_reward_all,frame_range);
        reward_peak_diff = [reward_peak_diff; frames_reward_all(reward_peak)-peak_frames(curr_peak)];
        if sum(reward_peak)~=0
            reward_peak_event(curr_peak) = 1;
        end
        
    end
    
    % average all peaks
    peak_trace = [];
    for peak_num = 1:length(peak_frames);
        % only take peaks with frame_sample = surrounding frames
        if peak_frames(peak_num)-frame_sample > 1 && ...
                peak_frames(peak_num)+frame_sample < size(roi_trace_long,2)
            peak_trace(peak_num,:) = roi_trace_long...
                (curr_roi,peak_frames(peak_num)-frame_sample:peak_frames(peak_num)+frame_sample);
        end
    end
    
    if hist_flag == 1;
        % create stacked histogram for peak events
        
        % get the max peak for bin, round to nearest tenth
        peak_max_max = ceil(max(roi_trace_long(curr_roi,peak_frames))*10)/10;
        
%         [hist_all_peak_n hist_all_peak_x] = hist(spike_amplitudes{curr_roi},[1:0.1:peak_max_max]);
%         [hist_left_down_peak_n hist_left_down_peak_x] = ...
%             hist(spike_amplitudes{curr_roi}((left_down_peak_event == 1)),[1:0.1:peak_max_max]);
%         [hist_left_up_peak_n hist_left_up_peak_x] = ...
%             hist(spike_amplitudes{curr_roi}((left_up_peak_event == 1)),[1:0.1:peak_max_max]);
%         [hist_reward_peak_n hist_reward_peak_x] = ...
%             hist(spike_amplitudes{curr_roi}((reward_peak_event == 1)),[1:0.1:peak_max_max]);
%         
%         figure;
%         bar(hist_all_peak_x,[hist_all_peak_n; hist_reward_peak_n; hist_left_up_peak_n; hist_left_down_peak_n]');
%         colormap(summer)
%         legend('All','Reward','Left up','Left down');
%         title(['ROI ' num2str(curr_roi)]);
         
        % plot raw trace with ticks: down = black, up = red, reward = blue
        figure;hold on;
        plot(roi_trace_long(curr_roi,:))
        for i = 1:length(frames_left_down_all)
            line([frames_left_down_all(i) frames_left_down_all(i)], [0 0.2],'color','black')
        end
        for i = 1:length(frames_left_up_all)
            line([frames_left_up_all(i) frames_left_up_all(i)], [0.3 0.5],'color','red')
        end
        for i = 1:length(frames_reward_all)
            line([frames_reward_all(i) frames_reward_all(i)], [0.6 0.8],'color','blue')
        end
        plot(spike_frames{curr_roi},spike_amplitudes{curr_roi},'.r');
    end
    
    if reward_trace_flag == 1;
        figure;
        imagesc(reward_trace);
        caxis([0 2]);
        title(num2str(curr_roi));
        x_sec = round(((str2num(get(gca,'XTickLabel')) - frame_sample+1).*(1.24*128/1000)).*10)./10;
        set(gca,'XTickLabel',x_sec);
        xlabel('Time from reward (s)')
        ylabel('Reward number')
    end
    
    reward_trace_concat(:,:,curr_mod) = reward_trace;
    
end

% plot average response
sqr_num = ceil(sqrt(size(reward_trace_concat,3)));
figure;
for i = 1:size(reward_trace_concat,3)
    subplot(sqr_num,sqr_num,i);
    plot(nanmean(reward_trace_concat(:,:,i)),'linewidth',2,'color','k')
    hold on
    nanstd_reward_trace_u = nanmean(reward_trace_concat(:,:,i))+...
        nanstd(reward_trace_concat(:,:,i));
    nanstd_reward_trace_l = nanmean(reward_trace_concat(:,:,i))-...
        nanstd(reward_trace_concat(:,:,i));
    nansem_reward_trace_u = nanmean(reward_trace_concat(:,:,i)) +...
        ((nanstd(reward_trace_concat(:,:,i))/sqrt...
        (size(reward_trace_concat(:,:,i),1))));
    nansem_reward_trace_l = nanmean(reward_trace_concat(:,:,i)) -...
        ((nanstd(reward_trace_concat(:,:,i))/sqrt...
        (size(reward_trace_concat(:,:,i),1))));
    x_reward = [1:1:size(reward_trace_concat(:,:,i),2)];
    %jbfill(x_reward,nanstd_reward_trace_u,...
    %    nanstd_reward_trace_l,'red','none',0,0.5);
    jbfill(x_reward,nansem_reward_trace_u,...
        nansem_reward_trace_l,'blue','none',0,0.5);
    set(gca,'xlim',[0 2*frame_sample+1]);
    x_sec = round(((str2num(get(gca,'XTickLabel')) - frame_sample-1).*(1.24*128/1000)).*10)./10;
    set(gca,'XTickLabel',x_sec);
    ylabel('{\Delta}F/F_0')
    xlabel('Time from reward (s)')
    title(num2str(mod_cells(i)));
    line([frame_sample frame_sample],ylim,'linestyle','--','color','k');
end


%% Plot events with all peaks

figure;hold on;
for i = 1:size(roi_trace_long,1)
    for j = 1:size(spike_frames{i})
        % color lines by amplitude - darker is larger, -0.2 to make all vis
        if spike_amplitudes{i}(j) > 0; %FIX THIS! spike amplitudes should never be neg, but are sometimes
            color_amp = 1-(spike_amplitudes{i}(j)/(max(spike_amplitudes{i})+0.2));
            line([spike_frames{i}(j) spike_frames{i}(j)], ...
                [i i+0.9],'color',[color_amp color_amp color_amp]);
        end
    end
end

topline = length(spike_frames);
for i = 1:length(frames_left_down_all)
    line([frames_left_down_all(i) frames_left_down_all(i)], [topline+1 topline+1.9],'color','green')
end
for i = 1:length(frames_left_up_all)
    line([frames_left_up_all(i) frames_left_up_all(i)], [topline+2 topline+2.9],'color','red')
end
for i = 1:length(frames_reward_all)
    line([frames_reward_all(i) frames_reward_all(i)], [topline+3 topline+3.9],'color','blue')
end
for i = 1:length(frames_cue_onset_all)
    line([frames_cue_onset_all(i) frames_cue_onset_all(i)], [topline+4 topline+4.9],'color','cyan')
end
ylim([0 topline+5]);

%% Plot stack of traces
figure;
hold on;
for i = 1:size(mod_up_cells,1)
    %for i = 1:length(temp)
    plot(roi_trace_long(mod_up_cells(i),:) + i,'-k');
end
 
% % plot lever down as dotted black lines
% for i = 1:length(frames_left_down_all)
%     line([frames_left_down_all(i) frames_left_down_all(i)],ylim,'linestyle','--','color','k');
% end
% 
% % plot lever up as dotted blue lines
% for i = 1:length(frames_left_up_all)
%     line([frames_left_up_all(i) frames_left_up_all(i)],ylim,'linestyle','--','color','b');
% end

% plot rewards as dashed red lines
for i = 1:length(frames_reward_all)
    line([frames_reward_all(i) frames_reward_all(i)],ylim,'linestyle','--','color','r');
end

% plot non-rewards as dashed magenta lines
for i = 1:length(frames_unrewarded_all)
    line([frames_unrewarded_all(i) frames_unrewarded_all(i)],ylim,'linestyle','--','color','m');
end

% % plot air onset as dotted blue lines
% for i = 1:length(frames_cue_onset_all)
%     line([frames_cue_onset_all(i) frames_cue_onset_all(i)],ylim,'linestyle',':','color','b');
% end

%% Show which ROIs had events around given timeframe
frame_back = 10;
frame_forward = 10;
% events around reward/left up
surround_reward = zeros(size(roi_trace_long,1),length(frames_reward_all));
for curr_reward = 1:length(frames_reward_all)
    curr_frame = frames_reward_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_reward(curr_roi,curr_reward) = event_detect_amp;
        end
    end
end
figure;imagesc(surround_reward);
title('Events around reward/left up')

% get behavior-related cells:

% events around left down during air
surround_left_down_air = zeros(size(roi_trace_long,1),length(frames_left_down_air_all));
for curr_left_down_air = 1:length(frames_left_down_air_all)
    curr_frame = frames_left_down_air_all(curr_left_down_air);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_left_down_air(curr_roi,curr_left_down_air) = event_detect_amp;
        end
    end
end
figure;imagesc(surround_left_down_air);
title('Events around left down during air')

%% Find reward-related cells (requires last cell be run)
tic
n_sample = 1000;

frame_back = 10;
frame_forward = 10;

surround_random_all = zeros(size(roi_trace_long,1),n_sample);
frames_random = zeros(n_sample, length(frames_reward_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_long,2),length(frames_reward_all));
end

w = waitbar(0)
for curr_roi = 1:length(spike_frames)
    frames_random_spike = [];
    frames_random_spike = repmat(frames_random,[1 1 length(spike_frames{curr_roi})]) ;
    % get difference between random frames and spikes
    for curr_spike = 1:length(spike_frames{curr_roi});
        frames_random_spike(:,:,curr_spike) = frames_random_spike(:,:,curr_spike) - spike_frames{curr_roi}(curr_spike);
    end
    % find any of that matrix that are between frame_back and frame_forward
    frames_random_spike = frames_random_spike > -frame_back & frames_random_spike < frame_forward;
    surround_random_all(curr_roi,:) = sum(any(frames_random_spike,3),2);
    disp(['ROI' num2str(curr_roi)]);
    waitbar(curr_roi/length(spike_frames),w);
end
close(w);

% make rewarded events binary
surround_reward_binary = surround_reward;
surround_reward_binary(surround_reward_binary > 0) = 1;
surround_reward_binary = sum(surround_reward_binary,2);

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [surround_random_all(curr_roi,:) surround_reward_binary(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end

mod_up_cells = [];
mod_down_cells = [];

mod_up_cells = find(p > 0.95);
mod_down_cells = find(p < 0.05);
toc

%% show up/down/non-modulated cells

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[0.5 0.5 0.5]);
end

for i = 1:length(mod_up_cells)
    patch(polygon.ROI{mod_up_cells(i)}(:,1),-polygon.ROI{mod_up_cells(i)}(:,2),[0 1 0]);
end

for i = 1:length(mod_down_cells)
    patch(polygon.ROI{mod_down_cells(i)}(:,1),-polygon.ROI{mod_down_cells(i)}(:,2),[1 0 0]);
end

title('Cells correlated with trained task, p > 0.05')

figure
hold on
for i = 1:length(polygon.ROI)
    patch(polygon.ROI{i}(:,1),-polygon.ROI{i}(:,2),[1-p(i) p(i) 0]);
end
title('Cells correlated with trained task, graded')

%% Plot by kmeans grouping for general trace
[temp1 centers] = kmeans(roi_trace_long,2);
figure;hold on;
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

for i = 1:length(a)
    patch(polygon.ROI{a(i)}(:,1),polygon.ROI{a(i)}(:,2),[1 0 0]);
end
for i = 1:length(b)
    patch(polygon.ROI{b(i)}(:,1),polygon.ROI{b(i)}(:,2),[0 1 0]);
end
for i = 1:length(c)
    patch(polygon.ROI{c(i)}(:,1),polygon.ROI{c(i)}(:,2),[0 0 1]);
end
title('Cells sorted through K Means')

% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    temp_trace(i,:) = roi_trace_long(sort_temp1(i,2),:);
end
figure;imagesc(temp_trace);
figure;imagesc(corrcoef(temp_trace'))
figure;plot(centers(1,:));hold on;plot(centers(2,:),'--r');
title('K-means centers');
legend('Center 1','Center 2');

%% Plot by kmeans grouping for which rewards had events around them
surround_reward_bin = surround_reward;
surround_reward_bin(surround_reward > 0) = 1;

[temp1 centers] = kmeans(surround_reward_bin,2);
figure;hold on;
a = find(temp1 == 1);
b = find(temp1 == 2);
c = find(temp1 == 3);

for i = 1:length(a)
    patch(polygon.ROI{a(i)}(:,1),polygon.ROI{a(i)}(:,2),[1 0 0]);
end
for i = 1:length(b)
    patch(polygon.ROI{b(i)}(:,1),polygon.ROI{b(i)}(:,2),[0 1 0]);
end
for i = 1:length(c)
    patch(polygon.ROI{c(i)}(:,1),polygon.ROI{c(i)}(:,2),[0 0 1]);
end
title('Cells sorted through K Means')

% sort roi_trace_long by clusters
temp1(:,2) = [1:length(temp1)]';
sort_temp1 = sortrows(temp1,1);
temp_trace = zeros(size(surround_reward));
for i = 1:size(surround_reward,1);
    temp_trace(i,:) = surround_reward(sort_temp1(i,2),:);
end
figure;imagesc(temp_trace);
figure;imagesc(corrcoef(temp_trace'))
figure;plot(centers(1,:));hold on;plot(centers(2,:),'--r');
title('K-means centers');
legend('Center 1','Center 2');

% plot the correlations between groups

corrcoef_1  = corrcoef(roi_trace_long(a,:)');
corrcoef_1(corrcoef_1 == 1) = NaN;
mean_corrcoef_1 = nanmean(corrcoef_1(:));
sem_corrcoef_1 = nanstd(corrcoef_1(:))./sqrt(size(roi_trace_long,1));

corrcoef_2  = corrcoef(roi_trace_long(b,:)');
corrcoef_2(corrcoef_2 == 1) = NaN;
mean_corrcoef_2 = nanmean(corrcoef_2(:));
sem_corrcoef_2 = nanstd(corrcoef_2(:))./sqrt(size(roi_trace_long,1));

corrcoef_all  = corrcoef(roi_trace_long');
corrcoef_all(corrcoef_all == 1) = NaN;
mean_corrcoef_all = nanmean(corrcoef_all(:));
sem_corrcoef_all = nanstd(corrcoef_all(:))./sqrt(size(roi_trace_long,1));

means = [mean_corrcoef_1 mean_corrcoef_2 mean_corrcoef_all];
sems = [sem_corrcoef_1 sem_corrcoef_2 sem_corrcoef_all];

figure; hold on;
h = bar(means)
% this can't be copied and pasted, very lame
%ch = get(h,'Children');
%set(ch,'FaceVertexCData',[1 0 0 ; 0 1 0; .5 .5 .5])
errorb(means,sems,'top')
colormap(gray)
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'Group 1','Group 2','All'});
colormap(gray)
ylim([0 0.3])
ylabel('Correlation')

%% Find unrewarded-related cells
% find which ROIs had events around given timeframe
frame_back = 0;
frame_forward = 20;
% events around reward/left up
surround_unrewarded = zeros(size(roi_trace_long,1),length(frames_unrewarded_all));
for curr_reward = 1:length(frames_unrewarded_all)
    curr_frame = frames_unrewarded_all(curr_reward);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_unrewarded(curr_roi,curr_reward) = event_detect_amp;
        end
    end
end
figure;imagesc(surround_unrewarded);
title('Events around reward/left up')

% get behavior-related cells:

% events around left down during air
surround_left_down_air = zeros(size(roi_trace_long,1),length(frames_left_down_air_all));
for curr_left_down_air = 1:length(frames_left_down_air_all)
    curr_frame = frames_left_down_air_all(curr_left_down_air);
    curr_frame_search = curr_frame-frame_back:curr_frame+frame_forward;
    for curr_roi = 1:length(spike_frames)
        event_detect = [];
        event_detect = ismember(spike_frames{curr_roi},curr_frame_search);
        if any(event_detect)
            % get the largest event which occured around the searched time
            event_detect_amp = max(spike_amplitudes{curr_roi}(event_detect));
            surround_left_down_air(curr_roi,curr_left_down_air) = event_detect_amp;
        end
    end
end


tic
n_sample = 10000;

surround_random_all = zeros(size(roi_trace_long,1),n_sample);
frames_random = zeros(n_sample, length(frames_unrewarded_all));
% create matrix for random reward times
for num_rand = 1:n_sample;
    frames_random(num_rand,:) = ...
        randsample(size(roi_trace_long,2),length(frames_unrewarded_all));
end

w = waitbar(0)
for curr_roi = 1:length(spike_frames)
    frames_random_spike = [];
    frames_random_spike = repmat(frames_random,[1 1 length(spike_frames{curr_roi})]) ;
    % get difference between random frames and spikes
    for curr_spike = 1:length(spike_frames{curr_roi});
        frames_random_spike(:,:,curr_spike) = frames_random_spike(:,:,curr_spike) - spike_frames{curr_roi}(curr_spike);
    end
    % find any of that matrix that are between frame_back and frame_forward
    frames_random_spike = frames_random_spike > -frame_back & frames_random_spike < frame_forward;
    surround_random_all(curr_roi,:) = sum(any(frames_random_spike,3),2);
    disp(['ROI' num2str(curr_roi)]);
    waitbar(curr_roi/length(spike_frames),w);
end
close(w);

% make rewarded events binary
surround_unrewarded_binary = surround_unrewarded;
surround_unrewarded_binary(surround_unrewarded_binary > 0) = 1;
surround_unrewarded_binary = sum(surround_unrewarded_binary,2);

p = zeros(size(roi_trace_long,1),1);
% loop through each ROI, get percentile of events around reward from random
for curr_roi = 1:size(roi_trace_long,1)
    % tack on test value to end of random distribution
    distr_test = [surround_random_all(curr_roi,:) surround_unrewarded_binary(curr_roi)];
    % rank all values, nonparametric way to get percentile
    distr_rank = tiedrank(distr_test);
    % percentile is the ratio of the test value to all total ranks
    p(curr_roi) = distr_rank(end)/length(distr_rank);
end

mod_up_cells_unrewarded = [];
mod_down_cells_unrewarded = [];

mod_up_cells_unrewarded = find(p > 0.95);
mod_down_cells_unrewarded = find(p < 0.05);
toc

%% Divide caltrace conts y by 1/4 and call them polygon.ROI
polygon.ROI = cell(size(CONTS));
for i = 1:length(CONTS)
    polygon.ROI{i}(:,1) = CONTS{i}(:,1);
    polygon.ROI{i}(:,2) = CONTS{i}(:,2)/4;
end

%% run through with peakdet, as a prelim for finding events
w = waitbar(0)
for i = 1:size(roi_trace_long,1);
    delta = abs(2*std(roi_trace_long(i,:)));
    if delta == 0
        continue
    end
    [maxtab{i} mintab{i}] = peakdet(roi_trace_long(i,:),delta);
    if isempty(maxtab{i})
        maxtab{i} = [0 0];
    end
    waitbar(i/size(roi_trace_long,1),w);
end
close(w);

spike_frames = cell(size(maxtab));
spike_frames = cellfun(@(x) x(:,1), maxtab,'UniformOutput',0);

