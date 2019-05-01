%% Get correlation between movie pixels and lever movement (rewarded + all)
tic
% Pick xsg files first
[xsg_filename,xsg_path]=uigetfile('*.xsg','Choose XSG files','Multiselect','on');
cd(xsg_path)
cd ..
if ~iscell(xsg_filename)
    xsg_filename = mat2cell(xsg_filename)
end

% pick tiffs
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');

% pick behavior file, parse relevant behavior
[bhv_filename bhv_path] = uigetfile('*','Behavior file');
bhv = load([bhv_path bhv_filename],'-MAT');
parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;

num_trials = length(parsed_events);

% Get timing from behavior file: rewards, cues, itis, trial starts, trial
% ends, which trials were rewarded, rewarded lever presses
all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);

% parallel script
try
    matlabpool
catch me
end

t_date = [];
t_mouse = [];
t_loop = [];
t_subloop = [];
t_correction = [];
for i = 1:length(tiff_filename)
    [t_date_curr t_mouse_curr t_loop_curr t_subloop_curr t_correction] = strread(tiff_filename{i},'%d %s %d %d %s.tif', 'delimiter','_');
    t_loop = [t_loop t_loop_curr];
    t_subloop = [t_subloop t_subloop_curr];
end

for curr_xsg = 1:length(xsg_filename)
    
    % get loop number, find the tiff files, check for all files present
    curr_loop = str2num(xsg_filename{curr_xsg}(end-7:end-4));
    tiff_curr_loop = find(t_loop == curr_loop);
    if all(t_subloop(tiff_curr_loop) ~= [1:8])
        disp(['Don''t have 1-8 for loop' num2str(curr_loop) ',skipping']);
        continue
    end
    tiff_filename_curr = tiff_filename(tiff_curr_loop);
    
    % load in tiff images, columns are frames
    im_temp = [];
    im = [];
    
    for i = 1:length(tiff_filename_curr);
        img_filename = [tiff_path tiff_filename_curr{i}];
        imageinfo=imfinfo(img_filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        clear im_temp_long
        clear im_temp
        
        disp('Loading file (native format)....')
        parfor loadframe = 1:numframes
            im_temp = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
            im_temp_long(:,loadframe) = im_temp(:);
        end
        
        im = [im im_temp_long];
        disp(['Done ' num2str(i) '/' num2str(length(tiff_filename_curr))])
    end
    clear im_temp_long
    
    % load corresponding xsg
    xsg = load([xsg_path xsg_filename{curr_xsg}],'-MAT');
    lever_force = xsg.data.acquirer.trace_2;
    
    X = M;
    Y = N;
    numframes = size(im,2);
    % to make resampling work: numframes must be even number
    % if it's not, cut out one frame at the end
    if mod(numframes,2) ~= 0
        im(:,end) = [];
        numframes = size(im,2);
    end
    lever_force_resample = resample(lever_force,numframes,length(lever_force));
    
    % Find frame for each reward, grab surrounding frames and resampled force
    
    % Parse header, get framerate
    img_info = imageinfo(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    % sample rate for xsg
    xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
    
    % Get trials in seconds since started
    curr_trial_list = [];
    trial_list = [];
    trial_stitch = []; % unused at the moment since looking only at single xsg
    curr_trial_list = AP_ReadBitCode([xsg_path xsg_filename{curr_xsg}]);
    if ~isempty(curr_trial_list)
        % this is commented out for now b/c looking at 1 xsg at a time
        %         curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
        trial_list = [trial_list; curr_trial_list];
    end
    
    % Get all trial offsets in seconds
    xsg_bhv_offset = zeros(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            xsg_bhv_offset(curr_trial) = NaN;
        end
    end
    
    % Get lever presses when rewarded
    reward_trials_num = find(reward_trials == 1);
    lever_surround_sec = 1;
    lever_surround_sample = round((lever_surround_sec*framerate)/2);
    lever_force_resample_reward = [];
    lever_force_resample_reward_parsed = {};
    force_point = nan(num_trials,1);
    im_reward = [];
    for curr_reward_trial = reward_trials_num';
        % Skip this if trial start wasn't recorded
        if ~isnan(xsg_bhv_offset(curr_reward_trial)) & ...
                ~ismember(curr_reward_trial,trial_stitch)
            % Current xsg trial time
            xsg_trial_indx = find(trial_list(:,2) == curr_reward_trial,1);
            lever_frame_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(framerate);
            % xsg lever time = relative lever sample trial start + xsg trial start
            lever_frame_xsg = round(lever_frame_bhv + trial_list(xsg_trial_indx,1)*framerate);
            % check that the frames in the range exist, take what's there
            clear frame_range
            frame_range = lever_frame_xsg-lever_surround_sample:lever_frame_xsg+lever_surround_sample;
            frame_range(frame_range < 1) = [];
            frame_range(frame_range > size(im,2)) = [];
            % grab surrounding frames
            lever_force_resample_reward = [lever_force_resample_reward ...
                lever_force_resample(frame_range)'];
            curr_num_rewards = length(lever_force_resample_reward_parsed);
            lever_force_resample_reward_parsed{curr_num_rewards+1} = ...
                lever_force_resample(frame_range)';
            im_reward = [im_reward ...
                im(:,frame_range)];
            force_point(curr_reward_trial) = lever_frame_xsg;
        end
    end
    
    cc_img = zeros(Y*X,1);
    cc_img_reward = zeros(Y*X,1);
    parfor px_curr = 1:(X*Y);
        % slice the image variable for parallel
        im_slice = im(px_curr,:);
        im_reward_slice = im_reward(px_curr,:);
        
        % skip first and last 10 frames to avoid resampling edge effects
        im_slice = double(im_slice(10:end-10));
        im_reward_slice = double(im_reward_slice(10:end-10));
        
        curr_cc = corrcoef(im_slice,lever_force_resample(10:end-10));
        curr_cc_reward = corrcoef(im_reward_slice,lever_force_resample_reward(10:end-10));
        cc_img(px_curr) = curr_cc(2);
        cc_img_reward(px_curr) = curr_cc_reward(2);
        disp([num2str(px_curr) '/' num2str(X*Y)]);
    end
    
    curr_date = tiff_filename{1}(1:6);
    savename = ['/usr/local/lab/People/Andy/Data/CorrData/' curr_date '_' xsg_filename{curr_xsg}(1:end-4) '_reward'];
    save(savename,'cc_img','cc_img_reward','im_reward','lever_force_resample_reward','xsg_filename','tiff_filename');
    disp(['Saved ' savename]);
end
toc
disp(['Finished analysis']);


%% plot those results
plot_all_flag = 1;

[cc_filename,cc_path]=uigetfile('*.mat','Choose CC files','Multiselect','on');

if ~iscell(cc_filename)
    cc_filename = {cc_filename};
end
cc_dates = {};
cc_dates = cellfun(@(x) x(1:6),cc_filename,'UniformOutput',0);

unique_cc_dates = unique(cc_dates);

for curr_date = 1:length(unique_cc_dates)
    curr_date_indx = strcmp(cc_dates,unique_cc_dates(curr_date));
    cc_filename_curr_date = cc_filename(curr_date_indx);
    
    cc_all = [];
    cc_reward_all = [];
    for i = 1:length(cc_filename_curr_date)
        load([cc_path cc_filename_curr_date{i}], 'cc_img','cc_img_reward','lever_force_resample_reward');
        cc_im_reshape = reshape(cc_img,512,512);
        cc_im_reward_reshape = reshape(cc_img_reward,512,512);
        if plot_all_flag == 1;
            figure;imagesc(cc_im_reshape);
            caxis([-0.2 0.2])
            title([unique_cc_dates{curr_date} ' ' num2str(i)])
            figure;imagesc(cc_im_reward_reshape);
            caxis([-0.2 0.2])
            title([unique_cc_dates{curr_date} ': Reward ' num2str(i)])
            if mod(length(lever_force_resample_reward),29) == 0
                lever_force_resample_reward_reshape = reshape(...
                    lever_force_resample_reward,29,length(lever_force_resample_reward)/29);
                figure;plot(lever_force_resample_reward_reshape)
                title([unique_cc_dates{curr_date} ' ' num2str(i)])
            end
        end
        
        cc_all = [cc_all cc_img];
        cc_reward_all = [cc_reward_all cc_img_reward];
    end
    clear cc_mean cc_mean_reward
    cc_mean = nanmean(cc_all,2);
    cc_mean_reward = nanmean(cc_reward_all,2);
    clear cc_mean_reshape cc_mean_reward_reshape
    cc_mean_reshape = reshape(cc_mean,512,512);
    cc_mean_reward_reshape = reshape(cc_mean_reward,512,512);
    
    figure;imagesc(cc_mean_reshape);
    title(['Mean CC: ' unique_cc_dates{curr_date}]);
    caxis([-0.2 0.2]);
    
    figure;imagesc(cc_mean_reward_reshape);
    title(['Mean CC: ' unique_cc_dates{curr_date} ': Reward']);
    caxis([-0.2 0.2]);
    
end

%% Get correlation between movie pixels and lever movement (rewarded + all)
%%% WHOLE MOVIE

tic

% initialize concat vars
im_concat = [];
im_reward_concat = [];
lever_force_resample_concat = [];
lever_force_resample_reward_concat = [];
    
% Pick xsg files first
[xsg_filename,xsg_path]=uigetfile('*.xsg','Choose XSG files','Multiselect','on');
cd(xsg_path)
cd ..
if ~iscell(xsg_filename)
    xsg_filename = mat2cell(xsg_filename)
end

% pick tiffs
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');

% pick behavior file, parse relevant behavior
[bhv_filename bhv_path] = uigetfile('*','Behavior file');
bhv = load([bhv_path bhv_filename],'-MAT');
parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;

num_trials = length(parsed_events);

% Get timing from behavior file: rewards, cues, itis, trial starts, trial
% ends, which trials were rewarded, rewarded lever presses
all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);

% parallel script
try
    matlabpool
catch me
end

t_date = [];
t_mouse = [];
t_loop = [];
t_subloop = [];
t_correction = [];
for i = 1:length(tiff_filename)
    [t_date_curr t_mouse_curr t_loop_curr t_subloop_curr t_correction] = strread(tiff_filename{i},'%d %s %d %d %s.tif', 'delimiter','_');
    t_loop = [t_loop t_loop_curr];
    t_subloop = [t_subloop t_subloop_curr];
end

for curr_xsg = 1:length(xsg_filename)
    
    % get loop number, find the tiff files, check for all files present
    curr_loop = str2num(xsg_filename{curr_xsg}(end-7:end-4));
    tiff_curr_loop = find(t_loop == curr_loop);
    if all(t_subloop(tiff_curr_loop) ~= [1:8])
        disp(['Don''t have 1-8 for loop' num2str(curr_loop) ',skipping']);
        continue
    end
    tiff_filename_curr = tiff_filename(tiff_curr_loop);
    
    % load in tiff images, columns are frames
    im_temp = [];
    im = [];
    
    for i = 1:length(tiff_filename_curr);
        img_filename = [tiff_path tiff_filename_curr{i}];
        imageinfo=imfinfo(img_filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        clear im_temp_long
        clear im_temp
        
        disp('Loading file (native format)....')
        parfor loadframe = 1:numframes
            im_temp = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
            im_temp_long(:,loadframe) = im_temp(:);
        end
        
        im = [im im_temp_long];
        disp(['Done ' num2str(i) '/' num2str(length(tiff_filename_curr))])
    end
    clear im_temp_long
    
    % load corresponding xsg
    xsg = load([xsg_path xsg_filename{curr_xsg}],'-MAT');
    lever_force = xsg.data.acquirer.trace_2;
    
    X = M;
    Y = N;
    numframes = size(im,2);
    % to make resampling work: numframes must be even number
    % if it's not, cut out one frame at the end
    if mod(numframes,2) ~= 0
        im(:,end) = [];
        numframes = size(im,2);
    end
    lever_force_resample = resample(lever_force,numframes,length(lever_force));
    
    % Find frame for each reward, grab surrounding frames and resampled force
    
    % Parse header, get framerate
    img_info = imageinfo(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    % sample rate for xsg
    xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
    
    % Get trials in seconds since started
    curr_trial_list = [];
    trial_list = [];
    trial_stitch = []; % unused at the moment since looking only at single xsg
    curr_trial_list = AP_ReadBitCode([xsg_path xsg_filename{curr_xsg}]);
    if ~isempty(curr_trial_list)
        % this is commented out for now b/c looking at 1 xsg at a time
        %         curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
        trial_list = [trial_list; curr_trial_list];
    end
    
    % Get all trial offsets in seconds
    xsg_bhv_offset = zeros(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            xsg_bhv_offset(curr_trial) = NaN;
        end
    end
    
    % Get lever presses when rewarded
    reward_trials_num = find(reward_trials == 1);
    lever_surround_sec = 1;
    lever_surround_sample = round((lever_surround_sec*framerate)/2);
    lever_force_resample_reward = [];
    lever_force_resample_reward_parsed = {};
    force_point = nan(num_trials,1);
    im_reward = [];
    for curr_reward_trial = reward_trials_num';
        % Skip this if trial start wasn't recorded
        if ~isnan(xsg_bhv_offset(curr_reward_trial)) & ...
                ~ismember(curr_reward_trial,trial_stitch)
            % Current xsg trial time
            xsg_trial_indx = find(trial_list(:,2) == curr_reward_trial,1);
            lever_frame_bhv = (parsed_events{curr_reward_trial}.pokes.L(lever_cue_indx{curr_reward_trial},1) - all_trial_start{curr_reward_trial}(1))*(framerate);
            % xsg lever time = relative lever sample trial start + xsg trial start
            lever_frame_xsg = round(lever_frame_bhv + trial_list(xsg_trial_indx,1)*framerate);
            % check that the frames in the range exist, take what's there
            clear frame_range
            frame_range = lever_frame_xsg-lever_surround_sample:lever_frame_xsg+lever_surround_sample;
            frame_range(frame_range < 1) = [];
            frame_range(frame_range > size(im,2)) = [];
            % grab surrounding frames
            lever_force_resample_reward = [lever_force_resample_reward ...
                lever_force_resample(frame_range)'];
            curr_num_rewards = length(lever_force_resample_reward_parsed);
            lever_force_resample_reward_parsed{curr_num_rewards+1} = ...
                lever_force_resample(frame_range)';
            im_reward = [im_reward ...
                im(:,frame_range)];
            force_point(curr_reward_trial) = lever_frame_xsg;
        end
    end
    
    cc_img = zeros(Y*X,1);
    cc_img_reward = zeros(Y*X,1);
    
    % concatenate everything so far
    im_concat = [im_concat im];
    im_reward_concat = [im_reward_concat im_reward];
    lever_force_resample_concat = [lever_force_resample_concat ...
        lever_force_resample];
    lever_force_resample_reward_concat = [lever_force_resample_reward_concat ...
        lever_force_resample_reward];
end

lever_force_resample_concat = lever_force_resample_concat(:);

% get full correlation for day
parfor px_curr = 1:(X*Y);
    % slice the image variable for parallel
    im_slice = im_concat(px_curr,:);
    im_reward_slice = im_reward_concat(px_curr,:);
    
    % skip first and last 10 frames to avoid resampling edge effects
    im_slice = double(im_slice(10:end-10));
    im_reward_slice = double(im_reward_slice(10:end-10));
    
    curr_cc = corrcoef(im_slice,lever_force_resample_concat(10:end-10));
    curr_cc_reward = corrcoef(im_reward_slice,lever_force_resample_reward_concat(10:end-10));
    cc_img(px_curr) = curr_cc(2);
    cc_img_reward(px_curr) = curr_cc_reward(2);
    disp([num2str(px_curr) '/' num2str(X*Y)]);
end

curr_date = tiff_filename{1}(1:6);
savename = ['/usr/local/lab/People/Andy/Data/CorrData/' curr_date '_' xsg_filename{curr_xsg}(1:end-4) '_FullDay'];
save(savename,'cc_img','cc_img_reward','im_reward','lever_force_resample_reward','xsg_filename','tiff_filename');
disp(['Saved ' savename]);

toc
disp(['Finished analysis']);


%% ROIs based on correlations
% 1) Crop edges for motion 2) Gaussian filter, 3) threshold,
% 4) largest to smallest cc, include adjoining 95 (?) corr and exclude 
% 95>x>50(?) for processing 5) exclude if under pixel size limit

cc_mean_reward_reshape = reshape(cc_img_reward,512,512);

% crop edges - motion makes artificial correlations
edge_crop = 15;
cc_mean_reward_reshape_cropped = cc_mean_reward_reshape ...
    (edge_crop:end-edge_crop,edge_crop:end-edge_crop);

% gaussian smooth correlated images
hsize = 10;
sigma = 5;
h = fspecial('gaussian', hsize, sigma);
clear cc_smooth;
cc_smooth = imfilter(cc_mean_reward_reshape_cropped,h);

% threshold
clear cc_smooth_thresh;
cc_smooth_thresh = abs(cc_smooth) > 0.05;

% get values from threshold, put border back on
clear small_cc_thresh_value cc_thresh_value
small_cc_thresh_value = cc_mean_reward_reshape_cropped;
small_cc_thresh_value(~cc_smooth_thresh) = NaN;
cc_thresh_value = nan(size(small_cc_thresh_value) + 2*edge_crop - 1);
cc_thresh_value(edge_crop:end-edge_crop,edge_crop:end-edge_crop) = ...
    small_cc_thresh_value;

% largest to smallest cc, find 95% correlated

% first, index eligible pixles
px_eligible = find(~isnan(cc_thresh_value));
% attach absolute cc values to those
px_eligible(:,2) = abs(cc_thresh_value(px_eligible));
% sort by absolute value
px_eligible_sorted = sortrows(px_eligible,2);
% make it descending
px_eligible_sorted = flipud(px_eligible_sorted);

% loop through all eligible pixels to make
curr_px_eligible_indx = 1;
while curr_px_eligible_indx < length(px_eligible_sorted+1)
    curr_px_eligible = px_eligible_sorted(curr_px_eligible_indx,1);
    curr_px_trace = double(im_concat(curr_px_eligible,:));
    % correlate that pixel with all other pixels
    curr_cc_img = [];
    parfor curr_cc_px = 1:size(im_concat,1)
        curr_cc_px_trace = double(im_concat(curr_cc_px,:));
        curr_cc = corrcoef(curr_px_trace, curr_cc_px_trace);
        curr_cc_img(curr_cc_px) = curr_cc(2);
        disp(num2str(curr_cc_px));
    end
    curr_cc_img_thresh = find (curr_cc_img > 0.95);
end








