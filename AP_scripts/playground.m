%% post to alyx

onLoad; % initializes missing-http

myAlyx = alyx.getToken([], 'andy', 'xpr1mnt');

use_date = datestr(now,'yyyy-mm-dd');
use_time = datestr(now,'HH:MM');

clear d
d.user = 'andy';

animals = {'AP031'};

water = [1.2];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.water_administered = water(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWater = alyx.postData(myAlyx, 'water-administrations', d);
end

% weight = [26.7];
% for curr_animal = 1:length(animals)
%     d.subject = animals{curr_animal};
%     d.weight = weight(curr_animal);
%     d.date_time = sprintf('%sT%s',use_date,use_time);
%     newWeight = alyx.postData(myAlyx, 'weighings', d);
% end



%% loop through days, do regression

animal = 'AP017';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

experiment = 1;
for curr_day_idx = 9:4:length(days)
    
    day = days{curr_day_idx};
    AP_load_experiment;
    
    % Skip the first n seconds to do this
    skip_seconds = 10;
    use_frames = (frame_t > skip_seconds);
    
    % Make choiceworld event traces
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    choiceworld_event_trace = [];
    
    %use_contrasts = [0,0.125,0.25,0.5,1];
    use_contrasts = [0.5,1];
    use_trialSides = [-1,1];
    use_hitValues = [1];
    
    for trialSide_idx = 1:length(use_trialSides)
        
        curr_trialSide = use_trialSides(trialSide_idx);
        
        use_trials = ismember(signals_events.trialContrastValues,use_contrasts) &  ...
            signals_events.trialSideValues == curr_trialSide & ...
            ismember(signals_events.hitValues,use_hitValues);
        align_times = signals_events.interactiveOnTimes(use_trials(1:length(signals_events.interactiveOnTimes)))';
        choiceworld_event_trace = [choiceworld_event_trace;histcounts(align_times,frame_edges)];
        
    end
    
    %choiceworld_event_trace = [choiceworld_event_trace;licking_trace];
    %choiceworld_event_trace = [choiceworld_event_trace;wheel_speed];
    
    %choiceworld_event_trace = zscore(choiceworld_event_trace,[],2);
    
    use_svs = 1:50;
    kernel_frames = -7:35;
    lambda = 0;
    zs = false;
    cvfold = 5;
    
    [k,predicted_fluor,explained_var] = ...
        AP_regresskernel(choiceworld_event_trace(:,use_frames), ...
        fVdf(use_svs,use_frames), ...
        kernel_frames,lambda,zs,cvfold);
    
    % Reshape kernel and convert to pixel space
    k_r = permute(reshape(k,size(choiceworld_event_trace,1),length(kernel_frames),length(use_svs)),[3,2,1]);
    
    r_px = zeros(size(U,1),size(U,2),size(k_r,2),size(k_r,3),'single');
    for curr_event = 1:size(k_r,3);
        r_px(:,:,:,curr_event) = svdFrameReconstruct(Udf(:,:,use_svs),k_r(:,:,curr_event));
    end
    
    a = diff(r_px,[],4);
    
    AP_image_scroll(a,kernel_frames/samplerate);
    caxis([prctile(a(:),[1,99])]*4);
    truesize
    
end

%% Regress position stim to fluorescence

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Make choiceworld event traces
frame_edges = [frame_t,frame_t(end)+1/samplerate];
choiceworld_event_trace = [];

use_azimuths = unique(signals_events.trialAzimuthValues);

for curr_condition = 1:length(use_azimuths)
    use_trials = ismember(signals_events.trialAzimuthValues,use_azimuths(curr_condition));
    align_times = signals_events.stimOnTimes(use_trials(1:length(signals_events.stimOnTimes)))';
    choiceworld_event_trace = [choiceworld_event_trace;histcounts(align_times,frame_edges)];
end

use_svs = 1:50;
kernel_frames = -7:35;
lambda = 0;
zs = false;
cvfold = 5;

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel(choiceworld_event_trace(:,use_frames), ...
    fVdf(use_svs,use_frames), ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
k_r = permute(reshape(k,size(choiceworld_event_trace,1),length(kernel_frames),length(use_svs)),[3,2,1]);

r_px = zeros(size(U,1),size(U,2),size(k_r,2),size(k_r,3),'single');
for curr_event = 1:size(k_r,3);
    r_px(:,:,:,curr_event) = svdFrameReconstruct(Udf(:,:,use_svs),k_r(:,:,curr_event));
end

AP_image_scroll(r_px,kernel_frames/samplerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize


%% Decode stim position

% Get orientation decoding as a function of SVs
window_length = 0.5;
start_window = 0.2;

use_svs = 200;

stim_matrix_distance = nan(length(use_svs),1);
correct_decoding = nan(length(use_svs),1);

use_azimuths = unique(signals_events.trialAzimuthValues);

confusion_mat = nan(length(use_azimuths),length(use_azimuths),length(use_svs));

samplerate = 1./median(diff(frame_t));
surround_samplerate = 1/(samplerate*1);
surround_time = start_window:surround_samplerate: ...
    start_window+window_length;

align_surround_times = bsxfun(@plus, signals_events.stimOnTimes', surround_time);

% don't use stim that are out of frame time bounds
use_trials = ~any(align_surround_times > frame_t(end),2);

avg_response_v = permute(nanmean(interp1(frame_t,fV',align_surround_times),2),[3,1,2]);

avg_response_v = avg_response_v(1:use_svs,:);

stim_order = zeros(length(use_azimuths),length(signals_events.stimOnTimes));
for curr_condition = 1:length(use_azimuths)
    stim_order(curr_condition,signals_events.trialAzimuthValues == use_azimuths(curr_condition)) = 1;
end

lambda = 1e7;
zs = false;
cvfold = 5;

%stim_order = shake(stim_order,1);
stim_order = stim_order(:,use_trials);
[~,stim_order_idx] = max(stim_order,[],1);
avg_response_v = avg_response_v(:,use_trials);

[k,predicted_stim,explained_var] = ...
    AP_regresskernel(avg_response_v,stim_order,0,lambda,zs,cvfold);

max_likely_stim = bsxfun(@eq,predicted_stim,max(predicted_stim,[],1));
correct_decoding = sum(max_likely_stim(:) & stim_order(:))./sum(stim_order(:));

[~,predicted_stims] = max(predicted_stim,[],1);
confusion_mat(:,:) = bsxfun(@rdivide,confusionmat(stim_order_idx,predicted_stims),sum(stim_order,2));

stim_matrix_distance = pdist2(predicted_stim(:)',stim_order(:)');


% Get orientation decoding as a function of SVs
window_length = 0.5;
start_window = [0:0.5:3];

use_svs = 1000;

stim_matrix_distance = nan(length(start_window),1);
correct_decoding = nan(length(start_window),1);
confusion_mat = nan(length(use_azimuths),length(use_azimuths),length(start_window));

for curr_start_window = 1:length(start_window)
    
    use_azimuths = unique(signals_events.trialAzimuthValues);
    
    samplerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(samplerate*1);
    surround_time = start_window(curr_start_window):surround_samplerate: ...
        start_window(curr_start_window)+window_length;
    
    align_surround_times = bsxfun(@plus, signals_events.stimOnTimes', surround_time);
    
    % don't use stim that are out of frame time bounds
    use_trials = ~any(align_surround_times > frame_t(end),2);
    
    avg_response_v = permute(nanmean(interp1(frame_t,fV',align_surround_times),2),[3,1,2]);
    
    avg_response_v = avg_response_v(1:use_svs,:);
    
    stim_order = zeros(length(use_azimuths),length(signals_events.stimOnTimes));
    for curr_condition = 1:length(use_azimuths)
        stim_order(curr_condition,signals_events.trialAzimuthValues == use_azimuths(curr_condition)) = 1;
    end
    
    lambda = 0;
    zs = false;
    cvfold = 5;
    
    %stim_order = shake(stim_order(:,use_trials),1);
    [~,stim_order_idx] = max(stim_order,[],1);
    avg_response_v = avg_response_v(:,use_trials);
    
    [k,predicted_stim,explained_var] = ...
        AP_regresskernel(avg_response_v,stim_order,0,lambda,zs,cvfold);
    
    max_likely_stim = bsxfun(@eq,predicted_stim,max(predicted_stim,[],1));
    correct_decoding(curr_start_window) = sum(max_likely_stim(:) & stim_order(:))./sum(stim_order(:));
    
    [~,predicted_stims] = max(predicted_stim,[],1);
    confusion_mat(:,:,curr_start_window) = bsxfun(@rdivide,confusionmat(stim_order_idx,predicted_stims),sum(stim_order,2));
    
    stim_matrix_distance(curr_start_window) = pdist2(predicted_stim(:)',stim_order(:)');
    
    disp(curr_start_window)
    
end


%% Load AP017 choiceworld days, regress stim and decisions

% Get choiceworld days
animal = 'AP017';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};
vanillaChoiceworld_days = false(size(days));
for curr_day = 1:length(days)
    
    day = days{curr_day};
    % In the event of multiple experiments, use the last one
    expDay_dir = dir([expInfo_path filesep days{curr_day}]);
    exp_nums = [expDay_dir(3:end).name];
    
    [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(end),'block');
    
    if ~block_exists
        continue
    end
    
    % Load the block file
    load(block_filename)
    
    [~,expDef] = fileparts(block.expDef);
    vanillaChoiceworld_days(curr_day) = strcmp(expDef,'vanillaChoiceworld');
    
end

use_days = days(vanillaChoiceworld_days);
k_px_diff_choice = cell(length(use_days),1);
k_px_diff_stim = cell(length(use_days),1);
for curr_day = 1:length(use_days);
    
    animal = 'AP017';
    day = use_days{curr_day};
    
    expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
    expDay_dir = dir([expInfo_path filesep day]);
    exp_nums = [expDay_dir(3:end).name];
    experiment = exp_nums(end);
    
    AP_load_experiment
    
    if ~exist('fV','var')
        continue
    end
    
    %%%%%%%%%%%%%%%%%
    % CHOICE DECODING
    %%%%%%%%%%%%%%%%%
    
    % Get left/right choice trials (of chosen contrasts)
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    
    left_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.125,0.06,0]) &  ...
        ((ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[1])) | ...
        (ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[0]))) & ...
        ~signals_events.repeatTrialValues;
    
    right_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.125,0.06,0]) &  ...
        ((ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[0])) | ...
        (ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[1]))) & ...
        ~signals_events.repeatTrialValues;
    
    % Get photodiode flips closest to stim presentations
    photodiode_name = 'photoDiode';
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
    photodiode_flip_times = Timeline.rawDAQTimestamps( ...
        find((Timeline.rawDAQData(1:end-1,photodiode_idx) <= 2) ~= ...
        (Timeline.rawDAQData(2:end,photodiode_idx) <= 2)) + 1);
    
    [~,closest_stimOn_photodiode] = ...
        arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
        photodiode_flip_times)), ...
        1:length(signals_events.stimOnTimes));
    stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
    
    % Get time to first movement after stim onset
    wheel_move_thresh = 1;
    wheel_move_start_times = frame_t( ...
        find((wheel_velocity(1:end-1) <= wheel_move_thresh) & ...
        (wheel_velocity(2:end) > wheel_move_thresh)) + 1);
    stimOn_move_times = ...
        arrayfun(@(x) ...
        wheel_move_start_times(find(wheel_move_start_times > stimOn_times(x),1)) - ...
        stimOn_times(x),1:length(stimOn_times),'uni',false);
    stimOn_move_times(cellfun(@isempty,stimOn_move_times)) = {NaN};
    stimOn_move_times = cell2mat(stimOn_move_times);
    
    % % Restrict trials to non-move before interactive on
    left_trials = left_trials & stimOn_move_times >= 0.5;
    right_trials = right_trials & stimOn_move_times >= 0.5;
    
    % Fix the parameters
    use_window = [-0.2,0.5];
    
    samplerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(samplerate*1);
    surround_time = use_window(1):surround_samplerate: ...
        use_window(2);
    
    align_surround_times_left = bsxfun(@plus, stimOn_times(left_trials)', surround_time);
    align_surround_times_right = bsxfun(@plus, stimOn_times(right_trials)', surround_time);
    
    fV_align_left = interp1(frame_t,fV',align_surround_times_left);
    fV_align_right = interp1(frame_t,fV',align_surround_times_right);
    
    stim_order = [[ones(1,sum(left_trials)),zeros(1,sum(right_trials))]; ...
        [zeros(1,sum(left_trials)),ones(1,sum(right_trials))]];
    
    use_svs = 100;
    lambda = 1e6;
    zs = false;
    cvfold = 5;
    
    fV_align_all = [reshape(permute(fV_align_left(:,:,1:use_svs),[2,3,1]),[],sum(left_trials)), ...
        reshape(permute(fV_align_right(:,:,1:use_svs),[2,3,1]),[],sum(right_trials))];
    
    [k,predicted_stim,explained_var] = ...
        AP_regresskernel(fV_align_all,stim_order,0,lambda,zs,cvfold);
    
    max_likely_stim = bsxfun(@eq,predicted_stim,max(predicted_stim,[],1));
    correct_decoding = sum(max_likely_stim(:) & stim_order(:))./sum(stim_order(:));
    
    k2 = reshape(k,[],use_svs,2);
    k_px_l = svdFrameReconstruct(U(:,:,1:use_svs),k2(:,:,1)');
    k_px_r = svdFrameReconstruct(U(:,:,1:use_svs),k2(:,:,2)');
    
    baseline_time = find(surround_time < 0,1,'last');
    k_px_l_norm = bsxfun(@minus,k_px_l,nanmean(k_px_l(:,:,1:baseline_time),3));
    k_px_r_norm = bsxfun(@minus,k_px_r,nanmean(k_px_r(:,:,1:baseline_time),3));
    k_px_diff = k_px_l_norm - k_px_r_norm;
    
    k_px_diff_choice{curr_day} = k_px_diff;
    
    %%%%%%%%%%%%%%%%%
    % STIMULUS DECODING
    %%%%%%%%%%%%%%%%%
    
    
    % Get left/right choice trials (of chosen contrasts)
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    
    left_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.125,0.06,0]) &  ...
        ((ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[1])) | ...
        (ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[0]))) & ...
        ~signals_events.repeatTrialValues;
    
    right_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.125,0.06,0]) &  ...
        ((ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[0])) | ...
        (ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[1]))) & ...
        ~signals_events.repeatTrialValues;
    
    % Get photodiode flips closest to stim presentations
    photodiode_name = 'photoDiode';
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
    photodiode_flip_times = Timeline.rawDAQTimestamps( ...
        find((Timeline.rawDAQData(1:end-1,photodiode_idx) <= 2) ~= ...
        (Timeline.rawDAQData(2:end,photodiode_idx) <= 2)) + 1);
    
    [~,closest_stimOn_photodiode] = ...
        arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
        photodiode_flip_times)), ...
        1:length(signals_events.stimOnTimes));
    stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
    
    % Get time to first movement after stim onset
    wheel_move_thresh = 1;
    wheel_move_start_times = frame_t( ...
        find((wheel_velocity(1:end-1) <= wheel_move_thresh) & ...
        (wheel_velocity(2:end) > wheel_move_thresh)) + 1);
    stimOn_move_times = ...
        arrayfun(@(x) ...
        wheel_move_start_times(find(wheel_move_start_times > stimOn_times(x),1)) - ...
        stimOn_times(x),1:length(stimOn_times),'uni',false);
    stimOn_move_times(cellfun(@isempty,stimOn_move_times)) = {NaN};
    stimOn_move_times = cell2mat(stimOn_move_times);
    
    % % Restrict trials to non-move before interactive on
    left_trials = left_trials & stimOn_move_times >= 0.5;
    right_trials = right_trials & stimOn_move_times >= 0.5;
    
    % Fix the parameters
    use_window = [-0.2,0.5];
    
    samplerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(samplerate*1);
    surround_time = use_window(1):surround_samplerate: ...
        use_window(2);
    
    align_surround_times_left = bsxfun(@plus, stimOn_times(left_trials)', surround_time);
    align_surround_times_right = bsxfun(@plus, stimOn_times(right_trials)', surround_time);
    
    fV_align_left = interp1(frame_t,fV',align_surround_times_left);
    fV_align_right = interp1(frame_t,fV',align_surround_times_right);
    
    stim_order = [[ones(1,sum(left_trials)),zeros(1,sum(right_trials))]; ...
        [zeros(1,sum(left_trials)),ones(1,sum(right_trials))]];
    
    use_svs = 100;
    lambda = 1e6;
    zs = false;
    cvfold = 5;
    
    fV_align_all = [reshape(permute(fV_align_left(:,:,1:use_svs),[2,3,1]),[],sum(left_trials)), ...
        reshape(permute(fV_align_right(:,:,1:use_svs),[2,3,1]),[],sum(right_trials))];
    
    [k,predicted_stim,explained_var] = ...
        AP_regresskernel(fV_align_all,stim_order,0,lambda,zs,cvfold);
    
    max_likely_stim = bsxfun(@eq,predicted_stim,max(predicted_stim,[],1));
    correct_decoding = sum(max_likely_stim(:) & stim_order(:))./sum(stim_order(:));
    
    k2 = reshape(k,[],use_svs,2);
    k_px_l = svdFrameReconstruct(U(:,:,1:use_svs),k2(:,:,1)');
    k_px_r = svdFrameReconstruct(U(:,:,1:use_svs),k2(:,:,2)');
    
    baseline_time = find(surround_time < 0,1,'last');
    k_px_l_norm = bsxfun(@minus,k_px_l,nanmean(k_px_l(:,:,1:baseline_time),3));
    k_px_r_norm = bsxfun(@minus,k_px_r,nanmean(k_px_r(:,:,1:baseline_time),3));
    k_px_diff = k_px_l_norm - k_px_r_norm;
    
    k_px_diff_stim{curr_day} = k_px_diff;
    
    disp(curr_day)
    clearvars -except curr_day use_days k_px_diff_choice k_px_diff_stim
    
end


%% After alignment: align the kernels

k_px_diff_stim_reg = cell(length(use_days),1);
k_px_diff_choice_reg = cell(length(use_days),1);

k_px_diff_stim_reg{ref_im_num} = k_px_diff_stim{ref_im_num};
k_px_diff_choice_reg{ref_im_num} = k_px_diff_choice{ref_im_num};

for curr_day = setdiff(1:length(use_days),[4,ref_im_num]);
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    k_px_diff_stim_reg{curr_day} = imwarp(k_px_diff_stim{curr_day},tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
    
    k_px_diff_choice_reg{curr_day} = imwarp(k_px_diff_choice{curr_day},tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
end


%% Plot correlation of kernel with average kernel

avg_stim_k = nanmean(cat(4,k_px_diff_stim_reg{:}),4);
avg_choice_k = nanmean(cat(4,k_px_diff_choice_reg{:}),4);

stim_k_corr = cellfun(@(x) corrcoef(x(:),avg_stim_k(:)),k_px_diff_stim_reg(setdiff(1:length(use_days),4)),'uni',false);
stim_k_corr = cellfun(@(x) x(2),stim_k_corr);

choice_k_corr = cellfun(@(x) corrcoef(x(:),avg_choice_k(:)),k_px_diff_choice_reg(setdiff(1:length(use_days),4)),'uni',false);
choice_k_corr = cellfun(@(x) x(2),choice_k_corr);

choice_stim_k_corr = cellfun(@(x) corrcoef(x(:),avg_stim_k(:)),k_px_diff_choice_reg(setdiff(1:length(use_days),4)),'uni',false);
choice_stim_k_corr = cellfun(@(x) x(2),choice_stim_k_corr);

figure;
subplot(1,2,1); hold on;
plot(setdiff(1:length(use_days),4),stim_k_corr,'k','linewidth',2);
plot(setdiff(1:length(use_days),4),choice_k_corr,'r','linewidth',2);
plot(setdiff(1:length(use_days),4),choice_stim_k_corr,'b','linewidth',2);

legend({'Stim kernel','Choice kernel','Choice -> Stim'})
xlabel('Day');
ylabel('Correlation of kernel with average kernel');

% Performance of 0.5/1
go_right_correct = 1-(sum(bhv.go_left_trials(:,1:2),2)./sum(bhv.n_trials_condition(:,1:2),2));
go_left_correct = sum(bhv.go_left_trials(:,end-1:end),2)./sum(bhv.n_trials_condition(:,end-1:end),2);
subplot(1,2,2); hold on;
plot(go_right_correct,'k','linewidth',2);
plot(go_left_correct,'r','linewidth',2);
legend({'Go right','Go left'})
xlabel('Day');
ylabel('Fraction correct (0.5/1 contrast)');

%% %%%%%%%%%%%%%%%%% LEVER

%% Load AP015 lever days, regress stim/lever/reward

% Get choiceworld days
animal = 'AP015';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};
lever_days = days(end-5:end);

use_days = lever_days;
k_px = cell(length(use_days),1);
for curr_day = 1:length(use_days);
    
    animal = 'AP015';
    day = use_days{curr_day};
    
    expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
    expDay_dir = dir([expInfo_path filesep day]);
    exp_nums = [expDay_dir(3:end).name];
    experiment = exp_nums(end);
    
    AP_load_experiment
    
    % Skip the first n seconds to do this
    skip_seconds = 10;
    use_frames = (frame_t > skip_seconds);
    
    % Make choiceworld event traces
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    signals_event_trace = [];
    
    align_times = signals_events.stimOnTimes';
    signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];
    
    align_times = signals_events.lever_r_flipTimes(signals_events.lever_r_flipValues == -1)';
    signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];
    
    align_times = signals_events.lever_r_flipTimes(signals_events.lever_r_flipValues == 1)';
    signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];
    
    align_times = signals_events.totalWaterTimes';
    signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];
    
    use_svs = 1:50;
    kernel_frames = -35*5:35*5;
    lambda = 0;
    zs = false;
    cvfold = 5;
    
    [k,predicted_fluor,explained_var] = ...
        AP_regresskernel(signals_event_trace(:,use_frames), ...
        fV(use_svs,use_frames), ...
        kernel_frames,lambda,zs,cvfold);
    
    % Reshape kernel and convert to pixel space
    k_r = permute(reshape(k,size(signals_event_trace,1),length(kernel_frames),length(use_svs)),[3,2,1]);
    
    r_px = zeros(size(U,1),size(U,2),size(k_r,2),size(k_r,3),'single');
    for curr_event = 1:size(k_r,3);
        r_px(:,:,:,curr_event) = svdFrameReconstruct(U(:,:,use_svs),k_r(:,:,curr_event));
    end
    
    k_px{curr_day} = r_px;
    
    disp(curr_day)
    clearvars -except curr_day use_days k_px
    
end














%%

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 1300 & templateDepths < 2000)));

align_times = stim_times_grid{4,30};

raster_window = [-1,2];

psth_bin_size = 0.001;
[psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
    use_spikes,align_times, ...
    raster_window, psth_bin_size);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(psth,smWin,'same');
figure; hold on;
plot(bins(20:end-20),psth_smooth(20:end-20)','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Population spikes');
xlabel('Time from stim onset')


% PSTH by stim condition
a = cell(numel(stim_times_grid),1);
stim_psth = nan(size(stim_times_grid));
for curr_stim_idx = 1:numel(stim_times_grid);
    
    [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
        use_spikes,stim_times_grid{curr_stim_idx}, ...
        raster_window, psth_bin_size);
    
    stim_psth(curr_stim_idx) = max(psth);
    a{curr_stim_idx} = psth;
    
end



use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 1300 & templateDepths < 2000)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/samplerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

% Interpolate stim screen for all times
stim_screen_interp = abs(single(interp1(stim_times,reshape(stim_screen,[],length(stim_times))',frame_t,'nearest')'));
stim_screen_interp(isnan(stim_screen_interp)) = 0;

% Use only the onsets of stim
stim_screen_interp = single([zeros(size(stim_screen_interp,1),1),diff(stim_screen_interp,[],2)] == 1);

fluor_kernel_frames = -8:-2;
lambda = 1000;

k = AP_regresskernel(frame_spikes, ...
    stim_screen_interp,fluor_kernel_frames,lambda,true,5);


k = AP_regresskernel(stim_screen_interp,frame_spikes, ...
    [-2:8],lambda,true,5);



%%



% load dataset
dataset = '1';
calcium_train = csvread([dataset '.train.calcium.csv']);
spike_train = csvread([dataset '.train.spikes.csv']);

% plot example neuron
figure
n = 5;         % index of neuron
t = (0:length(calcium_train(:,n))-1)/100;     % 100Hz sampling rate, each bin 10 ms
plot(t,zscore(calcium_train(:,n)),'k')
hold on
plot(t,spike_train(:,5)-2,'r')

xlim([t(1) 400])
ylim([-2.5 7])

xlabel('Time (s)')
ylabel('Fluorescence / Spike rate')

legend('calcium','spikes')

% saving predictions
spike_pred = calcium_train;     % fill spike_pred with your own predictions
csvwrite([dataset '.train.pred.csv'], spike_pred);


%% %%%%%%%%%%%%%%%%%%%%%%%%% TEMPLATE TO DO STUFF IN BATCH
% NOTE! use this now: experiments = AP_find_experiments(animal,protocol);


animal = 'AP026';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};


batch_vars = struct;
for curr_day = 1:length(use_days);
    
    day = use_days{curr_day};
    
    expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
    expDay_dir = dir([expInfo_path filesep day]);
    experiment = 1;
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    stimIDs = signals_events.trialSideValues(1:end-1).*signals_events.trialContrastValues(1:end-1);
    stim_onsets = stimOn_times(1:length(stimIDs));
    
    % Average (single frame) responses to stimuli
    surround_window = [0.2,0.3];
    samplerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(samplerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    
    % Make choiceworld event trace
    use_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25]) &  ...
        ismember(signals_events.trialSideValues,[1]) & ...
        ismember(signals_events.hitValues,[1]);
    align_times = signals_events.stimOnTimes(use_trials(1:length(signals_events.stimOnTimes)))';
    
    align_times(align_times + surround_time(1) < frame_t(2) | ...
        align_times + surround_time(2) > frame_t(end)) = [];
    
    align_surround_times = bsxfun(@plus, align_times, surround_time);
    peri_stim_v = nanmean(permute(nanmean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]),2);
    
    surround_im = svdFrameReconstruct(U,peri_stim_v);
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    disp(curr_day)
    clearvars -except use_days curr_day animal batch_vars
    
end


%% After alignment: align the kernels

batch_vars_reg = batch_vars;

for curr_day = setdiff(1:length(use_days),ref_im_num);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.surround_im{curr_day};
    curr_im(isnan(curr_im)) = 0;
    batch_vars_reg.surround_im{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
    
end


a = cat(5,batch_vars_reg.im_stim{:});
b = nanmean(a,5);
AP_image_scroll(b,kernel_frames/35)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Behavior in batch


%% lick task over days

animal = 'AP022';
% just get all days for now (eventually choose, supply date range, etc)
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

% Initialize the behavior structure
bhv = struct;

for curr_day = 1:length(days)
    day = days{curr_day};
    experiment = 1;
    AP_load_experiment;
    
    lick_idx = strcmp({Timeline.hw.inputs.name}, 'beamLickDetector');
    lick_trace = Timeline.rawDAQData(:,lick_idx) > 2.5;
    lick_times = Timeline.rawDAQTimestamps(find(lick_trace(2:end) & ~lick_trace(1:end-1))+1);
    
    epoch_times = reshape([stimOn_times,stimOff_times]',[],1);
    epoch_hit = histcounts(lick_times,epoch_times) > 0;
    
    orientations = signals_events.trialOrientationValues;
    stim_hit = epoch_hit(1:2:end)';
    
    stim_hit_licktime_cell = arrayfun(@(x) lick_times(find(...
        lick_times >= stimOn_times(x),1)),1:length(stimOn_times),'uni',false);
    stim_hit_licktime = nan(size(stim_hit));
    stim_hit_licktime(stim_hit) = [stim_hit_licktime_cell{stim_hit}];
    stim_lick_delay = stim_hit_licktime - stimOn_times;
    
    % Get lick aligned to stim hit
    surround_interval = [-2,2];
    surround_time = surround_interval(1):Timeline.hw.samplingInterval:surround_interval(2);
    water_surround_lick = bsxfun(@plus,stim_hit_licktime(~isnan(stim_hit_licktime)), ...
        surround_time);
    water_aligned_lick = interp1(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,lick_idx), ...
        water_surround_lick);
    
    bhv(curr_day).orientations = orientations;
    bhv(curr_day).stim_hit = stim_hit;
    bhv(curr_day).stim_lick_delay = stim_lick_delay;
    bhv(curr_day).water_aligned_lick = water_aligned_lick;
    
end

figure;

subplot(1,3,1); hold on;
plotSpread({bhv.stim_lick_delay});
errorbar(cellfun(@nanmedian,{bhv.stim_lick_delay}), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),{bhv.stim_lick_delay}),'r')
xlabel('Day');
ylabel('Time from stim to lick')

subplot(2,3,2); hold on; set(gca,'ColorOrder',copper(length(days)));
frac_stim = arrayfun(@(x) smooth(+bhv(x).stim_hit,20),1:length(days),'uni',false);
for i = 1:length(days)
    plot(frac_stim{i},'linewidth',2);
end
ylim([0 1]);
xlabel('Trial');
ylabel('Fraction hit trials');

subplot(2,3,5);
plot(cellfun(@median,frac_stim),'k','linewidth',2);
xlabel('Day');
ylabel('Median smoothed fraction hit');
ylim([0 1]);

subplot(1,3,3)
day_trials_plotted = arrayfun(@(x) size(bhv(x).water_aligned_lick,1),1:length(days));
imagesc(surround_time,1:sum(day_trials_plotted),1-vertcat(bhv.water_aligned_lick));
colormap(gray);
cumulative_trials = cumsum(day_trials_plotted);
for i = 1:length(days)
    line(xlim,[cumulative_trials(i),cumulative_trials(i)],'color','r');
end
ylabel('Trial');
xlabel('Time from rewarded lick');



%% Testing predictability measures

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 2000 & templateDepths < 3000)));
frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/samplerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

frame_spikes = smooth(frame_spikes,35)';
%frame_spikes = frame_spikes - smooth(frame_spikes,5000,'loess')';

% Thin spikes at random a bunch of times
num_reps = 10;
thin_spikes_frac = 0.5;

frame_spikes_thinned = zeros(num_reps,length(frame_t),'single');
for curr_rep = 1:num_reps
    curr_thinned_spikes = use_spikes;
    curr_thinned_spikes( ...
        randperm(length(curr_thinned_spikes), ...
        round(length(curr_thinned_spikes)*thin_spikes_frac))) = [];
    
    frame_spikes_thinned(curr_rep,:) = histcounts(curr_thinned_spikes,frame_edges);
end

frame_spikes_use = [frame_spikes;frame_spikes_thinned];



use_svs = 1:50;
kernel_frames = -20:20;
lambda = 1e6;
zs = [false,false];
cvfold = 5;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fV(use_svs,use_frames), ...
    frame_spikes_use(:,use_frames),kernel_frames,lambda,zs,cvfold);




% Evaluating prediction with Q (or L?)
% Q = -1/2(integral)f(t)^2 dt + (sum over s) f(ts)
% L = (sum over s) log f(ts) - (integral) f df

%f_real = zscore(frame_spikes(:,use_frames),[],2);
f_real = frame_spikes(:,use_frames);
f_pred = predicted_spikes;
T = frame_t(find(use_frames,1,'last')) - frame_t(find(use_frames,1,'first'));

Q = (-1/2)*(sum(f_pred.^2./T,2)) + sum(f_pred,2);
Q_norm = Q./(T*sum(f_pred./T,2).^2);











%% Regression from MUA by depth to fluorescence: Thin spikes

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(1000,3500,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [400,1500,2000,2300,3000,4000];
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

samplerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    %         curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline(depth_group == curr_depth & ...
        ismember(spike_templates,find(msn)));
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

spike_thin_frac = 1-(min(sum(frame_spikes,2))./sum(frame_spikes,2));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    %     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline(depth_group == curr_depth & ...
        ismember(spike_templates,find(msn)));
    
    curr_thinned_spikes = curr_spike_times;
    curr_thinned_spikes( ...
        randperm(length(curr_thinned_spikes), ...
        round(length(curr_thinned_spikes)*spike_thin_frac(curr_depth)))) = [];
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    frame_spikes(curr_depth,:) = histcounts(curr_thinned_spikes,frame_edges);
    
end

use_svs = 1:50;
kernel_frames = -35*3:35*3;
downsample_factor = 1;
lambda = 0;
zs = [true,false];
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel( ...
    downsample(frame_spikes(:,use_frames)',downsample_factor)', ...
    downsample(fV(use_svs,use_frames)',downsample_factor)', ...
    kernel_frames_downsample,lambda,zs,cvfold);

k = permute(reshape(k,size(frame_spikes,1),length(kernel_frames_downsample),[]),[3,2,1]);

fluor_kernel = arrayfun(@(x) svdFrameReconstruct(U(:,:,use_svs),k(:,:,x)),1:size(k,3),'uni',false);
AP_image_scroll(cat(4,fluor_kernel{:}),kernel_frames_downsample*downsample_factor/samplerate);
truesize;

% Get map of explained variance
downsample_factor = 10;
spatial_explained_var = AP_spatial_explained_var(U(:,:,use_svs), ...
    fV(use_svs,use_frames),predicted_fluor,downsample_factor);
figure;imagesc(spatial_explained_var);
caxis([-max(abs(spatial_explained_var(:))),max(abs(spatial_explained_var(:)))]);
colormap(colormap_BlueWhiteRed);
colorbar;



%% Testing normalizing rate by template

%% STA for multiunit

use_spikes_idx = ismember(spike_templates,find(templateDepths > 0 & templateDepths < 1500));

use_spikes = spike_times_timeline(use_spikes_idx);
use_templates = spike_templates(use_spikes_idx);

use_templates_unique = unique(use_templates);

frame_spikes_templates = zeros(length(use_templates_unique),length(frame_t));
for curr_template_idx = 1:length(use_templates_unique)
    
    curr_template = use_templates_unique(curr_template_idx);
    
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);
    curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
        curr_spike_times + surround_times(2) > frame_t(end)) = [];
    
    if isempty(curr_spike_times)
        continue
    end
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/samplerate];
    frame_spikes(curr_template_idx,:) = histcounts(curr_spike_times,frame_edges);
    
end

frame_spikes_zs = mean(zscore(frame_spikes,[],2),1);

surround_times = [-3,3];
samplerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*samplerate):round(surround_times(2)*samplerate);
sta_t = surround_frames./samplerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames));

% Prepare weighting matrix for each frame
frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames);
    frames_w(:,curr_sta_frame) = ...
        circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_v = fV*frames_w;
sta_im = svdFrameReconstruct(U,sta_v);

% Normalize the image
%sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,px_mean),px_std);

% Draw the movie
AP_image_scroll(sta_im,sta_t);








%% Map of explained variance spikes -> fluorescence

sample_rate = (1/median(diff(frame_t)))*1;

% Skip the first/last n seconds to do this
skip_seconds = 60;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = linspace(700,3500,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [3500 Inf];

[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

binned_spikes = zeros(length(depth_group_edges_use)-1,length(time_bins)-1);
for curr_depth = 1:length(depth_group_edges_use)-1
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    %     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
    %         ismember(spike_templates,find(msn)));
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
    
end

use_svs = 1:50;
kernel_frames = -50:50;
downsample_factor = 1;
lambda = 0;
zs = [true,false];
cvfold = 5;

fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';

% TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf_resample, ...
%     binned_spikes,kernel_frames_downsample,lambda,zs,cvfold);
% TO USE dfV
[k,predicted_fluor,explained_var] = ...
    AP_regresskernel(binned_spikes(:,2:end-1), ...
    conv2(diff(fVdf_resample,[],2),[1,1]/2,'valid'), ...
    kernel_frames,lambda,zs,cvfold);

k = permute(reshape(k,size(binned_spikes,1),length(kernel_frames),[]),[3,2,1]);

fluor_kernel = arrayfun(@(x) svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,x)),1:size(k,3),'uni',false);
AP_image_scroll(cat(4,fluor_kernel{:}),kernel_frames/framerate);
truesize;

% Get map of explained variance
downsample_factor = 10;
spatial_explained_var = AP_spatial_explained_var(Udf(:,:,use_svs), ...
    conv2(diff(fVdf_resample,[],2),[1,1]/2,'valid'),predicted_fluor,downsample_factor);
figure;imagesc(spatial_explained_var);
caxis([-max(abs(spatial_explained_var(:))),max(abs(spatial_explained_var(:)))]);
colormap(colormap_BlueWhiteRed);
colorbar;







%%

raster_window = [-2,2];
use_spikes = spike_templates == 420;
psthViewer(spike_times(use_spikes),spike_templates(use_spikes), ...
    reward_t_timeline,raster_window,ones(size(reward_t_timeline)));


%% Plot spike properties

plot_templates = templateDepths < 1500;

% Plot the waveforms and spike statistics
figure;

subplot(1,3,1); hold on;
plot(waveform_t,nanmean(waveforms(uin,:),1),'c')
plot(waveform_t,nanmean(waveforms(tan,:),1),'g')
plot(waveform_t,nanmean(waveforms(fsi,:),1),'b')
plot(waveform_t,nanmean(waveforms(msn,:),1),'m')
xlabel('Time (ms)')
legend({'MSN','FSI','TAN','UIN'});

subplot(1,3,2); hold on;
stem3( ...
    templateDuration_us(msn & plot_templates)/1000, ...
    prop_long_isi(msn & plot_templates), ...
    spike_rate(msn & plot_templates),'m');

stem3( ...
    templateDuration_us(fsi & plot_templates)/1000, ...
    prop_long_isi(fsi & plot_templates), ...
    spike_rate(fsi & plot_templates),'b');

stem3( ...
    templateDuration_us(tan & plot_templates)/1000, ...
    prop_long_isi(tan & plot_templates), ...
    spike_rate(tan & plot_templates),'g');

stem3( ...
    templateDuration_us(uin & plot_templates)/1000, ...
    prop_long_isi(uin & plot_templates), ...
    spike_rate(uin & plot_templates),'c');

xlabel('waveform duration (ms)')
ylabel('frac long ISI')
zlabel('spike rate')

set(gca,'YDir','reverse')
set(gca,'XDir','reverse')
view(3);
grid on;
axis vis3d;

isi_edges = [0:0.001:0.3];
isi_centers = conv(isi_edges,[1,1]/2,'valid');

isi_hist = nan(max(spike_templates),length(isi_centers));
for curr_template = 1:max(spike_templates)
    curr_spikes = spike_times_timeline(spike_templates == curr_template);
    curr_isi = diff(curr_spikes);
    isi_hist(curr_template,:) = histcounts(curr_isi,isi_edges);
end

isi_hist_norm = bsxfun(@rdivide,bsxfun(@minus,isi_hist,min(isi_hist,[],2)),max(isi_hist,[],2) - min(isi_hist,[],2));

msn_isi = nanmedian(isi_hist_norm(plot_templates & msn,:),1);
fsi_isi = nanmedian(isi_hist_norm(plot_templates & fsi,:),1);
tan_isi = nanmedian(isi_hist_norm(plot_templates & tan,:),1);
uin_isi = nanmedian(isi_hist_norm(plot_templates & uin,:),1);

subplot(1,3,3); hold on;
plot(isi_centers,msn_isi,'m','linewidth',2);
plot(isi_centers,fsi_isi,'b','linewidth',2);
plot(isi_centers,tan_isi,'g','linewidth',2);
plot(isi_centers,uin_isi,'c','linewidth',2);

xlabel('Inter-spike interval');
ylabel('Normalized frequency');

%% Transfer all basket data to zserver
% I don't have the guts to movefile, so I'm copyfile and then will check

basket_path = '\\basket.cortexlab.net\data\ajpeters';
zserver_path = '\\zserver.cortexlab.net\Data\Subjects';

% Get all animal folders in basket
basket_dir = dir([basket_path filesep 'AP*']);
animals = {basket_dir.name};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    animal_dir = dir([basket_path filesep animal]);
    days = {animal_dir(3:end).name};
    
    for curr_day = 1:length(days)
        day = days{curr_day};
        
        curr_basket_path = [basket_path filesep animal filesep day filesep 'ephys' filesep '*'];
        curr_zserver_path = [zserver_path filesep animal filesep day filesep 'ephys' filesep 'kilosort'];
        
        curr_basket_files = dir(curr_basket_path);
        if length(curr_basket_files) > 2
            disp([curr_basket_path ' ---> ' curr_zserver_path])
            if ~exist(curr_zserver_path,'dir')
                mkdir(curr_zserver_path)
            end
            copyfile(curr_basket_path,curr_zserver_path);
        else
            disp(['Skipping: ' curr_basket_path])
        end
    end
    
end


%% Check that the contents and sizes of all transferred files are ok

basket_path = '\\basket.cortexlab.net\data\ajpeters';
zserver_path = '\\zserver.cortexlab.net\Data\Subjects';

% Get all animal folders in basket
basket_dir = dir([basket_path filesep 'AP*']);
animals = {basket_dir.name};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    animal_dir = dir([basket_path filesep animal]);
    days = {animal_dir(3:end).name};
    
    for curr_day = 1:length(days)
        day = days{curr_day};
        
        curr_basket_path = [basket_path filesep animal filesep day filesep 'ephys' filesep '*'];
        curr_zserver_path = [zserver_path filesep animal filesep day filesep 'ephys' filesep 'kilosort'];
        
        curr_basket_files = dir(curr_basket_path);
        curr_zserver_files = dir(curr_zserver_path);
        
        if length(curr_basket_files) > 2
            same_files = isequaln(curr_basket_files(~[curr_basket_files.isdir]),curr_zserver_files(~[curr_zserver_files.isdir]));
            if ~same_files
                error(['Bad: ' curr_basket_path])
            elseif same_files
                disp(['Good: ' curr_basket_path])
            end
        else
            disp(['Skipping: ' curr_basket_path])
        end
    end
end






%% Testing move onset-aligned MUA
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        use_trials = true(1,min(length(signals_events.trialSideValues),length(stimOn_times)));
        stimIDs = signals_events.trialSideValues(use_trials).*signals_events.trialContrastValues(use_trials);
        stim_onsets = stimOn_times(use_trials);
        
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-5.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_move_hit = nan(6,length(t_bins),length(conditions));
        mua_move_miss = nan(6,length(t_bins),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);
                
                use_trials = find(stimIDs == curr_condition & signals_events.hitValues(use_trials) == 1 & ~isnan(wheel_move_time));
                use_move_onsets = wheel_move_time(use_trials);
                if length(use_move_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                use_trials = find(stimIDs == curr_condition & signals_events.missValues(use_trials) == 1 & ~isnan(wheel_move_time));
                use_move_onsets = wheel_move_time(use_trials);
                if length(use_move_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_move_hit(:,:,:,curr_day) = mua_move_hit;
        batch_vars(curr_animal).mua_move_miss(:,:,:,curr_day) = mua_move_miss;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_move_choiceworld'],'batch_vars');

%% Testing feedback onset-aligned MUA
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        use_trials = true(1,min(length(signals_events.trialSideValues),length(stimOn_times)));
        stimIDs = signals_events.trialSideValues(use_trials).*signals_events.trialContrastValues(use_trials);
        stim_onsets = stimOn_times(use_trials);
        
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-5.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_feedback_hit = nan(6,length(t_bins),length(conditions));
        mua_feedback_miss = nan(6,length(t_bins),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);
                
                curr_trials = find(stimIDs == curr_condition & signals_events.hitValues(use_trials) == 1);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_feedback_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                curr_trials = find(stimIDs == curr_condition & signals_events.missValues(use_trials) == 1);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_feedback_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_feedback_hit(:,:,:,curr_day) = mua_feedback_hit;
        batch_vars(curr_animal).mua_feedback_miss(:,:,:,curr_day) = mua_feedback_miss;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_feedback_choiceworld'],'batch_vars');




%% trying shuffle data correlation with GPU instead of for loop
% blargh this wasn't any faster

tic

corr_mua_mua = cell(size(event_aligned_spikes,3),size(event_aligned_spikes,3));
for curr_mua1 = 1:size(event_aligned_spikes,3)
    for curr_mua2 = 1:size(event_aligned_spikes,3)
        
        curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua1),[],1);
        curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua2),[],1);
        curr_data2_shuff = curr_data2(shuff_idx);
        
        corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
        
        corr_shuff = ...
            gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
            gpuArray(curr_data2_shuff)))./(length(align_times)-1);
        
        corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
        corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
        
        corr_mua_mua{curr_mua1,curr_mua2} = corr_real;
        corr_mua_mua{curr_mua1,curr_mua2}(~corr_sig) = 0;
        
    end
end

toc

tic

curr_data1 = zscore(event_aligned_spikes,[],1);
curr_data2 = zscore(event_aligned_spikes,[],1);

curr_data2_shuff = zscore(shake(repmat(event_aligned_spikes,[1,1,1,n_shuff]),1),[],1);

corr_real = gather(pagefun(@mtimes,gpuArray(permute(curr_data1,[2,1,4,3])), ...
    gpuArray(curr_data2)));

corr_sig = nan(size(corr_real));
for curr_page = 1:size(curr_data1,3); % too big to do without FOR
    corr_shuff = ...
        gather(pagefun(@mtimes,gpuArray(permute(curr_data1(:,:,curr_page),[2,1])), ...
        gpuArray(curr_data2_shuff)));
    corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],4);
    corr_sig(:,:,:,curr_page) = corr_real(:,:,:,curr_page) < corr_shuff_cutoff(:,:,:,1) | ...
        corr_real(:,:,:,curr_page) > corr_shuff_cutoff(:,:,:,2);
end

corr_real(~corr_sig) = 0;

toc

% this is the memory step that breaks it
corr_shuff = ...
    gather(pagefun(@mtimes,gpuArray(permute(curr_data1,[2,1,4,5,3])), ...
    gpuArray(curr_data2_shuff)));


%% (fixing file name mistake)

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_dir = dir(data_path);
a = cellfun(@(x) any(strfind(x,'_rly')),{data_dir.name});
for fn = {data_dir(a).name}
    use_fn = fn{:};
    idx = strfind(use_fn,'_rly');
    use_fn(idx+3) = [];
    movefile([data_path filesep fn{:}],[data_path filesep use_fn]);
end


%% logistic regression using Vs instead of ROIs

aUdf = single(AP_align_widefield(animal,day,Udf));

% Get aligned V's
% Set times to pull
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

event_aligned_V = nan(length(stimOn_times),length(t),size(U,3),2);
for curr_align = 1:2
    
    % Align by stim or movement
    switch curr_align
        case 1
            use_align = stimOn_times;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    
    % Get event-aligned Vs
    t_peri_event = bsxfun(@plus,use_align,t);
    event_aligned_V(:,:,:,curr_align) = ...
        interp1(conv2(frame_t,[1,1]/2,'valid'),diff(fVdf,[],2)',t_peri_event);
    
end

% Pick trials to use
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5 & ...
    stim_to_move < 0.5;

cvfold = 10;

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = signals_events.trialSideValues == -1;
R_trials = signals_events.trialSideValues == 1;

D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);

D.response = ((trial_choice(use_trials)'+1)/2)+1;
D.repeatNum = ones(sum(use_trials),1);

% Fit stim model
g_stim_all = GLM(D).setModel('AP_test_stim').fit;
behavParameterFit = g_stim_all.parameterFits;

% Fit cross-validated stim model
g_stim = GLM(D).setModel('AP_test_stim').fitCV(cvfold);
pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
likelihood = pL.*(g_stim.data.response==1) + pR.*(g_stim.data.response==2);
loglik_bpt_stim = mean(log2(likelihood));

% Get "empirical model"
trial_partition_ordered = round(linspace(1,cvfold,sum(use_trials)))';
trial_partition = trial_partition_ordered(randperm(length(trial_partition_ordered)));

contrasts = [0,0.06,0.125,0.25,0.5,1];
contrast_sides = unique([-1*contrasts,contrasts]);
contrast_side_trial = diff(D.stimulus,[],2);
[~,contrast_sides_id] = ismember(contrast_side_trial,contrast_sides','rows');

pL_empirical = nan(sum(use_trials),1);
for curr_cv = 1:cvfold
    
    train_trials = trial_partition ~= curr_cv;
    test_trials = trial_partition == curr_cv;
    
    go_left_n_train = accumarray(contrast_sides_id(train_trials), ...
        D.response(train_trials) == 1,[length(contrast_sides),1]);
    contrast_sides_n_train = accumarray(contrast_sides_id(train_trials),1, ...
        [length(contrast_sides),1]);
    go_left_frac_train = go_left_n_train./contrast_sides_n_train;
    
    pL_empirical(test_trials) = go_left_frac_train(contrast_sides_id(test_trials));
    
end

likelihood = pL_empirical.*(D.response==1) + (1-pL_empirical).*(D.response==2);
loglik_bpt_empirical = mean(log2(likelihood(likelihood ~= 0)));

% Fit fluor model (time averaged)
use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.2 & t < -0.02;

use_align = 1;
use_t = use_t_stim;

use_Vs = 1:100;

D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
% D.offset_ZL = log(pL_empirical./(1-pL_empirical));
D.neur = double(squeeze(nanmean(event_aligned_V(use_trials,use_t,use_Vs,use_align),2)));
g_neur = GLM(D).setModel('AP_test_V_stim').fitCV(cvfold);
pL = g_neur.p_hat(:,1);
pR = g_neur.p_hat(:,2);
likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
loglik_bpt_fluor = mean(log2(likelihood(~isnan(likelihood))));

a = nanmean(g_neur.parameterFits(:,2:end),1);
b = svdFrameReconstruct(aUdf(:,:,use_Vs),a');
figure;imagesc(b);
axis image off
caxis([-1e-11,1e-11])
colormap(colormap_BlueWhiteRed);

% Fit fluor model (across all timepoints)
use_Vs = 1:100;

pV = nan(size(aUdf,1),size(aUdf,2),length(t),2);
loglik_bpt_fluor = nan(size(t));
for curr_align = 1
    for curr_t = 1:length(t)
        
        D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
        D.neur = double(squeeze(event_aligned_V(use_trials,curr_t,use_Vs,curr_align)));
        g_neur = GLM(D).setModel('AP_test_V_stim').fitCV(cvfold);
        
        pL = g_neur.p_hat(:,1);
        pR = g_neur.p_hat(:,2);
        likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
        loglik_bpt_fluor(curr_t) = mean(log2(likelihood));
        
        pV_mean_cv = nanmean(g_neur.parameterFits(:,2:end),1);
        pV(:,:,curr_t,curr_align) = svdFrameReconstruct(aUdf(:,:,use_Vs),pV_mean_cv');
        
        AP_print_progress_fraction(curr_t,length(t));
    end
end
AP_image_scroll(pV,t);
figure;
plot(t,loglik_bpt_fluor - loglik_bpt_stim,'k','linewidth',2)
xlabel('Time (s)');
ylabel('Log likelihood increase');

%% Plot rate over time (for data club)

n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spikeDepths,depth_group_edges);

spiking_stat_window = 60*0.5; % seconds
spiking_stat_bins = min(spike_times_timeline):spiking_stat_window: ...
    max(spike_times_timeline);
spiking_stat_centers = conv2(spiking_stat_bins,[1,1]/2,'valid');

% Get firing rate across the session
bin_spikes = nan(n_depths, ...
    length(spiking_stat_bins)-1);
for curr_depth = 1:n_depths
    bin_spikes(curr_depth,:) = ...
        histcounts(spike_times_timeline(depth_group == curr_depth), ...
        spiking_stat_bins);
end

figure; hold on
set(gca,'ColorOrder',copper(n_depths));
plot(spiking_stat_centers/60,bsxfun(@rdivide,bin_spikes,bin_spikes(:,1))','linewidth',2);
% plot(spiking_stat_centers/60,bin_spikes','linewidth',2);
line(repmat(acqLive_timeline(1),1,2)/60,ylim,'color','k');
line(repmat(acqLive_timeline(2),1,2)/60,ylim,'color','k');
xlabel('Time from experiment start (min)');
ylabel('Spikes normalized to bin 1')




n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spikeDepths,depth_group_edges);

spiking_stat_window = 0.005; % seconds
spiking_stat_bins = acqLive_timeline(1):spiking_stat_window: ...
    acqLive_timeline(2);
spiking_stat_centers = conv2(spiking_stat_bins,[1,1]/2,'valid');

% Get firing rate across the session
bin_spikes = nan(n_depths, ...
    length(spiking_stat_bins)-1);
for curr_depth = 1:n_depths
    bin_spikes(curr_depth,:) = ...
        histcounts(spike_times_timeline(depth_group == curr_depth), ...
        spiking_stat_bins);
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_spikes = conv2(bin_spikes,smWin,'same');


%% Scatter move L/R trials

early_condition = find(ismember(conditions(:,4),[1],'rows'));
curr_trials = ismember(trial_id,early_condition) & use_trials';

use_align = stimOn_times;
use_align(isnan(use_align)) = 0;
t_peri_event = bsxfun(@plus,use_align,t);
event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);

use_t = t > 0.05 & t < 0.15;
ddf_mean = squeeze(nanmean(event_aligned_ddf(:,use_t,:),2));
ddf_group = reshape(bsxfun(@plus,(trial_choice(curr_trials)'+1)/4,1:n_rois),[],1);

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

plot_conditions = conditions(:,4) == 1;
plot_contrasts = [unique(conditions(plot_conditions,1))*sides(1); ...
    unique(conditions(plot_conditions,1))*sides(2)];
[~,sort_idx] = sort(plot_contrasts);

plot_roi = 10;
figure; hold on;
L_trials = trial_conditions(:,3) == -1;
R_trials = trial_conditions(:,3) == 1;

scatter(trial_conditions(curr_trials & L_trials,1).*trial_conditions(curr_trials & L_trials,2), ...
    ddf_mean(curr_trials & L_trials,plot_roi),100,'k','filled','MarkerEdgeColor','r','linewidth',3);
scatter(trial_conditions(curr_trials & R_trials,1).*trial_conditions(curr_trials & R_trials,2), ...
    ddf_mean(curr_trials & R_trials,plot_roi),100,'k','filled','MarkerEdgeColor','b','linewidth',3);

%% Checking that all face/eyecam files are copied over

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocols = {'vanillaChoiceworld', ...
    'stimSparseNoiseUncorrAsync', ...
    'stimKalatsky', ...
    'AP_choiceWorldStimPassive'};

for curr_protocol = 1:length(protocols)
    
    protocol = protocols{curr_protocol};
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        experiments = AP_find_experiments(animal,protocol);
        
        experiments = experiments([experiments.imaging]);
        
        for curr_day = 1:length(experiments);
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment;
            [~,eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
            [~,facecam_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam');
            
            if ~eyecam_exists
                disp([animal filesep day filesep protocol filesep 'eyecam']);
            end
            if ~facecam_exists
                disp([animal filesep day filesep protocol filesep 'facecam']);
            end
            
        end
    end
end


%% testing peter model mean time (on concatenated data)

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < -0.02;
use_t_align = [use_t_stim;use_t_move];

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = 6;

cvfold = 20;

% Concatenate fluorescence
fluor_cat = cat(1,fluor_all{:});

% HERE: I started to estimate "gain" parameters for all days to equalize,
% I'm just going to include that in the model
% [avg_fluor,gname] = grpstats(fluor_cat(:,:,1,1),[diff(D.stimulus,[],2),trial_day],{'mean','gname'});
% grp_contrast = cellfun(@str2num,gname(:,1));
% grp_day = cellfun(@str2num,gname(:,2));
%
% contrasts = [0,0.06,0.125,0.25,0.5,1];
% sides = [-1,1];
% contrast_sides = sort(unique(bsxfun(@times,contrasts',sides)));
%
% day_gain = ones(length(fluor_all),1);
%
% model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
%
% start_params = [0,0,1,0];
%
% options = struct;
% options.MaxFunEvals = 1000;
% options.Display = 'off';
%
% params_fit = fminsearch(model_sse,start_params,options);
% activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:) = ...
%     params_fit;

% Concatenate MUA and normalize
mua_cat_raw = cat(1,mua_all{:});
t_baseline = t < 0;
softnorm = 5;
mua_baseline = nanmean(mua_cat_raw(:,t_baseline,:,1),2);
mua_cat = bsxfun(@rdivide,mua_cat_raw,mua_baseline + softnorm);

% Get L-R activity
fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
    fluor_cat(:,:,n_rois/2+1:end,:));

% Concatenate behavioural data
D = struct;
D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all,'uni',false));
D.response = cell2mat(cellfun(@(x) x.response,D_all,'uni',false));
D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all,'uni',false));

trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
    num2cell(1:length(fluor_all))',fluor_all,'uni',false));
D.day = trial_day;

% % (replicate original contrast-shuffling typo)
% D.stimulus(D.stimulus(:,1) ~= 0,1) = shake(D.stimulus(D.stimulus(:,1) ~= 0,1));
% D.stimulus(D.stimulus(:,2) ~= 0,2) = shake(D.stimulus(D.stimulus(:,2) ~= 0,2));

% Get zero-contrasts as a subset
zero_contrasts = D.stimulus(:,1) == 0 & D.stimulus(:,2) == 0;
D_zero = struct;
D_zero.stimulus = D.stimulus(zero_contrasts,:);
D_zero.response = D.response(zero_contrasts);
D_zero.repeatNum = D.repeatNum(zero_contrasts);

% Fit stim all
use_model = 'AP_test_stim';
g_stim_all = GLM(D).setModel(use_model).fit;
behavParameterFit = g_stim_all.parameterFits;

D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));

% Fit stim cross-validated
use_model = 'AP_test_stim';

g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
likelihood = pL.*(g_stim.data.response == 1) + pR.*(g_stim.data.response == 2);

loglik_bpt_stim = nanmean(log2(likelihood));

% Fit stim + activity (all time)
use_model = 'AP_test_neur_stimoffset_dayscale';

loglik_bpt_fluor = nan(length(t),n_rois,2);
for curr_align = 1:2
    for curr_area = 1:n_rois
        for curr_t = 1:length(t)
            
            D.neur = reshape(fluor_cat_hemidiff(:,curr_t,curr_area,curr_align),[],1);
            
            clear g_act
            g_act = GLM(D).setModel(use_model).fitCV(cvfold);
            pL = g_act.p_hat(:,1);
            pR = g_act.p_hat(:,2);
            likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
            
            loglik_bpt_fluor(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
            AP_print_progress_fraction(curr_t,length(t));
        end
    end
end
loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;

% Fit stim + activity (mean time)
use_model = 'AP_test_neur_stimoffset_dayscale';

loglik_bpt_fluor = nan(n_rois,2);
for curr_align = 1:2
    for curr_area = 1:n_rois
        
        D.neur = reshape(nanmean(fluor_cat_hemidiff(:,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
        
        clear g_act
        g_act = GLM(D).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
        
        loglik_bpt_fluor(curr_area,curr_align) = nanmean(log2(likelihood));
        
        AP_print_progress_fraction(curr_area,n_rois);
    end
end
loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;

loglik_bpt_mua = nan(n_depths,2);
for curr_align = 1:2
    for curr_area = 1:n_depths
        
        D.neur = reshape(nanmean(mua_cat(:,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
        
        clear g_act
        g_act = GLM(D).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
        
        loglik_bpt_mua(curr_area,curr_align) = nanmean(log2(likelihood));
        
        AP_print_progress_fraction(curr_area,n_depths);
    end
end
loglik_increase_mua = loglik_bpt_mua - loglik_bpt_stim;

% Fit stim + activity (ZERO CONTRASTS)
use_model = 'AP_test_neur_nostim';

loglik_bpt_fluor_zerocontrast = nan(n_rois,2);
loglik_bpt_guess = nan(n_rois,2);
for curr_align = 1:2
    for curr_area = 1:n_rois
        
        D_zero.neur = reshape(nanmean(fluor_cat_hemidiff(zero_contrasts,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
        
        clear g_act
        g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
        
        loglik_bpt_fluor_zerocontrast(curr_area,curr_align) = nanmean(log2(likelihood));
        loglik_bpt_guess(curr_area,curr_align) = g_act.guess_bpt;
        
        AP_print_progress_fraction(curr_area,n_rois);
    end
end
loglik_increase_fluor_zerocontrast = bsxfun(@minus,loglik_bpt_fluor_zerocontrast,nanmean(loglik_bpt_guess,1));

loglik_bpt_mua_zerocontrast = nan(n_depths,2);
loglik_bpt_guess = nan(n_depths,2);
for curr_align = 1:2
    for curr_area = 1:n_depths
        
        D_zero.neur = reshape(nanmean(mua_cat(zero_contrasts,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
        
        clear g_act
        g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
        
        loglik_bpt_mua_zerocontrast(curr_area,curr_align) = nanmean(log2(likelihood));
        loglik_bpt_guess(curr_area,curr_align) = g_act.guess_bpt;
        
        AP_print_progress_fraction(curr_area,n_depths);
    end
end
loglik_increase_mua_zerocontrast = bsxfun(@minus,loglik_bpt_mua_zerocontrast,nanmean(loglik_bpt_guess,1));

% Fit stim + activity (all L/R ROIs together)
use_model = 'AP_test_roi_stimoffset';
roi_params = nan(n_rois,2);
for curr_align = 1:2
    
    D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
    D.neur = squeeze(nanmean(fluor_cat(:,use_t_align(curr_align,:),:,curr_align),2));
    
    clear g_act
    g_act = GLM(D).setModel(use_model).fit;
    roi_params(:,curr_align) = g_act.parameterFits(2:end);
    
    AP_print_progress_fraction(curr_align,2);
end




%% testing peter model mean time (standalone)

% Prepare fluorescence
% (load widefield ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

aUdf = single(AP_align_widefield(animal,day,Udf));
roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));

% (get ddf)
roi_traces_derivative = diff(roi_trace,[],2);
roi_traces_derivative(roi_traces_derivative < 0) = 0;

% (make L-R traces)
roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
    roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);

% Get event-aligned activity
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
for curr_align = 1:2
    switch curr_align
        case 1
            use_align = stimOn_times;
            use_align(isnan(use_align)) = 0;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    
    t_peri_event = bsxfun(@plus,use_align,t);
    
    % Fluorescence
    event_aligned_ddf(:,:,:,curr_align) = ...
        interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
    
end

% Pick trials to use
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5 & ...
    stim_to_move < 0.5;

cvfold = 10;

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = signals_events.trialSideValues(1:n_trials) == -1;
R_trials = signals_events.trialSideValues(1:n_trials) == 1;

D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);

% D.response = (trial_choice(use_trials)'+1)/2+1;
D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
D.repeatNum = ones(sum(use_trials),1);

% Fit stim crossvalidated
use_model = 'AP_test_stim';

g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
likelihood = pL.*(g_stim.data.response==1) + pR.*(g_stim.data.response==2);

stim_params = nanmean(g_stim.parameterFits,1);
loglik_bpt_stim = nanmean(log2(likelihood(likelihood ~= 0)));

% Fit stim all
use_model = 'AP_test_stim';
g_stim_all = GLM(D).setModel(use_model).fit;
behavParameterFit = g_stim_all.parameterFits;

% (temp: double check stim is correct)
[psych,grp] = grpstats(trial_choice(use_trials) == 1,trial_conditions(use_trials,1).*trial_conditions(use_trials,2),{'mean','gname'});
grp = cellfun(@str2num,grp);

% Fit stim w/ neural
use_model = 'AP_test_neur_stimoffset';

use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < 0.02;

curr_align = 1;
use_t = use_t_stim;

D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
loglik_bpt = nan(n_rois,1);
for curr_area = 1:n_rois
    
    D.neur = reshape(nanmean(event_aligned_ddf(use_trials,use_t,curr_area,curr_align),2),[],1);
    
    clear g_fluor
    g_act = GLM(D).setModel(use_model).fitCV(cvfold);
    pL = g_act.p_hat(:,1);
    pR = g_act.p_hat(:,2);
    likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
    
    loglik_bpt_fluor(curr_area) = nanmean(log2(likelihood(likelihood ~= 0)));
    
end
loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;
figure;plot(loglik_increase_fluor,'k')




%% attempt to binarize ctx>str kernels

n_experiments = arrayfun(@(x) length(batch_vars(x).r_px),1:length(animals));
ctx_str_k_binary = false(437,416,6,sum(n_experiments));

curr_expt = 1;
for curr_animal = 1:length(animals)    
    for curr_day = 1:length(batch_vars(curr_animal).r_px);
        
        r_px = batch_vars(curr_animal).r_px{curr_day};
        
        r_px_max = squeeze(max(r_px,[],3)).^3;
        
        r_px_max(isnan(r_px_max)) = 0;
        for i = 1:n_depths
            r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[20,20]);
        end
        
        r_px_binary = false(size(r_px_max));
        for curr_depth = 1:size(r_px_max,3)
            r_px_binary(:,:,curr_depth) = imbinarize(r_px_max(:,:,curr_depth), ...
                prctile(reshape(r_px_max(:,:,curr_depth),[],1),90));
        end
                
        r_px_com = sum(bsxfun(@times,r_px_binary,permute(1:n_depths,[1,3,2])),3)./sum(r_px_binary,3);
        com_colored = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
        
        figure;image(com_colored);
        AP_reference_outline('ccf_aligned','w');AP_reference_outline('retinotopy','m');
        axis image off;
        drawnow;
        
        ctx_str_k_binary(:,:,:,curr_expt) = r_px_binary;
        curr_expt = curr_expt + 1;
    end
end

disp('done')


%% try to fix hemo correction
% somehow I just keep making it worse??

animal = 'AP029'; day = '2017-12-14'; experiment = 3; verbose = true; load_parts.imaging = true; AP_load_experiment;

% Trying in V space? but they didn't do this at first?
% (actually this is hemocorrectnonlocal)
a = Vh_Un(:,use_hemo_frames)'\Vn_th(:,use_hemo_frames)';
b = transpose(Vh_Un'*a);
c = Vn_th - b;

highpassCutoff = 0.01; % Hz
[b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
fVn_hemo = single(filtfilt(b100s,a100s,double(c)')');

[Udf,fVdf] = dffFromSVD(Un,c,avg_im);


V2 = HemoCorrectNonlocal(Vn_th',Vh_Un',framerate,[0.01,13])';
fV2 = single(filtfilt(b100s,a100s,double(V2)')');
[Udf,fVdf] = dffFromSVD(Un,fV2,avg_im);


% try instead: df/f first, then hemo correct
% no this looks like crap too
[Undf,Vndf] = dffFromSVD(Un,Vn,avg_im_n);

Vn_th = SubSampleShift(Vndf,1,2);
Vh_Undf = ChangeU(Uh,Vh,Undf);
hemo_freq = [7,13];
[~,hemo_tform] = HemoCorrectLocal(Un,Vn_th(:,use_hemo_frames),Vh_Un(:,use_hemo_frames),framerate,hemo_freq,3);
zVh_Undf = bsxfun(@minus, Vh_Undf, mean(Vh_Undf(:,use_hemo_frames),2));
Vn_hemo = transpose(Vn_th' - zVh_Undf'*hemo_tform');



if verbose; disp('Correcting hemodynamics...'); end

skip_seconds = 60*1;
use_hemo_frames = frame_t > skip_seconds;

min_frames = min(size(Vn,2),size(Vh,2));
Vn = Vn(:,1:min_frames);
Vh = Vh(:,1:min_frames);

Vn_th = SubSampleShift(Vn,1,2);
Vh_Un = ChangeU(Uh,Vh,Un);

hemo_tform_fn = [data_path filesep 'hemo_tform.mat'];
if exist(hemo_tform_fn,'file')
    % If the hemo tform matrix has been computed, load and fix
    if verbose; disp('Using old hemo tform...'); end;
    load(hemo_tform_fn)
    zVh_Un = bsxfun(@minus, Vh_Un, mean(Vh_Un(:,use_hemo_frames),2));
    Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
else
    % If no p hemo tform matrix, compute and save
    if verbose; disp('Computing hemo tform...'); end
    %     hemo_freq = [0.1,2];
    hemo_freq = [7,13];
    [Vn_hemo,hemo_tform] = AP_hemocorrectlocal(Un,Vn_th(:,use_hemo_frames),Vh_Un(:,use_hemo_frames),framerate,hemo_freq,3);
    [Vn_hemo,hemo_tform] = HemoCorrectLocal(Un,Vn_th(:,use_hemo_frames),Vh_Un(:,use_hemo_frames),framerate,hemo_freq,3);
    save(hemo_tform_fn,'hemo_tform');
    % Close the figures (hacky - but function isn't mine)
    close(gcf)
    close(gcf)
end

if verbose; disp('Filtering...'); end;
% Don't bother filtering heartbeat, just detrend and highpass
% fVn_hemo = detrendAndFilt(Vn_hemo, framerate);
highpassCutoff = 0.01; % Hz
[b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');

dVn_hemo = detrend(Vn_hemo', 'linear')';
% wasn't zero-lag filtered before? why not?
%fVn_hemo = filter(b100s,a100s,dVn_hemo,[],2);
fVn_hemo = single(filtfilt(b100s,a100s,double(dVn_hemo)')');

% Do this for the colors individually, in case they're used
dVn = detrend(Vn', 'linear')';
fVn = single(filtfilt(b100s,a100s,double(dVn)')');

dVh = detrend(Vh', 'linear')';
fVh = single(filtfilt(b100s,a100s,double(dVh)')');

% set final U/V to use
fV = fVn_hemo;
U = Un;
avg_im = avg_im_n;
frame_t = th; % shifted to use hemo color times


if verbose; disp('Done.'); end

% Make dF/F
[Udf,fVdf] = dffFromSVD(U,fV,avg_im);




% (trying medfilting U's first for crappy data)
% (doesn't make hemo correction any better)
Un_m = Un;
Uh_m = Uh;
gUn = gpuArray(Un);
gUh = gpuArray(Uh);
for curr_U = 1:size(Un,3)
    Un_m(:,:,curr_U) = gather(medfilt2(gUn(:,:,curr_U),[3,3]));
    Uh_m(:,:,curr_U) = gather(medfilt2(gUh(:,:,curr_U),[3,3]));
    AP_print_progress_fraction(curr_U,size(U,3));
end

Vn_m = ChangeU(Un,Vn,Un_m);
Vh_m = ChangeU(Uh,Vh,Uh_m);

Vn_m_th = SubSampleShift(Vn_m,1,2);
Vh_m_Un = ChangeU(Uh_m,Vh_m,Un_m);

skip_seconds = 60*1;
use_hemo_frames = frame_t > skip_seconds;

hemo_freq = [7,13];
[~,hemo_tform] = HemoCorrectLocal(Un_m,Vn_m_th(:,use_hemo_frames),Vh_m_Un(:,use_hemo_frames),framerate,hemo_freq,3);
zVh_m_Un = bsxfun(@minus, Vh_m_Un, mean(Vh_m_Un(:,use_hemo_frames),2));
Vn_m_hemo = transpose(Vn_m_th' - zVh_m_Un'*hemo_tform');

highpassCutoff = 0.01; % Hz
[b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
dVn_m_hemo = detrend(Vn_m_hemo', 'linear')';
fVn_m_hemo = single(filtfilt(b100s,a100s,double(dVn_m_hemo)')');

[Udf,fVdf] = dffFromSVD(Un_m,fVn_m_hemo,avg_im);




% try medfilt U and then changing V?
% why does this totally break it? that really doesn't make sense
Udf_medfilt = gpuArray(double(Udf));
for curr_U = 1:size(U,3);
    Udf_medfilt(:,:,curr_U) = medfilt2(Udf(:,:,curr_U),[5,5]);
    AP_print_progress_fraction(curr_U,size(U,3));
end
Udf_medfilt = gather(Udf_medfilt);

fVdf_medfilt = ChangeU(Udf,fVdf,Udf_medfilt);


% Medfilt before df/f? this looks like it helps
U_m = U;
gU = gpuArray(U);
for curr_U = 1:size(U,3);
    U_m(:,:,curr_U) = gather(medfilt2(gU(:,:,curr_U),[5,5]));
    AP_print_progress_fraction(curr_U,size(U,3));
end

fV_medfilt = ChangeU(U,fV,U_m);

[Udf,fVdf] = dffFromSVD(U_m,fV_medfilt,medfilt2(avg_im,[5,5]));



% Linear fitting for each pixel - way better
AP_hemo = AP_hemo_correct(Un,Vn_th,Vh_Un,framerate,3);

hemo_freq = [7,13];
KH_hemo = HemoCorrectLocal(Un,Vn_th(:,use_hemo_frames),Vh_Un(:,use_hemo_frames),framerate,hemo_freq,3);
   



highpassCutoff = 0.01; % Hz
[b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');

fVn = single(filtfilt(b100s,a100s,double(Vn)')');
fVh = single(filtfilt(b100s,a100s,double(Vh)')');

[Un_df,fVn_df] = dffFromSVD(Un,fVn,avg_im_n);
[Uh_df,fVh_df] = dffFromSVD(Uh,fVh,avg_im_h);

fVn_df_th = SubSampleShift(fVn_df,1,2);        
fVh_df_Un = ChangeU(Uh_df,fVh_df,Un_df);       

skip_seconds = 60*1;
use_hemo_frames = frame_t > skip_seconds;

hemo_freq = [7,13];
[Vn_hemo,hemo_tform] = HemoCorrectLocal(Un_df,fVn_df_th,fVh_df_Un,framerate,hemo_freq,3);




%%

% Plot activity difference of conditions
plot_act = mua_allcat(:,:,3,1);

rxn_times = linspace(0,0.5,6);
rxn_time_centers = rxn_times(1:end-1)+diff(rxn_times)/2;
n_rxn_times = length(rxn_times)-1;

activity_split_0 = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

activity_split_6 = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_side_allcat == 1 & ...
    trial_contrast_allcat == 0.06 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

n_boot = 1000;
activity_6_0_max = cell2mat(cellfun(@ (x0,x6) ...
    prctile(bootstrp(n_boot,@(x) max(nanmean(x,1)),x6) - ...
    bootstrp(n_boot,@(x) max(nanmean(x,1)),x0),[50,2.5,97.5]), ...
    activity_split_0,activity_split_6,'uni',false)');

figure;
AP_errorfill(rxn_time_centers,activity_6_0_max(:,1), ...
    bsxfun(@minus,activity_6_0_max(:,2:3),activity_6_0_max(:,1)),'k',0.5);
line(xlim,[0,0],'linewidth',2);
ylabel('Activity stim - no stim')
xlabel('Reaction time (center)')


%%

% Plot activity difference of conditions
plot_act = mua_allcat(:,:,3,1);

rxn_times = linspace(0,0.5,6);
rxn_time_centers = rxn_times(1:end-1)+diff(rxn_times)/2;
n_rxn_times = length(rxn_times)-1;

activity_split_nostim = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

activity_split_stim = cell(5,n_rxn_times);
for curr_contrast_idx = 1:length(contrasts)-1
    activity_split_stim(curr_contrast_idx,:) = arrayfun(@(x) ...
        plot_act( ...
        move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
        trial_side_allcat == 1 & ...
        trial_contrast_allcat == contrasts(curr_contrast_idx+1) & ...
        trial_choice_allcat == -1,:), ...
        1:length(rxn_times)-1,'uni',false);
end

n_boot = 1000;
activity_diff_max = cellfun(@ (nostim,stim) ...
    prctile(bootstrp(n_boot,@(x) max(nanmean(x,1)),stim) - ...
    bootstrp(n_boot,@(x) max(nanmean(x,1)),nostim),[50,2.5,97.5]), ...
    repmat(activity_split_nostim,5,1),activity_split_stim,'uni',false)';
activity_diff_max_plot = arrayfun(@(x) vertcat(activity_diff_max{:,x}),1:5,'uni',false);

figure;
for i = 1:5    
    subplot(5,1,i);
    AP_errorfill(rxn_time_centers,activity_diff_max_plot{i}(:,1)', ...
        bsxfun(@minus,activity_diff_max_plot{i}(:,2:3), ...
        activity_diff_max_plot{i}(:,1))','k',0.5);
    line(xlim,[0,0],'linewidth',2);
    ylabel('Activity stim - no stim')
    xlabel('Reaction time (center)')
end














