clear all

%% Get leverpress-triggered average image

[roi_filename, roi_pathname] = uigetfile('*.*','ROI file');
[acq_filename, acq_pathname] = uigetfile('*.xsg','Acq file');
[bhv_filename, bhv_pathname] = uigetfile('*.mat','Behavior file');
[img_filename, img_pathname] = uigetfile('*.tif','Image file');
n_channels = 1; % in the final
% Get average frame n before and n after presses
n_frames = 20;

% load ROI data
load([roi_pathname roi_filename],'-MAT')

% pull out time for each for each frame from header, messy
img_info = imfinfo([img_pathname img_filename]);
msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
lines = img_info(1).Height;
fpms = msPerLine*lines; % time for single frame capture
% Get trials in s since started
TrialNumList = TK_ReadBitCode([acq_pathname acq_filename]);
% put trial time in ms
TrialNumList(2,:) = TrialNumList(2,:)*1000;

bhv = load([bhv_pathname bhv_filename],'-MAT');
lever_L = [];
lever_R = [];
reward_L = [];
reward_R = [];

for trial = TrialNumList(1,:)
    % get time points for all events in this trial
    bhv_current = bhv.saved_history.RewardsSection_LastTrialEvents{trial,1};
    % fixes possible bug of putting in old timepoints at the end
    LastGoodTrial = find(diff(bhv_current(:,3)) < 0,1);
    if ~isempty(LastGoodTrial)
        bhv_current = bhv_current(1:LastGoodTrial,:);
    end
    % get time lag of behavior to acq (t_bhv + timeDiff = t_acq) in ms
    trialStart_bhv = bhv_current(find(bhv_current(:,1) == 101 & bhv_current(:,2) == 0,1,'last'),3)*1000;
    timeDiff = TrialNumList(2,find(TrialNumList(1,:) == trial)) - trialStart_bhv;
    % time for first leverpress (left = 5 or right = 3) in ms, NaN for none
    % only while waiting for lever press (state 43)
    lever_L = [lever_L; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff]; % if first: & bhv_current(:,1) == 43
    lever_L(isempty(lever_L)) = NaN;
    lever_R = [lever_R; bhv_current(find(bhv_current(:,2) == 3),3)*1000 + timeDiff]; % if first: & bhv_current(:,1) == 43
    lever_R(isempty(lever_R)) = NaN;
    % get rewarded leverpress times
    reward_L = [reward_L; bhv_current(find(bhv_current(:,2) == 5 & bhv_current(:,1) == 43),3)*1000 + timeDiff];
    reward_R = [reward_R; bhv_current(find(bhv_current(:,2) == 3 & bhv_current(:,1) == 43),3)*1000 + timeDiff];

end
    
% Get frame for each leverpress and rewarded leverpress
frames_lever_L = ceil(lever_L(~isnan(lever_L))/fpms)';
frames_lever_R = ceil(lever_R(~isnan(lever_R))/fpms)';
frames_reward_L = ceil(reward_L(~isnan(reward_L))/fpms)';
frames_reward_R = ceil(reward_R(~isnan(reward_R))/fpms)';

% round up for number of channels
frames_lever_L = (1-mod(frames_lever_L,n_channels))+frames_lever_L;
frames_lever_R = (1-mod(frames_lever_R,n_channels))+frames_lever_R;
frames_reward_L = (1-mod(frames_reward_L,n_channels))+frames_reward_L;
frames_reward_R = (1-mod(frames_reward_R,n_channels))+frames_reward_R;

numframes=size(roi_trace,2); % number of frames from roi_trace
numrois=size(roi_trace,1); % number or ROIs

%% Event-triggered avg LEFT

% Get rid of leverpress events that don't have n_frames surrounding
frames_lever_L(frames_lever_L < (n_frames*n_channels) | frames_lever_L > (numframes+1-n_frames*n_channels)) = [];
num_leverpress = length(frames_lever_L);

w = waitbar(0,'Calculating trace ETA for LEFT');
for curr_roi = 1:size(roi_trace,1)
    tempTrace_L=zeros(num_leverpress,n_frames*2+1);
    for leverpress = 1:num_leverpress;
        curr_frames = frames_lever_L(leverpress)-n_frames*n_channels:n_channels:frames_lever_L(leverpress)+n_frames*n_channels;
        tempTrace_L(leverpress,:) = roi_trace(curr_roi,curr_frames);
    end
    avgTrace_L(curr_roi,:) = mean(tempTrace_L,1);
    SEMTrace_L(curr_roi,:) = std(tempTrace_L,1)./sqrt(num_leverpress);
    allTrace_L{curr_roi} = tempTrace_L;
    waitbar(curr_roi/size(roi_trace,1),w);
end
close(w);

% for i = 1:size(roi_trace,1);
%     figure
%     plot(avgTrace_L(i,:),'k')
%     hold on
%     plot(avgTrace_L(i,:)+SEMTrace_L(i,:).*1.96,'k--')
%     plot(avgTrace_L(i,:)-SEMTrace_L(i,:).*1.96,'k--')
%     title(['Event Average Left, ROI = ' num2str(i) ', n = ' num2str(num_leverpress)])
% end

%% Event-triggered avg RIGHT

% Get rid of leverpress events that don't have n_frames surrounding
frames_lever_R(frames_lever_R < (n_frames*n_channels) | frames_lever_R > (numframes+1-n_frames*n_channels)) = [];
num_leverpress = length(frames_lever_R);

w = waitbar(0,'Calculating trace ETA for RIGHT');
for curr_roi = 1:size(roi_trace,1)
    tempTrace_R=zeros(num_leverpress,n_frames*2+1);
    for leverpress = 1:num_leverpress;
        curr_frames = frames_lever_R(leverpress)-n_frames*n_channels:n_channels:frames_lever_R(leverpress)+n_frames*n_channels;
        tempTrace_R(leverpress,:) = roi_trace(curr_roi,curr_frames);
    end
    avgTrace_R(curr_roi,:) = mean(tempTrace_R,1);
    SEMTrace_R(curr_roi,:) = std(tempTrace_R,1)./sqrt(num_leverpress);
    allTrace_R{curr_roi} = tempTrace_R;
    waitbar(curr_roi/size(roi_trace,1),w);
end
    close(w);

% i = 2;
% figure
% plot(avgTrace_R(i,:),'k')
% hold on
% plot(avgTrace_R(i,:)+SEMTrace_R(i,:).*1.96,'k--')
% plot(avgTrace_R(i,:)-SEMTrace_R(i,:).*1.96,'k--')
% title(['Event Average Left, n = ' num2str(num_leverpress)])