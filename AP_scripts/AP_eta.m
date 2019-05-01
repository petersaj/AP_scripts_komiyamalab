clear all

%% Get leverpress-triggered average image

[acq_filename, acq_pathname] = uigetfile('*.xsg','Acq file');
[bhv_filename, bhv_pathname] = uigetfile('*.mat','Behavior file');
[img_filename, img_pathname] = uigetfile('*.tif','Image file');
n_channels = 2;
% Get average frame n before and n after presses
n_frames = 20;

% pull out time for each for each frame from header, messy
img_info = imfinfo([img_pathname img_filename]);
msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
lines = img_info(1).Height;
fpms = msPerLine*lines;
% Get trials in s since started
TrialNumList = TK_ReadBitCode([acq_pathname acq_filename]);
% put trial time in ms
TrialNumList(2,:) = TrialNumList(2,:)*1000;

bhv = load([bhv_pathname bhv_filename],'-MAT');
lever_L = [];
lever_R = [];

for trial = TrialNumList(1,:)
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
    lever_L = [lever_L; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff]; %bhv_current(:,1) == 43 & if first
    %lever_L(isempty(lever_L)) = NaN;
    lever_R = [lever_R; bhv_current(find(bhv_current(:,2) == 3),3)*1000 + timeDiff]; %bhv_current(:,1) == 43 & if first
    %lever_R(isempty(lever_R)) = NaN;
    % store leverpress times in TrialNumList, row3 = L, row4 = R
    %TrialNumList(3,find(TrialNumList(1,:) == trial)) = [lever_L + timeDiff];
    %TrialNumList(4,find(TrialNumList(1,:) == trial)) = [lever_R + timeDiff];
end
    
% Get frame for each leverpress
frames_lever_L = ceil(lever_L(~isnan(lever_L))/fpms)';
frames_lever_R = ceil(lever_R(~isnan(lever_R))/fpms)';

% round up for number of channels
frames_lever_L = (1-mod(frames_lever_L,n_channels))+frames_lever_L;
frames_lever_R = (1-mod(frames_lever_R,n_channels))+frames_lever_R;


numframes=length(img_info);
M=img_info(1).Width;
N=img_info(1).Height;

%% Event-triggered avg LEFT
im_L=zeros(N,M,n_frames*2+1);

% Get rid of leverpress events that don't have n_frames surrounding
frames_lever_L(frames_lever_L < (n_frames*n_channels) | frames_lever_L > (numframes+1-n_frames*n_channels)) = [];
num_leverpress = length(frames_lever_L);

w = waitbar(0,'Calculating ETA for LEFT');
for leverpress = 1:num_leverpress;
    for frame_current = frames_lever_L(leverpress)-n_frames*n_channels:n_channels:frames_lever_L(leverpress)+n_frames*n_channels;
        ETA_frame = ceil((frame_current-(frames_lever_L(leverpress)-n_frames*n_channels)+1)/n_channels);
        im_L(:,:,ETA_frame)=im_L(:,:,ETA_frame)+double(imread([img_pathname img_filename],'tiff',frame_current))./num_leverpress;
    end    
    waitbar(leverpress/num_leverpress,w,'Calculating ETA for LEFT');
end
close(w);

% Commented out all right lever-related analyses

% %% Event-triggered avg RIGHT
% im_R=zeros(N,M,n_frames*2+1);
% 
% % Get rid of leverpress events that don't have n_frames surrounding
% frames_lever_R(frames_lever_R < (n_frames*n_channels) | frames_lever_R > (numframes+1-n_frames*n_channels)) = [];
% num_leverpress = length(frames_lever_R);
% 
% w = waitbar(0,'Calculating ETA for RIGHT');
% for leverpress = 1:num_leverpress;
%     for frame_current = frames_lever_R(leverpress)-n_frames*n_channels:n_channels:frames_lever_R(leverpress)+n_frames*n_channels;
%         ETA_frame = ceil((frame_current-(frames_lever_R(leverpress)-n_frames*n_channels)+1)/n_channels);
%         im_R(:,:,ETA_frame)=im_R(:,:,ETA_frame)+double(imread([img_pathname img_filename],'tiff',frame_current))./num_leverpress;
%     end
%     waitbar(leverpress/num_leverpress,w,'Calculating ETA for RIGHT');
% end
% close(w);
% 
previewmov = figure;
h_waitbar = waitbar(0,'Frames');

% make tiff file
[img_filepath, img_filebase] = fileparts(img_filename);

for tif_frameL = 1:size(im_L,3)
    frameL = uint16(im_L(:,:,tif_frameL));
    imwrite(frameL,[img_pathname filesep img_filebase  '_etaL.tif'],'tif','Compression','none','WriteMode','append');
end
% for tif_frameR = 1:size(im_R,3)
%     frameR = uint16(im_R(:,:,tif_frameR));
%     imwrite(frameR,[img_pathname filesep img_filebase  '_etaR.tif'],'tif','Compression','none','WriteMode','append');
% end

% make movie for matlab
for u = 1:size(im_L,3)
    imagesc([im_L(:,:,u)]); %im_R(:,:,u)]);
    colormap(gray);
    im_mov(u) = getframe(gcf);
    waitbar(u/size(im_L,3), h_waitbar, ['Frames: ' num2str(u) '/' num2str(numframes)]);
end
close(h_waitbar, previewmov);

close all