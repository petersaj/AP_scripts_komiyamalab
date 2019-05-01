function [peak_frames peak_amplitudes] = AP_trace_peak(roi_trace_df)
% [peak_frames peak_amplitudes] = AP_trace_peak(roi_trace_df)
% Get peaks in a calcium trace

%peakstd = [5 5]; % definitions for max trigger, in stds - NOT USED atm
peak_back = 30;
local_max = 30;
drop_percentage = .8;
peak_frames = {};
peak_amplitudes = {};

num_rois = size(roi_trace_df,1);

% this is to catch a dumb bug: if ROI is offscreen, just make noise
for i = 1:num_rois
    if all(roi_trace_df(i,:) == roi_trace_df(i,1))
        roi_trace_df(i,:) = rand(size(roi_trace_df(i,:)));
    end
end

% Smooth trace for easier detection
size_smooth = 30;
smooth_type = 'loess';

roi_trace_df_smooth = zeros(size(roi_trace_df));
for i = 1:size(roi_trace_df,1)
    if any(isnan(roi_trace_df(i,:)))
        continue
    end
    roi_trace_df_smooth(i,:) = smooth(roi_trace_df(i,:),size_smooth,smooth_type);
end

smooth_std = zeros(size(roi_trace_df,1),1);
smooth_std = sqrt(sum((abs(roi_trace_df - roi_trace_df_smooth).^2)/size(roi_trace_df_smooth,2),2));

disp('Finding calcium trace peaks...')
fprintf('%2d',0);
tic
for i = 1:num_rois
    % Estimate noise through smoothing
    %curr_trace_std = smooth_std(i);
    %minpeak = curr_trace_std*peakstd;
    minpeak = 0.2;
    % find peaks in trace
    [peak_frames{i} peak_amplitudes{i}] = AP_schmittDetect2(roi_trace_df_smooth(i,:),minpeak,peak_back,local_max,drop_percentage);
    fprintf('%c%c%2d',8,8,round(100*i/num_rois));
end
toc

% % eliminate peaks at loop junctions
% for curr_roi = 1:num_rois;
%     file_frames = cellfun(@length,roi_trace_long_split);
%     file_frames_start = [1;cumsum(file_frames)+1];
%     loop_frames_start = file_frames_start(1:8:end);
%     for i = 1:length(loop_frames_start);
%         if i == 1
%             spike_train(:,loop_frames_start(i):loop_frames_start(i)+20) = 0;
%         else
%             spike_train(:,loop_frames_start(i)-20:loop_frames_start(i)+20) = 0;
%         end
%     end
% end

% % eliminate peaks at file junctions because they're mistimed
% for i = 1:num_rois
%     if ~isempty(spike_frames{i})
%         file_stitch = [];
%         file_stitch = find(mod(spike_frames{i}-1,movie_frames) == 0);
%         peak_frames(file_stitch,:) = [];
%     end
% end

%%%% what's been kind of working
%[peak_frames] = ap_schmittDetect(roi_trace_long(24,:),minpeak,5,20,0.5);
%%%%