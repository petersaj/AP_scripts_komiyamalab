% Plots baseline estimation parameters through smoothing, requires
% roi_trace_long

%% Find optimal solution for fitting via smoothing
size_file = 4000;

smooth_types = {'sgolay','loess'}; %'rlowess','rloess',
smooth_values = [70:2:120];

results_mean = cell(length(smooth_types),length(smooth_values));
results_std = cell(length(smooth_types),length(smooth_values));
w = waitbar(0);
for curr_smooth_type_num = 1:length(smooth_types);
    curr_smooth_type = smooth_types{curr_smooth_type_num};
    for curr_smooth_value_num = 1:length(smooth_values);
        curr_smooth_value = smooth_values(curr_smooth_value_num);
        roi_trace_long_smooth = zeros(size(roi_trace_long));
        for i = 1:size(roi_trace_long,1);
            for j = 1:size(roi_trace_long,2)/size_file;
                roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),curr_smooth_value,curr_smooth_type);
            end
        end
        roi_trace_diff = abs(roi_trace_long - roi_trace_long_smooth);
        results_mean{curr_smooth_type_num,curr_smooth_value_num} = mean(roi_trace_diff,2);
        results_std{curr_smooth_type_num,curr_smooth_value_num} = sqrt(sum((abs(roi_trace_long - roi_trace_long_smooth).^2)/size(roi_trace_long_smooth,2),2));
        waitbar(curr_smooth_value_num/length(smooth_values),w,['Smooth value: ' num2str(curr_smooth_type_num),'/' num2str(length(smooth_types))]);
    end
end
close(w);

%% Plot optimal value test
figure;
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;
color_type = {'k','r','b','g'};
smooth_types_plot = {'sgolay','loess'};

for curr_smooth_type_num = 1:length(smooth_types_plot);
    curr_smooth_type = smooth_types_plot{curr_smooth_type_num};
    smooth_type_match = find(strcmp(smooth_types_plot{curr_smooth_type_num},smooth_types) == 1);
    
    subplot(2,1,1)
    all_std = zeros(size(roi_trace_long,1),size(results_std,2));
    for i = 1:size(roi_trace_long,1)
        all_std(i,:) = cellfun(@(x) x(i), results_std(smooth_type_match,:));
    end
    errorbar(mean(all_std),std(all_std),color_type{curr_smooth_type_num},'linewidth',2)
    set(gca,'XTick',[1:length(smooth_values)]);
    set(gca,'XTickLabel',smooth_values);
    xlabel('Smoothing span','FontSize',20)
    ylabel('Estimated noise std','FontSize',20)
    
    subplot(2,1,2)
    all_mean = zeros(size(roi_trace_long,1),size(results_mean,2));
    for i = 1:size(roi_trace_long,1)
        all_mean(i,:) = cellfun(@(x) x(i), results_mean(smooth_type_match,:));
    end
    errorbar(mean(all_mean),std(all_mean),color_type{curr_smooth_type_num},'linewidth',2)
    set(gca,'XTick',[1:length(smooth_values)]);
    set(gca,'XTickLabel',smooth_values);
    xlabel('Smoothing span','FontSize',20)
    ylabel('Abs. difference','FontSize',20)
end
subplot(2,1,1);
legend('location','NW',smooth_types_plot,'FontSize',15)

%% Compare individual traces of loess and sgolay

size_file = 4000;
size_smooth = 30;
roi = 10;

figure;

smooth_type = 'lowess';
roi_trace_long_smooth = zeros(size(roi_trace_long));

for i = roi;
    for j = 1:size(roi_trace_long,2)/size_file;
        roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end

subplot(2,1,1); hold on;
plot(roi_trace_long(roi,:),'k','linewidth',2)
plot(roi_trace_long_smooth(roi,:),'r','linewidth',2)
title([smooth_type ': ' num2str(size_smooth)],'FontSize',20);

smooth_type = 'loess';
roi_trace_long_smooth = zeros(size(roi_trace_long));
for i = roi;
    for j = 1:size(roi_trace_long,2)/size_file;
        roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end

subplot(2,1,2); hold on;
plot(roi_trace_long(roi,:),'k','linewidth',2)
plot(roi_trace_long_smooth(roi,:),'r','linewidth',2)
title([smooth_type ': ' num2str(size_smooth)],'FontSize',20);

linkaxes;

%% Compare std of gaussian method to std of smoothing method
size_file = 4000;
size_smooth = 90;
smooth_type = 'loess';

% noise through smoothing
roi_trace_long_smooth = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    for j = 1:size(roi_trace_long,2)/size_file;
        roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end
smooth_std = zeros(size(roi_trace_long,1),1);
smooth_std = sqrt(sum((abs(roi_trace_long - roi_trace_long_smooth).^2)/size(roi_trace_long_smooth,2),2));

smooth_baseline = zeros(size(roi_trace_long,1),1);
for i = 1:size(roi_trace_long,1);
    options = optimset('Algorithm','interior-point');
    lb = min(roi_trace_long(i,:));
    ub = max(roi_trace_long(i,:));
    starting = lb+(ub-lb)/2;
    
    estimate_baseline_cutoff = fminsearch(@(x) abs(smooth_std(i)-std(roi_trace_long(i,(roi_trace_long(i,:) < x)))),starting);
    estimate_baseline = mean(roi_trace_long(i,(roi_trace_long(i,:) < estimate_baseline_cutoff)));
    smooth_baseline(i) = estimate_baseline;
end

% noise through gaussian estimation
gauss_std = zeros(size(roi_trace_long,1),1);
gauss_baseline = zeros(size(roi_trace_long,1),1);
for curr_roi = 1:size(roi_trace_long,1)
    cellTrace = roi_trace_long(curr_roi,:);
    if ~any(isfinite(cellTrace));
        continue
    end
    try
        norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',2); % there are two hidden curves
    catch me
        norm_fit_obj = gmdistribution.fit(roi_trace_long(curr_roi,:)',1);
    end
    baseline_curve = find(norm_fit_obj.mu == min(norm_fit_obj.mu));
    gauss_std(curr_roi) = sqrt(norm_fit_obj.Sigma(baseline_curve(1)));
    dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
    x = min(cellTrace):dx:max(cellTrace);
    x = x(1:end-1); % it adds an extra, not sure why but can fix later
    try
        norm_fit = pdf(norm_fit_obj,x');
    catch me
        continue
    end
    gauss_baseline(curr_roi) = x(norm_fit == max(norm_fit));
end

figure;hold on
plot(smooth_std,'k','linewidth',2)
plot(gauss_std,'--r','linewidth',2)
legend({'Smoothing','Gaussian-mixture'});
xlabel('ROI','FontSize',20);
ylabel('Baseline STD.','FontSize',20)

figure;hold on
plot(smooth_baseline,'k','linewidth',2)
plot(gauss_baseline,'--r','linewidth',2)
legend({'Smoothing','Gaussian-mixture'});
xlabel('ROI','FontSize',20);
ylabel('Baseline','FontSize',20)

%% Calculate baseline by finding gaussian fit
tic
size_file = 4000;
size_smooth = 30;
smooth_type = 'loess';

% Estimate noise through smoothing
roi_trace_long_smooth = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    for j = 1:size(roi_trace_long,2)/size_file;
        roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end
smooth_std = zeros(size(roi_trace_long,1),1);
smooth_std = sqrt(sum((abs(roi_trace_long - roi_trace_long_smooth).^2)/size(roi_trace_long_smooth,2),2));

% Estimate baseline: threshold trace, mirror below threshold to above
% threshold, match standard deviation with estimated noise
roi_trace_df = zeros(size(roi_trace_long));
estimate_baseline_cutoff_all = zeros(size(roi_trace_long,1),1);
estimate_baseline_all = zeros(size(roi_trace_long,1),1);
for roi = 1:size(roi_trace_long,1);
    options = optimset('Algorithm','interior-point');
    lb = min(roi_trace_long(roi,:));
    ub = max(roi_trace_long(roi,:));
    starting = lb+(ub-lb)/2;
    
    %     % this function minimizes the absolute difference between the estimated
    %     % standard deviation, and the standard deviation of the mirrored trace
    %     % below threshold
    %     estimate_baseline = fminsearch(@(x) ...
    %         abs(smooth_std(i)- ...
    %         std([(roi_trace_long(i,(roi_trace_long(i,:) < x))-x) ...
    %         -(roi_trace_long(i,(roi_trace_long(i,:) < x))-x)])),starting);
    %
    %%% THIS IS FOR THE MODE RESTIMATION
    % this function minimizes the absolute difference between the estimated
    % standard deviation, and the standard deviation of the mirrored trace
    % below threshold
    estimate_baseline = fminsearch(@(x) ...
        AP_fitBaselineGaussian(x,smooth_std(roi),roi_trace_long(roi,:)),starting);
    
    % correct the baseline by looking for the mode of the cutoff assumed to
    % be just baseline (get better mode estimate by rounding baseline_trace
    % to the nearest interger)
    baseline_ub = estimate_baseline + 2*smooth_std(roi);
    baseline_lb = estimate_baseline - 2*smooth_std(roi);
    baseline_trace = roi_trace_long(roi,(roi_trace_long(roi,:) > baseline_lb &...
        roi_trace_long(roi,:) < baseline_ub));
    baseline_mode = mode(round(baseline_trace));
    estimate_baseline_all(roi) = baseline_mode;
    %%%
    %     estimate_baseline_all(i) = estimate_baseline;
    roi_trace_df(roi,:) = (roi_trace_long(roi,:)-estimate_baseline)/estimate_baseline;
end
toc
%% Check histogram for baseline
for roi = 20:30
    figure;
    subplot(2,1,1); hold on;
    plot(roi_trace_long(roi,:),'k');
    plot(roi_trace_long_smooth(roi,:),'r');
    line(xlim, [estimate_baseline_all(roi) estimate_baseline_all(roi)],'color','b', 'linestyle','-','linewidth',2);
    %plot(spike_frames{roi},roi_trace_long(roi,spike_frames{roi}),'.m','MarkerSize',30)
    % plot dotted lines for STDs
    line(xlim, [estimate_baseline_all(roi)+10*smooth_std(i) estimate_baseline_all(roi)+10*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [estimate_baseline_all(roi)+20*smooth_std(i) estimate_baseline_all(roi)+20*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [estimate_baseline_all(roi)+30*smooth_std(i) estimate_baseline_all(roi)+30*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [estimate_baseline_all(roi)+40*smooth_std(i) estimate_baseline_all(roi)+40*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [estimate_baseline_all(roi)+50*smooth_std(i) estimate_baseline_all(roi)+50*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    
    % plot lines for file separations
    x_plot = xlim;
    a = (estimate_baseline_all(roi)+smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
    b = (estimate_baseline_all(roi)-smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
    jbfill([x_plot(1)+1:x_plot(2)],a,b,'b','none',0,0.5);
    a = (estimate_baseline_all(roi)+2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
    b = (estimate_baseline_all(roi)-2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
    jbfill([x_plot(1)+1:x_plot(2)],a,b,'g','none',0,0.5);
    
    subplot(2,1,2); hold on;
    baseline_ub = estimate_baseline_all(roi) + 3*smooth_std(roi);
    baseline_lb = estimate_baseline_all(roi) - 3*smooth_std(roi);
    baseline_trace = roi_trace_long(roi,(roi_trace_long(roi,:) > baseline_lb &...
        roi_trace_long(roi,:) < baseline_ub));
    hist(baseline_trace,100);
    line([estimate_baseline_all(roi) estimate_baseline_all(roi)],ylim,'color','k','linewidth',2)
    line([mode(round(baseline_trace)) mode(round(baseline_trace))],ylim,'color','r','linewidth',2)
end

%% Get an estimate of (large) calcium transients

%peakstd = [1.5*ones(1,5) 1*ones(1,10)]; % definitions for max trigger, in stds
peak_std = [5*ones(1,10)]; % definitions for max trigger, in stds
spike_frames = {zeros(size(roi_trace_long,1))};
movie_frames = 4000;
peak_back = 10;
local_max = 30;
drop_percentage = .8;

parfor roi = 1:size(roi_trace_long,1);
    minpeak = smooth_std(roi)*peak_std;
    spike_frames{roi} = ap_schmittDetect(roi_trace_long_smooth(roi,:),minpeak,peak_back,local_max,drop_percentage);
    % eliminate peaks at file junctions because they're mistimed
    if ~isempty(peak_frames)
        file_stitch = [];
        file_stitch = find(mod(peak_frames-1,movie_frames) == 0);
        spike_frames{roi}(file_stitch) = [];
    end
    disp([num2str(roi) '/' num2str(size(roi_trace_long,1))]);
end

%% Calculate moving baseline by finding gaussian fit within window
tic

baseline_window = 6000;
baseline_window_step = 2000;
min_baseline_length = 300;
max_consec = 20;
size_file = 4000;

size_smooth = 90;
smooth_type = 'loess';

% Estimate noise through smoothing
roi_trace_long_smooth = zeros(size(roi_trace_long));
for i = 1:size(roi_trace_long,1);
    for j = 1:size(roi_trace_long,2)/size_file;
        roi_trace_long_smooth(i,size_file*(j-1)+1:size_file*(j-1)+size_file) = smooth(roi_trace_long(i,size_file*(j-1)+1:size_file*(j-1)+size_file),size_smooth,smooth_type);
    end
end
smooth_std = zeros(size(roi_trace_long,1),1);
smooth_std = sqrt(sum((abs(roi_trace_long - roi_trace_long_smooth).^2)/size(roi_trace_long_smooth,2),2));

% Estimate baseline: threshold trace, mirror below threshold to above
% threshold, match standard deviation with estimated noise
roi_trace_df = zeros(size(roi_trace_long));

% make a baseline matrix: all points not estimated are NaN
estimate_baseline_all = nan(size(roi_trace_long));
baseline_mode_all = nan(size(roi_trace_long));
baseline_trace_indx_all = cell(size(roi_trace_long,1),1);

for roi = 1:size(roi_trace_long,1);
    
    curr_frame = baseline_window/2+1;
    
    while curr_frame < size(roi_trace_long,2)-baseline_window/2;
        
        curr_frame_range = [curr_frame-baseline_window/2:curr_frame+baseline_window/2];
        curr_trace = roi_trace_long(roi,curr_frame_range);
        
        options = optimset('Algorithm','interior-point');
        options = optimset('MaxFunEvals',2000);
        lb = min(roi_trace_long(roi,:));
        ub = max(roi_trace_long(roi,:));
        starting = lb+(ub-lb)/2;
        
        %     % this function minimizes the absolute difference between the estimated
        %     % standard deviation, and the standard deviation of the mirrored trace
        %     % below threshold
        %     estimate_baseline = fminsearch(@(x) ...
        %         abs(smooth_std(i)- ...
        %         std([(roi_trace_long(i,(roi_trace_long(i,:) < x))-x) ...
        %         -(roi_trace_long(i,(roi_trace_long(i,:) < x))-x)])),starting);
        %
        %%% THIS IS FOR THE MODE RESTIMATION
        % this function minimizes the absolute difference between the estimated
        % standard deviation, and the standard deviation of the mirrored trace
        % below threshold
        estimate_baseline = fminsearch(@(x) ...
            AP_fitBaselineGaussian(x,smooth_std(roi),curr_trace),starting);
        
        % correct the baseline by looking for the mode of the cutoff assumed to
        % be just baseline (get better mode estimate by rounding baseline_trace
        % to the nearest interger)
        baseline_ub = estimate_baseline + 2*smooth_std(roi);
        baseline_lb = estimate_baseline - 2*smooth_std(roi);
        baseline_trace = curr_trace(curr_trace > baseline_lb &...
            curr_trace < baseline_ub);
        baseline_trace_indx = find(curr_trace > baseline_lb &...
            curr_trace < baseline_ub) + curr_frame-baseline_window/2;
        baseline_mode = mode(round(baseline_trace));
        
        baseline_trace_below_indx = find(curr_trace < baseline_mode);
        
        % find longest consecutive run
        consec_below = [0 diff(baseline_trace_below_indx) == 1 0];
        p = find(~consec_below);%find 0
        longest_consec = max(diff(p)-1);
        
        % 1) ignore if number of baseline points below minimum
        % 2) ignore if based off of too many consecutive samples below thresh
        if length(baseline_trace) > min_baseline_length &&...
                longest_consec < max_consec
            baseline_mode_all(roi,curr_frame) = baseline_mode;
            estimate_baseline_all(roi,curr_frame) = estimate_baseline;
            baseline_trace_indx_all{roi} = [baseline_trace_indx_all{roi} baseline_trace_indx];
        end
                
        curr_frame = curr_frame + baseline_window_step;
        
    end
end

% if no good estimate, then just try to estimate by the whole trace
for roi = 1:size(roi_trace_long,1);
    if ~all(isnan(estimate_baseline_all(roi,:)))
        continue
    end
    disp(['WARNING: No good baseline for ' num2str(roi) ', estimating from whole trace']);
     options = optimset('Algorithm','interior-point');
        options = optimset('MaxFunEvals',2000);
        lb = min(roi_trace_long(roi,:));
        ub = max(roi_trace_long(roi,:));
        starting = lb+(ub-lb)/2;
        
        curr_trace = roi_trace_long(roi,:);
        curr_frame = round(size(roi_trace_long,2)/2);
        %%% THIS IS FOR THE MODE RESTIMATION
        % this function minimizes the absolute difference between the estimated
        % standard deviation, and the standard deviation of the mirrored trace
        % below threshold
        estimate_baseline = fminsearch(@(x) ...
            AP_fitBaselineGaussian(x,smooth_std(roi),curr_trace),starting);
        
        % correct the baseline by looking for the mode of the cutoff assumed to
        % be just baseline (get better mode estimate by rounding baseline_trace
        % to the nearest interger)
        baseline_ub = estimate_baseline + 2*smooth_std(roi);
        baseline_lb = estimate_baseline - 2*smooth_std(roi);
        baseline_trace = curr_trace(curr_trace > baseline_lb &...
            curr_trace < baseline_ub);
        baseline_trace_indx = find(curr_trace > baseline_lb &...
            curr_trace < baseline_ub);
        baseline_trace_indx_all{roi} = [baseline_trace_indx_all{roi} baseline_trace_indx];
        baseline_mode = mode(round(baseline_trace));
        baseline_mode_all(roi,curr_frame) = baseline_mode;
        estimate_baseline_all(roi,curr_frame) = estimate_baseline;
end

% interpolate baseline between estimated baseline points
estimate_baseline_discrete = estimate_baseline_all;
for roi = 1:size(roi_trace_long,1);
    if all(isnan(estimate_baseline_all(roi,:)))
        disp(['WARNING: No good estimate for ROI ' num2str(roi)]);
        continue
    end
    frame_estimates_x = [];
    frame_estimates_x = find(isnan(estimate_baseline_all(roi,:)) == 0);
    frame_estimates_y = [];
    frame_estimates_y = estimate_baseline_all(roi,frame_estimates_x);
    % assume baseline start and end are the same as the first and last
    % estimates of baseline
    frame_estimates_x = [1 frame_estimates_x size(roi_trace_long,2)];
    frame_estimates_y = [frame_estimates_y(1) frame_estimates_y frame_estimates_y(end)];
    estimate_baseline_all(roi,:) = ...
        interp1(frame_estimates_x,frame_estimates_y,1:size(roi_trace_long,2));
end

% interpolate baseline between estimated baseline mode points
baseline_mode_discrete = baseline_mode_all;
for roi = 1:size(roi_trace_long,1);
    if all(isnan(baseline_mode_all(roi,:)))
        continue
    end
    frame_estimates_x = [];
    frame_estimates_x = find(isnan(baseline_mode_all(roi,:)) == 0);
    frame_estimates_y = [];
    frame_estimates_y = baseline_mode_all(roi,frame_estimates_x);
    % assume baseline start and end are the same as the first and last
    % estimates of baseline
    frame_estimates_x = [1 frame_estimates_x size(roi_trace_long,2)];
    frame_estimates_y = [frame_estimates_y(1) frame_estimates_y frame_estimates_y(end)];
    baseline_mode_all(roi,:) = ...
        interp1(frame_estimates_x,frame_estimates_y,1:size(roi_trace_long,2));
end

% this is to actually correct the raw trace
roi_trace_df = (roi_trace_long-baseline_mode_all)./baseline_mode_all;

toc

%% Check histogram for moving
for roi = 5
    figure;
    %subplot(2,1,1); 
    hold on;
    plot(baseline_trace_indx_all{roi},roi_trace_long(roi,baseline_trace_indx_all{roi}),'.b');
    plot(roi_trace_long(roi,:),'k');
    plot(roi_trace_long_smooth(roi,:),'r');
    plot(estimate_baseline_all(roi,:),'-b','linewidth',2);
    plot(estimate_baseline_discrete(roi,:),'.m','MarkerSize',10);
    
    plot(baseline_mode_all(roi,:),'-g','linewidth',2);
    plot(baseline_mode_discrete(roi,:),'.c','MarkerSize',10)

    % plot dotted lines for STDs (from avg baseline)
    curr_mean_baseline = mean(estimate_baseline_all(roi,:));
    line(xlim, [curr_mean_baseline+10*smooth_std(i) curr_mean_baseline+10*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+20*smooth_std(i) curr_mean_baseline+20*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+30*smooth_std(i) curr_mean_baseline+30*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+40*smooth_std(i) curr_mean_baseline+40*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+50*smooth_std(i) curr_mean_baseline+50*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    
%     % fill area of baseline for std
%     x_plot = xlim;
%     a = (estimate_baseline_all(roi)+smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     b = (estimate_baseline_all(roi)-smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     jbfill([x_plot(1)+1:x_plot(2)],a,b,'b','none',0,0.5);
%     a = (estimate_baseline_all(roi)+2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     b = (estimate_baseline_all(roi)-2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     jbfill([x_plot(1)+1:x_plot(2)],a,b,'g','none',0,0.5);
    
%     subplot(2,1,2); hold on;
%     baseline_ub = estimate_baseline_all(roi) + 2*smooth_std(roi);
%     baseline_lb = estimate_baseline_all(roi) - 2*smooth_std(roi);
%     baseline_trace = roi_trace_long(roi,(roi_trace_long(roi,:) > baseline_lb &...
%         roi_trace_long(roi,:) < baseline_ub));
%     hist(baseline_trace,100);
%     line([estimate_baseline_all(roi) estimate_baseline_all(roi)],ylim,'color','k','linewidth',2)
%     line([mode(round(baseline_trace)) mode(round(baseline_trace))],ylim,'color','r','linewidth',2)
end

%% Plot baseline estimates
for roi = 67
    figure;
    %subplot(2,1,1); 
    hold on;
    plot(baseline_trace_indx_all{roi},roi_trace_long(roi,baseline_trace_indx_all{roi}),'.b');
    plot(roi_trace_long(roi,:),'k');
    plot(roi_trace_long_smooth(roi,:),'r');
    plot(estimate_baseline_all(roi,:),'-b','linewidth',2);
    plot(estimate_baseline_discrete(roi,:),'.m','MarkerSize',10);
    
    plot(baseline_mode_all(roi,:),'-g','linewidth',2);
    plot(baseline_mode_discrete(roi,:),'.c','MarkerSize',10)

    % plot dotted lines for STDs (from avg baseline)
    curr_mean_baseline = mean(estimate_baseline_all(roi,:));
    line(xlim, [curr_mean_baseline+10*smooth_std(i) curr_mean_baseline+10*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+20*smooth_std(i) curr_mean_baseline+20*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+30*smooth_std(i) curr_mean_baseline+30*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+40*smooth_std(i) curr_mean_baseline+40*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    line(xlim, [curr_mean_baseline+50*smooth_std(i) curr_mean_baseline+50*smooth_std(i)],'color','k', 'linestyle','--','linewidth',1);
    
%     % fill area of baseline for std
%     x_plot = xlim;
%     a = (estimate_baseline_all(roi)+smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     b = (estimate_baseline_all(roi)-smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     jbfill([x_plot(1)+1:x_plot(2)],a,b,'b','none',0,0.5);
%     a = (estimate_baseline_all(roi)+2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     b = (estimate_baseline_all(roi)-2*smooth_std(roi))*ones(1,round(x_plot(2)-x_plot(1)));
%     jbfill([x_plot(1)+1:x_plot(2)],a,b,'g','none',0,0.5);
    
%     subplot(2,1,2); hold on;
%     baseline_ub = estimate_baseline_all(roi) + 2*smooth_std(roi);
%     baseline_lb = estimate_baseline_all(roi) - 2*smooth_std(roi);
%     baseline_trace = roi_trace_long(roi,(roi_trace_long(roi,:) > baseline_lb &...
%         roi_trace_long(roi,:) < baseline_ub));
%     hist(baseline_trace,100);
%     line([estimate_baseline_all(roi) estimate_baseline_all(roi)],ylim,'color','k','linewidth',2)
%     line([mode(round(baseline_trace)) mode(round(baseline_trace))],ylim,'color','r','linewidth',2)
end
