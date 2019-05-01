frame_search = 10;
event_decay = cell(size(roi_trace_long,1),1);
for curr_roi = 1:size(roi_trace_long,1);
    if ~isempty(roi_peak_max{curr_roi});
    peak_frames = roi_peak_max{curr_roi}(:,1);
    event_decay{curr_roi} = zeros(length(peak_frames),frame_search+1);
    for curr_peak = 1:length(peak_frames)
        if peak_frames(curr_peak)-frame_search > 1 && ...
            peak_frames(curr_peak)+frame_search < size(roi_trace_long,2)
        event_decay_range = peak_frames(curr_peak):peak_frames(curr_peak)+frame_search;
        event_decay{curr_roi}(curr_peak,:) = (roi_trace_long(curr_roi,event_decay_range) ...
             ./roi_trace_long(curr_roi,peak_frames(curr_peak))) + 1;
        end
    end
    end
end
x_ms = [0:frame_search].*(1.24*128);
figure;hold on;
median_event_decay = [];
for i = 1:length(event_decay)
    if ~isempty(event_decay{i});
    plot(x_ms,median(event_decay{i}',2));
    median_event_decay = [median_event_decay; median(event_decay{i}',2)'];
    end
end
line([xlim],[1/exp(1)+1 1/exp(1)+1],'linestyle','--','color','k')

tau = [];
for i = 1:size(median_event_decay,1);
    [p] = polyfit(x_ms,log(abs(median_event_decay(i,:)-1)),1);
    tau(i) = -1/(p(1));
end
tau2 = [];
[p2] = polyfit(x_ms,log(abs(mean(median_event_decay,1)-1)),1)
tau2 = -1/(p2(1));

figure; hold on;
good_curves = size(median_event_decay,1);
plot(x_ms,mean(median_event_decay));
plot(x_ms,mean(median_event_decay) + (std(median_event_decay)./sqrt(good_curves)).*1.98,'--');
plot(x_ms,mean(median_event_decay) - (std(median_event_decay)./sqrt(good_curves)).*1.98,'--');
line([mean(tau) mean(tau)],[ylim],'color','k')
line([mean(tau)+(std(tau)/sqrt(good_curves))*1.98 mean(tau)+(std(tau)/sqrt(good_curves))*1.98],[ylim],'color','k','linestyle','--')
line([mean(tau)-(std(tau)/sqrt(good_curves))*1.98 mean(tau)-(std(tau)/sqrt(good_curves))*1.98],[ylim],'color','k','linestyle','--')
line([tau2 tau2],[ylim],'color','r')
line([xlim],[1/exp(1)+1 1/exp(1)+1],'color','k')


    
    