figure
n_rois = size(roi_trace,1);
for i = 1:n_rois
    hold on
    subplot(n_rois,1,i)
    plot(roi_trace(i,:)')
end

roi_corrcoef = corrcoef(roi_trace');
imagesc(roi_corrcoef);
colormap(gray);
