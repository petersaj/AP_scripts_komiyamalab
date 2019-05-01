function [roi_bgcorrected_dFoF, intermediates] = AM_get_dFoF_from_img(img, roi, bg, framerate, bg_cutoff)
% selecting bg_pixels is adapted from AP_cell_ROI_pro.m.
% (names of the variables are modified)
%                       AM 5/16/2013
    if(nargin<4)
        framerate = 6.3004;
    end
    
    if(nargin<5)
        bg_cutoff = 3;  %% different from the default value of AP_cell_ROI_pro
    end


    % ROI = bsxfun(@minus, ROI_pol, topleft - 1);
    % bgROI = bsxfun(@minus, bgROI_pol, topleft - 1);

    N = size(img,1);
    M = size(img,2);
    L = size(img,3);

    im = reshape(img,[N*M, L]);

        % get ROI trace
    roi_mask = poly2mask(roi(:,1),roi(:,2),N,M);
    roi_trace = mean(im(roi_mask(:),:));
    % roi_pixel_trace = im(roi_mask(:),:);

    % get background trace
    bg_mask = poly2mask(bg(:,1),bg(:,2),N,M);
    
    
    % orig_bgROI_mask = poly2mask(bgROI(:,1),bgROI(:,2),N,M);
    bg_pixel_trace = im(bg_mask(:),:);

    bg_pixel_diff = ...
        bsxfun(@minus,double(bg_pixel_trace),double(roi_trace));

    % define outliers by the distribution of differences. There is a
    % gaussian (sometimes a little buried) - so just estimate the mean with
    % mode and the std

    % get mode rounded to nearest integer
    bg_diff_mode = mode(round(bg_pixel_diff(:)*1))/1;
    bg_diff_std = std(bg_pixel_diff(:));
    
%     bg_trace_outliers = ...
%         bsxfun(@gt,bg_pixel_diff, bg_diff_mode + bg_cutoff*bg_diff_std);
    bg_trace_max = max(bg_pixel_diff,[],2);
    
    bg_mask_index = find(bg_mask);
%    bg_pixel_outliers = bg_mask_index(any(bg_trace_outliers,2));
%     bg_pixel_outliers = bg_mask_index(bg_trace_max >  bg_diff_mode + bg_cutoff*bg_diff_std);
    bg_pixel_fixed = bg_mask_index(bg_trace_max <=  bg_diff_mode + bg_cutoff*bg_diff_std);
    
%     bg_mask(bg_pixel_outliers) = false;
    
    bg_mask_fixed = false(N,M);
    bg_mask_fixed(bg_pixel_fixed) = true;
    
    
    disp(['Selecting ' num2str(sum(sum(bg_mask_fixed))) '/' num2str(sum(sum(bg_mask)))]);
    
    MIN_N_PIXEL = 5;
    
    % if there are < 10 pixels, use whole bg
    if sum(bg_mask_fixed(:)) < MIN_N_PIXEL
%         bg_mask(bg_pixel_outliers) = true;
        
        if sum(bg_mask(:)) < MIN_N_PIXEL
            bg_mask_fixed = bg_mask;
            disp(['Warning: bgROI has less than ' num2str(MIN_N_PIXEL) ' pixels, using the whole ROI']);
        else
            [~, sorted_index] = sort(bg_trace_max);
            bg_pixel_fixed = bg_mask_index(sorted_index(1:MIN_N_PIXEL)); % 10 rois with smallest max response
            bg_mask_fixed = false(N,M);
            bg_mask_fixed(bg_pixel_fixed) = true;
            disp(['Warning: bgROI has less than ' num2str(MIN_N_PIXEL) ' pixels, increased threshold to ' num2str( ...
                (bg_trace_max(sorted_index(MIN_N_PIXEL))-bg_diff_mode)/bg_diff_std)]);
        end
        
        
    %     set(polygon.bghandles{roi},'color','r');
% 
%     % give a warning if < 50 pixels
%     elseif sum(bg_mask(:)) < 50
%         disp(['Warning: bgROI has less than 50 pixels ' ]);
%     %     set(polygon.bghandles{roi},'color',[1 0.5 0]);
% 
%     else
%     %     disp('ok')
%     %     set(polygon.bghandles{roi},'color','g');  

    end

    bg_trace = mean(im(bg_mask_fixed(:),:));

    [bg_dFoF, bg_baseline] = AP_baselineEstimation(bg_trace, framerate);
    bg_diff = bg_trace - bg_baseline;
    roi_bgcorrected_trace = roi_trace - bg_diff; 
    [roi_bgcorrected_dFoF, roi_bgcorrected_baseline] = AP_baselineEstimation(roi_bgcorrected_trace, framerate);
    
    intermediates.bg_mask = bg_mask_fixed;
    
    intermediates.bg_dFoF = bg_dFoF;
    intermediates.bg_baseline = bg_baseline;
    intermediates.bg_diff = bg_diff;
    intermediates.bg_trace = bg_trace;

    intermediates.roi_trace = roi_trace;
    
    intermediates.roi_bgcorrected_trace = roi_bgcorrected_trace;
    intermediates.roi_bgcorrected_dFoF = roi_bgcorrected_dFoF;
    intermediates.roi_bgcorrected_baseline = roi_bgcorrected_baseline;
    
end

%%
% clf
% subplot(211);
% plot([roi_trace' bg_trace' roi_bgcorrect' roi_corrected_baseline' roi_bg_baseline']);
% legend('ROI','BG ROI','Corrected','ROI baseline', 'BG ROI baseline')
% box off;
% legend('boxoff')
% subplot(212);
% plot([roi_corrected_dFoF' roi_bg_dFoF']);
% legend('ROI DFOF','BG ROI DFOF');
% box off;
% legend('boxoff')
% % saveas(gcf,fullfile(mousename,[filename '.png']),'png');
% toc
% close gcf
%%
% [u, s, v] = svd(im,'econ');
% 
% figure;
% for i = 1:3
%     subplot(2,4,i);
%     imagesc(reshape(u(:,i),size(img,1),size(img,2)));
%     title(num2str(s(i,i)));
% end
% subplot(2,4,4);
% plot(diag(s(1:50,1:50)),'.');
% subplot(2,1,2);
% plot(v(:,1:3));
% 
% 
% 
% 
% 
% 
