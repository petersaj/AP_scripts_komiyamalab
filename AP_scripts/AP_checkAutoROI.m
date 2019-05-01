
%% visually run through autosort ROIs, get which ones are good

f = figure;

% go through to get rid of ROIs that look clearly off
bad_roi = [];
i = 1;

fig_dim = 5;
while i <= ceil(size(polygon.autosort,3)/5);
    
    for j = 1:fig_dim^2;
        if ((i-1)*fig_dim^2+j) <= size(polygon.autosort,3)
            subplot(fig_dim,fig_dim,j);
            imagesc(polygon.autosort(:,:,(i-1)*fig_dim^2+j));
            colormap(gray);
            axis off
            title(num2str((i-1)*fig_dim^2+j));         
        else
            subplot(fig_dim,fig_dim,j);
            cla
            title('');
        end
    end
    r = input('Bad ROI? or b/f = scroll, q = quit: ','s');
    if strcmp(r,'b');
        i = i-1;
    elseif strcmp(r,'f');
        i = i+1;
    elseif strcmp(r,'q');
        break
        close(f);
    else
        bad_roi = [bad_roi;str2num(r)];
    end
end

roi_trace_long(bad_roi,:) = [];
polygon.handles(bad_roi) = [];
polygon.ROI(bad_roi) = [];
polygon.autosort(:,:,bad_roi) = [];

%% Check corrcoef, check traces over 0.5 correlation


% !! maybe not worry about it if the ROIs don't touch eachother

% get half of the correlation coefficients between traces
trace_corrcoef = tril(corrcoef(roi_trace_long'),-1);
[check_row check_col] = find(trace_corrcoef > 0.5);
unique_row = unique(check_row);
unique_col = unique(check_col);
% loop through, check out roi placement and trace and pick good/bad
h = figure;
for i = 1:length(unique_row)
    clf
    check_indx = find(check_row == unique_row(i));
    % plot all the traces
    subplot(length(check_indx)+1,2,1);
    plot(roi_trace_long(unique_row(i),:)); % plot the first, from check_rows
    title(num2str(unique_row(i)))
    axis off
    for num_check = 1:length(check_indx) % plot all the others
        subplot(length(check_indx)+1,2,(num_check*2+1));
        plot(roi_trace_long(check_col(check_indx(num_check)),:))
        title(num2str(check_col(check_indx(num_check))));
        axis off
    end
    linkaxes;
    ylim([-1 2])
    
    % plot all ROIs
    subplot(length(check_indx)+1,2,2);
    imagesc(polygon.autosort(:,:,unique_row(i))); % plot the first, from check_rows
    for num_check = 1:length(check_indx) % plot all the others
        subplot(length(check_indx)+1,2,(num_check*2+2));
        imagesc(polygon.autosort(:,:,check_col(check_indx(num_check))))
    end
    
    % input for which ones are bad
    r = input('Bad ROI? q = quit: ','s');
    if strcmp(r,'q');
        break
        close(f);
    else
        bad_roi = [bad_roi;str2num(r)'];
    end
    
    % replace bad roi_trace_long with NaNs to keep ROI numbers
    roi_trace_long(str2num(r)',:) = NaN;
    
end
close(h);
%%
roi_trace_long(bad_roi,:) = [];
polygon.handles(bad_roi) = [];
polygon.ROI(bad_roi) = [];
polygon.autosort(:,:,bad_roi) = [];
