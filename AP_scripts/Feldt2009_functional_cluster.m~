% this functionally clusters cells as described in Feldt et al 2009


animals = {'AP71' 'AP72' 'AP73' 'AP74' 'AP75' 'AP76' 'AP78' 'AP79'};

curr_animal = 5;

animal = animals{curr_animal};

% load ROI labels and identify cell type
analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_roi_template'];
roilabel_name = [animal '_roilabels.roilabel'];
load([analysis_path filesep roilabel_name],'-MAT')
cells = 1:length(roi_labels);
gad_cells = cellfun(@(x) any(strcmp('gad',x)),roi_labels);
pyr_cells = ~gad_cells;
filled_cells = cellfun(@(x) any(strcmp('filled',x)),roi_labels);

pyr_unfilled = find(pyr_cells & ~filled_cells);
gad_unfilled = find(gad_cells & ~filled_cells);

peak_matrix = AP_caEvents(im.roi_trace_df,pyr_unfilled,gad_unfilled);

pyr_peaks = peak_matrix(pyr_unfilled,:);



%% Loop functional clustering algorithm until no significant similarities


% 1) Create matrix of pairwise similarity

pyr_peaks_idx = arrayfun(@(x) find(pyr_peaks(x,~isnan(pyr_peaks(x,:)))),1:size(pyr_peaks,1),'uni',false);
pyr_peaks_active = cellfun(@(x) ~isempty(x),pyr_peaks_idx);

pyr_peaks_amd = cellfun(@(y) cellfun(@(x) ...
    nanmean(min(abs(bsxfun(@minus,y,x')),[],2))/(size(x,2)/(length(x)+1)), ...
    pyr_peaks_idx(pyr_peaks_active)), ...
    pyr_peaks_idx(pyr_peaks_active),'uni',false);


gauss_filt = fspecial('gaussian',[1 100],20);
pyr_peaks_filt = conv2(pyr_peaks,gauss_filt,'same');
pyr_peaks_gauss_similarity = AP_itril(corrcoef(pyr_peaks_filt'),-1);

pyr_peaks_amd_norm = AP_itril((vertcat(pyr_peaks_amd{:})+vertcat(pyr_peaks_amd{:}))/2,-1);

% 2) Calculate 95% CI by jitter

% jitter and repeat, store 500 highest values for 95% conf
reps = 100;
similarity_95 = -Inf(length(pyr_peaks_gauss_similarity),round((reps*0.05)));
for rep = 1:reps
    % jitter data
    pyr_peaks_jitter_idx = ...
        cellfun(@(x) x+round(normrnd(0,30,1,length(x))),pyr_peaks_idx,'uni',false);
    
    pyr_peaks_jitter = zeros(size(pyr_peaks));
    for i = 1:size(pyr_peaks_jitter_idx)
        pyr_peaks_jitter(i,pyr_peaks_jitter_idx{i}) = 1;
    end
    pyr_peaks_jitter_filt = conv2(pyr_peaks_jitter,gauss_filt,'same');
    pyr_peaks_jitter_gauss_similarity = AP_itril(corrcoef(pyr_peaks_jitter_filt'),-1);
    
%     pyr_peaks_jitter_amd = cellfun(@(y) cellfun(@(x) ...
%     nanmean(min(abs(bsxfun(@minus,y,x')),[],2))/(size(x,2)/(length(x)+1)), ...
%     pyr_peaks_jitter(pyr_peaks_active)), ...
%     pyr_peaks_jitter(pyr_peaks_active),'uni',false);
% 
%     pyr_peaks_jitter_amd_norm = AP_itril(...
%         (vertcat(pyr_peaks_jitter_amd{:})+vertcat(pyr_peaks_jitter_amd{:}))/2,-1);

    similarity_sort = sort([pyr_peaks_jitter_gauss_similarity similarity_95],2,'ascend');    
    similarity_95 = similarity_sort(:,1:round((reps*0.05)));
    
    rep
    
end

% 3) Group significantly similar cells, record significance


% 4) Extract grouped cells
























