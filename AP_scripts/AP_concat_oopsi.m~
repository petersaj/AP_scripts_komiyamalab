function AP_concat_oopsi(animal)

% set oopsi parameters
oopsi_param.tau = 1;
oopsi_param.pyr_lambda = 0.1;
oopsi_param.gad_lambda = 0.3;

% find days
data_path = ['/usr/local/lab/People/Andy/Data/' animal];
data_dir = dir(data_path);
data_names = {data_dir.name};
data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
    'UniformOutput',false);
data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
data_folders = data_fullnames(data_isfolders);
data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);

% load labels file
label_path = ['/usr/local/lab/People/Andy/Data/' animal filesep ...
    animal '_roi_template'];
labels_name = [animal '_roilabels.roilabel'];
load([label_path filesep labels_name],'-MAT');

% get cell labels
cells = 1:length(roi_labels);
gad_cells = find(cellfun(@(x) any(strcmp('gad',x)),roi_labels))';
pyr_cells = cells(~ismember(cells,gad_cells));

% initialize results cell arrays
concat_oopsi = cell(length(cells),length(days));
master_roi_trace_df = cell(length(cells),length(days));
master_index = cell(length(cells),length(days));
master_daylength = nan(size(days));

concat_days = 1:length(days);
% AP71 day 1+2 bad
if strcmp(animal,'AP71')
    concat_days = 3:length(days);
end

% Get full concatenated trace over days to do deconvolution
for curr_day = concat_days
    clearvars -except days animal oopsi_param ...
        pyr_cells gad_cells concat_oopsi ...
        curr_day concat_days master_roi_trace_df master_index master_daylength ...
        roi_df_oopsi
    
    day = num2str(days{curr_day});
    
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_name = [day '_' animal '_bgROI_analysis.mat'];
    load([analysis_path filesep analysis_name]);
    
    inactive_thresh = 0; % threshold for stitching
    num_loops = length(im.roi_trace_split)/8;
    file_lengths = cellfun(@length,im.roi_trace_split);
    file_lengths_loops = reshape(file_lengths,8,num_loops);
    loop_lengths = sum(file_lengths_loops);
    
    % Break up baseline corrected trace into loops
    roi_trace_df_split = mat2cell(im.roi_trace_df,[size(im.roi_trace_df,1)], ...
        loop_lengths);
    
    split_indx = mat2cell(1:size(im.roi_trace_df,2),1, ...
        loop_lengths);
    
    curr_loopcut_df = cell(size(im.roi_trace_df,1),1);
    curr_loopcut_stitch = cell(size(im.roi_trace_df,1));
    for curr_cell = 1:size(im.roi_trace_df,1)
        
        % Stitch together loops at inactive parts
        % Skip the first 10 frames in case of weird things
        loop_start_stitches = cellfun(@(x) find(x(curr_cell,11:end) ...
            < inactive_thresh,1)+10, ...
            roi_trace_df_split,'UniformOutput',false);
        loop_stop_stitches = cellfun(@(x) find(x(curr_cell,:) ...
            < inactive_thresh,1,'last'), ...
            roi_trace_df_split,'UniformOutput',false);
        F_stitched = cellfun(@(x,y,z) x(curr_cell,y:z),roi_trace_df_split,...
            loop_start_stitches,loop_stop_stitches,'UniformOutput',false);
        stitching_indx = cellfun(@(x,y,z) x(y:z), split_indx, ...
            loop_start_stitches, loop_stop_stitches,'UniformOutput',false);
        
        F = [F_stitched{:}];
        
        if isempty(F)
            continue
        end
        
        master_roi_trace_df{curr_cell,curr_day} = F;
        master_index{curr_cell,curr_day} = [stitching_indx{:}];
    end
    
    % get the total number of frames for each day (to index back in)
    master_daylength(curr_day) = size(im.roi_trace_df,2);
    
    disp(['Loaded ' day])
    
end

roi_df_oopsi = cell(size(master_roi_trace_df));

tic
disp('Deconvolving...')
fprintf('%2d',0)
for curr_cell = 1:size(im.roi_trace_df,1)
    
    % choose appropriate lambda
    if ismember(curr_cell,pyr_cells)
        lambda = oopsi_param.pyr_lambda;
    elseif ismember(curr_cell,gad_cells)
        lambda = oopsi_param.gad_lambda;
    end
    
    F = cell2mat(master_roi_trace_df(curr_cell,:));
    
    F(isnan(F)) = 0;
    
    % set options/parameters
    clear V
    V.dt = 1/bhv.framerate;
    V.est_sig = 1;
    V.est_lam = 1;
    V.est_gam = 1;
    V.est_b = 1;
    V.est_a = 1;
    
    clear P
    P.gam = 1-(V.dt/oopsi_param.tau);
    P.lam = lambda;
    
    % do oopsi deconvolution
    [n_out P_out V_out] = fast_oopsi(F,V,P);
    norm_oopsi = n_out/max(n_out);
    
    % split deconvolution back into days
    oopsi_lengths = cellfun(@length,master_roi_trace_df(curr_cell,:));
    norm_oopsi_split = mat2cell(norm_oopsi',1,oopsi_lengths);
    
    for curr_day = concat_days
        roi_df_oopsi{curr_cell,curr_day} = zeros(1,master_daylength(curr_day));
        roi_df_oopsi{curr_cell,curr_day}(master_index{curr_cell,curr_day}) = ...
            norm_oopsi_split{curr_day};
    end
    
    fprintf('%c%c%2d',8,8,floor(100*curr_cell/size(im.roi_trace_df,1)));
end

disp('Finished deconvolution')
toc

clearvars -except animal oopsi_param roi_df_oopsi days concat days

%%%% Save as a whole in a seperate file
disp('Saving as whole...')
savename = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI/' animal '_concat_oopsi.mat'];
save(savename,'oopsi_param','roi_df_oopsi','-v7.3')

%%%% Save within each analysis file seperately
disp('Saving within analysis files...')
for curr_day = concat_days
    clear im bhv
    day = days{curr_day};
    analysis_path = ['/usr/local/lab/People/Andy/Data/' animal filesep animal '_bgROI'];
    analysis_name = [day '_' animal '_bgROI_analysis.mat'];
    analysis_fullname = [analysis_path filesep analysis_name];
    load(analysis_fullname);
    
    im.roi_concat_oopsi = cell2mat(roi_df_oopsi(:,curr_day));
    
    save(analysis_fullname,'bhv','im');
end

disp(['Finished concat oopsi: ' animal]);















