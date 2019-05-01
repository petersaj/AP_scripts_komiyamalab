function AP_getConcatTrace_batch(animal,roi_path)
% AP_getConcatTrace_batch(animal,roi_path)

disp('CURRENTLY SPECIFIC TO ANDY')

% get roi files and corresponding days
roi_dir = dir(roi_path);
roi_dir_filenames = {roi_dir.name};
roi_files_indx = cellfun(@(x) ~isempty(strfind(x,'.roi')), roi_dir_filenames);
roi_filenames = roi_dir_filenames(roi_files_indx);
roi_filenames = sort(roi_filenames);

days = cellfun(@(x) x(1:6), roi_filenames,'UniformOutput',0);

for curr_day = 1:length(days);
    
    tiff_path = ['/usr/local/lab/People/Andy/Data/' ...
        animal filesep days{curr_day}];
    
    dir_currfolder = dir(tiff_path);
    im_dir_filenames = {dir_currfolder.name};
    tiff_files_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), im_dir_filenames);
    tiff_filenames = im_dir_filenames(tiff_files_indx);
    
    % find only tiff files have been turboreg'd
    tiff_corrected_filenames_indx = find(cellfun(@(x) strcmp(x(end-3:end),'.tif'),tiff_filenames));
    tiff_corrected_filenames = tiff_filenames(tiff_corrected_filenames_indx);
    tiff_corrected_filenames = sort(tiff_corrected_filenames);
    % display them - to know they're in order
    disp(vertcat(tiff_corrected_filenames{:}));
    
    % append directory to tiff filename
    tiff_corrected_filenames = cellfun(@(x) [tiff_path filesep x], ...
        tiff_corrected_filenames,'UniformOutput',false);
    
    % load ROIs
    load([roi_path filesep roi_filenames{curr_day}],'-MAT')

    disp('Getting and concatenating traces');
    % convert tiff names to cell, in case there's only one file
    if ~iscell(tiff_corrected_filenames);
        tiff_corrected_filenames = {tiff_corrected_filenames};
    end
    
    % get traces from current tiff file
    imageinfo=imfinfo(tiff_corrected_filenames{1},'tiff');
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % create masks from ROIs
    cellMask = zeros(N*M,length(polygon.ROI));
    for n_polygon = 1:length(polygon.ROI);
        if ~isfield(polygon,'autosort');
            temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            cellMask(:,n_polygon) = temp_mask(:);
        elseif isfield(polygon,'autosort');
            temp_mask = polygon.autosort(:,:,n_polygon);
            cellMask(:,n_polygon) = temp_mask(:);
        end
    end
    
    % Create mask for background ROIs, if available
    if isfield(polygon,'bgROI')
        for n_polygon = 1:length(polygon.bgROI);
            cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
        end
    end
    
    % get the traces for the ROIs
    roi_trace_long_split = cell(length(tiff_corrected_filenames),1);
    roi_trace_long_split_bg = cell(length(tiff_corrected_filenames),1);
    
    for i = 1:length(tiff_corrected_filenames);
        img_filename = [tiff_corrected_filenames{i}];
        
        % get traces from current tiff file
        imageinfo=imfinfo(img_filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        im = [];
        im = zeros(N*M,numframes);
        
        tic
        disp('Loading file....')
        for loadframe = 1:numframes
            im_temp = [];
            im_temp = imread(img_filename,'tiff',loadframe);
            % image: columns are frames
            im(:,loadframe) = im_temp(:);
        end
        disp('Done')
        toc
   
        disp(['Getting activity from file ' num2str(i) '/' num2str(length(tiff_corrected_filenames))]);
        
        roi_trace = zeros(size(cellMask,2),numframes);
        traceVal = zeros(size(cellMask,2),numframes);
        
        roi_trace_bg = zeros(size(cellMask,2),numframes);
        traceVal_bg = zeros(size(cellMask,2),numframes);
        
        if ~isfield(polygon,'autosort');
            
            % treat zeros as missing values
            im(im == 0) = NaN;
            
            for n_polygon = 1:size(cellMask,2)
                % sum of pixel brighness / number of pixels, avg fluor in ROI
                traceVal(n_polygon,:) = nanmean(im(cellMask(:,n_polygon),:));
                traceVal_bg(n_polygon,:) = nanmean(im(cellMask_bg(:,n_polygon),:));
            end
            
        elseif isfield(polygon,'autosort');
            % weighted fluorescence
            traceVal = (cellMask'*im)./repmat(sum(cellMask',2),1,size(im,2));
        end
        roi_trace = traceVal;
        roi_trace_bg = traceVal_bg;
        
        % save current file's traces
        roi_trace_long_split{i} = roi_trace;
        roi_trace_long_split_bg{i} = roi_trace_bg;
    
        clear roi_trace roi_trace_bg
    end
    
    % concatenate traces into one continuous trace for the day
    roi_trace_long = [];
    roi_trace_long = horzcat(roi_trace_long_split{:});
    
    roi_trace_long_bg = [];
    roi_trace_long_bg = horzcat(roi_trace_long_split_bg{:});
    
    % clear to save memory
    im = [];
    
    savename = [];
    savename = [roi_path filesep roi_filenames{curr_day}(1:end-4) ...
        '_analysis.mat'];
    % save roi_trace_long and roi_trace_long_split in the ROI path
    if isfield(polygon,'bgROI')
        save(savename,'roi_trace_long','roi_trace_long_split', ...
            'roi_trace_long_bg','roi_trace_long_split_bg');
    else
        save(savename,'roi_trace_long','roi_trace_long_split');
    end
    
end