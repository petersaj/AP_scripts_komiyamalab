function [roi_trace_long roi_trace_long_split roi_trace_long_bg ...
    roi_trace_long_split_bg] = AP_getConcatTrace()
% [roi_trace_long roi_trace_long_split roi_trace_long_bg ...
%    roi_trace_long_split_bg] = AP_getConcatTrace()
% 
% roi_trace_long is the concatenated trace
% roi_trace_long split is the trace as a cell array, split up by file
% currently, zero values are treated as missing data points and ignored

% try
%     matlabpool
% catch me
% end

[roi_file roi_path] = uigetfile('*.roi','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
disp('Getting and concatenating traces');
if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end
roi_trace_all = [];
roi_trace_all_bg = [];
roi_trace_long_split = cell(length(tiff_filename),1);
roi_trace_long_split_bg = cell(length(tiff_filename),1);
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    % Create mask for ROIs
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
            %temp_mask = poly2mask(polygon.bgROI{n_polygon}(:,1)',polygon.bgROI{n_polygon}(:,2)',N,M);
            cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
        end
    end
    
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
    
    disp(['Getting activity from file ' num2str(i) '/' num2str(length(tiff_filename))]);
 
    roi_trace = zeros(size(cellMask,2),numframes);
    traceVal = zeros(size(cellMask,2),numframes);
    
    roi_trace_bg = zeros(size(cellMask,2),numframes);
    traceVal_bg = zeros(size(cellMask,2),numframes);
    if ~isfield(polygon,'autosort');
        
        % CURRENTLY: TREAT ZEROS AS MISSING DATA!
        % NOTE: sometimes for the b-scope zero points occur, but sparsely
        im(im == 0) = NaN;
        
        for n_polygon = 1:size(cellMask,2)
            % sum of pixel brighness / number of pixels, avg fluor in ROI
            traceVal(n_polygon,:) = nanmean(im(cellMask(:,n_polygon),:));
            traceVal_bg(n_polygon,:) = nanmean(im(cellMask_bg(:,n_polygon),:));
        end
        
    elseif isfield(polygon,'autosort');
        % weighted fluorescence
        traceVal = cellMask'*im;
    end
    roi_trace = traceVal;
    roi_trace_bg = traceVal_bg;
    
    % concat in case different number of frames
    roi_trace_long_split{i} = roi_trace;
    roi_trace_all = [roi_trace_all roi_trace];
    
    roi_trace_long_split_bg{i} = roi_trace_bg;
    roi_trace_all_bg = [roi_trace_all_bg roi_trace_bg];
    
    clear roi_trace roi_trace_bg
end

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = roi_trace_all;
roi_trace_long_bg = roi_trace_all_bg;
%roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);
% clear to save memory
im = [];