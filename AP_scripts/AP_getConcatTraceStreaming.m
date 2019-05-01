function roi_trace_long = AP_getConcatTraceStreaming()

try
    matlabpool
catch me
end

[roi_file roi_path] = uigetfile('*.roi','ROI file');
load([roi_path roi_file],'-MAT')
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
w = waitbar(0,'Getting and concatenating traces');
if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end
roi_trace_all = [];
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    
    % get info from current tiff file
    tsStack = TIFFStack(img_filename);
    numframes = size(tsStack,3);
    M = size(tsStack,2);
    N = size(tsStack,1);
    
    autosort_flag = 0;
    
    % Create mask for selected ROI, get trace
    for n_polygon = 1:length(polygon.ROI);
        if ~isfield(polygon,'autosort');
            temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
            cellMask(:,n_polygon) = temp_mask(:);
            autosort_flag = 0;
        elseif isfield(polygon,'autosort');
            temp_mask = polygon.autosort(:,:,n_polygon);
            cellMask(:,n_polygon) = temp_mask(:);
            autosort_flag = 1;
        end
    end
    
    disp('Getting trace....')
    roi_trace = [];
    roi_trace = zeros(size(cellMask,2),numframes);
    parfor loadframe = 1:numframes
        im = [];        
        im = tsStack(:,:,loadframe);
        % image: columns are frames
        im = im(:);
        
        traceVal = [];
        if autosort_flag == 0;
            for n_polygon = 1:size(cellMask,2)
                % sum of pixel brighness / number of pixels, avg fluor in ROI
                traceVal(n_polygon,:) = sum(im(cellMask(:,n_polygon),:)./sum(sum(cellMask(:,n_polygon))));
            end
        elseif autosort_flag == 1;
            for n_polygon = 1:size(cellMask,2)
                % sum of pixel brighness / number of pixels, avg fluor in ROI
                traceVal(n_polygon,:) = cellMask(cellMask(:,n_polygon)~=0,n_polygon)'*im(cellMask(:,n_polygon)~=0,:)...
                    ./sum(sum(cellMask(:,n_polygon)~=0));
            end
        end
        roi_trace(:,loadframe) = traceVal;
        disp(num2str(loadframe/numframes));
    end
    disp('Done');

    % concat in case different number of frames
    roi_trace_all = [roi_trace_all roi_trace];
    clear roi_trace
    waitbar(i/length(tiff_filename),w,'Getting and concatenating traces');
    
end
close(w);

% reshape traces from stacks to long form
num_roi = size(roi_trace_all,1);
num_frames = size(roi_trace_all,2);
num_img = size(roi_trace_all,3);
roi_trace_long = roi_trace_all;
%roi_trace_long = reshape(roi_trace_all,num_roi,num_frames*num_img);
% clear to save memory
im = [];