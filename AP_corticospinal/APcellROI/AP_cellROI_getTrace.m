
% Get image information
tiff_fullfilename = get(handles.tiffFile,'string');
imageinfo=imfinfo(tiff_fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% Create mask for selected ROI, get trace
if ~isfield(polygon,'autosort');
    cellMask = poly2mask(polygon.ROI{listselect}(:,1)',polygon.ROI{listselect}(:,2)',N,M);
elseif isfield(polygon,'autosort');
    cellMask = polygon.autosort(:,:,listselect);
end

w = waitbar(0,'Creating ROI Trace');
for frame=1:numframes
    frame_reshape = reshape(handles.im(:,frame),N,M);
    normactivityMap = (double(frame_reshape).*cellMask)./sum(cellMask(:));
    traceVal = sum(normactivityMap(:));
    cellTrace(frame) = traceVal;
    waitbar(frame/numframes,w,'Creating ROI Trace'); % not with parfor
end
waitbar(100,w,'Normalizing Trace')
norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
x = min(cellTrace):dx:max(cellTrace);
x = x(1:end-1); % it adds an extra, not sure why but can fix later
norm_fit = pdf(norm_fit_obj,x');
baseline_firing = x(norm_fit == max(norm_fit));

% change to df/f make peakin relative to df
cellTrace = cellTrace./baseline_firing;
close(w);

    

