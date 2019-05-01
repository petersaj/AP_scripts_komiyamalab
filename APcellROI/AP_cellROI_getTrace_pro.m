
% Get image information
tiff_fullfilename = get(handles.tiffFile,'string');
imageinfo=imfinfo(tiff_fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% Create mask for selected ROI, get trace
if ~isfield(polygon,'autosort');
    cellMask = poly2mask(polygon.ROI{listselect}(:,1)',polygon.ROI{listselect}(:,2)',N,M);
    cellMask_1d = cellMask(:);
    cellTrace = mean(handles.im(cellMask_1d,:));
elseif isfield(polygon,'autosort');
    cellMask = polygon.autosort(:,:,listselect);
    cellMask_1d = cellMask(:);
    % careful here: matrix math requires at least single so this can be a
    % memory sink on a small computer
    cellTrace = cellMask_1d'*single(handles.im);
end



% norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
% dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
% x = min(cellTrace):dx:max(cellTrace);
% x = x(1:end-1); % it adds an extra, not sure why but can fix later
% norm_fit = pdf(norm_fit_obj,x');
% baseline_firing = x(norm_fit == max(norm_fit));

%estimate baseline with median
baseline_firing = median(cellTrace);

% change to df/f make peakin relative to df
%cellTrace = cellTrace./baseline_firing;

% NOT ESTIMATING DF/F FOR NOW!
    

