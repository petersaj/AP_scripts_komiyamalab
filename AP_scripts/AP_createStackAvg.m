function AP_createStackAvg();
% Create an average image out of selected stacks

[filelist,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','on');

if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end

% Create parsed reference

savedir = tiff_path;

% get image size from first file
imageinfo=imfinfo(filelist{1},'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

disp('Creating average of stacks')
ref_numframes = 0;
im_t = zeros(N,M);

for pfiles = 1:length(filelist);
    ref_iminfoCurrent=imfinfo(filelist{pfiles},'tiff');
    ref_numframesCurrent=length(ref_iminfoCurrent);
    for ref_frame = 1:ref_numframesCurrent
        im_t(:,:)=im_t(:,:)+double(imread(filelist{pfiles},'tiff',ref_frame))./ref_numframes;
    end
    disp(['Creating average image: ' ...
        num2str(100*(pfiles/length(filelist))) '%']);
end

% save temporary combined average
im_t = uint16(im_t);
for windows_lock = 1:100
    try
        imwrite(im_t,[savedir filesep filelist{1} '_AVG.tif'],'tif','compression','none','writemode','overwrite');
        break;
    catch me
        pause(0.2);
    end
end
