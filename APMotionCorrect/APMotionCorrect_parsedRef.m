% Create parsed reference

savedir = filelist{1}{1,1};

% get image size from first file
imageinfo=imfinfo(filelist{1}{2,1},'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

disp('Creating Reference from Parsed')
ref_numframes = 0;
im_t = zeros(N,M);
% get total number of frames from all files
for pfiles2 = 1:size(filelist,1);
    ref_iminfo=imfinfo(filelist{pfiles2}{2,1},'tiff');
    ref_numframes_long=length(ref_iminfo);
    ref_numframes=ref_numframes + ref_numframes_long/channels;
end
for pfiles = 1:size(filelist,1);
    ref_iminfoCurrent=imfinfo(filelist{pfiles}{2,1},'tiff');
    ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
    for ref_frame = 1:channels:ref_numframesCurrent*channels
        im_t(:,:)=im_t(:,:)+double(imread(filelist{pfiles}{2,1},'tiff',ref_frame))./ref_numframes;
    end
    disp(['Creating Reference from Parsed: ' ...
        num2str(100*(pfiles/size(filelist,1))) '%']);
end
% % save temporary combined average
% im_t = uint16(im_t);
% for windows_lock = 1:100
%     try
%         imwrite(im_t,[savedir filesep 'parsedAverage.tif'],'tif','Compression','none','WriteMode','overwrite');
%         break;
%     catch me
%         pause(0.2);
%     end
% end
% clear im_t

