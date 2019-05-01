% Create parsed reference

% savedir = filelist{1}{1,1};

% Use second-to-last file
ref_filename = filelist{end-1};

% get image size from first file
imageinfo=imfinfo(ref_filename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

disp('Creating Reference: Using average of penultimate file')
%ref_numframes = 0;
im_t = zeros(N,M);

% % get total number of frames from all files
% for pfiles2 = 1:length(filelist);
%     ref_iminfo=imfinfo(filelist{pfiles2},'tiff');
%     ref_numframes_long=length(ref_iminfo);
%     ref_numframes=ref_numframes + ref_numframes_long/channels;
% end


% New reference: find most stable 100 frames of last movie since most
% likely to not have movement

ref_numframes = length(imageinfo);
im_ref = zeros(N*M,ref_numframes);
for ref_frame = 1:ref_numframes
   curr_frame = double(imread(ref_filename,'tiff',ref_frame));
   im_ref(:,ref_frame) = curr_frame(:);
end

im_diff = abs(diff(im_ref,1,2));
im_diff = im_diff(:,50:end-50);
[im_diff_min im_diff_min_indx] = min(smooth(sum(im_diff),100));
im_diff_min_indx = im_diff_min_indx + 50;
%im_t = reshape(mean(im_ref(:,im_diff_min_indx-49:im_diff_min_indx+49),2), ...
%    [512 512]);
im_t = reshape(mean(im_ref,2), ...
    [512 512]);

if all(isnan(im_t(:)))
    error('Problem creating reference image');
end
% 
% im_diff = zeros(ref_numframes-1,1);
% curr_diff = 1;
% for pfiles = length(filelist);
%     ref_iminfoCurrent=imfinfo(filelist{pfiles},'tiff');
%     ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
%     for ref_frame = 1:channels:ref_numframesCurrent*channels      
%         curr_im = double(imread(filelist{pfiles},'tiff',ref_frame));     
%         im_t(:,:) = im_t(:,:) + curr_im;
%         im_diff(curr_diff) = sum(abs(curr_im(:) - last_im(:)));
%         curr_diff = curr_diff + 1;
%         last_im = curr_im;
%     end
%     disp(['Creating Reference from Parsed: ' ...
%         num2str(100*(pfiles/length(filelist))) '%']);
%  end

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

