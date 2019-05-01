function Turboreg_AP_pro(targetfilename, sourcefilename, savefilename, channels)
%
% Call ImageJ class and Turboreg class to do motion correction on input
% image data file. Only works with translation of Turboreg. For other
% transformation algorithms, further modification is needed.
%
% squre_flag -- whether to scale the image to sqare before run Turboreg.
%
% Based on TurboregTK by TK.
%
% NX - May 09
%

% Load reference(tiff), otherwise assume im_t from gui
if ischar(targetfilename)
    im_t(:,:) = imread(targetfilename,'tiff');
else
    im_t = targetfilename;
end

target = im_t;
target_ij = array2ijStackAK(target);


% Load in data to be corrected(tiff), otherwise just assume im_s from gui
% if ischar(sourcefilename)
% [sourcepath, sourcename, sourceext] = fileparts(sourcefilename);
%    n_frame = length(imfinfo(sourcefilename,'tiff'))/channels;
%     for u = 1:channels:n_frame*channels
%         im_s(:,:,ceil(u/channels))=imread(sourcefilename,'tiff',u);
%     end
% else
%     im_s = sourcefilename;
%     n_frame = size(im_s,3);
% end

% Get information about source image
imageinfo=imfinfo(sourcefilename,'tiff');
if isfield(imageinfo(1),'ImageDescription')
    image_description = imageinfo(1).ImageDescription;
else
    disp('Warning! No scanimage header!')
end
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

%Do turboreg

clear im_s

disp('Done.')
clear im_registered
loadframes_channels = 1:channels:numframes;
% NOTE! ScanImage records uint/int depending on the version and can change
% without warning!! Load first frame, check type, initialize accordingly
frame_sample = imread(sourcefilename,'tiff',1,'Info',imageinfo);
integer_type = 'uint16'; % AP 140417: complicated tiff writing otherwise
%integer_type = class(frame_sample);
im_registered = zeros(N,M,length(loadframes_channels),integer_type);
disp('Turboreg correcting...')

%Cropping= ['0 0 ' num2str(size(target,2)) ' ' num2str(size(target,1))];
% Crop 10 frames on all sides: this makes a significant difference
crop_border = 150;
Cropping = [num2str(crop_border) ' ' num2str(crop_border) ' '...
    num2str(size(target,2)-crop_border) ' ' num2str(size(target,1)-crop_border)];
Transformation='-translation';
center_x = num2str(fix(size(target,2)/2)-1-crop_border);
center_y = num2str(fix(size(target,1)/2)-1-crop_border);
landmarks = [center_x,' ',center_y,' ',center_x,' ',center_y];
cmdstr=['-align -window s ',Cropping,' -window t ', Cropping,' ', ...
    Transformation, ' ', landmarks,' -hideOutput'];

al=IJAlign_AK;
for curr_frame = 1:length(loadframes_channels);
    ii = loadframes_channels(curr_frame);
    im_s = imread(sourcefilename,'tiff',ii,'Info',imageinfo);
    source = im_s;
    source_ij=array2ijStackAK(source);
    registered_ij = al.doAlign(cmdstr, source_ij, target_ij);
    registered = ij2arrayAK(registered_ij);
    a=eval([integer_type '(round(registered))']);
    im_registered(:,:,curr_frame) = a;
  
    disp(['Corrected frame (' savefilename '): ' num2str(curr_frame) '/' num2str(length(loadframes_channels))]);

    clear source_ij registered_ij
end

% find/save offsets
% calculate offsets based on rows/columns = 0 in one direction
clear zeros_y zeros_x turboreg_offsets
zeros_y = permute(any(im_registered,2),[1 3 2]);
zeros_y_offset = zeros_y-flipud(zeros_y);
turboreg_offsets.y = (sum(abs(zeros_y_offset)).*sign(zeros_y_offset(1,:)))/2;
zeros_x = permute(any(im_registered,1),[2 3 1]);
zeros_x_offset = zeros_x-flipud(zeros_x);
turboreg_offsets.x = (sum(abs(zeros_x_offset)).*sign(zeros_x_offset(1,:)))/2;
save([savefilename '_offsets.mat'],'turboreg_offsets')

% something about the java code might leak memory, so clear the shit out of
% everything related to java
java.lang.Runtime.getRuntime.gc % java garbage collector
clear source_ij registered_ij al mex

disp('Done. Saving...')
% You might get an error here about not being able to write to file
% because of permissions: SOLUTION IS TO NOT HAVE WINDOWS EXPLORER OPEN
% SIMULTANEOUSLY (but this windows lock solution might help fix)
% tic
% saveastiff(im_registered,[savefilename '.tif'],1,1,1,0,image_description);
% toc
for curr_frame = 1:length(loadframes_channels);
    if curr_frame == 1,
        for windows_lock = 1:100
            try
                if exist('image_description','var')
                    imwrite(im_registered(:,:,curr_frame),[savefilename '.tif'],'tif','Compression','none','WriteMode','overwrite', ...
                        'Description', imageinfo(curr_frame).ImageDescription);
                else
                    imwrite(im_registered(:,:,curr_frame),[savefilename '.tif'],'tif','Compression','none','WriteMode','overwrite');
                end
                break;
            catch me
                pause(0.2);
                continue
            end
            disp('Error! Didn''t Write')
            keyboard
        end
    else
        for windows_lock = 1:100
            try
                if exist('image_description','var')
                    imwrite(im_registered(:,:,curr_frame),[savefilename '.tif'],'tif','Compression','none','WriteMode','append', ...
                        'Description', imageinfo(curr_frame).ImageDescription);
                else
                    imwrite(im_registered(:,:,curr_frame),[savefilename '.tif'],'tif','Compression','none','WriteMode','append');
                end
                break;
            catch me
                pause(0.2);
                continue
            end
            disp('Error! Didn''t Write')
            keyboard
        end
    end;
    disp(['Frame written (' savefilename '): ' num2str(curr_frame) '/' num2str(length(loadframes_channels))]);
end

disp('Done.')