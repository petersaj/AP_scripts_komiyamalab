function Turboreg_AP_pro_2ch_parallel(targetfilename, sourcefilename, savefilename, channels, channel_correct)
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
spath = javaclasspath('-static');
spath = cell2mat(spath(1));
javafolder = strfind(spath,['java' filesep]);
javafolder = spath(1:javafolder+3);

javaaddpath([javafolder filesep 'jar' filesep 'ij.jar'])
javaaddpath([javafolder filesep 'jar' filesep 'mij.jar'])
javaaddpath([javafolder])

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
loadframes_channels = channel_correct:channels:numframes;
% NOTE! ScanImage records uint/int depending on the version and can change
% without warning!! Load first frame, check type, initialize accordingly
frame_sample = imread(sourcefilename,'tiff',1,'Info',imageinfo);
integer_type = 'uint16'; % AP 140417: complicated tiff writing otherwise
%integer_type = class(frame_sample);

red_registered = zeros(N,M,length(loadframes_channels),integer_type);
disp('Turboreg correcting...')

%Cropping= ['0 0 ' num2str(size(target,2)) ' ' num2str(size(target,1))];
% Crop 10 frames on all sides: this makes a significant difference
crop_border = 10;
Cropping = [num2str(crop_border) ' ' num2str(crop_border) ' '...
    num2str(size(target,2)-crop_border) ' ' num2str(size(target,1)-crop_border)];
Transformation='-translation';
center_x = num2str(fix(size(txarget,2)/2)-1-crop_border);
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
    red_registered(:,:,curr_frame) = a;
  
%     disp(['Corrected frame (' savefilename '): ' num2str(curr_frame) '/' num2str(length(loadframes_channels))]);

    clear source_ij registered_ij
end

% find/save offsets
% calculate offsets based on rows/columns = 0 in one direction
clear zeros_y zeros_x turboreg_offsets
zeros_y = permute(any(red_registered,2),[1 3 2]);
zeros_y_offset = zeros_y-flipud(zeros_y);
turboreg_offsets.y = (sum(abs(zeros_y_offset)).*sign(zeros_y_offset(1,:)))/2;
zeros_x = permute(any(red_registered,1),[2 3 1]);
zeros_x_offset = zeros_x-flipud(zeros_x);
turboreg_offsets.x = (sum(abs(zeros_x_offset)).*sign(zeros_x_offset(1,:)))/2;
save([savefilename '_offsets.mat'],'turboreg_offsets')

% load in green channel
disp('Loading in green channel, applying offsets from red...')
greenframes_channels = setdiff(1:numframes,loadframes_channels);
green_raw = zeros(size(red_registered));
for curr_frame = 1:length(greenframes_channels);
    ii = greenframes_channels(curr_frame);
    green_raw(:,:,curr_frame) = imread(sourcefilename,'tiff',ii,'Info',imageinfo);
end

% correct green by red offsets
green_registered = zeros(size(green_raw),'uint16');
for i = 1:size(green_raw,3)
    green_registered(zeros_y(:,i),zeros_x(:,i),i) = ...
        green_raw(flipud(zeros_y(:,i)),flipud(zeros_x(:,i)),i);
end

% error out if different size red and green channels
if any(size(red_registered) ~= size(green_registered))
    error('Red and green registered images not the same size')
end

% interlace frames for saving
im_registered = zeros(size(red_registered,1),size(red_registered,2), ...
    size(red_registered,3)*2,'uint16');
im_registered(:,:,greenframes_channels) = green_registered;
im_registered(:,:,loadframes_channels) = red_registered;

% something about the java code might leak memory, so clear the shit out of
% everything related to java
java.lang.Runtime.getRuntime.gc % java garbage collector
clear source_ij registered_ij al
    
disp('Done. Saving...')
% You might get an error here about not being able to write to file
% because of permissions: SOLUTION IS TO NOT HAVE WINDOWS EXPLORER OPEN
% SIMULTANEOUSLY (but this windows lock solution might help fix)
% tic
% saveastiff(im_registered,[savefilename '.tif'],1,1,1,0,image_description);
% toc
if isunix
    SlashInds = strfind(sourcefilename,'/');
    ShortFileName = sourcefilename((SlashInds(end)+1):end);
    tempfilename = ['/dev/shm/' ShortFileName];
else
    SlashInds = strfind(sourcefilename,'\');
    ShortFileName = sourcefilename((SlashInds(end)+1):end);
    tempfilename = ['G:\' ShortFileName];
end 


for curr_frame = 1:size(im_registered,3);
    if curr_frame == 1,
        writemode = 'overwrite';
    else
        writemode = 'append';
    end
    for windows_lock = 1:100
        try
%             if exist('image_description','var')
                imwrite(im_registered(:,:,curr_frame),tempfilename,'tif','Compression','none','WriteMode',writemode, ...
                    'Description',imageinfo(curr_frame).ImageDescription);
%             else
%                 imwrite(im_registered(:,:,curr_frame),tempfilename,'tif','Compression','none','WriteMode',writemode);
%             end
            break;
        catch me
            pause(0.2);
            continue
        end
        disp('Error! Didn''t Write')
        keyboard
    end    
%     disp(['Frame written (' savefilename '): ' num2str(curr_frame) '/' num2str(size(im_registered,3))]);
end
movefile(tempfilename,[savefilename '.tif'],'f')


disp('Done.')