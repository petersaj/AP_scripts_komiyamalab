% Make AVI movie with imaging on top and lever trace on bottom

%% load in image file

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','off');

img_filename = [tiff_path tiff_filename];

% get traces from current tiff file
imageinfo=imfinfo(img_filename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% try
%     matlabpool
% catch me
% end

im = zeros(N,M,numframes,'uint16');
disp('Loading file....')
for loadframe = 1:numframes
    im(:,:,loadframe) = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
    disp(['Frame ' num2str(loadframe) '/' num2str(numframes)]);
end
disp('Done')

%% make a movie with cells on top and trace on bottom
use_lever = bhv.lever_force_resample(3000:end);
use_frames = round((length(use_lever)/length(bhv.lever_force_resample)) * ...
    size(im,3));
im_use = im(:,:,end-use_frames:end);


%f = figure;
s1 = subplot(3,4,[2 3 6 7]); colormap(gray)
s2 = subplot(3,1,3);

plot(-use_lever,'k','linewidth',2)
xlim([0 length(use_lever)+1])


winsize = get(f,'Position');

force = length(use_lever);
force_inc = force/size(im_use,3);

subplot(s2)
ylim([-2 -1])
tempy = [-1.95 -1.05];

curr = patch([(i-1)*force_inc i*force_inc ...
       i*force_inc (i-1)*force_inc], ...
       [tempy(1) tempy(1) tempy(2) tempy(2)],...
       'g');
set(gca,'Children',flipud(get(gca,'Children')));
set(gca,'Visible','off')
   
for i = 1:use_frames
   subplot(s1)
   imagesc(im_use(:,:,i));
   caxis([0 1000])
   set(gca,'Visible','off')
   
   subplot(s2)
   set(curr,'XData', ...
       [(i-1)*force_inc+2 i*force_inc ...
       i*force_inc (i-1)*force_inc]);
   
   A(i) = getframe(f);
end

%% Save movie
filename = 'AP75_day1.avi';
writerObj = VideoWriter(filename);
writerObj.FrameRate = 20;
open(writerObj);
writeVideo(writerObj,A);
close(writerObj);






%% 
fs1 = 10; % Original sampling frequency in Hz
t1 = 0:0.01:2*pi; % Time vector
x = cos(t1); % Define a linear sequence
xpad = [repmat(x(1), 1, 100), x, repmat(x(end), 1, 100)];
tpad = [-1/fs1*10 : 1/fs1: 0-1/fs1, t1, 1+1/fs1:1/fs1:1+1/fs1*10];
ypad = resample(xpad,2,3); % Now resample it
t2 = (0:(length(ypad)-1))*2/(3*fs1) - 1; % New time vector
plot(t1,x,'*',t2,ypad,'o',(-0.5:0.01:1.5),(-0.5:0.01:1.5)+100,':')
legend('original','resampled'); xlabel('Time')



