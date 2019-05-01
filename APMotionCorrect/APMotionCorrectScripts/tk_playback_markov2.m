function [fixeddata,countdata]=tk_playback_markov(fullfilename,savefile,offsets,edgebuffer,playback,lambda,numchannels)
%fixeddata,countdata]=playback_markov(imagedata,offsets,edgebuffer,playback)

%pick out the number of frames and size of images
% numframes=size(imagedata,1);
imageinfo=imfinfo(fullfilename,'tiff');
[filepath, filename] = fileparts(fullfilename);
numframes=length(imageinfo);
if ~exist('numchannels')
    %automated way of determining whether this is a one channel or two channel
    %movie... look at the std of frame 1 vs 2 and frame 1 vs 3.
    first=double(imread(fullfilename,'tiff',1));
    second=double(imread(fullfilename,'tiff',2));
    third=double(imread(fullfilename,'tiff',3));
    junka=first-second;
    junkb=first-third;
    if (std(junka(:))>1.25*std(junkb(:)))
        numchannels=2;
    else
        numchannels=1;
    end
end
if numchannels==2
        numframes=numframes/2;
end    
N=imageinfo(1).Height;
M=imageinfo(1).Width;
%load the first three images just for the code
imagedata=zeros(3,N,M);
for i = 1:3
    if numchannels==2
        imagedata(i,:,:)=imread(fullfilename,'tiff',2*i-1);
    else
        imagedata(i,:,:)=imread(fullfilename,'tiff',i);
    end
end

%dynamically set the dynamic range for playback by taking the 99.99% pixel
%removes outliers. can adjust downward for automated gain control
sortpix=imagedata(1:3,:,:);
sortpix=sortpix(:);
sortpix=sort(sortpix);
k=round(length(sortpix)*.9999);
maxpixel=sortpix(k)
minpixel=sortpix(1)


%initialize the data 
countdata=int8(zeros(N,M));
fixeddata=zeros(N,M);

%in case you just wanted to do a few frames
frames=1:numframes;

%index for offsets
m=0;

%for saving fixed image
h = waitbar(0,'');
imagedata=zeros(1,N,M);
for i=frames
    if numchannels==2
        imagedata(1,:,:)=imread(fullfilename,'tiff',2*i-1);
    else
        imagedata(1,:,:)=imread(fullfilename,'tiff',i);
    end
    %initialize a matrix just for this frame for simplicity
    correctimage=zeros(N,M);
    countimage=zeros(N,M);
    
    %loop over the lines we are considering placing
    for j=1+edgebuffer:N-edgebuffer
        %increment the offset counter
        m=m+1;
        
        %pick out the line of the data we are placing
        lineextract=double(squeeze(imagedata(1,j,:)));
 
        %pick out the line number where it is going
        linenumber=offsets(1,m)+j;
        %pick out the relative shift in X within that line
        shift=offsets(2,m);
   
        %need different bounds for shifts left and shifts right
        %but in general increments correctimage and countimage correctly
        if shift<0    
            correctimage(linenumber,1:end+shift)=correctimage(linenumber,1:end+shift)+lineextract(1-shift:end)';
            countimage(linenumber,1:end+shift)=countimage(linenumber,1:end+shift)+ones(1,M+shift);
        
        else
            correctimage(linenumber,shift+1:end)=correctimage(linenumber,shift+1:end)+lineextract(1:end-shift)';
            countimage(linenumber,shift+1:end)=countimage(linenumber,shift+1:end)+ones(1,M-shift);
        end
      
    end
    %show the framenumber you have just completed
    waitbar(i/frames(end),h,'Correcting and saving images');
    %store the results into the data structure
    fixeddata(:,:)=correctimage;
    countdata(:,:)=countimage;
    im = uint16(correctimage./(countimage));
    for windows_lock = 1:100
        try
            imwrite(im,[savefile '.tif'],'tif','Compression','none','WriteMode','append', ...
                'Description', imageinfo(1).ImageDescription);
            break;
        catch me
            pause(0.2);
        end
    end
    %if you want to visualize the results do so
    if playback

        
        figure(3);
        clf;
        set(gcf,'Position',[50 314 560 420]);
        %finds those pixels which were not sampled and sets their counts to
        %infinity for display purposes
        badones=find(countimage==0);
        countimage(badones)=Inf;
        
        %display, normalizing for multiple samples
        imagesc(correctimage./(countimage));
        colormap(gray);
        %set the dynamic range
        caxis([minpixel maxpixel]);
        %title it by frame number
        title(num2str(i));
        
%         figure(4);
%         clf;
%         set(gcf,'Position',[50+560 314 560 420]);
%         
%         imagesc(squeeze(imagedata(i,:,:)));
%         colormap gray;
%         caxis([minpixel maxpixel]);
%         title(num2str(i));
        pause(.1);
        
    end

end
close(h);



