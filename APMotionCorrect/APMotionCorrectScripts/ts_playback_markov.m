function [fixeddata,countdata]=playback_markov(maxdx, maxdy, lambdax, lambday, gain, fullfilename, correctimagefilename, combinedimagefilename, numchannels, ImageDescription, ...
    frameanalysismode, imagedata, offsets, edgebuffer)
%fixeddata,countdata]=playback_markov(imagedata,offsets,edgebuffer,playback,totvel)
%        playback_markov(fullfilename, correctimagefilename, numchannels, 1, chone, offsets,edgebuffer,1);


%pick out the number of frames and size of images
playback = 1;
imageinfo = imfinfo(fullfilename, 'tiff');
numframes = length(imageinfo)/numchannels;
N=imageinfo(1).Height;%size(imagedata,2);
M=imageinfo(1).Width;%size(imagedata,3);


%add info on ImageDescription / change numchannels
ImageDescription = [ImageDescription, sprintf('state.analysis.maxdx=%d', maxdx), 13, ...
    sprintf('state.analysis.maxdy=%d', maxdy), 13, ...
    sprintf('state.analysis.lambdax=%2.5f', lambdax), 13, ...
    sprintf('state.analysis.lambday=%2.5f', lambday), 13, ...
    sprintf('state.analysis.gain=%2.5f', gain), 13];

%change numberOfChannelsAcquire
BP = find(double(ImageDescription)==13);%Return
NumOfLine = length(BP)+1;
BP = [0 BP length(ImageDescription)+1];
for line=1:NumOfLine
    Comment{line} = ImageDescription(BP(line)+1 : BP(line+1)-1);
    if length(Comment{line})>34 && isequal(Comment{line}(1:34), 'state.acq.numberOfChannelsAcquire=')
        ImageDescription = [ImageDescription(1:BP(line)), 'state.acq.numberOfChannelsAcquire=1', ...
            ImageDescription(BP(line+1):end)];
    end
end

%change numberOfChannelsSave
BP = find(double(ImageDescription)==13);%Return
NumOfLine = length(BP)+1;
BP = [0 BP length(ImageDescription)+1];
for line=1:NumOfLine
    Comment{line} = ImageDescription(BP(line)+1 : BP(line+1)-1);
    if length(Comment{line})>31 && isequal(Comment{line}(1:31), 'state.acq.numberOfChannelsSave=')
        ImageDescription = [ImageDescription(1:BP(line)), 'state.acq.numberOfChannelsSave=1', ...
            ImageDescription(BP(line+1):end)];
    end
end




%dynamically set the dynamic range for playback by taking the 99.99% pixel
%removes outliers. can adjust downward for automated gain control
sortpix=imread(fullfilename,'tiff',1);
sortpix=sortpix(:);
sortpix=sort(sortpix);
k=round(length(sortpix)*.9999);
maxpixel=sortpix(k)
minpixel=sortpix(1)


%initialize the data 
% countdata=int8(zeros(numframes,N,M));
% fixeddata=zeros(numframes,N,M);

%in case you just wanted to do a few frames
frames=1:numframes;

%index for offsets
m=0;


for i=frames
    %initialize a matrix just for this frame for simplicity
    correctimage=zeros(N,M);
    countimage=zeros(N,M);
    %loop over the lines we are considering placing
    for j=1+edgebuffer:N-edgebuffer
        %increment the offset counter
        m=m+1;

        %pick out the line of the data we are placing
        if frameanalysismode
            if numchannels==2
                originalimage = imread(fullfilename,'tiff',2*i-1);
            else
                originalimage = imread(fullfilename,'tiff',i);
            end
            lineextract = double(squeeze(originalimage(j,:)))';
        else
            lineextract = squeeze(imagedata(i,j,:));
        end

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
    %disp(i);
    %store the results into the data structure
%     fixeddata(i,:,:)=correctimage;
%     countdata(i,:,:)=countimage;

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
        %pause(.1);
    end

    im = uint16(round(correctimage./(countimage)));
    combinedimage = [originalimage zeros(size(originalimage,1), 5) im];
    
    if i==min(frames)
        imwrite(im,correctimagefilename,'tif','Compression','none','WriteMode','overwrite', 'Description', ImageDescription);
        imwrite(combinedimage,combinedimagefilename,'tif','Compression','none','WriteMode','overwrite', 'Description', ImageDescription);
    else
        imwrite(im,correctimagefilename,'tif','Compression','none','WriteMode','append');
        imwrite(combinedimage,combinedimagefilename,'tif','Compression','none','WriteMode','append');
    end
end
disp(['corrected image saved in ', correctimagefilename]);
disp(['combined image saved in ', combinedimagefilename]);


