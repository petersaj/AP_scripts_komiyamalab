function [fixeddata,countdata]=playback_markov(imagedata,offsets,edgebuffer,playback,totvel)
%fixeddata,countdata]=playback_markov(imagedata,offsets,edgebuffer,playback,totvel)

%pick out the number of frames and size of images
numframes=size(imagedata,1);
N=size(imagedata,2);
M=size(imagedata,3);


%dynamically set the dynamic range for playback by taking the 99.99% pixel
%removes outliers. can adjust downward for automated gain control
sortpix=imagedata(1:3,:,:);
sortpix=sortpix(:);
sortpix=sort(sortpix);
k=round(length(sortpix)*.9999);
maxpixel=sortpix(k)
minpixel=sortpix(1)


%initialize the data 
countdata=int8(zeros(numframes,N,M));
fixeddata=zeros(numframes,N,M);

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
        lineextract=double(squeeze(imagedata(i,j,:)));
 
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
    disp(i);
    %store the results into the data structure
    fixeddata(i,:,:)=correctimage;
    countdata(i,:,:)=countimage;

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



