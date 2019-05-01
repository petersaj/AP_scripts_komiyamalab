function playback_crosscorrect(chone,stilloffsets)

maxpixel=max(chone(:));
minpixel=min(chone(:));


for i=1:size(chone,1)
    
    thisframe=zeros(size(chone,2),size(chone,3));
    
    starty=max(1,1-stilloffsets(i,1));
    endy=min(size(chone,2),size(chone,2)-stilloffsets(i,1));
    startx=max(1,1-stilloffsets(i,2));
    endx=min(size(chone,3),size(chone,3)-stilloffsets(i,2));
    
    starty2=max(1,1+stilloffsets(i,1));
    endy2=min(size(chone,2),size(chone,2)+stilloffsets(i,1));
    startx2=max(1,1+stilloffsets(i,2));
    endx2=min(size(chone,3),size(chone,3)+stilloffsets(i,2));
      
    thisframe(starty:endy,startx:endx)=squeeze(chone(i,starty2:endy2,startx2:endx2));
    figure(2);
    imagesc(thisframe);
    colormap gray;
    caxis([minpixel maxpixel]);
    axis image
    pause(.05);
end
