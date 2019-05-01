function playback_crosscorrect2(chone,stilloffsets)

maxpixel=max(chone(:));
minpixel=min(chone(:));
stilloffsets2=floor(stilloffsets/2);


for i=1:size(chone,1)
    i
    thisframe=zeros(2*size(chone,2),2*size(chone,3));
    
    starty=max(1+mod(stilloffsets(i,1),2),1+stilloffsets(i,1));
    endy=min(2*size(chone,2),2*size(chone,2)+stilloffsets(i,1));
    startx=max(1+mod(stilloffsets(i,2),2),1+stilloffsets(i,2));
    endx=min(2*size(chone,3),2*size(chone,3)+stilloffsets(i,2));
    
    starty2=max(1,1-stilloffsets2(i,1));
    endy2=min(size(chone,2),size(chone,2)-stilloffsets2(i,1));
    startx2=max(1,1-stilloffsets2(i,2));
    endx2=min(size(chone,3),size(chone,3)-stilloffsets2(i,2));
    
     if i==102
         startx
         starty
   end

    thisframe(starty:2:endy,startx:2:endx)=squeeze(chone(i,starty2:endy2,startx2:endx2));
    thisframe(starty+1:2:endy,startx:2:endx)=squeeze(chone(i,starty2:endy2-mod(stilloffsets(i,1),2),startx2:endx2));
    thisframe(starty+1:2:endy,startx+1:2:endx)=squeeze(chone(i,starty2:endy2-mod(stilloffsets(i,1),2),startx2:endx2-mod(stilloffsets(i,2),2)));
   % thisframe(starty:2:2*endy,startx+1:2:2*endx)=squeeze(chone(i,starty2:endy2,startx2:endx2));
%                
    figure(2);
    imagesc(thisframe);
    colormap gray;
    caxis([minpixel maxpixel]);
    axis image
    pause(.05);
end
