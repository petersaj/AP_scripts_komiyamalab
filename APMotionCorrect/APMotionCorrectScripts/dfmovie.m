avgimage=squeeze(sum(fixeddata)./sum(countdata));


for i=1:size(fixeddata,1)
    
    figure(1);
    clf;
    
    thisframe=squeeze(fixeddata(i,:,:))./double(squeeze(countdata(i,:,:)));
    thisframe=thisframe./avgimage;
    thisframe=thisframe.*(squeeze(countdata(i,:,:)>0));
    
    imagesc(thisframe);
    caxis([1 2]);
    colormap gray;
    pause(.1);
end

