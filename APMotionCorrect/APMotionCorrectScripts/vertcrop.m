clear all

if ~exist('filename')
    [filename,path]=uigetfile('*.tif','pick your tif file');
else
    path='../';
end

fullfilename=[path filename];
imageinfo=imfinfo(fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

frameone=imread(fullfilename,'tiff',1);

figure(1);
clf;
imagesc(frameone);
axis equal;
axis tight;
colormap gray;

[x,y]=ginput(1);

x=floor(x);

chone=zeros(numframes,N,M);
for i=1:numframes
       
        chone(i,:,:)=imread(fullfilename,'tiff',i);
        
end

chone=chone(:,:,x:M);

fullfilename=[fullfilename(1:end-4) '_crop.tif'];

for i=1:numframes
    imwrite(squeeze(uint16(chone(i,:,:))),fullfilename,'tiff','Compression','none','WriteMode','append');
end
    