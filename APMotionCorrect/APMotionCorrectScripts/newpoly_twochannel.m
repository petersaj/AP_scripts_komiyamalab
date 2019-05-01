clear all
edgebuffer=7;

[offsetfile,offsetpath]=uigetfile('*.mat','pick the offset file');
fulloffsetfile=[offsetpath offsetfile];
load(fulloffsetfile);

figure(60);
plot(offsets');
drawnow;
[tifffilename,tiffpath]=uigetfile('*.tif','pick your tiff file');
fullfilename=[tiffpath tifffilename];
imageinfo=imfinfo(fullfilename);

channels=[1 1 1 1 1 1];
loadchannels;

imageinfo=imfinfo(fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

chone=zeros(numframes/2,N,M);
chtwo=zeros(numframes/2,N,M);
for i=1:2:numframes
  if mod(i,10)==1
    disp(i);
  end
  chone((i+1)/2,:,:)=imread(fullfilename,'tiff',i);
  chtwo((i+1)/2,:,:)=imread(fullfilename,'tiff',i+1);
end
numframes=numframes/2;

% chone=zeros(numframes,N,M);
% for i=1:numframes
%   if mod(i,10)==1
%     disp(i);
%   end
%   chone(i,:,:)=imread(fullfilename,'tiff',i);
% end



[fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,0);
[fixeddata2,countdata2]=playback_markov(chtwo,zeros(size(offsets)),edgebuffer,0);
%[F,Fr,colors,X,Y]=extract_F(fixeddata,countdata,squeeze(mean(fixeddata2)./mean(countdata2)));
[F,Fr,masks,colors]=extract_F_poly(fixeddata,countdata,squeeze(mean(fixeddata2)./mean(countdata2)));
numcolors=size(colors,2);

Fc=zeros(size(F));
for j=1:numcolors

    junk=F(:,j);
    junk=junk;

    window=round(110*(64/N));
    junk2=zeros(size(junk));
    for k=1:length(junk)
        cut=junk(max(1,k-window):min(numframes,k+window));
        cutsort=sort(cut);
        a=round(length(cut)*.08);
        junk2(k)=cutsort(a);
    end
    Fc(:,j)=(junk-junk2);
  
end


samplingrate=1; %in Khz
binwidth=.025; %bin velocity/stimulus data in seconds
framestarts=find(diff(diff(data(:,5))<-.1)==1);
IFI=round(mean(diff(framestarts)));
if length(framestarts)==numframes
  framestarts=[framestarts(1)-IFI;framestarts];
end

%time=(1:length(datacut))/(samplingrate*1000);



medians=median(data(:,1:4))
medmatrix=repmat(medians,length(data),1);
datanew=round((data(:,1:4)-medmatrix)/.025);
figure(11);
clf;
hold on;
plot(datanew(:,2),'b');
   
for j=1:4
    a=datanew(:,j);
    position=cumsum(a);
    newposition=spline(1:75:length(position),position(1:75:end),1:length(position));
    figure(j);
    clf;
    hold on;
    plot(position);
    plot(newposition,'r');
    drawnow;
    pause(.1);
    
    anew=diff(newposition);
    datanew(:,j)=[0 anew];
end
figure(11);
plot(datanew(:,2),'g');

data(:,1:4)=datanew;

samplingrate=1; %in Khz
framestarts=find(diff(diff(data(:,5))<-.1)==1);
IFI=round(mean(diff(framestarts)));
if length(framestarts)==numframes
    framestarts=[framestarts(1)-IFI;framestarts];
end

data=data(framestarts(1):framestarts(end)+IFI,:);
framestarts=framestarts-framestarts(1)+1;

totvel=sqrt(data(:,1).^2+data(:,2).^2+data(:,4).^2);
totvel=totvel/137;
figure(1);
clf;
hold on;
plot(totvel);

window=750;
samplingrate=1;
time=(1:length(data))/(samplingrate*1000);
acctime=[1+window:window:length(data)-window];
acctime=time(acctime);
lt=[2*(1+edgebuffer):2:2*(N-edgebuffer)];
linetimes=[];
for i=1:numframes
    linetimes=[linetimes lt+framestarts(i)-framestarts(1)];
end
linetimes=linetimes/(1000);
j=0
totacc=zeros(1,length(totvel)/(window));
for i=1+window:window:length(totvel)-window
%for i=5E4
    j=j+1;
    if mod(j,10)==0
        disp([i length(totvel)]);
    end
    cut=totvel(i-window:i+window);    
    fitmodel=fit([-window:window]',cut,'a*x+b','StartPoint',[(cut(1+window)-cut(1))/window cut(1)]);
    totacc(j)=fitmodel.a;
    figure(1);
    x=[-window:window];
    plot(x+i,fitmodel.a*x+fitmodel.b,'r');
    drawnow;
end
totacc=totacc(1:j);

totaccsplined=spline(acctime,totacc,linetimes);
endtime=find(linetimes>acctime(end));
endtime=endtime(1);


velbin=zeros(1,numframes);
for i=1:numframes
    velbin(i)=mean(totvel(framestarts(i):framestarts(i+1)));
end
displace=(1:size(F,2));
displace=repmat(displace,numframes,1);

order=[1:numcolors];
figure(2);
clf 
hold on;
for j=1:length(order)

  subplot('position',[.05 .38 .93 .6]);

  hold on;
  
  H=plot((1:numframes),5*(Fc(:,order(j))-1)+2*displace(:,j));
  
  set(H,'Color',colors(:,order(j)));
  xlim([0 numframes]);
end
axis tight;
subplot('position',[.05 .05 .93 .06]);
plot(linetimes*1000/(2*N),offsets(1,:));
xlim([0 numframes]);
ylabel('Y');

subplot('position',[.05 .12 .93 .06]);
plot(linetimes*1000/(2*N),offsets(2,:));
xlim([0 numframes]);
set(gca,'Xtick',[]);
ylabel('X');

subplot('position',[.05 .20 .93 .06]);
plot(time*1000/(2*N),data(:,6));
set(gca,'Xtick',[]);
xlim([0 numframes]);

subplot('position',[.05 .26 .93 .1]);
plot(time*1000/(2*N),totvel);
xlim([0 numframes]);
axis tight;
set(gca,'Xtick',[]);


% 
% 
 Fname=sprintf('%s_Fpoly.mat',tifffilename(1:end-4))
 fullFname=[tiffpath Fname];
 save(fullFname,'F','Fr','Fc','N','M','masks','colors','data','framestarts','IFI','velbin','totacc','totaccsplined','linetimes','totvel','numframes','numcolors');
% 
% cellFname=sprintf('%s_cells_poly.jpg',tifffilename(1:end-4));
% fullcellFname=[tiffpath cellFname];
% figure(45);
% saveas(gcf,fullcellFname,'jpg');