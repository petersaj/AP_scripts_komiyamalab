function split_createPI_ManRef(fullfilename,fullPIsavefilename,fullreferencefilename,frameanalysismode, maxdx,maxdy,numchannels,lambda)
%split_createPI(filename,maxdx,maxdy,numchannels)
%FILENAME is the file located in the data directory you want to correct
%MAXDX is the maximum offset in pixels you estimate is observed in the X
%direction (across columns)
%MAXDY is the maximum offset in pixels you estimate is observed in the Y
%direction (across rows)
%these two parameters should be set high enough so that the maximal offset
%is never reached, but as small as possible as the running time increases
%and
%accuracy of placement decreases as these values get larger.
%NUMCHANNELS is the number of channels in the movie
tic;


%these two parameters are the maximum offsets to consider in pixels across
%the whole movie, these should be set high enough that that offset is never
%reached, but as small as possible as the running time increases and
%accuracy of placement decreases as these values get larger
%maxdx=input('max delta X: ');
%maxdy=input('max delta Y: ');
% maxdx=5;
% maxdy=5;

Nx=2*maxdx+1;
Ny=2*maxdy+1;

%creates the full filename with the path, and extracts the number of frames
%etc from the tiff file information
imageinfo=imfinfo(fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;


msecperframe=1.24*N; %fixed from 2 to 1.24, AP 8/16/11
framesper30sec=floor(30000/msecperframe);%check Takashi

stop=false;

% SECTION COMMENTED TO LOADING IN SINGLE FRAMES
% h = waitbar(0.0,'Loading Data...');
% set(h,'Name','Loading Data...');
% 
% 
% % if its a two channel movie the data must be loaded appropriately
% if numchannels==2
%     chone=zeros(numframes/2,N,M);
%     %chtwo=zeros(numframes/2,N,M);
%     for i=1:2:numframes
%         if mod(i,10)==1
%             waitbar(i/numframes,h,['Loading Frame... ' num2str(i)]);
%         end
%         chone((i+1)/2,:,:)=imread(fullfilename,'tiff',i);
%         %chtwo((i+1)/2,:,:)=imread(fullfilename,'tiff',i+1);
%     end
%     if exist('h') delete(h); end
%     %make sure there are no zeros (ln(0) is bad...)
%     %this was put in to keep crib data from breaking the program
%     %chtwo(find(chtwo<1))=1;
%     chone(find(chone<1))=1;
% 
%     numframes=numframes/2;
% 
%     %find_gain takes a movie and calculates the gain by comparing the std and
%     %mean of pixels over a "stillest" section of movie
%     %in this case 50 frames as indicated by the second parameter passed.
%     %see find_gain for more details...
%     %if you have a manaul calibration you could substitute that here..
% 
     for im_sample = 1:numchannels:27*numchannels
         im_gain(ceil(im_sample/numchannels),:,:) = double(imread(fullfilename,'tiff',im_sample));
     end
     gain=find_gain(im_gain,25);%regardless of N of channels, use green channel.
     disp(['gain = ', num2str(gain)]);
%     %gain = 10.3;
% 
%     %16.88
% 
%     %I used to normalize the data here so it is in units of photons
%     %however, it makes more sense to store the raw data and do the
%     %normalization after when you do the HMM in case you want to put in a
%     %manual gain factor from a calibration experiment.
%     %chtwo=chtwo/gain;
%     %break the movie up into 30 second intervals
% 
%     
%     %numstills=floor(numframes/framesper30sec)+1;
%     %stillimages are the reference frames which we will be aligning too
%     %     for i=1:numstills
%     %         [minval,minimage]=min(mean(mean(abs(diff(chtwo(1+(i-1)*framesper30sec:min(i*framesper30sec,numframes),:,:),1)),2),3));
%     %         stillimages(i,:,:)=chtwo(minimage+(i-1)*framesper30sec,:,:);
%     %     end
%     %     %if its a one channel movie it needs to be loaded another way
%     %but the rest of the steps are very similar.
% else
%     chone=zeros(numframes,N,M);
%     for i=1:numframes
%         if mod(i,10)==0
%             waitbar(i/numframes,h,['Loading Frame... ' num2str(i)]);
%         end
%         chone(i,:,:)=imread(fullfilename,'tiff',i);
%         %chtwo((i+1)/2,:,:)=imread(filename,'tiff',i+1);
%     end
%     if exist('h') delete(h); end
%     chone(find(chone<1))=1;
%     gain=find_gain(chone,25);%10.3
%     disp(['gain = ', num2str(gain)]);
%     %make sure there are no zeros .. ln(0) is bad
% end

if frameanalysismode
    clear chone;
end

h = waitbar(0.0,'Getting reference image');
set(h,'Name','Getting reference image');
frames=1:numframes;

if ischar(fullreferencefilename)
    if ~exist(fullreferencefilename,'file')
        AverageFrame = zeros(N,M);
        if frameanalysismode
            for i=frames
                waitbar(i/numframes,h,['Calculating reference image...']);
                if numchannels==2
                    temp = imread(fullfilename,'tiff',2*i-1);
                else
                    temp = imread(fullfilename,'tiff',i);
                end
                AverageFrame = AverageFrame+double(temp);
            end
            close(h);
            AverageFrame = AverageFrame/numframes;
        else
            AverageFrame = squeeze(mean(chone,1));
        end
        refimage = round(AverageFrame); clear AverageFrame
    else
        refimage = double(imread(fullreferencefilename,'tiff'));
    end
elseif isnumeric(fullreferencefilename)
    refimage = fullreferencefilename;
end


edgebuffer=maxdy;

%this function aligns all the reference frames to another using maximal
%cross correlation looking over a range of offsets (in this case 5)
%hopefully these values should all be very small as the reference frames
%have basically the same overall shape.. we found this to be true in our
%data but it is not gaurunteed.

waitbar(0.0,h, 'Reference Images...');
set(h,'Name','Calculating PI');

waitbar(0.1,h,'Aligning Reference Images...');
%stilloffsets=refs_align(stillimages,5)


%this is the total number of lines we will consider placing since we don't
%try to place the top/bottom maxdy lines
numlines=(N-2*edgebuffer+1)*numframes;

%this is the matrix for which we will store the fits for each offset.
PIsave=zeros(numlines,Ny,Nx);
%m is the index we will use to order the lines for which we have placed
m=0;
%this stores which line of the data is at which index
msave=zeros(numlines,2);
%this stores the times it took to do each frame for my crude progress bar
tocs=zeros(numframes,1);

waitbar(0.2,h,'Calculating PI...');
for i=1:numchannels:numframes
    %note the time at the start of the frame
    tic;
    
    % load in one frame at a time for memory issues
    chone(1,:,:)=double(imread(fullfilename,'tiff',i));
    
    %which reference frame are we aligning to
    %stillindex=floor(i/framesper30sec)+1;
    %pull out that reference frame
    %refimage=squeeze(stillimages(stillindex,:,:));

    %scan over all the lines in that frame that we are considering
    for j=1+edgebuffer:N-edgebuffer
        %increment our index
        m=m+1;
        %note which line it is
        msave(m,:)=[i j];

        %pull out the line we are placing

        if frameanalysismode
            if numchannels==1
                temp = imread(fullfilename,'tiff',i);
            elseif numchannels==2
                temp = imread(fullfilename,'tiff',i*2-1);
            end
            lineextract = double(squeeze(temp(j,:)))';
        else
            lineextract=squeeze(chone(1,j,:));
        end

        %this function takes the line we are trying to fit, the reference
        %image.. the maximum offsets, and the expected position of the line
        %and returns the log fit probabilities into the PIsave matrix
        %see create_fit_markov for more details
        PIsave(m,:,:)=create_PI_markov(lineextract,refimage,maxdx,maxdy,j);

        %this is if you wanted to visualize the fits as they are being made
        %useful for debugging, but slow.
        %       figure(1);
        %       imagesc(squeeze(PIsave(m,:,:)));
        %       %colorbar;
        %       pause(.05);
    end
    %note the time at the end of the frame
    tocs(i)=toc;

    %estimate the average time it took for frames to be processed
    delframes=mean(tocs(max(i-9,1):i));
    %use that to estimate the time remaining
    minremain=(numframes-i)*delframes/60;
    %show what frame you are on, how much time is remaining and how much
    %time that frame took
    waitbar(0.2+.7*(i/length(frames)),h,['PI: Frame ' num2str(i) ', ETA ' num2str(minremain) ' min']);
    %disp([i minremain tocs(i)]);

end



%find where in the fullfilename the filename is, and cut out the string of
%the filename without the .tif
%use that to save every variable into a _PI.mat file in one directory up,
%then down to matfiles.

junk=strfind(fullfilename,'\');
if ~isempty(junk)
    filebase=fullfilename(junk(end):end-4);
else
    filebase=fullfilename(1:end-4);
end
disp('Saving results of PI calculation');
waitbar(0.9,h,'Saving results of PI');
save([fullPIsavefilename '_PI.mat']);
waitbar(1,h,'Done!');
if exist('h') delete(h); end
