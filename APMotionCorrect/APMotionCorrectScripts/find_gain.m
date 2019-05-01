function gain=find_gain(data,streak)
%takes a movie (DATA) that is numframes x N x N and calculates the GAIN of
%the optical system which took the movie by looking at the relationship
%between the standard deviation and the mean of those pixels.
%this of course only works if there is no movements, so this function look
%for the "stillest" section of frames (STREAK frames long).  "stillest" is
%judged by having the smallest mean absolute difference between frames.

h = waitbar(0.0,'Calculating Mean Absolute Difference Between Frames...');
set(h,'Name','Auto Calculating Gain...');
set(h,'Position',[50 0 360 72]);

%pull out the number of frames and size of each frame
%assumed to be square.
 numframes=size(data,1);
 N=size(data,2);

 %calculate the mean absolute difference between frames
 diffvector=mean(mean(abs(diff(data)),2),3);
 
 %initialize the vector which stores the mean absolute difference between
 %frames for all the possible streaks
 diffavg=zeros(1,length(diffvector)-streak);
 
 %fill in that vector
 for i=1:numframes-streak-1
     waitbar(.1+.8*(i/(numframes-streak-1)),h,['Calculating Mean Diff. Starting at Frame ' num2str(i)]);
     diffavg(i)=mean(diffvector(i:i+streak));
 end
 
 %find the best streak
 [junk,minframe]=min(diffavg);
 set(h,'Name',['AutoCalculating Gain...Streak at ' num2str(minframe)]);
 %cut out the data across that streak
 datacut=data(minframe:minframe+streak,:,:);
 
 waitbar(.9,h,['Calculating Means..']);
 %find the mean of every pixel across that streak
 means=mean(datacut);
 
 %turn it into a vector...
 means=means(:);
 
 waitbar(.925,h,['Calculating Stds..']);
%calculate the stds of every pixel across that streak and turn it into a
%vector
 stds=std(datacut);
 stds=stds(:);
 
 waitbar(.95,h,['Fitting to Sqrt Function..']);
 %fit a sqrt function to the relationship between the mean and stds
 curve=fit(means,stds,'a*x^.5');
 
% the gain is related to a^2
gain=curve.a.^2;
waitbar(1,h,['Done!']);

if exist('h') delete(h); end
  
 