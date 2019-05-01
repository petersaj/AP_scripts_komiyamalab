%this program is meant to create a crib data set which simulates the effect
%of aquiring laser scanning microscopy while a sample is moving
%the result is two files, test.tiff and realmotion.mat. test.tiff is a
%multiframe tiff dataset containing the simulated two photon movie and
%contains not explicit information about the movements.
%realmotion.mat contains four vectors (two of importance).  yofn, xofn
%contain the location of the beam at every pixel aquired in the dataset.
%The motion is simulated with pixel resolution, acknowledging that
%approximating the position of the sample to be constant throughout a
%single line is an assumption of the model.
%feel free to change the various parameters of this simulation in orer to
%explore the behavior of the model in other regimes.

clear all;
%we will use the yfp dataset as a sample that we will simulate aquiring
filename='../data/41307_006.tif';
%a frame from the movie that I judged to be relatively still
scaled=imread(filename,'tiff',18);

%image the frame
figure(5);
imagesc(scaled);
colormap gray;
axis equal
axis tight

%the number of pixels in the sample
N=size(scaled,1);

%we need to simulate motion, so we will be aquiring only over the center of
%the image.  scaled2 is thus our reference image, and we will shift the
%portion of the image scaled2 is sampling from around wiht the motion.
scaled2=scaled(N/4:3*N/4,N/4:3*N/4);


%this is the size of the images we will actually be acquiring
N2=size(scaled2,1);


%we will simulate the motion such that it is gaurunteed to be continuous,
%taux and tauy are the characteristic time constant (in pixels) of the
%first order differential equation driven by a combination of random and
%deterministic variables.
taux=200;
tauy=200;
%tx, and ty are the corresponding fractional elements used in simulating
%the motion.
tx=1-(1/taux);
ty=1-(1/tauy);

%we first create a drive vector and a corresponding isdriving vector
%containing the deterministic portions of the motion.  the drive consists
%of alternating segments of no motion, and segments of sinusoidal motion.
%i've created structures for different timecourses of x and y motion, but
%only make use of one as the motion we see in the brain appears to be
%fairly linear in nature.
drive=[];
isdriving=[];
%loop over 5 periods of no drive
for i=1:5
    %Inter Drive Interval is poisson distributed with a mean of 5 frames, 
    %with a refractory period of 3 frames
    %NOTE: this really should be longer to match our mouse running behavior
    %physiological, but this simulation takes so long as it is i've made it shorter.
    %The consequence of this is that the automatic gain calculation when run
    %with streak=20 frames doesn't work so well and the increased
    %variance created by looking over frames with motion causes the gain to
    %be overestimated (2.7 when it should be 1).  You could decrease this
    %factor in the code before running the HMM, but I wanted to leave it the default
    %for the real data cases.  Also, note however that this does not dramatically
    %affect the results.
    
    IDI=round(-5*N2*N2*log(rand(1)))+3*N2*N2;
    %add on a pad of zeros of the length of the Inter Drive Interval
    drive=[drive [zeros(2,IDI)]];
    %note there is no drive during this period
    isdriving=[isdriving [zeros(1,IDI)]];
    %create a timebase to oscillate over the course of 4 frames with a
    %period of a little less than a frame (1/1.2).
    wt=(2*pi*[1:4*N2*N2]/(N2*N2/1.2));
    %create a sinusodial drive
    ymove=sin(wt);
    %create a different drive for x if you want.. this won't end up being
    %used
    xmove=sin(.75*wt);
    %add the sinusodial drive onto the vector
    drive=[drive [ymove;xmove]];
    %note that there is drive during this period
    isdriving=[isdriving ones(1,length(ymove))];
end
%pad on another non drive period at the end
IDI=round(-5*N2*N2*log(rand(1)))+3*N2*N2;
drive=[drive [zeros(2,IDI)]];
%pad on even more zeros to make the offsets an integer number of frames
%long
intframelength=ceil(length(drive)/(N2*N2))*N2*N2;
drive(:,length(drive):(ceil(length(drive)/(N2*N2))*N2*N2))=0;
isdriving(length(isdriving):length(drive))=0;

%plot out the deterministic portion of the drive
figure(1);
clf;
hold on;
plot(drive');

%initialize the real motion vectors
yofn=zeros(1,length(drive));
xofn=zeros(1,length(drive));
for i=1:intframelength
    %each component of the position obeys the same dynamics with slightly
    %different parameters.
    %the drive contains 3 components
    %1) an always on random component (6*randn(1))
    %2) a random component which is on only during drive (30 or 10)*randn(1)
    %3) a deterministic component, where the x,y components are perfectly
    %correlated but of different magnitudes.. corresponding to motion along
    %a line (7 or .7)*drive(1,i)
    yofn(i+1)=ty*yofn(i)+(1-ty)*(6*randn(1)+(30*randn(1)*isdriving(i))+7*drive(1,i));
    xofn(i+1)=tx*xofn(i)+(1-tx)*(6*randn(1)+(10*randn(1)*isdriving(i))-.7*drive(1,i));
end
%round the positions off to the nearest pixel
%it would be more realistic to leave them floats an interpolate the sample,
%but this program already takes forever, and since we are adding noise
%already i'm not sure how big a difference it would make anyway.
xofn=round(xofn);
yofn=round(yofn);
%plot out the real motion on top of the deterministic part
plot(xofn);
plot(yofn,'r');

%initialize two vectors just to keep track of the actual location of the
%beam, instead of just the relative position of the sample to the optical
%axis.
xoft=zeros(1,length(drive));
yoft=zeros(1,length(drive));
%indices to keep track of line and pixel number
n=0;
m=0;

%calculate the number of frames in the motion sequence
numframes=length(drive)/(N2*N2);

%initialize the matrix for storing the simulated data
data=zeros(numframes,N2,N2);

%loop over frames
for fn=1:numframes
    %loop over lines within a frame
    for j=1:N2
        %increment line counter
        n=n+1
        %loop over pixels within a line
        for i=1:N2
            %increment pixel counter
            n=n+1;
            
            %pick out scaled2 from the portion of scaled corresponding to
            %entire frame if the sample was held at this position
            scaled2=scaled(N/4+yofn(n):3*N/4+yofn(n),N/4+xofn(n):3*N/4+xofn(n));
       
            %note the actual location of the beam on the sample
            yoft(n)=yofn(n)+j;
            xoft(n)=xofn(n)+i;
        
            %the datat that goes in this pixel spot is scaled2(j,i) but
            %sampled with poisson statistics.  5.7 is the gain i calculated
            %for this image, so this should create accurate noise for this
            %image, however you can make 5.7 larger if you want to make the
            %image noiser.  The resulting values should result in an image which was
            %acquired with a gain of 1.  you could rescale it back to add
            %in additional gain if you wanted.
            %NOTE: this is the step that takes the longest in the program
            data(fn,j,i)=random('poiss',double(scaled2(j,i))/5.7);
        end
   
    end
end

%save the crib data to a tiff file
filename='test.tif';
%loop over frames
for i=1:numframes
    %first image the crib data for that frame
    figure(3);
    clf;
    imagesc(squeeze(data(i,:,:)));
    colormap gray;
    pause(.1);   
    
    %write the frame to a tiff file without compression, appending it to
    %the end, so make sure that test.tiff in the scripts directory doesn't
    %exist.
    imwrite(uint16(squeeze(data(i,:,:))),filename, 'TIFF','Compression','none','WriteMode','append');
end

%save the real motion data to a .mat file for future reference.
save realmotion yofn xofn yoft xoft

