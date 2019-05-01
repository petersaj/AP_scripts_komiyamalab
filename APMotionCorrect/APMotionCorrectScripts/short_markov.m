
%double checks that you loaded a file and asks you to load one if you
%haven't
%you also need to have lambda set before you run this... normally run
%maxlambda that calls this.. but you can set lambda by hand if you want.
if ~exist('PIsave')
    if ~exist('filename')
        [filename,pathname]=uigetfile('*Pi.mat','pick your file');
        junk=strfind(filename,'/');
        filebase=filename(1:end-4);
    else
        pathname='../';
        junk=strfind(filename,'/');
        filebase=filename(junk(end)+1:end-4);
        filename=['matfiles/' filebase '_PI.mat'];
    end
    load([pathname filename]);
end


%creates the basic exponential model for the transistion probabilities,
%this in terms of relative change in offset, normalizing appropriately
%we make sure its big enough so that it covers all possible differences in
%offsets
[xx,yy]=meshgrid(-2*maxdx:2*maxdx,-2*maxdy:2*maxdy);
rr=sqrt(xx.^2+yy.^2);
rel_trans_prob=exp(-(abs(rr)./lambda));
rel_trans_prob=rel_trans_prob/(sum(sum(rel_trans_prob)));

%to plot it out if you wanted
%figure(9);
%imagesc(rel_trans_prob);

%now build up the entire transition probability matrix where you index a
%pair of offsets as a single hashed value by making use of the reshape
%function in matlab.  an offset pair will now be refered to a state.
trans_prob=zeros(Nx*Ny,Nx*Ny);
for i=-maxdx:maxdx
    for j=-maxdy:maxdy
        trans_prob((i+maxdx)*(maxdy*2+1)+j+maxdy+1,:)=reshape(rel_trans_prob((maxdy-j+1):(3*maxdy-j+1),(maxdx-i+1):(3*maxdx-i+1)),Nx*Ny,1);
    end
end
%there are some alternative ways you could consider renormalizing this
%matrix, but we decided not to use any of these..
%trans_prob=trans_prob./repmat(sum(trans_prob),Nx*Ny,1);
%trans_prob=trans_prob./repmat(sum(trans_prob,2),1,Nx*Ny);

%translate it into a log probability
trans_prob=log(trans_prob);

%image the whole matrix if you want
% figure(20);
% imagesc(trans_prob);
% colorbar;

frames=1:20;
%this is going to be the matrix which keeps track of the transition which
%was the most probable way to get to a particular state from the previous
%state.. we will fill this matrix up, and the backtrack the most probable
%path as is the standard way of the veterbi algorithm.
savemax=zeros(Nx*Ny,length(frames)*N,'int16');

%P is the probability vector which describes the maximal probability of
%being in a particular state at the current time step as we march forward
%we will start with a uniform distribution across offsets
P=ones(1,Nx*Ny);

%vector to save the time it took to process each frame
tocs=zeros(numframes,1);

%index which will march over lines considered
m=0;

%loop over all frames
for i=frames
    %note the time
    tic
    %loop over all the lines considered in that frame
    for j=1+edgebuffer:N-edgebuffer

        %increment our index for lines considered
        m=m+1;

        %replicate the starting probabilities at the time step before for
        %all the possible offsets
        Prep=repmat(P,Nx*Ny,1);
        %calculate the matrix of probabilities of being in one state at
        %the previous time point and then transitioning to another
        Pnew=Prep+trans_prob;

        %this was my C implementation of the previous 2 lines.. it speeds
        %things up but i have disabled it here... if you want to compile
        %the C for your platform you can.. the file there.
        % Pnew = makepifast(trans_prob,P);

        %calculate the most probable way to wind up in a given state, and
        %what the probability is.. save which path you took to get that
        %value, and update P.
        [P,savemax(:,m)]= max(Pnew',[],1);

        %now add in the fit data to adjust the probability of being in a
        %given state.

        %pull out the relevant matrix of values
        PI=squeeze(PIsave(m,:,:));

        %          figure(10);
        %          imagesc(PI);
        %          drawnow;
        %          pause(.1);

        %reshape them into a vector
        PIvec=reshape(PI,Nx*Ny,1);

        %i shift all the probabilities by the mean just to keep things
        %from hitting round off errors.. this doesn't affect the
        %calculation, just keeps things reasonable.
        PIvec=PIvec/gain;
        PIvec=PIvec-mean(PIvec);



        %add on the fits to the probabilities

        P=P+PIvec';

        %plot out the maximal probabilities and the fits if you want
        % figure(60);
        % imagesc(reshape(P,Nx,Ny));
        % figure(61);
        % imagesc(reshape(PIvec,Nx,Ny));
    end

    %note the time
    tocs(i)=toc;

    %calculate the average time per frame so far
    delframes=mean(tocs(max(i-9,1):i));
    %use that to estimate the time remaining
    minremain=(length(frames)-i)*delframes/60;
    %display what frame you are on, how long it should take, and how long
    %the current frame took
    waitbar((jj-1)/length(lambdas)+dl*(i/length(frames)),h,['Lambda=' num2str(lambda) ' Frame ' num2str(i)]);
    %disp([i minremain tocs(i)]);

end



clear offsets;
numlines=m;

%find the state that was the most probable ending point
%and what the total probablity was for this value of lambda
[totprob,mprob]=max(P);

%initialize the path of most probable states
thepath=zeros(1,numlines);

%for my interest i save the fits along the path
PIpath=zeros(1,numlines);
%calculate the total fit for this path without considering transition
%probabilities.
Ptot=0;

%march backward from the last line considered to the first
for k=numlines:-1:1
    %pull out the fits from the current line
    PI=squeeze(PIsave(k,:,:));
    %turn it into a vector
    PIvec=reshape(PI,Nx*Ny,1);

    %what is the fit for the most probable state at this timepoint
    PIpath(k)=PIvec(mprob);
    %add that to the total fit
    Ptot=Ptot+PIvec(mprob);

    %save that point along the path
    thepath(k)=mprob;
    %remember what was the most probable way was to get to that state.. update mprob
    mprob=savemax(mprob,k);
end

%unhash the path in terms of state index into a pair of offsets
offsets(1,:)=mod(thepath,Ny);
offsets(1,find(offsets(1,:)==0))=Ny;
offsets(1,:)=offsets(1,:)-maxdy-1;;
offsets(2,:)=((thepath-mod(thepath,Ny))/Ny)-maxdx;

%adjust for the alignments between reference frames
% offsetfix=1:numlines;
% offsetfix=floor(offsetfix/(N-2*edgebuffer));
% offsetfix=floor(offsetfix/framesper30sec)+1;
% offsetfix=stilloffsets(offsetfix,:)';
% offsets=offsets-offsetfix;

%plot out the offsets
figure(2);
set(gcf,'Position',[50 200 750 500]);
plot(offsets');
title(['\lambda' sprintf(' = %6.5f',lambda)]);
xlabel('line number');
ylabel('offset (relative pixels)');
