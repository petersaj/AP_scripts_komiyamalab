%double checks that you loaded a _PI file and asks you to load one if you
%haven't

%you also need to have lambda set before you run this... normally run
%maxlambda, and have it call this script.. but you can set lambda by hand if you want.
if ~exist('PIsave')
    if ~exist('filename')
      [filename,pathname]=uigetfile('*PI.mat','pick your file');
      junk=strfind(filename,'\');
      filebase=filename(1:end-4);
    else
      %pathname='../';
        junk=strfind(filename,'\');
        filebase=filename(junk(end)+1:end-4);
        filename=['matfiles\' filebase '_PI.mat'];
    end
     load([pathname filename]);
end

if ~exist('h')
    h = waitbar(0.0,['Run with Lambda=' num2str(lambda)]);
    set(h,'Name',['HMM Lambda=' num2str(lambda)]);
end
%creates the basic exponential model for the transistion probabilities,
%this in terms of relative change in offset.  We normalize appropriately
%and make sure its big enough so that it covers all possible differences in
%offsets
[xx,yy]=meshgrid(-2*maxdx:2*maxdx,-2*maxdy:2*maxdy);
%TK090224 correct for uneven sampling
ratio = M/N;
%TK090224 assume X is rostral-caudal and Y is medial-lateral, and X has
%twice as much motion and correct for it
switch motionBias
    case 1 % no bias
        ratio; % do nothing
    case 2 % motion bias for X
        ratio = ratio*2; % more X than Y
    case 3 % motion bias for Y
        ratio = ratio/2; % more Y than X
end
rr=sqrt(xx.^2+(yy*ratio).^2);
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
%figure(20);
%imagesc(trans_prob);

frames=1:numframes;

%this is going to be the matrix which keeps track of the transition which
%was the most probable way to get to a particular state from the previous
%state.. we will fill this matrix up, and the backtrack the most probable
%path as is the standard way of the veterbi algorithm.
savemax=zeros(Nx*Ny,numframes*N,'int16');

%P is the probability vector which describes the maximal probability of
%being in a particular state at the current time step as we march forward
%we will start with a uniform distribution across offsets
P=ones(1,Nx*Ny);

%vector to save the time it took to process each frame
tocs=zeros(numframes,1);

%index which will march over lines considered
m=0;

%loop over all frames
h = waitbar(0,'Running HMM');
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
        %Pnew = makepifast(trans_prob,P);
        %Pnew = Pnew';

        %calculate the most probable way to wind up in a given state, and
        %what the probability is.. save which path you took to get that
        %value, and update P.
        [P,savemax(:,m)]= max(Pnew',[],1);

        %now add in the fit data to adjust the probability of being in a
        %given state.

        %pull out the relevant matrix of values
        PI=squeeze(PIsave(m,:,:));
        %reshape them into a vector
        PIvec=reshape(PI,Nx*Ny,1);

        %i shift all the probabilities by the mean just to keep things
        %from hitting round off errors.. this doesn't affect the
        %calculation, just keeps things reasonable.
        PIvec=PIvec/gain;
        PIvec=PIvec-mean(PIvec);


        %add on the fits to the probabilities

        P=P+PIvec';

        % plot out the maximal probabilities and the fits if you want
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
    minremain=(numframes-i)*delframes/60;
    %display what frame you are on, how long it should take, and how long
    %the current frame took
     waitbar((i/length(frames)),h,['Running HMM.. Frame ' num2str(i) ' ETA:' num2str(minremain) ' min']);
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
PItot=0;

%march backward from the last line considered to the first
for k=numlines:-1:1
    if(mod(k,N)==0)
%     waitbar((k/numlines),h,['Backtracing path.. line ' num2str(k)]);
    end
    %pull out the fits from the current line
    PI=squeeze(PIsave(k,:,:));
    %turn it into a vector
    PIvec=reshape(PI,Nx*Ny,1);

    %what is the fit for the most probable state at this timepoint
    PIpath(k)=PIvec(mprob);
    %add that to the total fit
    PItot=PItot+PIvec(mprob);

    %save that point along the path
    thepath(k)=mprob;
    %remember what was the most probable way was to get to that state.. update mprob
    mprob=savemax(mprob,k);
end

%unhash the path in terms of state index into a pair of offsets
offsets(1,:)=mod(thepath,Ny);
offsets(1,find(offsets(1,:)==0))=Ny;
offsets(1,:)=offsets(1,:)-maxdy-1;
offsets(2,:)=((thepath-mod(thepath,Ny))/Ny)-maxdx;

%adjust for the alignments between reference frames
% offsetfix=1:numlines;
% offsetfix=floor(offsetfix/(N-2*edgebuffer));
% offsetfix=floor(offsetfix/framesper30sec)+1;
% offsetfix=stilloffsets(offsetfix,:)';
% offsets=offsets-offsetfix;

%plot out the offsets
figure(2);
plot(offsets');
title(['\lambda' num2str(lambda)]);
xlabel('line number');
ylabel('offset (relative pixels)');

%save the results to an F file which is named according to the value of
%lambda used.
%pathname='../';
[markov_save_path, markov_save_filebase] = fileparts([filename '.mat']);
savename=[savedir filesep markov_save_filebase '.mat'];
save(savename,'offsets','thepath','P','PItot','PIpath','totprob');

close(h);
