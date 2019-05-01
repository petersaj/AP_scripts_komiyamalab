function [final_trace onset] = AH_baselineEstimation(Raw_trace,framerate,detectEvents,compareplot)

% Made into a function by AP
% [final_trace onset] = AH_baselineEstimation(Raw_trace,framerate,detectEvents,compareplot)
% this code is currently kind of a mess


%%% Command file
%%% Set Parameters and March through the procedure
%%% Alex Heitman   Winter Rotation Komiyama Lab

%%% !!!! if the trace is too noisy  there can be some problems !!!! %%%%
%%% !!!! Sometimes we miss a start value that is too high    !!!! %%%%


%%% Input = Raw.mat with [Row=ROIS , Columns=Trace values]
%%% Output = Final_BigEvents.mat and Final_Trace.mat
%%% Final_Trace.mat has final_trace (normalized and straightened flour
%%% trace)
%%% Final_BigEvents.mat has abovethresh (binary matrix of conservative
%%% above/below threhsold )

%%% Still need to do template matching for small events
%%% Still need to convert abovethresh to a true event number
%%% Still need to construct a more aggresive "Accurate" detector
%%% Still need to work on graded events




%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%
%percentile for initial baseline (not used currently)
params.init_baseline=.3; 
% Normed initial trace that we classify as "On State"
params.toohigh=.5;   
% time around "On State" that counts as "On state" in front and behind
% these first two aren't used
params.upticklag=2*round(framerate); 
params.downticklag=4*round(framerate);  
% this one is used symmetrically
params.restsurround = 4*round(framerate);
% minutes of baseline smoothing
params.baseline_window=1;  
% framerate
params.framerate = framerate;

%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1

%%%%%%%%%%%%%%%%%%%%%%%%%



display('First Pass Normalization');
%Final_A1_Norm_Init;

%%% Final_A1_Norm_Init
%%% Input = Raw.mat with a Raw_trace
%%% Uses init_baseline from Params.mat
%%% Raw_trace matrix = (rows are ROIs, columns are flour tr)
%%% does a rough normalization at init_baseline  (shouldn't matter too much)
%%% baseline hack of 30th percentile of trace  (trying to be robust)

%%% Output = Norm.mat  which has Norm_Trace

%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%
plots=0;
trace=Raw_trace;

%%%%%%%%%%%%%%% Identifiy Init_Baseline  %%%%%%%%%

[ROInum,frames]=size(trace);
time_steps=linspace(1,frames,frames);
time_secs=(1/params.framerate)*time_steps;
one_vec=ones(1,frames);
boundary=.1;
backframe=floor((1-boundary)*frames);
frontframe=floor(boundary*frames);
base_matrix=zeros(size(trace));

% This is used to find initial on/off states, doesn't makes sense to be a
% single number, so do by large-scale moving average
smoothtime=params.baseline_window*params.framerate*60; 
for i = 1:ROInum
    base_matrix(i,:) = smooth(trace(i,:),smoothtime);
end
% compensate for edge smoothing by setting edges to first/last reliable point
smooth_edge = ceil(smoothtime/2);
start_smooth = base_matrix(:,smooth_edge);
end_smooth = base_matrix(:,size(base_matrix,2)-smooth_edge);

base_matrix(:,1:smooth_edge) = repmat(start_smooth,1,smooth_edge);
base_matrix(:,end-(smooth_edge-1):end) = repmat(end_smooth,1,smooth_edge);

% for i=1:ROInum
%         values=sort(trace(i,frontframe:backframe));
%         baseline_vec(i)=values(floor(init_baseline*frames));  % the baseline
% end

% % AP EDIT: estimate baseline by mode
% for i = 1:ROInum
%     baseline_vec(i) = mode(round(trace(i,frontframe:backframe)));
%     % sometimes this is zero by error? if so, median
%     if baseline_vec(i) == 0
%        baseline_vec(i) = median(trace(i,frontframe:backframe));
%     end
% end
% base_matrix=baseline_vec*one_vec;

%%%%%%%%%%% Normalize %%%%%%%%%%%%%%
Norm_Trace=( trace - base_matrix )./( base_matrix );

%%%%%%%% Plotting initial baseline %%%%%%%%%%
if plots
    for i=1:16
        r=rem(i,4);
        if r==0
            r=4;
        end
        if r==1
            figure; hold on
            title('Black=Init_Baseline, Red=Raw_Trace')
        end
        j=ceil(rand*ROInum);
        subplot(2,2,r); hold on
        plot(time_steps,Raw_trace(j,:),'r');
        plot(time_steps,one_vec*baseline_vec(j),'k');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 2

%%%%%%%%%%%%%%%%%%%%%%%%%

display('Initial On/Off States');
%Final_A2a_Rest_Init;

%%% Final_A2a_Rest
%%% Input = Norm.mat with a Norm_Trace
%%% "toohigh" from Params.mat
%%% "upticklag" and "downticklag" from Params.mat
%%% Computes an initial estimate of on and off states (using toohigh)
%%% Output = Init_Rest.mat  which has On off state segregation

%%% NOTE If Data is really noisy (noisier than "toohigh") this section may
%produce some bad traces


%%%%%%%%% Intialization %%%%%%%%
plots=0;
% get the on/off states based on smoothed trace
% this initial smoothing value is a bit arbitrary:
% for now based off of looking at decent ones from b-scope (one second)
smooth_size = round(framerate);
smooth_trace = zeros(size(Norm_Trace));
for i = 1:size(Raw_trace,1);
    smooth_trace(i,:) = smooth(Norm_Trace(i,:),smooth_size,'loess');
end

%%%%%%%%%%% Finding Off/On  On/Off Transitions %%%%%%%%
rest_time=ones(size(Raw_trace));
% find where the trace goes over the limit to define an "on" state
overtop=smooth_trace>params.toohigh;
rest_time(overtop)=0;

% set points surrounding over-thresholds to an on state with convolution
rest_surround_kernel = ones(1,round(params.restsurround/2));
on_states = conv2((1-rest_time),rest_surround_kernel,'same');
rest_time(on_states > 0) = 0;
rest_timepts=rest_time;

% Alex used to do this kind of a complicated way

% fbdifference=rest_time(:,1:end-1)-rest_time(:,2:end);
% zerovec=zeros(ROInum,1);
% fbdifference1=[fbdifference zerovec];
% fbdifference2=[zerovec fbdifference];
% [uptickrow, uptickcolumn]=find(fbdifference1==1);
% [downtickrow,downtickcolumn]=find(fbdifference2==-1);
% uptickmatrix=ones(size(trace));
% downtickmatrix=ones(size(trace));
% duration=size(trace,2);
% 
% %%%% Construct On/Off Matrices  then Rest_Timepts matrix %%%%
% for i=1:length(uptickcolumn)
%     firstcolumn=max(1,uptickcolumn(i)-params.upticklag);
%     uptickmatrix(uptickrow(i),firstcolumn:uptickcolumn(i))=0;
% end
% for i=1:length(downtickcolumn)
%     lastcolumn=min(downtickcolumn(i)+params.downticklag,size(trace,2));
%     downtickmatrix(downtickrow(i),downtickcolumn(i):lastcolumn)=0;
% end
% rest_timepts=rest_time.*uptickmatrix.*downtickmatrix;

%%%% Optional Plotting %%%%
if plots==1
    for k=1:12
        j=1+floor( (ROInum-1)*rand(1));
        figure; hold on
        upperbound=max(Norm_Trace(j,1:4000));
        plot(linspace(1,4000,4000),Norm_Trace(j,1:4000),'k',...
            linspace(1,4000,4000),toohigh*(-rest_timepts(j,1:4000)+1),'o')
        axis([0,4000,min(Norm_Trace(j,1:4000)),upperbound])
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 3

%%%%%%%%%%%%%%%%%%%%%%%%%

display('Initial Moving Baseline');
%Final_A2b_MovingBaseline_Init;

%%% Final_A2b_MovingBaseline_Init
%%% Input = Init_Rest.mat , Norm.mat
%%% "baseline_window" from  Params.mat
%%% smoothes out the off-states for a moving baseline
%%% subtracts out the baseline to give us inter_trace

%%% Output = Inter_Traces.mat  with Inter_trace
%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%


%%%%%%%%%%%% Calculate Moving Baseline %%%%%%%%%%
smoothtime=params.baseline_window*params.framerate*60; 
trace_matrix=Norm_Trace;
resttime=rest_timepts;
full_movingmean=zeros(ROInum,frames);
rest_std=zeros(ROInum,1);
for i=1:ROInum
    trace=trace_matrix(i,:);
    rest=find(resttime(i,:));
    resttrace=trace(rest);
    rest_movingmean=( smooth(resttrace,smoothtime) )';
    rest_std(i)=std(resttrace-rest_movingmean);
    j=1;
    for k=1:frames
        if resttime(i,k)==1
            ind=min(frames,j);
            full_movingmean(i,k)=rest_movingmean(ind);
            j=j+1;
        end
        if resttime(i,k)==0
            ind=max(1,j-1);
            try
            full_movingmean(i,k)=rest_movingmean(ind);
            catch me
                disp(['Cell ' num2str(i) ' too noisy, median''ing']);
                full_movingmean(i,:) = median(round(trace));
                break
            end
            j=j+0;
        end
    end
end


%%%%%%%%%%% Create New Normalized Trace %%%%%%%%%%%
inter_trace=trace_matrix-full_movingmean;
movingmean=full_movingmean;


%%%%%%%%%%% Plotting %%%%%%%%%%
if plots
    bins=300;
    for i=1:ROInum
        j=i;
        r=rem(i,4);
        if r==0
            r=4;
        end
        if r==1
            figure;
        end
        
        
        subplot(2,2,r); hold on
        plot(time_steps,Norm_Trace(j,:),'r');
        plot(time_steps,-.5*resttime(j,:)+.5,'y.')   ;
        plot(time_steps,movingmean(j,:),'k');
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 4

%%%%%%%%%%%%%%%%%%%%%%%%%

display('Second Iteration of Baseline Procedure');
%Final_A3a_SecondIter;

%%% Final_A3a_SecondIter
%%% Input = Inter_Traces.mat  and Norm.mat
%%% Uses baseline window
%%% Runs the Crank a second time
%%% On-States defined as Positive Points exceeding the most Negative Point
%%% Smooth out the off-state for approx of the Moving Baseline
%%% Plots diagnostics of moving baseline

%%% Output = Final_Traces.mat,'final_trace','movingmean','rest_std','notrest');

%%%% Initialize %%%%%%%%
stds=rest_std;

trace=inter_trace;

%%%%%%%%%% Find above threshold On-States %%%%%%%%%%%%%
onevec=ones(1,size(trace,2));
stds_matrix=stds*onevec;
normed_trace=trace./stds_matrix;
thresh=-min(normed_trace,[],2);      %%%% Bigger than the most negative
abovethresh=zeros(size(trace));
onset=zeros(size(trace,1),size(trace,2));
offset=zeros(size(trace,1),size(trace,2));
for i=1:size(trace,1)
    smooth1=smooth(normed_trace(i,:),9/frames,'loess');
    above=find(smooth1>thresh(i));
    abovethresh(i,above)=1;
    
    above1=abovethresh(i,1:end-1);
    above2=abovethresh(i,2:end);
    above_backdiff=[above2-above1,0];
    upticks=find(above_backdiff==1);
    onset(i,upticks)=1;
    downticks=find(above_backdiff==-1);
    offset(i,downticks)=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Second Part %%%%%%%%%%%%%%%%%%%%


%%%%%% Initialize %%%%%%%%%%%
smoothtime=params.baseline_window*params.framerate*60;
notrest=abovethresh;
trace_matrix=Norm_Trace;

%%%%% Compute moving_mean and the Final Trace %%%%%%%%
full_movingmean=zeros(ROInum,frames);
rest_std=zeros(ROInum,1);
for i=1:ROInum
    trace=trace_matrix(i,:);
    
    starts=find(onset(i,:));
    if ~isempty(starts)
        for j=1:length(starts)
            zed=max(1,starts(j)-params.upticklag);
            notrest(i,zed:starts(j))=1;
        end
    end
    
    ends=find(offset(i,:));
    if ~isempty(ends)
        for j=1:length(ends)
            zed=min(frames,ends(j)+params.downticklag);
            notrest(i,ends(j):zed)=1;
        end
    end
    
    
    rest=find(notrest(i,:)==0);
    resttrace=trace(rest);
    rest_movingmean=( smooth(resttrace,smoothtime) )';
    rest_std(i)=std(resttrace-rest_movingmean);
    
    j=1;
    for k=1:frames
        
        if notrest(i,k)==0
            ind=min(frames,j);
            full_movingmean(i,k)=rest_movingmean(ind);
            j=j+1;
        end
        
        if notrest(i,k)==1
            ind=max(1,j-1);
            try
            full_movingmean(i,k)=rest_movingmean(ind);
            catch me
                disp(['Cell ' num2str(i) ' too noisy, median''ing']);
                full_movingmean(i,:) = median(round(trace));
                break
            end
            j=j+0;
        end
    end
    
end

final_trace=trace_matrix-full_movingmean;
movingmean=full_movingmean;

%%%%%% Plotting Original Trace, Computed Baseline, Final Trace %%%%%%%%
if exist('compareplot','var')
    if compareplot
        bins=300;
        for i=1:min(100,ROInum)
            j=ceil(ROInum*rand(1));
            r=rem(i,2);
            if r==0
                r=2;
            end
            if r==1
                figure;
            end
            
            
            subplot(2,2,r); hold on
            plot(time_steps,Norm_Trace(j,:),'r');
            % plot(time_steps,notrest(j,:)*thresh(j)*rest_std(j),'y.')   ;
            plot(time_steps,movingmean(j,:),'k');
            
            subplot(2,2,r+2); hold on
            plot(time_steps,final_trace(j,:),'b');
            plot(time_steps,notrest(j,:)*thresh(j)*rest_std(j),'y.')   ;
            plot(time_steps,zeros(1,length(time_steps)),'k');
            
        end
    end
end


% %%%%%% Diagnostic Plots of the Moving Mean %%%%%%%%%%%
% basecorr=corrcoef(movingmean');
% figure; hold on; imagesc(basecorr); colorbar; title('Correlation of Moving Baselines');
%
% basecorrvec=[];
% for i=1:ROInum-1
%     for j=1:i
%         basecorr(i,j)=0;
%     end
%     basecorrvec=[basecorrvec basecorr(i,i+1:end)];
% end
% figure; hold on; hist(basecorrvec);
% title('Histogram of Correlation Structure of Moving Baseline')



%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 5 (optional)

%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('detectEvents','var');
    if detectEvents
        display('Finding Big Events');
        
        %%% Final_A4a_BigEvents  %%%%
        %%% Input = Norm.mat  Final_Trace.mat  Params.mat
        %%% Output= Final_BigEvents.mat
        %%% Thresh= One standard deviation above the most negative standard
        %%% deviation
        
        plots=0;
        stds=rest_std;
        trace=final_trace;
        
        
        onevec=ones(1,size(trace,2));
        stds_matrix=stds*onevec;
        normed_trace=trace./stds_matrix;
        thresh=-min(normed_trace,[],2)+1;      %%% Big Difference
        abovethresh=zeros(size(trace));
        onset=zeros(size(trace,1),size(trace,2));
        offset=zeros(size(trace,1),size(trace,2));
        smoothed1=zeros(size(trace,1),size(trace,2));
        for i=1:size(trace,1)
            smooth1=smooth(normed_trace(i,:),9/frames,'loess');
            smoothed1(i,:)=smooth1;
            above=find(smooth1>thresh(i));
            abovethresh(i,above)=1;
            
            above1=abovethresh(i,1:end-1);
            above2=abovethresh(i,2:end);
            above_backdiff=[above2-above1,0];
            upticks=find(above_backdiff==1);
            onset(i,upticks)=1;
            downticks=find(above_backdiff==-1);
            offset(i,downticks)=1;
        end
        
        final_renormed=normed_trace;
%         save('Final_BigEvents.mat',...
%             'abovethresh','thresh','onset','offset','final_renormed','smoothed1');
        plotnum=40;
        plotroi=ceil((ROInum-1)*rand(1,plotnum));
        if plots==1;
            for i=1:plotnum
                j=plotroi(i);
                r=rem(i,4);
                if r==0
                    r=4;
                end
                if r==1
                    figure;
                end
                subplot(2,2,r); hold on
                plot((time_steps),final_trace(j ,:),'k');
                %plot((time_secs),smooth(final_Ketamine(j,:),9/size(final_Ketamine,2),'loess'),'k');
                if ~isempty( find(abovethresh(j,:)>0) )
                    plot(find(abovethresh(j,:)>0),thresh(j)*rest_std(j),'r.');
                end
                if ~isempty( find(onset(j,:)==1) )
                    plot(find(onset(j,:)==1),1.1*thresh(j)*rest_std(j),'g*');
                end
                if ~isempty( find(offset(j,:)==1) )
                    plot(find(offset(j,:)==1),.9*thresh(j)*rest_std(j),'y*');
                end
                
            end
        end
        
        
%         figure;
%         subplot(1,2,1); hold on
%         imagesc(onset+abovethresh); colormap(gray);
%         xlabel('BigEvents');
%         subplot(1,2,2); hold on
%         imagesc(final_trace); colormap('default');
%         xlabel('Straightened Flourescence Trace')
        
    end
end

disp('Finished baseline normalization')
