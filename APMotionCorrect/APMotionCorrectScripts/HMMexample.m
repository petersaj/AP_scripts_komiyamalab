
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Two channel example (chone=Calcium Green AM, chtwo=SR101
%This is the data from Figure 4 of the Neuron paper)
% uncomment lines 9-16
% OR uncomment lines 19-24 if you have previously calculated PI
% OR lines 27-31 if you have calculated both PI and the offsets 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %example of start to finish
clear all
%close all
pathname = 'C:\Users\Andy\Documents\MATLAB\AnalysisCodeForAndy\MotionCorrection\testfiles\';
[filename,path]=uigetfile('*.tif','pick your tif file');
split_createPI(filename,3,3,2);
load([pathname filename(1:end-4) '_PI.mat']);
maxlambda;
[fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);

% %example of running the HMM part after calculating the fits
% clear all
% close all
% load('../matfiles/31007_012_PI.mat');
% maxlambda;
% [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);

% %example of playing back the result without doing any further calculation
% clear all
% close all
% load('../matfiles/31007_012_PI.mat');
% load('../matfiles/31007_012_L0.32.mat');
% [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%One channel example (chone=YFP under Thy1 in cortical neurons)
%(uncomment lines 35-40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% split_createPI('41307_006_crop.tif',4,17,1);
% load('../matfiles/41307_006_crop_PI.mat');
% maxlambda;
% [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);
% 
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %crib data example (uncomment lines 45-81)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% split_createPI('test.tif',4,13,1);
% load('../matfiles/test_PI.mat');
% maxlambda;
% [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);
% 
% %load the actual motion that was used to create this crib dataset
% %see makecribdata.m for details on how exactly it was created
% load realmotion
% 
% %the offsets don't contain predictions for the lines at the edges
% %to avoid any problems with a lack of a proper place within the reference 
% %image to place a line, so we need to create a vector which contains the
% %absolute line number for the predictions we've made
% 
% %the absolute line numbers for every line scanned
% linenumbers=1:(length(yofn)/N);
% %their relative line number within a frame
% linemods=mod(linenumbers,N);
% linemods(find(linemods==0))=N;
% %return the indices which have relative line numbers in the center of the
% %frame.. the ones we actually have predictions for
% linenumbers=find(and(linemods>edgebuffer,linemods<N-edgebuffer+1));
% 
% %coplot the real motion (with pixel resolution) with the predicted motion
% %(with line resolution)
% %not the best way to plot it, but simple enough
% %note the sections of motion that have no prediction this method of visualization
% %just draws a line between them, we have simply chosen not to use these
% %lines, though linear interpolation wouldn't be a terrible choice.
% figure(10);
% clf;
% hold on;
% plot((1:length(yofn))/size(chone,2),yofn,'r');
% plot((1:length(yofn))/size(chone,2),xofn,'k');
% plot(linenumbers,offsets');