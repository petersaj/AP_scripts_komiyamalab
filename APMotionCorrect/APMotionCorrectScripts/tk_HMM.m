function tk_HHM(filename, referencefilename, loadmode,maxdx,maxdy,numchannels, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Two channel example (chone=Calcium Green AM, chtwo=SR101
%This is the data from Figure 4 of the Neuron paper)
% uncomment lines 9-16
% OR uncomment lines 19-24 if you have previously calculated PI
% OR lines 27-31 if you have calculated both PI and the offsets 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %example of start to finish
% clear all
% close all
if isempty(loadmode)
loadmode = 0;
end

% filename = 'TK_090205_JF26230_beh_004_Green.tif';
% referencefilename='ALLAVG_reg_TK_090205_JF26230_beh_Green.tif';
%filename to be analyzed, path from one level up from where the script is
% filename=['data\' filename];
% referencefilename=['data\' referencefilename];
if ~exist('filename')
    [filename,pathname]=uigetfile('*.tif','pick your tif file');
else
        [pathname, filename] = fileparts(filename);
end
fullfilename=[pathname filesep filename];
if ~exist('referencefilename')
    [referencefilename,referencepathname]=uigetfile('*.tif','pick your reference file');
else
        [referencepathname, referencefilename] = fileparts(referencefilename);
end
fullreferencefilename=[referencepathname filesep referencefilename];

junk=strfind(filename,'\');
if ~isempty(junk)
    filebase=filename(junk(end)+1:end-4);
else
    filebase=filename(1:end-4);
end

if     strcmp(referencepathname, 'C:\Program Files\MATLAB\R2007b\work\HMM\')
fullsavefilename = [pathname 'matfiles\' filebase '_PI.mat'];
else
    fullsavefilename = [pathname filesep filebase '_PI.mat'];
end


switch loadmode
    case 0,
        if ~exist('maxdx')
            maxdx = 15;
        end
        if ~exist('maxdy')
            maxdy = 4;
        end
        if ~exist('numchannels')
            numchannels = 1;
        end
        tk_split_createPI(fullfilename,fullreferencefilename,fullsavefilename,maxdx,maxdy,numchannels);
        load(fullsavefilename);
        maxlambda;
        [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);
    case 1,
        if ~exist('maxdx')
            maxdx = 15;
        end
        if ~exist('maxdy')
            maxdy = 4;
        end
        if ~exist('numchannels')
            numchannels = 1;
        end
        tk_split_createPI(fullfilename,fullreferencefilename,fullsavefilename,maxdx,maxdy,numchannels);
    case 2,
        load(fullsavefilename);
        load('TK_090205_JF26230_beh_004-20G_L0.40.mat');
        [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,1);
    case 3,
        load(fullsavefilename);
        if ~exist('lambda')
            lambda = 0.050;
        end
        h = waitbar(0.0,['Run with Lambda=' num2str(lambda)]);
        set(h,'Name',['HMM Lambda=' num2str(lambda)]);
        tk_markov_on_PIsave;
        [fixeddata,countdata]=tk_playback_markov(fullfilename,offsets,edgebuffer,0);
    case 4,
        if ~exist('maxdx')
            maxdx = 15;
        end
        if ~exist('maxdy')
            maxdy = 4;
        end
        if ~exist('numchannels')
            numchannels = 1;
        end
        tk_split_createPINoAvg(fullfilename,fullreferencefilename,fullsavefilename,maxdx,maxdy,numchannels);
    case 5,
        if ~exist('maxdx')
            maxdx = 15;
        end
        if ~exist('maxdy')
            maxdy = 4;
        end
        if ~exist('numchannels')
            numchannels = 1;
        end
        tk_split_createPINoAvg(fullfilename,fullreferencefilename,fullsavefilename,maxdx,maxdy,numchannels);
        load(fullsavefilename);
        if ~exist('lambda')
            lambda = 0.050;
        end
        h = waitbar(0.0,['Run with Lambda=' num2str(lambda)]);
        set(h,'Name',['HMM Lambda=' num2str(lambda)]);
        tk_markov_on_PIsave;
        [fixeddata,countdata]=tk_playback_markov(fullfilename,offsets,edgebuffer,0);
end




