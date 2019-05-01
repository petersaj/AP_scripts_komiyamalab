function tk_HHMLoop(filedir, referencefilename, loadmode,maxdx,maxdy,numchannels, lambda)
if ~exist ('filedir')
    filedir  = uigetdir('C:\Program Files\MATLAB\R2007b\work\HMM','Pick a source directory');
else
    if ~isdir(filedir)
        error('Second argument must be a valid directory');
    end
end
if ~exist('referencefilename')
    [referencefilename,referencepathname]=uigetfile('C:\Program Files\MATLAB\R2007b\work\HMM','pick your reference file','*.tif');
else
        [referencepathname, referencefilename] = fileparts(referencefilename);
end
fullreferencefilename=[referencepathname referencefilename];
if ~exist('loadmode')
loadmode = 3;
end
if ~exist('maxdx')
    maxdx = 15;
end
if ~exist('maxdy')
    maxdy = 4;
end
if ~exist('numchannels')
    numchannels = 1;
end
if ~exist('lambda')
    lambda = 0.050;
end

files = dir(filedir);
fileNames = {files.name};

for i=3:length(files)
    [PATHSTR,NAME,EXT,VERSN] = fileparts(files(i).name);
    if isequal(EXT, '.tif')
        file = [filedir filesep files(i).name];
        tk_HMM(file, fullreferencefilename, loadmode,maxdx,maxdy,numchannels, lambda);
    end
end
