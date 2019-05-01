function [cor] = tk_imCorr(varargin)
if isempty(varargin)
    [sourcefilename, sourcefilepath]  = uigetfile('*.tif','Pick source image','E:\ImagingData\');
    sourcefile = fullfile(sourcefilepath, sourcefilename);
    [targetfilename, targetfilepath]  = uigetfile('*.tif','Pick target image','E:\ImagingData\');
    targetfile = fullfile(targetfilepath, targetfilename);
else
    sourcefile = varargin{1};
    [sourcefilepath, sourcefilename] = fileparts(sourcefile);
    targetfile = varargin{2};
    [targetfilepath, targetfilename] = fileparts(targetfile);
end
[PATHSTR,NAME,EXT,VERSN] = fileparts(sourcefile);

if size(imfinfo(sourcefile),2)~=1
    error('sourcefile needs to be a single frame')
end

sourceim = double(imread(sourcefile));


cor = [];
for i = 1:size(imfinfo(targetfile),2)
    targetim = double(imread(targetfile,i));
    correlation = corrcoef(sourceim, targetim);
    cor=[cor, correlation(2)];
end
  