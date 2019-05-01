function tk_RemoveRed(varargin)
if isempty(varargin)
    [sourcefilename, sourcefilepath]  = uigetfile('Pick a source file');
    sourcefile = fullfile(sourcefilepath, sourcefilename);
else
    sourcefile = varargin{1};
    [sourcefilepath, sourcefilename] = fileparts(sourcefile);
end

info=imfinfo(sourcefile);

for i=1:2:size(info,1)
    im=imread(sourcefile,i);
    filename = [sourcefilename,'_Green.tif'];
    imwrite(im,fullfile(sourcefilepath, filename),'tif','Compression','none','WriteMode','append');
end
