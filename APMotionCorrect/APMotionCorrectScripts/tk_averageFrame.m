function tk_averageFrame(varargin)
if isempty(varargin)
    [sourcefilename, sourcefilepath]  = uigetfile('*.tif','Pick an image');
    sourcefile = fullfile(sourcefilepath, sourcefilename);
else
    sourcefile = varargin{1};
    [sourcefilepath, sourcefilename] = fileparts(sourcefile);
end
[PATHSTR,NAME,EXT,VERSN] = fileparts(sourcefile);
info=imfinfo(sourcefile);
cdata = [];
for i = 1 : size(info,1)
    if isempty(cdata)
        cdata = double(imread(sourcefile,i));
    else
        cdata = cdata + double(imread(sourcefile,i));
    end
end
cdata = uint16(cdata / size(info,2));
    filename = ['AVG_',sourcefilename,'.tif'];
imwrite(cdata,fullfile(sourcefilepath, filename),'tif','Compression','none','WriteMode','append');