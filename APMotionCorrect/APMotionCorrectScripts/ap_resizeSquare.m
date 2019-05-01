function ap_resizeSquare(varargin)
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
    im_square = reshape([im(:) im(:) im(:) im(:)]',size(im,1)*4,[]);
    max_intensity = max(max(im_square));
    im_square = im_square.*(2^16/max_intensity);
    filename = ['SQ_' sourcefilename '.png'];
    imwrite(im_square,fullfile(sourcefilepath, filename),'Bitdepth',16);
end