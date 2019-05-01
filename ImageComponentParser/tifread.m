function image = tifread(filename)
% image = tifread(filename): reads all frames from a tif file.

info = imfinfo(filename);
image = zeros(info(1).Height, info(1).Width, length(info));

for i = 1:length(info)
    image(:,:,i) = imread(filename,'index',i,'Info',info);
end
