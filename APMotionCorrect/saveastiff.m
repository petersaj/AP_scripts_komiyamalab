function saveastiff(data, path, samples_per_pixel, compression_type, always_yes, no_message, image_description)
% AP 120327: added image description
%
% Examples :
% 
% [X,Y,Z] = peaks(100);
% Z_index = uint8((Z - min(Z(:))) * (255 / (max(Z(:)) - min(Z(:)))));
% 
% % 8-bit, grayscale image
% saveastiff(uint8(Z_index), 'Z_uint8.tif', 1);
% 
% % 8-bit, grayscale image (using default parameter)
% saveastiff(uint8(Z_index), 'Z_uint8.tif');
% 
% % Lossless compression
% % Compression. 1 : Uncompressed, 5 : lossless LZW, 7 : lossy JPEG, 8 : lossless Adobe-style
% saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', 1, 5);
% 
% % Overwrite with no question (default is 0 : ask question)
% saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', 1, 5, 1);
% 
% % All Messages are skipped. (default is 0 : display all message)
% saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', 1, 5, 1, 1);
% 
% % 16-bit, grayscale image, no compression, no quesion, no message
% saveastiff(uint16(Z_index), 'Z_uint16.tif', 1, 1, 1, 1);
% 
% % 32-bit single, grayscale image
% saveastiff(Z, 'Z_single.tif', 1, 1, 1, 1);
% 
% % Generate a RGB color image from the original gray scale image.
% Z_color = uint8(ind2rgb(Z_index, hsv(256)*256));
% 
% % RGB color image
% saveastiff(Z_color, 'Z_rgb.tif', 3, 1, 1, 1);
% 
% % Save each R, G and B chanels of the color image, separately.
% saveastiff(Z_color, 'Z_rgb_channel.tif', 1, 1, 1, 1);
% 
% % Generate a multi-frame RGB color image.
% Z_color_multiframe = reshape([Z_color(:)*0.2 Z_color(:)*0.6 Z_color(:)], 100, 100, 3, 3);
% 
% % Save the multi-frame RGB color image
% saveastiff(Z_color_multiframe, 'Z_rgb_multiframe.tif', 3, 1, 1, 1);
% 
% % Add white noise to the RGB color image
% Z_color_noisy = uint8(single(Z_color) + rand(100, 100, 3).*50);
% 
% % Save the noise-added RGB color image
% saveastiff(Z_color_noisy, 'Z_rgb_noisy.tif', 3, 1, 1, 1);
% 
% % 32-bit single, 50x50x50 volume data
% saveastiff(single(rand(50, 50, 50)), 'volume_50x50x50.tif', 1, 1, 1, 1);
% 
% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

errcode = 0;
try
if isreal(data) == false
    errcode = 8; assert(0);
end
if nargin < 6
    no_message = false;
end
if nargin < 5
    always_yes = false;
end

if isempty(dir(path)) == false && always_yes == false
    reply = input('File already exist, do you want continue? Y/n: ', 's');
    if isempty(reply), reply = 'Y'; end
    if any(upper(reply) ~= 'Y')
        errcode = 7; assert(0);
    end
else
    fid = fopen(path, 'w');
    if fid == -1
        errcode = 1; assert(0);
    end
    
    fclose(fid);
    delete(path);
end

if nargin < 4
    compression_type = 1;
end
if nargin < 3
    samples_per_pixel = 1;
end
if (samples_per_pixel == 1 && ndims(data) > 3) ...
    || (samples_per_pixel == 3 && ndims(data) > 4)
    errcode = 10; assert(0);
end


if isempty(data)
    errcode = 2; assert(0);
end

if samples_per_pixel == 1
    if ndims(data) >= 4, errcode = 3; assert(0); end;
    [height width depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
elseif samples_per_pixel == 3
    if ndims(data) >= 5, errcode = 3; assert(0); end;
    [height width rgb depth] = size(data);
    if rgb ~= 3, errcode = 4; assert(0); end;
    if strcmp(class(data), 'uint8'), data = uint8(data); end;
    tagstruct.Photometric = Tiff.Photometric.RGB;
else
    errcode = 5; assert(0);
end

tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.SamplesPerPixel = samples_per_pixel;
tagstruct.RowsPerStrip = height;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.ImageDescription = image_description;

% http://en.wikipedia.org/wiki/Tagged_Image_File_Format
% 1 : Uncompressed, 5 : lossless LZW, 7 : lossy JPEG, 8 : lossless Adobe-style
switch compression_type
    case {1, 5, 7, 8}
        tagstruct.Compression = compression_type;
    otherwise
        if always_yes == false
            reply = input('Unsupported compression type : Data will not be compressed.\nDo you wish to continue? Y/n: ', 's');
            if isempty(reply), reply = 'Y'; end
            if reply == 'Y' || reply =='y'
                tagstruct.Compression = 1;
            else
                errcode = 7; assert(0);
            end
        else
            tagstruct.Compression = 1;
        end
end

switch class(data)
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = 1;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = 2;
        if samples_per_pixel == 3
            errcode = 8; assert(0);
        end
    case {'single'}
        tagstruct.SampleFormat = 3;
    otherwise
        errcode = 9; assert(0);
end

switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    otherwise
        errcode = 9; assert(0);
end

tStart = tic;
tempfile = '~$temporal.tif';
t = Tiff(tempfile, 'w');
fileattrib(tempfile, '+w', '', 's');
for d = 1:depth
    t.setTag(tagstruct);
    if samples_per_pixel == 1
        t.write(data(:, :, d));
    else
        t.write(data(:, :, :, d));
    end
    if d ~= depth
        t.writeDirectory();
    end
end

t.close();
fileattrib(tempfile, '-w');
movefile(tempfile, path, 'f');

tElapsed = toc(tStart);
if no_message == false
    display(sprintf('File saved successfully. Elapsed time : %.3f s.', tElapsed));
end

catch exception
    if exist('t', 'var')
        t.close();
        delete(tempfile);
    end
    
    if no_message == false
        switch errcode
            case 1
                error 'Invalide path.';
            case 2
                error 'No image provided in data.';
            case 3
                error 'Data dimension is too large.';
            case 4
                error 'Third dimesion (color depth) is not 3.';
            case 5
                error '''samples_per_pixel'' must be 1 or 3. (1 : grayscale, 3 : rgb color)';
            case 6
                error 'int8, int16 or int32 format for rgb image is not supported.';
            case 7
                error 'File save canceled.';
            case 8
                error 'Complex number is not supported.';
            case 9
                error 'Unsupported data type.';
            case 10
                error 'Dimension of source data is too large.'
            otherwise
                rethrow(exception);
        end
    end
end
end