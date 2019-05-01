function oimg = loadtiff(path)
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

s = warning('off', 'all'); % To ignore unknown TIFF tag.

% Frame number
tiff = Tiff(path, 'r');
n = tiff.getTag('ImageWidth');
m = tiff.getTag('ImageLength');
p = 1;
while ~tiff.lastDirectory()
    p = p + 1;
    tiff.nextDirectory();
end

% Grayscale or color
samples_per_pixel = tiff.getTag('SamplesPerPixel');

% Pixel format
switch tiff.getTag('SampleFormat')
    case 1
        switch tiff.getTag('BitsPerSample')
            case 8
                data_type = 'uint8';
            case 16
                data_type = 'uint16';
            case 32
                data_type = 'uint32';
        end
    case 2
        switch tiff.getTag('BitsPerSample')
            case 8
                data_type = 'int8';
            case 16
                data_type = 'int16';
            case 32
                data_type = 'int32';
        end
    case 3
        data_type = 'single';
end


if samples_per_pixel == 1
    oimg = zeros(m, n, p, data_type); % grayscle
else
    oimg = zeros(m, n, samples_per_pixel, p, data_type); % color
end

% read image data and save to 'oimg'.
for k=1:p
    tiff.setDirectory(k);
    
    if samples_per_pixel == 1
        oimg(:, :, k) = tiff.read();
    else
        oimg(:, :, :, k) = tiff.read();
    end
end
tiff.close();

warning(s);
end
