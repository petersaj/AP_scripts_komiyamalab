function corrcoef_diff = AP_affineAlign(params,im_orig,im_align_z);
% align images x,y,z

% define parameters
dx = params(1);
dy = params(2);
theta = params(3);

% fix for resolution for now - later don't hardcode?
im_orig; % this is the original image (one slice)
im_align_z; % run through all slices of stack to align

% center image and constrain output to image size
[N M] = size(im_orig);

udata = [1 N] - (N+1)/2;% input image x
vdata = [1 M] - (M+1)/2;% input image y

xdata = udata;% output image x
ydata = vdata;% output image y

% define affine transform matrix
A = ...
    [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
     dx          dy          1];
 
tform_translate = maketform('affine',A);

% transform image
im_align_tform = imtransform(im_align_z,tform_translate,...
    'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);

% get 1-correlation of image (function to minimize, ideally 0)
corrcoef_tform = corrcoef(im_orig,im_align_tform);
corrcoef_diff = 1-corrcoef_tform(2);
