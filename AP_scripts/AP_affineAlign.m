function corrcoef_diff = AP_affineAlign(params,im_orig,im_align,xy_ratio);
% align images x,y,theta
% fix for resolution for now - later don't hardcode?
im_orig = im_orig(repmat(1:size(im_orig,1),xy_ratio,1),:);
im_align= im_align(repmat(1:size(im_align,1),xy_ratio,1),:);

% define parameters
dx = params(1);
dy = params(2);
theta = params(3);

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
 
tform_rotate = maketform('affine',A);

% transform image
im_align_tform = imtransform(im_align,tform_rotate,...
    'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);

% get 1-correlation of image (function to minimize, ideally 0)
corrcoef_tform = corrcoef(im_orig,im_align_tform);
corrcoef_diff = 1-corrcoef_tform(2);
