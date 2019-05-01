% Align stack to single slice

% pick original image (single slice)
[orig_file orig_path] = uigetfile('*.tif','Pick original image');

% pick stack image (multiple slices)
[stack_file stack_path] = uigetfile('*.tif','Pick stack to align');

% load in original image
im_orig = double(imread([orig_path orig_file],'tiff'));

% load in stack
stack_info = imfinfo([stack_path stack_file],'tiff');
M=stack_info(1).Width;
N=stack_info(1).Height;
stack_num = length(stack_info);
im_stack = zeros(N,M,stack_num);
for i = 1:stack_num
    im_stack(:,:,i) = double(imread([stack_path stack_file],'tiff',i));
end    

% input parameters and minimize 1-corrcoef
align_fun = zeros(stack_num,1);
estimates = zeros(stack_num,2);

tic
w = waitbar(0,'Aligning stack to original')
for j = 1:stack_num
    im_align_z = im_stack(:,:,j);
    starting = [0 0 0]; % dx dy theta
    options = optimset('Algorithm','interior-point');
    lb = [-30 -30 -20];
    ub = [30 30 20];
    [estimates(j,:) align_fun(j)] = fmincon(@(x) AP_3DaffineAlign(x,im_orig,im_align_z),starting,[],[],[],[],lb,ub,[],options);
    waitbar(j/stack_num,w,'Aligning stack to original');
end
close(w);
toc
