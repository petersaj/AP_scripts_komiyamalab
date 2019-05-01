%% Testing multiple ROI 

[summed_file summed_path] = uigetfile('*.tif','Select summed movies');

imageinfo = imfinfo([summed_path filesep summed_file]);
M = imageinfo(1).Height;
N = imageinfo(1).Width;
numframes = length(imageinfo);

%%
im = zeros(ceil(M),ceil(N),numframes);
for curr_frame = 1:numframes
    
    curr_im = imread([summed_path filesep summed_file], ...
        'tiff',curr_frame,'Info',imageinfo);
 
    im(:,:,curr_frame) = curr_im;
    
    disp(curr_frame);
end

border = 50;
im = im(50:end-border+1,border:end-border+1,:);
[M_new,N_new,numframes] = size(im);

%%
ds = 2;
h = fspecial('disk',2);

M_ds = ceil(M_new/ds);
N_ds = ceil(N_new/ds);

im_smoothds = zeros(M_ds,N_ds);

for curr_frame = 1:numframes
    curr_frame_smooth = conv2(im(:,:,curr_frame),h,'same');
    im_smoothds(:,:,curr_frame) = curr_frame_smooth(1:ds:end,1:ds:end);
end


%%

split_num = [3 3];

pc_kurtosis_split = cell(split_num(1),split_num(2));

M_split = floor(M_ds/split_num(1));
N_split = floor(N_ds/split_num(2));

for x_split = 1:split_num(2)
    for y_split = 1:split_num(1)
                
        [coeff score latent] = princomp(zscore(reshape( ...
            im_smoothds(1+M_split*(y_split-1):M_split*y_split, ...
            1+M_split*(x_split-1):M_split*x_split,:),[],numframes),[],2)');
        
        coeff_reshape = reshape(coeff,M_split,N_split,size(coeff,2));

        var_thresh = find(cumsum(latent)/sum(latent) > 0.95,1);
        curr_kurtosis = kurtosis(coeff_reshape(:,:,1:var_thresh),[],3);
        pc_kurtosis_split{y_split,x_split} = (curr_kurtosis- ...
            min(curr_kurtosis(:)))./max(curr_kurtosis(:)-min(curr_kurtosis(:)));
        
    end
end





%% None of that splitting is necessary: PCA is the same regardless of orientation

[summed_file summed_path] = uigetfile('*.tif','Select summed movies');

imageinfo = imfinfo([summed_path filesep summed_file]);
M = imageinfo(1).Height;
N = imageinfo(1).Width;
numframes = length(imageinfo);

im = zeros(ceil(M),ceil(N),numframes);
for curr_frame = 1:numframes
    
    curr_im = imread([summed_path filesep summed_file], ...
        'tiff',curr_frame,'Info',imageinfo);
 
    im(:,:,curr_frame) = curr_im;
    
    disp(curr_frame);
end

border = 50;
im = im(50:end-border+1,border:end-border+1,:);
[M_new,N_new,numframes] = size(im);

im = reshape(im,[],numframes);
%%
[coeff,score,latent] = princomp(zscore(zscore(im,[],1),[],2));




















