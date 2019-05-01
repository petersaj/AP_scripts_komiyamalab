function [jPC_scores] = AP_jPCA(jPC_data,trial_lengths,type,handle)
% jPC_scores = AP_jPCA(jPC_data,trial_lengths,type)
%
% Perform jPCA/rotation or symmPCA/expansion analysis a la 
% Churchland & Shenoy 2012. 
%
% Input: 
%
% jPC_data - the data to do analysis on. This must have time in rows
% and cells in columns. Concatenate trials vertically
%
% trial_lengths - matrix of frame length for each trial
%
% type (optional) - choose 'rotation' for jPCA and 'expansion' for symmPCA.
% if not included, jPCA will be performed
%
% handle (optional) - handle of axis to plot on, make empty if new figure,
% don't include if don't want plot
%
% Output:
% 
% A figure plotting the data in jPC space. Color represents time (from
% blue to red).


% using the top 6 PCs

num_trials = length(trial_lengths);
if sum(trial_lengths) ~= size(jPC_data,1);
    error('Trial lengths don''t add to total ')
end

trial_frames = [0 cumsum(trial_lengths(1:end-1)) sum(trial_lengths)-1];

% Choose 'rotation' or 'expansion' for the type of jPCA/symmPCA
if nargin == 3 && ~isempty(type)
    pca_split = type;
else
    pca_split = 'rotation';
end

% perform PCA on the data, take the first six components
[coeff score latent] = princomp(jPC_data);
pcs = jPC_data*coeff(:,1:6);

% dynamics: find transformation matrix from PCs to derivative of PCs
pcs_diff = diff(pcs);
pcs_compare = pcs(1:end-1,:);

X = pcs_compare;
X_dot = pcs_diff;
M = X\X_dot;

X_tilde = blkdiag(X,X,X,X,X,X);

% find the vector k which ultimately gives M_skew
k_start = ones(15,1)*0.001;
options = optimset('MaxFunEvals',300000);

switch pca_split
    case 'rotation'
        k = fminsearch(@(x) jPCA_rotation_fminsearch ...
            (x,X_tilde,X_dot),k_start,options);
    case 'expansion'
        k = fminsearch(@(x) jPCA_expansion_fminsearch ...
            (x,X_tilde,X_dot),k_start,options);
end

% construct final M_skew
idx = ones(6);
tril_idx = logical(tril(idx,-1));
triu_idx = logical(triu(idx,1));
M_skew = zeros(6);
M_skew(tril_idx) = k;
M_skew_t = M_skew';
switch pca_split
    case 'rotation'
        M_skew_t(tril_idx) = -k;
    case 'expansion'
        M_skew_t(tril_idx) = k;
end
M_skew = M_skew_t';

% find the largest eigenvalues and pull out the top 2 (complex conj pair)
clear j
[V D] = eig(M_skew);
D = diag(D);
switch pca_split
    case 'rotation'
        jPCA_idx = find(abs(D) == max(abs(D)));
        jPCA_1 = V(:,jPCA_idx(1)) + V(:,jPCA_idx(2));
        jPCA_2 = j*(V(:,jPCA_idx(1)) - V(:,jPCA_idx(2)));
    case 'expansion'
        jPCA_idx = [D [1:length(D)]'];
        jPCA_idx = sortrows(jPCA_idx,1);
        jPCA_idx = jPCA_idx(end:-1:end-1,2);
        jPCA_1 = V(:,jPCA_idx(1)) + V(:,jPCA_idx(2));
        jPCA_2 = V(:,jPCA_idx(1)) - V(:,jPCA_idx(2));
end

% get the final jPCs, project data onto jPC space
pop_jPCs = X*[jPCA_1 jPCA_2];

if exist('handle','var')
    if isempty(handle);
        figure;
        handle = axes;
        hold on;
    end
    title('jPCs')
    for i = 1:num_trials
        curr_frames = trial_frames(i)+1:trial_frames(i+1);
        z = zeros(size(curr_frames));
        x = pop_jPCs(curr_frames,1)';
        y = pop_jPCs(curr_frames,2)';
        col = curr_frames - curr_frames(1);
        axis(handle);
        surface([x;x],[y;y],[z;z],[col;col],'facecol','no', ...
            'edgecol','interp','linew',2);
    end
end

jPC_scores = pop_jPCs;

end

function sse = jPCA_rotation_fminsearch(k,X_tilde,X_dot)

% construct Hk
idx = ones(6);
tril_idx = logical(tril(idx,-1));
triu_idx = logical(triu(idx,1));
Hk = zeros(6);
Hk(tril_idx) = k;
Hk_t = Hk';
Hk_t(tril_idx) = -k;
Hk = Hk_t';

sse = sqrt(sum((X_dot(:) - X_tilde*Hk(:)).^2));

end

function sse = jPCA_expansion_fminsearch(k,X_tilde,X_dot)

% construct Hk
idx = ones(6);
tril_idx = logical(tril(idx,-1));
triu_idx = logical(triu(idx,1));
Hk = zeros(6);
Hk(tril_idx) = k;
Hk_t = Hk';
Hk_t(tril_idx) = k;
Hk = Hk_t';

sse = sqrt(sum((X_dot(:) - X_tilde*Hk(:)).^2));

end