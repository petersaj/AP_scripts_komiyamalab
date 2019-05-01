function sse = Churchland2012_jPCA_expansion_fminsearch(k,X_tilde,X_dot)

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