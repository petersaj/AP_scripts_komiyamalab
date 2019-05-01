function Churchland2012_jPCA_fminsearch

idx = ones(6);
tril_idx = logical(tril(idx,-1));
triu_idx = logical(triu(idx,1));
% construct m
m = zeros(6);
m(tril_idx) = k;
m_t = m';
m_t(tril_idx) = -k;
m = m_t(:)';
