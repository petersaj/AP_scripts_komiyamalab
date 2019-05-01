function [x_rangenorm] = AP_rangenorm(x,dim)
% x_rangenorm = AP_rangenorm(x,dim)
%
% Normalize the range of rows or columns of a matrix from 0 to 1
% dim = dimension, if not entered assumes dim = 1

if nargin == 1
    dim = 1;
end

x_minsub = bsxfun(@minus,x,min(x,[],dim));
x_rangenorm = bsxfun(@times,x_minsub,1./max(x_minsub,[],dim));


end

