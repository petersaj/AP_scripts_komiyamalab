function p = AP_signrank_matrix(x,y)
% p = AP_signrank_matrix(x,y);
% Matlab signrank function, but on matricies
% rows are observations, columns are dimensions

if nargin == 1
    
    p = nan(1,size(x,2));
    for i = 1:size(x,2)
        p(i) = signrank(x(:,i));
    end
    
elseif nargin == 2
    
    if size(x) ~= size(y)
        error('Sizes of x,y don''t match')
    end
    
    p = nan(1,size(x,2));
    for i = 1:size(x,2)
        if all(isnan(x(:,i))) || all(isnan(y(:,i)))
            continue
        end
        p(i) = signrank(x(:,i),y(:,i));
    end
    
end


