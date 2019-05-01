function p = AP_ranksum_matrix(x,y)
% p = AP_ranksum_matrix(x,y);
% Matlab ranksum function, but on matricies
% rows are observations, columns are dimensions

if size(x) ~= size(y)
    error('Sizes of x,y don''t match')
end

p = nan(1,size(x,2));
for i = 1:size(x,2)
    if all(isnan(x(:,i))) || all(isnan(y(:,i)))
        continue
    end
    p(i) = ranksum(x(:,i),y(:,i));

end