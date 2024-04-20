function [B] = median_disparity_filter(A)
% A -disparity map - invalid disparity should be infinity 
idx_inf = isinf(A);
fun = @(x) medifilt(x);

C = padarray(A,[1,1],'replicate','both');
B = C;
while sum(isinf(B(:)))>0
    
    B = colfilt(B,[3 3],'sliding',fun);
end
B(1,:) = []; B(:,1) = []; B(:,end) = []; B(end,:) = [];
B(~idx_inf) = A(~idx_inf);



end

function [y] = medifilt(x)

x(x==0) = inf;
num_vals  = 9;
num_inf = sum(x == inf);
x_sort = sort(x);
num_valid_vals = num_vals-num_inf;
num_valid_vals2 = ceil(num_valid_vals/2);
num_valid_vals2(num_valid_vals2==0) = 1;
sz = size(x);
row = num_valid_vals2;
col = 1:size(x,2);
yind = sub2ind(sz,row,col);

% x(isnan(x)) = [];
y = x_sort(yind);

end