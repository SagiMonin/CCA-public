tmpdisparity_map = disparity_map+1;
tmpdisparity_map(tmpdisparity_map==74) = 73;
% tmpdisparity_map(tmpdisparity_map==73) = 72;
tmpdisparity_map(tmpdisparity_map==1) = 2;
% tmpdisparity_map(tmpdisparity_map==0) = 1;

% Only need to pass costs of max+2 neighboors
[xs,ys]=meshgrid(1:size(Cost_SGM,2),1:size(Cost_SGM,1));
cost_neig = zeros(size(Cost_SGM,1),size(Cost_SGM,2),3);
cost_neig(:,:,1) = Cost_SGM(sub2ind(size(Cost_SGM),ys,xs,tmpdisparity_map-1));
cost_neig(:,:,2) = Cost_SGM(sub2ind(size(Cost_SGM),ys,xs,tmpdisparity_map));
cost_neig(:,:,3) = Cost_SGM(sub2ind(size(Cost_SGM),ys,xs,tmpdisparity_map+1));


% Calculate left, middle and right scores
left_val = cost_neig(:,:,1);
mid_val = cost_neig(:,:,2);
right_val = cost_neig(:,:,3);

% Shift values by middle values
left_val = -(left_val-mid_val); % Shift value and flip 
right_val = -(right_val-mid_val);  % Shift value and flip 
mid_val = mid_val-mid_val; % Zero middle value

% Need to remove negative values - this can happen on the border of
% disparitys
left_val(left_val<0)=0;
right_val(right_val<0)=0;

left_smaller=left_val<=right_val;
params.interpolant = 'parabola';
% Choose interpolaton method of data
switch params.interpolant
    case 'parabola'
        % Parabola model
        f=@(x_) x_./(x_+1);
    case 'f1'
        % from: "Design of Interpolation Functions for Sub-Pixel Accuracy Stereo-Vision Systems"
        f=@(x_) 0.25.*(x_+x_.^2);
    case 'f2'
        % from: "New Sub-Pixel Interpolation Function for Accurate Real-Time Stereo-Matching Algorithms"
        f=@(x_) 0.5.*(sin(x_.*pi/2-pi/2)+1);
    case 'f3'
        % from: "New Sub-Pixel Interpolation Function for Accurate Real-Time Stereo-Matching Algorithms"
        f=@(x_) 0.25.*(x_+x_.^4);
    case 'f4'
        % from: "New Sub-Pixel Interpolation Function for Accurate Real-Time Stereo-Matching Algorithms"
        f=@(x_) max(0.25.*(x_+x_.^4),1-cos(x_.*pi/3));
    case 'f5'
        % from: "Design of Interpolation Functions for Sub-Pixel Accuracy Stereo-Vision Systems"
        f=@(x_) 0.5-0.5.*cos(x_*pi/2);
end

% Define functions 
d_final_1 = @(leftVal_,rightVal_) -0.5+f(leftVal_./rightVal_);
d_final_2 = @(leftVal_,rightVal_) 0.5-f(rightVal_./leftVal_);

% Calculate values x values - interpolation of data to sub-pixel values
x1 = d_final_1(left_val,right_val);
x2 = d_final_2(left_val,right_val);
x = x2;
x(left_smaller) = x1(left_smaller);
x(left_val==0 & right_val==0)=0; % avoid nans

% Shift data according to offset+biases pre-calibrated.
biases = interp1(params.offset, params.bias, x(:));
% x = x - reshape(biases,size(x)); % Locations of minimum 


[good05,good1,good2,good4,rmse] = metrics_stereo(GT,disparity_map+x,mask,num_disp,0);
display([good05, good1,good2,good4 rmse])

[good05,good1,good2,good4,rmse] = metrics_stereo(GT,disparity_map-x,mask,num_disp,0);
display([good05, good1,good2,good4 rmse])
