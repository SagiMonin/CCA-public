function parab = subPixInterp(score_neig, params)
%% This is based on ?Design of Interpolation Functions for Sub-Pixel Accuracy Stereo-Vision Systems?,2012
% Output: parab - output parabolas

% Calculate left, middle and right scores
left_val = score_neig(:,:,1);
mid_val = score_neig(:,:,2);
right_val = score_neig(:,:,3);

% Shift values by middle values
left_val = -(left_val-mid_val); % Shift value and flip 
right_val = -(right_val-mid_val);  % Shift value and flip 
mid_val = mid_val-mid_val; % Zero middle value

% Need to remove negative values - this can happen on the border of
% disparitys
left_val(left_val<0)=0;
right_val(right_val<0)=0;

left_smaller=left_val<=right_val;

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
x = x - reshape(biases,size(x)); % Locations of minimum 
% x = zeros(size(x));

% Build parabolas with minimum at the desired place
left_val = -left_val; % Reverat values to original sign
right_val = -right_val; % Reverat values to original sign

% The parabolas are calculated as approximations
parab.a =(left_val+right_val)./2; % a - calculated with the 3 original cost values (ignoring sub-pixel refinment), so just take second derivative 
parab.a(isnan(x(:))) = -1e-50; % Not sure when nan happens
parab.a(parab.a==0) = -1e-50; % In case 3 points are at a line

parab.b = x.*(-2.*parab.a); % b - is adjusted so x is the minimum x = b/-2a -> b = -2ax
parab.c = mid_val; % Don't really need this value

