function [parab, invalid_mask] = estSubPixENCC(score_neig, dispar_int_val, im_L, im_R, params)
%% This is based on ?Design of Interpolation Functions for Sub-Pixel Accuracy Stereo-Vision Systems?,2012
% Output: parab - output parabolas
left_val = score_neig(:,:,1);
mid_val = score_neig(:,:,2);
right_val = score_neig(:,:,3);

right_larger=left_val<=right_val;

%%
%Need to shift disparity values to be positive.
% min_dispar_val = min(dispar_int_val(:));
% dispar_int_val = dispar_int_val - min_dispar_val;

% Need to shift image 
% im_R = imtranslate(im_R,[-min_dispar_val,0],'linear','OutputView' ,'same','FillValues',0);
% 
% tic;
% [tau_clean, invalid_mask] = subPixelENCC(im_L, im_R, dispar_int_val, params, right_larger);
% toc 
% 
% dispar_int_val=dispar_int_val+min_dispar_val;
% tau_clean = (tau_clean+1).*double(tau_clean<-0.5)+tau_clean.*double(tau_clean>=-0.5); 
% x = tau_clean;
% 
% est = tau_clean.*(right_larger<0.5)-tau_clean.*(right_larger>0.5)+dispar_int_val;
% x = tau_clean.*(right_larger<0.5)-tau_clean.*(right_larger>0.5);
% x(x<-0.5) = -0.5; x(x>0.5) = 0.5;
%%+
%% This is how to compensate if the disparity is positive
max_disp = max(dispar_int_val(:))+1;
dispar_map = dispar_int_val - max_disp;
im_R_shift = imtranslate(im_R, [-max_disp,0],'linear','OutputView' ,'same','FillValues',0);
% tic
[ T3, W3,max_map] = encc_subpixel_alter2(im_L, im_R_shift, [params.encc_window,params.encc_window], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);
% [ T3, W3,max_map] = encc_subpixel_alter3(im_L, im_R_shift, [params.encc_window2,params.encc_window], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);
% toc
% x = double(max_map==2).*T3+double(max_map==3);%.*(T3+1)+double(max_map==1).*(T3-1);


% tic
% [ T4,max_map2] = encc_subpixel_alter3_2(im_L, im_R_shift, [params.encc_window2,params.encc_window], max(abs(dispar_map(:)))+1, 1, 0, dispar_map);

% [ x] = encc_subpixel_alter2_parallal(im_L, im_R_shift, [params.encc_window2,params.encc_window], max(abs(dispar_map(:)))+1, 0, dispar_map);
% toc
x = double(max_map==2).*T3+double(max_map==3).*(T3+1);%+double(max_map==1).*(T3-1);
% x = double(max_map2==1).*T4+double(max_map2==2).*(T4+1);%+double(max_map==1).*(T3-1);
% x = (T3+1);%+double(max_map==1).*(T3-1);

invalid_mask = zeros(size(x));
invalid_mask = isnan(x);
% invalid_mask(max_map~=2 | max_map~=3) = 1;
x(isnan(x)) = 0;
% final_dispar = dispar_map;
% final_dispar = zeros(size(max_map));
% final_dispar = dispar_map+double(max_map==2).*T3+double(max_map==3).*(T3+1) + max_disp;

%%
%


% Shift values by middle values
left_val = -(left_val-mid_val); % Shift value and flip 
right_val = -(right_val-mid_val);  % Shift value and flip 
mid_val = mid_val-mid_val; % Zero middle value

% Need to remove negative values - this can happen on the border of
% disparitys
left_val(left_val<0)=0;
right_val(right_val<0)=0;

x(left_val==0 & right_val==0)=0; % avoid nans
% Build parabolas with minimum at the desired place
left_val = -left_val; % Reverat values to original sign
right_val = -right_val; % Reverat values to original sign


% % The parabolas are calculated as approximations
parab.a =(left_val+right_val)./2; % a - calculated with the 3 original cost values (ignoring sub-pixel refinment), so just take second derivative 
parab.b = x.*(-2.*parab.a); % b - is adjusted so x is the minimum x = b/-2a -> b = -2ax
parab.c = mid_val; % Don't really need this value
parab.a(isnan(x(:))) = 0; % Not sure when nan happens

% 
% % Shift parabaols by disparaity values so all parabolas are centered around zero.
aNew = parab.a;
bNew = parab.b-2*parab.a.*dispar_int_val;
cNew = parab.a.*dispar_int_val.^2-parab.b.*dispar_int_val+parab.c;
parab.a = aNew;
parab.b = bNew;
parab.c = cNew;

