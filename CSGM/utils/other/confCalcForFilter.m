function [conf_score] = confCalcForFilter(dispar,im_guide)
filt_size = 9;
filt = fspecial('average', filt_size); % Avg

im_guide_gray = rgb2gray(im_guide/255).*255;

%% Horizontal grad - Reduce confidence where we 
[Gx,~] = imgradientxy(im_guide_gray);
c1 = exp(-abs(Gx)/50);
%% Neig variance - Reduce confidence where we don't have texture
w_v = 100;
eps_v = 0.5;

var_guide = (imgaussfilt(im_guide_gray.^2,9) - imgaussfilt(im_guide_gray,9).^2);
c2 = exp(-max(0,w_v./var_guide - eps_v));

%% Neigh disparity
sig_u = 0.1;

d_i1 = dispar;

d_i1_shift_left = imtranslate(d_i1,[-1,0]);
d_i1_shift_right = imtranslate(d_i1,[1,0]);
in1 = min( ((d_i1-d_i1_shift_left)./sig_u.^2 ).^2 , ((d_i1-d_i1_shift_right)./sig_u.^2 ).^2);
c3 = exp(-in1);
% conf = conf.*exp(-in1);

conf_score = c1.*c2.*c3;


end

