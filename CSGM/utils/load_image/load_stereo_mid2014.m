function [im_L,im_R,im_guide_L,im_guide_R,GT_X_L,GT_X_R,GT_Y_L,GT_Y_R] = load_stereo_mid2014(path,reduce_fact,y_flag)
% Load data of 2014 from stereo.

GT_Y_L = 0; GT_Y_R = 0;
% Load y
if y_flag == 1
    GT_Y_L = readpfm([path,'\disp0y.pfm']);   
    GT_Y_L = imresize(GT_Y_L, 1/reduce_fact)/(reduce_fact);
    
    GT_Y_R = readpfm([path,'\disp1y.pfm']);
    GT_Y_R = imresize(GT_Y_R, 1/reduce_fact)/(reduce_fact);
end

% Load x
GT_X_L = readpfm([path,'\disp0.pfm']);
GT_X_L = imresize(GT_X_L, 1/reduce_fact)/(reduce_fact);

GT_X_R = readpfm([path,'\disp1.pfm']);
GT_X_R = imresize(GT_X_R, 1/reduce_fact)/(reduce_fact);

% Load image L+R
im_L = double(imread([path,'\im0.png']))/255;
im_R = double(imread([path,'\im1.png']))/255;

lvl_pyr = log2(reduce_fact)+1;
% [im_blur_L, im_diff_L, im_pyr_L] = pyrImg(im_L, lvl_pyr); 
% [im_blur_R, im_diff_R, im_pyr_R] = pyrImg(im_R, lvl_pyr);
% im_L = im_pyr_L{1};
% im_R = im_pyr_R{1};
im_L = imresize(im_L, 1/reduce_fact);
im_R = imresize(im_R, 1/reduce_fact);

im_guide_L = im_L*255;
im_guide_R = im_R*255;


im_L = rgb2gray(im_L)*255;
im_R = rgb2gray(im_R)*255;






end

