function [conf] = conf_sobel(im_guide)
h =  ones(5);
thr = 200;
rat = 100;
%% Look at gradients
xsobel = [-1 0 1;
          -2 0 2;
          -1 0 1];
  
im_guide = rgb2gray(im_guide/255)*255;
guide_filt_x = imfilter(abs(imfilter(im_guide,xsobel,'same')),h,'same');

guide_filt_x_thr = max(rat./guide_filt_x,0);
conf = exp(-guide_filt_x_thr);
end

