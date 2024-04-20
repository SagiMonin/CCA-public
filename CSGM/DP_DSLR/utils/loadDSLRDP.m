function [imL,imR, im_guide] = loadDSLRDP(pathData,imname, imgExt)
% Load data - form DSLR images, pre-processing if needed.
% imL = double(imread([pathData,imname,'_L',imgExt]));
% imR = double(imread([pathData,imname,'_R',imgExt]));
% im_guide = imL;
% imL=sum(imL,3)./3;
% imR=sum(imR,3)./3;

% This works better
imL = im2double(imread([pathData,imname,'_L',imgExt]));
imR = im2double(imread([pathData,imname,'_R',imgExt]));


im_guide = imL*255;
imL=rgb2gray(imL)*255;
imR=rgb2gray(imR)*255;


end

