function [imL,imR, im_guide] = loadDSLRDP_Qual(pathData,imname, imgExt)
% Load data - form DSLR images, pre-processing if needed.
% imL = double(imread([pathData,imname,'_L',imgExt]));
% imR = double(imread([pathData,imname,'_R',imgExt]));
% im_guide = imL;
% imL=sum(imL,3)./3;
% imR=sum(imR,3)./3;

% This works better
imL = im2double(imread([pathData,imname,'_L',imgExt]));
imR = im2double(imread([pathData,imname,'_R',imgExt]));
clip = 100;

im_guide = imL(clip:end-clip+1,clip:end-clip+1,:)*255;
imL=rgb2gray(imL(clip:end-clip+1,clip:end-clip+1,:))*255;
imR=rgb2gray(imR(clip:end-clip+1,clip:end-clip+1,:))*255;


end

