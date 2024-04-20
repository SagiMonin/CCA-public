clc;
clear all;

addpath('C:\source\depthPlatformGIT\Matlab\opencvMatlab\mexopencv')
addpath('C:\source\depthPlatformGIT\Matlab\opencvMatlab\mexopencv\opencv_contrib')

img = rgb2gray(imread('img.png'));
depthImg=  imread('depthImg.png');


fastSmoother = cv.FastGlobalSmootherFilter(img,'Lambda',30,'SigmaColor',7.5,'LambdaAttenuation',0.19,'NumIter',2);
depthImgFiltered = fastSmoother.filter(depthImg);

figure(1)
imshow(depthImg,[]);impixelinfo;
figure(2)
imshow(depthImgFiltered,[]);impixelinfo;

