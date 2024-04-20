%% CSGM - wrapper
% This is opitmized for DP - data set from "Modeling Defocus-Disparity in Dual-Pixel Sensors"

% Variables - debug + data location
pltFlag = 1;
dataFlag = 'DSLR';
imname = '0001';
imgExt = '.png';
pathData = 'C:\Users\local_admin\CSGM\DP_data\ICCP2020\Qualitative\';

% Parameters for algorithm 
crop_flag = 1;

if strcmp(dataFlag, 'DSLR')
    [im_L,im_R,im_guide] = loadDSLRDP(pathData,imname, imgExt);
    if crop_flag == 1
        x = 1401:2400;
        y = 1401:2400;
        im_L = im_L(y,x);
        im_R = im_R(y,x);
        im_guide = im_guide(y,x,:);
    end
    if pltFlag == 1
        figure(101); 
        ax(1) = subplot(121); imshow(im_L,[]); title('Left image');
        ax(2) = subplot(122); imshow(im_R,[]); title('Right image');
        linkaxes(ax);
    end
end

imL = im_L;
imR = imtranslate(im_L, [-2.3,0],'linear','OutputView' ,'same','FillValues',0);
% imtranslate(im_R,[-idx_dispar,0],'linear','OutputView' ,'same','FillValues',0))
wSize = [15 21];


zmFlag = 1;
swFlag = 0; % keep it disabled (no much sense for zero disparity range, too slow for mid/high resolution)
dispRange = 3;

% due to zero disp range, imR can be swapped with imL
% [D0, T0, W0] = encc(imL, imR, wSize, dispRange, zmFlag, swFlag);

[D, T, W] = enccMEX(imL, imR, wSize, dispRange, zmFlag, swFlag);

% T>1.5
% D = D-T;
% D(abs(D)>dispRange) = nan;
% D(W<0.97) = nan; % optional: keep confident values

figure;imagesc(D);colormap(bone)
D0 = D; D0(isnan(D))=0;
