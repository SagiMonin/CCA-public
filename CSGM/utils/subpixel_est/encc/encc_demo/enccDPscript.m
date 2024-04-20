L = imread('20210705_161427_778_left.pgm');
R = imread('20210705_161427_778_right.pgm');
L = double(L)/255;
R = double(R)/255;

resLev = 1; % 1; 2; % resLev = 0 => original resolution
if resLev==0
    imL = L;
    imR = R;
    wSize = [15 21];
elseif resLev==1 % half resolution
    imL = imresize(L,.5);
    imR = imresize(R,.5);
    wSize = [11 15];
else % resLev=2 (quarter resolution)
    imL = imresize(L,.25);
    imR = imresize(R,.25);
    wSize = [5 9];
end

zmFlag = 1;
swFlag = 0; % keep it disabled (no much sense for zero disparity range, too slow for mid/high resolution)
dispRange = 0;

% due to zero disp range, imR can be swapped with imL
% [D0, T0, W0] = encc(imL, imR, wSize, dispRange, zmFlag, swFlag);

[D, T, W] = enccMEX(imL, imR, wSize, dispRange, zmFlag, swFlag);
D = D-T;
D(abs(D)>1) = nan;
D(W<0.97) = nan; % optional: keep confident values

figure;imagesc(D);colormap(bone)
D0 = D; D0(isnan(D))=0;

