function [imL,imR, im_guide] = loadPhone2021(pathData,imname, imgExt)
% Load data - form DSLR images, pre-processing if needed.
% imL = double(imread([pathData,imname,'_L',imgExt]));
% imR = double(imread([pathData,imname,'_R',imgExt]));
% im_guide = imL;
% imL=sum(imL,3)./3;
% imR=sum(imR,3)./3;

flag_opt = 0;

cropX = 336-20; cropY = 252-20;
normKerSize=7; % Sagi - spesific value for dual-PD (used for normalization)

% nrm=@(img,ch) ((img(:,:,ch)-min(min(img(:,:,ch))))./(max(max(img(:,:,ch)))-min(min(img(:,:,ch))))).*255;
% for ii=1:3
%     rgb1(:,:,ii)=nrm(rgb1,ii);
%     rgb2(:,:,ii)=nrm(rgb2,ii);
% end
bias = 1024; 
maxval = 2^14-1;

% Load calibraiton images
im_l_calib = double(imread([pathData,'calibration\white_sheet_left.png']));
im_r_calib = double(imread([pathData,'calibration\white_sheet_right.png']));

% Load images
imL = double(imread([pathData,imname,'_left',imgExt]));
imR = double(imread([pathData,imname,'_right',imgExt]));

% Bilateral filter substraction
diam = 15;
sigSpace = 3;
sigColor = 20;

if flag_opt == 0 % Just normalize
%     imL = sqrt(imL-bias);
%     imR = sqrt(imR-bias);
    imL = sqrt((imL-bias)/maxval)*255;
    imR = sqrt((imR-bias)/maxval)*255;
    % Normalize DP image to [0,1] (after remove bias since that is dark lvl
%     imL = (imL-bias)/maxval * 255;
%     imR = (imR-bias)/maxval * 255;
    
    
    % Crop - this is what they use
    imL = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    imR = imR(cropY:end-cropY-1,cropX:end-cropX-1);
    im_guide = imL;
    
elseif flag_opt == 1 %This is what google does
    % Normalize DP image to [0,1] (after remove bias since that is dark lvl
    imL = (imL-bias)/maxval;
    imR = (imR-bias)/maxval;
    
    % Normalize calibration  image to [0,1] (after remove bias since that is dark lvl
    im_l_calib = (im_l_calib-bias)/maxval;
    im_r_calib = (im_r_calib-bias)/maxval;
    
    % Vignneting correction
    window_size = 101;
    calib_right = imfilter(im_l_calib./im_r_calib , ones(window_size)/window_size^2);
    imR = imR.*calib_right;

    % Crop - this is what they use
    imL = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    imR = imR(cropY:end-cropY-1,cropX:end-cropX-1);

    % Scale by max
%     max_val = max([imR(:);imL(:)]);
%     imL = imL/max_val;
%     imR = imR/max_val;

    im_guide = imL;
elseif flag_opt == 2 % Try bilateral filter (according to Hirschmuller)
    % Normalize DP image to [0,1] (after remove bias since that is dark lvl
    imL = (imL-bias)/maxval;
    imR = (imR-bias)/maxval;
    imL = (imL-bias);%/maxval;
    imR = (imR-bias);%/maxval;
    
    im_guide = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    
    imL = imL - cv.bilateralFilter(imL ,'Diameter',diam, 'SigmaSpace', sigSpace, 'SigmaColor', sigColor);
    imR = imR - cv.bilateralFilter(imR ,'Diameter',diam, 'SigmaSpace', sigSpace, 'SigmaColor', sigColor);
    % Crop - this is what they use
    imL = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    imR = imR(cropY:end-cropY-1,cropX:end-cropX-1);
elseif flag_opt == 3 % Try bilateral filter (according to Hirschmuller) with vigneeting correction
    % Normalize DP image to [0,1] (after remove bias since that is dark lvl
    imL = (imL-bias)/maxval;
    imR = (imR-bias)/maxval;
    
    
    % Normalize calibration  image to [0,1] (after remove bias since that is dark lvl
    im_l_calib = (im_l_calib-bias)/maxval;
    im_r_calib = (im_r_calib-bias)/maxval;
    
    % Vignneting correction
    window_size = 101;
    calib_right = imfilter(im_l_calib./im_r_calib , ones(window_size)/window_size^2);
    imR = imR.*calib_right;
    
    
    im_guide = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    
    imL = imL - cv.bilateralFilter(imL ,'Diameter',diam, 'SigmaSpace', sigSpace, 'SigmaColor', sigColor);
    imR = imR - cv.bilateralFilter(imR ,'Diameter',diam, 'SigmaSpace', sigSpace, 'SigmaColor', sigColor);
    % Crop - this is what they use
    imL = imL(cropY:end-cropY-1,cropX:end-cropX-1);
    imR = imR(cropY:end-cropY-1,cropX:end-cropX-1);

    
end

end


