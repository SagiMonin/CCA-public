% Choose data set - parameters don't change
imgset = 'training';
imgsize = 'Q';



%% Parameters to choose 
params = param_stereo; % This holds all relevant parameters
im_num = 1; % Numeber of image to use 


%% Load data + adjust images
[im_L,im_R,GT,mask,ndisp] = load_midbury(im_num,imgsize,imgset);

im_L = double(im_L);
im_R = double(im_R);

im_guide_L = im_L;
im_guide_R = im_R;

im_L = rgb2gray(im_L/255)*255;
im_R = rgb2gray(im_R/255)*255;

im_L_orig = im_L;
im_R_orig = im_R;

% Adjust the range of disparities to the chosen resolution
if imgsize == 'Q'
    DisparityRange = [-round(ndisp(im_num)/4),-1];
elseif imgsize == 'H'
    DisparityRange = [-round(ndisp(im_num)/2),-1];
else
    DisparityRange = [-round(ndisp(im_num)),-1];
end
%     
 im_L = im_L - cv.bilateralFilter(im_L ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
 im_R = im_R - cv.bilateralFilter(im_R ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);

%% Calculate integer disparity 

% Left image disparity 
params.dispar_vals = DisparityRange(1)-1:DisparityRange(2)+1;
[cost_neig, conf_score, dispar_int_val] = disparCostAltConf(im_L, im_R, params, im_guide_L);

% Right image disparity 
params.dispar_vals = -fliplr(params.dispar_vals); 
[cost_neig_R, conf_score_R, dispar_int_val_R] = disparCostAltConf(im_R,im_L, params, im_guide_R);


%% C-SGM 
% Refine Condidence scores - Speckle reduction
val_reduce = 1e-10;
if 1
    speckle_filter = cv.filterSpeckles(dispar_int_val, 1000, 128, 8);
    idxs_filterLR = (speckle_filter == 1000);
    conf_score(idxs_filterLR) = val_reduce;

    speckle_filter = cv.filterSpeckles(dispar_int_val_R, 1000, 128, 8);
    idxs_filterRL = (speckle_filter == 1000);
    conf_score_R(idxs_filterRL) = val_reduce;
end

% Refine Condidence scores - L-R check

if 1
    LeftRightThresh = 1;
    [LR_filt,RL_filt]=leftRightConsistency(dispar_int_val,dispar_int_val_R,LeftRightThresh);
    idxs_filter = isnan(LR_filt);
    LR_filt(isnan(LR_filt)) = inf;
    LR_filt(idxs_filterLR) = inf;
    [LR_filt] = median_disparity_filter(LR_filt);
    conf_score(idxs_filter) = val_reduce;

    idxs_filter = isnan(RL_filt);
    conf_score_R(idxs_filter) = val_reduce;
    RL_filt(isnan(RL_filt)) = inf;
    RL_filt(idxs_filterRL) = inf;
    [RL_filt] = median_disparity_filter(RL_filt);
end

% Run CSGM
[dispar_map, sum_parab] = wrapperCSGMStereo(im_L_orig, im_R_orig, im_guide_L, params ,cost_neig, double(conf_score), dispar_int_val);
[dispar_map_R, sum_parab_R] = wrapperCSGMStereo(im_R_orig, im_L_orig, im_guide_R, params, cost_neig_R, double(conf_score_R), dispar_int_val_R );

% Calculate raw score before filter
[good05_CSGM,good1_CSGM,good2_CSGM,good4_CSGM,rmse_CSGM] = metrics_stereo(GT,-dispar_map,mask,0,0);


% Prepare post - filter: L-R check + speckle reduction                      
dispar_LR = double(dispar_map);
dispar_RL = double(dispar_map_R);

% L-R check
LeftRightThresh=1.0;
[dispar_LR_filter,~]=leftRightConsistency(dispar_LR,dispar_RL,LeftRightThresh);
dispar_LR_filter(isnan(dispar_LR_filter)) = 0;

% Speckle check
dispar_LR_filter = double(cv.filterSpeckles(double(dispar_LR_filter*16), 0, 140, 8))/16;

masknan = ~(dispar_LR_filter==0);
mask2 = masknan & mask;
                   
% Create files to save - confidence, disparity and guide image
conf_score2 = double(mask2);
conf_score = conf_score2;
conf_score(conf_score<=0) = 1e-150;

% Fill values with median filter
dispar_LR_filter5 = -dispar_LR_filter;
dispar_LR_filter5(dispar_LR_filter5==0) = inf;
[curr_dispar] = median_disparity_filter(dispar_LR_filter5);


im_guide = uint8(im_guide_L);
path_save = 'C:\Users\sagim\OneDrive - Technion\CCA_PUBLISH\demo_Stereo\';
save([path_save,'data_im',num2str(im_num),'.mat'], 'im_guide','curr_dispar', 'conf_score');


%% Load result and show graph

output_solver = load('result_filter1.mat', 'output_solver'); output_solver = output_solver.output_solver;
output_solver = medfilt2(output_solver,[3,3]);

[good05_filter,good1_filter,good2_filter,good4_filter,rmse_filter] = metrics_stereo(GT,output_solver,mask,0,0);

figure(5);
subplot(221); imagesc(im_guide); title('Guide image');
subplot(222); imagesc(GT); title('gt');
subplot(223); imagesc(-dispar_map); title('C-SGM result');
subplot(224); imagesc(output_solver); title('C-SGM + filter');