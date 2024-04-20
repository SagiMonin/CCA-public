% Choose data set - parameters don't change
imgset = 'training';
imgsize = 'Q';


%% Parameters to choose 
im_num = 1; % Numeber of image to use 


% Block size for raw cost aggergation. simple block matching
params.block_size = 5; 
% Penalise disparity different than 1 for neighbor 
params.P1 = (1/16)*8*params.block_size.^2; 
% Penalise for more than 1 disparity different from the neighbor
params.P2 = (1/16)*32*params.block_size.^2; 
% Number of SGM path direction
params.directions_num = 8; 
params.gaussKerSigma = params.block_size; % Window size
params.cost = 'BT';
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
min_disp = 0;
params.disparity_range = [min_disp ceil(ndisp(im_num)+1)]; 


%     
 im_L = im_L - cv.bilateralFilter(im_L ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);
 im_R = im_R - cv.bilateralFilter(im_R ,'Diameter',15, 'SigmaSpace', 3, 'SigmaColor', 20);

%% Calculate integer disparity 

% Left image disparity 
params.dispar_vals = DisparityRange(1)-1:DisparityRange(2)+1;
[~, ~, dispar_int_val] = disparCostAltConf(im_L, im_R, params, im_guide_L);

% Right image disparity 
params.dispar_vals = -fliplr(params.dispar_vals); 
[~, ~, dispar_int_val_R] = disparCostAltConf(im_R,im_L, params, im_guide_R);

%% SGM 
% Run SGM
[ disparity_map,Cost_SGM ] = SGMWrapper( im_L, im_R, params );
Cost_SGM_L = fliplr(Cost_SGM);
disparity_map_L = fliplr(disparity_map);

[ disparity_map_R,Cost_SGM_R ] = SGMWrapper( fliplr(im_R), fliplr(im_L), params );

% Sub-pixel - refinment - L
max_idx_L = disparity_map_L+1;
[xs,ys]=meshgrid(1:size(Cost_SGM_L,2),1:size(Cost_SGM_L,1));
cost_neig_L = zeros(size(Cost_SGM_L,1),size(Cost_SGM_L,2),3);
cost_neig_L(:,:,1) = Cost_SGM_L(sub2ind(size(Cost_SGM_L),ys,xs,max_idx_L(:,:,1)-1));
cost_neig_L(:,:,2) = Cost_SGM_L(sub2ind(size(Cost_SGM_L),ys,xs,max_idx_L(:,:,1)));
cost_neig_L(:,:,3) = Cost_SGM_L(sub2ind(size(Cost_SGM_L),ys,xs,max_idx_L(:,:,1)+1));

params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];
params.bias = zeros(21,1);
params.interpolant = 'parabola';
parab_L = subPixInterp(-cost_neig_L, params);

x_L = -parab_L.b./(2*parab_L.a);
x_L(isnan(x_L)) = 0;
disparity_map_L = disparity_map_L+x_L;

% Sub-pixel - refinment - R
max_idx_R = disparity_map_R+1;
[xs,ys]=meshgrid(1:size(Cost_SGM_R,2),1:size(Cost_SGM_R,1));
cost_neig_R = zeros(size(Cost_SGM_R,1),size(Cost_SGM_R,2),3);
cost_neig_R(:,:,1) = Cost_SGM_R(sub2ind(size(Cost_SGM_R),ys,xs,max_idx_R(:,:,1)-1));
cost_neig_R(:,:,2) = Cost_SGM_R(sub2ind(size(Cost_SGM_R),ys,xs,max_idx_R(:,:,1)));
cost_neig_R(:,:,3) = Cost_SGM_R(sub2ind(size(Cost_SGM_R),ys,xs,max_idx_R(:,:,1)+1));

params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];
params.bias = zeros(21,1);
params.interpolant = 'parabola';
parab_R = subPixInterp(-cost_neig_R, params);

x_R = -parab_R.b./(2*parab_R.a);
x_R(isnan(x_R)) = 0;
disparity_map_R = disparity_map_R+x_R;


% Calculate raw score before filter
[good05_SGM,good1_SGM,good2_SGM,good4_SGM,rmse_SGM] = metrics_stereo(GT,disparity_map_L,mask,0,0);

% Prepare post - filter: L-R check + speckle reduction                      
dispar_LR = -double(disparity_map_L);
dispar_RL = double(disparity_map_R);

% L-R check
LeftRightThresh=1.0;
[dispar_LR_filter,~]=leftRightConsistency(dispar_LR,dispar_RL,LeftRightThresh);
dispar_LR_filter(isnan(dispar_LR_filter)) = 0;

% Speckle check
dispar_LR_filter2 = double(cv.filterSpeckles(double(-dispar_LR_filter*16), 0, 140, 8))/16;

mask2 = dispar_LR_filter2>0;

                   
% Create files to save - confidence, disparity and guide image
conf_score = double(mask2);
conf_score(conf_score==0) = 1e-150;

% Fill values with median filter
dispar_LR_filter5 = dispar_LR_filter2;
dispar_LR_filter5(dispar_LR_filter5==0) = inf;
[curr_dispar] = median_disparity_filter(dispar_LR_filter5);


im_guide = uint8(im_guide_L);
path_save = 'C:\Users\local_admin\Documents\GitHub\C-SGM\demo_Stereo\';
save([path_save,'data_im_sgm',num2str(im_num),'.mat'], 'im_guide','curr_dispar', 'conf_score');


%% Load result and show graph

output_solver = load('result_filter_SGM1.mat', 'output_solver'); output_solver = output_solver.output_solver;
output_solver = medfilt2(output_solver,[3,3]);

[good05_filter,good1_filter,good2_filter,good4_filter,rmse_filter] = metrics_stereo(GT,output_solver,mask,0,0);

figure(5);
subplot(221); imagesc(im_guide); title('Guide image');
subplot(222); imagesc(GT); title('gt');
subplot(223); imagesc(-dispar_map); title('C-SGM result');
subplot(224); imagesc(output_solver); title('C-SGM + filter');