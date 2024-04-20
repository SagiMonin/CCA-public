function params = param_stereo()
% Configurable parameters for algorithm 
params.dispRange = [-1 ,1]; % This is updated later
params.numIter = 4; % Number of iterations 

params.cost = 'BT'; % Cost function type - BT/SAD
params.gaussKerSigma = 5;% Window size for cost calculation
params.normKerSize = 7;
params.confidenceThresh = -1e-3; % Threshold for flatting parabolas 

% Border - refine
params.border_len = params.gaussKerSigma*2+max(abs(params.dispRange(:))); % Sagi - how much border is reliable 
params.penalty_border = 10; % Reduce confidence on border

% Param - penalty on edges and large-diff
params.P1param = 1;  % Weight of averaging pixels with disparity close to each other
params.P2param = 0.05;  % Weight of averaging pixels with large disparity difference
params.thr_dis = [2,2,2,2]; % If disparity is smaller than this threshold we use P1 - penalty
params.thr_dis2 = [2,2,180,180]; % If disparity is smaller than this threshold we use P2 - penalty (if it is equal to thr_dis - then use P1). Above this thrreshold don't propagate disparity but choose.
params.thr_rat = 10e2; % Ratio between confidence score of two adjacent pixels to deceide if to propagate parabola or stop.
params.sigmaEdges = 4; % Used to identify edges to give less penalty for smoothness
params.thr_edge = 5e-1; % This is used to detect edges for which we don't propagate disparity 


% Pyramid values - not used in stereo with large disparity
params.plt_pyr = 0; % Plot each level of the pyramid
params.levels = 1;% SGM pyramid levels
params.priorW = 1; % This is prior when jumping from low level to higher level in the pyramid.
params.idx_pyr = 1; % First level of pyramid - this stays 1
params.pyr_max = params.levels;

%% ENNC - parameter
params.encc_window = 5;
params.encc_window2 = params.encc_window;
%% Sub-pixel interpolation - This is based on simulations Sagi.K did to check pixel-locking
params.interpolant = 'ENCC';% ['f1', 'ENCC'] % Interpolation kernel for sub-pixel estimation
params.bias = [3.34275500790682e-05;0.00973533187061548;0.0220652986317873;0.0296091139316559;0.0247546602040529;2.39057189901359e-05;-0.0246540009975433;-0.0294469539076090;-0.0218638572841883;-0.00956899579614401;2.35768729908159e-05;0.00971676409244537;0.0220222137868404;0.0295578259974718;0.0247263666242361;3.81323552574031e-05;-0.0245587956160307;-0.0292639490216970;-0.0216563567519188;-0.00947358552366495;0.000167622143635526];
params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];
params.bias = zeros(size(params.offset));

