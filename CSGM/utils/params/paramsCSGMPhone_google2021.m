function params = paramsCSGMPhone_google2021()
% Configurable parameters for algorithm 
flagLRC = 1; % LRC processing - for phone-cameras
flagPost = 0; % Post processing 

params.dispRange = [-4 ,4]; % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
params.numIter = [6,2,1,1,1] ; % Number of iterations 

% params.subpix = ''; 'ENCC'; 

params.cost = 'SAD'; % Cost function type
params.gaussKerSigma = 3;% STD of gaussian for average-weighting of score
params.normKerSize = 7; % Used for image normalization

params.confidenceThresh=-12e-2; % Threshold for flatting parabolas %%% Tune-flag

% Border - refine
params.border_len = params.gaussKerSigma*2+max(abs(params.dispRange(:))); % Sagi - how much border is reliable 
params.penalty_border = 1000; % Reduce confidence on border

% Param - penalty on edges and large-diff
params.P1param = 5;  % Penalty for large diff in dispairty
params.sigmaEdges = 2; % Used to identify edges to give less penalty for smoothness
% params.currentToPrevWeight=0.2; % Not-used % SAgi - look at code % A large number means that we consider more the confidence of the previous sample so a "noiser" result with be created

% Pyramid values
params.plt_pyr = 0; % Plot each level of the pyramid
params.levels = [4];% SGM pyramid levels
params.priorW = 1e-2; % This is prior when jumping from low level to higher level in the pyramid.
params.idx_pyr = 1; % First level of pyramid - this stays 1
params.pyr_max = 3;
%% ENNC - parameter
params.encc_window = [15,15];

%% Sub-pixel interpolation - This is based on simulations Sagi.K did to check pixel-locking
params.interpolant = 'f1';% ['f1', 'ENCC'] % Interpolation kernel for sub-pixel estimation
params.bias = [3.34275500790682e-05;0.00973533187061548;0.0220652986317873;0.0296091139316559;0.0247546602040529;2.39057189901359e-05;-0.0246540009975433;-0.0294469539076090;-0.0218638572841883;-0.00956899579614401;2.35768729908159e-05;0.00971676409244537;0.0220222137868404;0.0295578259974718;0.0247263666242361;3.81323552574031e-05;-0.0245587956160307;-0.0292639490216970;-0.0216563567519188;-0.00947358552366495;0.000167622143635526];
params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];

%% LRC parameters ---- This might be needed for data from phones
if flagLRC == 1
    params.doLRC=true;
    params.LRC_level=1;
    params.LRC_pyr_lvl = 5;
end

%% Post-process filtering parameters 
% if flagPost == 1
params.applyPostFilter=false; %  Sagi - false (no pst process)
params.guidedIter=0;

% Filter parameters - post processing parameters filter
params.fgsLambda=50;
params.fgsSigmaColor=3.5;
params.guidedDegreeOfSmoothing=0.0001*256*256;
params.guidedNeighborhoodSize=[3 3];
% end

% Unused parameters
%params.patchBorderSize=[16 16]; %
%params.isNormalized=false;% are the images pre-normalized?
%params.useSAD=true;
%params.normKerSize=7; % Sagi - spesific value for dual-PD (used for normalization)
%% Downscale data ---- This might be needed for data from phones
params.downSampleSize=[4 8]./2; % Sagi - Ignore - should be [1,1] (For google might need to tune)

