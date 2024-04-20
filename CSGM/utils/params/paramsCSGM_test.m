function params = paramsCSGM_test(val_gauss, val_range, val_conf,  val_p1, val_sig, val_lvl, val_priorW)
%%% Parameters - tested on one image (DSLR with GT) based on metric
% gaussKerSigma = 7 - smoothess result 

% params.dispRange = [-5 ,5]; - Guess depends on image, but too large is bad

% params.numIter = [3,3,3] ; % More than 3 iterations starts to fail (Very
% fast) Also looks like its best to have same number of iterations for each
% pyramid lvl.

% params.confidenceThresh= - 7e-2; % Threshold for flatting parabolas %%% Tune-flag
% params.penalty_border = 1000; % This is not that sensitive - just needs to be large 
% params.border_len = params.gaussKerSigma*2+max(abs(params.dispRange(:))); Not sensitive at all
% params.P1param = 2;  % 1.5-5 gives same result
% params.sigmaEdges=5; % Identifity edges (3-5 similar result)
% params.levels = 2-3;% More than that seemed to over-smooth (need to check different number of iterations for each lvl
% params.priorW = 0.05; % Larger values increase smoothness 
% params.interpolant = 'f5'; % Small increase with f2-f5 over f1
% biases+offsets - don't matter for metrics of aiws(1+2) + spearman corr


% Configurable parameters for algorithm 
flagLRC = 0; % LRC processing - for phone-cameras
flagPost = 0; % Post processing 

params.dispRange = [-val_range ,val_range]; % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
params.numIter = [3,3,3] ; % Number of iterations 

% params.subpix = ''; 'ENCC'; 

params.cost = 'SAD'; % Cost function type
params.gaussKerSigma = val_gauss;% STD of gaussian for average-weighting of score

params.confidenceThresh = val_conf; % Threshold for flatting parabolas %%% Tune-flag

% Border - refine
params.border_len = params.gaussKerSigma*2+max(abs(params.dispRange(:))); % Sagi - how much border is reliable 
params.penalty_border = 1000; % Reduce confidence on border

% Param - penalty on edges and large-diff
params.P1param = val_p1;  % Penalty for large diff in dispairty
params.sigmaEdges = val_sig; % Used to identify edges to give less penalty for smoothness
% params.currentToPrevWeight=0.2; % Not-used % SAgi - look at code % A large number means that we consider more the confidence of the previous sample so a "noiser" result with be created

% Pyramid values
params.plt_pyr = 0; % Plot each level of the pyramid
params.levels = val_lvl;% SGM pyramid levels
params.priorW = val_priorW; % This is prior when jumping from low level to higher level in the pyramid.
params.idx_pyr = 1; % First level of pyramid - this stays 1
 
%% ENNC - parameter
params.encc_window = [15,15];

%% Sub-pixel interpolation - This is based on simulations Sagi.K did to check pixel-locking
params.interpolant = 'f5';% ['f1', 'ENCC'] % Interpolation kernel for sub-pixel estimation
params.bias = [3.34275500790682e-05;0.00973533187061548;0.0220652986317873;0.0296091139316559;0.0247546602040529;2.39057189901359e-05;-0.0246540009975433;-0.0294469539076090;-0.0218638572841883;-0.00956899579614401;2.35768729908159e-05;0.00971676409244537;0.0220222137868404;0.0295578259974718;0.0247263666242361;3.81323552574031e-05;-0.0245587956160307;-0.0292639490216970;-0.0216563567519188;-0.00947358552366495;0.000167622143635526];
params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];

%% LRC parameters ---- This might be needed for data from phones
if flagLRC == 1
    params.doLRC=true;
    params.LRCLevel=4;
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

