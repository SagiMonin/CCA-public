function params = paramsCSGMStereo()
%% These are not used
params.dispRange = [-12,0]; % Disparity range - (this is for first level of pyramid) for DP this is symmetric 
params.cost = 'SAD'; % Cost function type
params.gaussKerSigma = 3;% STD of gaussian for average-weighting of score
%%
% params.numIter = [3 2 2 ] ; % Number of iterations 

% Confidence score for a 
params.confidenceThresh=-1e-2; % Threshold for flatting parabolas %%% Tune-flag
% This needs additional tuning

% Border - refine
params.border_len = round(params.gaussKerSigma*2+max(abs(params.dispRange(:)))); % Sagi - how much border is reliable 
params.penalty_border = 10; % Reduce confidence on border

% Param - penalty on edges and large-diff

params.P1param = 1; %  % Penalty for large diff in dispairty
% P1param - This is a weighted average of pixels - higher means we smooth
% it more
params.sigmaEdges = 3; % Used to identify edges to give less penalty for smoothness 
% Sigma - (Higher means we will only identify real strong edges, lower
% means we will find texture as edges)
% Propagation stops at strong edges

% Pyramid values
params.plt_pyr = 0; % Plot each level of the pyramid
params.levels = [1];% SGM pyramid levels
params.priorW = 0.8;0.045; % This is prior when jumping from low level to higher level in the pyramid.
params.idx_pyr = 1; % First level of pyramid - this stays 1
params.numIter = [1] ; % Number of iterations 

%% ENNC - parameter
params.encc_window = 15;

%% Sub-pixel interpolation - This is based on simulations Sagi.K did to check pixel-locking
params.interpolant = 'ENCC';% ['f1', 'ENCC'] % Interpolation kernel for sub-pixel estimation
params.bias = [3.34275500790682e-05;0.00973533187061548;0.0220652986317873;0.0296091139316559;0.0247546602040529;2.39057189901359e-05;-0.0246540009975433;-0.0294469539076090;-0.0218638572841883;-0.00956899579614401;2.35768729908159e-05;0.00971676409244537;0.0220222137868404;0.0295578259974718;0.0247263666242361;3.81323552574031e-05;-0.0245587956160307;-0.0292639490216970;-0.0216563567519188;-0.00947358552366495;0.000167622143635526];
params.offset = [-1,-0.900000000000000,-0.800000000000000,-0.700000000000000,-0.600000000000000,-0.500000000000000,-0.400000000000000,-0.300000000000000,-0.200000000000000,-0.100000000000000,0,0.100000000000000,0.200000000000000,0.300000000000000,0.400000000000000,0.500000000000000,0.600000000000000,0.700000000000000,0.800000000000000,0.900000000000000,1];
params.bias = zeros(size(params.offset)); [3.34275500790682e-05;0.00973533187061548;0.0220652986317873;0.0296091139316559;0.0247546602040529;2.39057189901359e-05;-0.0246540009975433;-0.0294469539076090;-0.0218638572841883;-0.00956899579614401;2.35768729908159e-05;0.00971676409244537;0.0220222137868404;0.0295578259974718;0.0247263666242361;3.81323552574031e-05;-0.0245587956160307;-0.0292639490216970;-0.0216563567519188;-0.00947358552366495;0.000167622143635526];

