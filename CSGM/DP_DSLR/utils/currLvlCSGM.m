function [dispar_map, sum_parab,conf_score_no_suprress] = currLvlCSGM(im_L, im_R, im_guide, params, prior_lvl_parab, prior_lvl_dispar )
% Calculate current level C-SGM

% Mark which level we are
flag_lvl = 1; % Current level in upscale
if ~isempty(prior_lvl_parab)
    flag_lvl = 2; % This is not the first lvl
end

%% This calculates score per pixel per disparity -> generate parabolas after sub pixel shifting
if flag_lvl == 1 % This is the first level of scales
    params.dispar_vals = params.dispRange(1)-1:params.dispRange(end)+1; % Add -+1 values to disparity range for calculating parabolas    
else  % This is the second level of scales
    params.dispar_vals = round(min(prior_lvl_dispar(:)))-2:ceil(max(prior_lvl_dispar(:))+2);
end

if ~strcmp(params.interpolant, 'ENCC')
    % Calculate cost and neighboors + confidence score + disparity integer value.
    [cost_neig, conf_score, dispar_int_val,conf_score_no_suprress] = disparCost(im_L, im_R, params);

    % Generate parab after sub-pixel estimation
    parab = genParab(cost_neig, dispar_int_val, params);
elseif strcmp(params.interpolant, 'ENCC')
    % Calculate cost and neighboors + confidence score + disparity integer value.
    [cost_neig, conf_score, dispar_int_val,conf_score_no_suprress] = disparCost(im_L, im_R, params);
    [parab,invalid_mask] = estSubPixENCC(cost_neig, dispar_int_val, im_L, im_R, params);
    conf_score(invalid_mask == 1) = 0.01^2; % Invalid conf
end

if 0
    dispar_map =-parab.b./(2.*parab.a);
    figure; imagesc(dispar_map)
end

% Multiply with confidence scores
parab.a = parab.a.*conf_score;
parab.b = parab.b.*conf_score;
parab.c = parab.c.*conf_score;    

% If we have data from prior from previous level of pyramid sacle the
% parabolas
if flag_lvl > 1
    parab.a = parab.a + params.priorW * prior_lvl_parab.a;
    parab.b = parab.b + params.priorW * prior_lvl_parab.b;
    parab.c = parab.c + params.priorW * prior_lvl_parab.c;
end

% Refine parabaols - low confidence + at borders
parab = refineParab(parab, params);

%%% If use cost that reduces size of image - we might need to reduce
%%% im_guide size

% Cost prop 8-way
sum_parab = propCSGM(parab, im_guide, params);

% Post-filter
if params.applyPostFilter
    display('Not implemented');
%     [sgmDisparity]=postFilter(sumA,sumB,img1forEdgesSampled,params,plt);
else
    dispar_map =-sum_parab.b./(2.*sum_parab.a);
end

% if strcmp(params.interpolant, 'ENCC')
%     dispar_map = disparity_res;
% end



end